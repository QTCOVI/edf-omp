c
c-----------------------------------------------------------------------
c
      subroutine oqsloc (aominp,mxcent,covx,epsneg,maxb,nloc,nloctarget,
     & udat,wfnfile,maxcritov,critov,ixcent,largwr,warn,skiph,lw,lerr)
c
c-----PARAMETERS--------------------------------------------------------
c
c     aominp()                        INPUT/OUTPUT
c
c     Atomic Overlap Matrix (AOM) of Canonical MOs in all the centers.
c
c-----AOMO()                          OUTPUT
c
c     AOM between orthonormalized Localized MOs in all the centers in the
c     output
c
c.....c()                            OUTPUT 
c
c     Rotation matrix giving the nonorthonormal localized MOs in terms of 
c     the input MOs: Loc_i = SUM_j C_ij CANMO_j
c
c.....mxcent                         INPUT
c
c     Maximum number of N in the search of (Nc,2e) bonds.
c
c.....covx                           INPUT
c
c     Factor multiplying the sum of covalent radii of two atoms A and B  
c     used in the connect routine to determined whether A and B are 
c     linked or not from a purely geomtrical point of view.
c
c-----epsneg                         INPUT
c
c     Diagonalizations in this routine involve definite positive 
c     matrices that should have all of their eigvenvalues positive.
c     Numerical errors, however, make that this not necessarily happens 
c     and some of the eigenvalues can be marginally negative. The value 
c     of 'epsneg' is the largest value (in absolute value) allowed in order
c     that the method ends nicely. For instance, if epsneg=-1D-6
c     and all the eigenvalues are greater that this number the 
c     process continues up to the end. However, if an eigvenvalue
c     is smaller that -1D-6, for instance -1D-5, this routine stops.
c
c-----maxb                           INPUT
c
c     Maximum number of links between two atoms A and B for which
c     a (2c,2e) localized MOs between these atoms is searched. Atoms
c     directly connected in a geometrical sence have a 'geometrical
c     distance' or link of 1, if A is connected to C which, in turn,
c     is connected to B the 'geometrical distance' or link between
c     A and B is 2, and so on. 
c
c.....nloc                           OUTPUT
c
c     Number of localized MOs found. It should be equal to nloctarget
c     (see below) if the routine ends nicely
c
c.....nloctarget                      INPUT
c
c     Half the number of electrons of the molecule (assumed to have
c     a even number of electrons, i.e. nalpha = nbeta)
c
c.....udat                           INPUT
c
c     Unit number of original WFN file
c
c-----wfnfile                        INPUT
c
c     Name of the WFN file
c
c-----maxcritov                      INPUT
c
c     Maximum order 'n' in the search of (nc,2e) localized MOs
c
c.....critov(1)...critov(mxcent)     INPUT
c
c     Tolerance used in localization: If an eigvenvalue of
c     S^(1/2) rho S^{1/2} is greater that critov(i) it means that the
c     associated eigenvector corresponds to an MO assumed to be loca-
c     lized in the 'i' centers involved in the diagonalization.
c
c.....ixcent(1)...ixcent(mxcent)     INPUT
c
c     Logical variable: If ixcent(i)=.true. (i-center,2e) bonds are
c     searched, otherwise (i-center,2e) bonds are skipped.
c
c-----largwr                         INPUT
c
c     Logical variable: if largwr=.true. a large output is requested,
c     otherwise a short output is used.
c
c-----warn                           OUTPUT
c
c     Returned  .false. if the routine ends nicely. Otherwise warn is
c     returned .true.
c
c-----skiph                          INPUT
c
c     Logical variable: If .true. hydrogen atoms are skipped in one-
c     center localizations. Otherwise hydrogens are not skipped.
c
c-----lw                             INPUT
c
c     Logical output unit
c
c-----lerr                           INPUT
c
c     Logical output unit of error messages
c
c-----------------------------------------------------------------------
c                         
      USE      space_for_wfncoef
      USE      space_for_wfnbasis
      USE      space_for_rdm1
      include 'implicit.inc'
      include 'constants.inc'
      include 'param.inc'
      include 'wfn.inc'
      include 'error.inc'
      include 'mline.inc'
      real   (kind=8) c(nmo,nmo)
      real   (kind=8) co(nmo,nmo)
      real   (kind=8) aom(ncent,nmo,nmo),aominp(ncent,nmo,nmo)
      real   (kind=8) aomold(ncent,nmo,nmo)
      real   (kind=8) aomo(ncent,nmo,nmo)
      real   (kind=8) sumaom(nloctarget,nloctarget)
      real   (kind=8) sg(nmo,nmo),eigen(nmo),uvec(nmo,nmo)
      real   (kind=8) uvecx(nloctarget,nloctarget)
      integer(kind=4) iab(2,nmo*(nmo-1)/2)
      real   (kind=8) xloci(nmo)
      integer(kind=4) bond(nmo)
      integer(kind=4) enlaces(ncent,ncent)
      real   (kind=8),allocatable, dimension (:,:)   :: rho
      real   (kind=8),allocatable, dimension (:)     :: eigs,work
      real   (kind=8),allocatable, dimension (:)     :: eigloc
      real   (kind=8),allocatable, dimension (:)     :: uveceig
      real   (kind=8),allocatable, dimension (:)     :: ocup
      real   (kind=8),allocatable, dimension (:,:)   :: v
      real   (kind=8),allocatable, dimension (:,:)   :: sloc,sloc0
      real   (kind=8),allocatable, dimension (:,:)   :: ctmp
      real   (kind=8),allocatable, dimension (:,:,:) :: aomx
      integer(kind=8),allocatable, dimension (:)     :: ntuples
      integer(kind=4),allocatable, dimension (:)     :: inside
      integer(kind=4),allocatable, dimension (:,:)   :: wh
      integer(kind=4),allocatable, dimension (:)     :: coord
      integer(kind=4),allocatable, dimension (:,:)   :: madis
      integer(kind=4),allocatable, dimension (:)     :: nmoc
      logical warn,warno,is_hydrogen,skiph
      parameter (maxcoord = 20)
      real   (kind=8), parameter :: fifty_percent = 50d0
c
c-----Data used in the orthonormalization process
c
      real(kind=8),    allocatable,dimension (:,:)   :: cc
      real(kind=8),    allocatable,dimension (:,:)   :: cco
      real(kind=8),    allocatable,dimension (:,:)   :: ornor
c
      real   (kind=8)       epsneg
      real   (kind=8)       critov(maxcritov)
      logical               ixcent(maxcritov)
      integer (kind=4)      p,q
      logical               largwr,dow,first
      integer(kind=4)       udat,udatnw
      character(len=*)      wfnfile
      character (len=mline) wfnloc
      character (len=mline) locform(nmo)
      character (len=30)    typebond

      character*1    digs(0:9)
      data digs / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
c
c-----------------------------------------------------------------------
c
      call timer (2,iopenloc,'_oqsloc   ',-1)

       aom = aominp
Cc
Cc-----Read 1-RDM in the canonical MO basis
Cc
C      lu1 = 55
C      open (lu1,file=trim(wfnfile)//".1rdm",form='unformatted',
C     &      iostat=ios,status='old')
C      if (ios.ne.0) write (lerr,100) trim(wfnfile)//".1rdm"
C      if (ios.ne.0) stop
C      read (lu1) nmox
C      if (nmox.ne.nmo) then
C        write (lerr,1123) trim(wfnfile)//".1rdm",nmox,nmo
C        stop
C      endif
C      allocate (rho(nmo,nmo))
C      rho = 0d0
C      do
C        read (lu1,iostat=ios) m,n,value
C        if (ios.ne.0) exit
C        rho(m,n) = value
C        rho(n,m) = value
C      enddo
C
      allocate (rho(nmo,nmo))
      rho(1:nmo,1:nmo) = c1et(1:nmo,1:nmo)
c
c-----This order means that the routine is only applicable when the
c     number of ALPHA and BETA MOs is the same.
c
      rho = rho * 0.5d0
c
c.....Allocate arrays
c
      allocate (eigs(nmo)          )
      allocate (eigloc(nmo)        )
      allocate (ntuples(mxcent)    )
      allocate (nmoc(mxcent)       )
      allocate (inside(ncent)      )
      allocate (work(3*nmo-1)      )
      allocate (v(nmo,nmo)         )
      allocate (ctmp(nmo,nmo)      )
      allocate (madis(ncent,ncent) )
      allocate (coord(ncent)       )
      allocate (wh(ncent,maxcoord) )
      allocate (uveceig(nmo)       )
      ncenteff = min(mxcent,ncent)
      write (lw,10) ' B E G I N  '
      write (lw,21) ncenteff,(critov(j),j=1,ncenteff)
      if (     skiph) write (lw,51) '     '
      if (.not.skiph) write (lw,51) ' NOT '
      write (lw,52) maxb
c
c-----Obtain coordinations, connectivities, and distance arrays
c
      call connect (lw,covx,ncent,maxcoord,largwr,nbonds,coord,wh,madis)
c
c-----Largest number of bonds between two 'connected' atoms
c
      longest_chain = maxval(madis)
c
c-----Print the input AOM arrays
c
      if (largwr) then
        write (lw,*) '# Atomic Overlap matrix (AOM)'
        write (lw,*) '#'
        do i=1,ncent
          write (lw,*) '     ATOM ',i
          do j=1,nmo
            write (lw,133) (aom(i,j,k),k=1,j)
          enddo
        enddo
        write (lw,*) '#'
        write (lw,*) '#    TOTAL'
        do j=1,nmo
          write (lw,133) (sum(aom(:,j,k)),k=1,j)
        enddo
      endif
c
c.....Start the process.
c
      natoms = 1
      do while (natoms.le.ncenteff) 
        xtuples = 1d0
        do i = 0,natoms-1
          xtuples = xtuples*dble(ncent-i)/dble(natoms-i)
        enddo
        ntuples(natoms) = nint(xtuples)
        natoms = natoms+1
      enddo     
      allocate (sloc(nmo,nmo)         )
      allocate (sloc0(nmo,nmo)        )
      sloc    = 0d0
      sloc0   = 0d0
      ctmp    = 0d0
      n3      = 3*nmo-1
      nloc    = 0
      natoms  = 1
      nmoc    = 0
      warn    = .false.
      nbab    = 0
      enlaces = 0
      write (lw,28) epsneg
      do while (natoms.le.ncenteff.and.nloc.lt.nmo)
        if (ixcent(natoms)) then
          write (lw,11) natoms
          first = .true.
          lis = 1
          dow = .true.
          do while (lis.le.ntuples(natoms).and.nloc.lt.nmo)
            call ngen (natoms,ncent,ncent,inside,first)
            if (natoms.eq.1) then
              is_hydrogen = nint(charge(inside(1))).eq.1
              if (is_hydrogen.and.skiph) goto 111
            elseif (natoms.eq.2) then
              iat1=inside(1)
              iat2=inside(2)
              m12=madis(iat1,iat2)
              if (m12.lt.1.and.m12.gt.maxb) goto 111
            elseif (natoms.ge.3) then
              ic1=coord(inside(1))
              ic2=coord(inside(2))
              ic3=coord(inside(3))
              if (ic1.le.2.and.ic2.le.2.and.ic3.le.2) goto 111
            endif
c
c-----------Group overlap integrals are obtained from the AOM
c
            sg = 0d0
            do i=1,natoms
              sg = sg+aom(inside(i),:,:)
            enddo
            call dsyev ('V','U',nmo,sg,nmo,eigen(1:nmo),work,n3,info)
            call errdsyev (info,'oqsloc.f')
            uvec = sg
c
c-----------Test that sg() is definite positive
c
            do m=1,nmo
              if (eigen(m).gt.0d0) then
              elseif (eigen(m).gt.epsneg) then
                write (lerr,222) eigen(m),inside(1:natoms)
                eigen(m)=0d0
              else
                write (lerr,1124) m,eigen(m)
                stop
              endif
            enddo
c
c-----------Compute S^{1/2} (stored in sg())
c
            sg = 0d0 
            do k=1,nmo
              do l=1,nmo
                forall (m=1:nmo) uveceig(m)=uvec(l,m)*sqrt(eigen(m))
                sg(k,l) = dot_product(uvec(k,:),uveceig(:))
              enddo
            enddo
c
c-----------This is S^{1/2} * 1-RDM * S^{1/2}. The previous "localized"
c           density is substracted
c
            v = matmul(sg,matmul(rho,transpose(sg))) - sloc0
            eigs = 0d0
            call dsyev ('V','U',nmo,v,nmo,eigs(1:nmo),work,n3,info)
            call errdsyev (info,'oqsloc.f')
            if (largwr) then
              write (lw,*) '#'
              write (lw,432) (inside(i),i=1,natoms)
              write (lw,*) '#'
              write (lw,133) (eigs(k),k=1,nmo)
            endif
c
c-----------Be carefull with the following loop. It only works if the 
c           used diagonalization routine (in this case, 'dsyev') 
c           returns the eigenvalues in ascencing order. Otherwise
c           before this loop eigenvalues and eigvectors should be
c           reordered such that eigs(1) is the smallest eigvenvalue
c           eigs(2) the following, and so on.
c
            do i = nmo,1,-1
              if (eigs(i).lt.critov(natoms)) exit
              nloc = nloc + 1
              if (natoms.eq.2) then
                nbab=nbab+1
                if (nbab.gt.nmo*(nmo-1)/2) then
                  stop '# uhfoqs.f: Increase the first dimension of iab'
                endif
                iab(1:2,nbab)=inside(1:2)
                bond(nbab)=nloc
              endif
              nmoc(natoms) = nmoc(natoms)+1
              if (nloc.gt.nmo) write (lw,433)
              if (nloc.gt.nmo) stop
              eigloc(nloc) = eigs(i)
              ctmp(:,nloc) = v(:,i)
              do j=1,nmo
                do k=1,nmo
                  prod = ctmp(j,nloc)*ctmp(k,nloc)*eigloc(nloc)
                  sloc(j,k) = sloc(j,k)+prod
                enddo
              enddo
              if (dow) then
                write (lw,113) natoms,critov(natoms)
                dow = .false.
              endif
              if (nloc.le.nloctarget) then
                write (lw,116) nloc,eigs(i),
     &            (atnam(inside(j))(3:4),inside(j),j=1,natoms)
                if (natoms.eq.2) then
                  i1=min(inside(1),inside(2))
                  i2=max(inside(1),inside(2))
                  enlaces(i1,i2)=enlaces(i1,i2)+1
                endif
              else
                write (locform(nloc-nloctarget),116) nloc,eigs(i),
     &            (atnam(inside(j))(3:4),inside(j),j=1,natoms)
              endif
            enddo
 111        continue
            lis=lis+1
          enddo
        endif
        sloc0 = sloc
        natoms = natoms+1 
      enddo
c
      if (nloc.lt.nloctarget) then
        write (lw,434)
        warn = .true.
        return
      elseif (nloc.eq.nloctarget) then
        write (lw,300)
        write (lw,301) (nmoc(i),i=1,mxcent)
      else
        write (lw,435)
        do i=nloctarget+1,nloc
          write (lw,'(a)') trim(locform(i-nloctarget))
        enddo
        warn = .true.
        return
      endif

      c(1:nloc,:) = transpose(ctmp(:,1:nloc))
c
c-----Each OPENLOC localized MO (LMO) is normalized to 1.0 in the
c     fragment it was obtained, but is normalized to eigloc^{-1} in R^3. 
c     For this reason each LMO will be multiplied by sqrt(eigloc) such
c     that it will be normalized in R^3 instead of the fragment.
c
***** forall (i=1:nloc) c(i,1:n)=c(i,1:n)*sqrt(eigloc(i))
c
c-----The above is true but the diagonalization routine already returns
c     normalized eigenvectors. This is why we comment the line 'forall'
c
      if (largwr) then
        write (lw,20)
        do i=1,nloc
          write (lw,33) i
          write (lw,30) (c(i,j),j=1,nmo)
        enddo
      endif
c
c-----------------------------------------------------------------------
c                                                                      !
c-----Compute AOM between the new nonorthonormal localized MOs         !
c                                                                      !
      aomold = aom                                                     !
      do i=1,ncent                                                     !
        aom(i,1:nloc,1:nloc) = matmul(c(1:nloc,:),                     !
     &   matmul(aomold(i,:,:),transpose(c(1:nloc,:))))                 !
      enddo                                                            !
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----Write diagonal AOM elements in each atom or fragment.            !
c                                                                      !
      write (lw,68) 'non-orthonormal'                                  !
      write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'   !
      do ip=1,nloc                                                     !
        write (lw,410,advance='no') ip                                 !
        if (ncent.gt.10) then                                          !
          write (lw,41 ) (100d0*aom(k,ip,ip),k=1,10)                   !
          write (lw,411) (100d0*aom(k,ip,ip),k=11,ncent)               !
          write (lw,*)                                                 !
        else                                                           !
          write (lw,41) (100d0*aom(k,ip,ip),k=1,ncent)                 !
        endif                                                          !
      enddo                                                            !
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----effective number of atoms expanded by each localized MO          !
c                                                                      !
      do ip=1,nloc                                                     !
        xloci(ip)=0d0                                                  !
        do k = 1, ncent                                                !
          spp = aom(k,ip,ip)                                           !
          xloci(ip) = xloci(ip) + spp * spp                            !
        enddo                                                          !
      enddo                                                            !
      write (lw,520) 'non-orthonormal',(1d0/xloci(ip),ip=1,nloc)       !
c
c-----Print type of bond between each pair of atoms
c
      write (lw,320)
 320  format (' # Type of bond between atoms',/,' # ',70('-'))
      do i=1,ncent-1
        do j=i+1,ncent
          deltais=0d0
          if (enlaces(i,j).gt.0) then
            do k=1,nbab
              iabmin=min(iab(1,k),iab(2,k))
              iabmax=max(iab(1,k),iab(2,k))
              ib=bond(k)
              di=4d0*aom(iabmin,ib,ib)*aom(iabmax,ib,ib)
              if (iabmin.eq.i.and.iabmax.eq.j) then
                deltais=deltais+di
              endif
            enddo
          endif
          if (enlaces(i,j).eq.0) then
            typebond = 'NO BOND                       '
          elseif (enlaces(i,j).eq.1) then
            typebond = 'SINGLE BOND                   '
          elseif (enlaces(i,j).eq.2) then
            typebond = 'DOUBLE BOND                   '
          elseif (enlaces(i,j).eq.3) then
            typebond = 'TRIPLE BOND                   '
          elseif (enlaces(i,j).eq.4) then
            typebond = 'QUADRUPLE BOND                '
          elseif (enlaces(i,j).eq.5) then
            typebond = 'QUINTUPLE BOND                '
          else
            typebond = 'SEXTUPLE or higher order BOND '
          endif
          write (lw,321) i,j,typebond,deltais
        enddo
      enddo
 321  format (' # Bond type atoms ',2I4,' : ',a16,'  Delta = ',F15.8)
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----Write approximate delocalization indices                         !
c                                                                      !
      write (lw,331)                                                   !
      do i=1,nbab                                                      !
        ib=bond(i)                                                     !
        di = 4d0 * aom(iab(1,i),ib,ib) * aom(iab(2,i),ib,ib)           !
        write (lw,65) iab(1:2,i),di                                    !
      enddo                                                            !
 65   format (' # Atoms ',2I3,'   Delta = ',F15.8)                     !
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----Orthogonalize the nonorthonormal Localized MOs by Lowdin. This   !
c     is equivalent to the maximum overlap method. First compute the   !
c     overlap in R^3 of the nonorthonormal Localized MOs.              !
c                                                                      !
c-----First, obtain S, where S is the overlap matrix in R^3 betwen     !
c     the just found nonorthonormal Localized MOs.                     !
c                                                                      !
c-----Overlaps in R^3 of the nonorthormal OQS MOs                      !
c                                                                      !
      sumaom=0d0                                                       !
      do k=1,ncent                                                     !
        sumaom(:,:)=sumaom(:,:)+aom(k,1:nloc,1:nloc)                   !
      enddo                                                            !
      if (largwr) then                                                 !
        write (lw,*) '#'                                               !
        write (lw,*) '# Overlaps in R^3 of the nonorthormal OQS MOs'   !
        do i=1,nloc                                                    !
          write (lw,133) (sumaom(i,j),j=1,i)                           !
        enddo                                                          !
      endif                                                            !
c                                                                      !
      n3=3*nloc-1                                                      !
      uvecx = sumaom                                                   !
c                                                                      !
c-----Diagonalize S                                                    !
c                                                                      !
      call dsyev ('V','U',nloc,uvecx,nloc,eigen(1:nloc),work,n3,info)  !
      call errdsyev (info,'oqsloc.f')                                  !
c                                                                      !
c-----Obtain  S^{-1/2} (stored in sumaom())                            !
c                                                                      !
      sumaom = 0d0                                                     !
      do i=1,nloc                                                      !
        do j=1,nloc                                                    !
          value=0d0                                                    !
          do k=1,nloc                                                  !
            value=value+uvecx(i,k)*uvecx(j,k)/sqrt(eigen(k))           !
          enddo                                                        !
          sumaom(i,j) = value                                          !
        enddo                                                          !
      enddo                                                            !
      co(1:nmo,1:nloc)=matmul(transpose(c(1:nloc,1:nmo)),sumaom(:,:))  !
      co(1:nloc,1:nmo)=transpose(co(1:nmo,1:nloc))                     !

      do i=1,ncent                                                     !
        aomo(i,1:nloc,1:nloc)=matmul(transpose(sumaom),                !
     &  matmul(aom(i,1:nloc,1:nloc),sumaom))                           !
      enddo                                                            !
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c                                                                      !
c-----Write orthonormalized MOs in terms of the canonical MOs          !
c                                                                      !
      if (largwr) then                                                 !
        write (lw,2001)                                                !
        do i=1,nloc                                                    !
          write (lw,33) i                                              !
          write (lw,30) (co(i,j),j=1,nmo)                              !
        enddo                                                          !
      endif                                                            !
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----Write diagonal AOM elements in each atom or fragment between     !
c     orthonormalized localized MOs                                    !
c                                                                      !
      write (lw,68) 'orthonormal'                                      !
      write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'   !
      do ip=1,nloc                                                     !
        write (lw,410,advance='no') ip                                 !
        if (ncent.gt.10) then                                          !
          write (lw,41)  (100d0*aomo(k,ip,ip),k=1,10)                  !
          write (lw,411) (100d0*aomo(k,ip,ip),k=11,ncent)              !
          write (lw,*)                                                 !
        else                                                           !
          write (lw,41) (100d0*aomo(k,ip,ip),k=1,ncent)                !
        endif                                                          !
      enddo                                                            !
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----effective number of atoms expanded by each orthonormal MO        !
c                                                                      !
      do ip=1,nloc                                                     !
        xloci(ip)=0d0                                                  !
        do k = 1, ncent                                                !
          spp = aomo(k,ip,ip)                                          !
          xloci(ip) = xloci(ip) + spp * spp                            !
        enddo                                                          !
      enddo                                                            !
      write (lw,520) 'orthonormalized',(1d0/xloci(ip),ip=1,nloc)       !
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----Overlaps in R^3 of the orthonormalized OQS MOs                   !
c                                                                      !
      sumaom=0d0                                                       !
      do k=1,ncent                                                     !
        sumaom=sumaom+aomo(k,1:nloc,1:nloc)                            !
      enddo                                                            !
      if (largwr) then                                                 !
        write (lw,*) '#'                                               !
        write (lw,*) '# Overlaps in R^3 of the orthonormalized OQS MOs'!
        do i=1,nloc                                                    !
          write (lw,133) (sumaom(i,j),j=1,i)                           !
        enddo                                                          !
      endif                                                            !
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----Find localized MOs in terms of the primitive basis               !
c                                                                      !
      allocate (cc(nloc,nprims))                                       !
      allocate (cco(nloc,nprims))                                      !
      do i=1,nloc                                                      !
        do j=1,nprims                                                  !
          value1=0d0                                                   !
          value2=0d0                                                   !
          do k=1,nmo                                                   !
            ik=nmo+k                                                   !
            value1=value1+coef(ik,j)*c(i,k)                            !
            value2=value2+coef(ik,j)*co(i,k)                           !
          enddo                                                        !
          cc(i,j)  = value1                                            !
          cco(i,j) = value2                                            !
        enddo                                                          !
      enddo                                                            !
c                                                                      !
c-----------------------------------------------------------------------
c
c-----Write AOM between the non-orthonormal localized MOs
c
      if (largwr) then
        write (lw,115) 'Non-orthonormal'
        do k=1,ncent
          write (lw,*) '# Center',k
          do i=1,nloc
            write (lw,101) (aom(k,i,j),j=1,i)
          enddo
        enddo
      endif

      if (largwr) then
        write (lw,115) 'orthonormal'
        do k=1,ncent
          write (lw,*) '# Center',k
          do i=1,nloc
            write (lw,101) (aomo(k,i,j),j=1,i)
          enddo
        enddo
      endif
c
c-----Compute the overlaps in R^3 between non-orthonormal and
c     orthonormal localized MOs
c
      allocate (ornor(nloc,nloc))
      ornor=matmul(c(1:nloc,1:nmo),transpose(co(1:nloc,1:nmo)))
      write (lw,200) 'non-orthonormal and orthonormal MOs'
      do i=1,nloc
        write (lw,101)(ornor(i,j),j=1,nloc)
      enddo
c
c-----write a WFN file with non-orthonormal localized MOs
c
      allocate (ocup(nloc),stat=ier)
      if (ier.ne.0) stop '# oqsloc.f: Cannot allocate ocup()'
      ocup(1:nloc) = 2d0
      inpr = nprims
      wfnloc = trim(wfnfile)//"-open"
      udatnw = udat+10
      open (unit=udatnw,file=trim(wfnloc),status='unknown')
      write (lw,*)
      write (lw,351) trim(wfnloc)
      call cdafhmos (udat,udatnw,0,ocup,cc,nloc,inpr)
c
c-----write a WFN file with orthonormal localized MOs
c
      inpr = nprims
      wfnloc = trim(wfnfile)//"-open-ortho"
      udatnw = udat+10
      open (unit=udatnw,file=trim(wfnloc),status='unknown')
      write (lw,352) trim(wfnloc)
      call cdafhmos (udat,udatnw,0,ocup,cco,nloc,inpr)
      deallocate (ocup,stat=ier)
      if (ier.ne.0) stop '# oqsloc.f: Cannot deallocate ocup()'
c
c-----Write file with AOM between the non-orthonormal localized MOs
c     and between the orthonormal localized MOs
c
      lu18 = udat+11
      wfnloc = trim(wfnfile)//"-open.aom"
      open (unit=lu18,file=trim(wfnloc),status='unknown')
      write (lw,353) trim(wfnloc)
      do k=1,ncent
        write (lu18,'(I4,a)') k,' <=== AOM within this center'
        write (lu18,80) ((aom(k,m,j),m=1,j),j=1,nloc)
      enddo
      close (unit=lu18)
      wfnloc = trim(wfnfile)//"-open-ortho.aom"
      open (unit=lu18,file=trim(wfnloc),status='unknown')
      write (lw,354) trim(wfnloc)
      do k=1,ncent
        write (lu18,'(I4,a)') k,' <--- AOM in this center'
        write (lu18,80) ((aomo(k,m,j),m=1,j),j=1,nloc)
      enddo
      close (unit=lu18)
c
c-----Analyze the MOs and their localizations in all the centers
c
      critics=0.02d0
      allocate (aomx(ncent,nloc,nloc))
      aomx=aom(1:ncent,1:nloc,1:nloc)
      call aomtest (aomx,nloc,ncent,critics,lw,largwr)
      write (lw,3340) 'Open quantum-systems localized MOs'
      call newdrvmult (nloc,nprims,ncent,cc,lw)

      write (lw,10) '  E N D  '

      deallocate (rho      )
      deallocate (eigs     )
      deallocate (ntuples  )
      deallocate (work     )
      deallocate (v        )
      deallocate (uveceig  )
      deallocate (sloc     )
      deallocate (sloc0    )
      deallocate (ctmp     )
      deallocate (inside   )
      deallocate (eigloc   )
      deallocate (madis    )
      deallocate (nmoc     )
      deallocate (coord    )
      deallocate (wh       )
      deallocate (cc       )
      deallocate (cco      )
      deallocate (ornor    )
      deallocate (aomx     )
c
      return
c
      call timer (4,iopenloc,'_oqsloc   ',-1)
c
c.....Formats
c
 410  format (' # ',I4)
 20   format (' # Rotation Matrix (READ BY ROWS)',6x,
     & 'OQS MO_i =  SUM_j CANMO_j C_ij',/,' # ',
     & 'OQS MOs have been renormalized such that <i|i>_R^3=1')
 2001 format (' # Rotation Matrix (READ BY ROWS)',6x,
     & 'Orthonormal OQS MO_i =  SUM_j CANMO_j C_ij')
 30   format (5(1x,E15.8))
 41   format (1000(10x,10(1x,F5.1,'%')))
 411  format (1000(17x,10(1x,F5.1,'%'),/))
 300  format (' #',/,' # ',70('-'),/,
     &  ' # CONVERGED (1--N)-CENTER OPEN SYSTEM ANALYSIS',/,' #')
 66   format (/' # Diagonal overlaps of each MO in each atom')
 68   format (/' # Localization of each ',a,
     &  ' OQS MO in each atom')
 33   format (' # MO ',I3)
 113  format (/,1x,'# ',78('-'),/,2x,I2,
     & '-center OQS analysis, Critical overlap = ',F12.5,/,
     & ' # ',78('-'))
 116  format (4x,'MO ',I4,' with eigenvalue = ',F12.5,
     & ') localized on atom(s)',10(1x,'(',A2,I2,')'))
 21   format (' # CRITOV (1 ..',I2,')  = ',10(1x,F5.3))
 10   format (/,' # ',70('-'),/,' # ',a,
     & '(1--N)-CENTER OPEN SYSTEM LOCALIZATION ANALYSIS',/, 
     & ' # ',70('-'))
 11   format (' #',/,
     & ' # OQS localization analysis on ',I2,' centers',/,' #')
 69   format (' #',/,' # SUM of eigenvalues = ',F18.8,/,' #')
 70   format (' # ',7x,I4,'-center = ',F18.8)
 71   format (' # ',/,
     &  ' # Electrons associated to successful diagonalizations')
 710  format (' # Atom ',I4,'   Electrons = ',F18.8)
 711  format (' # SUM  ',17x,'= ',F18.8)
 334  format (/,' # Second & fourth moment orbital spread of ',a,/,3x,
     &'I.M. Hoyvik & P. Jorgensen [Chem. Rev. 116, 3306-3326 (2016)]',
     &/,3x,
     & 'I.M. Hoyvik, B. Jansik & P. Jorgensen [JCP 137, 2224114 (2012)'
     & /,3x,61('-'),/,4x,'Orbital',1x,'sigma_2^p',5x,'sigma_4^p',
     & 5x,'beta_p',8x,'beta_p^4',/,3x,61('-'))
 351  format (1x,
     &   "# Writing file '",a,"'",1x,"with non-orthonormal MOs")
 352  format (1x,
     &   "# Writing file '",a,"'",1x,"with orthonormal MOs")
 353  format (1x,"# Writing file '",a,"'",1x,
     &  "with AOM for non-orthonormal MOs")
 354  format (1x,"# Writing file '",a,"'",1x,
     &  "with AOM for orthonormal MOs")
 1000 format (
     & ' # OQS Localized MOs are orthonormalized with Lowdin')
 112  format (2(' #',/),' #',90('-'),/,' # ',a,
     & ' OF RESULTS WITH ORTHONORMALIZED OPTIMIZED MOs',/,' #',90('-'))
 51   format (' # In single atom analyses, H atoms are',a,'skipped')
 52   format (
     & " # In 'atom1-atom2' searches, maximum distance index is ",I2)
 432  format (
     & ' # Eigenvalues of the OQS of atoms ',20(1x,I3,'+'))
 433  format (' #',/,
     & ' # Too many OQS Localized MOs',/,
     & ' # Try again decreasing the critov() values',/,' #')
 434  format (' #',/,' # oqsloc.f: Unable to localize all the MOs',
     & /,' #',/,' # NON CONVERGED (1--N)-CENTER OPEN SYSTEM ANALYSIS',
     & /,' #')
 435  format (' #',
     & /,' # oqsloc.f: ! WARNING: Too many localized MOs !',/,' #')
 100  format (' # oqsloc.f: File ',a,' NOT FOUND')
 1123 format (' # NMO in file ',a,'(',I4,') # Current value = ',I4)
 1124 format (' # Eigenvalue ',I4,' of sg() is negative, ',E17.10)
 28   format (' # The value of EPSNEG is ',E15.8)
 222  format (' # Negative eigenvalue ',E15.8,
     &    ' set to 0.0 diagonalizing fragment ',100I4)
 301  format (' # Number of MOs involving 1, 2, 3, ... centers',/,
     &        ' # ',20(1x,I3))
 115  format (' # AOM between ',a,' localized MOs')
 101  format (1000(' # ',5(1x,F15.8),/))
 133  format (5(1x,F15.8))
 200  format (' # Overlaps in R^3 between ',a,/,
     & ' # Non-orthonormal(CHANGING DOWN)\Orthonormal(CHANGING RIGHT)')
 201  format (' # Oxidation state of atom (',a4,I4,' ) is ',I3)
 80   format (6(1x,e16.10))
 330  format (' #',/,
     &        ' # -------------------------------------------',/,
     &        ' # Oxidation states of atoms',/,
     &        ' # -------------------------------------------',/,
     &        ' #')
 331  format (/,' # Counting bonds and their delocalization indices',/,
     &          ' # -----------------------------------------------')
 3340 format (/,' # Second & fourth moment orbital spread of ',a)
 520  format (' # Atoms expanded by each ',a,' OQS MO',/,
     &    1000(5(3x,F12.5),/))
      end
