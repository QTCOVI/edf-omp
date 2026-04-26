c
c-----------------------------------------------------------------------
c
      subroutine uhfoqsloc (aom,aomo,rho,c,co,icode,nab,
     & nmotarget,mxcent,covx,epsneg,maxb,nloc,udat,wfnfile,
     & maxcritov,critov,ixcent,largwr,warn,skiph,lw,lerr)
c
c-----PARAMETERS--------------------------------------------------------
c
c-----AOM()                          INPUT/OUTPUT
c
c     Atomic Overlap Matrix (AOM) between Canonical MOs in all the centers
c     in the input and AOM matrix betweel Localized MOs in the output.
c
c-----AOMO()                          OUTPUT
c
c     AOM between orthonormalized Localized MOs in all the centers in the
c     output
c
c-----RHO()                          INPUT
c
c     First order RDM of alpha (icode=1) or beta (icode=2) electrons
c
c-----C()                            OUTPUT 
c
c     Rotation matrix giving the nonorthonormal localized MOs in terms of 
c     the input MOs: Loc_i = SUM_j C_ij CANMO_j
c
c-----CO()                           OUTPUT 
c
c     Rotation matrix giving the orthonormal localized MOs in terms of 
c     the input MOs: Loc_i = SUM_j C_ij CANMO_j
c
c-----ICODE                          INPUT 
c
c     = 1 ---> alpha block
c     = 2 ---> beat  block
c
c-----NAB                            INPUT
c
c     Number of MOs involved in the localization  (equal to the number
c     of alpha or beta MOs in the full set of NMO orbitals of the WFN
c     file)
c
c-----NMOTARGET                      INPUT
c
c     Expected number of Localized MOs that should be found in case 
c     that the method works find
c
c.....MXCENT                         INPUT
c
c     Maximum number of N in the search of (Nc,2e) bonds.
c
c.....COVX                           INPUT
c
c     Factor multiplying the sum of covalent radii of two atoms A and B  
c     used in the connect routine to determined whether A and B are 
c     linked or not from a purely geomtrical point of view.
c
c-----EPSNEG                         INPUT
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
c-----MAXB                           INPUT
c
c     Maximum number of links between two atoms A and B for which
c     a (2c,2e) localized MOs between these atoms is searched. Atoms
c     directly connected in a geometrical sence have a 'geometrical
c     distance' or link of 1, if A is connected to C which, in turn,
c     is connected to B the 'geometrical distance' or link between
c     A and B is 2, and so on. 
c
c.....NLOC                           OUTPUT
c
c     Number of localized MOs found. It should be equal to nmotarget
c     (see below) if the routine ends nicely
c
c.....UDAT                           INPUT
c
c     Unit number of original WFN file
c
c-----WFNFILE                        INPUT
c
c     Name of the WFN file
c
c-----MAXCRITOV                      INPUT
c
c     Maximum coordination number assumed for an atom in the 'connect'
c     routine
c
c.....CRITOV(1)...CRITOV(MXCENT)     INPUT
c
c     Tolerance used in localization: If an eigvenvalue of
c     S^(1/2) rho S^{1/2} is greater that critov(i) it means that the
c     associated eigenvector corresponds to an MO assumed to be loca-
c     lized in the 'i' centers involved in the diagonalization.
c
c.....IXCENT(1)...IXCENT(MXCENT)     INPUT
c
c     Logical variable: If ixcent(i)=.true. (i-center,2e) bonds are
c     searched, otherwise (i-center,2e) bonds are skipped.
c
c-----LARGWR                         INPUT
c
c     Logical variable: if largwr=.true. a large output is requested,
c     otherwise a short output is used.
c
c-----WARN                           OUTPUT
c
c     Returned  .false. if the routine ends nicely. Otherwise warn is
c     returned .true.
c
c-----skiph                          INPUT
c
c     Logical variable: If .true. hydrogen atoms are skipped in one-
c     center localizations. Otherwise hydrogens are not skipped.
c
c-----LW                             INPUT
c
c     Logical output unit
c
c-----LERR                           INPUT
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
      include 'corr.inc'
      include 'error.inc'
      include 'mline.inc'
      real   (kind=8) c(nab,nab)
      real   (kind=8) co(nab,nab)
      real   (kind=8) rho(nab,nab)
      real   (kind=8) aom(ncent,nab,nab)
      real   (kind=8) aomold(ncent,nab,nab)
      real   (kind=8) aomo(ncent,nab,nab)
      real   (kind=8) sumaom(nmotarget,nmotarget)
      real   (kind=8) sg(nab,nab),eigen(nab),uvec(nab,nab)
      integer(kind=4) iab(2,nmo*(nmo-1)/2)
      real   (kind=8) xloci(nmo)
      integer(kind=4) bond(nmo)
      integer(kind=4) enlaces(ncent,ncent)
      real   (kind=8),allocatable, dimension (:)     :: eigs,work
      real   (kind=8),allocatable, dimension (:)     :: eigloc
      real   (kind=8),allocatable, dimension (:)     :: uveceig
      real   (kind=8),allocatable, dimension (:,:)   :: v
      real   (kind=8),allocatable, dimension (:,:)   :: sloc,sloc0
      real   (kind=8),allocatable, dimension (:,:)   :: ctmp
      integer(kind=8),allocatable, dimension (:)     :: ntuples
      integer(kind=4),allocatable, dimension (:)     :: inside
      integer(kind=4),allocatable, dimension (:,:)   :: wh
      integer(kind=4),allocatable, dimension (:)     :: coord
      integer(kind=4),allocatable, dimension (:,:)   :: madis
      integer(kind=4),allocatable, dimension (:)     :: nmoc
      logical warn,warno,is_hydrogen,skiph
      parameter (maxcoord = 20)
c
      real   (kind=8)       epsneg
      real   (kind=8)       critov(maxcritov)
      logical               ixcent(maxcritov)
      integer (kind=4)      p,q
      logical               largwr,dow,esta,first
      integer(kind=4)       udat,udatnw
      character(len=*)      wfnfile
      character (len=mline) wfnloc
      character (len=mline) locform(nab)
      character (len=30)    typebond
      character (len=5)     spin

      character*1    digs(0:9)
      data digs / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
c
c-----------------------------------------------------------------------
c
      call timer (2,iopenloc,'_uhfoqsloc',-1)
c
c.....Allocate arrays
c
      n = nab
      allocate (eigs(nab)          )
      allocate (eigloc(nab)        )
      allocate (ntuples(mxcent)    )
      allocate (nmoc(mxcent)       )
      allocate (inside(ncent)      )
      allocate (work(3*nab-1)      )
      allocate (v(nab,nab)         )
      allocate (ctmp(nab,nab)      )
      allocate (madis(ncent,ncent) )
      allocate (coord(ncent)       )
      allocate (wh(ncent,maxcoord))
      allocate (uveceig(nab)       )
      ncenteff = min(mxcent,ncent)
      if (icode.eq.1) write (lw,10) ' B E G I N  ','ALPHA block'
      if (icode.eq.2) write (lw,10) ' B E G I N  ','BETA  block'
      write (lw,21) ncenteff,(critov(j),j=1,ncenteff)
      if (     skiph) write (lw,51) '     '
      if (.not.skiph) write (lw,51) ' NOT '
      write (lw,52) maxb
c
c-----Obtain coordinations, connectivities, and distance arrays
c
      call connect (lw,covx,ncent,maxcoord,largwr,nbonds,
     &  coord,wh,madis)
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
      allocate (sloc(nab,nab)         )
      allocate (sloc0(nab,nab)        )
      sloc    = zero
      sloc0   = zero
      ctmp    = zero
      n3      = 3*nab-1
      nloc    = 0
      natoms  = 1
      nmoc    = izero
      warn    = .false.
      nbab    = 0
      enlaces = 0
      write (lw,28) epsneg
      do while (natoms.le.ncenteff.and.nloc.lt.nab)
        if (ixcent(natoms)) then
          write (lw,11) natoms
          first = .true.
          lis = 1
          dow = .true.
          do while (lis.le.ntuples(natoms).and.nloc.lt.nab)
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
            sg = zero
            do i=1,natoms
              sg = sg+aom(inside(i),:,:)
            enddo
            call dsyev ('V','U',nab,sg,nab,eigen(1:nab),work,n3,info)
            call errdsyev (info,'uhfoqsloc.f')
            uvec=sg
c
c-----------Test that sg() is definite positive
c
            do m=1,nab
              if (eigen(m).gt.0d0) then
              elseif (eigen(m).gt.epsneg) then
                write (lerr,222) eigen(m),inside(1:natoms)
                eigen(m)=0d0
              else
                write (lerr,1124) m,eigen(m)
                stop
              endif
            enddo
            sg = 0d0   ! This is really S^{1/2}
            do k=1,nab
              do l=1,nab
                forall (m=1:nab) uveceig(m)=uvec(l,m)*sqrt(eigen(m))
                sg(k,l) = dot_product(uvec(k,:),uveceig(:))
              enddo
            enddo
c
c-----------This is S^{1/2} * 1-RDM * S^{1/2}. The previous "localized"
c           density is substracted
c
            v = matmul(sg,matmul(rho,transpose(sg))) - sloc0
            eigs = 0d0
            call dsyev ('V','U',nab,v,nab,eigs(1:nab),work,n3,info)
            call errdsyev (info,'uhfoqsloc.f')
            if (largwr) then
              write (lw,*) '#'
              write (lw,432) (inside(i),i=1,natoms)
              write (lw,*) '#'
              write (lw,133) (eigs(k),k=1,nab)
            endif
c
c-------------Be carefull with the following loop. It only works if the 
c             used diagonalization routine (in this case, 'dsyev') 
c             returns the eigenvalues in ascencing order. Otherwise
c             before this loop eigenvalues and eigvectors should be
c             reordered such that eigs(1) is the smallest eigvenvalue
c             eigs(2) the following, and so on.
c
            do i = nab,1,-1
              if (eigs(i).lt.critov(natoms)) exit
              nloc = nloc + 1
              if (natoms.eq.2) then
                nbab=nbab+1
                if (nbab.gt.nmo*(nmo-1)/2) then
                  stop '# uhfoqsloc.f.f:Increase first dimension of iab'
                endif
                iab(1:2,nbab)=inside(1:2)
                bond(nbab)=nloc
              endif
              nmoc(natoms) = nmoc(natoms)+1
              if (nloc.gt.nab) write (lw,433)
              if (nloc.gt.nab) stop
              eigloc(nloc)    = eigs(i)
              ctmp(:,nloc)    = v(:,i)
              do j=1,nab
                do k=1,nab
                  prod=ctmp(j,nloc)*ctmp(k,nloc)*eigloc(nloc)
                  sloc(j,k)=sloc(j,k)+prod
                enddo
              enddo
              if (dow) then
                write (lw,113) natoms,critov(natoms)
                dow = .false.
              endif
              if (nloc.le.nmotarget) then
                write (lw,116) nloc,eigs(i),
     &            (atnam(inside(j))(3:4),inside(j),j=1,natoms)
                if (natoms.eq.2) then
                  i1=min(inside(1),inside(2))
                  i2=max(inside(1),inside(2))
                  enlaces(i1,i2)=enlaces(i1,i2)+1
                endif
              else
                write (locform(nloc-nmotarget),116) nloc,eigs(i),
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
      if (nloc.lt.nmotarget) then
        write (lw,434)
        warn = .true.
        return
      elseif (nloc.eq.nmotarget) then
        write (lw,300)
        write (lw,301) (nmoc(i),i=1,mxcent)
      else
        write (lw,435)
        do i=nmotarget+1,nloc
          write (lw,'(a)') trim(locform(i-nmotarget))
        enddo
        return
      endif
c
c-----Print type of bond between each pair of atoms
c
      write (lw,320)
 320  format (' # Type of bond between pairs of atoms',/,' # ',70('-'))
      do i=1,ncent-1
        do j=i+1,ncent
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
          write (lw,321) i,j,typebond
        enddo
      enddo
 321  format (' # Type of bond betweeen atoms ',2I4,' : ',a30)
      c(1:nloc,:) = transpose(ctmp(:,1:nloc))
c
c-----Each OPENLOC localized MO (LMO) is normalized to 1.0 in the
c     fragment it was obtained, but is normalized to eigloc^{-1} in R^3. 
c     For this reason each LMO will be multiplied by sqrt(eigloc) such
c     that it will be normalized in R^3 instead of the fragment.
c
******forall (i=1:nloc) c(i,1:n)=c(i,1:n)*sqrt(eigloc(i))
c
c-----The above is true but the diagonalization routine already returns
c     normalized eigenvectors. That is why we comment the line 'forall'
c

      if (largwr) then
        write (lw,20)
        do i=1,nloc
          write (lw,33) i
          write (lw,30) (c(i,j),j=1,nab)
        enddo
      endif
c
c-----Compute AOM between the new nonorthonormal localized MOs
c
      aomold = aom
      do i=1,ncent
        aom(i,1:nloc,1:nloc) = matmul(c(1:nloc,:),
     &  matmul(aomold(i,:,:),transpose(c(1:nloc,:))))
      enddo
c
c-----Write diagonal AOM elements in each atom or fragment.
c
      write (lw,68) 'non-orthonormal'
      write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
      do ip=1,nloc
        write (lw,410,advance='no') ip
        if (ncent.gt.10) then
          write (lw,41 ) (hundred*aom(k,ip,ip),k=1,10)
          write (lw,411) (hundred*aom(k,ip,ip),k=11,ncent)
          write (lw,*) 
        else
          write (lw,41) (hundred*aom(k,ip,ip),k=1,ncent)
        endif
      enddo
c
c-----------------------------------------------------------------------
c
c-----effective number of atoms expanded by each localized MO
c
      do ip=1,nloc
        xloci(ip)=0d0
        do k = 1, ncent
          spp = aom(k,ip,ip)
          xloci(ip) = xloci(ip) + spp * spp
        enddo
      enddo
      write (lw,520) 'non-orthonormal',(1d0/xloci(ip),ip=1,nloc)
c
c-----Approximate two-center delocalization indices
c
      write (lw,*)
      write (lw,*) '# Approximate two-center delocalization indices'
      write (lw,*) '# ---------------------------------------------'
      do i=1,nbab
        i1=iab(1,i)
        i2=iab(2,i)
        ib=bond(i)
        di = 2d0 * aom(i1,ib,ib) * aom(i2,ib,ib)
        write (lw,65) iab(1:2,i),di
      enddo
 65   format (' # Atoms ',2I3,'   Delta = ',F15.8)
c
c-----Orthogonalize the nonorthonormal Localized MOs by Lowdin. This is
c     equivalent to the maximum overlap method. First compute the overlap
c     in R^3 of the nonorthonormal Localized MOs.
c
c-----First, obtain S, where S is the overlap matrix in R^3 betwen
c     the just found nonorthonormal Localized MOs.
c
      sumaom=0d0
      do k=1,ncent
        sumaom(:,:)=sumaom(:,:)+aom(k,1:nloc,1:nloc)
      enddo
      if (largwr) then
        write (lw,*) '#'
        write (lw,*) '# Overlaps in R^3 of OQS MOs'
        do i=1,nloc
          write (lw,133) (sumaom(i,j),j=1,i)
        enddo
      endif
      n3=3*nloc-1
      uvec(1:nloc,1:nloc) = sumaom
c
c-----Diagonalize S
c
      call dsyev ('V','U',nloc,uvec(1:nloc,1:nloc),
     &      nloc,eigen(1:nloc),work,n3,info)
      call errdsyev (info,'uhfoqsloc.f')
c
c-----Obtain  S^{-1/2} (stored in sumaom())
c
      sumaom = 0d0
      do i=1,nloc
        do j=1,nloc
          value=0d0
          do k=1,nloc
            value=value+uvec(i,k)*uvec(j,k)/sqrt(eigen(k))
          enddo
          sumaom(i,j) = value
        enddo
      enddo
      co(1:nab,1:nloc)=matmul(transpose(c(1:nloc,1:nab)),sumaom(:,:))
      co(1:nloc,1:nab)=transpose(co(1:nab,1:nloc))
      do i=1,ncent
        aomo(i,1:nloc,1:nloc)=matmul(transpose(sumaom),
     &  matmul(aom(i,1:nloc,1:nloc),sumaom))
      enddo
c
c-----Write orthonormalized MOs in terms of the canonical MOs
c
      if (largwr) then
        write (lw,2001)
        do i=1,nloc
          write (lw,33) i
          write (lw,30) (co(i,j),j=1,nab)
        enddo
      endif
c
c-----------------------------------------------------------------------
c
c-----Write diagonal AOM elements in each atom or fragment between
c     orthonormalized localized MOs
c
      write (lw,68) 'orthonormal'
      write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
      do ip=1,nloc
        write (lw,410,advance='no') ip
        if (ncent.gt.10) then
          write (lw,41)  (hundred*aomo(k,ip,ip),k=1,10)
          write (lw,411) (hundred*aomo(k,ip,ip),k=11,ncent)
          write (lw,*) 
        else
          write (lw,41) (hundred*aomo(k,ip,ip),k=1,ncent)
        endif
      enddo
c
c-----------------------------------------------------------------------
c
c-----effective number of atoms expanded by each orthonormal MO
c
      do ip=1,nloc
        xloci(ip)=0d0
        do k = 1, ncent
          spp = aomo(k,ip,ip)
          xloci(ip) = xloci(ip) + spp * spp
        enddo
      enddo
      write (lw,520) 'orthonormalized ',(1d0/xloci(ip),ip=1,nloc)
c
c-----------------------------------------------------------------------
c
c-----Overlaps in R^3 of the orthonormalized OQS MOs
c
      sumaom=0d0
      do k=1,ncent
        sumaom=sumaom+aomo(k,1:nloc,1:nloc)
      enddo
      if (largwr) then
        write (lw,*) '#'
        write (lw,*) '# Overlaps in R^3 of the orthonormalized OQS MOs'
        do i=1,nloc
          write (lw,133) (sumaom(i,j),j=1,i)
        enddo
      endif
c
c-----------------------------------------------------------------------
c
      deallocate (eigs     )
      deallocate (ntuples  )
      deallocate (work     )
      deallocate (v        )
      deallocate (ctmp     )
      deallocate (sloc     )
      deallocate (sloc0    )
      deallocate (uveceig  )
      deallocate (inside   )
      deallocate (eigloc   )
      deallocate (madis    )
      deallocate (coord    )
      deallocate (nmoc     )
      deallocate (wh       )
c
      write (lw,10) '  E N D  ','     '
      return
      call timer (4,iopenloc,'_uhfoqsloc   ',-1)
c
c.....Formats
c
 410  format (' # ',I4)
 20   format (/' # Rotation Matrix (READ BY ROWS)',6x,
     & 'Open System MO_i =  SUM_j CANMO_j C_ij',/,' # ',
     & 'Open System MOs have been renormalized such that <i|i>_R^3=1')
 2001 format (/' # Rotation Matrix (READ BY ROWS)',6x,
     & 'Orthonormal OQS MO_i =  SUM_j CANMO_j C_ij')
 30   format (5(1x,E15.8))
 41   format (1000(10x,10(1x,F5.1,'%')))
 411  format (1000(17x,10(1x,F5.1,'%'),/))
 300  format (' #',/,' # ',70('-'),/,
     &  ' # CONVERGED (1--N)-CENTER OPEN SYSTEM ANALYSIS',/,' #')
 66   format (/' # Diagonal overlaps of each MO in each atom')
 68   format (/' # Localization of each ',a,
     &  ' Open System MO in each atom')
 33   format (' # MO ',I3)
 113  format (/,1x,'# ',78('-'),/,2x,I2,
     & '-center Open System analysis, Critical overlap = ',F12.5,/,
     & ' # ',78('-'))
 116  format (4x,'MO ',I4,' with eigenvalue = ',F12.5,
     & ') localized on atom(s)',10(1x,'(',A2,I2,')'))
 21   format (' # CRITOV (1 ..',I2,')  = ',10(1x,F5.3))
 10   format (/,' # ',90('-'),/,' # ',a,
     & '(1--N)-CENTER OPEN SYSTEM LOCALIZATION ANALYSIS',3x,a,
     & /,' # ',90('-'))
 11   format (' #',/,
     & ' # Open System localization analysis on ',I2,' centers',/,' #')
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
 351  format (/1x,"# Writing file '",a,"'",14x,
     &  "with orthonormalized MOs")
 1000 format (
     & ' # Open System Localized MOs are orthonormalized with Lowdin')
 112  format (2(' #',/),' #',90('-'),/,' # ',a,
     & ' OF RESULTS WITH ORTHONORMALIZED OPTIMIZED MOs',/,' #',90('-'))
 51   format (' # In single atom analyses, H atoms are',a,'skipped')
 52   format (
     & " # In 'atom1-atom2' searches, maximum distance index is ",I2)
 432  format (
     & ' # Eigenvalues of the Open System of atoms ',20(1x,I3,'+'))
 433  format (' #',/,
     & ' # Too many Open Systems Localized MOs',/,
     & ' # Try again decreasing the critov() values',/,' #')
 434  format (' #',/,' # uhfoqsloc.f: Unable to localize all the MOs',
     & /,' #',/,' # NON CONVERGED (1--N)-CENTER OPEN SYSTEM ANALYSIS',
     & /,' #')
 435  format (' #',
     & /,' # uhfoqsloc.f: ! WARNING: Too many localized MOs !',/,' #')
 100  format (' # uhfoqsloc.f: File ',a,' NOT FOUND')
 1123 format (' # NMO in file ',a,'(',I4,') # Current value = ',I4)
 1124 format (' # Eigenvalue ',I4,' of sg() is negative, ',E17.10)
 28   format (' # The value of EPSNEG is ',E15.8)
 222  format (' # Negative eigenvalue ',E15.8,
     &    ' set to 0.0 diagonalizing fragment ',100I4)
 301  format (' # Number of MOs involving 1, 2, 3, ... centers',/,
     &        ' # ',20(1x,I3))
 352  format (" # Writing file '",a,"' with ",I2,
     &  "-center ortho-normalized MOs")
 353  format (" # Writing file '",a,"' with ",I2,
     &  "-center non-orthonormal MOs")
 520  format (' # Atoms expanded by each ',a,' OQS MO',/,
     &    1000(5(3x,F12.5),/))
 133  format (5(1x,F15.8))
      end
