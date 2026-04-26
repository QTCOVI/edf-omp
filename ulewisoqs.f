c
c-----------------------------------------------------------------------
c
      subroutine ulewisoqs (aom,mal,mbe,nfrag,atoms_in_frag,
     & ifrag,mxfrag,covx,epsneg,maxb,nloc,udat,wfnfile,maxcritov,
     & critov,ixfrag,verbose,warn,skiph,lw,lerr)
c
c-----PARAMETERS--------------------------------------------------------
c
c-----aom()                            INPUT/OUTPUT
c
c     Atomic Overlap Matrix (AOM) of MOs.
c
c-----mal,mbe
c
c     Number of ALPHA and BETA electrons of the molecule
c
c-----nfrag
c
c     Number of fragments
c
c-----atoms_in_frag(1..nfrag)
c
c     Number of atoms in each fragment
c
c-----ifrag(1..ncent,1..nfrag)
c
c     Labels of the atoms that belong to each fragment. If there are F
c     atoms in fragment i (i.e. atoms_in_frag(i)=F), only the first
c     'F' indices ifrag(1..F,i) have meaning.
c
c-----mxfrag                         INPUT
c
c     Maximum number of N in the search of (nc,2e) bonds.
c
c-----covx                           INPUT
c
c     Factor multiplying the sum of covalent radii of two atoms A and B  
c     used in the connect routine to determined whether A and B are 
c     linked or not from a purely geomtrical point of view. It is only
c     used in the 'confrag.f' routine.
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
c-----maxb                          INPUT
c
c     Maximum number of links between two fragments A and B for which
c     a (2c,2e) localized MOs between these fragments is searched. 
c     Fragments directly connected in a geometrical sence have a 
c     'geometrical distance' or link of 1, if A is connected to C which, 
c     in turn, is connected to B the 'geometrical distance' or link 
c     between A and B is 2, and so on. 
c
c-----nloc                           OUTPUT
c
c     Number of localized MOs found. It should be equal to nloctarget
c     (see below) if the routine ends nicely
c
c-----udat                           INPUT
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
c-----critov(1)...critov(mxfrag)     INPUT
c
c     Tolerance used in localization: If an eigvenvalue of
c     S^(1/2) rho S^{1/2} is greater that critov(i) it means that the
c     associated eigenvector corresponds to an MO assumed to be loca-
c     lized in the 'i' fragments involved in the diagonalization.
c
c-----ixfrag(1)...ixfrag(mxfrag)     INPUT
c
c     Logical variable: If ixfrag(i)=.true. (i-fragment,2e) bonds are
c     searched, otherwise (i-fragment,2e) bonds are skipped.
c
c
c-----verbose                        INPUT
c
c     Logical variable: if verbose=.true. a large output is requested,
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
c     fragment localizations. Otherwise hydrogens are not skipped.
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
      USE      space_for_sym
      include 'implicit.inc'
      include 'constants.inc'
      include 'param.inc'
      include 'wfn.inc'
      include 'corr.inc'
      include 'error.inc'
      include 'mline.inc'
      integer(kind=4) atoms_in_frag(nfrag)
      integer(kind=4) ifrag(ncent,nfrag)
      real   (kind=8),allocatable, dimension (:,:) :: c
      real   (kind=8),allocatable, dimension (:,:) :: co
      real   (kind=8) aom(ncent,nmo,nmo)
      real   (kind=8) eos(nfrag),eosord(nfrag)
      integer(kind=4) imost(nmo)
      integer(kind=4) ieos(nfrag)
      integer(kind=4) ntofrag(nfrag)
      real   (kind=8), allocatable, dimension (:,:,:) :: fom
      real   (kind=8), allocatable, dimension (:,:,:) :: fomold
      real   (kind=8), allocatable, dimension (:,:,:) :: fomo
      real   (kind=8), allocatable, dimension (:,:)   :: sumfom
      real   (kind=8) fominp(nfrag,nmo,nmo)

      real   (kind=8), allocatable, dimension (:,:) :: sg
      real   (kind=8), allocatable, dimension (:)   :: eigen
      real   (kind=8), allocatable, dimension (:,:) :: uvec
      real   (kind=8), allocatable, dimension (:,:) :: uvecx

      integer(kind=4),allocatable, dimension (:)     :: nfrags
      integer(kind=4),allocatable, dimension (:)     :: infrag
      integer(kind=4),allocatable, dimension (:)     :: iordi,nidxi
      integer(kind=4),allocatable, dimension (:)     :: iordj,nidxj
      integer(kind=4),allocatable, dimension (:)     :: ndiffe
      integer(kind=4),allocatable, dimension (:,:)   :: idiffe
      logical eqfrags, conflicting

      integer(kind=4), allocatable, dimension (:,:) :: iab
      real   (kind=8), allocatable, dimension (:)   :: xloci
      integer(kind=4) bond(nmo)
      integer(kind=4) bonds(nfrag,nfrag)
      integer(kind=4),allocatable, dimension (:)     :: isigma
      real   (kind=8),allocatable, dimension (:,:)   :: onerdm
      real   (kind=8),allocatable, dimension (:)     :: eigs,work
      real   (kind=8),allocatable, dimension (:)     :: eigdeg
      real   (kind=8),allocatable, dimension (:)     :: eigloc
      real   (kind=8),allocatable, dimension (:)     :: uveceig
      real   (kind=8),allocatable, dimension (:)     :: ocup
      real   (kind=8),allocatable, dimension (:,:)   :: v
      real   (kind=8),allocatable, dimension (:,:)   :: sloc,sloc0

      real   (kind=8),allocatable, dimension (:,:)   :: xsumfom

      real   (kind=8),allocatable, dimension (:,:)   :: ctmp
      real   (kind=8),allocatable, dimension (:,:,:) :: fomx
      integer(kind=8),allocatable, dimension (:)     :: ntuples
      integer(kind=4),allocatable, dimension (:)     :: inside
      integer(kind=4),allocatable, dimension (:,:)   :: wh
      integer(kind=4),allocatable, dimension (:)     :: coord
      integer(kind=4),allocatable, dimension (:,:)   :: madis
      integer(kind=4),allocatable, dimension (:)     :: nmoc
      real   (kind=8),allocatable, dimension (:,:)   :: ccprev,ccprevo
      real   (kind=8),allocatable, dimension (:,:)   :: cret,coret
      real   (kind=8),allocatable, dimension (:,:)   :: dcoef
      logical warn,warno,is_hydrogen,only_one_atom,skiph,degsys
      parameter (maxcoord = 20)
      real   (kind=8), parameter :: fifty_percent = 50.0d+00
c
c-----Data used in the orthonormalization process
c
      real(kind=8),    allocatable,dimension (:,:)   :: cc
      real(kind=8),    allocatable,dimension (:,:)   :: cco
      real(kind=8),    allocatable,dimension (:,:)   :: ornor
c
      real   (kind=8)       epsneg
      real   (kind=8)       critov(maxcritov)
      logical               ixfrag(maxcritov)
      integer (kind=4)      p,q
      integer (kind=4)      FragsLewis
      logical               verbose,dow,first,valid
      integer(kind=4)       udat,udatnw
      character(len=*)      wfnfile
      character (len=mline) wfnloc
      character (len=mline) uppcase
      character (len=mline) locform(nmo+nmo)
      character (len=30)    typebond

      character (len=6) spin
      character (len=1) digs(0:9)
      data digs / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
c
c-----------------------------------------------------------------------
c
      call timer (2,iulewisoqs,'_ulewisoqs',-1)
c
c-----Obtain the input Fragment Overlap Matrices (FOM)
c
      fominp = 0.0d+00
      eos(1:nfrag) = 0.0d+00
      do ngr=1,nfrag
        do i=1,atoms_in_frag(ngr)
          fominp(ngr,:,:) = fominp(ngr,:,:) + aom(ifrag(i,ngr),:,:)
          eos(ngr) = eos(ngr) + charge(ifrag(i,ngr))
        enddo
      enddo


cAQUI
C     allocate (nfrags(ncent))
C     allocate (infrag(ncent))
C     nfrags(1:ncent) = 0
C     do i=1,nfrag
C       do k=1,ncent
C         if (atoms_in_frag(i).eq.k) then
C           nfrags(k) = nfrags(k) + 1
C           infrag(nfrags(k)) = i
C         endif
C       enddo
C     enddo
C     write (6,*) nfrags(1:ncent)
C     STOP



      

c
c-----Test if two fragments with the same number of atoms are comparable
c
C     allocate (iordi(ncent))
C     allocate (nidxi(ncent))
C     allocate (ndiffe(nfrag))
C     allocate (idiffe(nfrag,nfrag))
c
c-----Order the number of atoms of each fragment
c
C     do i=1,nfrag
C       iordi(i) = i
C       nidxi(i) = atoms_in_frag(i)
C     enddo
C     call iqcksort(nidxi,iordi,nfrag,1,nfrag)
c
c-----Determine the number of fragments with a given size.
c
C     nsizes = 1
C     ndiffe(nsizes) = 1
C     idiffe(nsizes,ndiffe(nsizes)) = 1
C     do i=2,nfrag
C       if (nidxi(iordi(i)).eq.nidxi(iordi(i-1))) then
C         ndiffe(nsizes) = ndiffe(nsizes) + 1
C       else
C         nsizes = nsizes + 1
C         ndiffe(nsizes) = 1
C       endif
C     enddo
C     write (lw,*) '#'
C     write (lw,9823) nsizes
C     do i=1,nsizes
C       write (lw,*) '# Size with ',nidxi(1),' atoms = ',ndiffe(i)
C     enddo
C     deallocate (iordi)
C     deallocate (nidxi)
C  
C     STOP
c
c-----Write composition of the fragments
c
      do i=1,nfrag
        write (lw,912) i,atoms_in_frag(i),ifrag(1:atoms_in_frag(i),i)
      enddo
c
c-----Obtain coordinations, connectivities, and distance arrays
c
      allocate (madis(nfrag,nfrag) )
      allocate (coord(nfrag)       )
      allocate (wh(nfrag,maxcoord) )
      call confrag (lw,covx,ncent,nfrag,atoms_in_frag,ifrag,maxcoord,
     &  verbose,nbonds,coord,wh,madis)
      deallocate (wh)
c
c-----Largest number of bonds between two 'connected' atoms
c
      longest_chain = maxval(madis)
c
c-----Analyze symmetry
c
      call allocate_space_for_sym (nmo,ncent)
      call sym (lw,.false.)

      ispin = 1
 1    continue

      if (ispin.eq.1) then
        spin = '-alpha'
        nab = nalpha
        nloctarget = mal
        allocate (onerdm(nab,nab))
        allocate (isigma(nab))
        isigma(1:nab) = ialpha(1:nab)
        onerdm(1:nab,1:nab) = c1ea(isigma(1:nab),isigma(1:nab))
      else
        spin = '-beta'
        nab = nbeta
        nloctarget = mbe
        allocate (onerdm(nab,nab))
        allocate (isigma(nab))
        isigma(1:nab) = ibeta(1:nab)
        onerdm(1:nab,1:nab) = c1ea(isigma(1:nab),isigma(1:nab))
      endif
c
c.....Allocate arrays
c
      allocate (eigdeg(nab)        )
      allocate (eigs(nab)          )
      allocate (eigloc(nab)        )
      allocate (ntuples(mxfrag)    )
      allocate (nmoc(mxfrag)       )
      allocate (inside(nfrag)      )
      allocate (work(3*nab-1)      )
      allocate (v(nab,nab)         )
      allocate (ctmp(nab,nab)      )
      allocate (uveceig(nab)       )
      nfrageff = min(mxfrag,nfrag)
      write (lw,10) ' B E G I N  ', trim(spin(2:leng(uppcase(spin))))
      write (lw,21) nfrageff,(critov(j),j=1,nfrageff)
      if (     skiph) write (lw,51) '     '
      if (.not.skiph) write (lw,51) ' NOT '
      write (lw,52) maxb
c
c-----Print the input AOM arrays
c
      allocate (fom(nfrag,nab,nab))
      allocate (eigen(nab))
      allocate (uvec(nab,nab))
      allocate (sg(nab,nab))
      allocate (iab(2,nab*(nab-1)/2))

      forall (i=1:nfrag) 
        fom(i,1:nab,1:nab) = fominp(i,isigma(1:nab),isigma(1:nab))
      endforall

      if (verbose) then
        write (lw,*) '# Fragment Overlap matrix (FOM)'
        write (lw,*) '#'
        do i=1,nfrag
          write (lw,*) ' Fragment ',i
          do j=1,nab
            write (lw,133) (fom(i,j,k),k=1,j)
          enddo
        enddo
        write (lw,*) '#'
        write (lw,*) '#    Total'
        do j=1,nab
          write (lw,133) (sum(fom(:,j,k)),k=1,j)
        enddo
      endif
c
c.....Start the process.
c
      FragsLewis = 1
      do while (FragsLewis.le.nfrageff) 
        xtuples = 1.0d+00
        do i = 0,FragsLewis-1
          xtuples = xtuples*dble(nfrag-i)/dble(FragsLewis-i)
        enddo
        ntuples(FragsLewis) = nint(xtuples)
        FragsLewis = FragsLewis+1
      enddo     
      allocate (sloc(nab,nab)         )
      allocate (sloc0(nab,nab)        )
      sloc    = 0.0d+00
      sloc0   = 0.0d+00
      ctmp    = 0.0d+00
      n3      = 3*nab-1
      nloc    = 0
      nmoc    = 0
      warn    = .false.
      degsys  = .false.
      nbab    = 0
      bonds   = 0
      FragsLewis  = 1
      write (lw,28) epsneg
      do while (FragsLewis.le.nfrageff.and.nloc.lt.nloctarget)
        if (ixfrag(FragsLewis)) then
          first = .true.
          lis = 1
          dow = .true.
          ntuf = ntuples(FragsLewis)
          do while (lis.le.ntuf.and.nloc.lt.nloctarget)
            call ngen (FragsLewis,nfrag,nfrag,inside,first)
            if (FragsLewis.eq.1) then
              is_hydrogen = nint(charge(inside(1))).eq.1
              only_one_atom = atoms_in_frag(inside(1)).eq.1
              is_hydrogen = is_hydrogen.and.only_one_atom
              if (is_hydrogen.and.skiph) goto 111
            elseif (FragsLewis.eq.2) then
              iat1=inside(1)
              iat2=inside(2)
              m12=madis(iat1,iat2)
              if (m12.lt.1.and.m12.gt.maxb) goto 111
            elseif (FragsLewis.ge.3) then
              valid = .true.
              do i=1,FragsLewis
                if (coord(inside(i)).lt.2) valid = .false.
              enddo
              if (.not.valid) goto 111
            endif
c
c-----------Group overlap integrals are obtained from the AOM
c
            sg = 0.0d+00
            do i=1,FragsLewis
              sg = sg+fom(inside(i),:,:)
            enddo
            call dsyev ('V','U',nab,sg,nab,eigen,work,n3,info)
            call errdsyev (info,'ulewisoqs.f')
            uvec = sg
c
c-----------Test that sg() is definite positive
c
            do m=1,nab
              if (eigen(m).gt.0.0d+00) then
              elseif (eigen(m).gt.epsneg) then
                write (lerr,222) eigen(m),epsneg,inside(1:FragsLewis)
                eigen(m) = 0.0d+00
              else
                write (lerr,1124) m,eigen(m)
              endif
            enddo
c
c-----------Compute S^{1/2} (stored in sg())
c
            sg = 0.0d+00 
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
            v = matmul(sg,matmul(onerdm,transpose(sg))) - sloc0
            eigs = 0.0d+00
            call dsyev ('V','U',nab,v,nab,eigs,work,n3,info)
            call errdsyev (info,'ulewisoqs.f')
            if (verbose) then
              write (lw,*) '#'
              write (lw,432) (inside(i),i=1,FragsLewis)
              write (lw,*) '#'
              write (lw,133) (eigs(k),k=1,nab)
            endif
c
c-----------Be carefull with the following loop. It only works if the 
c           used diagonalization routine (in this case, 'dsyev') 
c           returns the eigenvalues in ascencing order. Otherwise,
c           before this loop, eigenvalues and eigenvectors should be
c           reordered such that eigs(1) is the smallest eigvenvalue
c           eigs(2) the following, and so on.
c
            do i = nab,1,-1
              if (eigs(i).lt.critov(FragsLewis)) exit
              nloc = nloc + 1
              if (FragsLewis.eq.2) then
                nbab=nbab+1
                if (nbab.gt.nab*(nab-1)/2)
     &            stop '# ulewisoqs.f: Increase first dimension of iab'
                iab(1:2,nbab)=inside(1:2)
                bond(nbab)=nloc
              endif
              nmoc(FragsLewis) = nmoc(FragsLewis)+1
              if (nloc.gt.nloctarget) then
                write (lw,436) FragsLewis
                degsys  = .true.
                write (lw,671) (eigloc(k),k=1,nloctarget)
                write (lw,672) eigs(i)
                if (nloc-nloctarget.le.nmo+nmo) then
                  write (locform(nloc-nloctarget),1160) nloc,eigs(i),
     &                (inside(j),j=1,FragsLewis)
                  write (lw,'(a)') trim(locform(nloc-nloctarget))
                endif
                cycle
              endif
              eigloc(nloc) = eigs(i)
              ctmp(:,nloc) = v(:,i)
              do j=1,nab
                do k=1,nab
                  prod = ctmp(j,nloc)*ctmp(k,nloc)*eigloc(nloc)
                  sloc(j,k) = sloc(j,k)+prod
                enddo
              enddo
              if (dow) then
                write (lw,113) FragsLewis,critov(FragsLewis)
                dow = .false.
              endif
              if (nloc.le.nloctarget) then
                write (lw,1160) nloc,eigs(i),(inside(j),j=1,FragsLewis)
                if (FragsLewis.eq.2) then
                  i1=min(inside(1),inside(2))
                  i2=max(inside(1),inside(2))
                  bonds(i1,i2)=bonds(i1,i2)+1
                endif
              endif
            enddo
 111        continue
            lis=lis+1
          enddo
        endif
        sloc0 = sloc
        FragsLewis = FragsLewis+1 
      enddo
      deallocate (work)
c
 3344 continue
      if (nloc.lt.nloctarget) then
        write (lw,434)
        warn = .true.
        return
      elseif (nloc.eq.nloctarget) then
        write (lw,300)
        write (lw,301) (nmoc(i),i=1,mxfrag)
      else
        warn = .true.
        nloc = nloctarget 
      endif

      allocate (c(nloc,nab))
      c(1:nloc,:) = transpose(ctmp(:,1:nloc))
c
c-----The above is true but the diagonalization routine already returns
c     normalized eigenvectors. This is why we comment the line 'forall'
c
      if (verbose) then
        write (lw,20)
        do i=1,nloc
          write (lw,33) i
          write (lw,30) (c(i,j),j=1,nab)
        enddo
      endif
c
c-----------------------------------------------------------------------
c                                                                      !
c-----Compute AOM between the new nonorthonormal localized MOs         !
c                                                                      !
      allocate (fomold(nfrag,nab,nab))
      fomold = fom                                                     !
      do i=1,nfrag                                                     !
        fom(i,1:nloc,1:nloc) = matmul(c(1:nloc,:),                     !
     &  matmul(fomold(i,:,:),transpose(c(1:nloc,:))))                  !
      enddo                                                            !
c                                                                      !
c-----------------------------------------------------------------------
c
c-----Write diagonal AOM elements in each atom or fragment. 
c
      write (lw,68) 'non-orthonormal'
      write (lw,690) (i,i=1,min(nfrag,15))
      do ip=1,nloc
        write (lw,410,advance='no') ip
        if (nfrag.gt.15) then
          write (lw,41 ) (100.0d+00*fom(k,ip,ip),k=1,15)
          write (lw,411) (100.0d+00*fom(k,ip,ip),k=16,nfrag)
        else
          write (lw,41) (100.0d+00*fom(k,ip,ip),k=1,nfrag)
        endif
        previousMO = -1.0d+00
        do k=1,nfrag
          if (fom(k,ip,ip).gt.previousMO) then
            previousMO = fom(k,ip,ip)
            imost(ip) = k
          endif
          ieos(k) = k
          eosord(k) = 100.0d+00*fom(k,ip,ip)
        enddo
c
c-------Order Localization by increasing order.
c
        call qcksort (eosord,ieos,1,nfrag)
c
c-------Determine NEQUAL, the number of Maximal Localizations that differ
c       by less than 0.1.
c
        nequal = 1
        do k=nfrag-1,1,-1
          k1 = ieos(k)
          k2 = ieos(k+1)
          if (abs(eosord(k1)-eosord(k2)).lt.1.0D-01) then
            nequal = nequal + 1
          else
c
c-----------There are no other fragments with a Maximal
c           Localization. Hence, skip the remaining values.
c
            exit
          endif
        enddo
c
c-------Add 1/nequal to the oxidation state of each fragment with a 
c       maximal delocalization.
c
        do k=1,nequal
          kis = ieos(nfrag-k+1)
          eos(kis) = eos(kis) - 1.0d+00/dble(nequal)
        enddo
      enddo
      ntofrag(1:nfrag) = 0
      do ip=1,nloc
        ntofrag(imost(ip)) = ntofrag(imost(ip)) + 1
*       write (lw,171) ip,imost(ip)
      enddo
cAQUI
c
c-----Test if two fragments with the same number of atoms are comparable
c
      allocate (iordi(ncent))
      allocate (iordj(ncent))
      allocate (nidxi(ncent))
      allocate (nidxj(ncent))
      do i=1,nfrag-1
        ini = atoms_in_frag(i)
        do j=i+1,nfrag
          inj = atoms_in_frag(j)
          if (ini.eq.inj) then
            zi = 0.0d+00
            zj = 0.0d+00
            do k=1,ini
              zi = zi + charge(ifrag(k,i))
              zj = zj + charge(ifrag(k,j))
              iordi(k) = k
              iordj(k) = k
              nidxi(k) = idx(ifrag(k,i),1)
              nidxj(k) = idx(ifrag(k,j),1)
            enddo
            call iqcksort(nidxi,iordi,ini,1,ini)
            call iqcksort(nidxj,iordj,ini,1,inj)
            eqfrags = .true.
            do k=1,ini
              if (nidxi(iordi(k)).ne.nidxj(iordj(k))) then
                eqfrags = .false.
              endif
            enddo
            conflicting = eqfrags .and. (abs(zi-zj).lt.1.0d-04)
            if (conflicting.and.ntofrag(i).ne.ntofrag(j)) then
              ntotot = ntofrag(i)+ntofrag(j)
              ntofrag(i) = ntotot/2
              ntofrag(j) = ntotot/2
            endif
          endif
        enddo
      enddo
      write (lw,*)
      do k=1,nfrag
        write (lw,172) trim(spin(2:)),k,ntofrag(k)
      enddo
      deallocate (iordi)
      deallocate (iordj)
      deallocate (nidxi)
      deallocate (nidxj)
cAQUI
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----effective number of atoms expanded by each localized MO          !
c                                                                      !
      allocate (xloci(nloc))
      forall (ip=1:nloc) 
        xloci(ip)=dot_product(fom(:,ip,ip),fom(:,ip,ip))
      endforall
      write (lw,520) 'non-orthonormal',(1.0d+00/xloci(ip),ip=1,nloc)   !
c
c-----Print type of bond between each pair of atoms
c
      write (lw,320)
      kini = 0
      do i=1,nfrag-1
        do j=i+1,nfrag
          if (bonds(i,j).gt.0) then
            deltais=0.0d+00
            do k=1,bonds(i,j)
              iabmin = min(iab(1,k+kini),iab(2,k+kini))
              iabmax = max(iab(1,k+kini),iab(2,k+kini))
              ib = bond(k+kini)
              di = 2.0d+00*fom(iabmin,ib,ib)*fom(iabmax,ib,ib)
              deltais = deltais + di
              if (k.eq.1) write (lw,321) i,j,bonds(i,j),ib,di
              if (k.ne.1) write (lw,322) ib,di
            enddo
            if (abs(deltais).gt.1.0d-06) write (lw,323) deltais
          endif
          kini=kini+bonds(i,j)
        enddo
      enddo
c
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
      allocate (sumfom(nloc,nloc))
      sumfom=0.0d+00                                                   !
      do i=1,nloc
        do j=1,i
          do k=1,nfrag                                                 !
            sumfom(i,j)=sumfom(i,j)+fom(k,i,j)                         !
          enddo                                                        !
          sumfom(j,i)=sumfom(i,j)
        enddo
      enddo
      if (verbose) then                                                !
        write (lw,*) '#'                                               !
        write (lw,*) '# Overlaps in R^3 of the nonorthormal OQS MOs'   !
        do i=1,nloc                                                    !
          write (lw,133) (sumfom(i,j),j=1,i)                        !
        enddo                                                          !
      endif                                                            !
c                                                                      !
      n3=3*nloc-1                                                      !
      allocate (work(n3))
      allocate (uvecx(nloc,nloc))
      uvecx(1:nloc,1:nloc) = sumfom(1:nloc,1:nloc)
c                                                                      !
c-----Diagonalize S                                                    !
c                                                                      !
      call dsyev ('V','U',nloc,uvecx,nloc,eigen(1:nloc),work,n3,info)  !
      call errdsyev (info,'ulewisoqs.f')                                !
c                                                                      !
c-----Obtain  S^{-1/2} and store it in sumfom()                        !
c                                                                      !
      sumfom = 0.0d+00                                                 !
      do i=1,nloc                                                      !
        do j=1,nloc                                                    !
          value=0.0d+00                                                !
          do k=1,nloc                                                  !
            xeig = max(1.0d-20,eigen(k))
            value=value+uvecx(i,k)*uvecx(j,k)/sqrt(xeig)               !
          enddo                                                        !
          sumfom(i,j) = value                                          !
        enddo                                                          !
      enddo                                                            !

      allocate (co(1:nloc,1:nloc))
      co(1:nloc,1:nloc) = matmul(transpose(sumfom),c(1:nloc,1:nloc))

      allocate (fomo(1:nfrag,nloc,nloc))
      do i=1,nfrag                                                     !
        fomo(i,1:nloc,1:nloc)=matmul(transpose(sumfom),                !
     &  matmul(fom(i,1:nloc,1:nloc),sumfom))                           !
      enddo                                                            !
c
c-----Write orthonormalized MOs in terms of the canonical MOs
c
      if (verbose) then
        write (lw,2001)
        do i=1,nloc
          write (lw,33) i
          write (lw,30) (co(i,j),j=1,nloc)
        enddo
      endif
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----Write diagonal AOM elements in each atom or fragment between     !
c     orthonormalized localized MOs                                    !
c                                                                      !
c                                                                      !
      write (lw,68) 'orthonormal'                                      !
      write (lw,690) (i,i=1,min(nfrag,15))
      do ip=1,nloc                                                     !
        write (lw,410,advance='no') ip                                 !
        if (nfrag.gt.15) then                                          !
          write (lw,41)  (100.0d+00*fomo(k,ip,ip),k=1,15)              !
          write (lw,411) (100.0d+00*fomo(k,ip,ip),k=16,nfrag)          !
        else                                                           !
          write (lw,41) (100.0d+00*fomo(k,ip,ip),k=1,nfrag)            !
        endif                                                          !
      enddo                                                            !
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----effective number of atoms expanded by each orthonormal MO        !
c                                                                      !
      forall (ip=1:nloc)                                               !
        xloci(ip)=dot_product(fomo(:,ip,ip),fomo(:,ip,ip))             !
      endforall                                                        !
      write (lw,520) 'orthonormalized',(1.0d+00/xloci(ip),ip=1,nloc)   !
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----Overlaps in R^3 of the orthonormalized OQS MOs                   !
c                                                                      !
      sumfom(1:nloc,1:nloc) = sum(fomo(1:nfrag,1:nloc,1:nloc))         !
      do i=1,nloc
        do j=1,i
          sumfom(i,j) = sum(fomo(1:nfrag,i,j))
          sumfom(i,i) = sumfom(i,j)
        enddo
      enddo
      if (verbose) then                                                !
        write (lw,*) '#'                                               !
        write (lw,*) '# Overlaps in R^3 of the orthonormalized OQS MOs'!
        do i=1,nloc                                                    !
          write (lw,133) (sumfom(i,j),j=1,i)                           !
        enddo                                                          !
      endif                                                            !
c                                                                      !
c-----------------------------------------------------------------------
c                                                                      !
c-----Find localized MOs in terms of the primitive basis               !
c                                                                      !
      allocate (cc(nloc,nprims))                                       !
      allocate (cco(nloc,nprims))                                      !
      allocate (cret(nloc,nloc))
      allocate (coret(nloc,nloc))
      allocate (ccprev(nloc,nprims))
      allocate (ccprevo(nloc,nprims))
      allocate (dcoef(nloc,nprims))

      cret=c(1:nloc,1:nloc)
      coret=co(1:nloc,1:nloc)
      do i=1,nloc
        dcoef(i,1:nprims) = coef(nmo+isigma(i),1:nprims)
      enddo
      ccprev  = matmul(cret, dcoef)
      ccprevo = matmul(coret,dcoef)
      do i=1,nloc
        cc (i,:)  = ccprev (i,:)
        cco(i,:)  = ccprevo(i,:)
      enddo
      deallocate (ccprev )
      deallocate (ccprevo)
      deallocate (dcoef  )
c                                                                      !
c-----------------------------------------------------------------------
c
c-----Write AOM between the non-orthonormal localized MOs
c
      if (verbose) then
        write (lw,115) 'Non-orthonormal'
        do k=1,nfrag
          write (lw,*) '# Fragment',k
          do i=1,nloc
            write (lw,101) (fom(k,i,j),j=1,i)
          enddo
        enddo
      endif

      if (verbose) then
        write (lw,115) 'Orthonormal'
        do k=1,nfrag
          write (lw,*) '# Center',k
          do i=1,nloc
            write (lw,101) (fomo(k,i,j),j=1,i)
          enddo
        enddo
      endif
c
c-----Compute the overlaps in R^3 between non-orthonormal and
c     orthonormal localized MOs
c
      allocate (ornor(nloc,nloc))
      ornor=matmul(c(1:nloc,1:nloc),transpose(co(1:nloc,1:nloc)))
      write (lw,200) 'non-orthonormal and orthonormal MOs'
      do i=1,nloc
        write (lw,101)(ornor(i,j),j=1,nloc)
      enddo
c
c-----write a WFN file with non-orthonormal localized MOs
c
      inpr = nprims
      wfnloc = trim(wfnfile)//"-open"//trim(spin)
      udatnw = udat+10
      open (unit=udatnw,file=trim(wfnloc),status='unknown')
      write (lw,*)
      write (lw,351) trim(wfnloc)
      allocate (ocup(nloc),stat=ier)
      ocup(1:nloc)=2d0
      if (ier.ne.0) stop '# ulewisoqs.f: Cannot allocate ocup()'
      call cdafhmos (udat,udatnw,0,ocup,cc,nloc,inpr)
c
c-----write a WFN file with orthonormal localized MOs
c
      inpr = nprims
      wfnloc = trim(wfnfile)//"-open"//trim(spin)//"-ortho"
      udatnw = udat+10
      open (unit=udatnw,file=trim(wfnloc),status='unknown')
      write (lw,352) trim(wfnloc)
      call cdafhmos (udat,udatnw,0,ocup,cco,nloc,inpr)
      deallocate (ocup,stat=ier)
      if (ier.ne.0) stop '# ulewisoqs.f: Cannot deallocate ocup()'
c
c-----Write file with AOM between the non-orthonormal localized MOs
c     and between the orthonormal localized MOs
c
      lu18 = udat+11
      wfnloc = trim(wfnfile)//"-open"//trim(spin)//".aom"
      open (unit=lu18,file=trim(wfnloc),status='unknown')
      write (lw,353) trim(wfnloc)
      do k=1,nfrag
        write (lu18,'(I4,a)') k,' <=== AOM within this fragment'
        write (lu18,80) ((fom(k,m,j),m=1,j),j=1,nloc)
      enddo
      close (unit=lu18)
      wfnloc = trim(wfnfile)//"-open"//trim(spin)//"-ortho.aom"
      open (unit=lu18,file=trim(wfnloc),status='unknown')
      write (lw,354) trim(wfnloc)
      do k=1,nfrag
        write (lu18,'(I4,a)') k,' <--- AOM in this fragment'
        write (lu18,80) ((fomo(k,m,j),m=1,j),j=1,nloc)
      enddo
      close (unit=lu18)
c
c-----Analyze the MOs and their localizations in all the fragments
c
      critics=0.02d0
      allocate (fomx(nfrag,nloc,nloc))
      fomx=fom(1:nfrag,1:nloc,1:nloc)
      call aomtest (fomx,nloc,nfrag,critics,lw,verbose)
      write (lw,3340) 'Open quantum-systems localized MOs'
      call newdrvmult (nloc,nprims,ncent,cc,lw)

      write (lw,10) '  E N D  ',trim(spin(2:leng(uppcase(spin))))

      deallocate (onerdm   )
      deallocate (fom      )
      deallocate (eigen    )
      deallocate (uvec     )
      deallocate (uvecx    )
      deallocate (co       )
      deallocate (cret     )
      deallocate (coret    )
      deallocate (sg       )
      deallocate (c        )
      deallocate (fomold   )
      deallocate (fomo     )
      deallocate (sumfom   )
      deallocate (xloci    )
      deallocate (isigma   )
      deallocate (eigs     )
      deallocate (eigdeg   )
      deallocate (iab      )
      deallocate (ntuples  )
      deallocate (work     )
      deallocate (v        )
      deallocate (uveceig  )
      deallocate (sloc     )
      deallocate (sloc0    )
      deallocate (ctmp     )
      deallocate (inside   )
      deallocate (eigloc   )
      deallocate (nmoc     )
      deallocate (cc       )
      deallocate (cco      )
      deallocate (ornor    )
      deallocate (fomx     )


      if (ispin.eq.2) then
        write (lw,985)
        sumeos = 0.0d+00
        do i=1,nfrag
          if (abs(eos(i)-anint(eos(i))).lt.1.0D-1) then
            eosprint = dble(anint(eos(i)))
            sumeos = sumeos + eosprint
          else
            eosprint = eos(i)
            sumeos = sumeos + eosprint
          endif
          write (lw,986) i,eosprint
        enddo
        write (lw,987) sumeos
        call timer (4,iulewisoqs,'_ulewisoqs',-1)
        return
      endif
      ispin = 2
      goto 1
      deallocate (madis    )
      deallocate (coord    )
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
 41   format (1000(10x,15(1x,F5.1,'%')))
 411  format (1000(17x,15(1x,F5.1,'%'),/))
 300  format (' #',/,' # ',70('-'),/,
     &  ' # CONVERGED (1--N)-FRAGMENT OPEN SYSTEM ANALYSIS',/,' #')
 66   format (/' # Diagonal overlaps of each MO in each fragment')
 68   format (/' # Localization of each ',a,' OQS MO in each fragment')
 33   format (' # MO ',I3)
 113  format (/,1x,'# ',78('-'),/,2x,I2,
     & '-fragment OQS analysis, Critical overlap = ',F12.5,/,
     & ' # ',78('-'))
 116  format (4x,'MO ',I4,' with eigenvalue = ',F12.5,
     & ') localized on atom(s)',10(1x,'(',A2,I2,')'))
 1160 format (4x,'MO ',I4,' with eigenvalue = ',F15.8,
     & ') localized on fragment(s)',10(1x,I2))
 21   format (' # CRITOV (1 ..',I2,')  = ',10(1x,F5.3))
 10   format (/,' # ',80('-'),/,' # ',a,
     & '(1--N)-CENTER OPEN SYSTEM LOCALIZATION ANALYSIS','  SPIN ',a,/,
     & ' # ',80('-'))
 11   format (' #',/,
     & ' # OQS localization analysis on ',I2,' fragments',/,' #')
 69   format (' #',/,' # SUM of eigenvalues = ',F18.8,/,' #')
 70   format (' # ',7x,I4,'-fragment = ',F18.8)
 71   format (' # ',/,
     &  ' # Electrons associated to successful diagonalizations')
 710  format (' # Fragment ',I4,'   Electrons = ',F18.8)
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
     & ' # Eigenvalues of the OQS of fragments ',20(1x,I3,'+'))
 433  format (' #',/,
     & ' # Too many OQS Localized MOs',/,
     & ' # Try again decreasing the critov() values',/,' #')
 434  format (' #',/,' # ulewisoqs.f: Unable to localize all the MOs',
     & /,' #',/,' # NON CONVERGED Open Quantum Systems (OQS) Analysis',
     & /,' # Try running again with a greater value of MXFRAG',
     & /,' # or decreasing the values of CRITOV(1..MXFRAG)',/,' #')
 435  format (' #',
     & /,' # ulewisoqs.f: !! WARNING: Too many localized MOs !',/,
     &   ' # ++++++++++++++++++++++++++',/,
     &   ' # Possible degenerate system',/,
     &   ' # ++++++++++++++++++++++++++')
*436  format (' #',
*    & /,' # ulewisoqs.f: !! WARNING: Too many localized MOs !',/,
*    &   ' # +++++++++++++++++++++++++++++++++++++++++++++++++',/,
*    &   ' # POSSIBLE DEGENERATE SYSTEM',/,
*    &   ' # ++++++++++++++++++++++++++',/,' # There are',I2,
*    &   ' Localized MOs with very close eigenvalues:',/,' # ',
 436  format (' #',
     & /,' # ulewisoqs.f: !! WARNING: Too many localized MOs !',/,
     &   ' # +++++++++++++++++++++++++++++++++++++++++++++++++',/,
     &   ' # POSSIBLE DEGENERATE SYSTEM in the search of',I2,
     & '-fragment MOs')
 100  format (' # ulewisoqs.f: File ',a,' NOT FOUND')
 1123 format (' # NMO in file ',a,'(',I4,') # Current value = ',I4)
 1124 format (' # Eigenvalue ',I4,' of sg() is negative, ',E17.10)
 28   format (' # The value of EPSNEG is ',E15.8)
 222  format (' # Negative eigenvalue ',E15.8,' > ',E15.8,
     &    ' is set to 0.0 diagonalizing fragment ',100I4)
 301  format (' # Number of MOs involving 1, 2, 3, ... fragments',/,
     &        ' # ',20(1x,I3))
 115  format (' # AOM between ',a,' localized MOs')
 101  format (1000(5(1x,F17.10),/))
 133  format (5(1x,F17.10))
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
 520  format (/' # Fragments expanded by each ',a,' OQS MO',/,
     &    1000(5(3x,F12.5),/))
 320  format (' # Type of bond between fragments',/,' # ',88('-'))
 321  format (' # Fragments ',2I4,' : ',I2,' Delocalized MOs'/, 
     &  ' # MO ',I2,' -> 2 * S_ii(A) * S_ii(B) = ',F15.8)
 322  format (3x,'MO ',I2,' -> 2 * S_ii(A) * S_ii(B) = ',F15.8)
 323  format (19x,'Approximate DI = ',F15.8,/,' #')
 65   format (' # Fragments ',2I3,'   Delta = ',F15.8)
 985  format (' #',/,
     & ' # +++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' # Effective Oxidation States (EOS) of the fragments',/,
     & ' # +++++++++++++++++++++++++++++++++++++++++++++++++')
 986  format (' # Fragment ',I2,'  EOS = ',F15.8)
 987  format (' # Total            = ',F15.8)
 171  format (' # Electron of spin-orbital ',I3,' goes to fragment ',I3)
 172  format (' # Number of ',a,' electrons in fragment ',I3,' is ',I3)
 690  format ('  # MO\FRAGMENT -> ',15(I3,4x),' ...')
 671  format (' # Eigenvalues of the accepted Localized MOs',/,
     &          1000(' # ',5(1x,F15.8),/))
 672  format (' # First neglected MO has eigenvalue',16x,F15.8)
9823  format (' # Found fragments with ',I2,' different sizes')
9824  format (' # Fragments with ',I2,' atoms =',I4)
 912  format (' # Fragment ',I2,' formed by ',I3, ' atoms: ',1000I4)
      end
