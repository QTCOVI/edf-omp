c
c-----------------------------------------------------------------------
c
      subroutine itdafh (sg,c,mxcent,covx,maxbond,udat,ifilc,wfnfile,
     &   maxcritov,critov,dafhselect,ndafhsel,idafhsel,
     &   ixcent,largwr,warndafh,skiph,lw)
c
c.....Performs the localization of a set of NMO molecular orbitals into 
c     NCENT atoms by diagonalizing the AOM in single atoms and pairs of
c     atoms. After each diagonalization, the MOs which are localized in 
c     the atom or pair of atoms where the AOM was computed are excluded. 
c     The proccess goes on (progressively decreasing the CRITOV 
c     parameter) until there are NOT non localized MOs.
c                            
c                            
c-----PARAMETERS--------------------------------------------------------
c
c.....sg()                           INPUT/OUTPUT
c
c     Group overlap matrix between input MO's when entering the 
c     routine or between output MO's at the end of the routine.
c
c.....c()                            OUTPUT 
c
c     Rotation matrix giving the transformed MO's from the input MO's. 
c
c.....mxcent                         INPUT
c
c     Maximum number of centers (n) explored in the search of localized 
c     (nc-2e) MOs
c
c.....covx                           INPUT
c
c     Variable that determines whether two atoms A and B are bonded or
c     not. A a B are assumed to be bonded if Rij < (RA+RB)*COVX, where
c     RA and RB are the covalent radius stored in the array covrad()
c     of 'connect.f' routine
c
c.....maxbond                        INPUT
c
c     'maxbond' represents the maximum number of bonds necessary to 
c     travel between atoms A and B which is actually explored; i.e. if 
c     going from 'A' to 'B' requires more than 'maxbond' bonds the pair
c     A-B is not explored. The default value of 'maxbond' is 1, i.e. 
c     only directly bonded atoms are considered. The maximum value of 
c     'maxbond' is 'maxmaxb = 4'. This variable is only relevant in the
c     search of 2c-2e bonds.
c
c.....udat                           INPUT
c
c     Logical unit of the WFN file
c
c.....ifilc                          INPUT
c
c     Logical unit of the file containing the coefficients of Slater
c     determinants multi-det WFNs.
c
c.....wfnfile                        INPUT
c
c     Name of the file containing the WFN.
c
c.....critov()                       INPUT
c
c     critov(n) is the threshold value used to discriminate whether a MO
c     is localized or not in a set of 'n' centers when searching localized 
c     (nc-2e) MOs.
c
c.....dafhselect() 
c
c     When dafhselect(n) is .true. the search of localized nc-2e MOs is
c     carried out only for those n-tuples that are given by the idafhsel()
c     integer values (see below). If ixcent(n) is .false. the value of
c     dafhselect(n) is irrelevant.
c
c.....ndafhsel()
c
c     When dafhselect(n) is .true. ndafhsel(n) is the number of n-tuples
c     for which the search of localized nc-2e MOs is carried out. When
c     dafhselect(n) is .false.  the value of ndafhsel(n) is irrelevant.
c
c-----idafhsel(n,k,1...n)
c     When afhselect(n) is .true. idafhsel(n,k,1...n) are the indices
c     of the 'n' atoms for which the search of localized nc-2e MOs is
c     carried out. The index 'k' goes from 1 to ndafhsel(n)
c
c.....ixcent()                       INPUT
c
c     If ixcent(n) = .true. searches of localized (nc-2e) MOs are
c     performed. Otherwise, these searches are skipped.
c
c.....largwr                         INPUT
c
c     When largwr = .true. a slightly larger output is written.
c
c.....warndafh                       OUTPUT
c
c     It will be .true. at the end if the maximum possible number of
c     localized MOs has not been found. Otherwise, it will be .false.
c
c.....skiph                          INPUT
c
c     If skiph = .true. one-center two-electron searches (1c-2e) in
c     hydrogen atoms are skipped. Otherwise, these searches are done. 
c
c.....lw                             INPUT
c
c     Logical output unit
c
c-----------------------------------------------------------------------
c                         
      USE      space_for_wfncoef
      USE      space_for_wfnbasis
      include 'implicit.inc'
      include 'constants.inc'
      include 'param.inc'
      include 'wfn.inc'
      include 'error.inc'
      include 'mline.inc'
      real(kind=8)    c(nmo,nmo)
      real(kind=8)    sg(ncent,nmo,nmo)
      real(kind=8),   allocatable, dimension (:)     :: eigs,work
      real(kind=8),   allocatable, dimension (:)     :: eigloc
      real(kind=8),   allocatable, dimension (:)     :: sumcent,elec
      real(kind=8),   allocatable, dimension (:,:)   :: v
      real(kind=8),   allocatable, dimension (:,:)   :: aomt
      real(kind=8),   allocatable, dimension (:,:)   :: eta
      real(kind=8),   allocatable, dimension (:,:,:) :: sg0,sgx
      real(kind=8),   allocatable, dimension (:,:)   :: sloc,sloc0
      real(kind=8),   allocatable, dimension (:,:)   :: ctmp
      integer(kind=4),allocatable, dimension (:)     :: iord
      integer(kind=4),allocatable, dimension (:)     :: ncenloc
      integer(kind=8),allocatable, dimension (:)     :: ntuples
      integer(kind=4),allocatable, dimension (:)     :: inside
c
      integer(kind=4),allocatable, dimension (:,:)   :: wh
      integer(kind=4),allocatable, dimension (:)     :: coord
      integer(kind=4),allocatable, dimension (:,:)   :: madis
      parameter (maxcoord = 20)
c
      logical warndafh,warnofmo,is_hydrogen,skiph,stilltuples,doslt
c
c-----Data used in the orthonormalization process
c
      real(kind=8),    allocatable,dimension (:,:)   :: cc
c
      real(kind=8) critov(maxcritov)
      logical ixcent(maxcritov)
      logical          dafhselect(maxcritov)
      integer(kind=4) ndafhsel   (maxcritov)
      integer(kind=4) idafhsel   (maxcritov,maxcritov,maxcritov)
      logical largwr,dow,exfil,first
      integer(kind=4) udat,udatnw
      character(len=*) wfnfile
      character(len=mline) wfnloc
c
      call timer (2,ipid,'_itdafh   ',-1)
c
c.....Allocate arrays
c
      n = nmo
      if (.not.allocated(eigs)) then
        allocate (eigs(n),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate eigs()'
      endif
      if (.not.allocated(eigloc)) then
        allocate (eigloc(n),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate eigloc()'
      endif
      if (.not.allocated(ntuples)) then
        allocate (ntuples(mxcent),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate ntuples()'
      endif
      if (.not.allocated(inside)) then
        allocate (inside(ncent),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate inside()'
      endif
      if (.not.allocated(work)) then
        allocate (work(3*n-1),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate work()'
      endif
      if (.not.allocated(v)) then
        allocate (v(n,n),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate v()'
      endif
      if (.not.allocated(ctmp)) then
        allocate (ctmp(n,n),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate ctmp()'
      endif
      if (.not.allocated(ncenloc)) then
        allocate (ncenloc(n),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate ncenloc()'
      endif
      if (.not.allocated(madis)) then
        allocate (madis(ncent,ncent),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate madis()'
      endif
      if (.not.allocated(coord)) then
        allocate (coord(ncent),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate coord()'
      endif
      if (.not.allocated(wh)) then
        allocate (wh(ncent,maxcoord),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate wh()'
      endif
      ncenteff = min(mxcent,ncent)
      write (lw,10) ' B E G I N  '
      write (lw,21) ncenteff,(critov(j),j=1,ncenteff)
      if (     skiph) write (lw,51) '     '
      if (.not.skiph) write (lw,51) ' NOT '
      write (lw,52) maxbond
c
c-----Obtain coordinations, connectivities, and distance arrays
c
      call connect (lw,covx,ncent,maxcoord,largwr,nbonds,coord,wh,madis)
c
c-----Largest number of bonds between two 'connected' atoms
c
      longest_chain=maxval(madis)
c
c-----Print the input AOM arrays
c
      if (largwr) then
        if (.not.allocated(aomt)) then
          allocate (aomt(n,n),stat=ier)
          if (ier.ne.0) stop 'itdafh.f: Cannot allocate aomt()'
        endif
        aomt=zero
        write (lw,*) '# Atomic Overlap matrix (AOM)'
        write (lw,*) '#'
        do i=1,ncent
          write (lw,*) '     ATOM ',i
          do j=1,n
            do k=1,n
              aomt(j,k)=aomt(j,k)+sg(i,j,k)
            enddo
            write (lw,'(6(1x,F15.8))') (sg(i,j,k),k=1,j)
          enddo
        enddo
        write (lw,*) '#'
        write (lw,*) '#    TOTAL'
        do j=1,n
          write (lw,'(6(1x,F15.8))') (aomt(j,k),k=1,j)
        enddo
        if (allocated(aomt)) then
          deallocate (aomt,stat=ier)
          if (ier.ne.0) stop 'itdafh.f: Cannot deallocate aomt()'
        endif
      endif
c
c.....Start the process.
c
      natoms=1
      do while (natoms.le.ncenteff) 
        xtuples=1.0d+00
        do i=0,natoms-1
          xtuples=xtuples*dble(ncent-i)/dble(natoms-i)
        enddo
        ntuples(natoms)=int(xtuples)
        natoms=natoms+1
      enddo     
      if (.not.allocated(iord)) then
        allocate (iord(ncent),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate iord()'
      endif
      if (.not.allocated(sg0)) then
        allocate (sg0(ncent,n,n),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate sg0()'
      endif
      if (.not.allocated(sgx)) then
        allocate (sgx(ncent,n,n),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate sgx()'
      endif
      if (.not.allocated(sloc)) then
        allocate (sloc(n,n),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate sloc()'
      endif
      if (.not.allocated(sloc0)) then
        allocate (sloc0(n,n),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate sloc0()'
      endif
      if (.not.allocated(sumcent)) then
        allocate (sumcent(ncenteff),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate sumcent()'
      endif
      if (.not.allocated(elec)) then
        allocate (elec(ncent),stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot allocate elec()'
      endif
      sg0      = sg 
      sloc     = zero
      sloc0    = zero
      ctmp     = zero
      n3       = 3*n-1
      nloc     = 0
      natoms   = 1
      sumcent  = zero
      warndafh = .false.
      do while (natoms.le.ncenteff.and.nloc.lt.nmo)
        doslt = dafhselect(natoms)
        if (ixcent(natoms)) then
          write (lw,11) natoms
          dow = .true.
          lis = 1
          if (doslt) then
            stilltuples = lis.le.ndafhsel(natoms)
            write (lw,*) '# n-tuples of atoms to be examined'
            do j=1,ndafhsel(natoms)
              write (lw,'(20I4)')(idafhsel(natoms,j,i),i=1,natoms)
            enddo
          else
            write (lw,*) '# All n-tuples of atoms will be examined'
            first=.true.
            forall (i=1:ncent) iord(i)=i
            stilltuples = lis.le.ntuples(natoms)
          endif
          do while (stilltuples.and.nloc.lt.nmo)
            if (doslt) then
              do i=1,natoms
                inside(i) = idafhsel(natoms,lis,i)
              enddo
            else
              call ngen (natoms,ncent,ncent,inside,first)
            endif
            if (natoms.eq.1) then
              is_hydrogen = nint(charge(inside(1))).eq.1
              if (is_hydrogen.and.skiph) goto 111
            elseif (natoms.eq.2) then
              iat1=inside(1)
              iat2=inside(2)
              m12=madis(iat1,iat2)
              if (m12.lt.1.and.m12.gt.maxbond) goto 111
            elseif (natoms.ge.3) then
            endif
            v=zero
            do i=1,natoms
              v=v+sg0(inside(i),:,:)
            enddo
            v = v-sloc0
            eigs = 0.0d+00
            call dsyev ('V','U',n,v,n,eigs(1:n),work,n3,info)
            call errdsyev (info,'itdafh.f')
            if (largwr) then
              write (lw,*) '#'
              write (lw,432) (inside(i),i=1,natoms)
              write (lw,*) '#'
              write (lw,'(6(1x,F15.8))') (eigs(k),k=1,n)
            endif
            do i=n,1,-1
              if (eigs(i).lt.critov(natoms)) then
                exit
              else
                nloc = nloc + 1
                if (nloc.gt.nmo) then
                  write (lw,*) '#'
                  write (lw,*) '# Too many localized DAFH MOs. Try'//
     &                         ' again decreasing the critov() values.'
                  write (lw,*) '#'
                  stop              
                endif
                ncenloc(nloc)   = natoms
                eigloc(nloc)    = eigs(i)
                ctmp(:,nloc)    = v(:,i)
                sumcent(natoms) = sumcent(natoms) + eigs(i)
                do j=1,n
                  do k=1,n
                    prod=ctmp(j,nloc)*ctmp(k,nloc)*eigloc(nloc)
                    sloc(j,k)=sloc(j,k)+prod
                  enddo
                enddo
                if (dow) then
                  write (lw,113) natoms,critov(natoms)
                  dow = .false.
                endif
                if (eigs(i).gt.1.0d+00) eigs(i)=1.0d+00  ! Avoid eigenvalues > 1.0
                write (lw,116) nloc,eigs(i),
     &            (atnam(inside(j))(3:4),inside(j),j=1,natoms)
              endif
            enddo
 111        continue
            lis=lis+1
            if (doslt) then
              stilltuples = lis.le.ndafhsel(natoms)
            else
              stilltuples = lis.le.ntuples(natoms)
            endif
          enddo
        endif
        sloc0=sloc
        natoms=natoms+1 
      enddo
      if (nloc.lt.n) then
        write (lw,*)'# '
        write (lw,*)'# itdafh.f: Unable to localize all the MOs'
        write (lw,*)'# '
        write (lw,302)
        warndafh = .true.
        goto 2222
      else
        write (lw,300)
        c=transpose(ctmp)
        if (largwr) then
          write (lw,20)
          do i=1,n
            write (lw,33) i
            write (lw,30) (c(i,j),j=1,n)
          enddo
        endif
      endif
c
c-----Write diagonal AOM elements in each atom or fragment
c
      if (nloc.eq.n) then
        sumeigen=sum(eigloc(1:n))
        do i=1,ncent
          sg0(i,:,:) = matmul(c,matmul(sg(i,:,:),transpose(c)))
        enddo
c
c-------Avoid by hand diagonal AOM elements negative of greater than 1.0 
c
        do i=1,ncent
          do j=1,n
            if (sg0(i,j,j).gt.1.0d+00) sg0(i,j,j)=1.0d+00
            if (sg0(i,j,j).lt.0.0d+00) sg0(i,j,j)=0.0d+00
          enddo
        enddo
c
c-------Print the output AOM arrays
c
        if (largwr) then
          if (.not.allocated(aomt)) then
            allocate (aomt(n,n),stat=ier)
            if (ier.ne.0) stop 'itdafh.f: Cannot allocate aomt()'
          endif
          aomt=zero
          write (lw,*) '# Atomic Overlap matrix (AOM)'
          write (lw,*) '#'
          do i=1,ncent
            write (lw,*) '     ATOM ',i
            do j=1,n
              do k=1,n
                aomt(j,k)=aomt(j,k)+sg0(i,j,k)
              enddo
              write (lw,'(6(1x,F15.8))') (sg0(i,j,k),k=1,j)
            enddo
          enddo
          write (lw,*) '#'
          write (lw,*) '#    TOTAL'
          do j=1,n
            write (lw,'(6(1x,F15.8))') (aomt(j,k),k=1,j)
          enddo
          if (allocated(aomt)) then
            deallocate (aomt,stat=ier)
            if (ier.ne.0) stop 'itdafh.f: Cannot deallocate aomt()'
          endif
        endif

        write (lw,69) sumeigen
        do i=1,ncenteff
          if (abs(sumcent(i)).gt.zero) write (lw,70) i,sumcent(i)
        enddo
        write (lw,68)
        write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
        elec = zero
        do ip=1,n
          do j=1,ncent
            elec(j) = elec(j) + sg0(j,ip,ip)
          enddo
          write (lw,410,advance='no') ip
          if (ncent.gt.10) then
            write (lw,41 ) (100.0d+00*sg0(k,ip,ip),k=1,10)
            write (lw,411) (100.0d+00*sg0(k,ip,ip),k=11,ncent)
            write (lw,*) 
          else
            write (lw,41) (100.0d+00*sg0(k,ip,ip),k=1,ncent)
          endif
        enddo
        sg=sg0
        write (lw,71) 
        do i=1,ncent
          write (lw,710) i,elec(i)
        enddo
      endif
c
c-----Writing output corresponding to the orthonormalized DAFH MOs
c     using the symmetric Lowdin method.
c
      if (nloc.eq.n) then
        write (lw,112) 'BEGINNING'
        if (.not.allocated(cc)) then
          allocate (cc(nmo+nmo,nprims),stat=ier)
          if (ier.ne.0) stop 'edf.f: Cannot allocate cc()'
        endif
        do i=1,nmo
          do ip=1,nprims
            tmp=zero
            do j=1,nmo
              tmp=tmp+coef(j+nmo,ip)*c(i,j)
            enddo
            cc(i,ip)=tmp
            cc(i+nmo,ip)=tmp
          enddo
        enddo
        write (lw,1000)
        if (.not.allocated(eta)) then
          allocate (eta(nmo,nmo),stat=ier)
          if (ier.ne.0) stop 'itdafh.f: Cannot allocate eta()'
        endif
        call ofmortho (cc,eta,nmo,nprims,warnofmo,lw)
        if (warnofmo) then
          write (lw,*) '# itdafh.f: Singular matrix in ofmortho.f'
          goto 123
        endif
        do i=1,ncent
          sgx(i,:,:) = matmul(eta,matmul(sg0(i,:,:),transpose(eta)))
        enddo
        write (lw,68)
        write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
        elec = zero
        do ip=1,n
          do j=1,ncent
            elec(j) = elec(j) + sgx(j,ip,ip)
          enddo
          write (lw,410,advance='no') ip
          if (ncent.gt.10) then
            write (lw,41) (100.0d+00*sgx(k,ip,ip),k=1,10)
            write (lw,411) (100.0d+00*sgx(k,ip,ip),k=11,ncent)
            write (lw,*) 
          else
            write (lw,41) (100.0d+00*sgx(k,ip,ip),k=1,ncent)
          endif
        enddo
        write (lw,334) 'orthonormalized Localized MOs'
        call newdrvmult (nmo,nprims,ncent,cc(1:nmo,:),lw)
        inmo=nmo+nmo
        inpr=nprims
CCCC    occ(1:nmo)=zero
        wfnloc=trim(wfnfile)//"dafh-ortho"    
        udatnw=udat+10
        open (unit=udatnw,file=trim(wfnloc),status='unknown')
        write (lw,351) trim(wfnloc)
        call wrtwfn (udat,udatnw,0,ifilc,cc,inmo,inpr)
c
c.......Write AOM file with the orthonormalized MOs
c
        lu19    =  19
        exfil   = .true.
        wfnloc  = trim(wfnloc )//'.aom'
        do while (exfil)
          inquire (unit=lu19,opened=exfil)
          if (.not.exfil) then
            open (unit=lu19,file=trim(wfnloc),iostat=ios)
            if (ios.ne.0) then
              write (0,*) '# itdafh.f: Error openning '//trim(wfnloc)
              stop
            endif
          else
            lu19=lu19+1
          endif
        enddo
        write (lw,3511) trim(wfnloc)
        do i=1,ncent
          write (lu19,'(I4,a)') i,' <--- AOMLOC-ortho in this center'
          write (lu19,80) ((sgx(i,m,j),m=1,j),j=1,nmo)
        enddo
 123    write (lw,112) 'END'
        if (allocated(eta)) then
          deallocate (eta,stat=ier)
          if (ier.ne.0) stop 'itdafh.f: Cannot deallocate eta()'
        endif
      endif 
c
 2222 continue
      if (allocated(sg0)) then
        deallocate (sg0,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate sg0()'
      endif
      if (allocated(sgx)) then
        deallocate (sgx,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate sgx()'
      endif
      if (allocated(eigs)) then
        deallocate (eigs,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate eigs()'
      endif
      if (allocated(iord)) then
        deallocate (iord,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate iord()'
      endif
      if (allocated(ntuples)) then
        deallocate (ntuples,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate ntuples()'
      endif
      if (allocated(work)) then
        deallocate (work,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate work()'
      endif
      if (allocated(v)) then
        deallocate (v,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate v()'
      endif
      if (allocated(sloc)) then
        deallocate (sloc,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate sloc()'
      endif
      if (allocated(sloc0)) then
        deallocate (sloc0,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate sloc0()'
      endif
      if (allocated (ctmp)) then
        deallocate (ctmp,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate ctmp()'
      endif
      if (allocated (inside)) then
        deallocate (inside,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate inside()'
      endif
      if (allocated(ncenloc)) then
        deallocate (ncenloc,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate ncenloc()'
      endif
      if (allocated(eigloc)) then
        deallocate (eigloc,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate eigloc()'
      endif
      if (allocated(sumcent)) then
        deallocate (sumcent,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate sumcent()'
      endif
      if (allocated(elec)) then
        deallocate (elec,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate elec()'
      endif
      if (allocated(madis)) then
        deallocate (madis,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate madis()'
      endif
      if (allocated(coord)) then
        deallocate (coord,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate coord()'
      endif
      if (allocated(wh)) then
        deallocate (wh,stat=ier)
        if (ier.ne.0) stop 'itdafh.f: Cannot deallocate wh()'
      endif
      write (lw,10) '  E N D  '
      return
      call timer (4,ipid,'_itdafh   ',-1)
c
c.....Formats
c
 410  format (' # ',I4)
 20   format (' # Rotation Matrix (READ BY ROWS)',
     &        '    DAFH MO_i =  SUM_j CANMO_j C_ji')
 30   format (6(1x,F12.6))
 41   format (1000(10x,10(1x,F5.1,'%')))
 411  format (1000(17x,10(1x,F5.1,'%'),/))
 300  format (' #',/,' # ',70('-'),/,
     &  ' # CONVERGED (1--N)-CENTER DAFH ANALYSIS',/,' #')
 302  format (' #',/,' # ',70('-'),/,
     &        ' # NON CONVERGED (1--N)-CENTER DAFH ANALYSIS',/,' #')
 68   format (/' # Degree of localization of each MO in each atom')
 33   format (' # MO ',I3)
 113  format (/,1x,'# ',78('-'),/,2x,I2,
     & '-center DAFH analysis, Critical overlap = ',F12.5,/,
     & ' # ',78('-'))
 116  format (4x,'MO ',I4,' with eigenvalue = ',F12.5,
     & ') localized on atom(s)',10(1x,'(',A2,I2,')'))
 21   format (' # CRITOV (1 ..',I2,')  = ',10(1x,F5.3))
 10   format (/,' # ',70('-'),/,' # ',a,
     & '(1--N)-C E N T E R   D A F H   A N A L Y S I S',/,
     & ' # ',70('-'))
 11   format (' #',/,' # DAFH analysis on ',I2,' centers',/,' #')
 69   format (' #',/,' # SUM of eigenvalues = ',F18.8,/,' #')
 70   format (' # ',7x,I4,'-center = ',F18.8)
 71   format (' # ',/,' # Electron charge of all the atoms')
 710  format (' # Atom ',I4,'   Electrons = ',F18.8)
 334  format (/,' # Second & fourth moment orbital spread of ',a,/,3x,
     &'I.M. Hoyvik & P. Jorgensen [Chem. Rev. 116, 3306-3326 (2016)]',
     &/,3x,
     & 'I.M. Hoyvik, B. Jansik & P. Jorgensen [JCP 137, 2224114 (2012)'
     & /,3x,61('-'),/,4x,'Orbital',1x,'sigma_2^p',5x,'sigma_4^p',
     & 5x,'beta_p',8x,'beta_p^4',/,3x,61('-'))
 351  format (/1x,"# Writing file ","'",a,"' with orthonormalized MOs")
 1000 format (' # DAFH MOs are orthonormalized using the Lowdin method')
 112  format (2(' #',/),' #',90('-'),/,' # ',a,
     & ' OF RESULTS WITH ORTHONORMALIZED OPTIMIZED MOs',/,' #',90('-'))
 51   format (' # In single atom analyses, H atoms are',a,'skipped')
 52   format (
     & " # In 'atom1-atom2' searches, maximum distance index is ",I2)
 432  format (' # Eigenvalues of the DAFH of atoms ',20(1x,I3,'+'))
 3511 format (1x,"#",/,1x,"# Writing file ","'",a,"' with",
     & " AOMs over DAFH Localized+orthogonalized MOs")
 80   format (6(1x,e16.10))
      end
