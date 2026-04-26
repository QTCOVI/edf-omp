c
c-----------------------------------------------------------------------
c
      subroutine itdafhc (sg,dafh,c,mxcent,covx,maxbond,nloc,nelhalf,
     &  udat,wfnfile,maxcritov,critov,ixcent,largwr,warndafh,
     &  skiph,lw)
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
c.....dafh()                           INPUT/OUTPUT
c
c     Group overlap matrix between input MO's when entering the 
c     routine or between output MO's at the end of the routine.
c
c.....c()                            OUTPUT 
c
c     Rotation matrix giving the transformed MO's from the input MO's. 
c
c.....nmo                            INPUT
c
c     Number of MO's.
c
c.....ncent                          INPUT
c
c     Number of fragments 
c
c.....critov()                       INPUT
c
c     If the overlap of an MO with itself in a n-center domain is greater 
c     than CRITOV(n), the routine considers that this MO is fully localized 
c     on that domain. 
c
c.....largwr
c
c     When LARGWR=.TRUE. a slightly larger output is written.
c
c.....lw                             INPUT
c
c     Output logical unit
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
      real(kind=8)    sg(ncent,nmo,nmo),dafh(ncent,nmo,nmo)
      real(kind=8),   allocatable, dimension (:)     :: eigs,work
      real(kind=8),   allocatable, dimension (:)     :: eigloc
      real(kind=8),   allocatable, dimension (:)     :: sumcent,elec
      real(kind=8),   allocatable, dimension (:,:)   :: v
      real(kind=8),   allocatable, dimension (:,:)   :: aomt,dafht
      real(kind=8),   allocatable, dimension (:,:)   :: eta
      real(kind=8),   allocatable, dimension (:,:,:) :: sg0,sgx,dafh0
      real(kind=8),   allocatable, dimension (:,:)   :: dafhloc,dafhloc0
      real(kind=8),   allocatable, dimension (:,:)   :: ctmp
      integer(kind=4),allocatable, dimension (:)     :: iord
      integer(kind=4),allocatable, dimension (:)     :: ncenloc
      integer(kind=8),allocatable, dimension (:)     :: ntuples
      integer(kind=4),allocatable, dimension (:)     :: inside,outside
      real   (kind=8),allocatable, dimension (:)     :: ocup
c
      integer(kind=4),allocatable, dimension (:,:)   :: wh
      integer(kind=4),allocatable, dimension (:)     :: coord
      integer(kind=4),allocatable, dimension (:,:)   :: madis
      parameter (maxcoord = 20)
c
      logical :: first
      logical warndafh,warnofmo,is_hydrogen,skiph
c
c-----Data used in the orthonormalization process
c
      real(kind=8),    allocatable,dimension (:,:)   :: cc
c
      real(kind=8) critov(maxcritov)
      logical ixcent(maxcritov)
      logical   largwr,dow
c
      integer(kind=4) udat,udatnw
      character(len=*) wfnfile
      character(len=mline) wfnloc
      character(len=mline) locform(nmo)
c
      call timer (2,ipid,'_itdafhc  ',-1)
c
c.....Allocate arrays
c
      n = nmo
      if (.not.allocated(eigs)) then
        allocate (eigs(n),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate eigs()'
      endif
      if (.not.allocated(eigloc)) then
        allocate (eigloc(n),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate eigloc()'
      endif
      if (.not.allocated(ntuples)) then
        allocate (ntuples(mxcent),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate ntuples()'
      endif
      if (.not.allocated(inside)) then
        allocate (inside(ncent),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate inside()'
      endif
      if (.not.allocated(outside)) then
        allocate (outside(ncent),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate outside()'
      endif
      if (.not.allocated(work)) then
        allocate (work(3*n-1),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate work()'
      endif
      if (.not.allocated(v)) then
        allocate (v(n,n),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate v()'
      endif
      if (.not.allocated(ctmp)) then
        allocate (ctmp(n,n),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate ctmp()'
      endif
      if (.not.allocated(ncenloc)) then
        allocate (ncenloc(n),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate ncenloc()'
      endif
      if (.not.allocated(madis)) then
        allocate (madis(ncent,ncent),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate madis()'
      endif
      if (.not.allocated(coord)) then
        allocate (coord(ncent),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate coord()'
      endif
      if (.not.allocated(wh)) then
        allocate (wh(ncent,maxcoord),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate wh()'
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
c-----Print the input DAFH and AOM arrays
c
      if (largwr) then
        if (.not.allocated(aomt)) then
          allocate (aomt(n,n),stat=ier)
          if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate aomt()'
        endif
        if (.not.allocated(dafht)) then
          allocate (dafht(n,n),stat=ier)
          if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate dafht()'
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

        write (lw,*) '#'
        dafht=zero
        write (lw,*) '#'
        write (lw,*) '# Domain Averaged Fermi Hole (DAFH)'
        write (lw,*) '#'
        do i=1,ncent
          write (lw,*) '     ATOM ',i
          do j=1,n
            do k=1,n
              dafht(j,k)=dafht(j,k)+dafh(i,j,k)
            enddo
            write (lw,'(6(1x,F15.8))') (dafh(i,j,k),k=1,j)
          enddo
        enddo
        write (lw,*) '#'
        write (lw,*) '#    TOTAL'
        do j=1,n
          write (lw,'(6(1x,F15.8))') (dafht(j,k),k=1,j)
        enddo
        if (allocated(aomt)) then
          deallocate (aomt,stat=ier)
          if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate aomt()'
        endif
        if (allocated(dafht)) then
          deallocate (dafht,stat=ier)
          if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate dafht()'
        endif
      endif
c
c.....Start the process.
c
      natoms=1
      do while (natoms.le.ncenteff) 
        xtuples=1d0
        do i=0,natoms-1
          xtuples=xtuples*dble(ncent-i)/dble(natoms-i)
        enddo
        ntuples(natoms)=int(xtuples)
        natoms=natoms+1
      enddo     
      if (.not.allocated(iord)) then
        allocate (iord(ncent),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate iord()'
      endif
      if (.not.allocated(dafh0)) then
        allocate (dafh0(ncent,n,n),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate dafh0()'
      endif
      if (.not.allocated(sg0)) then
        allocate (sg0(ncent,n,n),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate sg0()'
      endif
      if (.not.allocated(sgx)) then
        allocate (sgx(ncent,n,n),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate sgx()'
      endif
      if (.not.allocated(dafhloc)) then
        allocate (dafhloc(n,n),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate dafhloc()'
      endif
      if (.not.allocated(dafhloc0)) then
        allocate (dafhloc0(n,n),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate dafhloc0()'
      endif
      if (.not.allocated(sumcent)) then
        allocate (sumcent(ncenteff),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate sumcent()'
      endif
      if (.not.allocated(elec)) then
        allocate (elec(ncent),stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate elec()'
      endif
      dafh0    = dafh 
      dafhloc  = zero
      dafhloc0 = zero
      ctmp     = zero
      n3       = 3*n-1
      nloc     = 0
      natoms   = 1
      sumcent  = zero
      warndafh = .false.
      do while (natoms.le.ncenteff.and.nloc.lt.nmo)
        if (ixcent(natoms)) then
          write (lw,11) natoms
          first=.true.
          lis=1
          dow=.true.
          forall (i=1:ncent) iord(i)=i
          do while (lis.le.ntuples(natoms).and.nloc.lt.nmo)
            call ngen (natoms,ncent,ncent,inside,first)
            if (natoms.eq.1) then
              is_hydrogen = nint(charge(inside(1))).eq.1
              if (is_hydrogen.and.skiph) goto 111
            elseif (natoms.eq.2) then
              iat1=inside(1)
              iat2=inside(2)
              m12=madis(iat1,iat2)
              if (m12.lt.1.and.m12.gt.maxbond) goto 111
            elseif (natoms.ge.3) then
              ic1=coord(inside(1))
              ic2=coord(inside(2))
              ic3=coord(inside(3))
              if (ic1.le.2.and.ic2.le.2.and.ic3.le.2) goto 111
            endif
            v=zero
            do i=1,natoms
              v=v+dafh0(inside(i),:,:)
            enddo
            v = v-dafhloc0
            eigs = 0d0
            call dsyev ('V','U',n,v,n,eigs(1:n),work,n3,info)
            call errdsyev (info,'itdafhc.f')
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
*               if (nloc.gt.nelhalf) then
                if (nloc.gt.nmo) then
                  write (lw,*) '#'
                  write (lw,*) '# Too many localized DAFH MOs. Try'//
     &                         ' again decreasing the critov() values.'
                  write (lw,*) '#'
                  stop              
                endif
                ncenloc(nloc)=natoms
                eigloc(nloc)=eigs(i)
                ctmp(:,nloc) = v(:,i)
                sumcent(natoms) = sumcent(natoms) + eigs(i)
                do j=1,n
                  do k=1,n
                    prod=ctmp(j,nloc)*ctmp(k,nloc)*eigloc(nloc)
                    dafhloc(j,k)=dafhloc(j,k)+prod
                  enddo
                enddo
                if (dow) then
                  write (lw,113) natoms,critov(natoms)
                  dow = .false.
                endif
                if (nloc.le.nelhalf) then
                  write (lw,116) nloc,eigs(i),
     &              (atnam(inside(j))(3:4),inside(j),j=1,natoms)
                else
                  write (locform(nloc-nelhalf),116) nloc,eigs(i),
     &              (atnam(inside(j))(3:4),inside(j),j=1,natoms)
                endif
              endif
            enddo
 111        continue
            lis=lis+1
          enddo
        endif
        dafhloc0=dafhloc
        natoms=natoms+1 
      enddo
      if (nloc.lt.nelhalf) then
        write (lw,*)'# '
        write (lw,*)'# itdafhc.f: Unable to localize all the MOs'
        write (lw,*)'# '
        write (lw,302)
        warndafh = .true.
      elseif (nloc.eq.nelhalf) then
        write (lw,300)
      else
        write (lw,*)
        write (lw,*)'# itdafhc.f: ! WARNING: Too many localized MOs !'
        write (lw,*)'# '
        do i=nelhalf+1,nloc
          write (lw,'(a)') trim(locform(i-nelhalf))
        enddo
*       nloc=nelhalf
      endif
      c(1:nloc,:)=transpose(ctmp(:,1:nloc))
      if (largwr) then
        write (lw,20)
        do i=1,nloc
          write (lw,33) i
          write (lw,30) (c(i,j),j=1,n)
        enddo
      endif
c
c-----Write diagonal AOM elements in each atom or fragment
c
      sumeigen=sum(eigloc(1:nloc))
      do i=1,ncent
        sg0(i,1:nloc,1:nloc) = matmul(c(1:nloc,:),
     &  matmul(sg(i,:,:),transpose(c(1:nloc,:))))
      enddo
      write (lw,69) sumeigen
      do i=1,ncenteff
        if (abs(sumcent(i)).gt.zero) write (lw,70) i,sumcent(i)
      enddo
      write (lw,68)
      write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
      elec = zero
      do ip=1,nloc
        do j=1,ncent
          elec(j) = elec(j) + sg0(j,ip,ip)
        enddo
        write (lw,410,advance='no') ip
        if (ncent.gt.10) then
          write (lw,41 ) (100d0*sg0(k,ip,ip),k=1,10)
          write (lw,411) (100d0*sg0(k,ip,ip),k=11,ncent)
          write (lw,*) 
        else
          write (lw,41) (100d0*sg0(k,ip,ip),k=1,ncent)
        endif
      enddo
      sg(:,1:nloc,1:nloc)=sg0(:,1:nloc,1:nloc)
      write (lw,71) 
      do i=1,ncent
        write (lw,710) i,elec(i)
      enddo
      write (lw,711) sum(elec(1:ncent))
c
c-----Writing output corresponding to the orthonormalized DAFH MOs
c     using the symmetric Lowdin method.
c
CCCCC     if (nloc.eq.n) then
CCCCC       write (lw,112) 'BEGINNING'
CCCCC       if (.not.allocated(cc)) then
CCCCC         allocate (cc(nmo+nmo,nprims),stat=ier)
CCCCC         if (ier.ne.0) stop 'edf.f: Cannot allocate cc()'
CCCCC       endif
CCCCC       do i=1,nmo
CCCCC         do ip=1,nprims
CCCCC           tmp=zero
CCCCC           do j=1,nmo
CCCCC             tmp=tmp+coef(j+nmo,ip)*c(i,j)
CCCCC           enddo
CCCCC           cc(i,ip)=tmp
CCCCC           cc(i+nmo,ip)=tmp
CCCCC         enddo
CCCCC       enddo
CCCCC       write (lw,1000)
CCCCC       if (.not.allocated(eta)) then
CCCCC         allocate (eta(nmo,nmo),stat=ier)
CCCCC         if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate eta()'
CCCCC       endif
CCCCC       call ofmortho (cc,eta,nmo,nprims,warnofmo,lw)
CCCCC       if (warnofmo) stop '# itdafhc.f: Singular matrix in ofmortho.f'
CCCCC       do i=1,ncent
CCCCC         sgx(i,:,:) = matmul(eta,matmul(sg(i,:,:),transpose(eta)))
CCCCC       enddo
CCCCC       write (lw,68)
CCCCC       write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
CCCCC       elec = zero
CCCCC       do ip=1,n
CCCCC         do j=1,ncent
CCCCC           elec(j) = elec(j) + sgx(j,ip,ip)
CCCCC         enddo
CCCCC         write (lw,410,advance='no') ip
CCCCC         if (ncent.gt.10) then
CCCCC           write (lw,41) (100d0*sgx(k,ip,ip),k=1,10)
CCCCC           write (lw,411) (100d0*sgx(k,ip,ip),k=11,ncent)
CCCCC           write (lw,*) 
CCCCC         else
CCCCC           write (lw,41) (100d0*sgx(k,ip,ip),k=1,ncent)
CCCCC         endif
CCCCC       enddo
CCCCC       write (lw,334) 'orthonormalized Localized MOs'
CCCCC       call newdrvmult (nmo,nprims,ncent,cc(1:nmo,:),lw)
CCCCC       inmo=nmo+nmo
CCCCC       inpr=nprims
CCCCC       occ(1:nmo)=zero
CCCCC       wfnloc=trim(wfnfile)//"dafh-ortho"    
CCCCC       udatnw=udat+10
CCCCC       open (unit=udatnw,file=trim(wfnloc),status='unknown')
CCCCC       write (lw,351) trim(wfnloc)
CCCCC       call wrtwfn (udat,udatnw,0,ifilc,cc,inmo,inpr)
CCCCC       write (lw,112) 'END'
CCCCC       if (allocated(eta)) then
CCCCC         deallocate (eta,stat=ier)
CCCCC         if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate eta()'
CCCCC       endif
CCCCC     else
        write (lw,112) 'BEGINNING'
        if (.not.allocated(cc)) then
          allocate (cc(nloc,nprims),stat=ier)
          if (ier.ne.0) stop 'edf.f: Cannot allocate cc()'
        endif
        do i=1,nloc
          do ip=1,nprims
            tmp=zero
            do j=1,nmo
              tmp=tmp+coef(j+nmo,ip)*c(i,j)
            enddo
            cc(i,ip)=tmp
          enddo
        enddo
        write (lw,1000)
        if (.not.allocated(eta)) then
          allocate (eta(nloc,nloc),stat=ier)
          if (ier.ne.0) stop ' # itdafhc.f: Cannot allocate eta()'
        endif
        call ofmortho (cc,eta,nloc,nprims,warnofmo,lw)
        if (warnofmo) stop ' # itdafhc.f: Singular matrix in ofmortho.f'
        do i=1,ncent
          sgx(i,1:nloc,1:nloc) = matmul(eta,
     &    matmul(sg(i,1:nloc,1:nloc),transpose(eta)))
        enddo
        write (lw,68)
        write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
        elec = zero
        do ip=1,nloc
          do j=1,ncent
            elec(j) = elec(j) + sgx(j,ip,ip)
          enddo
          write (lw,410,advance='no') ip
          if (ncent.gt.10) then
            write (lw,41) (100d0*sgx(k,ip,ip),k=1,10)
            write (lw,411) (100d0*sgx(k,ip,ip),k=11,ncent)
            write (lw,*) 
          else
            write (lw,41) (100d0*sgx(k,ip,ip),k=1,ncent)
          endif
        enddo
        write (lw,334) 'orthonormalized Localized MOs'
        call newdrvmult (nloc,nprims,ncent,cc,lw)
        inpr=nprims
        occ(1:nloc)=zero
        wfnloc=trim(wfnfile)//"dafh-ortho"    
        udatnw=udat+10
        open (unit=udatnw,file=trim(wfnloc),status='unknown')
        write (lw,351) trim(wfnloc)
        allocate (ocup(nloc),stat=ier)
        if (ier.ne.0) stop '# itdafhc.f: Cannot allocate ocup()'
        ocup(1:nloc)=occ(1:nloc)
        call cdafhmos (udat,udatnw,0,ocup,ocup,cc,nloc,inpr)
        write (lw,112) 'END'
        deallocate (ocup,stat=ier)
        if (ier.ne.0) stop '# itdafhc.f: Cannot deallocate ocup()'
        if (allocated(eta)) then
          deallocate (eta,stat=ier)
          if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate eta()'
        endif
CCCCC     endif 
c
      if (allocated(sg0)) then
        deallocate (sg0,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate sg0()'
      endif
      if (allocated(sgx)) then
        deallocate (sgx,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate sgx()'
      endif
      if (allocated(dafh0)) then
        deallocate (dafh0,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate dafh0()'
      endif
      if (allocated(eigs)) then
        deallocate (eigs,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate eigs()'
      endif
      if (allocated(iord)) then
        deallocate (iord,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate iord()'
      endif
      if (allocated(ntuples)) then
        deallocate (ntuples,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate ntuples()'
      endif
      if (allocated(work)) then
        deallocate (work,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate work()'
      endif
      if (allocated(v)) then
        deallocate (v,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate v()'
      endif
      if (allocated(dafhloc)) then
        deallocate (dafhloc,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate dafhloc()'
      endif
      if (allocated(dafhloc0)) then
        deallocate (dafhloc0,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate dafhloc0()'
      endif
      if (allocated (ctmp)) then
        deallocate (ctmp,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate ctmp()'
      endif
      if (allocated (inside)) then
        deallocate (inside,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate inside()'
      endif
      if (allocated (outside)) then
        deallocate (outside,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate outside()'
      endif
      if (allocated(ncenloc)) then
        deallocate (ncenloc,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate ncenloc()'
      endif
      if (allocated(eigloc)) then
        deallocate (eigloc,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate eigloc()'
      endif
      if (allocated(sumcent)) then
        deallocate (sumcent,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate sumcent()'
      endif
      if (allocated(elec)) then
        deallocate (elec,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate elec()'
      endif
      if (allocated(madis)) then
        deallocate (madis,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate madis()'
      endif
      if (allocated(coord)) then
        deallocate (coord,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate coord()'
      endif
      if (allocated(wh)) then
        deallocate (wh,stat=ier)
        if (ier.ne.0) stop ' # itdafhc.f: Cannot deallocate wh()'
      endif
      write (lw,10) '  E N D  '
      return
      call timer (4,ipid,'_itdafhc  ',-1)
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
 71   format (' # ',/,
     &  ' # Electrons associated to successful DAFH diagonalizations')
 710  format (' # Atom ',I4,'   Electrons = ',F18.8)
 711  format (' # SUM  ',17x,'= ',F18.8)
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
      end
