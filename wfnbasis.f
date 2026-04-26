c
c-----------------------------------------------------------------------
c
      module space_for_wfnbasis
      save
      real   (kind=8)   tote,gamma
      real   (kind=8),  allocatable, dimension (:,:)  :: xyz
      real   (kind=8),  allocatable, dimension (:,:)  :: sprim
      real   (kind=8),  allocatable, dimension (:)    :: oexp
      real   (kind=8),  allocatable, dimension (:)    :: occ,eorb,charge
      integer(kind=4),  allocatable, dimension (:)    :: icen,ityp
      integer(kind=4),  allocatable, dimension (:)    :: nprimc
      integer(kind=4),  allocatable, dimension (:)    :: ialpha,ibeta
      character*8,      allocatable, dimension (:)    :: atnam
      logical           sprimok
      end module space_for_wfnbasis
c-----------------------------------------------------------------------
      subroutine allocate_space_for_wfnbasis (nprims,nmo,ncent)
      USE        space_for_wfnbasis
c
      allocate (xyz(ncent,3),stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already allocated xyz() array ?'
      endif
c
      allocate (sprim(nprims,nprims),stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already allocated sprim() array ?'
      endif
c
      allocate (oexp(nprims),stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already allocated oexp() array ?'
      endif
c
      allocate (occ(nmo+nmo),stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already allocated occ() array ?'
      endif
c
      allocate (eorb(nmo),stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already allocated eorb() array ?'
      endif
c
      allocate (charge(ncent),stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already allocated charge() array ?'
      endif
c
      allocate (icen(nprims),stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already allocated icen() array ?'
      endif
c
      allocate (ityp(nprims),stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already allocated ityp() array ?'
      endif
c
      allocate (atnam(ncent),stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already allocated atnam() array ?'
      endif
c
      allocate (nprimc(ncent),stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already allocated nprimc() array ?'
      endif
c
      allocate (ialpha(nmo),stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already allocated ialpha() array ?'
      endif
c
      allocate (ibeta(nmo),stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already allocated ibeta() array ?'
      endif
c
      end subroutine allocate_space_for_wfnbasis
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_wfnbasis
      USE        space_for_wfnbasis
c
      if (allocated (xyz)) then
        deallocate (xyz,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated xyz() array ?'
        endif
      endif
c
      if (allocated (sprim)) then
      deallocate (sprim,stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already deallocated sprim() array ?'
      endif
      endif
c
      if (allocated (oexp)) then
        deallocate (oexp,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated oexp() array ?'
        endif
      endif
c
      if (allocated (occ)) then
        deallocate (occ,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated occ() array ?'
        endif
      endif
c
      if (allocated (eorb)) then
        deallocate (eorb,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated eorb() array ?'
        endif
      endif
c
      if (allocated (charge)) then
        deallocate (charge,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated charge() array ?'
        endif
      endif
c
      if (allocated (icen)) then
        deallocate (icen,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated icen() array ?'
        endif
      endif
c
      if (allocated (ityp)) then
        deallocate (ityp,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated ityp() array ?'
        endif
      endif
c
      if (allocated (atnam)) then
        deallocate (atnam,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated atnam() array ?'
        endif
      endif
c
      if (allocated (nprimc)) then
        deallocate (nprimc,stat=ier)
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated nprimc() array ?'
        endif
      endif
c
      deallocate (ialpha,stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already deallocated ialpha() array ?'
      endif
c
      deallocate (ibeta,stat=ier) 
      if (ier.ne.0) then
        stop 'wfnbasis.f: Already deallocated ibeta() array ?'
      endif
c
      end subroutine deallocate_space_for_wfnbasis
