c
c-----------------------------------------------------------------------
c
      module space_for_dafh
      save
      integer(kind=4)  ndafh
      real(kind=8),    allocatable,dimension (:)   :: cutoff
      integer(kind=4), allocatable,dimension (:,:) :: ibcen
      integer(kind=4), allocatable,dimension (:)   :: nbcen
      logical,         allocatable,dimension (:)   :: deplet
      end module space_for_dafh
c
c-----------------------------------------------------------------------
c
      subroutine allocate_space_for_dafh (mxdafh,ncent)
      USE        space_for_dafh
      allocate (cutoff(mxdafh),stat=ier)
      if (ier.ne.0) stop 'dafhspace.f: Cannot allocate cutoff()'
      allocate (nbcen(ncent),stat=ier)
      if (ier.ne.0) stop 'dafhspace.f: Cannot allocate nbcen()'
      allocate (deplet(mxdafh),stat=ier)
      if (ier.ne.0) stop 'dafhspace.f: Cannot allocate deplet()'
      allocate (ibcen(mxdafh,ncent),stat=ier)
      if (ier.ne.0) stop 'dafhspace.f: Cannot allocate ibcen()'
      end subroutine allocate_space_for_dafh
c
c-----------------------------------------------------------------------
c
      subroutine deallocate_space_for_dafh
      USE        space_for_dafh
      deallocate (cutoff,stat=ier)
      if (ier.ne.0) stop 'dafhspace.f: Cannot deallocate cutoff()'
      deallocate (nbcen,stat=ier)
      if (ier.ne.0) stop 'dafhspace.f: Cannot deallocate nbcen()'
      deallocate (deplet,stat=ier)
      if (ier.ne.0) stop 'dafhspace.f: Cannot deallocate deplet()'
      deallocate (ibcen,stat=ier)
      if (ier.ne.0) stop 'dafhspace.f: Cannot deallocate ibcen()'
      end subroutine deallocate_space_for_dafh
