c
c-----------------------------------------------------------------------
c
      module space_for_conf
      save
      real(kind=8),  allocatable,dimension (:)   :: cdet
      integer, allocatable,dimension (:)   :: kwa,kwb,nsiga,nsigb
      integer, allocatable,dimension (:,:) :: nconfa,nconfb
      real(kind=8)   xnorm
      integer  nelal,nelbe,npa,npb
      end module space_for_conf
c-----------------------------------------------------------------------
      subroutine allocate_space_for_conf (ndets,ndetsa,ndetsb,nmo)
      USE        space_for_conf
      if (.not.allocated(cdet)) then
        allocate (cdet(ndets),stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'spaceconf.f: Cannot allocate cdet()'
        endif
      endif
      if (.not.allocated(kwa)) then
        allocate (kwa(ndets),stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'spaceconf.f: Cannot allocate kwa()'
        endif
      endif
      if (.not.allocated(kwb)) then
        allocate (kwb(ndets),stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'spaceconf.f: Cannot allocate kwb()'
        endif
      endif
      if (.not.allocated(nsiga)) then
        allocate (nsiga(ndets),stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'spaceconf.f: Cannot allocate nsiga()'
        endif
      endif
      if (.not.allocated(nsigb)) then
        allocate (nsigb(ndets),stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'spaceconf.f: Cannot allocate nsigb()'
        endif
      endif
      if (.not.allocated(nconfa)) then
        allocate (nconfa(ndetsa,nmo),stat=ier)
        if (ier.ne.0) stop 'spaceconf.f: Cannot allocate nconfa()'
      endif
      if (.not.allocated(nconfb)) then
        allocate (nconfb(ndetsb,nmo),stat=ier)
        if (ier.ne.0) stop 'spaceconf.f: Cannot allocate nconfb()'
      endif
      end subroutine allocate_space_for_conf
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_conf
      USE        space_for_conf
      if (allocated(cdet)) then
        deallocate (cdet,stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'spaceconf.f: Cannot deallocate cdet()'
        endif
      endif
      if (allocated(kwa)) then
        deallocate (kwa,stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'spaceconf.f: Cannot deallocate kwa()'
        endif
      endif
      if (allocated(kwb)) then
        deallocate (kwb,stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'spaceconf.f: Cannot deallocate kwb()'
        endif
      endif
      if (allocated(nsiga)) then
        deallocate (nsiga,stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'spaceconf.f: Cannot deallocate nsiga()'
        endif
      endif
      if (allocated(nsigb)) then
        deallocate (nsigb,stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'spaceconf.f: Cannot deallocate nsigb()'
        endif
      endif
      if (allocated(nconfa)) then
        deallocate (nconfa,stat=ier)
        if (ier.ne.0) stop 'spaceconf.f: Cannot deallocate nconfa()'
      endif
      if (allocated(nconfb)) then
        deallocate (nconfb,stat=ier)
        if (ier.ne.0) stop 'spaceconf.f: Cannot deallocate nconfb()'
      endif
      end subroutine deallocate_space_for_conf
