c
c-----------------------------------------------------------------------
c
      module space_for_rdm1
      save
      real(kind=8), allocatable,dimension (:,:)  :: c1et,c1ea,c1eb
      end module space_for_rdm1
c-----------------------------------------------------------------------
      subroutine allocate_space_for_rdm1 (nmo)
        USE space_for_rdm1
        integer(kind=4) :: nmo
        allocate (c1et(nmo,nmo)) 
        allocate (c1ea(nmo,nmo)) 
        allocate (c1eb(nmo,nmo)) 
      end subroutine allocate_space_for_rdm1
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_rdm1 ()
        USE space_for_rdm1
        deallocate (c1et) 
        deallocate (c1ea) 
        deallocate (c1eb) 
      end subroutine deallocate_space_for_rdm1
c-----------------------------------------------------------------------
