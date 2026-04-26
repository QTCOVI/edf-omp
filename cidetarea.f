c
c-----------------------------------------------------------------------
c
      module space_for_cidet
      save
      real(kind=8),  allocatable, dimension (:)  :: cidet
      end module space_for_cidet
c-----------------------------------------------------------------------
      subroutine allocate_space_for_cidet (nelact)
      USE        space_for_cidet
      allocate (cidet(0:nelact),stat=ier) 
      if (ier.ne.0) then
        write (0,*) 'cidetarea.f: Already allocated cidet() ?'
      endif
      end subroutine allocate_space_for_cidet
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_cidet
      USE        space_for_cidet
      if (allocated (cidet)) then
        deallocate (cidet,stat=ier)
        if (ier.ne.0) stop 'cidetarea.f: Cannot  deallocate cidet ()'
      endif
      end subroutine deallocate_space_for_cidet
