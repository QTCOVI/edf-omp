c
c-----------------------------------------------------------------------
c
      module space_for_wfncoef
      save
      real(kind=8),  allocatable, dimension (:,:)  :: coef
      end module space_for_wfncoef
c-----------------------------------------------------------------------
      subroutine allocate_space_for_wfncoef (nprims,nmo)
      USE        space_for_wfncoef
      allocate (coef(nmo+nmo,nprims),stat=ier) 
      if (ier.ne.0) then
        write (0,*) 'wfncoef.f: Already allocated coef ?'
      endif
      end subroutine allocate_space_for_wfncoef
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_wfncoef
      USE        space_for_wfncoef
      if (allocated (coef)) then
        deallocate (coef,stat=ier)
        if (ier.ne.0) stop 'wfncoef.f: Cannot  deallocate coef ()'
      endif
      end subroutine deallocate_space_for_wfncoef
