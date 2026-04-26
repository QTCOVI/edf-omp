      module space_for_bicen
      save
      real(kind=8), allocatable, dimension (:,:) :: cbrag
      real(kind=8), allocatable, dimension (:)   :: rmsurf,rbrag,charatm
      end module space_for_bicen
c-----------------------------------------------------------------------
      subroutine allocate_space_for_bicen (ncent)
      USE        space_for_bicen
      if (.not.allocated(cbrag)) then
        allocate (cbrag(ncent,ncent),stat=ier) 
        if (ier.ne.0) stop 'bicenarea.f: Cannot allocate cbrag()'
      endif
      if (.not.allocated(rmsurf)) then
        allocate (rmsurf(ncent),stat=ier) 
        if (ier.ne.0) stop 'bicenarea.f: Cannot allocate rmsurf()'
      endif
      if (.not.allocated(rbrag)) then
        allocate (rbrag(ncent),stat=ier) 
        if (ier.ne.0) stop 'bicenarea.f: Cannot allocate rbrag()'
      endif
      if (.not.allocated(charatm)) then
        allocate (charatm(ncent),stat=ier) 
        if (ier.ne.0) stop 'bicenarea.f: Cannot allocate charatm()'
      endif
      end subroutine allocate_space_for_bicen
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_bicen
      USE        space_for_bicen
      if (allocated(cbrag)) then
        deallocate (cbrag,stat=ier) 
        if (ier.ne.0) stop 'bicenarea.f: Cannot deallocate cbrag()'
      endif
      if (allocated(rmsurf)) then
        deallocate (rmsurf,stat=ier) 
        if (ier.ne.0) stop 'bicenarea.f: Cannot deallocate rmsurf()'
      endif
      if (allocated(rbrag)) then
        deallocate (rbrag,stat=ier) 
        if (ier.ne.0) stop 'bicenarea.f: Cannot deallocate rbrag()'
      endif
      if (allocated(charatm)) then
        deallocate (charatm,stat=ier) 
        if (ier.ne.0) stop 'bicenarea.f: Cannot deallocate charatm()'
      endif
      end subroutine deallocate_space_for_bicen
