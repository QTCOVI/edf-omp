c
c-----------------------------------------------------------------------
c
      module space_for_maxwfnarea
      save
      real(kind=8),  allocatable,dimension (:,:)   :: xyzel
      logical, allocatable,dimension (:,:)   :: iopt
      real(kind=8)   epsdet
      integer  nmo,ncore,mal,mbe,lw,lr,lerr
      end module space_for_maxwfnarea
c-----------------------------------------------------------------------
      subroutine allocate_space_for_maxwfnarea (nel)
      USE        space_for_maxwfnarea
      allocate (xyzel(nel,3),stat=ier) 
      if (ier.ne.0) then
        stop 'maxwfnarea.f: Cannot allocate xyzel() array'
      endif
      allocate (iopt(nel,3),stat=ier) 
      if (ier.ne.0) then
        stop 'maxwfnarea.f: Cannot allocate iopt() array'
      endif
      end subroutine allocate_space_for_maxwfnarea
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_maxwfnarea ()
      USE        space_for_maxwfnarea
      deallocate (xyzel,stat=ier)
      if (ier.ne.0) then
        stop 'maxwfnarea.f: Cannot  deallocate xyzel array ()'
      endif
      deallocate (iopt,stat=ier)
      if (ier.ne.0) then
        stop 'maxwfnarea.f: Cannot  deallocate iopt array ()'
      endif
      end subroutine deallocate_space_for_maxwfnarea
