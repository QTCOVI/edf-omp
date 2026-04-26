c
c-----------------------------------------------------------------------
c
      module space_for_aomspin
      save
      real(kind=8), allocatable,dimension (:,:,:) :: aomalpha,aombeta
      end module space_for_aomspin
c
c-----------------------------------------------------------------------
c
      subroutine aomspin (ncent,nmo,aom,nalpha,nbeta)
      USE space_for_wfnbasis
      USE space_for_aomspin
      real   (kind=8) aom(ncent,nmo,nmo)
      integer(kind=4) ncent,nmo,nalpha,nbeta,i
c
      if (.not.allocated(aomalpha)) then
        allocate (aomalpha(ncent,nalpha,nalpha))
      endif
      if (.not.allocated(aombeta)) then
        allocate (aombeta (ncent,nbeta ,nbeta ))
      endif
      do i=1,ncent
        aomalpha(i,:,:) = aom(i,ialpha(1:nalpha),ialpha(1:nalpha))
        aombeta (i,:,:) = aom(i,ibeta (1:nbeta ),ibeta (1:nbeta ))
      enddo
      end subroutine aomspin
