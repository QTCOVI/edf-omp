c
c.....Read complex AOM (upper side triangle), normalize it
c     and construct the lower side triangle ensuring full hermiticity. 
c
      subroutine readcmplx (aom,ncent,nmo,lr,lerr)
      include    'constants.inc'
      integer ncent,nmo,lr,lerr,i,m,j
      complex*16 aom(ncent,nmo,nmo)
      real(kind=8), allocatable,dimension (:,:,:)  :: aomr
      real(kind=8)  aom1
      complex*16  aom1cmplx
c
      allocate (aomr(ncent,nmo,nmo))
      do i=1,ncent
        read (lr,*) iaom
        if (iaom.gt.ncent.or.iaom.le.0) then
          stop 'readcmplx.f: Fatal error reading complex aom() array'
        endif
        read (lr,80) ((aomr(iaom,m,j),m=1,j),j=1,nmo)
      enddo
      aom(1:ncent,1:nmo,1:nmo)=cmplx(aomr(1:ncent,1:nmo,1:nmo),0d0)
      deallocate (aomr)
c
c.....Normalize diagonal elements, setting to 0.0 the imaginary parts.
c
      do m=1,nmo
        aom1=zero
        do i=2,ncent
          aom1=aom1+dble(aom(i,m,m))
        enddo
        aom(1,m,m)=cmplx(one-aom1,zero)
      enddo

      do j=2,nmo
        do m=1,j-1
          aom1cmplx=cmplx(zero,zero)
          do i=2,ncent
            aom1cmplx=aom1cmplx+aom(i,m,j)
          enddo
          aom(1,m,j)=-aom1cmplx
        enddo
      enddo
c
c.....Fill in the lower part of AOM to ensure hermiticity
c
      do m=1,nmo
        do j=1,m-1
          aom(1:ncent,m,j)=conjg(aom(1:ncent,j,m))
        enddo
      enddo
c
*     do m=1,nmo
*       do j=1,m
*         aom1=zero
*         do i=1,ncent
*           aom1=aom1+dble(aom(i,m,j))
*         enddo
*         write (6,'(2I4,2F26.20)') m,j,aom1
*       enddo
*     enddo
      return
 80   format (6(1x,e16.10))
      end
