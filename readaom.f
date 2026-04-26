c
c.....Read the complex AOM written by DGrid (lower side triangle), 
c     normalize it and construct the upper side triangle ensuring 
c     full hermiticity. 
c
      subroutine readaom (aom,naom,ncent,nmo,lr,lerr)
      include    'constants.inc'
      include    'mline.inc'
*     parameter   (mline = 200)  
      integer naom,nmo,lr,lerr,i,m,j
      complex*16 aom(naom,nmo,nmo)
      real(kind=8), allocatable,dimension (:,:,:)  :: aomr,aomi
      real(kind=8)  aom1
      complex*16  aom1cmplx
      character*(mline) line
c
      allocate (aomr(naom,nmo,nmo))
      allocate (aomi(naom,nmo,nmo))

      nmobloc=nmo/10
      nmorest=nmo-nmobloc*10
      ibfin=nmobloc
      if (mod(nmo,10).ne.0) ibfin=ibfin+1

      do ntimes=1,2
        do i=1,naom
          inic0=0
          inic=10
          do ib=1,ibfin
            read (lr,'(a)') line
            read (lr,'(a)') line
            read (lr,'(a)') line
            read (lr,'(a)') line
            mas=1
            do m=inic0+1,nmo
              jini=inic0+1
              jfin=inic0+mas
              if (ntimes.eq.1) then
                read (lr,*) k,(aomr(i,m,j),j=jini,jfin)
              else
                read (lr,*) k,(aomi(i,m,j),j=jini,jfin)
              endif
              mas=min(mas+1,10)
            enddo
            inic0=inic
            inic=inic+10
          enddo
          read (lr,'(a)') line
          read (lr,'(a)') line
          read (lr,'(a)') line
        enddo
        if (ntimes.eq.1) then
          read (lr,'(a)') line
          read (lr,'(a)') line
          read (lr,'(a)') line
          read (lr,'(a)') line
          read (lr,'(a)') line
        endif
      enddo

      do m=1,nmo
        do j=1,m
          aom(1:naom,m,j)=cmplx(aomr(1:naom,m,j),aomi(1:naom,m,j))
        enddo
      enddo 
      deallocate (aomr,aomi)
c
c     In case that NCENT.EQ.NAOM ===> Normalize diagonal elements, 
c     setting to 0.0 the imaginary parts.
c
      if (naom.eq.ncent) then
        do m=1,nmo
          aom1=0d0
          do i=2,naom
            aom1=aom1+dreal(aom(i,m,m))
          enddo
          aom(1,m,m)=cmplx(one-aom1,zero)
        enddo
  
        do j=2,nmo
          do m=1,j-1
            aom1cmplx=cmplx(zero,zero)
            do i=2,naom
              aom1cmplx=aom1cmplx+aom(i,j,m)
            enddo
            aom(1,j,m)=-aom1cmplx
          enddo
        enddo
      endif
c
c.....Fill in the lower part of AOM to ensure hermiticity
c
      do m=1,nmo
        do j=1,m-1
          aom(1:naom,j,m)=conjg(aom(1:naom,m,j))
        enddo
      enddo
c
      return
 80   format (6(1x,e16.10))
 6    format (30(3I3,5x))
 7    format (30(2x,F16.10))
      end
