      subroutine zgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(n)
      complex*16 a(lda,n)
      real(kind=8) z(n)
      real(kind=8) rcond
      real(kind=8) zlange
c
      real(kind=8) anorm
      complex*16 work(2*n)
      integer info
      real(kind=8) rwork(2*n)
c
      anorm = zlange('1',n,n,a,lda,z)
      call zgetrf (n,n,a,lda,ipvt,info)
      call zgecon ('1',n,a,lda,anorm,rcond,work,rwork,info)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine zgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(n),job
      complex*16 a(lda,n),b(n)
      integer info
      call zgetrs('N',n,1,a,lda,ipvt,b,lda,info)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine zgedi(a,lda,n,ipvt,det,work,job,determ)
      integer lda,n,ipvt(n),job
      real(kind=8) det(2),work(n)
      complex*16 a(lda,n)
      real(kind=8) ten
      integer i
      complex*16 determ

      determ=cmplx(1.0d+00,0.0d+00)
      do i=1,n
        if (ipvt(i) .ne. i) determ=-determ
        determ=determ*a(i,i)
      enddo
      return
c
c     compute determinant
c
*     if (job/10 .eq. 0) go to 70
*        det(1) = 1.0d0
*        det(2) = 0.0d0
*        ten = 10.0d0
*        do 50 i = 1, n
*           if (ipvt(i) .ne. i) det(1) = -det(1)
*           det(1) = dble(a(i,i))*det(1)
c        ...exit
*           if (det(1) .eq. 0.0d0) go to 60
*  10       if (dabs(det(1)) .ge. 1.0d0) go to 20
*              det(1) = ten*det(1)
*              det(2) = det(2) - 1.0d0
*           go to 10
*  20       continue
*  30       if (dabs(det(1)) .lt. ten) go to 40
*              det(1) = det(1)/ten
*              det(2) = det(2) + 1.0d0
*           go to 30
*  40       continue
*  50    continue
*  60    continue
*  70 continue
*     return
      end
c
c-----------------------------------------------------------------------
c
      subroutine zjacobi (a,b,work,w,n)
      complex*16 a(n,n),b(n,n),work(2*n-1)
      integer info,itype,lwork,n,i
      real(kind=8) rwork(3*n-2),w(n)

      itype = 1
      lwork = 2*n-1
      b(1:n,1:n)=(0.0d+00,0.0d+00)
      do i=1,n
        b(i,i)=(1.0d+00,0.0d+00)
      enddo
      call zhegv (itype,'V','U',N,A,N,B,N,W,work,lwork,rwork,info)
      if (info.ne.0) stop 'znetlib.f: INFO # 0 returned by zhegv'
      return
      end
