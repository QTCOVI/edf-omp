!
!-----Computes using LAPACK the inverse of the a(N,N) matrix
!
      subroutine detlapack (a,ind,n,info,det)
      implicit none
      integer (kind=4) n,info,j
      real    (kind=8) signo,det
      integer (kind=4) ind(n)
      real    (kind=8) a(n,n)
!
!-----LU decomposition using Lapack
!
      call dgetrf (n,n,a,n,ind,info)
      if (info.gt.0) then
!
!-------The following warning is not issued because we are not 
!       going to solve a linear system of equations in this routine;
!       i.e. matrix 'A' has been put succesfully in the form
!    
!       A = P * L * U
!
!       where 'P' is a permutation matrix, 'L' is a lower triangular 
!       matrix with unit diagonal elements, and U is an upper triangular
!       matrix
!     
!       write (0,*) 'detlapack.f: U(i,i), with i=',info,' = 0'
!       write (0,*) 'detlapack.f: a() is a singular matrix'
      else if (info.lt.0) then
         write (0,*) 'detlapack.f: Input parameter ',-info,' is illegal'
      else
      endif
!
!-----Determinant
! 
      det=1.0d+00
      do j=1,n
        signo = 1.0d+00
        if (ind(j).ne.j) signo=-signo
        det = det * a(j,j) * signo
      enddo
      return
      end

