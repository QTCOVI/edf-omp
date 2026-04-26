!
!-----Computes the LU decomposition of the a(N,N) matrix using LAPACK
!
!     DGETRF computes an LU factorization of a general M-by-N matrix A
!     using partial pivoting with row interchanges.

!     The factorization has the form
!        A = P * L * U
!     where P is a permutation matrix, L is lower triangular with unit
!     diagonal elements (lower trapezoidal if m > n), and U is upper
!     triangular (upper trapezoidal if m < n).
!
      subroutine lulapack (a,ind,n,info)
      implicit none
      integer (kind=4) n,info,j
      integer (kind=4) ind(n)
      real    (kind=8) a(n,n)
!
!-----LU decomposition using Lapack
!
      call dgetrf (n,n,a,n,ind,info)

      if (info.gt.0) then
         write (0,*) 'lulapack.f: U(i,i), with i=',info,', is exactly 0'
         write (0,*) 'lulapack.f: a() is a singular matrix'
      else if (info.lt.0) then
         write (0,*) 'lulapack.f: Input parameter ',-info,', is illegal'
      else
      endif
      return
      end
