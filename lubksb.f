c
c-----------------------------------------------------------------------
c
      subroutine lubksb (a, n, np, indx, b)
c
c.....lubksb: back substitution of vector b into the LU decomposed
c     matrix a
c
c.....Input parameters:
c
c     a(,) ....... LU decomposed matrix
c     n, np ...... real and declared dimensions of the matrix
c     indx() ..... index vector of the LU decomposed rows
c     b() ........ column vector to be substituted
c
c.....Output parameters:
c
c     b() ........ solution vector
c
      implicit real*8 (a-h,o-z)
c
      dimension         a(np,np), indx(n), b(n)
c
c.....skip zero elements on b, and do the forward substitution
c
      ii = 0
      do i = 1, n
         ll = indx(i)
         sum = b(ll)
         b(ll) = b(i)
         if (ii.ne.0) then
            do j = ii, i-1
               sum = sum - a(i,j) * b(j)
            enddo
         else if (sum.ne.0d0) then
            ii = i
         endif
         b(i) = sum
      enddo
c
c.....now, complete with the back substitution
c
      do i = n, 1, -1
         sum = b(i)
         if (i.lt.n) then
            do j = i+1, n
               sum = sum - a(i,j) * b(j)
            enddo
         endif
         b(i) = sum / a(i,i)
      enddo
      return
      end
