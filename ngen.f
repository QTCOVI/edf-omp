c
c-----------------------------------------------------------------------
c
c
c.....Extracts a subset of N elements from the set of the first M natural 
c     numbers. Each time the routine is called, it produces the next sub-
c     set of 'N' elements. If FIRST.EQ..TRUE. the output is 1,2,...,N.
c
      subroutine ngen (n,m,mmax,comb,first)
      integer(kind=4) comb(mmax),n,m,i,j,mmax
      logical first
c
      if (m.gt.mmax) then
        stop '# ngen.f: The value of MMAX must be increased'
      endif
      if (first) then
        first=.false.
        do i=1,n
          comb(i)=i
        enddo
      else
        i=n
 1      do while (i.ge.1)
          if (comb(i).lt.m) then
            comb(i)=comb(i)+1
            do j=i+1,n
              comb(j)=comb(j-1)+1
              if (comb(j).gt.m) then
                 i=i-1
                 goto 1
              endif
            enddo
            return
          else
            i=i-1
          endif
        enddo
      endif
      return
      end
