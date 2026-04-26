c
c-----Prints Intersection of ordered arrays of integer numbers a[] and b[]
c     The nc common elements are returned in c[]
c
      subroutine intersection (a,b,c,na,nb,nab,nc)
      implicit none
      integer a(na),b(nb),c(nab),na,nb,nab,nc,i,j,k

      i=1
      j=1
      k=0
      do while (i.le.na .and. j.le.nb)
        if (a(i) < b(j)) then
          i=i+1
        elseif (b(j) < a(i)) then
          j=j+1
        elseif (a(i) == b(j)) then
          k=k+1
          c(k)=a(i)
          i=i+1
        endif
      enddo
      nc=k
      return
      end
