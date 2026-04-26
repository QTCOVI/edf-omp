c
c-----------------------------------------------------------------------
c
      subroutine iqcksort(iarr, iord, nmax, first, last)
c
c.....qcksort - sort the elements of ra in ascending order, using the
c     quicksort algorithm.
c
c.....input parameters:
c     iarr() ....... data to be sorted.
c     iord() ..... initial order of data in ra().
c     first .......... first element in ra() that will be analyzed.
c     last ........... las t element in ra() that will be analyzed.
c
c.....output parameters:
c     iord() ..... final order of data in ra().
c
c.....Maximum number of elements to be sorted depends on nstack value:
c     nstack...... ~2log2(last-first+1)
c
      implicit real(kind=8)   (a-h,o-z)
      integer           iarr(nmax),iord(nmax)
      integer           first, last
c
      parameter         (m=7, nstack=100,fm=7875d0,fa=211d0,fc=1663d0,
     &                  fmi=1.2698413d-4, zero=0d0)
      integer           istack(nstack)
c
      jstack=0
      l=first
      ir=last
      fx=zero
10    if(ir-l.lt.m)then
        do j=l+1, ir
          na=iord(j)
          a=iarr(na)
          do i=j-1, first, -1
            if(iarr(iord(i)).le.a) go to 12
            iord(i+1)=iord(i)
          enddo
          i=first-1
12        iord(i+1)=na
        enddo
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        i=l
        j=ir
        fx=mod(fx*fa+fc,fm)
        iq=l+(ir-l+1)*(fx*fmi)
        na=iord(iq)
        a=iarr(na)
        iord(iq)=iord(l)
20      continue
21        if (j.ge.first) then 
             if (a.lt.iarr(iord(j))) then
                 j=j-1
                 goto 21
             endif
          endif
          if(j.le.i)then
            iord(i)=na
            go to 30
          endif
          iord(i)=iord(j)
          i=i+1
22        if (i.le.last) then 
             if (a.gt.iarr(iord(i))) then
                i=i+1
                goto 22
             endif
          endif
          if(j.le.i)then
            iord(j)=na
            i=j
            go to 30
          endif
          iord(j)=iord(i)
          j=j-1
        go to 20
30      jstack=jstack+2
        if(jstack.gt.nstack) stop 'iqcksort: Increase nstack'
        if(ir-i.ge.i-l)then
          istack(jstack)=ir
          istack(jstack-1)=i+1
          ir=i-1
        else
          istack(jstack)=i-1
          istack(jstack-1)=l
          l=i+1
        endif
      endif
      go to 10
      end
