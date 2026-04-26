c
c-----------------------------------------------------------------------
c
      subroutine chsort (arr, iord, long, first, last, nstack)
c
c     character(len=long) version of the qcksort.f routine.
c
      include           'implicit.inc'
      include           'constants.inc'
      include           'stderr.inc'
      parameter         (m=7)
      character(len=*)  arr(*)
      character(len=long)  a
      dimension         iord(*),istack(nstack)
      integer           first,last
c
      jstack=0
      l=first
      ir=last
      fx=zero
10    if(ir-l.lt.m)then
        do j=l+1,ir
          na=iord(j)
          a=arr(na)
          do i=j-1,first,-1
            if(arr(iord(i)).le.a) go to 12
            iord(i+1)=iord(i)
          enddo
          i=first-1
12        iord(i+1)=na
        enddo
        if(jstack.eq.0) then
          return
        endif
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        i=l
        j=ir
        fx=mod(fx*fa+fc,fm)
        iq=l+(ir-l+1)*(fx*fmi)
        na=iord(iq)
        a=arr(na)
        iord(iq)=iord(l)
20      continue
21        if (j.ge.first) then 
             if (a.lt.arr(iord(j))) then
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
             if (a.gt.arr(iord(i))) then
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
        if (jstack.gt.nstack) then
          write (stderr,*) 'qqsort.f: !! Increase nstack value !!'
          stop
        endif
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
