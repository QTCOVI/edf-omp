c
c-----This routine takes the string variable 'totname' and suppress from 
c     it the last characters until a '.' charater is found. The result 
c     is returned in cutname. For instance, when  totname='file.wfn' the
c     returned result is cutname='file'
c
      subroutine topoint (totname,cutname)
      character*(*) totname,cutname
      character*1 ch,point,blanco
      integer*2 l,i

      l = leng(totname)
      point = '.'
      blanco = ' '
      ch = blanco
      do while (ch.ne.point.and.l.ge.1)
       ch=totname(l:l)
       l=l-1
      enddo
      if (l.eq.0) then
        stop ' # topoint.f: Input string TOTNAME does not contain a "."'
      else
        cutname=trim(totname(1:l))
      endif
      return
      end
