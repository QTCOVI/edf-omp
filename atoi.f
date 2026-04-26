c
c-----------------------------------------------------------------------
c
      integer*4 function atoi (string,nchstr)
c
c.....atoi   Convert ascii string to integer
c
      character*(*)     string
      integer           i,sign,nchstr
      logical           isdigit
c
      atoi=0
      i=1
      do while (string(i:i).eq.' ')
         i=i+1
      enddo
      sign=1
      if(string(i:i) .eq. '+' .or. string(i:i) .eq. '-') then 
         if (string(i:i) .eq. '-') sign=-1
         i=i+1
      endif
      do while (isdigit(string(i:i)))
         atoi=10*atoi+ichar(string(i:i))-ichar('0')
         i=i+1
         if (i.gt.nchstr) goto 1
      end do
 1    atoi=atoi*sign
      return
      end
