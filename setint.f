c
c-----------------------------------------------------------------------
c
      logical function setint (intvar, line, lp)
c
c.....setint - sets intvar to the integer value of the next word in
c     line.
c
c     If the next word is not a proper integer, setint returns .false.
c
      include           'implicit.inc'
c
      integer           intvar, lp
      character*(*)     line
c
      integer           i, ivar
      logical           isdig, done, iread
c
      character*(1)     ch, blank, newline, null
*mdc*if g77
*mdc*else
*      parameter ( blank = char(32), newline = char(10), null = char(0))
*mdc*endif
c
      isdig (ch) = ch.ge.'0' .and. ch.le.'9'
*mdc*if g77
      blank = char(32)
      newline = char(10)
      null = char(0)
*mdc*endif
c
      setint = .false.
      ll = len(line)
      done = .false.
      do while (lp.le.ll .and. .not.done)
         done = line(lp:lp) .ne. blank
         lp = lp + 1
      enddo
      if (.not.done) return
      lp = lp - 1
      i = lp
      if (line(i:i) .eq. '+') then
         isign = 1
         i = i + 1
      else if (line(i:i) .eq. '-') then
         isign = -1
         i = i + 1
      else
         isign = 1
      endif
      if (i.eq.ll+1) return
      ivar = 0
      done = .false.
      iread = .false.
      do while (.not.done .and. i.le.ll)
         if (isdig(line(i:i))) then
            ivar = 10*ivar + ichar(line(i:i)) - ichar('0')
            iread = .true.
         else if (index(blank//','//newline//null,line(i:i)).gt.0) then
            done = .true.
         else
            return
         endif
         i = i + 1
      enddo
      if (done) i = i - 1
      if (iread) then
         intvar = isign * ivar
         lp = i
         setint = .true.
      endif
      return
      end
