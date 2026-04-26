c
c-----------------------------------------------------------------------
c
      character*(*) function lower (string)
c
c.....lower - convert string to lowercase except where quoted
c           - string and lower will return the same
c
c.....WARNING: case dependent code!
c
      character*(*)     string
      character*(1)     quote
      integer           i, iadd
      logical           inquote
c
      quote = char(39)
      iadd = ichar('A') - ichar('a')
      i = 1
      inquote = .false.
      do i=1,len(string)
        if (string(i:i) .eq. quote) then
           inquote = .not. inquote
        else if (.not. inquote) then
           if (string(i:i).ge.'A' .and. string(i:i).le.'Z') then
              string(i:i) = char (ichar(string(i:i)) - iadd)
           endif
        endif
      enddo
      lower = string
      return
      end
