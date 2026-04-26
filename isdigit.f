c
c-----------------------------------------------------------------------
c
      logical function isdigit (c)
c
c.....isdigit - return true if c is a digit
c
      character*(1)     c
c
      isdigit = c.ge.'0' .and. c.le.'9'
      return
      end
c...................................................................
