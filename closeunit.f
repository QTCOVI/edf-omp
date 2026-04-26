c-----------------------------------------------------------------------
c
c-----If the file with unit number LU is open, this routine closes it,
c     removing or not removing the file depending if remove=.true. or 
c     .false, respectively.
c
c-----If the file with unit number LU is not open, nothing is done.
c
c-----------------------------------------------------------------------
c
      subroutine closeunit (lu,remove)
      integer(kind=4) lu
      logical exfil,remove

      inquire (unit=lu,opened=exfil)
      if (exfil) then
        if (remove) then
          close (lu,status='delete')
        else
          close (lu)
        endif
      endif
      return
      end
