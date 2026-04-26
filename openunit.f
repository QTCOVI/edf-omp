c-----------------------------------------------------------------------
c
c-----This routine returns a 'save' value LU to open the file 'name' 
c     with that value for its logical unit number. The opening of the
c     file is also done in the routine.
c
      subroutine openunit (lu,name)
      integer(kind=4) lu
      character*(*) name
      logical exfil
c
      exfil=.true.
      name=trim(name)
      do while (exfil)
        inquire (unit=lu,opened=exfil)
        if (.not.exfil) then
          open (unit=lu,file=trim(name),iostat=ios)
          if (ios.ne.0) then
            write (0,*) ' # openunit.f: Error openning '//trim(name)
            stop
          endif
        else
          lu=lu+1
        endif
      enddo
      return
      end
