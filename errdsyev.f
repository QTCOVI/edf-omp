      subroutine errdsyev (info,routine)
      implicit none
      integer(kind=4) info
      character(len=*) routine
 
      if (info.eq.0) then
        return
      elseif (info.gt.0) then
        write (0,1) trim(routine),'Illegal argument calling DSYEV'
        stop 
      else
        write (0,1) trim(routine),'DSYEV routine failed to converge'
        stop 
      endif
 1    format (' # ',a,': ',a)
      end
