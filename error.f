c
c-----------------------------------------------------------------------
c
      subroutine error (routine,message,errortype)
c
      include             'stderr.inc'
      include             'error.inc'
c
      character*(*)       routine,message
      integer             errortype
      character*(20)      chtype
c
      if (errortype.eq.faterr) then
         chtype='fatal error'
      else if (errortype.eq.warning) then
         chtype='warning'
      else if (errortype.eq.noerr) then
         chtype='no error'
      else
         chtype='unknown error'
      endif
      write (stderr,100) trim(chtype),routine,message
      if (errortype.eq.faterr) stop
c
 100  format (/,1x,'*** '/,1x,'*** ',a,' at ',a,' :'/,
     .          1x,'*** ',a/,1x,'***')
c
      end
