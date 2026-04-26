c
c-----------------------------------------------------------------------
c
c.....This routine determines the size occupied by a real*8 variable in 
c     a direct access file. This size is saved in the integer variable
c     RecLength, returned in common /sizew/ of 'length_of_a_record.inc'.
c
c-----------------------------------------------------------------------
c
      subroutine lengrec ()
      implicit real(kind=8) (a-h,o-z)
      include 'lengrec.inc'
      real(kind=8)  zero
      logical tested_size_8
      integer lopen
c
      RecLength = 2
      bytespow  = 63
      zero=0d0
 2    continue
      lopen=59
 4    lopen=lopen+1
      open (unit     =  lopen,
     &      file     = 'XxYyAaBb',
     &      access   = 'direct',
     &      recl     =  RecLength,
     &      err      =  4,
     &      form     = 'unformatted')
      write (lopen,rec=1,err=1) zero
      goto 3
 1    close (lopen)
      if (RecLength.ne.8) then
         RecLength = 8
         bytespow  = 31
         goto 2
      endif
 3    continue
      if (RecLength.ne.2.and.RecLength.ne.8) then
         stop 'lengrec.f: !! Error computing RecLength !! '
      else
         close (lopen,status='delete')
      endif
      end
