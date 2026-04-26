c
c-----------------------------------------------------------------------
c
      program dumi
c
c.....semilla - prints out the date and time at the beginning/end of the c     run.
c
      parameter (nchstr=8)
      character(len=30) :: date
      character*(nchstr) idum
      integer*4          num,atoi,nchaux
c
      date = fdate()
      idum = date(15:16)//date(18:19)//date(18:19)//date(15:16)
      nchaux = nchstr
      num = atoi (idum,nchaux)
      write (6,*) num
      end
