c
c-----------------------------------------------------------------------
c
      subroutine semilla (num)
c
c.....semilla - prints out the date and time at the beginning/end of the
c     run.
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
      return
      end







*     program test_fdate
*     integer(8) :: i, j
*     character(len=30) :: date
*     call fdate(date)
*     print *, 'Program started on ', trim(date)
*     do i = 1, 100000000 ! Just a delay
*         j = i * i - i
*     end do
*     call fdate(date)
*     print *, 'Program ended on ', trim(date)
*     end  program test_fdate
*Program started on 
*Mon Jan 17 21:08:31 2022
