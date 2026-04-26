c
c-----Variables used in the maximization of PSI^2
c
      parameter (maxoptim = 200)
      real*8 xyzel
      integer iopt
      common /optim/ xyzel(maxoptim,3),iopt(maxoptim,3)
