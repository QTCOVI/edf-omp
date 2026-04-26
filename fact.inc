c
c.....Factorials
c
      parameter (maxfac = 54)
      real (kind=8) fact(0:maxfac)
      real (kind=8) facd(-maxfac:maxfac)
      common/fact/fact,facd
