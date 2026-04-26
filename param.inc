c
c.....global parameters
c
c     parameter (maxtype      =    35)  ! S,P,D,F, and G primitites allowed.
      parameter (maxtype      =    56)  ! S,P,D,F,G, and H primitites allowed.
      parameter (nstack       =  4000)  ! Maximum number of probabilites
c                                       | that allow its ordering by de-
c                                       | creasing value.
c
c.....maxndet = Maximum allowed value for the dimensions of cj(),
c     cjcja(), iord(), and icj() arrays in 'lmowfn.f' = Maximum number dets
c
      parameter (maxndet      = 50000) 
      parameter (pi           = 3.14159265358979324D+00)
      common    /datai/ nlm(maxtype,3)
      logical icorr,rhf,rohf,uhf,cciqa,onlynat,thereisEDF,woccup
      common /logic / icorr,rhf,rohf,uhf,cciqa,onlynat,thereisEDF,woccup
