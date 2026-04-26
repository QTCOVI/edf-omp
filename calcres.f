      logical function calcres (ngroup,occa,occb,maxa,mina,maxb,minb)
      implicit none
      integer(kind=4) :: ngroup,k
      integer(kind=4) :: occa(ngroup)
      integer(kind=4) :: occb(ngroup)
      integer(kind=4) :: maxa(ngroup)
      integer(kind=4) :: mina(ngroup)
      integer(kind=4) :: maxb(ngroup)
      integer(kind=4) :: minb(ngroup)
      logical         :: docalc
    
      calcres=.true.
      do k=1,ngroup
        if (occa(k).gt.maxa(k).or.occb(k).gt.maxb(k)) then
          calcres=.false.
          return
        endif
        if (occa(k).lt.mina(k).or.occb(k).lt.minb(k)) then
          calcres=.false.
          return
        endif
      enddo
      return
      end
