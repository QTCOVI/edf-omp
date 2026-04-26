
c
c-----------------------------------------------------------------------
c
      subroutine sdwdeta (proba,nta,wca,nmo)
c
      USE        space_for_conf
      USE        space_for_sgarea
      include   'implicit.inc'
      include   'constants.inc'
c
      real   (kind=8)  :: oveaa(nta,nta)
      integer(kind=4)  :: indxa(nta)
      integer(kind=4)  :: wca(nmo)
      real   (kind=8)  :: proba
c
      call timer (2,isumdet,'_sdwdeta  ',-1)
      do m=1,nta
        do n=1,nta
          oveaa(m,n)=sg(wca(n),nconfa(1,m),nconfa(1,n))
        enddo
      enddo
      call detlapack (oveaa,indxa,nta,info,proba)
      call timer (4,isumdet,'_sdwdeta  ',-1)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine sdwdetb (probb,ntb,wcb,nmo)
c
      USE        space_for_conf
      USE        space_for_sgarea
      include   'implicit.inc'
      include   'constants.inc'
c
      real   (kind=8)  :: ovebb(ntb,ntb)
      integer(kind=4)  :: indxb(ntb)
      integer(kind=4)  :: wcb(nmo)
      real   (kind=8)  :: probb
c
      call timer (2,isumdet,'_sdwdetb  ',-1)
      do m=1,ntb
        do n=1,ntb
          ovebb(m,n)=sg(wcb(n),nconfb(1,m),nconfb(1,n))
        enddo
      enddo
      call detlapack (ovebb,indxb,ntb,info,probb)
      call timer (4,isumdet,'_sdwdetb  ',-1)
      return
      end
