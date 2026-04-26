c
c-----------------------------------------------------------------------
c
      subroutine sdwsumdet (epsdet,probal,ntimesa,ntimesb,nta,ntb,
     &   wca,wcb,ngroup,nmo)
c
c.......................................................................
c     !!!! This routine should only be used when there are the same
c          number of ALPHA and BETA electrons.
c
      USE        space_for_conf
      USE        space_for_sgarea
      include   'implicit.inc'
      include   'constants.inc'
c
      real   (kind=8)  :: oveaa(nta,nta)
      real   (kind=8)  :: ovebb(ntb,ntb)
      real   (kind=8)  :: ww(nmo)
      integer(kind=4)  :: indxa(nta)
      integer(kind=4)  :: indxb(ntb)
      integer(kind=4)  :: wca(ntimesa,nmo)
      integer(kind=4)  :: wcb(ntimesb,nmo)
      real   (kind=8)  :: epsdet,probal
c
      call timer (2,isumdet,'_sdwsumdet',-1)
c
c.....run over all pairs of determinants.
c
      proba=0d0
      do ia=1,ntimesa
        do m=1,nta
          do n=1,nta
            idom=wca(ia,n)
            oveaa(m,n)=sg(idom,nconfa(1,m),nconfa(1,n))
          enddo
        enddo
        call detlapack (oveaa,indxa,nta,info,probax)
        proba=proba+probax
      enddo
      probb=0d0
      do ib=1,ntimesb
        do m=1,ntb
          do n=1,ntb
            idom=wcb(ib,n)
            ovebb(m,n)=sg(idom,nconfb(1,m),nconfb(1,n))
          enddo
        enddo
        call detlapack (ovebb,indxb,ntb,info,probbx)
        probb=probb+probbx
      enddo
      probal=proba*probb
      call timer (4,isumdet,'_sdwsumdet',-1)
      return
      end
