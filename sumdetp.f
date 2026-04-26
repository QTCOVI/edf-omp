c
c-----------------------------------------------------------------------
c
      subroutine sumdetp (epsdet,probal,ntimesa,ntimesb,nta,ntb,
     &   ncore,wca,wcb,ngroup,ndets,nmo,mulliken)
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
      integer(kind=4)  :: confai(nta)
      integer(kind=4)  :: confbi(ntb)
      integer(kind=4)  :: confaj(nta)
      integer(kind=4)  :: confbj(ntb)
      integer(kind=4)  :: wca(ntimesa,nmo)
      integer(kind=4)  :: wcb(ntimesb,nmo)
      real   (kind=8)  :: epsdet,probal
      logical    mulliken
c
      call timer (2,isumdet,'_sumdetp  ',-1)
c
c.....run over all pairs of determinants.
c
      do i=1,ncore
        confai(i)=i
        confbi(i)=i
        confaj(i)=i
        confbj(i)=i
      enddo
      probal=0d0
      do i=1,ndets
         nsabi=nsiga(i)*nsigb(i)
         cdprim=cdet(i)
         kwai=kwa(i)
         kwbi=kwb(i)
         confai(ncore+1:nta)=nconfa(kwai,1:nta-ncore)+ncore
         confbi(ncore+1:ntb)=nconfb(kwbi,1:ntb-ncore)+ncore
         jend=i
         if (mulliken) jend=ndets
         do j=1,jend
           nsabj=nsiga(j)*nsigb(j)
           cdsecond=cdet(j)
           cd=cdprim*cdsecond*nsabi*nsabj
           if (abs(cd).ge.abs(epsdet)) then
             if (i.ne.j.and.(.not.mulliken)) cd=cd+cd
             kwaj=kwa(j)
             kwbj=kwb(j)
             confaj(ncore+1:nta)=nconfa(kwaj,1:nta-ncore)+ncore
             confbj(ncore+1:ntb)=nconfb(kwbj,1:ntb-ncore)+ncore
             proba=0d0
             do ia=1,ntimesa
               do m=1,nta
                 do n=1,nta
                   idom=wca(ia,n)
                   oveaa(m,n)=sg(idom,confai(m),confaj(n))
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
                   ovebb(m,n)=sg(idom,confbi(m),confbj(n))
                 enddo
               enddo
               call detlapack (ovebb,indxb,ntb,info,probbx)
               probb=probb+probbx
             enddo
             probal=probal+cd*proba*probb
           endif
         enddo
      enddo
      call timer (4,isumdet,'_sumdetp  ',-1)
      return
      end
