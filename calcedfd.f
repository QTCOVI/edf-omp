c
c-----------------------------------------------------------------------
c
      subroutine calcedfd 
     &  ( epsdet,probcut,nproba,nprobb,nprob,ngroup,ndets,nmo,ncore,
     &    nelact,nel,moval,ival,mocogrp,nelv,malv,mbev,ifilc,stdout,
     &    stderr,wfnfile,sg,resnca,resncb,largwr,orderp)
c
c.......................................................................
c     !!!! This routine can be used when the number of ALPHA and BETA 
c     electrons is different.
c
      USE        space_for_cidet
      include   'implicit.inc'
      include   'param.inc'
      include   'constants.inc'
      include   'lengrec.inc'
      include   'mline.inc'
c
      real(kind=8), allocatable,dimension (:,:,:) :: ovea
      real(kind=8), allocatable,dimension (:,:)   :: am,bm,oveaa
      real(kind=8), allocatable,dimension (:,:)   :: tpowa,tpowb
      real(kind=8), allocatable,dimension (:,:)   :: xlis
      real(kind=8), allocatable,dimension (:)  :: proba,probb,probt,pnew
      real(kind=8), allocatable,dimension (:)  :: w,ww,xli
      real(kind=8), allocatable,dimension (:)  :: p1,p1a,p1b,p2,d1,p3,d2
      real(kind=8), allocatable,dimension (:)  :: p2aa,p2bb,p2ab,p2ba
      real(kind=8), allocatable,dimension (:)  :: d1aa,d1bb,d1ab,d1ba
      real(kind=8), allocatable,dimension (:)  :: probord
      integer, allocatable,dimension (:,:)   :: resnc
      integer, allocatable,dimension (:)     :: ipvta,ipvtb,indxa,ind
      integer, allocatable,dimension (:)     :: indxb
      integer, allocatable,dimension (:)     :: iorda,iordb
      integer, allocatable,dimension (:)     :: npop
      integer, allocatable,dimension (:)     :: ordia,ordib,ordja,ordjb
      integer, allocatable,dimension (:)     :: nordia,nordib
      integer, allocatable,dimension (:)     :: kwa,kwb,nsiga,nsigb
      integer, allocatable,dimension(:)      :: idifa,idifb
      integer, allocatable,dimension (:)     :: ioprob
      integer, allocatable,dimension (:,:)   :: resncord
c
      real(kind=8)    sg(ngroup,nmo,nmo)
      real(kind=8)    dumi,random,deter(2)
      integer    resnca(nproba,ngroup),resncb(nprobb,ngroup)
      integer    stdout,stderr
      integer    ig(100)
      integer    mocogrp(ngroup)
      integer*4  idum
      logical    inlist,orderp
      logical    largwr
      integer    ival(moval)
*     parameter  (mline=200)
      character*(mline) wfnfile
      real(kind=8)     indet
c
      call timer (2,ipid1,'_calcedfd ',-1)
c
c.....Analyze the wave function and determine the different ALPHA/BETA 
c     configurations which appear on it. Write these configurations in
c     files 'confsa.dat' and 'confsb.dat'. Besides this, the arrays 
c     kwa() and kwb() store the order number of the ALPHA and BETA 
c     configuration appearing in each Slater determinant.
c
      allocate (ordia(nmo),ordja(nmo),iorda(nmo))
      allocate (ordib(nmo),ordjb(nmo),iordb(nmo))
      allocate (nordia(nmo),nordib(nmo))
      allocate (kwa(ndets),kwb(ndets))
      allocate (nsiga(ndets),nsigb(ndets))
      allocate (idifa(nmo),idifb(nmo))
c
c.....Open file to store alpha or beta configurations.
c
      indconfa=58
      indconfb=60
      open (indconfa,file=wfnfile(1:leng(wfnfile))//'confsa.dat',
     &      access='direct',recl=RecLength*nmo,form='unformatted')
      open (indconfb,file=wfnfile(1:leng(wfnfile))//'confsb.dat',
     &      access='direct',recl=RecLength*nmo,form='unformatted')
      xnorm=zero
c
      if (largwr) then
        write (stdout,*) '#'
        write (stdout,*) '# ALPHA AND BETA ACTIVE CONFIGURATIONS'
      endif
      do i=1,ndets
         mal=0
         mbe=0
         read (ifilc,rec=i) (cidet(m),m=0,nelact)
         xnorm=xnorm+cidet(0)*cidet(0)
         do m=1,nelact
           if (int(cidet(m)).gt.0) then
             mal=mal+1
             iorda(mal)=mal
             ordia(mal)=abs(int(cidet(m)))
           else
             mbe=mbe+1
             iordb(mbe)=mbe
             ordib(mbe)=abs(int(cidet(m)))
           endif
         enddo
c
c........Alpha electrons
c
         call iqcksort (ordia,iorda,nmo,1,mal)
         do j=1,mal
           nordia(j)=ordia(iorda(j))
         enddo
         call colord (ordia,nordia,mal,nmo,npera,nsiga(i),ndifa,idifa)
         if (i.eq.1) then
           npa=1
           kwa(i)=npa
           write (indconfa,rec=npa) (nordia(j),j=1,mal)
           if (largwr) then
             write (stdout,'(a,30I4)') ' # ALPHA ',(nordia(j),j=1,mal)
           endif
         else
           do k=1,npa
             inlist=.true.
             read (indconfa,rec=k) (ordja(j),j=1,mal)
             do j=1,mal
               inlist=inlist.and.(ordja(j).eq.nordia(j))
             enddo
             kwa(i)=k
             if (inlist) goto 67
           enddo
           npa=npa+1
           kwa(i)=npa
           write (indconfa,rec=npa) (nordia(j),j=1,mal)
           if (largwr) then
             write (stdout,'(a,30I4)') ' # ALPHA ',(nordia(j),j=1,mal)
           endif
         endif
c
c........Beta electrons
c
 67      call iqcksort (ordib,iordb,nmo,1,mbe)
         do j=1,mbe
           nordib(j)=ordib(iordb(j))
         enddo
         call colord (ordib,nordib,mbe,nmo,nperb,nsigb(i),ndifb,idifb)
         if (i.eq.1) then
           npb=1
           kwb(i)=npb
           write (indconfb,rec=npb) (nordib(j),j=1,mbe)
           if (largwr) then
             write (stdout,'(a,30I4)') ' # BETA  ',(nordib(j),j=1,mbe)
           endif
         else
           do k=1,npb
             inlist=.true.
             read (indconfb,rec=k) (ordjb(j),j=1,mbe)
             do j=1,mbe
               inlist=inlist.and.(ordjb(j).eq.nordib(j))
             enddo
             kwb(i)=k
             if (inlist) goto 68
           enddo
           npb=npb+1
           kwb(i)=npb
           write (indconfb,rec=npb) (nordib(j),j=1,mbe)
           if (largwr) then
             write (stdout,'(a,30I4)') ' # BETA  ',(nordib(j),j=1,mbe)
           endif
         endif
 68   enddo
      write (stdout,*) '# NUMBER OF alpha CONFIGS = ',npa
      write (stdout,*) '# NUMBER OF beta  CONFIGS = ',npb
      deallocate (ordia,ordja,ordib,ordjb,iorda,iordb,nordia,nordib)
c
      nprobts=nproba*nprobb
c
c.....Random numbers generation for alpha block.
c
      write (stdout,44)
      call semilla (idum)
      idum=-idum
      dumi=random(idum)
      allocate (am(nproba,nproba))
      allocate (tpowa(nproba,ngroup))
      allocate (ipvta(nproba),w(nproba))
 4000 am=zero
      tpowa=zero
      do i=1,nproba
        do k=1,ngroup-1
          tpowa(i,k)=two*random(idum)-one
        enddo
      enddo
      do i=1,nproba
        do j=1,nproba
          aco=one
          do k=1,ngroup-1
            aco=aco*tpowa(i,k)**resnca(j,k)
          enddo
          am(i,j)=aco
        enddo
      enddo
      call timer (2,ipid2,'_dgeco    ',-1)
      call dgeco (am,nproba,nproba,ipvta,rcond,w)
      call timer (4,ipid2,'_dgeco    ',-1)
      write (stdout,444) rcond,'ALPHA'
*     if (one + rcond .eq. one) goto 4000
      if (abs(rcond).lt.1d-50) goto 4000
      deallocate (w)
c
c.....run over all pairs of different alpha configs
c
      allocate (ordia(nmo),ordja(nmo))
      iper=0
      shundred=dble(npa)*(dble(npa)+1)/twohundred
      xhundred=shundred
      indet=0d0
      indproba=59
      open (indproba,file=wfnfile(1:leng(wfnfile))//'probala.dat',
     &     access='direct',recl=RecLength*2*nproba,form='unformatted')
      mal=ncore+mal
      allocate (proba(nproba))
      allocate (ovea(ngroup,moval,moval),oveaa(moval,moval))
      allocate (indxa(moval),ww(moval))
      do m=1,ncore
        ordia(m)=m
        ordja(m)=m
      enddo
      do i=1,npa
        read (indconfa,rec=i) (ordia(m+ncore),m=1,mal-ncore)
        do m=1,mal-ncore
          ordia(m+ncore)=ordia(m+ncore)+ncore
        enddo
        do j=1,i
          read (indconfa,rec=j) (ordja(m+ncore),m=1,mal-ncore)
          do m=1,mal-ncore
            ordja(m+ncore)=ordja(m+ncore)+ncore
          enddo
          indet=indet+1d0
c
c.........Testing how the most costly part of the run is going on.
c
          if (indet.ge.xhundred) then
            xhundred=xhundred+shundred
            iper=iper+1
            if (npa.gt.10) write (stderr,'(1x,a,i12,a)')
     &         '# Linear System Solution (ALPHA):', iper,' % done'
            call flush(stdout)
          endif
c
c.........Overlap matrices in ALPHA block.
c
          do m=1,malv
            ioma=ordia(ival(m))
            do k=1,malv
              ioka=ordja(ival(k))
              do igr=1,ngroup
                 ovea(igr,m,k)=sg(igr,ioma,ioka)
              enddo
            enddo
          enddo
c
          if (malv.gt.izero) then
            call timer (2,ipid3,'_determ   ',-1)
            do n=1,nproba
              do m=1,malv
                do k=1,malv
                  dumi=ovea(ngroup,m,k)
                  do igr=1,ngroup-1
                     dumi=dumi+tpowa(n,igr)*ovea(igr,m,k)
                  enddo
                  oveaa(m,k)=dumi
                enddo
              enddo
              call dgeco (oveaa,moval,malv,indxa,rcond,ww)
              job=10
              deter=zero
              call dgedi (oveaa,moval,malv,indxa,deter,ww,job)
              proba(n)=deter(1)*tenp**deter(2)
            enddo
            call timer (4,ipid3,'_determ   ',-1)
            job=0
            call timer (2,ipid4,'_dgesl    ',-1)
            call dgesl (am,nproba,nproba,ipvta,proba,job)
            call timer (4,ipid4,'_dgesl    ',-1)
          else
            nproba=1
            proba(1)=one
          endif
          irec=i*(i-1)/2+j
          write (indproba,rec=irec)(proba(n),n=1,nproba)
        enddo
      enddo
      deallocate (tpowa,am,ipvta)
      deallocate (ordia,ordja,indxa)
c
c.....Random numbers generation for beta block.
c
      write (stdout,44)
      allocate (bm(nprobb,nprobb))
      allocate (tpowb(nprobb,ngroup))
      allocate (ipvtb(nprobb),w(nprobb))
 5000 bm=zero
      tpowb=zero
      do i=1,nprobb
        do k=1,ngroup-1
          tpowb(i,k)=two*random(idum)-one
        enddo
      enddo
      do i=1,nprobb
        do j=1,nprobb
          aco=one
          do k=1,ngroup-1
            aco=aco*tpowb(i,k)**resncb(j,k)
          enddo
          bm(i,j)=aco
        enddo
      enddo
      call timer (2,ipid2,'_dgeco    ',-1)
      call dgeco (bm,nprobb,nprobb,ipvtb,rcond,w)
      call timer (4,ipid2,'_dgeco    ',-1)
      write (stdout,444) rcond,'BETA '
      if (one + rcond .eq. one) goto 5000
      deallocate (w)
c
c.....run over all pairs of different beta configs
c
      allocate (ordib(nmo),ordjb(nmo))
      allocate (indxb(moval))
      iper=0
      shundred=dble(npb)*(dble(npb)+1)/twohundred
      xhundred=shundred
      indet=0d0
      indprobb=61
      open (indprobb,file=wfnfile(1:leng(wfnfile))//'probalb.dat',
     &     access='direct',recl=RecLength*2*nprobb,form='unformatted')
      mbe=ncore+mbe
      allocate (probb(nprobb))
      do m=1,ncore
        ordib(m)=m
        ordjb(m)=m
      enddo
      do i=1,npb
        read (indconfb,rec=i) (ordib(m+ncore),m=1,mbe-ncore)
        do m=1,mbe-ncore
          ordib(m+ncore)=ordib(m+ncore)+ncore
        enddo
        do j=1,i
          read (indconfb,rec=j) (ordjb(m+ncore),m=1,mbe-ncore)
          do m=1,mbe-ncore
            ordjb(m+ncore)=ordjb(m+ncore)+ncore
          enddo
          indet=indet+1d0
c
c..........Testing how the most costly part of the run is going on.
c
           if (indet.ge.xhundred) then
             xhundred=xhundred+shundred
             iper=iper+1
             if (npb.gt.10) write (stderr,'(1x,a,i12,a)')
     &          '# Linear System Solution (BETA ):', iper,' % done'
             call flush(stdout)
           endif
c
c..........Overlap matrices in BETA  block.
c
           do m=1,mbev
             iomb=ordib(ival(m))
             do k=1,mbev
               iokb=ordjb(ival(k))
               do igr=1,ngroup
                  ovea(igr,m,k)=sg(igr,iomb,iokb)
               enddo
             enddo
           enddo
           if (mbev.gt.izero) then
             call timer (2,ipid6,'_determ   ',-1)
             do n=1,nprobb
               do m=1,mbev
                 do k=1,mbev
                   dumi=ovea(ngroup,m,k)
                   do igr=1,ngroup-1
                      dumi=dumi+tpowb(n,igr)*ovea(igr,m,k)
                   enddo
                   oveaa(m,k)=dumi
                 enddo
               enddo
               call dgeco (oveaa,moval,mbev,indxb,rcond,ww)
               job=10
               deter=zero
               call dgedi (oveaa,moval,mbev,indxb,deter,ww,job)
               probb(n)=deter(1)*tenp**deter(2)
             enddo
             call timer (4,ipid6,'_determ   ',-1)
             job=0
             call timer (2,ipid7,'_dgesl    ',-1)
             call dgesl (bm,nprobb,nprobb,ipvtb,probb,job)
             call timer (4,ipid7,'_dgesl    ',-1)
           else
             nprobb=1
             probb(1)=one
           endif
           irec=i*(i-1)/2+j
           write (indprobb,rec=irec)(probb(n),n=1,nprobb)
        enddo
      enddo
      deallocate (tpowb,bm,ipvtb)
      deallocate (ordib,ordjb)
c
      deallocate (ovea,oveaa,indxb,ww)
c
      iper=0
      shundred=ndets*(ndets+1)/twohundred
      xhundred=shundred
      indet=0d0
c
c.....run over all pairs of determinants.
c
      allocate (probt(nprobts))
      probt=zero
      do i=1,ndets
         read (ifilc  ,rec=i) cidet(0)
         cdprim=cidet(0)
         do j=1,i
           indet=indet+1d0
c
c..........Testing how the most costly part of the run is going on.
c
           if (indet.ge.xhundred) then
             xhundred=xhundred+shundred
             iper=iper+1
             if (ndets.gt.10) write (stderr,'(1x,a,i12,a)')
     &       '# Probability_ALPHA x Probability_BETA :', iper,' % done'
             call flush(stdout)
           endif
           read (ifilc,rec=j) cidet(0)
           cdsecond=cidet(0)
           cd=cdprim*cdsecond/xnorm
           cd=cd*nsiga(i)*nsiga(j)*nsigb(i)*nsigb(j)
           if (abs(cd).le.abs(epsdet)) goto 1000
           if (i.ne.j) cd=cd+cd
           imina=min(kwa(i),kwa(j))
           iminb=min(kwb(i),kwb(j))
           imaxa=max(kwa(i),kwa(j))
           imaxb=max(kwb(i),kwb(j))
           ireca=imaxa*(imaxa-1)/2+imina
           irecb=imaxb*(imaxb-1)/2+iminb
           read (indproba,rec=ireca)(proba(k),k=1,nproba)
           read (indprobb,rec=irecb)(probb(k),k=1,nprobb)
c
c..........Total spin-splitted probabilities for this pair of dets.
c
           ij=0
           do ia=1,nproba
             do ib=1,nprobb
               ij=ij+1
               probt(ij)=probt(ij)+cd*proba(ia)*probb(ib)
             enddo
           enddo
 1000    enddo
      enddo
c
c.....Write probabilities.
c
      write (stdout,200) ngroup,nprobts
      nprobg=0
      sum=zero
      sumtot=zero
      ij=0
      do i=1,nproba
        do j=1,nprobb
          ij=ij+1
          sumtot=sumtot+probt(ij)
          if (probt(ij).ge.probcut) then
            nprobg=nprobg+1
            sum=sum+probt(ij)
            write (stdout,100) probt(ij),
     &      (resnca(i,igr)+mocogrp(igr),
     &       resncb(j,igr)+mocogrp(igr),igr=1,ngroup)
          endif
        enddo
      enddo
      write (stdout,106) sum, nprobg, probcut, sumtot
c
c.....Average population of alpha and beta electrons in each group.
c     Store also in proba() and probb() arrays the alpha and beta 
c     probabilities.
c     
      proba=zero
      probb=zero
      ik=0
      do i=1,nproba
        do k=1,nprobb
          ik=ik+1
          proba(i)=proba(i)+probt(ik)
          probb(k)=probb(k)+probt(ik)
        enddo
      enddo
c
c.....Write EDF for ALPHA and BETA electrons independently
c
      write (stdout,203) 'ALPHA','BETA ',ngroup,'ALPHA',nproba
      nprobg=0
      sum=zero
      sumtot=zero
      do i=1,nproba
        sumtot=sumtot+proba(i)
        if (proba(i).ge.probcut) then
          nprobg=nprobg+1
          sum=sum+proba(i)
          write (stdout,100) proba(i),
     &      (resnca(i,igr)+mocogrp(igr),igr=1,ngroup)
        endif
      enddo
      deallocate (proba)
      write (stdout,106) sum, nprobg, probcut, sumtot
      write (stdout,203) 'BETA ','ALPHA',ngroup,'BETA ',nprobb
      nprobg=0
      sum=zero
      sumtot=zero
      do i=1,nprobb
        sumtot=sumtot+probb(i)
        if (probb(i).ge.probcut) then
          nprobg=nprobg+1
          sum=sum+probb(i)
          write (stdout,100) probb(i),
     &      (resncb(i,igr)+mocogrp(igr),igr=1,ngroup)
        endif
      enddo
      deallocate (probb)
      write (stdout,106) sum, nprobg, probcut, sumtot
c
c.....Compute spin-splitted localization and delocatization indices.
c
      allocate (p1a(ngroup),p1b(ngroup))
      allocate (xlis(ngroup,2))
      if (ngroup.gt.1) then
        ngpair=ngroup*(ngroup-1)/2
        allocate (p2aa(ngpair))
        allocate (p2bb(ngpair))
        allocate (p2ab(ngpair))
        allocate (p2ba(ngpair))
        allocate (d1aa(ngpair))
        allocate (d1bb(ngpair))
        allocate (d1ab(ngpair))
        allocate (d1ba(ngpair))
      endif
      p1a =zero
      p1b =zero
      p2aa=zero
      p2ab=zero
      p2ba=zero
      p2bb=zero
      d1aa=zero
      d1ab=zero
      d1ba=zero
      d1bb=zero
      xlis=zero
      np=0
      do ia=1,nproba
        do ib=1,nprobb
          np=np+1
          p=probt(np)
          do i=1,ngroup
            p1a(i)=p1a(i)+p*resnca(ia,i)
            p1b(i)=p1b(i)+p*resncb(ib,i)
            xlis(i,1)=xlis(i,1)+p*resnca(ia,i)*resnca(ia,i)
            xlis(i,2)=xlis(i,2)+p*resncb(ib,i)*resncb(ib,i)
          enddo
          if (ngroup.gt.1) then
            ipair=0
            do i=1,ngroup
              do j=1,i-1
                ipair=ipair+1
                p2aa(ipair)=p2aa(ipair)+resnca(ia,i)*resnca(ia,j)*p
                p2ab(ipair)=p2ab(ipair)+resnca(ia,i)*resncb(ib,j)*p
                p2ba(ipair)=p2ba(ipair)+resncb(ib,i)*resnca(ia,j)*p
                p2bb(ipair)=p2bb(ipair)+resncb(ib,i)*resncb(ib,j)*p
              enddo
            enddo
          endif
        enddo
      enddo
c
      write (stdout,119)
      do i=1,ngroup
        write (stdout,19) i,p1a(i),i,p1b(i)
      enddo
      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (stdout,29) 
     &       i,j,p2aa(ipair),i,j,p2ab(ipair),
     &       i,j,p2ba(ipair),i,j,p2bb(ipair)
          enddo
        enddo
      endif
c 
      np=0
      do ia=1,nproba
        do ib=1,nprobb
          np=np+1
          p=probt(np)
          if (ngroup.gt.1) then
            ipair=0
            do i=1,ngroup
              do j=1,i-1
                ipair=ipair+1
                nai=resnca(ia,i)
                naj=resnca(ia,j)
                nbi=resncb(ib,i)
                nbj=resncb(ib,j)
                d1aa(ipair)=d1aa(ipair)-2*p*(p1a(i)-nai)*(p1a(j)-naj)
                d1ab(ipair)=d1ab(ipair)-2*p*(p1a(i)-nai)*(p1b(j)-nbj)
                d1ba(ipair)=d1ba(ipair)-2*p*(p1b(i)-nbi)*(p1a(j)-naj)
                d1bb(ipair)=d1bb(ipair)-2*p*(p1b(i)-nbi)*(p1b(j)-nbj)
              enddo
            enddo
          endif
        enddo
      enddo
c
      do i=1,ngroup
        xlis(i,1)=-(xlis(i,1)-p1a(i)*(p1a(i)+one))
        xlis(i,2)=-(xlis(i,2)-p1b(i)*(p1b(i)+one))
        terma=zero
        termb=zero
        if (p1a(i).gt.epsq8) terma=hundred*xlis(i,1)/p1a(i)
        if (p1b(i).gt.epsq8) termb=hundred*xlis(i,2)/p1b(i)
        write (stdout,69) i,i,xlis(i,1),terma,i,i,xlis(i,2),termb
      enddo
      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (stdout,49) i,j,d1aa(ipair),i,j,d1ab(ipair),
     &                        i,j,d1ba(ipair),i,j,d1bb(ipair)
          enddo
        enddo
      endif
c
      if (ngroup.gt.1) then
        deallocate (p2aa)
        deallocate (p2bb)
        deallocate (p2ab)
        deallocate (p2ba)
        deallocate (d1aa)
        deallocate (d1bb)
        deallocate (d1ab)
        deallocate (d1ba)
      endif
      deallocate (xlis)
c
c.....Obtain spinless multiple-fragment electron population covariances.
c
      call ndelta (probt,p1a,p1b,resnca,resncb,
     &    ngroup,nproba,nprobb,nprobts,stdout)
c
      deallocate (p1a,p1b)
c
c.....Computes spinless probabilities from spin resolved probabilities.
c.....Obtain the number of spinless real space resonant structures.
c
      combi=one
      do i=0,nel-1
        combi=combi*dble(nel+ngroup-1-i)/dble(nel-i)
      enddo
      nprobt=int(combi+epsq10)
c
      allocate (pnew(nprobt))
      allocate (resnc(nprobt,ngroup))
      allocate (ind(ngroup))
      pnew=zero
      n=0
      do ia=1,nproba
        do ib=1,nprobb
          n=n+1
          if (n.eq.1) then
            np=1
            do i=1,ngroup
              resnc(np,i)=resnca(ia,i)+resncb(ib,i)
            enddo
            pnew(np)=probt(n)
          else
            do i=1,ngroup
              ind(i)=resnca(ia,i)+resncb(ib,i)
            enddo
            do j=1,np
              inlist=.true.
              do k=1,ngroup
                inlist=inlist.and.(ind(k).eq.resnc(j,k))
              enddo
              if (inlist) then
                pnew(j)=pnew(j)+probt(n)
                goto 1
              endif
            enddo
            np=np+1
            pnew(np)=probt(n)
            do i=1,ngroup
              resnc(np,i)=resnca(ia,i)+resncb(ib,i)
            enddo
          endif
 1      enddo
      enddo
      deallocate (probt,ind)
c
c.....Determine if probabilities are ordered or not.
c
      write (stdout,202) ngroup,np
      npnew=0
      sumtot=zero
      sum=zero
      do n=1,np
        sumtot=sumtot+pnew(n)
        if (pnew(n).ge.probcut) then
          npnew=npnew+1
          sum=sum+pnew(n)
        endif
      enddo
c
      if (npnew.le.nstack.and.orderp) then
c
c.......Order by increasing value the highest probabilities.
c
        allocate (ioprob (npnew))
        allocate (resncord (npnew,ngroup))
        allocate (probord(npnew))
        npronew=0
        do n=1,np
          if (pnew(n).ge.probcut) then
            npronew=npronew+1
            probord(npronew)=pnew(n) 
            ioprob(npronew)=npronew
            do igr=1,ngroup
              resncord(npronew,igr) = resnc(n,igr)+2*mocogrp(igr)
            enddo
          endif
        enddo
        call qqsort (probord,ioprob,1,npnew,nstack)
        do n=npnew,1,-1
          m=ioprob(n)
          write (stdout,100) probord(m),(resncord(m,igr),igr=1,ngroup)
        enddo
        deallocate (ioprob,resncord,probord)
      else
c
c.......Probabilities are not ordered.
c
        do n=1,np
          if (pnew(n).ge.probcut) write (stdout,100) pnew(n),
     &        (resnc(n,igr)+2*mocogrp(igr),igr=1,ngroup)
        enddo
      endif
      write (stdout,106) sum, npnew, probcut, sumtot
c
c.....Compute spinless localization and delocalization indices
c
      allocate (p1(ngroup),xli(ngroup))
      p1=zero
      xli=zero
      allocate (npop(ngroup))
      if (ngroup.gt.1) then
        allocate (p2(ngroup*(ngroup-1)/2))
        allocate (d1(ngroup*(ngroup-1)/2))
        p2=zero
        d1=zero
      endif
      if (ngroup.gt.2) then
        allocate (p3(ngroup*(ngroup-1)*(ngroup-2)/6))
        allocate (d2(ngroup*(ngroup-1)*(ngroup-2)/6))
        p3=zero
        d2=zero
      endif
c
      nprob=np
      do np=1,nprob
        p=pnew(np)
        do i=1,ngroup
          npop(i)=resnc(np,i)
        enddo
        do i=1,ngroup
          p1(i)=p1(i)+npop(i)*p
          xli(i)=xli(i)+p*npop(i)*npop(i)
        enddo
        if (ngroup.gt.1) then
          ipair=0
          do i=1,ngroup
            do j=1,i-1
              ipair=ipair+1
              p2(ipair)=p2(ipair)+npop(i)*npop(j)*p
            enddo
          enddo
        endif
        if (ngroup.gt.2) then
          iter=0
          do i=1,ngroup
            do j=1,i-1
              do k=1,j-1
              iter=iter+1
              p3(iter)=p3(iter)+npop(i)*npop(j)*npop(k)*p
              enddo
            enddo
          enddo
        endif
      enddo
c
      write (stdout,11)
      do i=1,ngroup
        write (stdout,15) i,p1(i)
      enddo
      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (stdout,2) i,j,p2(ipair)
          enddo
        enddo
      endif
c
      if (ngroup.gt.2) then
        iter=0
        do i=1,ngroup
          do j=1,i-1
            do k=1,j-1
            iter=iter+1
            write (stdout,3) i,j,k,p3(iter)
            enddo
          enddo
        enddo
      endif
c
      do np=1,nprob
        p=pnew(np)
        do i=1,ngroup
          npop(i)=resnc(np,i)
        enddo
        if (ngroup.gt.1) then
          ipair=0
          do i=1,ngroup
            do j=1,i-1
              ipair=ipair+1
              d1(ipair)=d1(ipair)-2*p*(p1(i)-npop(i))*(p1(j)-npop(j))
            enddo
          enddo
        endif
        if (ngroup.gt.2) then
          iter=0
          do i=1,ngroup
            do j=1,i-1
              do k=1,j-1
                iter=iter+1
                d2(iter)=d2(iter)-
     &          2*p*(p1(i)-npop(i))*(p1(j)-npop(j))*(p1(k)-npop(k))
              enddo
            enddo
          enddo
        endif
      enddo  
      deallocate (pnew)
c
      do i=1,ngroup
        xli(i)=-(xli(i)-p1(i)*(p1(i)+one))
        write (stdout,6) i,i,xli(i),hundred*xli(i)/p1(i)
      enddo
      write (stdout,33)
      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (stdout,4) i,j,d1(ipair)
          enddo
        enddo
      endif
      if (ngroup.gt.2) then
        iter=0
        do i=1,ngroup
          do j=1,i-1
            do k=1,j-1
            iter=iter+1
               write (stdout,5) i,j,k,d2(iter)
            enddo
          enddo
        enddo
      endif
      deallocate (p1,xli)
      deallocate (npop)
      if (ngroup.gt.1) deallocate (p2)
      if (ngroup.gt.2) deallocate (p3)
      if (ngroup.gt.1) deallocate (d1)
      if (ngroup.gt.2) deallocate (d2)
      deallocate (resnc)
c
      close (indconfa,status='delete')
      close (indproba,status='delete')
      close (indconfb,status='delete')
      close (indprobb,status='delete')
      call timer (4,ipid1,'_calcedfd ',-1)
      return
c
c.....Formats.
c
 200  format (/,' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION',
     & ' INCLUDING SPIN', /,1x, '# ',72('-'),/,
     & ' # NUMBER OF GROUPS               = ',I8,/,
     & ' # TOTAL NUMBER OF PROBABILITIES  = ',I8,/,
     & ' # Gi(a) Gi(b) ARE THE NUMBER OF ALPHA AND BETA ELECTRONS ',
     & ' IN GROUP i',/,1x,'#',72('-'),/,
     & ' #     Probability',11x,'G1(a) G1(b) G2(a) G2(b) G3(a) G3(b)')
 203  format (/,' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION ',
     & 'FOR ',a,' ELECTRONS',/,
     & ' # FOR EACH VALUE, A SUM OVER ALL ',a,
     & ' RESONANT STRUCTURES HAS BEEN DONE',/,' #',72('-'),/,
     &  ' # NUMBER OF GROUPS',33x,' = ',I8,/,
     &  ' # TOTAL NUMBER OF PROBABILITIES FOR ',a,' ELECTRONS = ',I8,/,
     &   1x,'#',72('-'),/,
     & ' #     Probability            n1    n2    n3 ...')
 44   format (/,1x,'# Linear System Solver.',/)
 444  format (1x,'# RECIPROCAL OF RCOND = ',E18.12,' FOR ',a,
     &        1x,'ELECTRONS')
 100  format (' # ',F22.16,1x,12I6)
 106  format (1x,'#',72('-'),/,' #',F22.16,2x,'<-- SUM,',I8,
     & ' PROBABILITIES > ',E16.10,/,' #',F22.16,2x,'<--- TOTAL SUM',/,
     & 1x,'#',72('-'))
 112  format (1x,'#',/,1x,'# multiple-group delocalization indices',/,
     &        1x,'#')
 113  format (1x,'# DELTA = ',1x,F16.10,5x,20(I3))
 202  format (/,' # M-BASINS SPINLESS ELECTRON DISTRIBUTION FUNCTION',
     & /,1x,'#',72('-'),/,
     &  ' # NUMBER OF GROUPS               = ',I8,/,
     &  ' # TOTAL NUMBER OF PROBABILITIES  = ',I8,/,
     &   1x,'#',72('-'),/,
     & ' #     Probability            n1    n2    n3 ...')
 11   format (//1x,'Average populations and localization indices')
 33   format (//1x,'Delocalization indices,',
     &             ' Eq. (28) J. Chem. Phys.  126, 094102 (2007)')
 15   format (1x,'# <n(',I3,')>               = ',F16.10)
 2    format (1x,'# <n(',I3,') n(',I3,')>        = ',F16.10)
 3    format (1x,'# <n(',I3,') n(',I3,') n(',I3,')> = ',F16.10)
 4    format (1x,'# delta_(',I3,I3,')         = ',F16.10)
 5    format (1x,'# delta_(',I3,I3,I3,')      = ',F16.10)
 6    format (1x,'# delta_(',I3,I3,')         = ',F16.10,
     &        2x,'% Localization = ',F8.4)
 119  format (//1x,'Average populations and delocalization indices')
 19   format (1x,'<n(',I3,')_alpha>                  = ',F16.10,/,
     &        1x,'<n(',I3,')_beta>                   = ',F16.10)
 29   format (/1x,'<n(',I3,')_alpha n(',I3,')_alpha>     = ',F16.10,/,
     &         1x,'<n(',I3,')_alpha n(',I3,')_beta>      = ',F16.10,/,
     &         1x,'<n(',I3,')_beta  n(',I3,')_alpha>     = ',F16.10,/,
     &         1x,'<n(',I3,')_beta  n(',I3,')_beta>      = ',F16.10)
 49   format (/1x,'delta_(',I3,I3,')_{alpha,alpha)    = ',F16.10,/,
     &         1x,'delta_(',I3,I3,')_{alpha,beta )    = ',F16.10,/
     &         1x,'delta_(',I3,I3,')_{beta ,alpha)    = ',F16.10,/
     &         1x,'delta_(',I3,I3,')_{beta ,beta )    = ',F16.10)
 69   format (/1x,'delta_(',I3,I3,')_alpha            = ',F16.10,
     &         2x,'% Localization = ',F8.4,/,
     &         1x,'delta_(',I3,I3,')_beta             = ',F16.10,
     &         2x,'% Localization = ',F8.4)
 1120 format (3(1x,'#',/),1x,'# ',80('+'),/,
     &        1x,'# EXACT CALCULATION OF PROBABILITIES',/,1x,'#')
      end
