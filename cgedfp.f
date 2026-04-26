c
c-----------------------------------------------------------------------
c
      subroutine cgedfp 
     &   (line,wfnfile,probcut,ngroup,stdout,aom,nfugrp,ifugrp,orderp)
c
      include     'implicit.inc'
      include     'param.inc'
      include     'constants.inc'
      include     'wfn.inc'
      include     'corr.inc'
      include     'mline.inc'
*     parameter   (mline = 200)  

      real(kind=8), allocatable,dimension (:,:)   :: am,tpow,probs,oveaa
      real(kind=8), allocatable,dimension (:)     :: probsa,probsb
      real(kind=8), allocatable,dimension (:,:,:) :: ovea
      real(kind=8), allocatable,dimension (:)     :: ww,probt,pnew,w
      real(kind=8), allocatable,dimension (:)  :: p1,p1a,p1b,p2,d1,p3,d2
      real(kind=8), allocatable,dimension (:)  :: xli
      real(kind=8), allocatable,dimension (:)  :: probord
      real(kind=8), allocatable,dimension (:,:,:) :: sg
      integer, allocatable,dimension (:)     :: ipvt,indx,ind
      integer, allocatable,dimension (:)     :: irdela,icdela
      integer, allocatable,dimension (:)     :: irdelb,icdelb
      integer, allocatable,dimension (:)     :: irowcol,iord
      integer, allocatable,dimension (:,:)   :: resnc
      integer, allocatable,dimension (:,:)   :: resncord
      integer, allocatable,dimension (:)     :: npop
      integer, allocatable,dimension (:)     :: ioprob
      integer, allocatable,dimension (:,:)   :: eleca,elecb
      integer, allocatable,dimension (:,:)   :: resnca,resncb
      integer, allocatable,dimension (:)     :: ikeepar,ikeepac
      integer, allocatable,dimension (:)     :: ikeepbr,ikeepbc

      real(kind=8) aom(ncent,nmo,nmo)
      integer nfugrp(ngroup)
      integer ifugrp(ncent,ngroup)
      character*(mline) line2,fileout
      character*(*) line,wfnfile
      integer nmo,nlmoa,nlmob
      real(kind=8) dumi,random,deter(2)
      integer stdout,leng
      integer*4  idum
      logical    inlist,orderp,diaga,diagb,setint,cutedf,ok

      character*1    digs(10),ch1,ch2,ch3
      data digs / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
c
c-----------------------------------------------------------------------
c
      call getdate (line2)

      allocate (irdela(nmo))
      allocate (icdela(nmo))
      allocate (irdelb(nmo))
      allocate (icdelb(nmo))
      allocate (irowcol(nmo))
      allocate (iord(nmo))

      if (ngroup.eq.izero) then
        stop ' # cgedfp.f: !! No groups have been defined yet !!'
      endif

      cutedf=.false.
      lp=1
      ok=setint(nlmoa,line,lp)
      if (ok) then
        if (nlmoa.gt.nmo) return
        if (setint(nlmob,line,lp)) then
          if (nlmob.gt.nmo) return
          line=line(lp:)
          lp=1
          do j=1,nlmoa
            ok=setint(irdela(j),line,lp)
            if (irdela(j).le.0.or.irdela(j).gt.nmo) then
              stop ' # cgedfp.f: irdela(j).le.0.or.irdela(j).gt.nmo'
            endif
            if (j.gt.1) then
              do k=1,j-1
                if (irdela(k).eq.irdela(j)) then
                  stop ' # cgedfp.f: irdela(k).eq.irdela(j)'
                endif
              enddo
            endif
          enddo
          do j=1,nlmoa
            ok=setint(icdela(j),line,lp)
            if (icdela(j).le.0.or.icdela(j).gt.nmo) then
              stop ' # cgedfp.f: icdela(j).le.0.or.icdela(j).gt.nmo'
            endif
            if (j.gt.1) then
              do k=1,j-1
                if (icdela(k).eq.icdela(j)) then
                  stop ' # cgedfp.f: icdela(k).eq.icdela(j)'
                endif
              enddo
            endif
          enddo
          do j=1,nlmob
            ok=setint(irdelb(j),line,lp)
            if (irdelb(j).le.0.or.irdelb(j).gt.nmo) then
              stop ' # cgedfp.f: irdelb(j).le.0.or.irdelb(j).gt.nmo'
            endif
            if (j.gt.1) then
              do k=1,j-1
                if (irdelb(k).eq.irdelb(j)) then
                  stop ' # cgedfp.f: irdelb(k).eq.irdelb(j)'
                endif
              enddo
            endif
          enddo
          do j=1,nlmob
            ok=setint(icdelb(j),line,lp)
            if (icdelb(j).le.0.or.icdelb(j).gt.nmo) then
              stop ' # cgedfp.f: icdelb(j).le.0.or.icdelb(j).gt.nmo'
            endif
            if (j.gt.1) then
              do k=1,j-1
                if (icdelb(k).eq.icdelb(j)) then
                  stop ' # cgedfp.f: icdelb(k).eq.icdelb(j)'
                endif
              enddo
            endif
          enddo
          if (nlmoa.gt.0.or.nlmob.gt.0) cutedf=.true.
        endif
      endif
c
c-----Order the deleted rows and columns
c
      do i=1,nlmoa
        iord(i)=i
      enddo
      call iqcksort (irdela, iord, nlmoa, 1, nlmoa)
      do i=1,nlmoa
        irowcol(i)=irdela(iord(i))
      enddo
      irdela(1:nlmoa)=irowcol(1:nlmoa)
      do i=1,nlmoa
        iord(i)=i
      enddo
      call iqcksort (icdela, iord, nlmoa, 1, nlmoa)
      do i=1,nlmoa
        irowcol(i)=icdela(iord(i))
      enddo
      icdela(1:nlmoa)=irowcol(1:nlmoa)

      do i=1,nlmob
        iord(i)=i
      enddo
      call iqcksort (irdelb, iord, nlmob, 1, nlmob)
      do i=1,nlmob
        irowcol(i)=irdelb(iord(i))
      enddo
      irdelb(1:nlmob)=irowcol(1:nlmob)
      do i=1,nlmob
        iord(i)=i
      enddo
      call iqcksort (icdelb, iord, nlmob, 1, nlmob)
      do i=1,nlmob
        irowcol(i)=icdelb(iord(i))
      enddo
      icdelb(1:nlmob)=irowcol(1:nlmob)


      if (.not.cutedf) then
        deallocate (irdela,icdela,irdelb,icdelb,iord,irowcol)
        return
      endif
      fileout=wfnfile(1:leng(wfnfile))//'--CGEDF'
      if (nlmoa.gt.0) then
        fileout=fileout(1:leng(fileout))//'-rowsA'
        do j=1,nlmoa
          i1=mod(irdela(j),100)
          i2=irdela(j)-i1
          i3=mod(i1,10)
          i4 = i1-i3
          i2 = i2/100
          i4 = i4/10
          do i=1,10
            if (i2.eq.i-1) ch1=digs(i)
            if (i4.eq.i-1) ch2=digs(i)
            if (i3.eq.i-1) ch3=digs(i)
          enddo
          fileout=fileout(1:leng(fileout))//'-'//ch1//ch2//ch3
        enddo
        fileout=fileout(1:leng(fileout))//'-colsA'
        do j=1,nlmoa
          i1=mod(icdela(j),100)
          i2=icdela(j)-i1
          i3=mod(i1,10)
          i4 = i1-i3
          i2 = i2/100
          i4 = I4/10
          do i=1,10
            if (i2.eq.i-1) ch1=digs(i)
            if (i4.eq.i-1) ch2=digs(i)
            if (i3.eq.i-1) ch3=digs(i)
          enddo
          fileout=fileout(1:leng(fileout))//'-'//ch1//ch2//ch3
        enddo
      endif
      if (nlmob.gt.0) then
        fileout=fileout(1:leng(fileout))//'-rowsB'
        do j=1,nlmob
         i1=mod(irdelb(j),100)
         i2=irdelb(j)-i1
         i3=mod(i1,10)
         i4 = i1-i3
         i2 = i2/100
         i4 = I4/10
         do i=1,10
           if (i2.eq.i-1) ch1=digs(i)
           if (i4.eq.i-1) ch2=digs(i)
           if (i3.eq.i-1) ch3=digs(i)
         enddo
         fileout=fileout(1:leng(fileout))//'-'//ch1//ch2//ch3
        enddo
        fileout=fileout(1:leng(fileout))//'-colsB'
        do j=1,nlmob
          i1=mod(icdelb(j),100)
          i2=icdelb(j)-i1
          i3=mod(i1,10)
          i4 = i1-i3
          i2 = i2/100
          i4 = I4/10
          do i=1,10
            if (i2.eq.i-1) ch1=digs(i)
            if (i4.eq.i-1) ch2=digs(i)
            if (i3.eq.i-1) ch3=digs(i)
          enddo
          fileout=fileout(1:leng(fileout))//'-'//ch1//ch2//ch3
        enddo
      endif
      lu19=19
      open (file=fileout(1:leng(fileout)),unit=lu19)
      write (lu19,1000) line2(1:leng(line2))
      do i=1,ngroup
         write (lu19,500) i,nfugrp(i),(ifugrp(j,i),j=1,nfugrp(i))
      enddo
      write (lu19,*) 
      write (lu19,7) nlmoa
      write (lu19,9) 'ROWS:   ',(irdela(j),j=1,nlmoa)
      write (lu19,9) 'COLUMNS:',(icdela(j),j=1,nlmoa)
      write (lu19,8) nlmob
      write (lu19,9) 'ROWS:   ',(irdelb(j),j=1,nlmob)
      write (lu19,9) 'COLUMNS:',(icdelb(j),j=1,nlmob)

      allocate (ikeepar(nmo))
      allocate (ikeepac(nmo))
      allocate (ikeepbr(nmo))
      allocate (ikeepbc(nmo))

      allocate (sg(ngroup,nmo,nmo),stat=ier)
c
c.....Obtain group overlap integrals.
c
      sg=0D0
      if (ngroup.eq.1) then
        do m=1,nmo
          do j=1,nmo
            sg(1,m,j) = 0D0
          enddo
          sg(1,m,m)=1D0
        enddo
      else
        sg=0D0
        do k=1,ncent
          do ngr=1,ngroup
            do i=1,nfugrp(ngr)
              if (ifugrp(i,ngr).eq.k) then
                 do m=1,nmo
                   do j=1,nmo
                     sg(ngr,m,j)=sg(ngr,m,j)+aom(k,m,j)
                   enddo
                 enddo
              endif
            enddo
          enddo
        enddo
      endif
c
      k=0
      do i=1,nmo
        do j=1,nlmoa
          if (irdela(j).eq.i) goto 101
        enddo
        k=k+1
        ikeepar(k)=i
 101  enddo
      nmoar=k
      k=0
      do i=1,nmo
        do j=1,nlmoa
          if (icdela(j).eq.i) goto 201
        enddo
        k=k+1
        ikeepac(k)=i
 201  enddo
      nmoac=k
      if (nmoar.ne.nmoac) then
        stop ' # cgedfp.f: Error deleting alpha MOs'
      else
        nmoa=nmoar
        nlmoa=nmo-nmoa
        diaga=.true.
        do i=1,nlmoa
          if (irdela(i).ne.icdela(i)) diaga=.false.
        enddo
      endif

      k=0
      do i=1,nmo
        do j=1,nlmob
          if (irdelb(j).eq.i) goto 301
        enddo
        k=k+1
        ikeepbr(k)=i
 301  enddo
      nmobr=k
      k=0
      do i=1,nmo
        do j=1,nlmob
          if (icdelb(j).eq.i) goto 401
        enddo
        k=k+1
        ikeepbc(k)=i
 401  enddo
      nmobc=k
      if (nmobr.ne.nmobc) then
        stop ' # cgedfp.f: Error deleting beta MOs'
      else
        nmob=nmobr
        nlmob=nmo-nmob
        diagb=.true.
        do i=1,nlmob
          if (irdelb(i).ne.icdelb(i)) diagb=.false.
        enddo
      endif

      deallocate (irdela,icdela,irdelb,icdelb)
      allocate (eleca(2,ngroup))
      allocate (elecb(2,ngroup))

      malcut=nalpha-nlmoa
      mbecut=nbeta-nlmob

      eleca(1,1:ngroup)=izero
      elecb(1,1:ngroup)=izero
      eleca(2,1:ngroup)=malcut
      elecb(2,1:ngroup)=mbecut

      mmm=1
      nproba=1
      nprobb=1
      do i=1,ngroup-1
        nproba=nproba*(malcut+i)
        nprobb=nprobb*(mbecut+i)
        mmm=mmm*i
      enddo
      nproba=nproba/mmm
      nprobb=nprobb/mmm
      allocate (resnca(nproba,ngroup))
      allocate (resncb(nprobb,ngroup))

      write (lu19,975) nproba,nprobb
      call rnprobs (eleca,resnca,malcut,ngroup,nproba,lu19)
      call rnprobs (elecb,resncb,mbecut,ngroup,nprobb,lu19)

      nelv=nel-nlmoa-nlmob
      write (lu19,1130)
      if (ndets.ne.1) then
        stop ' # cgedfp.f: Only Single-Determinant-WFNs allowed'
      endif
c
      write (lu19,44)
      call semilla (idum)
      idum=-idum
      dumi=random(idum)
c
      npmax=max(nproba,nprobb)
c
c.....Construct the first member of the linear system.
c
      allocate (probsa(nproba))
      if (nmoa.gt.0) then
        allocate (am(nproba,nproba))
        allocate (tpow(nproba,ngroup))
        allocate (ipvt(nproba))
        allocate (ovea(ngroup,nmoa,nmoa))
        allocate (oveaa(nmoa,nmoa))
        allocate (ww(nmoa))
        allocate (w(nproba))
        allocate (indx(nmoa))
c
        do m=1,nmoa
          do k=1,nmoa
            do igr=1,ngroup
              ovea(igr,m,k)=sg(igr,ikeepar(m),ikeepac(k))
            enddo
          enddo
        enddo
 2000   am=zero
        tpow=zero
        do i=1,nproba
          do k=1,ngroup-1
            tpow(i,k)=two*random(idum)-one
          enddo
        enddo
        do i=1,nproba
          do j=1,nproba
            aco=one
            do k=1,ngroup-1
              nn=resnca(j,k)
              aco=aco*tpow(i,k)**nn
            enddo
            am(i,j)=aco
          enddo
        enddo
c
        call timer (2,ipid1,'_dgeco    ',-1)
        call dgeco (am,nproba,nproba,ipvt,rcond,w)
        write (lu19,444) rcond
        if (one + rcond .eq. one) goto 2000
        call timer (4,ipid1,'_dgeco    ',-1)
c
        call timer (2,ipid2,'_determ   ',-1)
c
        do n=1,nproba
          do m=1,nmoa
            do k=1,nmoa
              dumi=ovea(ngroup,m,k)
              do igr=1,ngroup-1
                 dumi=dumi+tpow(n,igr)*ovea(igr,m,k)
              enddo
              oveaa(m,k)=dumi
            enddo
          enddo
c
c.........Determinants calculation using Netlib DGEDI routine.
c
          call dgeco (oveaa,nmoa,nmoa,indx,rcond,ww)
          job=10
          call dgedi (oveaa,nmoa,nmoa,indx,deter,ww,job)
          probsa(n)=deter(1)*tenp**deter(2)
        enddo
c
        call timer (4,ipid2,'_determ   ',-1)
c
c.......Linear System Solver: Netlib DGECO routine.
c 
        call timer (2,ipid3,'_dgesl    ',-1)
        job=0
        call dgesl (am,nproba,nproba,ipvt,probsa,job)
        call timer (4,ipid3,'_dgesl    ',-1)
c
        deallocate (ovea,oveaa,ww,w,indx,am,tpow,ipvt)
      else
        probsa(1)=1d0
      endif
      write (lu19,200) 'ALPHA',ngroup,'ALPHA',nproba
      nprobg=0
      sum=zero
      sumtot=zero
      do n=1,nproba
        sumtot=sumtot+probsa(n)
        if (probsa(n).ge.probcut) then
          nprobg=nprobg+1
          sum=sum+probsa(n)
           write (lu19,100) probsa(n),(resnca(n,igr),igr=1,ngroup)
        endif
      enddo
      write (lu19,106) sum, nprobg, probcut, sumtot
c
c.....Construct the first member of the linear system.
c
      allocate (probsb(nprobb))
      if (nmob.gt.0) then
        allocate (am(nprobb,nprobb))
        allocate (tpow(nprobb,ngroup))
        allocate (ipvt(nprobb))
        allocate (ovea(ngroup,nmob,nmob))
        allocate (oveaa(nmob,nmob))
        allocate (ww(nmob))
        allocate (w(nprobb))
        allocate (indx(nmob))
c
c.....Overlap matrices
c
        do m=1,nmob
          do k=1,nmob
            do igr=1,ngroup
              ovea(igr,m,k)=sg(igr,ikeepbr(m),ikeepbc(k))
            enddo
          enddo
        enddo
        deallocate (sg)
c
 2001   am=zero
        tpow=zero
        do i=1,nprobb
          do k=1,ngroup-1
            tpow(i,k)=two*random(idum)-one
          enddo
        enddo
        do i=1,nprobb
          do j=1,nprobb
            aco=one
            do k=1,ngroup-1
              nn=resncb(j,k)
              aco=aco*tpow(i,k)**nn
            enddo
            am(i,j)=aco
          enddo
        enddo
c
        call timer (2,ipid1,'_dgeco    ',-1)
        call dgeco (am,nprobb,nprobb,ipvt,rcond,w)
        write (lu19,444) rcond
        if (one + rcond .eq. one) goto 2001
        call timer (4,ipid1,'_dgeco    ',-1)
c
        call timer (2,ipid2,'_determ   ',-1)
c
        do n=1,nprobb
          do m=1,nmob
            do k=1,nmob
              dumi=ovea(ngroup,m,k)
              do igr=1,ngroup-1
                 dumi=dumi+tpow(n,igr)*ovea(igr,m,k)
              enddo
              oveaa(m,k)=dumi
            enddo
          enddo
c
c.........Determinants calculation using Netlib DGEDI routine.
c
          call dgeco (oveaa,nmob,nmob,indx,rcond,ww)
          job=10
          call dgedi (oveaa,nmob,nmob,indx,deter,ww,job)
          probsb(n)=deter(1)*tenp**deter(2)
        enddo
        deallocate (ikeepar,ikeepac,ikeepbr,ikeepbc)
c
        call timer (4,ipid2,'_determ   ',-1)
c
c.......Linear System Solver: Netlib DGECO routine.
c 
        call timer (2,ipid3,'_dgesl    ',-1)
        job=0
        call dgesl (am,nprobb,nprobb,ipvt,probsb,job)
        call timer (4,ipid3,'_dgesl    ',-1)
        deallocate (ovea,oveaa,ww,w,indx,am,tpow,ipvt)
      else
        probsb(1)=1d0
      endif
      write (lu19,200) 'BETA ',ngroup,'BETA ',nprobb
      nprobg=0
      sum=zero
      sumtot=zero
      do n=1,nprobb
        sumtot=sumtot+probsb(n)
        if (probsb(n).ge.probcut) then
          nprobg=nprobg+1
          sum=sum+probsb(n)
           write (lu19,100) probsb(n),(resncb(n,igr),igr=1,ngroup)
        endif
      enddo
      write (lu19,106) sum, nprobg, probcut, sumtot
c
c.....Total spin-splitted probabilities 
c
      allocate (probt(nproba*nprobb))
      nprobg=0
      sum=zero
      sumtot=zero
      write (lu19,4001) ngroup,nproba*nprobb
      ij=0
      xmut=zero
      entab=zero
      do i=1,nproba
        do j=1,nprobb
          ij=ij+1
          probt(ij)=probsa(i)*probsb(j)
          pbtpp=probt(ij)/probsa(i)/probsb(j)
          if (pbtpp.gt.0d0) then
            xmut=xmut+probt(ij)*log(pbtpp)
          endif
          if (probt(ij).gt.0d0) then
            entab=entab-probt(ij)*log(probt(ij))
          endif
          sumtot=sumtot+probt(ij)
          if (probt(ij).ge.probcut) then
            nprobg=nprobg+1
            sum=sum+probt(ij)
            write (lu19,100) probt(ij),
     &      (resnca(i,igr),resncb(j,igr),igr=1,ngroup)
          endif
        enddo
      enddo
      write (lu19,106) sum, nprobg, probcut, sumtot
c
c.....Entropy of ALPHA and BETA distributions
c
      enta=zero
      do i=1,nproba
        if (probsa(i).gt.0d0) enta=enta-probsa(i)*log(probsa(i))
      enddo
      entb=zero
      do i=1,nprobb
        if (probsb(i).gt.0d0) entb=entb-probsb(i)*log(probsb(i))
      enddo
c
c.....Average population of alpha and beta electrons in each group.
c
      allocate (p1a(ngroup),p1b(ngroup))
      p1a=zero
      do i=1,nproba
        do j=1,ngroup
          p1a(j)=p1a(j)+resnca(i,j)*probsa(i)
        enddo
      enddo
      p1b=zero
      do i=1,nprobb
        do j=1,ngroup
          p1b(j)=p1b(j)+resncb(i,j)*probsb(i)
        enddo
      enddo
      if (diaga) then
        write (lu19,1110)
        do i=1,ngroup
          write (lu19,151) i,'ALPHA',p1a(i)
        enddo
      endif
      if (diagb) then
        write (lu19,1111)
        do i=1,ngroup
          write (lu19,151) i,'BETA ',p1b(i)
        enddo
      endif
c
c.....Obtain multiple-fragment electron population covariances.
c
      npmax=max(nproba,nprobb)
      allocate (probs(npmax,2))
      do i=1,nproba
        probs(i,1)=probsa(i)
      enddo
      do i=1,nprobb
        probs(i,2)=probsb(i)
      enddo
      if (diaga.and.diagb) then
        call sndelta (probs,p1a,p1b,resnca,resncb,ngroup,
     &       nproba,nprobb,npmax,lu19)
      endif
      deallocate (probsa,probsb,probs,p1a,p1b)
c
c.....Computes spinless probabilities from spin resolved probabilities.
c     Obtain the number of spinless real space resonant structures.
c
      combi=one
      do i=0,nelv-1
        combi=combi*dble(nelv+ngroup-1-i)/dble(nelv-i)
      enddo
      nprobt=int(combi+1q-10)
c
      allocate (pnew(nprobt))
      allocate (resnc(nprobt,ngroup))
      allocate (ind(ngroup))
c
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
      write (lu19,202) ngroup,np
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
      if ((npnew.le.nstack).and.orderp) then
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
              resncord(npronew,igr) = resnc(n,igr)
            enddo
          endif
        enddo
        call qqsort (probord,ioprob,1,npnew,nstack)
        do n=npnew,1,-1
          m=ioprob(n)
          write (lu19,100) probord(m),(resncord(m,igr),igr=1,ngroup)
        enddo
        deallocate (ioprob,resncord,probord)
      else
c
c.......Probabilities are not ordered.
c
        do n=1,np
          if (pnew(n).ge.probcut) write (lu19,100) pnew(n),
     &        (resnc(n,igr),igr=1,ngroup)
        enddo
      endif
      write (lu19,106) sum, npnew, probcut, sumtot
      xmut  = xmut/xlog2
      entab = entab/xlog2
      enta  = enta/xlog2
      entb  = entb/xlog2
      sarb  = enta - xmut
      sbra  = entb - xmut
      write (lu19,107) xmut,entab,enta,entb,sarb,sbra
c
      if (diaga.and.diagb) then
c
c.......Compute spinless localization and delocalization indices
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
        write (lu19,11)
        do i=1,ngroup
          write (lu19,15) i,p1(i)
        enddo
        if (ngroup.gt.1) then
          ipair=0
          do i=1,ngroup
            do j=1,i-1
              ipair=ipair+1
              write (lu19,2) i,j,p2(ipair)
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
              write (lu19,3) i,j,k,p3(iter)
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
     &            2*p*(p1(i)-npop(i))*(p1(j)-npop(j))*(p1(k)-npop(k))
                enddo
              enddo
            enddo
          endif
        enddo  
        deallocate (resnc,resnca,resncb,eleca,elecb)
        deallocate (pnew)
c
        do i=1,ngroup
          xli(i)=-(xli(i)-p1(i)*(p1(i)+one))
          if (abs(p1(i)).gt.zero) then
            write (lu19,6) i,i,xli(i),hundred*xli(i)/p1(i)
          else
            write (lu19,67) i,i,xli(i)
          endif
        enddo
        write (stdout,33)
        if (ngroup.gt.1) then
          ipair=0
          do i=1,ngroup
            do j=1,i-1
              ipair=ipair+1
              write (lu19,4) i,j,d1(ipair)
            enddo
          enddo
        endif
        if (ngroup.gt.2) then
          iter=0
          do i=1,ngroup
            do j=1,i-1
              do k=1,j-1
              iter=iter+1
                 write (lu19,5) i,j,k,d2(iter)
              enddo
            enddo
          enddo
        endif
        deallocate (p1,xli,npop)
        if (ngroup.gt.1) deallocate (p2,d1)
        if (ngroup.gt.2) deallocate (p3,d2)
      endif
c
      call getdate (line2)
      write (lu19,4000) line2(1:leng(line2))
      close (unit=lu19)
c
c.....Formats.
c
 500  format (' # FRAGMENT ',I2,' FORMED BY ',I3,' BASINS:',30(1x,I4))
 1000 format(' #',18(' '),
     & 'ELECTRON NUMBER DISTRIBUTION FUNCTIONS AND',/,
     & ' #',24(' '),'COARSE GRAINED DENSITY MATRICES',/,' #',25(' '),
     & '(c) Evelio Francisco 2012',/,' #',28(' '),
     & 'University of Oviedo',/,' #',72('+'),/,
     & ' # Calculation starts on ',a)
 7    format (1x,'# Deleting ',I2,
     &  ' ROWS & COLUMNS from the ALPHA Group Overlap Matrix')
 8    format (1x,'# Deleting ',I2,
     &  ' ROWS & COLUMNS from the BETA  Group Overlap Matrix')
 9    format (1x,'# ',a,2x,10I3)

 200  format (' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION ',
     & 'FOR ',a,' ELECTRONS',/,1x,'#',72('-'),/,
     &  ' # NUMBER OF GROUPS',33x,' = ',I8,/,
     &  ' # TOTAL NUMBER OF PROBABILITIES FOR ',a,' ELECTRONS = ',I8,/,
     &   1x,'#',72('-'),/,
     & ' #     Probability            n1    n2    n3 ...')
 4001 format (/,' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION',
     & ' INCLUDING SPIN', /,1x, '# ',72('-'),/,
     & ' # NUMBER OF GROUPS               = ',I8,/,
     & ' # TOTAL NUMBER OF PROBABILITIES  = ',I8,/,
     & ' # Gi(a) Gi(b) ARE THE NUMBER OF ALPHA AND BETA ELECTRONS ',
     & ' IN GROUP i',/,1x,'#',72('-'),/,
     & ' #     Probability',11x,'G1(a) G1(b) G2(a) G2(b) G3(a) G3(b)')
 202  format (/,' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION',
     & ' NOT INCLUDING SPIN',
     & /,1x,'#',72('-'),/,
     &  ' # NUMBER OF GROUPS               = ',I8,/,
     &  ' # TOTAL NUMBER OF PROBABILITIES  = ',I8,/,
     &   1x,'#',72('-'),/,
     & ' #     Probability            n1    n2    n3 ...')
 44   format (1x,'# Linear System Solver: Netlib DGECO routine')
 444  format (1x,'# RECIPROCAL OF RCOND VALUE IS ',E18.11)
 100  format (' # ',F22.15,1x,12I6)
 106  format (1x,'#',72('-'),/,' #',F22.15,2x,'<-- SUM,',I8,
     & ' PROBABILITIES > ',E16.9,/,' #',F22.15,2x,'<--- TOTAL SUM',/,
     & 1x,'#',72('-'))
 107  format (1x,'# Mutual ALPHA-BETA information, I(A,B) = ',F22.15,
     &  ' bits',/,1x,'# Joint ALPHA-BETA entropy,      S(A,B) = ',
     &   F22.15,' bits',/,1x,
     & '# Entropy of ALPHA electrons,      S(A) = ',F22.15,' bits',/
     &  ' # Entropy of BETA  electrons,      S(A) = ',F22.15,' bits',/,
     &  ' # Entropy(A|B) =           S(A)- I(A,B) = ',F22.15,' bits',/,
     &  ' # Entropy(B|A) =           S(B)- I(A,B) = ',F22.15,' bits')
 11   format (//1x,'Average populations and localization indices')
 33   format (//1x,'Delocalization indices,',
     &             ' Eq. (28) J. Chem. Phys.  126, 094102 (2007)')
 15   format (1x,'<n(',I3,')>               = ',F16.9)
 2    format (1x,'<n(',I3,') n(',I3,')>        = ',F16.9)
 3    format (1x,'<n(',I3,') n(',I3,') n(',I3,')> = ',F16.9)
 4    format (1x,'delta_(',I3,I3,')         = ',F16.9)
 5    format (1x,'delta_(',I3,I3,I3,')      = ',F16.9)
 6    format (1x,'delta_(',I3,I3,')         = ',F16.9,
     &        2x,'% Localization = ',F8.4)
 67   format (1x,'delta_(',I3,I3,')         = ',F16.9,
     &        2x,'% Localization = Undefined')
 1110 format (//1x,'Average ALPHA populations')
 1111 format (//1x,'Average  BETA populations')
 151  format (1x,'<n(',I3,')>_',a,' = ',F16.9)
 4000 format (' # Calculation ends on ',a)
 1130 format (1(1x,'#',/),1x,'#',72('+'),/,1x,
     & '# EXACT CALCULATION OF PROBABILITIES',/,1x,'#')
 642  format (1x,'# Isopycnic Localization of CORE MOs: ')
 975  format (' # Number of Resonance Structures for ALPHA e = ',I10,/,
     &        ' # Number of Resonance Structures for BETA  e = ',I10)
      end
