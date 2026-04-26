c
c-----------------------------------------------------------------------
c
      subroutine cindgr 
     &  ( epsdet,probcut,nproba,nprobb,ngroup,ndets,nmo,ncore,
     &    nelact,nel,moval,ival,mocogrp,naval,nbval,ifilc,lw,
     &    wfnfile,sg,resnca,resncb,orderp,wrlw,pnew,nprobt)
c
      USE        space_for_cidet
      USE        space_for_conf
      include   'implicit.inc'
      include   'param.inc'
      include   'constants.inc'
      include   'lengrec.inc'
      include   'mline.inc'
c
      real(kind=8), allocatable,dimension (:,:,:) :: ovea
      real(kind=8), allocatable,dimension (:,:)   :: am,bm,oveaa
      real(kind=8), allocatable,dimension (:,:)   :: tpowa,tpowb
      real(kind=8), allocatable,dimension (:)     :: proba,probb,probt
      real(kind=8), allocatable,dimension (:)     :: w,ww
      real(kind=8), allocatable,dimension (:)     :: probord
      integer, allocatable,dimension (:,:)   :: resnc
      integer, allocatable,dimension (:)     :: ipvta,ipvtb,indxa,ind
      integer, allocatable,dimension (:)     :: indxb
      integer, allocatable,dimension (:)     :: iorda,iordb
      integer, allocatable,dimension (:)     :: ordia,ordib,ordja,ordjb
      integer, allocatable,dimension (:)     :: nordia,nordib
      integer, allocatable,dimension (:)     :: ioprob
      integer, allocatable,dimension (:,:)   :: resncord
c
      real(kind=8)    sg(ngroup,nmo,nmo)
      real(kind=8)    dumi,random,deter(2)
      real(kind=8)    pnew(nprobt)
      integer    resnca(nproba,ngroup),resncb(nprobb,ngroup)
      integer    lw
      integer    mocogrp(ngroup)
      integer*4  idum
      logical    inlist,orderp,wrlw
      integer    ival(moval)
*     parameter  (mline=200)
      character*(mline) wfnfile
c
      call timer (2,icindgr,'_cindgr   ',-1)
      nprobts=nproba*nprobb
c
c.....Random numbers generation for alpha block.
c
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
      call dgeco (am,nproba,nproba,ipvta,rcond,w)
      if (one + rcond .eq. one) goto 4000
      deallocate (w)
c
c.....run over all pairs of different alpha configs
c
      allocate (ordia(nmo),ordja(nmo))
      iproba=59
      open (iproba,file=wfnfile(1:leng(wfnfile))//'probala.dat',
     &     access='direct',recl=RecLength*2*nproba,form='unformatted')
      mal=ncore+nelal
      allocate (proba(nproba))
      allocate (ovea(ngroup,moval,moval),oveaa(moval,moval))
      allocate (indxa(moval),ww(moval))
      do m=1,ncore
        ordia(m)=m
        ordja(m)=m
      enddo
      do i=1,npa
        ordia(ncore+1:mal)=nconfa(i,1:nelal)+ncore
        do j=1,i
          ordja(ncore+1:mal)=nconfa(j,1:nelal)+ncore
          do m=1,naval
            ioma=ordia(ival(m))
            do k=1,naval
              ioka=ordja(ival(k))
              do igr=1,ngroup
                 ovea(igr,m,k)=sg(igr,ioma,ioka)
              enddo
            enddo
          enddo
c
          if (naval.gt.izero) then
            do n=1,nproba
              do m=1,naval
                do k=1,naval
                  dumi=ovea(ngroup,m,k)
                  do igr=1,ngroup-1
                     dumi=dumi+tpowa(n,igr)*ovea(igr,m,k)
                  enddo
                  oveaa(m,k)=dumi
                enddo
              enddo
              call dgeco (oveaa,moval,naval,indxa,rcond,ww)
              job=10
              deter=zero
              call dgedi (oveaa,moval,naval,indxa,deter,ww,job)
              proba(n)=deter(1)*tenp**deter(2)
            enddo
            job=0
            call dgesl (am,nproba,nproba,ipvta,proba,job)
          else
            nproba=1
            proba(1)=one
          endif
          irec=i*(i-1)/2+j
          write (iproba,rec=irec)(proba(n),n=1,nproba)
        enddo
      enddo
      deallocate (tpowa,am,ipvta)
      deallocate (ordia,ordja,indxa)
c
c.....Random numbers generation for beta block.
c
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
      call dgeco (bm,nprobb,nprobb,ipvtb,rcond,w)
      if (one + rcond .eq. one) goto 5000
      deallocate (w)
c
c.....run over all pairs of different beta configs
c
      allocate (ordib(nmo),ordjb(nmo))
      allocate (indxb(moval))
      iprobb=61
      open (iprobb,file=wfnfile(1:leng(wfnfile))//'probalb.dat',
     &     access='direct',recl=RecLength*2*nprobb,form='unformatted')
      mbe=ncore+nelbe
      allocate (probb(nprobb))
      do m=1,ncore
        ordib(m)=m
        ordjb(m)=m
      enddo
      do i=1,npb
        ordib(ncore+1:mbe)=nconfb(i,1:nelbe)+ncore
        do j=1,i
          ordjb(ncore+1:mbe)=nconfb(j,1:nelbe)+ncore
           do m=1,nbval
             iomb=ordib(ival(m))
             do k=1,nbval
               iokb=ordjb(ival(k))
               do igr=1,ngroup
                  ovea(igr,m,k)=sg(igr,iomb,iokb)
               enddo
             enddo
           enddo
           if (nbval.gt.izero) then
             do n=1,nprobb
               do m=1,nbval
                 do k=1,nbval
                   dumi=ovea(ngroup,m,k)
                   do igr=1,ngroup-1
                      dumi=dumi+tpowb(n,igr)*ovea(igr,m,k)
                   enddo
                   oveaa(m,k)=dumi
                 enddo
               enddo
               call dgeco (oveaa,moval,nbval,indxb,rcond,ww)
               job=10
               deter=zero
               call dgedi (oveaa,moval,nbval,indxb,deter,ww,job)
               probb(n)=deter(1)*tenp**deter(2)
             enddo
             job=0
             call dgesl (bm,nprobb,nprobb,ipvtb,probb,job)
           else
             nprobb=1
             probb(1)=one
           endif
           irec=i*(i-1)/2+j
           write (iprobb,rec=irec)(probb(n),n=1,nprobb)
        enddo
      enddo
      deallocate (tpowb,bm,ipvtb,ordib,ordjb,ovea,oveaa,indxb,ww)
c
c.....run over all pairs of determinants.
c
      allocate (probt(nprobts))
      probt=zero
      do i=1,ndets
         cdprim=cdet(i)
         do j=1,i
           cdsecond=cdet(j)
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
           read (iproba,rec=ireca)(proba(k),k=1,nproba)
           read (iprobb,rec=irecb)(probb(k),k=1,nprobb)
           ij=0
           do ia=1,nproba
             do ib=1,nprobb
               ij=ij+1
               probt(ij)=probt(ij)+cd*proba(ia)*probb(ib)
             enddo
           enddo
 1000    enddo
      enddo
*     deallocate (kwa,kwb,nsiga,nsigb,proba,probb)
      deallocate (proba,probb)
c
c.....Computes spinless probabilities.
c
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
      if (wrlw) then
        write (lw,202) ngroup,np
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
c.........Order by increasing value the highest probabilities.
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
            write (lw,100) probord(m),(resncord(m,igr),igr=1,ngroup)
          enddo
          deallocate (ioprob,resncord,probord)
        else
c
c.........Probabilities are not ordered.
c
          do n=1,np
            if (pnew(n).ge.probcut) write (lw,100) pnew(n),
     &          (resnc(n,igr)+2*mocogrp(igr),igr=1,ngroup)
          enddo
        endif
        write (lw,106) sum, npnew, probcut, sumtot
      endif
      deallocate (resnc)
      close (iproba,status='delete')
      close (iprobb,status='delete')
      call timer (4,icindgr,'_cindgr   ',-1)
      return
c
c.....Formats.
c
 202  format (/,' # M-BASINS SPINLESS ELECTRON DISTRIBUTION FUNCTION',
     & /,1x,'#',72('-'),/,
     &  ' # NUMBER OF GROUPS               = ',I8,/,
     &  ' # TOTAL NUMBER OF PROBABILITIES  = ',I8,/,
     &   1x,'#',72('-'),/,
     & ' #     Probability            n1    n2    n3 ...')
 100  format (' # ',F22.16,1x,12I6)
 106  format (1x,'#',72('-'),/,' #',F22.16,2x,'<-- SUM,',I8,
     & ' PROBABILITIES > ',E16.10,/,' #',F22.16,2x,'<--- TOTAL SUM',/,
     & 1x,'#',72('-'))
      end
