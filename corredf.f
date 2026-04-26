c
c-----------------------------------------------------------------------
c
      subroutine corredf (epsdet,nproba,nprobb,ngroup,ndets,nmo,
     &  ncore,nelact,nel,moval,ival,naval,nbval,ifilc,lw,wfnfile,sg)
c
      USE        space_for_cidet
      USE        space_for_conf
      USE        space_for_rsrs

      include   'implicit.inc'
      include   'param.inc'
      include   'constants.inc'
      include   'mline.inc'
c
      real(kind=8),  allocatable,dimension (:,:,:) :: ovea
      real(kind=8),  allocatable,dimension (:,:)   :: am,bm,oveaa
      real(kind=8),  allocatable,dimension (:,:)   :: tpowa,tpowb
      real(kind=8),  allocatable,dimension (:)     :: proba,probb,probt
      real(kind=8),  allocatable,dimension (:)     :: w,ww
      real(kind=8),  allocatable,dimension (:)     :: probord
      integer, allocatable,dimension (:)     :: ipvta,ipvtb,indxa,ind
      integer, allocatable,dimension (:)     :: indxb
      integer, allocatable,dimension (:)     :: iorda,iordb
      integer, allocatable,dimension (:)     :: ordia,ordib,ordja,ordjb
      integer, allocatable,dimension (:)     :: nordia,nordib
c
      real(kind=8)     sg(ngroup,nmo,nmo)
      real(kind=8)     dumi,random,deter(2)
      integer    lw
      integer*4  idum
      integer    ival(moval)
*     parameter  (mline=200)
      character*(mline) wfnfile
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
     &     access='direct',recl=4*nproba,form='unformatted')
      mal=ncore+nelal

      if (.not.allocated(proba)) then
        allocate (proba(nproba),stat=ier)
        if (ier.ne.0) stop 'corredf.f: Cannot allocate proba()'
      endif
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
      deallocate (tpowa,am,ipvta,ordia,ordja,indxa)
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
     &     access='direct',recl=4*nprobb,form='unformatted')
      mbe=ncore+nelbe

      if (.not.allocated(probb)) then
        allocate (probb(nprobb),stat=ier)
        if (ier.ne.0) stop 'corredf.f: Cannot allocate probb()'
      endif
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
      if (.not.allocated(probt)) then
        allocate (probt(nproba*nprobb),stat=ier)
        if (ier.ne.0) stop 'corredf.f: Cannot allocate probt()'
      endif
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
c
c.....Probabilities, electron populations, and delocalization indices
c
      pexact(1:nprobt)=zero
      n=0
      do i=1,nproba
        do j=1,nprobb
          n=n+1
          nn=abtot(n)
          pexact(nn)=pexact(nn)+probt(n)
        enddo
      enddo

      p1(1:ngroup)=zero
      diexact(1:ngroup,1:ngroup)=zero
      do i=1,nn
        pex=pexact(i)
        do k=1,ngroup
          p1(k)=p1(k)+resnc(i,k)*pex
        enddo
      enddo
      do i=1,nn
        pex=pexact(i)
        do k=2,ngroup
          addi1=(p1(k)-resnc(i,k))
          do l=1,k-1
            addi2=(p1(l)-resnc(i,l))
            diexact(k,l)=diexact(k,l)-two*pex*addi1*addi2
          enddo
        enddo
      enddo

      if (allocated(probt)) then
        deallocate (probt,stat=ier)
        if (ier.ne.0) stop 'corredf.f: Cannot deallocate probt()'
      endif
      if (allocated(proba)) then
        deallocate (proba,stat=ier)
        if (ier.ne.0) stop 'corredf.f: Cannot deallocate proba()'
      endif
      if (allocated(probb)) then
        deallocate (probb,stat=ier)
        if (ier.ne.0) stop 'corredf.f: Cannot deallocate probb()'
      endif

      close (iproba,status='delete')
      close (iprobb,status='delete')
      return
      end
