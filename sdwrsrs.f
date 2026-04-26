c
c-----------------------------------------------------------------------
c
      subroutine sdwrsrs (nta,ntb,nmo,ngroup,lw,prob,noccup,na,nb)
      USE space_for_conf
      USE space_for_sgarea
      include     'implicit.inc'
      include     'fact.inc'
      include     'constants.inc'
      parameter (maxg=10)
      real*8     prob,proba,probb

      integer*4, allocatable,dimension (:)    :: iacta,iactb,iactt
      integer*4, allocatable,dimension (:)    :: in1,in2,in3,in4,in5,in6
      integer*4, allocatable,dimension (:)    :: ou1,ou2,ou3,ou4,ou5,ou6
      integer*4, allocatable,dimension (:)    :: io1,io2,io3,io4,io5,io6
      integer*4, allocatable,dimension (:)    :: wc

      integer*4 ngroup,lw,noccup(ngroup),na(ngroup),nb(ngroup)
      integer*4 da(maxg),db(maxg)
      integer*8 ntimesa(maxg)
      integer*8 ntimesb(maxg)
      integer*8 i1,i2,i3,i4,i5
      logical first(maxg)

      real   (kind=8) :: overa(nta,nta)
      integer(kind=4) :: indxa(nta)
 
      nel=nta+ntb
c
      allocate (iacta(nel))
      allocate (iactb(nel))
      allocate (iactt(nel))
      allocate (wc(nmo))
c
      allocate (in1(nel),ou1(nel),io1(nel))
      if (ngroup.gt.1) allocate (in2(nel),ou2(nel),io2(nel))
      if (ngroup.gt.2) allocate (in3(nel),ou3(nel),io3(nel))
      if (ngroup.gt.3) allocate (in4(nel),ou4(nel),io4(nel))
      if (ngroup.gt.4) allocate (in5(nel),ou5(nel),io5(nel))
      if (ngroup.gt.5) allocate (in6(nel),ou6(nel),io6(nel))
c
c     nta:   Number of ALPHA electrons
c     ntb:   Number of BETA  electrons
c     nacta: Number of groups with a non-zero number of ALPHA electrons
c     nactb: Number of groups with a non-zero number of BETA  electrons
c     na:    The ALPHA Real space resonance structure
c     nb:    The BETA  Real space resonance structure
c     noccup: The full spinless Real space resonance structure
c
      nacta=0
      nactb=0
      do i=1,ngroup
        nb(i)=noccup(i)-na(i)
        if (na(i).gt.0) then
          nacta=nacta+1
          iacta(nacta)=na(i)
          da(nacta)=i
        endif
        if (nb(i).gt.0) then
          nactb=nactb+1
          iactb(nactb)=nb(i)
          db(nactb)=i
        endif
      enddo

      if (nta.gt.maxfac) then
        write (0,*) 'NTA,MAXFAC = ',nta,maxfac
        stop 'sdwrsrs.f: NTA > MAXFAC. Increase MAXFAC in fact.inc'
      endif
      nrest=0
      do i=1,nacta-1
        num=nta-nrest
        n1=iacta(i)
        ntimesa(i)=nint(fact(num)/fact(n1)/fact(num-n1))
        nrest=nrest+n1
      enddo

      if (ntb.gt.maxfac) then
        write (0,*) 'NTB,MAXFAC = ',ntb,maxfac
        stop 'sdwrsrs.f: NTB > MAXFAC. Increase MAXFAC in fact.inc'
      endif
      nrest=0
      do i=1,nactb-1
        num=ntb-nrest
        n1=iactb(i)
        ntimesb(i)=nint(fact(num)/fact(n1)/fact(num-n1))
        nrest=nrest+n1
      enddo

      if (nacta.ge.7) stop 'sdwrsrs.f: Too many domains'
      if (nactb.ge.7) stop 'sdwrsrs.f: Too many domains'
      first(1:ngroup)=.true.
      if (nacta.eq.1) then
        proba=0.0d+00
        ia1=iacta(1)
        forall (k=1:ia1) wc(k)=da(1)
        call sdwdeta (probax,nta,wc,nmo)
        proba=proba+probax
      elseif (nacta.eq.2) then
        proba=0.0d+00
        first(1)=.true.
        ntt=nta
        forall (i=1:nta) io1(i)=i
        ia1=iacta(1)
        do i1=1,ntimesa(1)
          call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
          forall (k=1:ia1) wc(in1(k))=da(1)
          ia2=iacta(2)
          forall (k=1:ia2) wc(ou1(k))=da(2)
          call sdwdeta (probax,nta,wc,nmo)
          proba=proba+probax
        enddo
      elseif (nacta.eq.3) then
        proba=0.0d+00
        first(1)=.true.
        ntt=nta
        forall (i=1:nta) io1(i)=i
        ia1=iacta(1)
        do i1=1,ntimesa(1)
          call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
          forall (k=1:ia1) wc(in1(k))=da(1)
          first(2)=.true.
          ntt1=nta-ia1
          io2(1:ntt1)=ou1(1:ntt1)
          ia2=iacta(2)
          do i2=1,ntimesa(2)
            call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
            forall (k=1:ia2) wc(in2(k))=da(2)
            ia3=iacta(3)
            forall (k=1:ia3) wc(ou2(k))=da(3)
            call sdwdeta (probax,nta,wc,nmo)
            proba=proba+probax
          enddo
        enddo
      elseif (nacta.eq.4) then
        proba=0.0d+00
        first(1)=.true.
        ntt=nta
        forall (i=1:nta) io1(i)=i
        ia1=iacta(1)
        do i1=1,ntimesa(1)
          call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
          forall (k=1:ia1) wc(in1(k))=da(1)
          first(2)=.true.
          ntt1=nta-ia1
          io2(1:ntt1)=ou1(1:ntt1)
          ia2=iacta(2)
          do i2=1,ntimesa(2)
            call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
            forall (k=1:ia2) wc(in2(k))=da(2)
            first(3)=.true.
            ntt2=nta-ia1-ia2
            io3(1:ntt2)=ou2(1:ntt2)
            ia3=iacta(3)
            do i3=1,ntimesa(3)
              call gengen (3,ia3,ntt2,maxg,nel,io3,in3,ou3,first)
              forall (k=1:ia3) wc(in3(k))=da(3)
              ia4=iacta(4)
              forall (k=1:ia4) wc(ou3(k))=da(4)
              call sdwdeta (probax,nta,wc,nmo)
              proba=proba+probax
            enddo
          enddo
        enddo
      elseif (nacta.eq.5) then
        proba=0.0d+00
        first(1)=.true.
        ntt=nta
        forall (i=1:nta) io1(i)=i
        ia1=iacta(1)
        do i1=1,ntimesa(1)
          call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
          forall (k=1:ia1) wc(in1(k))=da(1)
          first(2)=.true.
          ntt1=nta-ia1
          io2(1:ntt1)=ou1(1:ntt1)
          ia2=iacta(2)
          do i2=1,ntimesa(2)
            call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
            forall (k=1:ia2) wc(in2(k))=da(2)
            first(3)=.true.
            ntt2=nta-ia1-ia2
            io3(1:ntt2)=ou2(1:ntt2)
            ia3=iacta(3)
            do i3=1,ntimesa(3)
              call gengen (3,ia3,ntt2,maxg,nel,io3,in3,ou3,first)
              forall (k=1:ia3) wc(in3(k))=da(3)
              first(4)=.true.
              ntt3=nta-ia1-ia2-ia3
              io4(1:ntt3)=ou3(1:ntt3)
              ia4=iacta(4)
              do i4=1,ntimesa(4)
                call gengen(4,ia4,ntt3,maxg,nel,io4,in4,ou4,first)
                forall (k=1:ia4) wc(in4(k))=da(4)
                ia5=iacta(5)
                forall (k=1:ia5) wc(ou4(k))=da(5)
                call sdwdeta (probax,nta,wc,nmo)
                proba=proba+probax
              enddo
            enddo
          enddo
        enddo
      elseif (nacta.eq.6) then
        proba=0.0d+00
        first(1)=.true.
        ntt=nta
        forall (i=1:nta) io1(i)=i
        ia1=iacta(1)
        do i1=1,ntimesa(1)
          call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
          forall (k=1:ia1) wc(in1(k))=da(1)
          first(2)=.true.
          ntt1=nta-ia1
          io2(1:ntt1)=ou1(1:ntt1)
          ia2=iacta(2)
          do i2=1,ntimesa(2)
            call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
            forall (k=1:ia2) wc(in2(k))=da(2)
            first(3)=.true.
            ntt2=nta-ia1-ia2
            io3(1:ntt2)=ou2(1:ntt2)
            ia3=iacta(3)
            do i3=1,ntimesa(3)
              call gengen (3,ia3,ntt2,maxg,nel,io3,in3,ou3,first)
              forall (k=1:ia3) wc(in3(k))=da(3)
              first(4)=.true.
              ntt3=nta-ia1-ia2-ia3
              io4(1:ntt3)=ou3(1:ntt3)
              ia4=iacta(4)
              do i4=1,ntimesa(4)
                call gengen(4,ia4,ntt3,maxg,nel,io4,in4,ou4,first)
                forall (k=1:ia4) wc(in4(k))=da(4)
                first(5)=.true.
                ntt4=nta-ia1-ia2-ia3-ia4
                io5(1:ntt4)=ou4(1:ntt4)
                ia5=iacta(5)
                do i5=1,ntimesa(5)
                  call gengen(5,ia5,ntt4,maxg,nel,io5,in5,ou5,first)
                  forall (k=1:ia5) wc(in5(k))=da(5)
                  ia6=iacta(6)
                  forall (k=1:ia6) wc(ou5(k))=da(6)
                  call sdwdeta (probax,nta,wc,nmo)
                  proba=proba+probax
                enddo
              enddo
            enddo
          enddo
        enddo
      endif

      if (nactb.eq.1) then
        probb=0.0d+00
        ia1=iactb(1)
        forall (k=1:ia1) wc(k)=db(1)
        call sdwdetb (probbx,ntb,wc,nmo)
        probb=probb+probbx
      elseif (nactb.eq.2) then
        probb=0.0d+00
        first(1)=.true.
        ntt=ntb
        forall (i=1:ntb) io1(i)=i
        ia1=iactb(1)
        do i1=1,ntimesb(1)
          call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
          forall (k=1:ia1) wc(in1(k))=db(1)
          ia2=iactb(2)
          forall (k=1:ia2) wc(ou1(k))=db(2)
          call sdwdetb (probbx,ntb,wc,nmo)
          probb=probb+probbx
        enddo
      elseif (nactb.eq.3) then
        probb=0.0d+00
        first(1)=.true.
        ntt=ntb
        forall (i=1:ntb) io1(i)=i
        ia1=iactb(1)
        do i1=1,ntimesb(1)
          call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
          forall (k=1:ia1) wc(in1(k))=db(1)
          first(2)=.true.
          ntt1=ntb-ia1
          io2(1:ntt1)=ou1(1:ntt1)
          ia2=iactb(2)
          do i2=1,ntimesb(2)
            call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
            forall (k=1:ia2) wc(in2(k))=db(2)
            ia3=iactb(3)
            forall (k=1:ia3) wc(ou2(k))=db(3)
            call sdwdetb (probbx,ntb,wc,nmo)
            probb=probb+probbx
          enddo
        enddo
      elseif (nactb.eq.4) then
        probb=0.0d+00
        first(1)=.true.
        ntt=ntb
        forall (i=1:ntb) io1(i)=i
        ia1=iactb(1)
        do i1=1,ntimesb(1)
          call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
          forall (k=1:ia1) wc(in1(k))=db(1)
          first(2)=.true.
          ntt1=ntb-ia1
          io2(1:ntt1)=ou1(1:ntt1)
          ia2=iactb(2)
          do i2=1,ntimesb(2)
            call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
            forall (k=1:ia2) wc(in2(k))=db(2)
            first(3)=.true.
            ntt2=ntb-ia1-ia2
            io3(1:ntt2)=ou2(1:ntt2)
            ia3=iactb(3)
            do i3=1,ntimesb(3)
              call gengen (3,ia3,ntt2,maxg,nel,io3,in3,ou3,first)
              forall (k=1:ia3) wc(in3(k))=db(3)
              ia4=iactb(4)
              forall (k=1:ia4) wc(ou3(k))=db(4)
              call sdwdetb (probbx,ntb,wc,nmo)
              probb=probb+probbx
            enddo
          enddo
        enddo
      elseif (nactb.eq.5) then
        probb=0.0d+00
        first(1)=.true.
        ntt=ntb
        forall (i=1:ntb) io1(i)=i
        ia1=iactb(1)
        do i1=1,ntimesb(1)
          call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
          forall (k=1:ia1) wc(in1(k))=db(1)
          first(2)=.true.
          ntt1=ntb-ia1
          io2(1:ntt1)=ou1(1:ntt1)
          ia2=iactb(2)
          do i2=1,ntimesb(2)
            call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
            forall (k=1:ia2) wc(in2(k))=db(2)
            first(3)=.true.
            ntt2=ntb-ia1-ia2
            io3(1:ntt2)=ou2(1:ntt2)
            ia3=iactb(3)
            do i3=1,ntimesb(3)
              call gengen (3,ia3,ntt2,maxg,nel,io3,in3,ou3,first)
              forall (k=1:ia3) wc(in3(k))=db(3)
              first(4)=.true.
              ntt3=ntb-ia1-ia2-ia3
              io4(1:ntt3)=ou3(1:ntt3)
              ia4=iactb(4)
              do i4=1,ntimesb(4)
                call gengen(4,ia4,ntt3,maxg,nel,io4,in4,ou4,first)
                forall (k=1:ia4) wc(in4(k))=db(4)
                ia5=iactb(5)
                forall (k=1:ia5) wc(ou4(k))=db(5)
                probb=probb+probbx
              enddo
            enddo
          enddo
        enddo
      elseif (nactb.eq.6) then
        probb=0.0d+00
        first(1)=.true.
        ntt=ntb
        forall (i=1:ntb) io1(i)=i
        ia1=iactb(1)
        do i1=1,ntimesb(1)
          call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
          forall (k=1:ia1) wc(in1(k))=db(1)
          first(2)=.true.
          ntt1=ntb-ia1
          io2(1:ntt1)=ou1(1:ntt1)
          ia2=iactb(2)
          do i2=1,ntimesb(2)
            call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
            forall (k=1:ia2) wc(in2(k))=db(2)
            first(3)=.true.
            ntt2=ntb-ia1-ia2
            io3(1:ntt2)=ou2(1:ntt2)
            ia3=iactb(3)
            do i3=1,ntimesb(3)
              call gengen (3,ia3,ntt2,maxg,nel,io3,in3,ou3,first)
              forall (k=1:ia3) wc(in3(k))=db(3)
              first(4)=.true.
              ntt3=ntb-ia1-ia2-ia3
              io4(1:ntt3)=ou3(1:ntt3)
              ia4=iactb(4)
              do i4=1,ntimesb(4)
                call gengen(4,ia4,ntt3,maxg,nel,io4,in4,ou4,first)
                forall (k=1:ia4) wc(in4(k))=db(4)
                first(5)=.true.
                ntt4=ntb-ia1-ia2-ia3-ia4
                io5(1:ntt4)=ou4(1:ntt4)
                ia5=iactb(5)
                do i5=1,ntimesb(5)
                  call gengen(5,ia5,ntt4,maxg,nel,io5,in5,ou5,first)
                  forall (k=1:ia5) wc(in5(k))=db(5)
                  ia6=iactb(6)
                  forall (k=1:ia6) wc(ou5(k))=db(6)
                  call sdwdetb (probbx,ntb,wc,nmo)
                  probb=probb+probbx
                enddo
              enddo
            enddo
          enddo
        enddo
      endif
      prob = proba * probb

      write (lw,6,advance='no') na(1:ngroup)
      write (lw,7,advance='no') nb(1:ngroup)
      write (lw,8) prob,proba,probb
 6    format (' # (alpha',20I3)
 7    format (' --- beta',20I3)
 8    format (' ) PROB = ',E17.10,' = ',E17.10,' x ',E17.10)
      deallocate (iacta,iactb,iactt,wc)
      deallocate (in1,ou1,io1)
      if (ngroup.gt.1) deallocate (in2,ou2,io2)
      if (ngroup.gt.2) deallocate (in3,ou3,io3)
      if (ngroup.gt.3) deallocate (in4,ou4,io4)
      if (ngroup.gt.4) deallocate (in5,ou5,io5)
      if (ngroup.gt.5) deallocate (in6,ou6,io6)
      return
      end
