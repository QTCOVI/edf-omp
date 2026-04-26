c
c-----------------------------------------------------------------------
c
      subroutine wrsrs 
     &  (epsdet,ndets,nta,ntb,nmo,ncore,ngroup,lw,prob,n,na,nb)
      USE space_for_sgarea
      include     'implicit.inc'
      include     'fact.inc'
      include     'constants.inc'
      parameter (maxg=10)
      real(kind=8) epsdet,prob,tottimes

      integer, allocatable,dimension (:)    :: iacta,iactb,iactt
      integer, allocatable,dimension (:)    :: in1,in2,in3,in4,in5,in6
      integer, allocatable,dimension (:)    :: ou1,ou2,ou3,ou4,ou5,ou6
      integer, allocatable,dimension (:)    :: io1,io2,io3,io4,io5,io6
      integer, allocatable,dimension (:)    :: wc
      integer, allocatable,dimension (:,:)  :: wca,wcb


      integer ngroup,lw,n(ngroup),na(ngroup),nb(ngroup)
      integer da(maxg),db(maxg),dd(maxg)
      integer*8 ntimes(maxg),nttim(2)
      logical first(maxg)
 
      if (ndets.eq.1) then
        call sdwrsrs (nta,ntb,nmo,ngroup,lw,prob,n,na,nb)
        return
      endif
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

      nacta=0
      nactb=0
      do i=1,ngroup
        nb(i)=n(i)-na(i)
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
      do isp=1,2
        if (isp.eq.1) then
          nactt=nacta
          ntab=nta
          iactt(1:nactt)=iacta(1:nacta)
          dd(1:nactt)=da(1:nacta)
        else
          nactt=nactb
          ntab=ntb
          iactt(1:nactt)=iactb(1:nactb)
          dd(1:nactt)=db(1:nactb)
        endif
        nrest=0
        tottimes=one
        if (ntab.gt.maxfac) then
          write (0,*) 'NTAB,MAXFAC = ',ntab,maxfac
          stop 'prsrs.f: !!! NTAB > MAXFAC. Increase MAXFAC in fact.inc'
        endif
        do i=1,nactt-1
          num=ntab-nrest
          n1=iactt(i)
          ntimes(i)=nint(fact(num)/fact(n1)/fact(num-n1))
          tottimes=tottimes*ntimes(i)
          nrest=nrest+n1
        enddo

        if (tottimes.gt.huge(l)) then
          write (0,123) tottimes
          write (lw,123) tottimes
          stop 
        endif
 123    format (3(' #',/),
     &    ' # prsrs.f: I am ver sorry.',/,
     &    ' # prsrs.f: A very HUGE integer number should be needed'/,
     &    ' # prsrs.f: nttim() array should need ',F22.0,' elements'/
     &    ' # --- STOP --- STOP --- STOP --- STOP --- STOP --- STOP',/,
     &   3(' #',/))
        nttim(isp)=idnint(tottimes)
        if (isp.eq.1) then
          allocate (wca(nttim(isp),nmo))
        else
          allocate (wcb(nttim(isp),nmo))
        endif

        first(1:ngroup)=.true.
        if (nactt.eq.1) then
          ia1=iactt(1)
          forall (k=1:ia1) wc(k)=dd(1)
          ncom=1
          if (isp.eq.1) then
            wca(ncom,1:ntab)=wc(1:ntab)
          else
            wcb(ncom,1:ntab)=wc(1:ntab)
          endif
        elseif (nactt.eq.2) then
          first(1)=.true.
          ntt=ntab
          forall (i=1:ntab) io1(i)=i
          ia1=iactt(1)
          ncom=0
          do i1=1,ntimes(1)
            ncom=ncom+1
            call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
            forall (k=1:ia1) wc(in1(k))=dd(1)
            ia2=iactt(2)
            forall (k=1:ia2) wc(ou1(k))=dd(2)
            if (isp.eq.1) then
              wca(ncom,1:ntab)=wc(1:ntab)
            else
              wcb(ncom,1:ntab)=wc(1:ntab)
            endif
          enddo
        elseif (nactt.eq.3) then
          first(1)=.true.
          ntt=ntab
          forall (i=1:ntab) io1(i)=i
          ia1=iactt(1)
          ncom=0
          do i1=1,ntimes(1)
            call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
            forall (k=1:ia1) wc(in1(k))=dd(1)
            first(2)=.true.
            ntt1=ntab-ia1
            io2(1:ntt1)=ou1(1:ntt1)
            ia2=iactt(2)
            do i2=1,ntimes(2)
              ncom=ncom+1
              call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
              forall (k=1:ia2) wc(in2(k))=dd(2)
              ia3=iactt(3)
              forall (k=1:ia3) wc(ou2(k))=dd(3)
              if (isp.eq.1) then
                wca(ncom,1:ntab)=wc(1:ntab)
              else
                wcb(ncom,1:ntab)=wc(1:ntab)
              endif
            enddo
          enddo
        elseif (nactt.eq.4) then
          first(1)=.true.
          ntt=ntab
          forall (i=1:ntab) io1(i)=i
          ia1=iactt(1)
          ncom=0
          do i1=1,ntimes(1)
            call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
            forall (k=1:ia1) wc(in1(k))=dd(1)
            first(2)=.true.
            ntt1=ntab-ia1
            io2(1:ntt1)=ou1(1:ntt1)
            ia2=iactt(2)
            do i2=1,ntimes(2)
              call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
              forall (k=1:ia2) wc(in2(k))=dd(2)
              first(3)=.true.
              ntt2=ntab-ia1-ia2
              io3(1:ntt2)=ou2(1:ntt2)
              ia3=iactt(3)
              do i3=1,ntimes(3)
                ncom=ncom+1
                call gengen (3,ia3,ntt2,maxg,nel,io3,in3,ou3,first)
                forall (k=1:ia3) wc(in3(k))=dd(3)
                ia4=iactt(4)
                forall (k=1:ia4) wc(ou3(k))=dd(4)
                if (isp.eq.1) then
                  wca(ncom,1:ntab)=wc(1:ntab)
                else
                  wcb(ncom,1:ntab)=wc(1:ntab)
                endif
              enddo
            enddo
          enddo
        elseif (nactt.eq.5) then
          first(1)=.true.
          ntt=ntab
          forall (i=1:ntab) io1(i)=i
          ia1=iactt(1)
          ncom=0
          do i1=1,ntimes(1)
            call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
            forall (k=1:ia1) wc(in1(k))=dd(1)
            first(2)=.true.
            ntt1=ntab-ia1
            io2(1:ntt1)=ou1(1:ntt1)
            ia2=iactt(2)
            do i2=1,ntimes(2)
              call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
              forall (k=1:ia2) wc(in2(k))=dd(2)
              first(3)=.true.
              ntt2=ntab-ia1-ia2
              io3(1:ntt2)=ou2(1:ntt2)
              ia3=iactt(3)
              do i3=1,ntimes(3)
                call gengen (3,ia3,ntt2,maxg,nel,io3,in3,ou3,first)
                forall (k=1:ia3) wc(in3(k))=dd(3)
                first(4)=.true.
                ntt3=ntab-ia1-ia2-ia3
                io4(1:ntt3)=ou3(1:ntt3)
                ia4=iactt(4)
                do i4=1,ntimes(4)
                  ncom=ncom+1
                  call gengen(4,ia4,ntt3,maxg,nel,io4,in4,ou4,first)
                  forall (k=1:ia4) wc(in4(k))=dd(4)
                  ia5=iactt(5)
                  forall (k=1:ia5) wc(ou4(k))=dd(5)
                  if (isp.eq.1) then
                    wca(ncom,1:ntab)=wc(1:ntab)
                  else
                    wcb(ncom,1:ntab)=wc(1:ntab)
                  endif
                enddo
              enddo
            enddo
          enddo
        elseif (nactt.eq.6) then
          first(1)=.true.
          ntt=ntab
          forall (i=1:ntab) io1(i)=i
          ia1=iactt(1)
          ncom=0
          do i1=1,ntimes(1)
            call gengen (1,ia1,ntt,maxg,nel,io1,in1,ou1,first)
            forall (k=1:ia1) wc(in1(k))=dd(1)
            first(2)=.true.
            ntt1=ntab-ia1
            io2(1:ntt1)=ou1(1:ntt1)
            ia2=iactt(2)
            do i2=1,ntimes(2)
              call gengen (2,ia2,ntt1,maxg,nel,io2,in2,ou2,first)
              forall (k=1:ia2) wc(in2(k))=dd(2)
              first(3)=.true.
              ntt2=ntab-ia1-ia2
              io3(1:ntt2)=ou2(1:ntt2)
              ia3=iactt(3)
              do i3=1,ntimes(3)
                call gengen (3,ia3,ntt2,maxg,nel,io3,in3,ou3,first)
                forall (k=1:ia3) wc(in3(k))=dd(3)
                first(4)=.true.
                ntt3=ntab-ia1-ia2-ia3
                io4(1:ntt3)=ou3(1:ntt3)
                ia4=iactt(4)
                do i4=1,ntimes(4)
                  call gengen(4,ia4,ntt3,maxg,nel,io4,in4,ou4,first)
                  forall (k=1:ia4) wc(in4(k))=dd(4)
                  first(5)=.true.
                  ntt4=ntab-ia1-ia2-ia3-ia4
                  io5(1:ntt4)=ou4(1:ntt4)
                  ia5=iactt(5)
                  do i5=1,ntimes(5)
                    ncom=ncom+1
                    call gengen(5,ia5,ntt4,maxg,nel,io5,in5,ou5,first)
                    forall (k=1:ia5) wc(in5(k))=dd(5)
                    ia6=iactt(6)
                    forall (k=1:ia6) wc(ou5(k))=dd(6)
                    if (isp.eq.1) then
                      wca(ncom,1:ntab)=wc(1:ntab)
                    else
                      wcb(ncom,1:ntab)=wc(1:ntab)
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        else
          stop 'prsrs.f: Too many domains'
        endif
      enddo
      call sumdetp (epsdet,prob,nttim(1),nttim(2),
     &   nta,ntb,ncore,wca,wcb,ngroup,ndets,nmo,.false.)

      deallocate (iacta,iactb,iactt,wca,wcb,wc)
      deallocate (in1,ou1,io1)
      if (ngroup.gt.1) deallocate (in2,ou2,io2)
      if (ngroup.gt.2) deallocate (in3,ou3,io3)
      if (ngroup.gt.3) deallocate (in4,ou4,io4)
      if (ngroup.gt.4) deallocate (in5,ou5,io5)
      if (ngroup.gt.5) deallocate (in6,ou6,io6)
      return
      end
