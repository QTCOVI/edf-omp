c
c-----------------------------------------------------------------------
c
      subroutine rhocond 
     &  (line,wfnfile,iwfn,ngroup,stdout,aom,nfugrp,ifugrp)
      USE        space_for_wfnbasis
      USE        space_for_wfncoef
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
      real(kind=8), allocatable,dimension (:,:,:) :: sg,sgnat
      real(kind=8), allocatable,dimension (:,:)   :: rhoconda,rhocondb
      real(kind=8), allocatable,dimension (:,:)   :: rhotot,rhovec
      real(kind=8), allocatable,dimension (:)     :: rhoeig
      real(kind=8), allocatable,dimension (:,:)   :: newcoef
      integer, allocatable,dimension (:)     :: ipvt,indx,ind
      integer, allocatable,dimension (:,:)   :: resnc
      integer, allocatable,dimension (:,:)   :: resncord
      integer, allocatable,dimension (:)     :: npop
      integer, allocatable,dimension (:)     :: condens
      integer, allocatable,dimension (:,:)   :: eleca,elecb
      integer, allocatable,dimension (:,:)   :: resnca,resncb
      integer, allocatable,dimension (:)     :: ikeepar,ikeepac
      integer, allocatable,dimension (:)     :: ikeepbr,ikeepbc

      real(kind=8) aom(ncent,nmo,nmo)
      integer nfugrp(ngroup)
      integer ifugrp(ncent,ngroup)
      character*(mline) line2,fileout,filewfn,wfncut
      character*(*) line,wfnfile
      character*80 wfnttl
      integer nmo
      real(kind=8) random,deter(2)
      integer stdout,leng
      integer*4  idum
      logical    inlist,setint,isthiscond,ok
      character (len=4) mode,fourchar

      character*1    digs(10),ch1,ch2,ch3
      data digs / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
c
c-----------------------------------------------------------------------
c
      call getdate (line2)

      if (ngroup.eq.izero) then
        stop '# rhocond.f: !! No groups defined yet !!'
      endif
      if (ndets.ne.1) then
        stop '# rhocond.f: Only Single-Det WFNs allowed'
      endif
      lw=20
c
c-----Read the condensation [C] to be computed.
c
      if (.not.allocated(condens)) allocate (condens(ngroup))
      lp=1
      nelread=0
      do i=1,ngroup
        ok=setint(condens(i),line,lp)
        if (ok) then
          condens(i)=abs(condens(i))
          if (condens(i).gt.nel) then
            write (0,*) '# rhocond.f: Too many e- in fragment ',i
            stop
          else
            nelread=nelread+condens(i)
          endif
        else
          write (0,*) '# rhocond.f: Bad format in RHOCOND keyword'
          stop
        endif
      enddo
c
c-----NELREAD must be equal to the number of electrons minus one.
c
      if (nelread.ne.nel-1) then
        write (0,*) '# rhocond.f: nelread = ',nelread
        write (0,*) '# rhocond.f: nel-1   = ',nel-1
        write (0,*) '# rhocond.f: Improper Number of electrons in [C]'
        stop
      endif
c
c.....Determine the name of the output file and WFN file for [C]
c
      filewfn=wfnfile(1:leng(wfnfile))
      call topoint (wfnfile,fileout)
      fileout=trim(fileout)//'.rhocond'
      do i=1,ngroup
        fileout=trim(fileout)//'-'//fourchar(condens(i))
        filewfn=trim(filewfn)//'-'//fourchar(condens(i))
      enddo
      open (file=fileout(1:leng(fileout)),unit=lw)
      write (lw,1000) line2(1:leng(line2))

      if (.not.allocated(ikeepar)) allocate (ikeepar(nmo))
      if (.not.allocated(ikeepac)) allocate (ikeepac(nmo))
      if (.not.allocated(ikeepbr)) allocate (ikeepbr(nmo))
      if (.not.allocated(ikeepbc)) allocate (ikeepbc(nmo))
      if (.not.allocated(sg)     ) allocate (sg(ngroup,nmo,nmo))
c
c.....Obtain group overlap integrals.
c
      sg=0D0
      if (ngroup.eq.1) then
        sg=0d0
        forall (m=1:nmo) sg(1,m,m)=1D0
      else
        sg=0D0
        do m=1,nmo
          do j=1,nmo
            do ngr=1,ngroup
              do i=1,nfugrp(ngr)
                k=ifugrp(i,ngr)
                sg(ngr,m,j)=sg(ngr,m,j)+aom(k,m,j)
              enddo
            enddo
          enddo
        enddo
      endif
c
      if (.not.allocated(eleca))    allocate (eleca(2,ngroup)  )
      if (.not.allocated(elecb))    allocate (elecb(2,ngroup)  )
      if (.not.allocated(rhoconda)) allocate (rhoconda(nmo,nmo))
      if (.not.allocated(rhocondb)) allocate (rhocondb(nmo,nmo))
      if (.not.allocated(rhotot))   allocate (rhotot(nmo,nmo)  )
c
c.....Obtain the resonance structures of [C] when an alpha MO is deleted
c
      nmobr=nmo
      nmobc=nmo
      nmob =nmo
      forall (i=1:nmo) ikeepbr(i)=i
      forall (i=1:nmo) ikeepbc(i)=i
      malcut=nalpha-1
      mbecut=nbeta
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
      nelv=nel-1
      combi=one
      do i=0,nelv-1
        combi=combi*dble(nelv+ngroup-1-i)/dble(nelv-i)
      enddo
      nprobt=int(combi+1D-10)
      if (.not.allocated(resnca)) allocate (resnca(nproba,ngroup))
      if (.not.allocated(resncb)) allocate (resncb(nprobb,ngroup))
      if (.not.allocated(probsa)) allocate (probsa(nproba))
      if (.not.allocated(probsb)) allocate (probsb(nprobb))
      if (.not.allocated(probt))  allocate (probt(nproba*nprobb))
      if (.not.allocated(pnew))   allocate (pnew(nprobt))
      if (.not.allocated(resnc))  allocate (resnc(nprobt,ngroup))
      if (.not.allocated(ind))    allocate (ind(ngroup))
      call rnprobs (eleca,resnca,malcut,ngroup,nproba,lw)
      call rnprobs (elecb,resncb,mbecut,ngroup,nprobb,lw)
      do irdel=1,nmo
        do icdel=1,irdel
          nmoar=0
          do i=1,nmo
            if (irdel.ne.i) then
              nmoar=nmoar+1
              ikeepar(nmoar)=i
            endif
          enddo
          nmoac=0
          do i=1,nmo
            if (icdel.ne.i) then
              nmoac=nmoac+1
              ikeepac(nmoac)=i
            endif
          enddo
          if (nmoar.ne.nmoac) then
            stop 'rhocond.f: Error deleting alpha MOs'
          endif
          nmoa=nmoar
          call semilla (idum)
          idum=-idum
          dumi=random(idum)
          if (nmoa.gt.0) then
            allocate (am(nproba,nproba))
            allocate (tpow(nproba,ngroup))
            allocate (ipvt(nproba))
            allocate (ovea(ngroup,nmoa,nmoa))
            allocate (oveaa(nmoa,nmoa))
            allocate (ww(nmoa))
            allocate (w(nproba))
            allocate (indx(nmoa))
            do m=1,nmoa
              do k=1,nmoa
                do igr=1,ngroup
                  ovea(igr,m,k)=sg(igr,ikeepar(m),ikeepac(k))
                enddo
              enddo
            enddo
            ok=.false.
            do while (.not.ok)
              am=zero
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
              call dgeco (am,nproba,nproba,ipvt,rcond,w)
              if (one + rcond .ne. one) ok=.true.
            enddo
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
              call dgeco (oveaa,nmoa,nmoa,indx,rcond,ww)
              job=10
              call dgedi (oveaa,nmoa,nmoa,indx,deter,ww,job)
              probsa(n)=deter(1)*tenp**deter(2)
            enddo
            job=0
            call dgesl (am,nproba,nproba,ipvt,probsa,job)
            deallocate (ovea,oveaa,ww,w,indx,am,tpow,ipvt)
          else
            probsa(1)=1d0
          endif
          if (nmob.gt.0) then
            allocate (am(nprobb,nprobb))
            allocate (tpow(nprobb,ngroup))
            allocate (ipvt(nprobb))
            allocate (ovea(ngroup,nmob,nmob))
            allocate (oveaa(nmob,nmob))
            allocate (ww(nmob))
            allocate (w(nprobb))
            allocate (indx(nmob))
            do m=1,nmob
              do k=1,nmob
                do igr=1,ngroup
                  ovea(igr,m,k)=sg(igr,ikeepbr(m),ikeepbc(k))
                enddo
              enddo
            enddo
            ok=.false.
            do while (.not.ok)
              am=zero
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
              call dgeco (am,nprobb,nprobb,ipvt,rcond,w)
              if (one + rcond .ne. one) ok=.true.
            enddo
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
              call dgeco (oveaa,nmob,nmob,indx,rcond,ww)
              job=10
              call dgedi (oveaa,nmob,nmob,indx,deter,ww,job)
              probsb(n)=deter(1)*tenp**deter(2)
            enddo
            job=0
            call dgesl (am,nprobb,nprobb,ipvt,probsb,job)
            deallocate (ovea,oveaa,ww,w,indx,am,tpow,ipvt)
          else
            probsb(1)=1d0
          endif
          ij=0
          do i=1,nproba
            do j=1,nprobb
              ij=ij+1
              probt(ij)=probsa(i)*probsb(j)
            enddo
          enddo
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
 1          enddo
          enddo
          do j=1,np
            isthiscond=.true.
            do i=1,ngroup
              if (resnc(j,i).ne.condens(i)) isthiscond=.false.
            enddo
            if (isthiscond) then
              jact=j
              goto 99
            endif
          enddo
 99       rhoconda(irdel,icdel)=(-1)**(irdel+icdel)*pnew(jact)
          rhoconda(icdel,irdel)=rhoconda(irdel,icdel)
        enddo
      enddo
      if (allocated(probsa)) deallocate (probsa)
      if (allocated(probsb)) deallocate (probsb)
      if (allocated(resnca)) deallocate (resnca)
      if (allocated(resncb)) deallocate (resncb)
      if (allocated(resnc )) deallocate (resnc )
      if (allocated(pnew  )) deallocate (pnew  )
      if (allocated(probt )) deallocate (probt )
      if (allocated(ind   )) deallocate (ind   )
c
c.....Obtain the resonance structures of [C] when a beta MO is deleted
c
      if (nalpha.eq.nbeta) then
        rhocondb(1:nmo,1:nmo)=rhoconda(1:nmo,1:nmo)
      endif
      if (nalpha.ne.nbeta) then
        nmoar=nmo
        nmoac=nmo
        nmoa= nmo
        do i=1,nmo
          ikeepar(i)=i
          ikeepac(i)=i
        enddo
        malcut=nalpha
        mbecut=nbeta-1
        do i=1,ngroup
          eleca(1,i)=izero
          elecb(1,i)=izero
          eleca(2,i)=malcut
          elecb(2,i)=mbecut
        enddo
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
        nelv=nel-1
        combi=one
        do i=0,nelv-1
          combi=combi*dble(nelv+ngroup-1-i)/dble(nelv-i)
        enddo
        nprobt=int(combi+1D-10)
        if (.not.allocated(resnca)) allocate (resnca(nproba,ngroup))
        if (.not.allocated(resncb)) allocate (resncb(nprobb,ngroup))
        if (.not.allocated(probsa)) allocate (probsa(nproba))
        if (.not.allocated(probsb)) allocate (probsb(nprobb))
        if (.not.allocated(probt))  allocate (probt(nproba*nprobb))
        if (.not.allocated(pnew))   allocate (pnew(nprobt))
        if (.not.allocated(resnc))  allocate (resnc(nprobt,ngroup))
        if (.not.allocated(ind))    allocate (ind(ngroup))
        call rnprobs (eleca,resnca,malcut,ngroup,nproba,lw)
        call rnprobs (elecb,resncb,mbecut,ngroup,nprobb,lw)
        nelv=nel-1
        do irdel=1,nmo
          do icdel=1,irdel
            nmobr=0
            do i=1,nmo
              if (irdel.ne.i) then
                nmobr=nmobr+1
                ikeepbr(nmobr)=i
              endif
            enddo
            nmobc=0
            do i=1,nmo
              if (icdel.ne.i) then
                nmobc=nmobc+1
                ikeepbc(nmobc)=i
              endif
            enddo
            if (nmobr.ne.nmobc) then
              stop 'rhocond.f: Error deleting beta MOs'
            endif
            nmob=nmobr
            call semilla (idum)
            idum=-idum
            dumi=random(idum)
            if (nmoa.gt.0) then
              allocate (am(nproba,nproba))
              allocate (tpow(nproba,ngroup))
              allocate (ipvt(nproba))
              allocate (ovea(ngroup,nmoa,nmoa))
              allocate (oveaa(nmoa,nmoa))
              allocate (ww(nmoa))
              allocate (w(nproba))
              allocate (indx(nmoa))
              do m=1,nmoa
                do k=1,nmoa
                  do igr=1,ngroup
                    ovea(igr,m,k)=sg(igr,ikeepar(m),ikeepac(k))
                  enddo
                enddo
              enddo
              ok=.false.
              do while (.not.ok) 
                am=zero
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
                call dgeco (am,nproba,nproba,ipvt,rcond,w)
                if (one + rcond .ne. one) ok=.true.
              enddo
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
                call dgeco (oveaa,nmoa,nmoa,indx,rcond,ww)
                job=10
                call dgedi (oveaa,nmoa,nmoa,indx,deter,ww,job)
                probsa(n)=deter(1)*tenp**deter(2)
              enddo
              job=0
              call dgesl (am,nproba,nproba,ipvt,probsa,job)
              deallocate (ovea,oveaa,ww,w,indx,am,tpow,ipvt)
            else
              probsa(1)=1d0
            endif
            if (nmob.gt.0) then
              allocate (am(nprobb,nprobb))
              allocate (tpow(nprobb,ngroup))
              allocate (ipvt(nprobb))
              allocate (ovea(ngroup,nmob,nmob))
              allocate (oveaa(nmob,nmob))
              allocate (ww(nmob))
              allocate (w(nprobb))
              allocate (indx(nmob))
              do m=1,nmob
                do k=1,nmob
                  do igr=1,ngroup
                    ovea(igr,m,k)=sg(igr,ikeepbr(m),ikeepbc(k))
                  enddo
                enddo
              enddo
              ok=.false.
              do while (.not.ok)
                am=zero
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
                call dgeco (am,nprobb,nprobb,ipvt,rcond,w)
                if (one + rcond .ne. one) ok=.true.
              enddo
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
                call dgeco (oveaa,nmob,nmob,indx,rcond,ww)
                job=10
                call dgedi (oveaa,nmob,nmob,indx,deter,ww,job)
                probsb(n)=deter(1)*tenp**deter(2)
              enddo
              job=0
              call dgesl (am,nprobb,nprobb,ipvt,probsb,job)
              deallocate (ovea,oveaa,ww,w,indx,am,tpow,ipvt)
            else
              probsb(1)=1d0
            endif
            ij=0
            do i=1,nproba
              do j=1,nprobb
                ij=ij+1
                probt(ij)=probsa(i)*probsb(j)
              enddo
            enddo
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
                      goto 1001
                    endif
                  enddo
                  np=np+1
                  pnew(np)=probt(n)
                  do i=1,ngroup
                    resnc(np,i)=resnca(ia,i)+resncb(ib,i)
                  enddo
                endif
 1001         enddo
            enddo
            do j=1,np
              isthiscond=.true.
              do i=1,ngroup
                if (resnc(j,i).ne.condens(i)) isthiscond=.false.
              enddo
              if (isthiscond) then
                jact=j
                goto 891
              endif
            enddo
 891        rhocondb(irdel,icdel)=(-1)**(irdel+icdel)*pnew(jact)
            rhocondb(icdel,irdel)=rhocondb(irdel,icdel)
          enddo
        enddo
      endif
c
c.....Analisis.
c
      rhotot(1:nmo,1:nmo)=rhoconda(1:nmo,1:nmo)+rhocondb(1:nmo,1:nmo)
      if (allocated(eleca   )) deallocate (eleca   )
      if (allocated(elecb   )) deallocate (elecb   )
      if (allocated(rhoconda)) deallocate (rhoconda)
      if (allocated(rhocondb)) deallocate (rhocondb)
      if (allocated(probsa  )) deallocate (probsa  )
      if (allocated(probsb  )) deallocate (probsb  )
      if (allocated(resnca  )) deallocate (resnca  )
      if (allocated(resncb  )) deallocate (resncb  )
      if (allocated(resnc   )) deallocate (resnc   )
      if (allocated(pnew    )) deallocate (pnew    )
      if (allocated(probt   )) deallocate (probt   )
      if (allocated(ind     )) deallocate (ind     )
      if (allocated(ikeepar )) deallocate (ikeepar )
      if (allocated(ikeepac )) deallocate (ikeepac )
      if (allocated(ikeepbr )) deallocate (ikeepbr )
      if (allocated(ikeepbc )) deallocate (ikeepbc )
      write (lw,87) (condens(i),i=1,ngroup)
      do i=1,nmo
       write (lw,88) (rhotot(i,j),j=1,i)
      enddo
      if (.not.allocated(rhoeig)) allocate (rhoeig(nmo))
      if (.not.allocated(rhovec)) allocate (rhovec(nmo,nmo))
      call jacobi (rhotot,nmo,nmo,rhoeig,rhovec,nrot)
      rhointcond=sum(rhoeig(1:nmo))
      write (lw,89) 
      write (lw,88) (rhoeig(i),i=1,nmo)
      write (lw,93) rhointcond
      do j=1,nmo
        write (lw,883) j
        write (lw,88) (rhovec(i,j),i=1,nmo)
      enddo
c
c.....Obtain group overlap integrals between natural MOs
c
      if (.not.allocated(sgnat)) allocate(sgnat(ngroup,nmo,nmo))
      sgnat=0D0
      do igr=1,ngroup
        do i=1,nmo
          do k=1,nmo
            tmpsg=0D0
            do m=1,nmo
              do j=1,nmo
                tmpsg=tmpsg+rhovec(m,i)*rhovec(j,k)*sg(igr,m,j)
              enddo
            enddo
            sgnat(igr,i,k)=tmpsg
          enddo
        enddo
      enddo
      do igr=1,ngroup
        write (lw,881) igr
        do i=1,nmo
          write (lw,88) (sgnat(igr,i,j),j=1,i)
        enddo
      enddo
c
c.....Probabilities of [S]
c
      if (.not.allocated(eleca))    allocate (eleca(2,ngroup)  )
      if (.not.allocated(elecb))    allocate (elecb(2,ngroup)  )
      malcut=nalpha
      mbecut=nbeta
      do i=1,ngroup
        eleca(1,i)=izero
        elecb(1,i)=izero
        eleca(2,i)=malcut
        elecb(2,i)=mbecut
      enddo
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
      nelv=nel
      combi=one
      do i=0,nelv-1
        combi=combi*dble(nelv+ngroup-1-i)/dble(nelv-i)
      enddo
      nprobt=int(combi+1D-10)
      if (.not.allocated(resnca)) allocate (resnca(nproba,ngroup))
      if (.not.allocated(resncb)) allocate (resncb(nprobb,ngroup))
      if (.not.allocated(probsa)) allocate (probsa(nproba))
      if (.not.allocated(probsb)) allocate (probsb(nprobb))
      if (.not.allocated(probt))  allocate (probt(nproba*nprobb))
      if (.not.allocated(pnew))   allocate (pnew(nprobt))
      if (.not.allocated(resnc))  allocate (resnc(nprobt,ngroup))
      if (.not.allocated(ind))    allocate (ind(ngroup))
      call rnprobs (eleca,resnca,malcut,ngroup,nproba,lw)
      call rnprobs (elecb,resncb,mbecut,ngroup,nprobb,lw)
      call semilla (idum)
      idum=-idum
      dumi=random(idum)
      allocate (am(nproba,nproba))
      allocate (tpow(nproba,ngroup))
      allocate (ipvt(nproba))
      allocate (ovea(ngroup,nmo,nmo))
      allocate (oveaa(nmo,nmo))
      allocate (ww(nmo))
      allocate (w(nproba))
      allocate (indx(nmo))
      do m=1,nmo
        do k=1,nmo
          do igr=1,ngroup
            ovea(igr,m,k)=sg(igr,m,k)
          enddo
        enddo
      enddo
      ok=.false.
      do while (.not.ok)
        am=zero
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
        call dgeco (am,nproba,nproba,ipvt,rcond,w)
        if (one + rcond .ne. one) ok=.true.
      enddo
      do n=1,nproba
        do m=1,nmo
          do k=1,nmo
            dumi=ovea(ngroup,m,k)
            do igr=1,ngroup-1
              dumi=dumi+tpow(n,igr)*ovea(igr,m,k)
            enddo
            oveaa(m,k)=dumi
          enddo
        enddo
        call dgeco (oveaa,nmo,nmo,indx,rcond,ww)
        job=10
        call dgedi (oveaa,nmo,nmo,indx,deter,ww,job)
        probsa(n)=deter(1)*tenp**deter(2)
      enddo
      job=0
      call dgesl (am,nproba,nproba,ipvt,probsa,job)
      deallocate (ovea,oveaa,ww,w,indx,am,tpow,ipvt)

      allocate (am(nprobb,nprobb))
      allocate (tpow(nprobb,ngroup))
      allocate (ipvt(nprobb))
      allocate (ovea(ngroup,nmo,nmo))
      allocate (oveaa(nmo,nmo))
      allocate (ww(nmo))
      allocate (w(nprobb))
      allocate (indx(nmo))
      do m=1,nmo
        do k=1,nmo
          do igr=1,ngroup
            ovea(igr,m,k)=sg(igr,m,k)
          enddo
        enddo
      enddo
      ok=.false.
      do while (.not.ok)
        am=zero
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
        call dgeco (am,nprobb,nprobb,ipvt,rcond,w)
        if (one + rcond .ne. one) ok=.true.
      enddo
      do n=1,nprobb
        do m=1,nmo
          do k=1,nmo
            dumi=ovea(ngroup,m,k)
            do igr=1,ngroup-1
               dumi=dumi+tpow(n,igr)*ovea(igr,m,k)
            enddo
            oveaa(m,k)=dumi
          enddo
        enddo
        call dgeco (oveaa,nmo,nmo,indx,rcond,ww)
        job=10
        call dgedi (oveaa,nmo,nmo,indx,deter,ww,job)
        probsb(n)=deter(1)*tenp**deter(2)
      enddo
      job=0
      call dgesl (am,nprobb,nprobb,ipvt,probsb,job)
      deallocate (ovea,oveaa,ww,w,indx,am,tpow,ipvt)

      ij=0
      do i=1,nproba
        do j=1,nprobb
          ij=ij+1
          probt(ij)=probsa(i)*probsb(j)
        enddo
      enddo

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
                goto 1007
              endif
            enddo
            np=np+1
            pnew(np)=probt(n)
            do i=1,ngroup
              resnc(np,i)=resnca(ia,i)+resncb(ib,i)
            enddo
          endif
 1007   enddo
      enddo

      do ig=1,ngroup
        do i=1,np
          resnc(i,ig)=resnc(i,ig)-1
        enddo
        popgrp=zero
        do i=1,nmo
          popgrp=popgrp+rhoeig(i)*sgnat(ig,i,i)
        enddo
        do i=1,np
          isthiscond=.true.
          do j=1,ngroup
            if (condens(j).ne.resnc(i,j)) isthiscond=.false.
          enddo
          if (isthiscond) then
            iact=i
          endif
        enddo
        do i=1,np
          resnc(i,ig)=resnc(i,ig)+1
        enddo
        popgrp=popgrp/pnew(iact)
        prsrs=pnew(iact)
        if (abs(prsrs).gt.zero) then
          write (lw,890) ig
          write (lw,88) (rhoeig(i)/prsrs,i=1,nmo)
        endif
        write (lw,772) ig,popgrp
      enddo
c
c.....Write the WFN file of the condensation [C]
c
      allocate (newcoef(nmo,nprims))
      do i=1,nmo
        do j=1,nprims
          temp = 0d0
          do k=1,nmo
            temp = temp + rhovec(k,i) * coef(k+nmo,j)
          enddo
          newcoef(i,j) = temp
        enddo
      enddo
      iwfnew=78
      write (0,112) filewfn(1:leng(filewfn))
      open (file=filewfn(1:leng(filewfn)),unit=iwfnew)
      rewind (iwfn)
      read (iwfn,1011) wfnttl
      write (iwfnew,'(a)') wfnttl(1:leng(wfnttl))//'      RHO[C]'
*     read (iwfn,1011) wfnttl
      read (iwfn,102) mode,nmodummy,nprimsdummy,ncentdummy
*     write (iwfnew,'(a)') wfnttl(1:leng(wfnttl))
      write (iwfnew,102) '    ',nmo,nprims,ncent
      do i=1,ncent
        read (iwfn,1011) wfnttl
        write (iwfnew,'(a)') wfnttl(1:LENG(wfnttl))
      enddo
*     read (iwfn,104) (icen(i),i=1,nprimsdummy)
      read (iwfn,104) (icendummy,i=1,nprimsdummy)
      nrec = nprims/20
      nres = mod (nprims,20)
      inic = 1
      ifin = 20
      do i=1,nrec
        write (iwfnew,1041) (icen(j),j=inic,ifin)
        inic=inic+20
        ifin=ifin+20
      enddo
      if (nres.gt.0) write (iwfnew,1041) (icen(i),i=20*nrec+1,nprims)
      read (iwfn,104) (itypdummy,i=1,nprimsdummy)
      inic = 1
      ifin = 20
      do i=1,nrec
        write (iwfnew,1042) (ityp(j),j=inic,ifin)
        inic=inic+20
        ifin=ifin+20
      enddo
      if (nres.gt.0) write (iwfnew,1042) (ityp(i),i=20*nrec+1,nprims)
      read (iwfn,105) (oexpdummy,i=1,nprimsdummy)
      nrec = nprims/5
      nres = mod (nprims,5)
      inic = 1
      ifin = 5
      do i=1,nrec
        write (iwfnew,1051) (oexp(j),j=inic,ifin)
        inic=inic+5
        ifin=ifin+5
      enddo
      if (nres.gt.0) write (iwfnew,1051) (oexp(i),i=5*nrec+1,nprims)
      do i=1,nmo      
        read (iwfn,106) occdummy,eorbdummy
        if (abs(prsrs).gt.zero) then
          write (iwfnew,1061) i,rhoeig(i),
     &      '   RENORM OCC =',rhoeig(i)/prsrs
        else
          write (iwfnew,1061) i,rhoeig(i),
     &      '  ORB. ENERGY =',0d0
        endif
        read (iwfn,107) (coefdummy,j=1,nprimsdummy)
        write (iwfnew,1071) (newcoef(i,j),j=1,nprims)
      enddo
      read (iwfn,108) check
      write (iwfnew,108) 'END DATA'
      read (iwfn,1011) wfnttl
      write (iwfnew,'(a)') wfnttl(1:leng(wfnttl))

      if (allocated(newcoef) ) deallocate (newcoef)
      if (allocated(condens) ) deallocate (condens)
      if (allocated(probsa ) ) deallocate (probsa )
      if (allocated(probsb ) ) deallocate (probsb )
      if (allocated(resnca ) ) deallocate (resnca )
      if (allocated(resncb ) ) deallocate (resncb )
      if (allocated(resnc  ) ) deallocate (resnc  )
      if (allocated(pnew   ) ) deallocate (pnew   )
      if (allocated(probt  ) ) deallocate (probt  )
      if (allocated(ind    ) ) deallocate (ind    )
      if (allocated(sg     ) ) deallocate (sg     )
      if (allocated(sgnat  ) ) deallocate (sgnat  )
      if (allocated(rhoeig ) ) deallocate (rhoeig )
      if (allocated(rhovec ) ) deallocate (rhovec )
c
      call getdate (line2)
      write (lw,4000) line2(1:leng(line2))
      close (unit=lw)
      close (unit=iwfnew)
      return
c
c.....Formats.
c
 100  format (' # ',F22.16,1x,12I6)
 1000 format(' # Calculation starts on ',a)
 4000 format (' # Calculation ends on ',a)
 87   format (' # CGDM rho[C] in the canonical basis for [C] = ',10I3)
 88   format (6(1x,E16.10))
 89   format (' # Eigenvalues (eta_i)')
 93   format (' # Integral to R^3 of rho(r)[C] = ',F16.10)
 991  format (' # Eigenvectors')
 890  format (' # Modified Eigenvalues eta_i/P(S) in FRAGMENT ',I3)
 77   format (' # Probabilities for all the the RSRS [S]')
 881  format (' # AOM from rho[C] in group ',I3)
 883  format (' # Eigenvector ',I3)
 772  format (' # Integral of Renormalized RHO[C] IN FRAGMENT ',I3,
     &        ' = ',F16.10)
 1011 format (A80)
 104  format (20x,20i3)
 1041 FORMAT ('CENTRE ASSIGNMENTS  ',20I3)
 1042 FORMAT ('TYPE ASSIGNMENTS    ',20I3)
 1051 FORMAT ('EXPONENTS ',5(1PE14.7))
 1061 FORMAT ('MO',10x,I3,11X,'OCC NO = ',F12.8,A15,E18.11)
 105  FORMAT (10X,5E14.7)
 106  FORMAT (35X,F12.8,15X,F12.8)
 107  FORMAT (5E16.8)
 108  FORMAT (A8)
 1071 FORMAT (5(1PE16.8))
 112  FORMAT (1X,'# Writing file ',a)
 102  FORMAT (4X,A4,10X,3(I5,15X))
      end
