c
c-----------------------------------------------------------------------
c
      subroutine critcplx (pcut,npais,npbis,ngroup,nmoa,nmob,lw,
     &     eleca,elecb,sga,sgb,rsrsa,rsrsb,mocore)
c
c-----Heavily modified version of 'nbasab.f' to deal with complex AOMs
c
c.......................................................................
c
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
c
      real(kind=8), allocatable,dimension (:,:) :: tpow
      complex*16,   allocatable,dimension (:,:) :: am,oveaa
      complex*16,   allocatable,dimension (:)   :: probs
      complex*16 determ

      real   (kind=8), allocatable,dimension (:,:) :: probr
      real   (kind=8), allocatable,dimension (:)   :: w,ww,probt,pnew
      real   (kind=8), allocatable,dimension (:)   :: p1,p1a,p2,d1,p3,d2
      real   (kind=8), allocatable,dimension (:)   :: p1b
      real   (kind=8), allocatable,dimension (:)   :: xli
      real   (kind=8), allocatable,dimension (:)   :: probord
      real   (kind=8), allocatable,dimension (:,:,:) :: sg
c
      integer(kind=4), allocatable,dimension (:)   :: ipvt,indx,ind
      integer(kind=4), allocatable,dimension (:,:) :: resnc
      integer(kind=4), allocatable,dimension (:,:) :: igtzero
      integer(kind=4), allocatable,dimension (:,:) :: resncord
      integer(kind=4), allocatable,dimension (:)   :: npop
      integer(kind=4), allocatable,dimension (:)   :: ioprob
      integer(kind=4), allocatable,dimension (:,:) :: rsrs
      integer(kind=4), allocatable,dimension (:,:,:) :: resalp
c
      complex*16       sga(ngroup,nmoa,nmoa)
      complex*16       sgb(ngroup,nmob,nmob)
      complex*16       dumaom
      real   (kind=8)  dumi,random,deter(2)
      real   (kind=8)  rcond
      integer(kind=4)  rsrsa(npais,ngroup)
      integer(kind=4)  rsrsb(npbis,ngroup)
      integer(kind=4)  mocore(ngroup)
      integer(kind=4)  eleca(2,ngroup)
      integer(kind=4)  elecb(2,ngroup)
      integer(kind=4)  idum
      integer(kind=4)  nnonz(2)
      logical          inlist,singular
      character(len=5) spin

      real(kind=8) :: numrandom(1000)
c
      call timer (2,ipid,'_critcplx ',-1)
c
c-----Generates 40 random numbers in the range [0,1] and stores them
c     in numrandom(1:40)
c     call rndnumber (numrandom,40)
c
c-----Compute the maximum of the number of alpha and beta probabilities
c
      npmax = max(npais,npbis)
      if (.not.allocated(probr)) then
        allocate (probr(npmax,2),stat=ier)
        if (ier.ne.0) stop 'critcplx.f: Cannot allocate probr()'
      endif
      if (.not.allocated(rsrs)) then
        allocate (rsrs(npmax,ngroup),stat=ier)
        if (ier.ne.0) stop 'critcplx.f: Cannot allocate rsrs()'
      endif
c
c-----We consider separately the ALPHA and BETA probabilities
c
      ispin = 1
 22   continue
      call semilla (idum)
      idum = -idum
      dumi = random(idum)
      if (ispin.eq.1) then
        nmois = nmoa
        npab = npais
        spin = 'ALPHA'
        rsrs(1:npais,1:ngroup) = rsrsa(1:npais,1:ngroup)
        allocate (sg(ngroup,nmois,nmois))
        sg(1:ngroup,1:nmois,1:nmois) = sga(1:ngroup,1:nmois,1:nmois)
      else
        nmois = nmob
        npab = npbis
        spin = 'BETA '
        rsrs(1:npbis,1:ngroup) = rsrsb(1:npbis,1:ngroup)
        allocate (sg(ngroup,nmois,nmois))
        sg(1:ngroup,1:nmois,1:nmois) = sgb(1:ngroup,1:nmois,1:nmois)
      endif

      if (.not.allocated(oveaa)) then
        allocate (oveaa(nmois,nmois),stat=ier)
        if (ier.ne.0) stop 'critcplx.f: Cannot allocate oveaa()'
      endif
      if (.not.allocated(ww)) then
        allocate (ww(nmois),stat=ier)
        if (ier.ne.0) stop 'critcplx.f: Cannot allocate ww()'
      endif
      if (.not.allocated(indx)) then
        allocate (indx(nmois),stat=ier)
        if (ier.ne.0) stop 'critcplx.f: Cannot allocate indx()'
      endif
      if (.not.allocated(am)) then
        allocate (am(npab,npab),stat=ier)
        if (ier.ne.0) stop 'critcplx.f: Cannot allocate am()'
      endif
      if (.not.allocated(tpow)) then
        allocate (tpow(npab,ngroup),stat=ier)
        if (ier.ne.0) stop 'critcplx.f: Cannot allocate tpow()'
      endif
      if (.not.allocated(ipvt)) then
        allocate (ipvt(npab),stat=ier)
        if (ier.ne.0) stop 'critcplx.f: Cannot allocate ipvt()'
      endif
      if (.not.allocated(w)) then
        allocate (w(npab),stat=ier)
        if (ier.ne.0) stop 'critcplx.f: Cannot allocate w()'
      endif
      if (.not.allocated(probs)) then
        allocate (probs(npab),stat=ier)
        if (ier.ne.0) stop 'critcplx.f: Cannot allocate probs()'
      endif
c
c.....Construct the first member of the linear system.
c
      epscond = 1.0d-30
      singular = .true.
      do while (singular)
        am(1:npab,1:npab)=(0.0d+00,0.0d+00)
        tpow=0.0d+00
        do i=1,npab
          do k=1,ngroup-1
            tpow(i,k)=2.0d+00*random(idum)-1.0d+00
          enddo
        enddo
        do i=1,npab
          do j=1,npab
            aco=1.0d+00
            do k=1,ngroup-1
              aco=aco*tpow(i,k)**rsrs(j,k)
            enddo
            am(i,j)=cmplx(aco,0.0d+00)
          enddo
        enddo
        call timer (2,ipid1,'_zgeco    ',-1)
        call zgeco (am,npab,npab,ipvt,rcond,w)
        call timer (4,ipid1,'_zgeco    ',-1)
        write (lw,444) rcond
        if (abs(rcond).ge.epscond) singular=.false.
      enddo
c
c.....Construct the second member of the linear system.
c
      call timer (2,ipid2,'_determ   ',-1)
      do n=1,npab
        do m=1,nmois
          do k=1,nmois
            dumaom=sg(ngroup,m,k)
            do igr=1,ngroup-1
               dumaom=dumaom+tpow(n,igr)*sg(igr,m,k)
            enddo
            oveaa(m,k)=dumaom
          enddo
        enddo
c
c.......Determinants calculation.
c
        forall (i=1:nmois) indx(i) = i
        call zgeco (oveaa,nmois,nmois,indx,rcond,ww)
        job=10
        call zgedi (oveaa,nmois,nmois,indx,deter,ww,job,determ)
        probs(n)=determ
      enddo
c
      call timer (4,ipid2,'_determ   ',-1)
c
c.....Linear System Solver: Netlib DGECO routine.
c 
      call timer (2,ipid3,'_zgesl    ',-1)
      job=0
      call zgesl (am,npab,npab,ipvt,probs,job)
      call timer (4,ipid3,'_zgesl    ',-1)
c
      do n=1,npab
        probr(n,ispin)=dreal(probs(n))
      enddo
      write (lw,200) trim(spin),ngroup,npab

      if (.not.allocated(resalp)) then
        allocate (resalp(ngroup,npab,2),stat=ier)
        if (ier.ne.0) stop 'critcplx.f: Cannot allocate resalp()'
      endif
      if (.not.allocated(igtzero)) then
        allocate (igtzero(npab,2),stat=ier)
        if (ier.ne.0) stop 'critcplx.f: Cannot allocate igtzero()'
      endif
      nnonz(ispin) = 0
      sumgt = 0.0d+00
      sumtot = sum(probr(1:npab,ispin))
      do n=1,npab
        if (probr(n,ispin).ge.pcut) then
          nnonz(ispin) = nnonz(ispin)+1
          igtzero(nnonz(ispin),ispin) = n
          nnz = nnonz(ispin)
          resalp(1:ngroup,nnz,ispin) = rsrs(n,1:ngroup)+mocore(1:ngroup)
          sumgt = sumgt+probr(n,ispin)
          write (lw,100) probr(n,ispin),
     &         (resalp(igr,nnz,ispin),igr=1,ngroup)
        endif
      enddo
      write (lw,1060) sumgt, nnonz(ispin), pcut, sumtot
      deallocate (tpow,am,ipvt,w,oveaa,ww,indx)
      if (ispin.eq.1.and.nspin.eq.2) then
        ispin = 2
        if (allocated(am)) then
          deallocate (am,stat=ier)
          if (ier.ne.0) stop 'critcplx.f: Cannot deallocate am()'
        endif
        if (allocated(tpow)) then
          deallocate (tpow,stat=ier)
          if (ier.ne.0) stop 'critcplx.f: Cannot deallocate tpow()'
        endif
        if (allocated(ipvt)) then
          deallocate (ipvt,stat=ier)
          if (ier.ne.0) stop 'critcplx.f: Cannot deallocate ipvt()'
        endif
        if (allocated(w)) then
          deallocate (w,stat=ier)
          if (ier.ne.0) stop 'critcplx.f: Cannot deallocate w()'
        endif
        if (allocated(oveaa)) then
          deallocate (oveaa,stat=ier)
          if (ier.ne.0) stop 'critcplx.f: Cannot deallocate oveaa()'
        endif
        if (allocated(ww)) then
          deallocate (ww,stat=ier)
          if (ier.ne.0) stop 'critcplx.f: Cannot deallocate ww()'
        endif
        if (allocated(indx)) then
          deallocate (indx,stat=ier)
          if (ier.ne.0) stop 'critcplx.f: Cannot deallocate indx()'
        endif
        if (allocated(probs)) then
          deallocate (probs,stat=ier)
          if (ier.ne.0) stop 'critcplx.f: Cannot deallocate probs()'
        endif
        goto 22
      else
        npbis = npais
        nnonz(2) = nnonz(1)
        igtzero(1:npab,2) = igtzero(1:npab,1)
        probr(:,2) = probr(:,1)
        rsrsb(1:npab,:) = rsrsa(1:npab,:)
        resalp(1:ngroup,1:npab,2) = resalp(1:ngroup,1:npab,1)
      endif
c
c.....Total spin-splitted probabilities 
c
      allocate (probt(npais*npbis))
      nprobg = 0
      sumgt = 0.0d+00
      sumtot = 0.0d+00
      write (lw,2001) ngroup,npais*npbis,nnonz(1)*nnonz(2)
      ij = 0
      do i=1,nnonz(1)
        do j=1,nnonz(2)
          ij = ij+1
          probt(ij)=probr(igtzero(i,1),1)*probr(igtzero(j,2),2)
          sumtot = sumtot+probt(ij)
          if (probt(ij).ge.pcut) then
            nprobg=nprobg+1
            sumgt = sumgt+probt(ij)
            write (lw,100) probt(ij),
     &          (resalp(igr,i,1),resalp(igr,j,2),igr=1,ngroup)
          endif
        enddo
      enddo
      write (lw,106) sumgt, nprobg, pcut, sumtot
c
c.....Average population of alpha and beta electrons in each group.
c
      allocate (p1a(ngroup))
      p1a = 0.0d+00
      do i=1,npais
        do j=1,ngroup
          p1a(j) = p1a(j) + rsrsa(i,j) * probr(i,1)
        enddo
      enddo
      allocate (p1b(ngroup))
      p1b = 0.0d+00
      do i=1,npbis
        do j=1,ngroup
          p1b(j) = p1b(j) + rsrsb(i,j) *probr(i,2)
        enddo
      enddo

      write (lw,111)
      do i=1,ngroup
        write (lw,151) i,'ALPHA',p1a(i)
        write (lw,151) i,'BETA ',p1b(i)
      enddo
c
c.....Obtain multiple-fragment electron population covariances.
c
      call sndelta (probr,p1a,p1b,rsrsa,rsrsb,ngroup,
     &     npais,npbis,npmax,lw)
c
c.....Computes spinless probabilities from spin resolved probabilities.
c     Obtain the number of spinless real space resonant structures.
c
      nprobt = 1
      do i=1,ngroup-1
        nprobt = nprobt*(eleca(2,i)+elecb(2,i)-eleca(1,i)-elecb(1,i)+1)
      enddo
      allocate (pnew(nprobt))
      allocate (resnc(nprobt,ngroup))
      allocate (ind(ngroup))
c
      pnew = 0.0d+00
      n = 0
      do ia=1,nnonz(1)
        cycleib: do ib=1,nnonz(2)
          n=n+1
          if (n.eq.1) then
            np=1
            do i=1,ngroup
              resnc(np,i)=resalp(i,ia,1)+resalp(i,ib,2)
            enddo
            pnew(np) = probt(n)
          else
            do i=1,ngroup
              ind(i)=resalp(i,ia,1)+resalp(i,ib,2)
            enddo
            do j=1,np
              inlist = .true.
              do k=1,ngroup
                inlist = inlist.and.(ind(k).eq.resnc(j,k))
              enddo
              if (inlist) then
                pnew(j) = pnew(j)+probt(n)
                cycle cycleib
              endif
            enddo
            np = np+1
            pnew(np) = probt(n)
            do i=1,ngroup
              resnc(np,i)=resalp(i,ia,1)+resalp(i,ib,2)
            enddo
          endif
        enddo cycleib
      enddo
      deallocate (probt,ind)
c
c.....Determine if probabilities are ordered or not.
c
      write (lw,202) ngroup,np
      npnew = 0
      sumtot = sum(pnew(1:np))
      sumgt = 0.0d+00
      do n=1,np
        if (pnew(n).ge.pcut) then
          npnew = npnew+1
          sumgt = sumgt+pnew(n)
        endif
      enddo
c
      if (npnew.le.nstack) then
c
c.......Order by increasing value the highest probabilities.
c
        allocate (ioprob (npnew))
        allocate (resncord (npnew,ngroup))
        allocate (probord(npnew))
        npronew = 0
        do n=1,np
          if (pnew(n).ge.pcut) then
            npronew = npronew+1
            probord(npronew) = pnew(n) 
            ioprob(npronew) = npronew
            do igr=1,ngroup
              resncord(npronew,igr) = resnc(n,igr)
            enddo
          endif
        enddo
        call qqsort (probord,ioprob,1,npnew,nstack)
        do n=npnew,1,-1
          m = ioprob(n)
          write (lw,100) probord(m),(resncord(m,igr),igr=1,ngroup)
        enddo
        deallocate (ioprob,resncord,probord)
      else
c
c.......Probabilities are not ordered.
c
        do n=1,np
          if (pnew(n).ge.pcut) then
             write (lw,100) pnew(n),(resnc(n,igr),igr=1,ngroup)
          endif
        enddo
      endif
      write (lw,106) sumgt, npnew, pcut, sumtot
c
c.....Compute spinless localization and delocalization indices
c
      allocate (p1(ngroup),xli(ngroup))
      p1 = 0.0d+00
      xli = 0.0d+00
      allocate (npop(ngroup))
      if (ngroup.gt.1) then
        allocate (p2(ngroup*(ngroup-1)/2))
        allocate (d1(ngroup*(ngroup-1)/2))
        p2 = 0.0d+00
        d1 = 0.0d+00
      endif
      if (ngroup.gt.2) then
        allocate (p3(ngroup*(ngroup-1)*(ngroup-2)/6))
        allocate (d2(ngroup*(ngroup-1)*(ngroup-2)/6))
        p3 = 0.0d+00
        d2 = 0.0d+00
      endif
c
      do npa=1,npais
        do npb=1,npbis
          p = probr(npa,1)*probr(npb,2)
          do i=1,ngroup
            npop(i) = rsrsa(npa,i)+rsrsa(npb,i) + 2*mocore(i)
            p1(i) = p1(i)+npop(i)*p
            xli(i) = xli(i)+p*npop(i)*npop(i)
          enddo
          if (ngroup.gt.1) then
            ipair = 0
            do i=1,ngroup
              do j=1,i-1
                ipair = ipair+1
                p2(ipair) = p2(ipair)+npop(i)*npop(j)*p
              enddo
            enddo
          endif
          if (ngroup.gt.2) then
            iter = 0
            do i=1,ngroup
              do j=1,i-1
                do k=1,j-1
                iter = iter+1
                p3(iter) = p3(iter)+npop(i)*npop(j)*npop(k)*p
                enddo
              enddo
            enddo
          endif
        enddo
      enddo
      write (lw,11)
      do i=1,ngroup
        write (lw,15) i,p1(i)
      enddo
      write (lw,150) sum(p1(1:ngroup))
      if (ngroup.gt.1) then
        ipair = 0
        do i=1,ngroup
          do j=1,i-1
            ipair = ipair+1
            write (lw,2) i,j,p2(ipair)
          enddo
        enddo
      endif
c
      if (ngroup.gt.2) then
        iter = 0
        do i=1,ngroup
          do j=1,i-1
            do k = 1,j-1
              iter = iter+1
              write (lw,3) i,j,k,p3(iter)
            enddo
          enddo
        enddo
      endif
c
      do npa=1,npais
        do npb=1,npbis
          do i=1,ngroup
            npop(i)=rsrsa(npa,i)+rsrsa(npb,i) + 2*mocore(i)
          enddo
          p = probr(npa,1)*probr(npb,2)
          if (ngroup.gt.1) then
            ipair = 0
            do i=1,ngroup
              do j=1,i-1
                ipair = ipair+1
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
      enddo
c
      do i=1,ngroup
        xli(i)=-(xli(i)-p1(i)*(p1(i)+1.0d+00))
        write (lw,6) i,i,xli(i),hundred*xli(i)/p1(i)
      enddo
      write (lw,33)
      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (lw,4) i,j,d1(ipair)
          enddo
        enddo
      endif
      if (ngroup.gt.2) then
        iter=0
        do i=1,ngroup
          do j=1,i-1
            do k=1,j-1
              iter=iter+1
              write (lw,5) i,j,k,d2(iter)
            enddo
          enddo
        enddo
      endif
      deallocate (p1,xli,npop)
      if (ngroup.gt.1) deallocate (p2,d1)
      if (ngroup.gt.2) deallocate (p3,d2)
c
      return
c
      call timer (4,ipid,'_critcplx ',-1)
      return
c
c.....Formats.
c
 200  format (/,' # M-BASINS EDFs (COMPLEX AOM VERSION) ',
     & 'FOR ',a,' ELECTRONS',/,1x,'#',72('-'),/,
     &  ' # NUMBER OF GROUPS',13x,' = ',I8,/,
     &  ' # TOTAL NUMBER OF PROBABILITIES = ',I8,/,
     &   1x,'#',72('-'),/,
     & ' #     Probability            n1    n2    n3 ...')
 2001 format (/,' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION',
     & ' INCLUDING SPIN', /,1x, '# ',72('-'),/,
     & ' # NUMBER OF GROUPS                              =  ', I8,/,
     & ' # TOTAL AND PROCESSED  NUMBER OF PROBABILITIES  = ',2(1x,I8),/,
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
 44   format (/,1x,'# Linear System Solver: Netlib DGECO routine',/)
 444  format (1x,'# RECIPROCAL OF RCOND VALUE IS ',E18.12)
 100  format (' # ',F22.15,1x,12I6)
 1060 format (1x,'#',72('-'),/,' # ',F22.15,2x,'<-- SUM,',I8,
     & ' PROBABILITIES > ',E16.9,/,' # ',F22.15,2x,
     &   '<-- TOTAL SUM',/,1x,'#',72('-'))
 106  format (1x,'#',72('-'),/,' # ',F22.15,2x,'<-- SUM,',I8,
     & ' PROBABILITIES > ',E16.9,/,' # ',F22.15,2x,
     &   '<-- TOTAL SUM (NON NECESSARILY EQUAL TO 1.0)',
     &   /,1x,'#',72('-'))
 11   format (//1x,'Average populations and localization indices')
 33   format (//1x,'Delocalization indices,',
     &             ' Eq. (28) J. Chem. Phys.  126, 094102 (2007)')
 15   format (1x,'<n(',I3,')>               = ',F16.9)
 150  format (1x,'Total population       = ',F16.9)
 2    format (1x,'<n(',I3,') n(',I3,')>        = ',F16.9)
 3    format (1x,'<n(',I3,') n(',I3,') n(',I3,')> = ',F16.9)
 4    format (1x,'delta_(',I3,I3,')         = ',F16.9)
 5    format (1x,'delta_(',I3,I3,I3,')      = ',F16.9)
 6    format (1x,'delta_(',I3,I3,')         = ',F16.9,
     &        2x,'% Localization = ',F16.9)
 111  format (/,' # ',35('+'),/,' # Average electron populations',/,
     & ' # ',35('+'))
 151  format (1x,'# <n(',I3,')>_{',a,'} = ',F16.9)
 112  format (1x,'#',/,1x,'# multiple-group delocalization indices',/,
     &        1x,'#')
 113  format (1x,'DELTA (alpha,beta,total) = ',3(1x,F16.9),5x,20(I3))
 1120 format (3(1x,'#',/),1x,'# ',80('+'),/,
     &        1x,'# EXACT CALCULATION OF PROBABILITIES',/,1x,'#')
 206  format (1x,'#',/,1x,'# THERE ARE NO ',A,' ELECTRONS',/,1x,'#')
      end
