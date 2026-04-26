c
c-----------------------------------------------------------------------
c
      subroutine compedfstd (pcut,npab,nprev,ngroup,nmo,lw,mimael,sg,
     &                       rsrs,mocore)
c
c.......................................................................
c
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
c
      real   (kind=8), allocatable,dimension (:,:) :: tpow
      real   (kind=8), allocatable,dimension (:,:) :: am,oveaa
      real   (kind=8), allocatable,dimension (:)   :: probs
      real   (kind=8), allocatable,dimension (:,:) :: probr
      real   (kind=8), allocatable,dimension (:)   :: w,ww,probt,pnew
      real   (kind=8), allocatable,dimension (:)   :: p1,p1a,p2,d1,p3
      real   (kind=8), allocatable,dimension (:)   :: d2,xli
      real   (kind=8), allocatable,dimension (:)   :: sp1,sxli,sp2,sp3
      real   (kind=8), allocatable,dimension (:)   :: sd1,sd2
      real   (kind=8), allocatable,dimension (:)   :: probord
      integer(kind=4), allocatable,dimension (:)   :: ipvt,indx,ind
      integer(kind=4), allocatable,dimension (:,:) :: resnc,resalp
      integer(kind=4), allocatable,dimension (:,:) :: rsrsx
      integer(kind=4), allocatable,dimension (:)   :: igtzero
      integer(kind=4), allocatable,dimension (:,:) :: resncord
      integer(kind=4), allocatable,dimension (:)   :: npop
      integer(kind=4), allocatable,dimension (:)   :: ioprob
      real   (kind=8)  sg(ngroup,nmo,nmo), dumaom
      real   (kind=8)  dumi,random,deter(2)
      real   (kind=8)  rcond
      integer(kind=4)  rsrs(nprev,ngroup),lw
      integer(kind=4)  mocore(ngroup)
      integer(kind=4)  mimael(2,ngroup)
      integer(kind=4)  idum
      logical          inlist,singular
c
      call timer (2,ipidmain,'_compedfst',-1)
c
c.....Set the number of electrons two twice the number of MOs
c
      call semilla (idum)
      idum=-idum
      dumi=random(idum)

      allocate ( am(npab,npab) )
      allocate ( tpow(npab,ngroup) )
      allocate ( ipvt(npab) )
      allocate ( w(npab) )
      allocate ( oveaa(nmo,nmo) )
      allocate ( ww(nmo) )
      allocate ( indx(nmo) )
      allocate ( probs(npab) )
      allocate ( probr(npab,2) )
c
c.....Construct the first member of the linear system.
c
      epscond=1d-30
      singular=.true.
      ii=0
      do while (singular)
        ii=ii+1
        if (ii.eq.100) then
          write (lw,*) '# compedfstd.f: Singular am() matrix. abort'
          stop
        endif
        am(1:npab,1:npab)=0.0d+00
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
            am(i,j)=aco
          enddo
        enddo
        call timer (2,ipidgeco,'_dgeco    ',-1)
        call dgeco (am,npab,npab,ipvt,rcond,w)
        call timer (4,ipidgeco,'_dgeco    ',-1)
*       write (lw,444) rcond
        if (abs(rcond).ge.epscond) singular=.false.
      enddo
c
c.....Construct the second member of the linear system.
c
      call timer (2,ipiddet,'_determ   ',-1)
!$omp parallel
!$omp& private(oveaa,indx,ww,deter,dumaom,rcond,m,k,igr,job)
      allocate(oveaa(nmo,nmo), indx(nmo), ww(nmo))
!$omp do schedule(dynamic)
      do n=1,npab
        do m=1,nmo
          do k=1,nmo
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
        call dgeco (oveaa,nmo,nmo,indx,rcond,ww)
        job=10
        call dgedi (oveaa,nmo,nmo,indx,deter,ww,job)
        probs(n)=deter(1)*tenp**deter(2)
      enddo
!$omp end do
      deallocate(oveaa, indx, ww)
!$omp end parallel
      call timer (4,ipiddet,'_determ   ',-1)
c
c.....Linear System Solver: Netlib DGECO routine.
c 
      call timer (2,ipidgesl,'_dgesl    ',-1)
      job=0
      call dgesl (am,npab,npab,ipvt,probs,job)
      call timer (4,ipidgesl,'_dgesl    ',-1)
c
      do n=1,npab
        probr(n,1)=probs(n)
      enddo
      write (lw,200) ngroup,npab

      if (.not.allocated(resalp)) then
        allocate (resalp(ngroup,npab),stat=ier)
        if (ier.ne.0) stop 'compedfstd.f: Cannot allocate resalp()'
      endif
      if (.not.allocated(igtzero)) then
        allocate (igtzero(npab),stat=ier)
        if (ier.ne.0) stop 'compedfstd.f: Cannot allocate igtzero()'
      endif
      nnonz=0
      sum=0.0d+00
      sumtot=0.0d+00
      do n=1,npab
        sumtot=sumtot+probr(n,1)
        if (probr(n,1).ge.pcut) then
          nnonz=nnonz+1
          igtzero(nnonz)=n
          resalp(1:ngroup,nnonz)=rsrs(n,1:ngroup)+mocore(1:ngroup)
          sum=sum+probr(n,1)
          write (lw,100) probr(n,1),(resalp(igr,nnonz),igr=1,ngroup)
        endif
      enddo
      write (lw,1060) sum, nnonz, pcut, sumtot
      deallocate (tpow,am,ipvt,w,oveaa,ww,indx)
c
c-----Closed shell, so that alpha and beta EDFs are equal
c
      probr(1:npab,2)=probr(1:npab,1)
c
c.....Total spin-splitted probabilities 
c
      allocate (probt(npab*npab))
      nprobg=0
      sum=0.0d+00
      sumtot=0.0d+00
      write (lw,2001) ngroup,npab*npab,nnonz*nnonz
      ij=0
      do i=1,nnonz
        do j=1,nnonz
          ij=ij+1
          probt(ij)=probr(igtzero(i),1)*probr(igtzero(j),2)
          sumtot=sumtot+probt(ij)
          if (probt(ij).ge.pcut) then
            nprobg=nprobg+1
            sum=sum+probt(ij)
            write (lw,100) probt(ij),
     &          (resalp(igr,i),resalp(igr,j),igr=1,ngroup)
          endif
        enddo
      enddo
      write (lw,106) sum, nprobg, pcut, sumtot
c
c.....Average population of alpha and beta electrons in each group.
c
      allocate (p1a(ngroup))
      p1a=0.0d+00
      do i=1,npab
        do j=1,ngroup
          p1a(j)=p1a(j)+rsrs(i,j)*probr(i,1)
        enddo
      enddo
      write (lw,111)
      do i=1,ngroup
        write (lw,151) i,'ALPHA or BETA',p1a(i)
      enddo
c
c.....Obtain multiple-fragment electron population covariances.
c
      allocate (rsrsx(npab,ngroup))
      rsrsx(1:npab,1:ngroup)=rsrs(1:npab,1:ngroup)
      call sndelta (probr,p1a,p1a,rsrsx,rsrsx,ngroup,npab,npab,npab,lw)
      deallocate (rsrsx)
c
c.....Computes spinless probabilities from spin resolved probabilities.
c     Obtain the number of spinless real space resonant structures.
c
      nprobt=1
      do i=1,ngroup-1
        nprobt=nprobt*(2*mimael(2,i)-2*mimael(1,i)+1)
      enddo
      allocate (pnew(nprobt))
      allocate (resnc(nprobt,ngroup))
      allocate (ind(ngroup))
c
      pnew=0.0d+00
      n=0
      do ia=1,nnonz
        do ib=1,nnonz
          n=n+1
          if (n.eq.1) then
            np=1
            do i=1,ngroup
              resnc(np,i)=resalp(i,ia)+resalp(i,ib)
            enddo
            pnew(np)=probt(n)
          else
            do i=1,ngroup
              ind(i)=resalp(i,ia)+resalp(i,ib)
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
              resnc(np,i)=resalp(i,ia)+resalp(i,ib)
            enddo
          endif
 1      enddo
      enddo
      deallocate (probt,ind)
c
c.....Determine whether probabilities are ordered or not.
c
      write (lw,202) ngroup,np
      npnew=0
      sumtot=0.0d+00
      sum=0.0d+00
      do n=1,np
        sumtot=sumtot+pnew(n)
        if (pnew(n).ge.pcut) then
          npnew=npnew+1
          sum=sum+pnew(n)
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
        npronew=0
        do n=1,np
          if (pnew(n).ge.pcut) then
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
      write (lw,106) sum, npnew, pcut, sumtot
c
c.....Compute spinless localization and delocalization indices
c
      allocate (p1(ngroup),xli(ngroup))
      p1=0.0d+00
      xli=0.0d+00
      allocate (npop(ngroup))
      if (ngroup.gt.1) then
        allocate (p2(ngroup*(ngroup-1)/2))
        allocate (d1(ngroup*(ngroup-1)/2))
        p2=0.0d+00
        d1=0.0d+00
      endif
      if (ngroup.gt.2) then
        allocate (p3(ngroup*(ngroup-1)*(ngroup-2)/6))
        allocate (d2(ngroup*(ngroup-1)*(ngroup-2)/6))
        p3=0.0d+00
        d2=0.0d+00
      endif
c
      ngpair = ngroup*(ngroup-1)/2
      ntripl = ngroup*(ngroup-1)*(ngroup-2)/6
!$omp parallel
!$omp& private(npop,sp1,sxli,sp2,sp3,p,i,j,k,npb,ipair,iter)
      allocate(npop(ngroup), sp1(ngroup), sxli(ngroup))
      sp1 = 0.0d+00
      sxli = 0.0d+00
      if (ngroup.gt.1) then
        allocate(sp2(ngpair))
        sp2 = 0.0d+00
      endif
      if (ngroup.gt.2) then
        allocate(sp3(ntripl))
        sp3 = 0.0d+00
      endif
!$omp do schedule(dynamic)
      do npa=1,npab
        do npb=1,npab
          p=probr(npa,1)*probr(npb,2)
          do i=1,ngroup
            npop(i)=rsrs(npa,i)+rsrs(npb,i)+2*mocore(i)
            sp1(i)=sp1(i)+npop(i)*p
            sxli(i)=sxli(i)+p*npop(i)*npop(i)
          enddo
          if (ngroup.gt.1) then
            ipair=0
            do i=1,ngroup
              do j=1,i-1
                ipair=ipair+1
                sp2(ipair)=sp2(ipair)+npop(i)*npop(j)*p
              enddo
            enddo
          endif
          if (ngroup.gt.2) then
            iter=0
            do i=1,ngroup
              do j=1,i-1
                do k=1,j-1
                iter=iter+1
                sp3(iter)=sp3(iter)+npop(i)*npop(j)*npop(k)*p
                enddo
              enddo
            enddo
          endif
        enddo
      enddo
!$omp end do
!$omp critical
      p1 = p1 + sp1
      xli = xli + sxli
      if (ngroup.gt.1) p2 = p2 + sp2
      if (ngroup.gt.2) p3 = p3 + sp3
!$omp end critical
      deallocate(npop, sp1, sxli)
      if (ngroup.gt.1) deallocate(sp2)
      if (ngroup.gt.2) deallocate(sp3)
!$omp end parallel
c
      write (lw,11)
      do i=1,ngroup
        write (lw,15) i,p1(i)
      enddo
      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (lw,2) i,j,p2(ipair)
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
            write (lw,3) i,j,k,p3(iter)
            enddo
          enddo
        enddo
      endif
c
!$omp parallel
!$omp& private(npop,sd1,sd2,p,i,j,k,npb,ipair,iter)
      allocate(npop(ngroup))
      if (ngroup.gt.1) then
        allocate(sd1(ngpair))
        sd1 = 0.0d+00
      endif
      if (ngroup.gt.2) then
        allocate(sd2(ntripl))
        sd2 = 0.0d+00
      endif
!$omp do schedule(dynamic)
      do npa=1,npab
        do npb=1,npab
          do i=1,ngroup
            npop(i)=rsrs(npa,i)+rsrs(npb,i)+2*mocore(i)
          enddo
          p=probr(npa,1)*probr(npb,2)
          if (ngroup.gt.1) then
            ipair=0
            do i=1,ngroup
              do j=1,i-1
                ipair=ipair+1
                sd1(ipair)=sd1(ipair)-2*p*(p1(i)-npop(i))*
     &                     (p1(j)-npop(j))
              enddo
            enddo
          endif
          if (ngroup.gt.2) then
            iter=0
            do i=1,ngroup
              do j=1,i-1
                do k=1,j-1
                  iter=iter+1
                  sd2(iter)=sd2(iter)-
     &            2*p*(p1(i)-npop(i))*(p1(j)-npop(j))*(p1(k)-npop(k))
                enddo
              enddo
            enddo
          endif
        enddo
      enddo
!$omp end do
!$omp critical
      if (ngroup.gt.1) d1 = d1 + sd1
      if (ngroup.gt.2) d2 = d2 + sd2
!$omp end critical
      deallocate(npop)
      if (ngroup.gt.1) deallocate(sd1)
      if (ngroup.gt.2) deallocate(sd2)
!$omp end parallel
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
      call timer (4,ipidmain,'_compedfst',-1)
      return
c
c.....Formats.
c
 200  format (/,' # M-BASINS EDFs FOR ALPHA OR BETA ELECTRONS',
     & /,1x,'#',72('-'),/,' # NUMBER OF GROUPS',13x,' = ',I8,/,
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
 100  format (' # ',F22.16,1x,12I6)
 1060 format (1x,'#',72('-'),/,' # ',F22.16,2x,'<-- SUM,',I8,
     & ' PROBABILITIES > ',F16.10,/,' # ',F22.16,2x,
     &   '<-- TOTAL SUM',/,1x,'#',72('-'))
 106  format (1x,'#',72('-'),/,' # ',F22.16,2x,'<-- SUM,',I8,
     & ' PROBABILITIES > ',F16.10,/,' # ',F22.16,2x,
     &   '<-- TOTAL SUM (NON NECESSARILY EQUAL TO 1.0)',
     &   /,1x,'#',72('-'))
 11   format (//1x,'Average populations and localization indices')
 33   format (//1x,'Delocalization indices,',
     &             ' Eq. (28) J. Chem. Phys.  126, 094102 (2007)')
 15   format (1x,'<n(',I3,')>               = ',F16.10)
 2    format (1x,'<n(',I3,') n(',I3,')>        = ',F16.10)
 3    format (1x,'<n(',I3,') n(',I3,') n(',I3,')> = ',F16.10)
 4    format (1x,'delta_(',I3,I3,')         = ',F16.10)
 5    format (1x,'delta_(',I3,I3,I3,')      = ',F16.10)
 6    format (1x,'delta_(',I3,I3,')         = ',F16.10,
     &        2x,'% Localization = ',F16.10)
 111  format (//1x,'Average ALPHA and BETA populations')
 151  format (1x,'<n(',I3,')>_{',a,'} = ',F16.10)
      end
