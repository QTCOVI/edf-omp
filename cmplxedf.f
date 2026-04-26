c
c-----------------------------------------------------------------------
c
      subroutine cmplxedf (pcut,npab,nprob,ngroup,nmo,lw,eleca,sg,
     &                     rsrs,mocore)
c
c-----Heavily modified version of 'nbasab.f' to deal with complex AOMs
c
c.......................................................................
c
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
c
      real(kind=8),    allocatable,dimension (:,:)    :: tpow
      complex*16, allocatable,dimension (:,:)   :: am,oveaa
      complex*16, allocatable,dimension (:)     :: probs
      complex*16 determ

      real(kind=8), allocatable,dimension (:,:)   :: probr
      real(kind=8), allocatable,dimension (:)     :: w,ww,probt,pnew
      real(kind=8), allocatable,dimension (:)     :: p1,p1a,p2,d1,p3,d2
      real(kind=8), allocatable,dimension (:)     :: xli
      real(kind=8), allocatable,dimension (:)     :: probord
      integer,   allocatable,dimension (:)        :: ipvt,indx,ind
      integer,   allocatable,dimension (:,:)      :: resnc,resalp
      integer,   allocatable,dimension (:)        :: igtzero
      integer,   allocatable,dimension (:,:)      :: resncord
      integer,   allocatable,dimension (:)        :: npop
      integer,   allocatable,dimension (:)        :: ioprob
c
      complex*16  sg(ngroup,nmo,nmo)
      complex*16  dumaom
      real(kind=8)     dumi,random,deter(2)
      real(kind=8)     rcond
      integer    rsrs(npab,ngroup),lw
      integer    mocore(ngroup)
      integer    eleca(2,ngroup)
      integer*4  idum
      logical    inlist,singular
c
      call timer (2,ipid,'_cmplxedf ',-1)
c
c.....Set the number of electrons two twice the number of MOs
c
      nel=nmo+nmo
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
      do while (singular)
        am(1:npab,1:npab)=(zero,zero)
        tpow=zero
        do i=1,npab
          do k=1,ngroup-1
            tpow(i,k)=two*random(idum)-one
          enddo
        enddo
        do i=1,npab
          do j=1,npab
            aco=one
            do k=1,ngroup-1
              aco=aco*tpow(i,k)**rsrs(j,k)
            enddo
            am(i,j)=cmplx(aco,zero)
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
        call zgeco (oveaa,nmo,nmo,indx,rcond,ww)
        job=10
        call zgedi (oveaa,nmo,nmo,indx,deter,ww,job,determ)
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
        probr(n,1)=dreal(probs(n))
      enddo
      write (lw,200) ngroup,npab

      if (.not.allocated(resalp)) then
        allocate (resalp(ngroup,npab),stat=ier)
        if (ier.ne.0) stop 'cmplxedf.f: Cannot allocate resalp()'
      endif
      if (.not.allocated(igtzero)) then
        allocate (igtzero(npab),stat=ier)
        if (ier.ne.0) stop 'cmplxedf.f: Cannot allocate igtzero()'
      endif
      nnonz=0
      sum=zero
      sumtot=zero
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
      do i=1,npab
        probr(i,2)=probr(i,1)
      enddo
c
c.....Total spin-splitted probabilities 
c
      allocate (probt(npab*npab))
      nprobg=0
      sum=zero
      sumtot=zero
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
      p1a=zero
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
      call sndelta (probr,p1a,p1a,rsrs,rsrs,ngroup,npab,npab,npab,lw)
c
c.....Computes spinless probabilities from spin resolved probabilities.
c     Obtain the number of spinless real space resonant structures.
c
      nprobt=1
      do i=1,ngroup-1
        nprobt=nprobt*(2*eleca(2,i)-2*eleca(1,i)+1)
      enddo
      allocate (pnew(nprobt))
      allocate (resnc(nprobt,ngroup))
      allocate (ind(ngroup))
c
      pnew=zero
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
c.....Determine if probabilities are ordered or not.
c
      write (lw,202) ngroup,np
      npnew=0
      sumtot=zero
      sum=zero
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
      do npa=1,npab
        do npb=1,npab
          p=probr(npa,1)*probr(npb,2)
          do i=1,ngroup
            npop(i)=rsrs(npa,i)+rsrs(npb,i)+2*mocore(i)
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
      enddo


*     nprob=np
*     do np=1,nprob
*       p=pnew(np)
*       do i=1,ngroup
*         npop(i)=resnc(np,i)
*       enddo
*       do i=1,ngroup
*         p1(i)=p1(i)+npop(i)*p
*         xli(i)=xli(i)+p*npop(i)*npop(i)
*       enddo
*       if (ngroup.gt.1) then
*         ipair=0
*         do i=1,ngroup
*           do j=1,i-1
*             ipair=ipair+1
*             p2(ipair)=p2(ipair)+npop(i)*npop(j)*p
*           enddo
*         enddo
*       endif
*       if (ngroup.gt.2) then
*         iter=0
*         do i=1,ngroup
*           do j=1,i-1
*             do k=1,j-1
*             iter=iter+1
*             p3(iter)=p3(iter)+npop(i)*npop(j)*npop(k)*p
*             enddo
*           enddo
*         enddo
*       endif
*     enddo


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
        xli(i)=-(xli(i)-p1(i)*(p1(i)+one))
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
      call timer (4,ipid,'_cmplxedf ',-1)
      return
c
c.....Formats.
c
 200  format (/,' # M-BASINS EDFs (COMPLEX AOM VERSION) ',
     & 'FOR ALPHA OR BETA ELECTRONS',/,1x,'#',72('-'),/,
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
 112  format (1x,'#',/,1x,'# multiple-group delocalization indices',/,
     &        1x,'#')
 113  format (1x,'DELTA (alpha,beta,total) = ',3(1x,F16.10),5x,20(I3))
 1120 format (3(1x,'#',/),1x,'# ',80('+'),/,
     &        1x,'# EXACT CALCULATION OF PROBABILITIES',/,1x,'#')
 206  format (1x,'#',/,1x,'# THERE ARE NO ',A,' ELECTRONS',/,1x,'#')
      end
