c
c-----------------------------------------------------------------------
c
      subroutine xnbasab (probcut,nproba,nprobb,nprob,ngroup,nmo,
     &   moval,ival,mocogrp,nel,nelv,malv,mbev,ifilc,stdout,sg,
     &   resnca,resncb,largwr,doentropy,orderp)

c
c.....Obtains the probabilities for all the real space resonant structures 
c     of the ALPHA and BETA electrons independently by solving two linear 
c     systems. After this, the overall probabilities are obtained. 
c
c.....!!!!! THIS ROUTINE IN ONLY FOR SINGLE DETERMINANT WAVE FUNCTIONS 
c
c.......................................................................
c
      USE         space_for_cidet
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
c
      real(kind=8), allocatable,dimension (:,:):: am,tpow,probs,ovea
      real(kind=8), allocatable,dimension (:)  :: w,ww,probt,pnew
      real(kind=8), allocatable,dimension (:)  :: p1,p1a,p1b,p2,d1,p3,d2
      real(kind=8), allocatable,dimension (:)  :: xli
      real(kind=8), allocatable,dimension (:)  :: probord
      integer, allocatable,dimension (:)       :: ipvt,indx,ind
      integer, allocatable,dimension (:)       :: ordia,ordib      
      integer, allocatable,dimension (:,:)     :: resnc
      integer, allocatable,dimension (:)       :: npop
      integer, allocatable,dimension (:,:)     :: resncord
      integer, allocatable,dimension (:)       :: ioprob
c
      real(kind=8)     sg(ngroup,nmo,nmo)
      real(kind=8)     dumi,random
      integer     resnca(nproba,ngroup),resncb(nprobb,ngroup),stdout
      integer     ig(100)
      integer     mocogrp(ngroup)
      integer*4   idum
      logical     inlist,orderp,doentropy,largwr
      integer     ival(moval)
c
      call timer (2,ipid,'_xnbasab  ',-1)
c
      allocate ( ordia(nmo) ) 
      allocate ( ordib(nmo) )
      mela=0
      melb=0
      read (ifilc,rec=1) (cidet(m),m=0,nel)
      do m=1,nel
        mm=int(cidet(m))
        if (mm.gt.0) then
          mela=mela+1
          ordia(mela)=iabs(mm)
        else
          melb=melb+1
          ordib(melb)=iabs(mm)
        endif
      enddo
c
      call semilla (idum)
      idum=-idum
      dumi=random(idum)
c
      npmax=max(nproba,nprobb)
      allocate ( am(npmax,npmax) )
      allocate ( tpow(npmax,ngroup) )
      allocate ( ipvt(npmax) )
      allocate ( w(npmax) )
      allocate ( ovea(ngroup,moval) )
      allocate ( ww(moval) )
      allocate ( indx(moval) )
      allocate ( probs(npmax,2) )

      write (stdout,1120) 
      do ii=1,2    ! 1 ==> alpha, 2 ==> beta
 2000   am=zero
        if (ii.eq.1) then
           nprobs=nproba
           malmbe=malv
        else
           nprobs=nprobb
           malmbe=mbev
        endif
        if (malmbe.gt.izero) then
c
c.........Construct the first member of the linear system.
c
          tpow=zero
          do i=1,nprobs
            do k=1,ngroup-1
              tpow(i,k)=two*random(idum)-one
            enddo
          enddo
          do i=1,nprobs
            do j=1,nprobs
              aco=one
              do k=1,ngroup-1
                if (ii.eq.1) nn=resnca(j,k)
                if (ii.eq.2) nn=resncb(j,k)
                aco=aco*tpow(i,k)**nn
              enddo
              am(i,j)=aco
            enddo
          enddo
c
          call timer (2,ipid1,'_dgeco    ',-1)
          call dgeco (am,npmax,nprobs,ipvt,rcond,w)
          write (stdout,444) rcond
*         if (one + rcond .eq. one) goto 2000
          if (abs(rcond).lt.1d-50) goto 2000
          call timer (4,ipid1,'_dgeco    ',-1)
c
c.........Overlap matrices
c
          do m=1,malmbe
            if (ii.eq.1) then
              ioma=ordia(ival(m))
            else
              ioma=ordib(ival(m))
            endif
            do igr=1,ngroup
               ovea(igr,m)=sg(igr,ioma,ioma)
            enddo
          enddo
c
          do n=1,nprobs
            det=one
            do m=1,malmbe
              dumi=ovea(ngroup,m)
              do igr=1,ngroup-1
                 dumi=dumi+tpow(n,igr)*ovea(igr,m)
              enddo
              det=det*dumi
            enddo
            probs(n,ii)=det
          enddo
c
c.........Linear System Solver.
c 
          call timer (2,ipid2,'_dgesl    ',-1)
          job=0
          call dgesl (am,npmax,nprobs,ipvt,probs(1,ii),job)
          call timer (4,ipid2,'_dgesl    ',-1)
c
          if (ii.eq.1) write (stdout,200) 'ALPHA',ngroup,'ALPHA',nprobs
          if (ii.eq.2) write (stdout,200) 'BETA ',ngroup,'BETA ',nprobs
          nprobg=0
          sum=zero
          sumtot=zero
          do n=1,nprobs
            sumtot=sumtot+probs(n,ii)
            if (probs(n,ii).ge.probcut) then
              nprobg=nprobg+1
              sum=sum+probs(n,ii)
              if (ii.eq.1) then
               write (stdout,100) probs(n,ii),
     &            (resnca(n,igr)+mocogrp(igr),igr=1,ngroup)
              else
               write (stdout,100) probs(n,ii),
     &            (resncb(n,igr)+mocogrp(igr),igr=1,ngroup)
              endif
            endif
          enddo
          write (stdout,106) sum, nprobg, probcut, sumtot
        else
          if (ii.eq.1) write (stdout,206) 'ALPHA'
          if (ii.eq.2) write (stdout,206) 'BETA'
          probs(1,ii)=one
        endif
      enddo
      deallocate (tpow,am,ipvt,ordia,ordib,w,ovea,ww,indx)
c
c.....Total spin-splitted probabilities 
c
      allocate (probt(nproba*nprobb))
      nprobg=0
      sum=zero
      sumtot=zero
      write (stdout,2001) ngroup,nproba*nprobb
      ij=0
      do i=1,nproba
        do j=1,nprobb
          ij=ij+1
          probt(ij)=probs(i,1)*probs(j,2)
          sumtot=sumtot+probt(ij)
*         if (probt(ij).ge.probcut) then
          if (probs(i,1).gt.probcut .and. probs(j,2).gt.probcut) then
            nprobg=nprobg+1
            sum=sum+probt(ij)
            write (stdout,100) probt(ij),
     &                         (resnca(i,igr)+mocogrp(igr),
     &                          resncb(j,igr)+mocogrp(igr),
     &                          igr=1,ngroup)
          endif
        enddo
      enddo
      write (stdout,1060) sum, nprobg, probcut, sumtot
c
c.....Average population of alpha and beta electrons in each group.
c
      allocate (p1a(ngroup),p1b(ngroup))
      p1a=zero
      do i=1,nproba
        do j=1,ngroup
          p1a(j)=p1a(j)+resnca(i,j)*probs(i,1)
        enddo
      enddo
      p1b=zero
      do i=1,nprobb
        do j=1,ngroup
          p1b(j)=p1b(j)+resncb(i,j)*probs(i,2)
        enddo
      enddo
      write (stdout,111)
      do i=1,ngroup
        write (stdout,151) i,'ALPHA',p1a(i)+mocogrp(i)
      enddo
      do i=1,ngroup
        write (stdout,151) i,'BETA ',p1b(i)+mocogrp(i)
      enddo
c
c.....Obtain multiple-fragment electron population covariances.
c
      call sndelta (probs,p1a,p1b,resnca,resncb,ngroup,
     &    nproba,nprobb,npmax,stdout)
c
c.....Computes spinless probabilities from spin resolved probabilities.
c     Obtain the number of spinless real space resonant structures.
c
      combi=one
      do i=0,nelv-1
        combi=combi*dble(nelv+ngroup-1-i)/dble(nelv-i)
      enddo
      nprobt=int(combi+epsq10)
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
c-----Computation of Mutual Entropy Information.
c
      if (doentropy) call mutent 
     &     (pnew,resnc,mocogrp,nprobt,np,ngroup,nel,largwr,stdout) 
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
        write (stdout,15) i,p1(i)+2*mocogrp(i)
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
      deallocate (p1,xli,npop)
      if (ngroup.gt.1) deallocate (p2,d1)
      if (ngroup.gt.2) deallocate (p3,d2)
      return
c
      close (ifilc,status='delete')
      call timer (4,ipid,'_xnbasab  ',-1)
      return
c
c.....Formats.
c
 1120 format (3(1x,'#',/),1x,'# ',80('+'),/,
     &        1x,'# APPROXIMATE CALCULATION OF PROBABILITIES',/,1x,'#')
 200  format (/,' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION ',
     & 'FOR ',a,' ELECTRONS',/,1x,
     & '# CORE ELECTRONS ARE INCLUDED IN THE RESONANCE STRUCTURES',
     &             /,1x,'#',72('-'),/,
     &  ' # NUMBER OF GROUPS',33x,' = ',I8,/,
     &  ' # TOTAL NUMBER OF PROBABILITIES FOR ',a,' ELECTRONS = ',I8,/,
     &   1x,'#',72('-'),/,
     & ' #     Probability            n1    n2    n3 ...')
 2001 format (/,' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION',
     & ' INCLUDING SPIN', /,1x, 
     & '# CORE ELECTRONS ARE INCLUDED IN THE RESONANCE STRUCTURES',
     &       /,1x,'# ',85('-'),/,
     & ' # NUMBER OF GROUPS               = ',I8,/,
     & ' # TOTAL NUMBER OF PROBABILITIES  = ',I8,/,
     & ' # Gi(a) Gi(b) ARE THE NUMBER OF ALPHA AND BETA ELECTRONS ',
     & ' IN GROUP i',/,1x,'# ',85('-'),/,' #     Probability',11x,
     & 'G1(a) G1(b) G2(a) G2(b) G3(a) G3(b) G4(a) G4(b) G5(a) G5(b)',/,
     & ' # ',85('-'))
 202  format (/,' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION',
     & ' NOT INCLUDING SPIN',/,1x,
     & '# CORE ELECTRONS ARE INCLUDED IN THE RESONANCE STRUCTURES',
     & /,1x,'#',85('-'),/,
     &  ' # NUMBER OF GROUPS               = ',I8,/,
     &  ' # TOTAL NUMBER OF PROBABILITIES  = ',I8,/,
     &   1x,'#',85('-'),/,
     & ' #     Probability            n1    n2    n3 ...')
 444  format (1x,'# RECIPROCAL OF RCOND VALUE IS ',E18.12)
 100  format (' # ',F22.16,1x,12I6)
 106  format (1x,'#',72('-'),/,' #',F22.16,2x,'<--- SUM,',I8,
     & ' PROBABILITIES > ',E16.10,/,' #',F22.16,2x,'<--- TOTAL SUM',/,
     & 1x,'#',72('-'))
 11   format (//1x,'Average populations and localization indices',/,
     &          1x,'(CORE ELECTRONS INCLUDED)')
 33   format (//1x,'Delocalization indices,',
     &             ' Eq. (28) J. Chem. Phys.  126, 094102 (2007)')
 15   format (1x,'<n(',I3,')>               = ',F16.10)
 2    format (1x,'<n(',I3,') n(',I3,')>        = ',F16.10)
 3    format (1x,'<n(',I3,') n(',I3,') n(',I3,')> = ',F16.10)
 4    format (1x,'delta_(',I3,I3,')         = ',F16.10)
 5    format (1x,'delta_(',I3,I3,I3,')      = ',F16.10)
 6    format (1x,'delta_(',I3,I3,')         = ',F16.10,
     &        2x,'% Localization = ',F8.4)
 111  format (//1x,
     & 'Average ALPHA and BETA populations (CORE ELECTRONS INCLUDED)')
 151  format (1x,'<n(',I3,')>_',a,' = ',F16.10)
 112  format (1x,'#',/,1x,'# multiple-group delocalization indices',/,
     &        1x,'#')
 113  format (1x,'DELTA (alpha,beta,total) = ',3(1x,F16.10),5x,20(I3))
 206  format (1x,'#',/,1x,'# THERE ARE NO ',A,' ELECTRONS',/,1x,'#')
 1060 format (1x,'#',85('-'),/,' # ',F22.16,2x,'<- SUM ( ',I8,
     . ' PROBS with P_{a} & P_{b} > ',E14.7,')',/,
     . ' # ',F22.16,2x,'<- TOTAL SUM',/,1x,'#',85('-'))
      end
