c
c-----------------------------------------------------------------------
c
      subroutine nbasabw (probcut,nproba,nprobb,ngroup,nmo,
     &   short,moval,ival,mocogrp,nel,nelv,malv,mbev,ifilc,stdout,sg,
     &   resnca,resncb,doentropy,orderp,pnew,nprobt)

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
      real*16, allocatable,dimension (:,:)   :: am,tpow,probs,oveaa
      real*16, allocatable,dimension (:,:,:) :: ovea
      real*16, allocatable,dimension (:)     :: w,ww,probt
      real*16, allocatable,dimension (:)     :: p1,p1a,p1b,p2,d1,p3,d2
      real*16, allocatable,dimension (:)     :: xli
      real*16, allocatable,dimension (:)     :: probord
      integer, allocatable,dimension (:)     :: ipvt,indx,ind
      integer, allocatable,dimension (:)     :: ordia,ordib      
      integer, allocatable,dimension (:,:)   :: resnc
      integer, allocatable,dimension (:,:)   :: resncord
      integer, allocatable,dimension (:)     :: npop
      integer, allocatable,dimension (:)     :: ioprob
c
      real*16     sg(ngroup,nmo,nmo)
      real*16     pnew(nprobt)
      real*16     dumi,random,deter(2)
      integer     resnca(nproba,ngroup),resncb(nprobb,ngroup),stdout
      integer     ig(100)
      integer     mocogrp(ngroup)
      integer*4   idum
      logical     inlist,doentropy,orderp,short
      integer     ival(moval)
c
      call timer (2,ipid,'_nbasabw  ',-1)
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
      allocate ( ovea(ngroup,moval,moval) )
      allocate ( oveaa(moval,moval) )
      allocate ( ww(moval) )
      allocate ( indx(moval) )
      allocate ( probs(npmax,2) )
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
          if (one + rcond .eq. one) goto 2000
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
            do k=1,malmbe
              if (ii.eq.1) then
                ioka=ordia(ival(k))
              else
                ioka=ordib(ival(k))
              endif
              do igr=1,ngroup
                 ovea(igr,m,k)=sg(igr,ioma,ioka)
              enddo
            enddo
          enddo
c
          call timer (2,ipid2,'_determ   ',-1)
c
          do n=1,nprobs
            do m=1,malmbe
              do k=1,malmbe
                dumi=ovea(ngroup,m,k)
                do igr=1,ngroup-1
                   dumi=dumi+tpow(n,igr)*ovea(igr,m,k)
                enddo
                oveaa(m,k)=dumi
              enddo
            enddo
c
c...........Determinants calculation using Netlib DGEDI routine.
c
            call dgeco (oveaa,moval,malmbe,indx,rcond,ww)
            job=10
            call dgedi (oveaa,moval,malmbe,indx,deter,ww,job)
            probs(n,ii)=deter(1)*tenp**deter(2)
          enddo
c
          call timer (4,ipid2,'_determ   ',-1)
c
c
c.........Linear System Solver: Netlib DGECO routine.
c 
          call timer (2,ipid3,'_dgesl    ',-1)
          job=0
          call dgesl (am,npmax,nprobs,ipvt,probs(1,ii),job)
          call timer (4,ipid3,'_dgesl    ',-1)
        else
          probs(1,ii)=one
        endif
      enddo
      deallocate (tpow,am,ipvt,ordia,ordib,w,ovea,oveaa,ww,indx)
c
c.....Total spin-splitted probabilities 
c
      allocate (probt(nproba*nprobb))
      ij=0
      do i=1,nproba
        do j=1,nprobb
          ij=ij+1
          probt(ij)=probs(i,1)*probs(j,2)
        enddo
      enddo
      deallocate (probs)
c
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
      if (.not.short) write (stdout,202) ngroup,np
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
              resncord(npronew,igr) = resnc(n,igr)+2*mocogrp(igr)
            enddo
          endif
        enddo
        call qqsort (probord,ioprob,1,npnew,nstack)
        do n=npnew,1,-1
          m=ioprob(n)
          if (.not.short) 
     &    write (stdout,100) probord(m),(resncord(m,igr),igr=1,ngroup)
        enddo
        deallocate (ioprob,resncord,probord)
      else
c
c.......Probabilities are not ordered.
c
        if (.not.short) then
          do n=1,np
            if (pnew(n).ge.probcut) write (stdout,100) pnew(n),
     &         (resnc(n,igr)+2*mocogrp(igr),igr=1,ngroup)
          enddo
        endif
      endif
      if (.not.short) write (stdout,106) sum, npnew, probcut, sumtot
      deallocate (resnc)
      call timer (4,ipid,'_nbasabw  ',-1)
      return
c
c.....Formats.
c
 202  format (' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION',
     & ' NOT INCLUDING SPIN',
     & /,1x,'#',72('-'),/,
     &  ' # NUMBER OF GROUPS               = ',I8,/,
     &  ' # TOTAL NUMBER OF PROBABILITIES  = ',I8,/,
     &   1x,'#',72('-'),/,
     & ' #     Probability            n1    n2    n3 ...')
 100  format (' # ',F22.16,1x,20I6)
 106  format (1x,'#',72('-'),/,' #',F22.16,2x,'<-- SUM,',I8,
     & ' PROBABILITIES > ',E16.10,/,' #',F22.16,2x,'<--- TOTAL SUM',/,
     & 1x,'#',72('-'))
      end
