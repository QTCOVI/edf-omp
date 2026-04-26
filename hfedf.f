c
c-----------------------------------------------------------------------
c
      subroutine hfedf (nproba,nprobb,ngroup,nmo,moval,ival,nel,
     &   malv,mbev,ifilc,sg)

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
      USE         space_for_rsrs

      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
c
      real(kind=8), allocatable,dimension (:,:)   :: am,tpow,probs,oveaa
      real(kind=8), allocatable,dimension (:,:,:) :: ovea
      real(kind=8), allocatable,dimension (:)     :: w,ww
      integer, allocatable,dimension (:)     :: ipvt,indx,ordia,ordib
c
      real(kind=8)      sg(ngroup,nmo,nmo)
      real(kind=8)      dumi,random,deter(2)
      integer*4   idum
      integer     ival(moval)
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
      call semilla (idum)
      idum=-idum
      dumi=random(idum)
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
          call dgeco (am,npmax,nprobs,ipvt,rcond,w)
          if (one + rcond .eq. one) goto 2000
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
            call dgeco (oveaa,moval,malmbe,indx,rcond,ww)
            job=10
            call dgedi (oveaa,moval,malmbe,indx,deter,ww,job)
            probs(n,ii)=deter(1)*tenp**deter(2)
          enddo
          job=0
          call dgesl (am,npmax,nprobs,ipvt,probs(1,ii),job)
        else
          probs(1,ii)=one
        endif
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
          pexact(nn)=pexact(nn)+probs(i,1)*probs(j,2)
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

      deallocate (tpow,am,ipvt,ordia,ordib,w,ovea,oveaa,ww,indx,probs)
      return
      end
