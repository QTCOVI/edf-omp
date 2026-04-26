
c-----------------------------------------------------------------------
c
      subroutine gendafhstd (nmo,ncent,naom,ngroup,nfugrp,ifugrp,pcut,
     &                   epsdafh,lr,lw,ler,lu18)
c
c.......................................................................
c
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
      real    (kind=8), allocatable,dimension (:,:)   :: sdiag,bdiag
      real    (kind=8), allocatable,dimension (:)     :: wdiag
      real    (kind=8), allocatable,dimension (:,:,:) :: sg,sgn
      real    (kind=8), allocatable,dimension (:,:,:) :: aom
      real    (kind=8), allocatable,dimension (:)     :: seigen
      real    (kind=8), allocatable,dimension (:)     :: xnorm
      real    (kind=8)  val1,val2,aom1
      integer(kind=4),  allocatable,dimension (:)     :: eignz
      integer(kind=4),  allocatable,dimension (:,:)   :: mimael
      integer(kind=4),  allocatable,dimension (:,:)   :: resnc
      integer(kind=4),  allocatable,dimension (:)     :: mocore
      integer(kind=4)   nzeig,zeig
      integer(kind=4)   nfugrp(ngroup),ifugrp(ncent,ngroup)
      logical           notingroup
c
      call timer (2,igdafh,'_gendafhst',-1)
c
c-----Read AOM
c
      rewind (lu18)
      allocate (aom(ncent,nmo,nmo))
      call readaomstd (aom,nfugrp,ifugrp,ncent,ngroup,nmo,naom,lu18,ler)
      rewind (lu18)
c
c.....Computes SUM(k=1,ngroup) AOM(k,i,i) for i's, and use these values 
c     to renormalize the MOs, as well as to recompute the full AOM.
c
      allocate (xnorm(nmo))
      ratio=dble(ncent)/dble(naom)
      do i=1,nmo
        aom1=0.0d+00
        do k=1,naom
          aom1=aom1+ratio*aom(k,i,i)
        enddo
        xnorm(i)=1.0d+00/sqrt(aom1)
      enddo
      do i=1,nmo
        do j=1,nmo
          aom(1:naom,i,j)=aom(1:naom,i,j)*xnorm(i)*xnorm(j)
        enddo
      enddo
      deallocate (xnorm)
c
c.....Determine the atoms in the last group
c
      nfulast=0
      do l=1,ncent
        notingroup=.true.
        do i=1,ngroup-1
          do j=1,nfugrp(i)
            if (l.eq.ifugrp(j,i)) notingroup=.false.
          enddo
        enddo
        if (notingroup) then
          nfulast=nfulast+1
          ifugrp(nfulast,ngroup)=l
        endif
      enddo
      nfugrp(ngroup)=nfulast
c
c.....Compute Group Overlap integrals of all groups but the last one.
c
      allocate (sg(ngroup,nmo,nmo))
      sg(1:ngroup,1:nmo,1:nmo)=0.0d+00
      do i=1,ngroup-1
        do j=1,nfugrp(i)
          k=ifugrp(j,i)
          sg(i,:,:)=sg(i,:,:)+aom(k,:,:)
        enddo
      enddo
c
c.....Compute Group Overlap integrals of the last group.
 
      do m=1,nmo
        aom1=0d0
        do i=1,ngroup-1
          aom1=aom1+sg(i,m,m)
        enddo
        sg(ngroup,m,m)=1.0d+00-aom1
      enddo

      do j=2,nmo
        do m=1,j-1
          aom1=0.0d+00
          do i=1,ngroup-1
            aom1=aom1+sg(i,j,m)
          enddo
          sg(ngroup,j,m)=-aom1
          sg(ngroup,m,j)=sg(ngroup,j,m)
        enddo
      enddo
c
c.....Compute DAFH for each group and diagonalize it.
c
      allocate (sdiag(nmo,nmo))
      allocate (bdiag(nmo,nmo))
      allocate (seigen(nmo))
      allocate (wdiag(nmo+nmo-1))
      allocate (eignz(nmo))
      allocate (mimael(2,ngroup))
      allocate (mocore(ngroup))
c
      write (lw,110)
      do m=1,ngroup
        do i=1,nmo
          do j=1,i
            sdiag(i,j)=sg(m,i,j)
            sdiag(j,i)=sdiag(i,j)
          enddo
        enddo
        call jacobi (sdiag,bdiag,wdiag,seigen,nmo)
        pop=0.0d+00
        do i=1,nmo
          pop=pop+seigen(i)
        enddo
c
c-------Determine the non-zero elements of seigen().
c
        nzeig=0
        pop=0.0d+00
        do i=1,nmo
          if (seigen(i).gt.epsdafh) then
            nzeig=nzeig+1
            pop=pop+seigen(i)
            eignz(nzeig)=i
          endif
        enddo
        zeig=nmo-nzeig
        write (lw,114) m,epsdafh,nzeig
        write (lw,222) (seigen(eignz(i)),i=1,nzeig)
        write (lw,113) pop
        mimael(1,m)=0
        mimael(2,m)=nzeig
        mocore(m)=0
      enddo
c
c.....Recompute mimael(1,ngroup)
c
      ielecmax=0.0d+00
      do m=1,ngroup-1
        ielecmax=ielecmax+mimael(2,m)
      enddo
      mimael(1,ngroup)=nmo-ielecmax

      write (lw,766) ngroup
      do i=1,ngroup
        ielmin=mimael(1,i)+mocore(i)
        ielmax=mimael(2,i)+mocore(i)
        write (lw,764) i,ielmin,ielmax
        write (lw,767) nfugrp(i)
        write (lw,765) (ifugrp(j,i),j=1,nfugrp(i))
      enddo
c
c-----Compute number of probabilities.
c
      npab=1
      do i=1,ngroup-1
        npab=npab*(mimael(2,i)-mimael(1,i)+1)
      enddo
      nprev=npab
c
      allocate (resnc(nprev,ngroup))
c
c.....Computation of resonance structures.
c
      call rnprobs (mimael,resnc,nzeig,ngroup,npab,lw)
      if (nprev.lt.npab) then
        write (ler,*) 'gendafhstd.f: NPREV<NPAB returned by rnprobs'
        stop
      endif
c
c-----Compute EDF.
c
      call compedfstd (pcut,npab,nprev,ngroup,nzeig,lw,mimael,sg,
     &                resnc,mocore)
c
c-----Deallocate arrays.
c
      deallocate (aom)
      deallocate (sg)
      deallocate (sdiag)
      deallocate (bdiag)
      deallocate (seigen)
      deallocate (wdiag)
      deallocate (eignz)
      deallocate (sgn)
      deallocate (mimael)
      deallocate (resnc)
      deallocate (mocore)
      call timer (4,igdafh,'_gendafhst',-1)
      return
c
 222  format (4(2x,F16.10))
 111  format (' # NUMBER OF EIGENVALUES < ',E16.10,3x,' = ',I5)
 766  format (/,' #',/,' # NUMBER OF GROUPS = ',I2)
 764  format (' # GROUP ',I2,
     &  ' MinElec and MaxElec (alpha or beta) = ',2I6)
 767  format (1x,'# Number of atoms in the group = ',I3)
 765  format (' # ',10I6)
 110  format (/,' # ',70('+'),/,' #',27x,'DAFH MODULE',/,' # ',70('+'))
 112  format (' # Group ',I2,5x,'Ordered DAFH Eigenvalues')
 113  format (' # SUM = ',F16.10)
 114  format (/,' # GROUP',I2,
     &          '   DAFH Eigenvalues > ',E16.10,3x,' = ',I5)
      end
