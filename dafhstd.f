
c-----------------------------------------------------------------------
c
      subroutine dafhstd (nmo,ncent,ngroup,nfugrp,ifugrp,pcut,
     &                   epsdafh,lr,lw,ler,lu18)
c
c.......................................................................
c
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
      real   (kind=8), allocatable,dimension (:,:)   :: sdiag,bdiag
      real   (kind=8), allocatable,dimension (:)     :: wdiag
      real   (kind=8), allocatable,dimension (:,:,:) :: sg,sgn
      real   (kind=8), allocatable,dimension (:,:,:) :: aom
      real   (kind=8), allocatable,dimension (:)     :: seigen
      integer(kind=4), allocatable,dimension (:)     :: eignz
      integer(kind=4), allocatable,dimension (:,:)   :: mimael
      integer(kind=4), allocatable,dimension (:,:)   :: resnc
      integer(kind=4), allocatable,dimension (:)     :: mocore
      real   (kind=8)  val1,val2,aom1
      integer(kind=4)  nzeig,zeig
      integer(kind=4)  nfugrp(ngroup),ifugrp(ncent,ngroup)
c
      call timer (2,idafhstd,'_dafhstd  ',-1)
c
c-----Read AOM
c
      rewind (lu18)
      allocate (aom(ncent,nmo,nmo))
      call readaomstd (aom,nfugrp,ifugrp,ncent,ngroup,nmo,naom,lu18,ler)
      rewind (lu18)
c
c.....Group Overlap integrals for all groups but the last one.
c
      allocate (sg(ngroup,nmo,nmo))
      sg(1:ngroup,1:nmo,1:nmo)=0.0d0+00
      do i=1,ngroup-1
        do j=1,nfugrp(i)
          k=ifugrp(j,i)
          sg(i,:,:)=sg(i,:,:)+aom(k,:,:)
        enddo
      enddo
c
c.....Compute DAFH and diagonalize it.
c
      allocate (sdiag(nmo,nmo))
      allocate (bdiag(nmo,nmo))
      allocate (seigen(nmo))
      allocate (wdiag(nmo+nmo-1))
c
      do i=1,nmo
        do j=1,i
          sdiag(i,j)=sg(1,i,j)
          sdiag(j,i)=sdiag(i,j)
        enddo
      enddo
      call jacobi (sdiag,bdiag,wdiag,seigen,nmo)
      pop=0.0d0+00
      do i=1,nmo
        pop=pop+seigen(i)
      enddo
      write (lw,112)
      write (lw,222) (seigen(i),i=1,nmo)
      write (lw,113) pop
c
c-----Determine the non-zero elements of seigen().
c
      allocate (eignz(nmo))
      nzeig=0
      pop=0.0d0+00
      do i=1,nmo
        if (seigen(i).gt.epsdafh) then
          nzeig=nzeig+1
          pop=pop+seigen(i)
          eignz(nzeig)=i
        endif
      enddo
      zeig=nmo-nzeig
      write (lw,111) epsdafh,zeig
      write (lw,114) epsdafh,nzeig
      write (lw,222) (seigen(eignz(i)),i=1,nzeig)
      write (lw,113) pop
c
c-----Reconstruct the Group Overlap Matrix.
c
      allocate (sgn(ngroup,nzeig,nzeig))
      do k=1,ngroup-1
        do i=1,nzeig
          i1=eignz(i)
          do j=1,i
            j1=eignz(j)
            val1=0.0d0+00
            do l=1,nmo
              do m=1,nmo
                val1=val1+sdiag(l,i1)*sdiag(m,j1)*sg(k,l,m)
              enddo
            enddo
            sgn(k,i,j)=val1
            sgn(k,j,i)=sgn(k,i,j)
          enddo
        enddo
      enddo
c
c.....Compute the Group Overlap Matrix for the last group.
c
      do i=1,nzeig
        aom1=0.0d+00
        do k=1,ngroup-1
          aom1=aom1+sgn(k,i,i)
        enddo
        sgn(ngroup,i,i)=1.0d+00-aom1
      enddo
      do i=2,nzeig
        do j=1,i-1
          aom1=0.0d0+00
          do k=1,ngroup-1
            aom1=aom1+sgn(k,i,j)
          enddo
          sgn(ngroup,i,j)=-aom1
          sgn(ngroup,j,i)=sgn(ngroup,i,j)
        enddo
      enddo
c
      allocate (mimael(2,ngroup))
      mimael(1,1:ngroup)=0
      mimael(2,1:ngroup)=nzeig
      allocate (mocore(ngroup))
      mocore(1:ngroup-1)=0
      mocore(ngroup)=zeig
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
      allocate (resnc(nprev,ngroup))
c
c.....Computation of resonance structures.
c
      call rnprobs (mimael,resnc,nzeig,ngroup,npab,lw)
      if (nprev.lt.npab) then
        write (ler,*) 'dafhstd.f: NPREV < NPAB returned by rnprobs.f'
        stop
      endif
c
c-----Compute EDF.
c
      call compedfstd (pcut,npab,nprev,ngroup,nzeig,lw,mimael,sgn,
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
      call timer (4,idafhstd,'_dafhstd  ',-1)
      return
c
 222  format (4(2x,F16.10))
 111  format (' # NUMBER OF EIGENVALUES < ',E16.10,3x,' = ',I5)
 766  format (/,' #',/,' # NUMBER OF GROUPS = ',I2)
 764  format (' # GROUP ',I2,
     &  ' MinElec and MaxElec (alpha or beta) = ',2I6)
 765  format (' # ',10I6)
 112  format (/,' # ',70('+'),/,' #',27x,'DAFH MODULE',/,
     & ' # ',70('+'),/,' # Ordered DAFH Eigenvalues')
 113  format (' # SUM = ',F16.10)
 114  format (' # DAFH Eigenvalues > ',E16.10,3x,' = ',I5)
 767  format (1x,'# Number of atoms in the group = ',I3)
      end
