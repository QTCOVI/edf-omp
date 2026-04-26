
c-----------------------------------------------------------------------
c
      subroutine gendafh (nmo,ncent,naom,ngroup,nfugrp,ifugrp,pcut,
     &                   epsdafh,lr,lw,lerr,lu18)
c
c.......................................................................
c
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
      include    'stderr.inc'
      complex*16, allocatable,dimension (:,:)   :: sdiag,bdiag
      complex*16, allocatable,dimension (:)     :: wdiag
      complex*16, allocatable,dimension (:,:,:) :: sg,sgn
      complex*16, allocatable,dimension (:,:,:) :: aom
      real(kind=8),     allocatable,dimension (:)     :: seigen
      real(kind=8),     allocatable,dimension (:)     :: xnorm
      integer,    allocatable,dimension (:)     :: eignz
      integer,    allocatable,dimension (:,:)   :: eleca,resnca,resncb
      integer,    allocatable,dimension (:)     :: mocore
      complex*16  val1,val2,aom1cmplx
      integer     nzeig,zeig
      integer     nfugrp(ngroup),ifugrp(ncent,ngroup)
      real(kind=8)      aom1
      logical     notingroup
c
      call timer (2,igdafh,'_gendafh  ',-1)
c
c-----Read AOM
c
      rewind (lu18)
      if (.not.allocated(aom)) then
        allocate (aom(naom,nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot allocate aom()'
      endif
      call readaom (aom,naom,ncent,nmo,lu18,lerr)
      rewind (lu18)
c
c.....Computes SUM(k=1,ngroup) AOM(k,i,i) for i's, and use these values 
c     to renormalize the MOs, as well as to recompute the full AOM.
c
      if (.not.allocated(xnorm)) then
        allocate (xnorm(nmo),stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot allocate xnorm()'
      endif
      ratio=dble(ncent)/dble(naom)
      do i=1,nmo
        aom1=cmplx(zero,zero)
        do k=1,naom
          aom1=aom1+ratio*dreal(aom(k,i,i))
        enddo
        xnorm(i)=one/sqrt(aom1)
      enddo
      do i=1,nmo
        do j=1,nmo
          aom(1:naom,i,j)=aom(1:naom,i,j)*xnorm(i)*xnorm(j)
        enddo
      enddo
      if (allocated(xnorm)) then
        deallocate (xnorm,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate xnorm()'
      endif
c
c.....Determine the atoms in the last group
c
      nfulast=izero
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
c.....Compute complex Group Overlap integrals of all but the last one.
c
      if (.not.allocated(sg)) then
        allocate (sg(ngroup,nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot allocate sg()'
      endif
      sg(1:ngroup,1:nmo,1:nmo)=cmplx(zero,zero)
      do i=1,ngroup-1
        do j=1,nfugrp(i)
          k=ifugrp(j,i)
          sg(i,1:nmo,1:nmo)=sg(i,1:nmo,1:nmo)+aom(k,1:nmo,1:nmo)
        enddo
      enddo
c
c.....Compute complex Group Overlap integrals of the last group.
 
      do m=1,nmo
        aom1=0d0
        do i=1,ngroup-1
          aom1=aom1+dreal(sg(i,m,m))
        enddo
        sg(ngroup,m,m)=cmplx(one-aom1,zero)
      enddo

      do j=2,nmo
        do m=1,j-1
          aom1cmplx=cmplx(zero,zero)
          do i=1,ngroup-1
            aom1cmplx=aom1cmplx+sg(i,j,m)
          enddo
          sg(ngroup,j,m)=-aom1cmplx
          sg(ngroup,m,j)=conjg(sg(ngroup,j,m))
        enddo
      enddo
c
c.....Compute DAFH for each group and diagonalize it.
c
      if (.not.allocated(sdiag)) then
        allocate (sdiag(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot allocate sdiag()'
      endif
      if (.not.allocated(bdiag)) then
        allocate (bdiag(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot allocate bdiag()'
      endif
      if (.not.allocated(seigen)) then
        allocate (seigen(nmo),stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot allocate seigen()'
      endif
      if (.not.allocated(wdiag)) then
        allocate (wdiag(nmo+nmo-1),stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot allocate wdiag()'
      endif
      if (.not.allocated(eignz)) then
        allocate (eignz(nmo),stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot allocate eignz()'
      endif
      if (.not.allocated(eleca)) then
        allocate (eleca(2,ngroup),stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot allocate eleca()'
      endif
      if (.not.allocated(mocore)) then
        allocate (mocore(ngroup),stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot allocate mocore()'
      endif
c
      write (lw,110)
      do m=1,ngroup
        do i=1,nmo
          do j=1,i
            sdiag(i,j)=sg(m,i,j)
            sdiag(j,i)=conjg(sdiag(i,j))
          enddo
        enddo
        call zjacobi (sdiag,bdiag,wdiag,seigen,nmo)
        pop=zero
        do i=1,nmo
          pop=pop+seigen(i)
        enddo
*       write (lw,112) m
*       write (lw,222) (seigen(i),i=1,nmo)
*       write (lw,113) pop
c
c-------Determine the non-zero elements of seigen().
c
        nzeig=izero
        pop=zero
        do i=1,nmo
          if (seigen(i).gt.epsdafh) then
            nzeig=nzeig+1
            pop=pop+seigen(i)
            eignz(nzeig)=i
          endif
        enddo
        zeig=nmo-nzeig
*       write (lw,111) epsdafh,zeig
        write (lw,114) m,epsdafh,nzeig
        write (lw,222) (seigen(eignz(i)),i=1,nzeig)
        write (lw,113) pop
        eleca(1,m)=izero
        eleca(2,m)=nzeig
*       mocore(m)=zeig
        mocore(m)=izero
      enddo
c
c.....Recompute eleca(1,ngroup)
c
      ielecmax=zero
      do m=1,ngroup-1
        ielecmax=ielecmax+eleca(2,m)
      enddo
      eleca(1,ngroup)=nmo-ielecmax

      write (lw,766) ngroup
      do i=1,ngroup
        ielmin=eleca(1,i)+mocore(i)
        ielmax=eleca(2,i)+mocore(i)
        write (lw,764) i,ielmin,ielmax
        write (lw,767) nfugrp(i)
        write (lw,765) (ifugrp(j,i),j=1,nfugrp(i))
      enddo
c
c-----Compute number of probabilities.
c
      npab=1
      do i=1,ngroup-1
        npab=npab*(eleca(2,i)-eleca(1,i)+1)
      enddo
      nprev=npab
c
      if (.not.allocated(resncb)) then
        allocate (resncb(nprev,ngroup),stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot allocate resncb()'
      endif
c
c.....Computation of resonance structures.
c
      call rnprobs (eleca,resncb,nzeig,ngroup,npab,lw)
      if (nprev.lt.npab) then
        write (stderr,*) 'gendafh.f: NPREV < NPAB returned by rnprobs.f'
        stop
      endif
c
      if (.not.allocated(resnca)) then
        allocate (resnca(npab,ngroup),stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot allocate resnca()'
      endif
      resnca(1:npab,1:ngroup)=resncb(1:npab,1:ngroup)
c
c-----Compute EDF.
c
      call cmplxedf (pcut,npab,nprob,ngroup,nzeig,lw,eleca,sg,
     &               resnca,mocore)
c
c-----Deallocate arrays.
c
      if (allocated(aom)) then
        deallocate (aom,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate aom()'
      endif
      if (allocated(sg)) then
        deallocate (sg,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate sg()'
      endif
      if (allocated(sdiag)) then
        deallocate (sdiag,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate sdiag()'
      endif
      if (allocated(bdiag)) then
        deallocate (bdiag,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate bdiag()'
      endif
      if (allocated(seigen)) then
        deallocate (seigen,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate seigen()'
      endif
      if (allocated(wdiag)) then
        deallocate (wdiag,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate wdiag()'
      endif
      if (allocated(eignz)) then
        deallocate (eignz,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate eignz()'
      endif
      if (allocated(sgn)) then
        deallocate (sgn,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate sgn()'
      endif
      if (allocated(eleca)) then
        deallocate (eleca,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate eleca()'
      endif
      if (allocated(resncb)) then
        deallocate (resncb,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate resncb()'
      endif
      if (allocated(resnca)) then
        deallocate (resnca,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate resnca()'
      endif
      if (allocated(mocore)) then
        deallocate (mocore,stat=ier)
        if (ier.ne.0) stop 'gendafh.f: Cannot deallocate mocore()'
      endif
c
      call timer (4,igdafh,'_gendafh  ',-1)
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
