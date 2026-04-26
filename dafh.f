
c-----------------------------------------------------------------------
c
      subroutine dafh (nmo,ncent,naom,ngroup,nfugrp,ifugrp,pcut,
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
      real(kind=8), allocatable,dimension (:)     :: seigen
      integer,    allocatable,dimension (:)     :: eignz
      integer,    allocatable,dimension (:,:)   :: eleca,resnca,resncb
      integer,    allocatable,dimension (:)     :: mocore
      complex*16  val1,val2,aom1cmplx
      integer     nzeig,zeig
      integer     nfugrp(ngroup),ifugrp(ncent,ngroup)
      real(kind=8)      aom1
c
      call timer (2,idafh,'_dafh     ',-1)
c
c-----Read AOM
c
      rewind (lu18)
      if (.not.allocated(aom)) then
        allocate (aom(naom,nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot allocate aom()'
      endif
      call readaom (aom,naom,ncent,nmo,lu18,lerr)
      rewind (lu18)
c
c.....Compute complex Group Overlap integrals of all but the last one.
c
      if (.not.allocated(sg)) then
        allocate (sg(ngroup,nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot allocate sg()'
      endif
      sg(1:ngroup,1:nmo,1:nmo)=cmplx(zero,zero)
      do i=1,ngroup-1
        do j=1,nfugrp(i)
          k=ifugrp(j,i)
          sg(i,1:nmo,1:nmo)=sg(i,1:nmo,1:nmo)+aom(k,1:nmo,1:nmo)
        enddo
      enddo
c
c.....Compute DAFH and diagonalize it.
c
      if (.not.allocated(sdiag)) then
        allocate (sdiag(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot allocate sdiag()'
      endif
      if (.not.allocated(bdiag)) then
        allocate (bdiag(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot allocate bdiag()'
      endif
      if (.not.allocated(seigen)) then
        allocate (seigen(nmo),stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot allocate seigen()'
      endif
      if (.not.allocated(wdiag)) then
        allocate (wdiag(nmo+nmo-1),stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot allocate wdiag()'
      endif
c
      do i=1,nmo
        do j=1,i
          sdiag(i,j)=sg(1,i,j)
          sdiag(j,i)=conjg(sdiag(i,j))
        enddo
      enddo
      call zjacobi (sdiag,bdiag,wdiag,seigen,nmo)
      pop=zero
      do i=1,nmo
        pop=pop+seigen(i)
      enddo
      write (lw,112)
      write (lw,222) (seigen(i),i=1,nmo)
      write (lw,113) pop
c
c-----Determine the non-zero elements of seigen().
c
      if (.not.allocated(eignz)) then
        allocate (eignz(nmo),stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot allocate eignz()'
      endif
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
      write (lw,111) epsdafh,zeig
      write (lw,114) epsdafh,nzeig
      write (lw,222) (seigen(eignz(i)),i=1,nzeig)
      write (lw,113) pop
c
c-----Reconstruct the Group Overlap Matrix.
c
      if (.not.allocated(sgn)) then
        allocate (sgn(ngroup,nzeig,nzeig),stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot allocate sgn()'
      endif
      do k=1,ngroup-1
        do i=1,nzeig
          i1=eignz(i)
          do j=1,i
            j1=eignz(j)
            val1=cmplx(zero,zero)
            do l=1,nmo
              do m=1,nmo
                val1=val1+conjg(sdiag(l,i1))*sdiag(m,j1)*sg(k,l,m)
              enddo
            enddo
            sgn(k,i,j)=val1
            sgn(k,j,i)=conjg(sgn(k,i,j))
          enddo
        enddo
      enddo
c
c.....Compute the Group Overlap Matrix for the last group.
c
      do i=1,nzeig
        aom1=0d0
        do k=1,ngroup-1
          aom1=aom1+dreal(sgn(k,i,i))
        enddo
        sgn(ngroup,i,i)=cmplx(one-aom1,zero)
      enddo
      do i=2,nzeig
        do j=1,i-1
          aom1cmplx=cmplx(zero,zero)
          do k=1,ngroup-1
            aom1cmplx=aom1cmplx+sgn(k,i,j)
          enddo
          sgn(ngroup,i,j)=-aom1cmplx
          sgn(ngroup,j,i)=conjg(sgn(ngroup,i,j))
        enddo
      enddo
c
      if (.not.allocated(eleca)) then
        allocate (eleca(2,ngroup),stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot allocate eleca()'
      endif
      eleca(1,1:ngroup)=izero
      eleca(2,1:ngroup)=nzeig
      if (.not.allocated(mocore)) then
        allocate (mocore(ngroup),stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot allocate mocore()'
      endif
      mocore(1:ngroup-1)=izero
      mocore(ngroup)=zeig
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
        if (ier.ne.0) stop 'dafh.f: Cannot allocate resncb()'
      endif
c
c.....Computation of resonance structures.
c
      call rnprobs (eleca,resncb,nzeig,ngroup,npab,lw)
      if (nprev.lt.npab) then
        write (stderr,*) 'dafh.f: NPREV < NPAB returned by rnprobs.f'
        stop
      endif
c
      if (.not.allocated(resnca)) then
        allocate (resnca(npab,ngroup),stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot allocate resnca()'
      endif
      resnca(1:npab,1:ngroup)=resncb(1:npab,1:ngroup)
c
c-----Compute EDF.
c
      call cmplxedf (pcut,npab,nprob,ngroup,nzeig,lw,eleca,sgn,
     &               resnca,mocore)
c
c-----Deallocate arrays.
c
      if (allocated(aom)) then
        deallocate (aom,stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot deallocate aom()'
      endif
      if (allocated(sg)) then
        deallocate (sg,stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot deallocate sg()'
      endif
      if (allocated(sdiag)) then
        deallocate (sdiag,stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot deallocate sdiag()'
      endif
      if (allocated(bdiag)) then
        deallocate (bdiag,stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot deallocate bdiag()'
      endif
      if (allocated(seigen)) then
        deallocate (seigen,stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot deallocate seigen()'
      endif
      if (allocated(wdiag)) then
        deallocate (wdiag,stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot deallocate wdiag()'
      endif
      if (allocated(eignz)) then
        deallocate (eignz,stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot deallocate eignz()'
      endif
      if (allocated(sgn)) then
        deallocate (sgn,stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot deallocate sgn()'
      endif
      if (allocated(eleca)) then
        deallocate (eleca,stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot deallocate eleca()'
      endif
      if (allocated(resncb)) then
        deallocate (resncb,stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot deallocate resncb()'
      endif
      if (allocated(resnca)) then
        deallocate (resnca,stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot deallocate resnca()'
      endif
      if (allocated(mocore)) then
        deallocate (mocore,stat=ier)
        if (ier.ne.0) stop 'dafh.f: Cannot deallocate mocore()'
      endif
c
      call timer (4,idafh,'_dafh     ',-1)
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
