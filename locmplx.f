c
c-----------------------------------------------------------------------
c
      subroutine locmplx (sg,nmo,ngroup,ncent,
     &         nfugrp,ifugrp,pcut,lw,epsloc,okw)
c
c.....ruedmis routine for complex sg() integrals. 
c-----------------------------------------------------------------------
c                         
      include 'implicit.inc'
      include 'constants.inc'
      include 'error.inc'
      include 'stderr.inc'
c
      complex*16, allocatable,dimension (:,:,:)  :: sgn
      real(kind=8), allocatable,dimension (:)      :: xloci
      real(kind=8), allocatable,dimension (:,:)    :: c
      integer,    allocatable,dimension (:)      :: iloc
      integer,    allocatable,dimension (:)      :: mocore,eignz
      integer,    allocatable,dimension (:,:)    :: eleca,resnca,resncb

      parameter   (closetone=0.001D0)
      complex*16  aom1cmplx
      real(kind=8)      aom1
      logical     noloc
      complex*16  sg(ngroup,nmo,nmo)
      integer     nfugrp(ngroup),ifugrp(ncent,ngroup)
c
      call timer (2,iloclx,'_locmplx  ',-1)
c
c.....Allocate arrays
c
      if (.not.allocated(c)) then
        allocate (c(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot allocate c()'
      endif
      if (.not.allocated(xloci)) then
        allocate (xloci(nmo),stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot allocate xloci()'
      endif
      if (.not.allocated(iloc)) then
        allocate (iloc(nmo),stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot allocate iloc()'
      endif

      write (lw,111)
      call cmplxlo (sg,c,xloci,nmo,ngroup,lw,okw)
      do i=1,nmo
        sgmax=0d0
        do k=1,ngroup
          if (dreal(sg(k,i,i)).gt.sgmax) then
            sgmax=dreal(sg(k,i,i))
            kis=k
          endif
        enddo
      enddo
c
c-----Determine the localized and non localized MOs and
c     determine to which center each localized MO belongs
c
      if (.not.allocated(eignz)) then
        allocate (eignz(nmo),stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot allocate eignz()'
      endif
      if (.not.allocated(mocore)) then
        allocate (mocore(ngroup),stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot allocate mocore()'
      endif
      nloc=izero
      nzeig=izero
      mocore(1:ngroup)=izero
      do i=1,nmo
        if (abs(xloci(i)-one).lt.epsloc) then 
          nloc=nloc+1
          iloc(nloc)=i
          smax=-1d0
          do k=1,ngroup
            this=abs(dreal(sg(k,i,i)))
            if (this.gt.smax) then
              smax=this
              kmax=k
            endif
          enddo
          mocore(kmax)=mocore(kmax)+1
          write (lw,112) i,kmax,this
        else
          nzeig=nzeig+1
          eignz(nzeig)=i
        endif
      enddo
      write (lw,*) nloc,' localized MOs,',nzeig,' NON localized MOs'
      write (lw,*) (iloc(i),i=1,nloc)
      write (lw,113) (mocore(i),i=1,ngroup)
c
c-----Reconstruct the Group Overlap Matrix.
c
      if (.not.allocated(sgn)) then
        allocate (sgn(ngroup,nzeig,nzeig),stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot allocate sgn()'
      endif
      do k=1,ngroup-1
        do i=1,nzeig
          i1=eignz(i)
          do j=1,i
            j1=eignz(j)
            sgn(k,i,j)=sg(k,i1,j1)
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
      if (.not.allocated(eleca)) then
        allocate (eleca(2,ngroup),stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot allocate eleca()'
      endif
      eleca(1,1:ngroup)=izero
      eleca(2,1:ngroup)=nzeig
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
        if (ier.ne.0) stop 'locmplx.f: Cannot allocate resncb()'
      endif
c
c.....Computation of resonance structures.
c
      call rnprobs (eleca,resncb,nzeig,ngroup,npab,lw)
      if (nprev.lt.npab) then
        write (stderr,*) 'locmplx.f: NPREV < NPAB returned by rnprobs.f'
        stop
      endif
c
      if (.not.allocated(resnca)) then
        allocate (resnca(npab,ngroup),stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot allocate resnca()'
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
      if (allocated(c)) then
        deallocate (c,stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot deallocate c()'
      endif
      if (allocated(xloci)) then
        deallocate (xloci,stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot deallocate xloci()'
      endif
      if (allocated(iloc)) then
        deallocate (iloc,stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot deallocate iloc()'
      endif
      if (allocated(mocore)) then
        deallocate (mocore,stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot deallocate mocore()'
      endif
      if (allocated(eignz)) then
        deallocate (eignz,stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot deallocate eignz()'
      endif
      if (allocated(sgn)) then
        deallocate (sgn,stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot deallocate sgn()'
      endif
      if (allocated(eleca)) then
        deallocate (eleca,stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot deallocate eleca()'
      endif
      if (allocated(resncb)) then
        deallocate (resncb,stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot deallocate resncb()'
      endif
      if (allocated(resnca)) then
        deallocate (resnca,stat=ier)
        if (ier.ne.0) stop 'locmplx.f: Cannot deallocate resnca()'
      endif
      call timer (4,iloclx,'_locmplx  ',-1)
      return
c
 111  format (' #',/,1x,'# ',75('-'),/,
     &  25x,'ISOPYCNIC LOCALIZATION MODULE')
 112  format (' MO ',I8,' : Localization on fragment ',I2,' = ',F16.10)
 113  format (' # MOs full-localized in each fragment:',20I6)
 766  format (' # NUMBER OF GROUPS = ',I2)
 764  format (' # GROUP ',I2,
     &  ' MinElec and MaxElec (alpha or beta) = ',2I6)
 767  format (1x,'# Number of atoms in the group = ',I3)
 765  format (' # ',10I6)
      end
