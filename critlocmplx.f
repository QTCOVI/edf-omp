c
c-----------------------------------------------------------------------
c
      subroutine critlocmplx (sga,sgb,nmo,ngroup,ntotattr,nattr,nspin,
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
      complex*16  sga(ngroup,nmo,nmo)
      complex*16  sgb(ngroup,nmo,nmo)
      complex*16,      allocatable,dimension (:,:,:) :: sgna,sgnb
      complex*16,      allocatable,dimension (:,:,:) :: sg,sgn
      real   (kind=8), allocatable,dimension (:)     :: xloci
      real   (kind=8), allocatable,dimension (:,:)   :: c
      integer(kind=4), allocatable,dimension (:)     :: iloc
      integer(kind=4), allocatable,dimension (:)     :: mocore,eignz
      integer(kind=4), allocatable,dimension (:,:)   :: eleca,elecb
      integer(kind=4), allocatable,dimension (:,:)   :: resnca,resncb

      parameter    (closetone = 0.001D0)
      complex*16   aom1cmplx
      real(kind=8) aom1
      logical      noloc
      integer(kind=4)  nfugrp(ngroup),ifugrp(ntotattr,ngroup)
      character(len=5) spin(2)
c
      call timer (2,iloclxcrit,'_critlocmp',-1)
c
c.....Allocate arrays
c
      if (.not.allocated(c)) then
        allocate (c(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate c()'
      endif
      if (.not.allocated(xloci)) then
        allocate (xloci(nmo),stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate xloci()'
      endif
      if (.not.allocated(iloc)) then
        allocate (iloc(nmo),stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate iloc()'
      endif
      if (.not.allocated(sg)) then
        allocate (sg(ngroup,nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate sg()'
      endif
      if (.not.allocated(sgn)) then
        allocate (sgn(ngroup,nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate sgn()'
      endif
      if (.not.allocated(eignz)) then
        allocate (eignz(nmo),stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate eignz()'
      endif
      if (.not.allocated(mocore)) then
        allocate (mocore(ngroup),stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate mocore()'
      endif
      if (.not.allocated(eleca)) then
        allocate (eleca(2,ngroup),stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate eleca()'
      endif
      if (.not.allocated(elecb)) then
        allocate (elecb(2,ngroup),stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate elecb()'
      endif
c
      spin(1) = 'ALPHA'
      spin(2) = 'BETA '
      ispin = 1
 22   continue
      if (ispin.eq.1) sg = sga
      if (ispin.eq.2) sg = sgb

      write (lw,111)
      call cmplxlo (sg,c,xloci,nmo,ngroup,lw,okw)
      do i=1,nmo
        sgmax = 0.0d+00
        do k=1,ngroup
          if (dreal(sg(k,i,i)).gt.sgmax) sgmax = dreal(sg(k,i,i))
        enddo
      enddo
c
c-----Determine the localized and non localized MOs and
c     determine to which center each localized MO belongs
c
      nloc = 0
      nzeig = 0
      mocore(1:ngroup)=0
      do i=1,nmo
        if (abs(xloci(i)-1.0d+00).lt.epsloc) then 
          nloc = nloc+1
          iloc(nloc) = i
          smax = -1.0d+00
          do k=1,ngroup
            this = abs(dreal(sg(k,i,i)))
            if (this.gt.smax) then
              smax = this
              kmax = k
            endif
          enddo
          mocore(kmax) = mocore(kmax)+1
          write (lw,112) i,kmax,this
        else
          nzeig = nzeig+1
          eignz(nzeig) = i
        endif
      enddo
      write (lw,33) nloc,nzeig
      write (lw,34) (iloc(i),i=1,nloc)
      write (lw,113) (mocore(i),i=1,ngroup)
c
c-----Reconstruct the Group Overlap Matrix.
c
      if (.not.allocated(sgn)) then
        allocate (sgn(ngroup,nzeig,nzeig),stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate sgn()'
      endif
      do k=1,ngroup-1
        do i=1,nzeig
          i1 = eignz(i)
          do j=1,i
            j1 = eignz(j)
            sgn(k,i,j) = sg(k,i1,j1)
            sgn(k,j,i) = conjg(sgn(k,i,j))
          enddo
        enddo
      enddo
c
c.....Compute the Group Overlap Matrix for the last group.
c
      do i=1,nzeig
        aom1 = 0.0d+00
        do k=1,ngroup-1
          aom1 = aom1+dreal(sgn(k,i,i))
        enddo
        sgn(ngroup,i,i) = cmplx(1.0d+00-aom1,0.0d+00)
      enddo
      do i=2,nzeig
        do j=1,i-1
          aom1cmplx = cmplx(0.0d+00,0.0d+00)
          do k=1,ngroup-1
            aom1cmplx = aom1cmplx+sgn(k,i,j)
          enddo
          sgn(ngroup,i,j) = -aom1cmplx
          sgn(ngroup,j,i) = conjg(sgn(ngroup,i,j))
        enddo
      enddo
      if (ispin.eq.1) then
        eleca(1,1:ngroup) = 0
        eleca(2,1:ngroup) = nzeig
        if (.not.allocated(sgna)) then
          allocate (sgna(ngroup,nzeig,nzeig),stat=ier)
          if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate sgna()'
        endif
        sgna(:,1:nzeig,1:nzeig) = sgn(:,1:nzeig,1:nzeig)
        write (lw,766) ngroup
        do i=1,ngroup
          ielmin = eleca(1,i)+mocore(i)
          ielmax = eleca(2,i)+mocore(i)
          write (lw,764) i,trim(spin(ispin)),ielmin,ielmax
          write (lw,767) nfugrp(i)
          write (lw,765) (ifugrp(j,i),j=1,nfugrp(i))
        enddo
      else
        elecb(1,1:ngroup) = 0
        elecb(2,1:ngroup) = nzeig
        if (.not.allocated(sgnb)) then
          allocate (sgnb(ngroup,nzeig,nzeig),stat=ier)
          if (ier.ne.0) stop 'dcritlocmplx.f: Cannot allocate sgnb()'
        endif
        sgnb(:,1:nzeig,1:nzeig) = sgn(:,1:nzeig,1:nzeig)
        write (lw,766) ngroup
        do i=1,ngroup
          ielmin = elecb(1,i)+mocore(i)
          ielmax = elecb(2,i)+mocore(i)
          write (lw,764) i,trim(spin(ispin)),ielmin,ielmax
          write (lw,767) nfugrp(i)
          write (lw,765) (ifugrp(j,i),j=1,nfugrp(i))
        enddo
      endif
      if (allocated(sgn)) then
        deallocate (sgn,stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot deallocate sgn()'
      endif
c
c-----Compute number of probabilities.
c
      if (ispin.eq.1) then
        npais=1
        do i=1,ngroup-1
          npais=npais*(eleca(2,i)-eleca(1,i)+1)
        enddo
        nprev=npais
c
        if (.not.allocated(resnca)) then
          allocate (resnca(nprev,ngroup),stat=ier)
          if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate resnca()'
        endif
c
c.......Computation of ALPHA resonance structures.
c
        call rnprobs (eleca,resnca,nzeig,ngroup,npais,lw)
        if (nprev.lt.npais) then
          write (stderr,*) 'critlocmplx.f: rnprobs returns NPREV<NPAIS'
          stop
        endif
      else
        npbis=1
        do i=1,ngroup-1
          npbis=npbis*(elecb(2,i)-elecb(1,i)+1)
        enddo
        nprev=npbis
c
        if (.not.allocated(resncb)) then
          allocate (resncb(nprev,ngroup),stat=ier)
          if (ier.ne.0) stop 'critlocmplx.f: Cannot allocate resncb()'
        endif
c
c.......Computation of BETA  resonance structures.
c
        call rnprobs (elecb,resncb,nzeig,ngroup,npbis,lw)
        if (nprev.lt.npbis) then
          write (stderr,*) 'critlocmplx.f: rnprobs returns NPREV<NPBIS'
          stop
        endif
      endif
      if (ispin.eq.1.and.nspin.eq.2) then
        ispin = 2
        goto 22
      else
        elecb(1,1:ngroup) = eleca(1,1:ngroup)
        elecb(2,1:ngroup) = eleca(2,1:ngroup)
        npbis = npais
        if (.not.allocated(resncb)) then
          allocate (resncb(nprev,ngroup),stat=ier)
          if (ier.ne.0) stop 'dnocrit.f: Cannot allocate resncb()'
        endif
        resncb(npbis,ngroup) = resnca(npais,ngroup)
        if (.not.allocated(sgnb)) then
          allocate (sgnb(ngroup,nzeig,nzeig),stat=ier)
          if (ier.ne.0) stop 'dnocrit.f: Cannot allocate sgnb()'
        endif
        sgnb(:,1:nzeig,1:nzeig) = sgna(:,1:nzeig,1:nzeig)
      endif
c
c-----Compute EDF.
c
      call xcritcplx (pcut,npais,npbis,ngroup,nmo,nmo,lw,
     &     eleca,elecb,sgna,sgnb,resnca,resncb,mocore)
c
c-----Deallocate arrays.
c
      if (allocated(c)) then
        deallocate (c,stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot deallocate c()'
      endif
      if (allocated(xloci)) then
        deallocate (xloci,stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot deallocate xloci()'
      endif
      if (allocated(iloc)) then
        deallocate (iloc,stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot deallocate iloc()'
      endif
      if (allocated(mocore)) then
        deallocate (mocore,stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot deallocate mocore()'
      endif
      if (allocated(sg)) then
        deallocate (sg,stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot deallocate sg()'
      endif
      if (allocated(sgn)) then
        deallocate (sgn,stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot deallocate sgn()'
      endif
      if (allocated(eignz)) then
        deallocate (eignz,stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot deallocate eignz()'
      endif
      if (allocated(sgn)) then
        deallocate (sgn,stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot deallocate sgn()'
      endif
      if (allocated(eleca)) then
        deallocate (eleca,stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot deallocate eleca()'
      endif
      if (allocated(resncb)) then
        deallocate (resncb,stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot deallocate resncb()'
      endif
      if (allocated(resnca)) then
        deallocate (resnca,stat=ier)
        if (ier.ne.0) stop 'critlocmplx.f: Cannot deallocate resnca()'
      endif
      call timer (4,iloclxcrit,'_critlocmp',-1)
      return
c
 111  format (' #',/,1x,'# ',75('-'),/,
     &  25x,'ISOPYCNIC LOCALIZATION MODULE')
 112  format (' MO ',I8,' : Localization on fragment ',I2,' = ',F17.10)
 113  format (' # MOs full-localized in each fragment:',/,
     &  1000(' # ',20(1x,I4),/))
 33   format (' # ',I4,' localized MOs,',I4,' NON localized MOs')
 34   format (' # Localized MOs are',/,1000(20(1x,I4),/))
 766  format (' # NUMBER OF GROUPS = ',I2)
 764  format (' # GROUP ',I2,' MinElec and MaxElec (',a,') = ',2I6)
 767  format (1x,'# Number of atoms in the group = ',I3)
 765  format (' # ',10I6)
      end
