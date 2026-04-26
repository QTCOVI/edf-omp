
c-----------------------------------------------------------------------
c
      subroutine dafhcrit (sga,sgb,eleca,elecb,nmo,ntotattr,nattr,nspin,
     &    ngroup,nfugrp,ifugrp,pcut,epsdafh,lr,lw,lerr)
c
c.......................................................................
c
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
      include    'stderr.inc'
      complex*16  sga(ngroup,nmo,nmo)
      complex*16  sgb(ngroup,nmo,nmo)

      complex*16, allocatable,dimension (:,:)      :: sdiag,bdiag
      complex*16, allocatable,dimension (:)        :: wdiag
      complex*16, allocatable,dimension (:,:,:)    :: sg,sgna,sgnb,sgn
      real(kind=8), allocatable,dimension (:)      :: seigen
      integer,    allocatable,dimension (:)        :: eignz
      integer(kind=4)                              :: eleca(2,ngroup)
      integer(kind=4)                              :: elecb(2,ngroup)
      integer(kind=4), allocatable,dimension (:,:) :: resnca,resncb
      integer(kind=4), allocatable,dimension (:,:) :: resncx
      integer,    allocatable,dimension (:)        :: mocore
      complex*16  val1,val2
      integer     nzeig,zeig
      integer     nfugrp(ngroup),ifugrp(ntotattr,ngroup)
      character(len=5) spin(2)
c
      call timer (2,idafhcrit,'_dafhcrit ',-1)
c
c.....Compute complex Group Overlap integrals of all but the last one.
c
      if (.not.allocated(sg)) then
        allocate (sg(ngroup,nmo,nmo),stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot allocate sg()'
      endif
c
      spin(1) = 'ALPHA'
      spin(2) = 'BETA '
      ispin = 1
 22   continue
      if (ispin.eq.1) sg = sga
      if (ispin.eq.2) sg = sgb
c
      if (.not.allocated(sdiag)) then
        allocate (sdiag(nmo,nmo),stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot allocate sdiag()'
      endif
      if (.not.allocated(bdiag)) then
        allocate (bdiag(nmo,nmo),stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot allocate bdiag()'
      endif
      if (.not.allocated(seigen)) then
        allocate (seigen(nmo),stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot allocate seigen()'
      endif
      if (.not.allocated(wdiag)) then
        allocate (wdiag(nmo+nmo-1),stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot allocate wdiag()'
      endif
      if (.not.allocated(mocore)) then
        allocate (mocore(ngroup),stat=ier)
        if (ier.ne.0) stop '# dnocrit.f: Cannot allocate mocore()'
      endif
c
c.....Compute Group Overlap matrix of the sum of groups 1 to NGROUP-1
c
      sdiag = cmplx(0.0d+00,0.0d+00)
      do i=1,ngroup-1
        sdiag(:,:)=sdiag(:,:)+sg(i,:,:)
      enddo
      call zjacobi (sdiag,bdiag,wdiag,seigen,nmo)
      pop = sum(seigen(1:nmo))
      write (lw,112) trim(spin(ispin))
      write (lw,222) (seigen(i),i=1,nmo)
      write (lw,113) pop
c
c-----Determine the non-zero elements of seigen().
c
      if (.not.allocated(eignz)) then
        allocate (eignz(nmo),stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot allocate eignz()'
      endif
      nzeig = 0
      pop = 0.0d+00
      do i=1,nmo
        if (seigen(i).gt.epsdafh) then
          nzeig = nzeig+1
          pop = pop+seigen(i)
          eignz(nzeig) = i
        endif
      enddo
      zeig = nmo-nzeig
      write (lw,111) epsdafh,zeig
      write (lw,114) epsdafh,nzeig
      write (lw,222) (seigen(eignz(i)),i=1,nzeig)
      write (lw,113) pop
c
c-----Reconstruct the GOM in terms of MOs from the diagonalization
c
      if (.not.allocated(sgn)) then
        allocate (sgn(ngroup,nzeig,nzeig),stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot allocate sgn()'
      endif
      do k=1,ngroup-1
        do i=1,nzeig
          i1 = eignz(i)
          do j=1,i
            j1 = eignz(j)
            val1 = cmplx(0.0d+00,0.0d+00)
            do l=1,nmo
              do m=1,nmo
                val1=val1+conjg(sdiag(l,i1))*sdiag(m,j1)*sg(k,l,m)
              enddo
            enddo
            sgn(k,i,j) = val1
            sgn(k,j,i) = conjg(sgn(k,i,j))
          enddo
        enddo
      enddo
c
c.....Compute the Group Overlap Matrix for the last group.
c
      sgn(ngroup,1:nzeig,1:nzeig) = (0.0d+00,0.0d+00)
      do m=1,nzeig
        sgn(ngroup,m,m) = (1.0d+00,0.0d+00)
      enddo
      do i=1,ngroup-1
        sgn(ngroup,:,:) = sgn(ngroup,:,:) - sgn(i,:,:)
      enddo
      if (ispin.eq.1) then
        eleca(1,ngroup) = max(0,nzeig - sum(eleca(2,1:ngroup-1)))
        eleca(2,ngroup) = nzeig - sum(eleca(1,1:ngroup-1))
        if (.not.allocated(sgna)) then
          allocate (sgna(ngroup,nzeig,nzeig),stat=ier)
          if (ier.ne.0) stop '# dnocrit.f: Cannot allocate sgna()'
        endif
        sgna(:,1:nzeig,1:nzeig) = sgn(:,1:nzeig,1:nzeig)
      else
        elecb(1,ngroup) = max(0,nzeig - sum(elecb(2,1:ngroup-1)))
        elecb(2,ngroup) = nzeig - sum(elecb(1,1:ngroup-1))
        if (.not.allocated(sgnb)) then
          allocate (sgnb(ngroup,nzeig,nzeig),stat=ier)
          if (ier.ne.0) stop '# dnocrit.f: Cannot allocate sgnb()'
        endif
        sgnb(:,1:nzeig,1:nzeig) = sgn(:,1:nzeig,1:nzeig)
      endif
      if (allocated(sgn)) then
        deallocate (sgn,stat=ier)
        if (ier.ne.0) stop '# dnocrit.f: Cannot deallocate sgn()'
      endif
      mocore(1:ngroup) = 0
      mocore(ngroup)   = zeig
      write (lw,766) ngroup
      do i=1,ngroup
        if (ispin.eq.1) then
          ielmin = eleca(1,i)+mocore(i)
          ielmax = eleca(2,i)+mocore(i)
        else
          ielmin = elecb(1,i)+mocore(i)
          ielmax = elecb(2,i)+mocore(i)
        endif
        write (lw,764) i,trim(spin(ispin)),ielmin,ielmax
        write (lw,767) nfugrp(i)
        write (lw,765) (ifugrp(j,i),j=1,nfugrp(i))
      enddo
c
c-----Compute number of probabilities.
c
      if (ispin.eq.1) then
        npais = 1
        do i=1,ngroup-1
          npais = npais*(eleca(2,i)-eleca(1,i)+1)
        enddo
        nprev = npais
        if (.not.allocated(resncx)) then
          allocate (resncx(nprev,ngroup),stat=ier)
          if (ier.ne.0) stop '# dnocrit.f: Cannot allocate resncx()'
        endif
c
c-------Computation of ALPHA resonance structures.
c
        call rnprobs (eleca,resncx,nzeig,ngroup,npais,lw)
        if (nprev.lt.npais) then
          write (lerr,*) '# dnocrit.f: rnprobs.f returned NPREV < NPAIS'
          stop
        endif
        if (.not.allocated(resnca)) then
          allocate (resnca(npais,ngroup),stat=ier)
          if (ier.ne.0) stop '# dnocrit.f: Cannot allocate resncx()'
        endif
        resnca(1:npais,1:ngroup) = resncx(1:npais,1:ngroup)
        if (allocated(resncx)) then
          deallocate (resncx,stat=ier)
          if (ier.ne.0) stop '# dnocrit.f: Cannot deallocate resncx()'
        endif
        nzeiga = nzeig
      else
        npbis = 1
        do i=1,ngroup-1
          npbis = npbis*(elecb(2,i)-elecb(1,i)+1)
        enddo
        nprev = npbis
        if (.not.allocated(resncx)) then
          allocate (resncx(nprev,ngroup),stat=ier)
          if (ier.ne.0) stop '# dnocrit.f: Cannot allocate resncx()'
        endif
c
c-------Computation of ALPHA resonance structures.
c
        call rnprobs (elecb,resncx,nzeig,ngroup,npbis,lw)
        if (nprev.lt.npbis) then
          write (lerr,*) '# dnocrit.f: rnprobs.f returned NPREV < NPBIS'
          stop
        endif
        if (.not.allocated(resncb)) then
          allocate (resncb(npbis,ngroup),stat=ier)
          if (ier.ne.0) stop '# dnocrit.f: Cannot allocate resncb()'
        endif
        resncb(1:npbis,1:ngroup) = resncx(1:npbis,1:ngroup)
        if (allocated(resncx)) then
          deallocate (resncx,stat=ier)
          if (ier.ne.0) stop '# dnocrit.f: Cannot deallocate resncx()'
        endif
        nzeigb = nzeig
      endif
      if (ispin.eq.1.and.nspin.eq.2) then
        ispin = 2
        goto 22
      else
        elecb(1,1:ngroup) = eleca(1,1:ngroup)
        elecb(2,1:ngroup) = eleca(2,1:ngroup)
        npbis = npais
        if (.not.allocated(resncb)) then
          allocate (resncb(npbis,ngroup),stat=ier)
          if (ier.ne.0) stop '# dnocrit.f: Cannot allocate resncb()'
        endif
        resncb(1:npbis,1:ngroup) = resnca(1:npais,1:ngroup)
        if (.not.allocated(sgnb)) then
          allocate (sgnb(ngroup,nzeig,nzeig),stat=ier)
          if (ier.ne.0) stop '# dnocrit.f: Cannot allocate sgnb()'
        endif
        sgnb(:,1:nzeig,1:nzeig) = sgna(:,1:nzeig,1:nzeig)
        nzeigb = nzeiga
      endif
      call xcritcplx (pcut,npais,npbis,ngroup,nzeiga,nzeigb,lw,
     &     eleca,elecb,sgna,sgnb,resnca,resncb,mocore)
c
c-----Deallocate arrays.
c
      if (allocated(sg)) then
        deallocate (sg,stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot deallocate sg()'
      endif
      if (allocated(sdiag)) then
        deallocate (sdiag,stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot deallocate sdiag()'
      endif
      if (allocated(bdiag)) then
        deallocate (bdiag,stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot deallocate bdiag()'
      endif
      if (allocated(seigen)) then
        deallocate (seigen,stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot deallocate seigen()'
      endif
      if (allocated(wdiag)) then
        deallocate (wdiag,stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot deallocate wdiag()'
      endif
      if (allocated(eignz)) then
        deallocate (eignz,stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot deallocate eignz()'
      endif
      if (allocated(sgn)) then
        deallocate (sgn,stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot deallocate sgn()'
      endif
      if (allocated(resncb)) then
        deallocate (resncb,stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot deallocate resncb()'
      endif
      if (allocated(resnca)) then
        deallocate (resnca,stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot deallocate resnca()'
      endif
      if (allocated(mocore)) then
        deallocate (mocore,stat=ier)
        if (ier.ne.0) stop '# dafhcrit.f: Cannot deallocate mocore()'
      endif
c
      call timer (4,idafhcrit,'_dafhcrit ',-1)
      return
c
 222  format (4(1x,F17.10))
 111  format (' # NUMBER OF EIGENVALUES < ',E17.10,3x,' = ',I5)
 766  format (/,' #',/,' # NUMBER OF GROUPS = ',I2)
 764  format (' # GROUP ',I2,' MinElec and MaxElec (',a,') = ',2I6)
 765  format (' # ',10I6)
 112  format (/,' # ',70('+'),/,' # Spin Block ',a,/,
     & ' # Ordered Eigenvalues of the GOM of groups 1...NGROUP-1',/,
     & ' # (GOM = Group Overlap Matrix)',/,1x,'#',70('+'))
 113  format (' # SUM = ',F17.10)
 114  format (' # DAFH Eigenvalues > ',E17.10,3x,' = ',I5)
 767  format (1x,'# Number of atoms in the group = ',I3)
      end
