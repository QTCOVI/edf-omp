c
c-----------------------------------------------------------------------
c
      subroutine fnofrag (aominp,inside,cciqa,mal,mbe,natoms,nprims,
     & nmo,ncent,udat,wfnfile,lw,lerr,verbose)
c
c-----PARAMETERS--------------------------------------------------------
c
c-----aominp()                        INPUT
c
c     Atomic Overlap Matrix (AOM) of Canonical MOs in all the centers.
c
c-----inside (1..natoms)              INPUT 
c
c     Indices of the NATOMS that define the fragment.
c
c-----cciqa        INPUT
c
c     Each of them is .TRUE. or .FALSE. depending on the type of WFN used
c
c-----mal,mbe                         INPUT
c
c     Number of ALPHA and BETA electrons     
c
c-----natoms                          INPUT
c
c     Number of atoms of the fragment (fragment is called A from now on)
c
c-----nprims                          INPUT
c
c     Number of primitives of the WFN
c
c-----nmo                             INPUT 
c
c     Number of canonical MOs in the WFN
c
c-----ncent                           INPUT
c
c     Number of atoms of the molecule
c
c-----udat                            INPUT
c
c     Unit number of original WFN file
c
c-----wfnfile                         INPUT
c
c     Name of the WFN file
c
c-----lw                              INPUT
c
c     Output unit
c
c-----lerr
c
c     Output unit for errors
c
c-----verbose                         INPUT
c
c     Logical variable: if verbose=.true. a large output is requested,
c     otherwise a short output is used.
c
c-----------------------------------------------------------------------
c                         
      USE      space_for_wfncoef
      USE      space_for_wfnbasis
      USE      space_for_rdm1
      include 'implicit.inc'
      include 'corr.inc'
      include 'mline.inc'
      real   (kind=8) aominp(ncent,nmo,nmo)
      integer(kind=4) inside(ncent)
      integer(kind=4) insidex(ncent)
      real   (kind=8) rho(nmo,nmo)
      real   (kind=8) occwfn(nmo)
      integer(kind=4) natoms
      integer(kind=4) mal,mbe
      integer(kind=4) iordatoms(ncent)
      integer(kind=4) typewfn
      logical cciqa
      character(len=4) fourchar
      logical          verbose
      integer(kind=4)  udat,udatnw
      character(len=*) wfnfile
      character (len=mline) wfnloc

      real   (kind=8),allocatable, dimension (:,:,:) :: aom
      real   (kind=8),allocatable, dimension (:,:,:) :: aomab
      real   (kind=8),allocatable, dimension (:,:)   :: sg
      real   (kind=8),allocatable, dimension (:,:)   :: cmat
      real   (kind=8),allocatable, dimension (:,:)   :: cmatord
      real   (kind=8),allocatable, dimension (:,:)   :: wvec
      real   (kind=8),allocatable, dimension (:,:)   :: uvec
      real   (kind=8),allocatable, dimension (:,:)   :: ueu
      real   (kind=8),allocatable, dimension (:,:)   :: onerdm
      real   (kind=8),allocatable, dimension (:,:)   :: rhoa
      real   (kind=8),allocatable, dimension (:,:)   :: dcoef
      real   (kind=8),allocatable, dimension (:,:)   :: cc
      real   (kind=8),allocatable, dimension (:)     :: epsi
      real   (kind=8),allocatable, dimension (:)     :: lambda
      real   (kind=8),allocatable, dimension (:)     :: xloci
      real   (kind=8),allocatable, dimension (:)     :: diagord
      real   (kind=8),allocatable, dimension (:)     :: diagaom
      real   (kind=8),allocatable, dimension (:)     :: ocup
      real   (kind=8),allocatable, dimension (:,:)   :: elec
      real   (kind=8),allocatable, dimension (:)     :: elecord
      integer(kind=4),allocatable, dimension (:)     :: ieos
      integer(kind=4),allocatable, dimension (:)     :: iaom
      integer(kind=4) isigma(nmo)
      integer(kind=4) nant(2)
      character (len=1) spin

c
c-----------------------------------------------------------------------
c
      call timer (2,ifnofrag,'_fnofrag  ',-1)
!
!-----CCIQA wavefunctions with Nalpha # Nbeta can not be used yet.
!
      if (cciqa.and.mal.ne.mbe) then
        write (lerr,*) '# fnofrag.f: Open-shell CCWFN not allowed yet'
        return
      endif 

      nspin = 2
c
c-----Ordering the atoms of the fragment.
c
      insidex = inside
      forall (i=1:natoms) iordatoms(i)=i
      call iqcksort (insidex, iordatoms, ncent, 1, natoms)
      forall (i=1:natoms) inside(i)=insidex(iordatoms(i))      
      do ispin=1,nspin
        if (ispin.eq.1) then
          spin = 'a'
          nab = nalpha
          rho = c1ea
          isigma(1:nab) = ialpha(1:nab)
        else
          spin = 'b'
          nab = nbeta
          rho = c1eb
          isigma(1:nab) = ibeta(1:nab)
        endif
        allocate (aom(1:ncent,1:nab,1:nab))
        allocate (aomab(1:ncent,1:nab,1:nab))
        allocate (sg(1:nab,1:nab))
        allocate (rhoa(1:nab,1:nab))
        allocate (wvec(1:nab,1:nab))
        allocate (uvec(1:nab,1:nab))
        allocate (ueu(1:nab,1:nab))
        allocate (onerdm(1:nab,1:nab))
        allocate (epsi(1:nab))
        allocate (lambda(1:nab))
        allocate (xloci(1:nab))
        allocate (cmat(1:nab,1:nab))
        allocate (diagaom(1:nab))
        allocate (diagord(1:nab))
        allocate (iaom(1:nab))
        allocate (ocup(nab))
        allocate (dcoef(1:nab,1:nprims))
c
c-------Extract from the current AOM matrix the part associated
c       to electrons with spin ALPHA or BETA
c
        aomab(:,1:nab,1:nab)=aominp(:,isigma(1:nab),isigma(1:nab))
c
c-------Extract from the current 1RDM matrix the part associated
c       to electrons with spin ALPHA or BETA
c
        onerdm(1:nab,1:nab)=rho(isigma(1:nab),isigma(1:nab))
c
c-------Group overlap matrix (S^A) in the fragment
c
        sg = 0.0d+00
        do i=1,natoms
          sg = sg + aomab(inside(i),:,:)
        enddo
        write (lw,987) (inside(i),i=1,natoms)
        if (ispin.eq.1) write (lw,34) 'ALPHA'
        if (ispin.eq.2) write (lw,34) 'BETA '
        if (verbose) then
          write (lw,*) '# Matrix S^A'
          do i=1,nab
            write (lw,4444) (sg(i,j),j=1,i)
          enddo
        endif
c
c-------Diagonalize S^A
c
        call jacobi (sg,nab,nab,epsi,wvec,nrot)
        write (lw,988) (epsi(i),i=1,nab)
        if (verbose) write (lw,989)
        do i=1,nab
          if (verbose) write (lw,990) i,(wvec(j,i),j=1,nab)
        enddo
        do j=1,nab
          wvec(:,j) = wvec(:,j) * sqrt(max(epsi(j),1.0D-20))
        enddo
        rhoa = matmul(transpose(wvec),matmul(onerdm,wvec))
        if (verbose) then
          write (lw,*) 
     &     '# Matrix \bar \rho_A = €^{1/2}*W^t * 1RDM * W*€^{1/2}'
          do i=1,nab
            write (lw,4444) (rhoa(i,j),j=1,i)
          enddo
        endif
        call jacobi (rhoa,nab,nab,lambda,uvec,nrot)
        do j=1,nab
          wvec(:,j) = wvec(:,j) / max(epsi(j),1.0D-20)
        enddo
        cmat = matmul(wvec,uvec)
        if (verbose) then
          write (lw,*) '# Transformation matrix from CANMOs to FNOs'
          write (lw,*) '# C = W * €^{-1/2} * U'
          do k=1,nab
            write (lw,*) '# FNO ',k
            write (lw,'(5(1x,E17.10))') (cmat(l,k),l=1,nab)
          enddo
        endif
c
c-------AOM in R^3 of FNOs
c
        if (verbose) then
          write (lw,*) '# U^t * €^{-1} * U = AOM of FNOs in R^3'
        endif
        do i=1,nab
          do j=1,nab
            value=0.0d+00
            do k=1,nab
              value=value+uvec(k,i)*uvec(k,j)/max(1.0D-20,epsi(k))
            enddo
            ueu(i,j)=value
          enddo
          if (verbose) then
            write (lw,'(5(1x,E17.10))') (ueu(i,j),j=1,i)
          endif
        enddo
c
c-------Normalized FNOs in R^3
c
        forall (j=1:nab) cmat(:,j)=cmat(:,j)/sqrt(ueu(j,j))
c
c-------cmat() contains now the Normalized FNOs in R^3
c
c
c-------Overlaps in R^3 of normalized in R^3 FNOs
c
        do i=1,nab
          do j=1,nab
            ueu(i,j)=ueu(i,j)/sqrt(ueu(i,i)*ueu(j,j))
          enddo
        enddo
        if (verbose) then
          write (lw,37)
          do k=1,nab
            write (lw,'(5(1x,F17.10))') (ueu(k,l),l=1,k)
          enddo
        endif
c
c-------Compute AOM between normalized in R^3 FNOs in all the atoms
c
        do i=1,ncent
          aom(i,:,:)=matmul(transpose(cmat),matmul(aomab(i,:,:),cmat))
        enddo
c
c-------Ordering FNOs by decreasing value of the diagonal overlap in the
c       group. 
c
        diagaom=0.0d+00
        do ip=1,nab
          iaom(ip)=ip
          do i=1,natoms
            diagaom(ip) = diagaom(ip) + aom(inside(i),ip,ip)
          enddo
        enddo
        call qcksort (diagaom, iaom, 1, nab)
        forall (i=1:nab) diagord(i)=diagaom(iaom(nab-i+1))
        write (lw,332) (diagord(i),i=1,nab)

        write (lw,68) 'normalized in R^3 FNO'
        nrepcent = ncent/10
        if (mod(ncent,10).gt.0) nrepcent=nrepcent+1
        inic=1
        ifin=min(10,ncent)
        do m=1,nrepcent
          write (lw,161) (i,i=inic,ifin)
          do ip=1,nab
            i=iaom(nab-ip+1)
            xloci(ip) = 0.0d+00  ! Effective number of centers of each FNO
            do k = 1, ncent
              spp = aom(k,i,i)
              xloci(ip) = xloci(ip) + spp * spp
            enddo
            write (lw,412) ip,(100.0d+00*aom(k,i,i),k=inic,ifin)
          enddo
          inic=ifin+1
          ifin=min(ifin+10,ncent)
        enddo

        forall (ip=1:nab) ocup(ip)=lambda(iaom(nab-ip+1))
        write (lw,991) (ocup(ip),ip=1,nab)
        write (lw,4445) sum(lambda(1:nab))
        write (lw,520) (1d0/xloci(ip),ip=1,nab)
c
c-------Find FNOs in the primitive basis, cc(). 
c
        allocate (cmatord(1:nab,1:nab))
        do i=1,nab
          cmatord(:,i) = cmat(:,iaom(nab-i+1))
        enddo
        if (verbose) then
          write (lw,*) '# Normalized in R^3 FNOs from canonical MOs'
          do i=1,nab
            write (lw,*) '  Normalized in R^3 FNO number ',i
            write (lw,'(5(1x,F17.10))') (cmatord(j,i),j=1,nab)
          enddo
        endif
        do i=1,nab
          dcoef(i,1:nprims) = coef(nmo+isigma(i),1:nprims)
        enddo
        cc = matmul(transpose(cmatord),dcoef)
c
c-------write a WFN file with FNOs
c
        inpr = nprims
        wfnloc = trim(wfnfile)//"-fno"
        if (nspin.eq.2) wfnloc = trim(wfnfile)//"-fno"//spin(1:1)
        do i=1,natoms
          wfnloc = trim(wfnloc)//"-"//fourchar(inside(i))
        enddo
        udatnw = udat+10
        open (unit=udatnw,file=trim(wfnloc),status='unknown')
        write (lw,*)
        write (lw,351) trim(wfnloc)
        call cdafhmos (udat,udatnw,0,ocup,cc,nab,inpr)
        deallocate (ocup,stat=ier)
        if (ier.ne.0) stop '# fnofrag.f: Cannot deallocate ocup()'

        if (ispin.eq.1) malmbe = mal
        if (ispin.eq.2) malmbe = mbe

        deallocate (aom)
        deallocate (aomab)
        deallocate (sg)
        deallocate (cmat)
        deallocate (cmatord)
        deallocate (wvec)
        deallocate (uvec)
        deallocate (onerdm)
        deallocate (rhoa)
        deallocate (ueu)
        deallocate (dcoef)
        deallocate (cc)
        deallocate (epsi)
        deallocate (lambda)
        deallocate (xloci)
        deallocate (diagaom)
        deallocate (diagord)
        deallocate (iaom)
      enddo
      call timer (4,ifnofrag,'_fnofrag  ',-1)
      return
c
c.....Formats
c
 410  format (' # ',I4)
 411  format (1000(17x,10(1x,F6.2,'%'),/))
 68   format (/' # Localization of each ',a,' in each atom')
 710  format (' # Atom ',I4,'   Electrons = ',F18.8)
 711  format (' # SUM  ',17x,'= ',F18.8)
 351  format (1x,80('+'),/,1x,"+ File '",a,"'",
     & " contains the Normalized in R^3 FNOs",/,1x,80('+'))
 353  format (1x,"# Writing file '",a,"'",1x,"with AOM for FNOs")
 101  format (1000(' # ',5(1x,F15.8),/))
 80   format (6(1x,e16.10))
 520  format (' # Atoms expanded by each FNO',/,1000(5(1x,F17.10),/))
 987  format (//' # ',77('+'),/,
     & ' # FNO ANALYSIS - FNO ANALYSIS - FNO ANALYSIS - FNO ANALYSIS',/
     & ' # FNOs are not normalized in R^3',/,
     & ' # NFNOs ARE normalized in R^3',/,' # ',77('+'),/,
     & ' # CURRENT Fragment (A in the following) formed by atoms:',/,
     & 1000(' # ',20(1x,I3),/))
 988  format (' # Eigenvalues of S^A are',/,1000(5(1x,F17.10),/))
 994  format (' # Full set of FNOs TOTAL POPULATION = ',F15.8,/,
     & ' # ',77('-'))
 991  format (' # ',77('-'),/,
     & ' # Fragment Natural Orbitals (FNO) occupations',/,' # ',
     & 77('-'),/,1000(5(1x,F17.10),/))
 989  format (' # Eigenvectors of S^A')
 990  format (' # Eigenvector ',I3,/,1000(5(1x,F17.10),/))
 42   format (' # Eigenvalues    > ',1PE15.8,' are ',I4)
 332  format (' # Diagonal Overlaps in the fragment'
     & ' of normalized in R^3 FNOs',/,1000(5(1x,F17.10),/))
 33   format (' # ALPHA Fragment Natural MOs (FNO) = BETA FNOs'/,
     &        ' # ALPHA = BETA FNOs are determined')
 34   format (' # ',a,' FNOs are determined')
 37   format (' # AOM in R^3 of the normalized in R^3 FNOs (NFNOs)',/,
     &        ' # i.e. AOM in R^3 of NFNOs')
 4444 format (5(1x,F17.10))
 4445 format (1x,F17.10,
     &  '  <-- Total Electron Population',/,' # ',117('-'))
 77   format (1000(2x,5(1x,F17.10),/))
 78   format (1000(2x,10(1x,I4),/))
 79   format (10x,F17.10,4x,I4,3x,I4)
 161  format (' #  MO\ATOM   ',10I9)
 412  format (' # ',I4,10x,10(1x,F7.2,'%'))
 770  format (1000(2x,5(1x,F17.10,2x,' ( ',I2,' ) '),/))
 24   format (I7,2x,F17.10,2x,I4,6x,I4)
      end
