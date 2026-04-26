c
c-----------------------------------------------------------------------
c
      subroutine neweos (aominp,atoms_in_frag,ifrag,cciqa,mal,mbe,
     & nprims,nmo,ncent,nfrag,critoverlap,udat,wfnfile,lw,lerr,verbose)
c
c-----PARAMETERS--------------------------------------------------------
c
c-----aominp()                        INPUT
c
c     Atomic Overlap Matrix (AOM) of Canonical MOs in all the centers.
c
c-----atoms_in_frag(1..nfrag)
c
c     Number of atoms included in each fragment
c
c-----ifrag()
c
c     Labels of the atoms included in each of the 'nfrag' fragments
c
c-----cciqa                           INPUT
c
c     .TRUE. or .FALSE. depending on the type of WFN used
c
c-----mal,mbe                         INPUT
c
c     Number of ALPHA and BETA electrons     
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
c-----nfrag
c
c     Number of fragments in which the molecule has been divided
c
c
c-----critoverlap                     INPUT
c
c
c-----Critical value of the overlap between two normalized in R^3 
c     FNOs that determines whether one considers that they almost
c     the same or not.
c
c
c.....udat                            INPUT
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
*     real   (kind=8) occwfn(nmo)
      integer(kind=4) atoms_in_frag(nfrag)
      integer(kind=4) ifrag(ncent,nfrag)
      integer(kind=4) natoms
      integer(kind=4) mal,mbe
      integer(kind=4) iordatoms(ncent)
      integer(kind=4) typewfn
      logical cciqa
      character(len=4) fourchar
      logical          verbose
      logical          samepopandz
      logical          samez
      integer(kind=4)  udat,udatnw
      character(len=*) wfnfile
      character (len=mline) wfnloc

      real   (kind=8),parameter                      :: epspop = 1.0d-08
      real   (kind=8),parameter                      :: epseig = 1.0d-08
      real   (kind=8),allocatable, dimension (:,:,:) :: aom
      real   (kind=8),allocatable, dimension (:,:,:) :: aomab
      real   (kind=8),allocatable, dimension (:,:)   :: sg,sgo
      real   (kind=8),allocatable, dimension (:,:)   :: SArhoSA
      real   (kind=8),allocatable, dimension (:,:)   :: cmat
      real   (kind=8),allocatable, dimension (:,:,:) :: fno
      real   (kind=8),allocatable, dimension (:,:)   :: cmatord
      real   (kind=8),allocatable, dimension (:,:)   :: wvec,wvecx
      real   (kind=8),allocatable, dimension (:,:)   :: uvec
      real   (kind=8),allocatable, dimension (:,:)   :: ueu
      real   (kind=8),allocatable, dimension (:,:)   :: onerdm
      real   (kind=8),allocatable, dimension (:,:)   :: rhoa
      real   (kind=8),allocatable, dimension (:,:)   :: dcoef
      real   (kind=8),allocatable, dimension (:,:)   :: cc
      real   (kind=8),allocatable, dimension (:)     :: epsi,epsix
      real   (kind=8),allocatable, dimension (:)     :: dswork
      real   (kind=8),allocatable, dimension (:)     :: lambda
      real   (kind=8),allocatable, dimension (:)     :: xloci
      real   (kind=8),allocatable, dimension (:)     :: diagord
      real   (kind=8),allocatable, dimension (:)     :: diagaom
      real   (kind=8),allocatable, dimension (:)     :: ocup
      real   (kind=8),allocatable, dimension (:,:)   :: ocupt

      real   (kind=8),allocatable, dimension (:)     :: ocupdeg
      integer(kind=8),allocatable, dimension (:)     :: iset
      integer(kind=8),allocatable, dimension (:,:)   :: nindex
      integer(kind=8),allocatable, dimension (:)     :: MOinG
      integer(kind=8),allocatable, dimension (:,:)   :: iMOinG
      integer(kind=8),allocatable, dimension (:)     :: listG

      real   (kind=8),allocatable, dimension (:,:)   :: elec

      real   (kind=8),allocatable, dimension (:,:)   :: sumepsi
      real   (kind=8),allocatable, dimension (:)     :: epord
      integer(kind=4),allocatable, dimension (:)     :: ipord

      real   (kind=8),allocatable, dimension (:)     :: elecord
      real   (kind=8),allocatable, dimension (:,:)   :: elecsave
      integer(kind=4),allocatable, dimension (:)     :: ieos
      integer(kind=4),allocatable, dimension (:,:)   :: iatom
      integer(kind=4),allocatable, dimension (:)     :: iaom
      real   (kind=8),allocatable, dimension (:)     :: zcharge
      integer(kind=4),allocatable, dimension (:,:)   :: nelec
      logical,        allocatable, dimension (:,:)   :: savemo
      integer(kind=4),allocatable, dimension (:)     :: similarFNO
      integer(kind=4) isigma(nmo)
      integer(kind=4) nant(2)
      real   (kind=8) eos(nfrag),eleca(nfrag),elecb(nfrag)
      character (len=1) spin
      real   (kind=8) rsigma(2)
      real   (kind=8) sigma(nfrag,2)
      real   (kind=8) ov(nmo,nmo)

c
c-----------------------------------------------------------------------
c
      call timer (2,ineweos,'_neweos   ',-1)
!
!-----CCIQA wavefunctions with Nalpha # Nbeta can not be used yet.
!
      if (cciqa.and.mal.ne.mbe) then
        write (lerr,*) '# neweos.f: Open-shell CCWFN not allowed yet'
        return
      endif 

      write (lw,100)
      nspin = 2
c
c-----LOOP over all the atoms of the molecule
c
      allocate (elec(nmo*nfrag,nspin))
      allocate (sumepsi(nmo*nfrag,nspin))
      allocate (iatom(nmo*nfrag,nspin))
      allocate (fno(nmo,nmo*nfrag,nspin))
      allocate (ocupt(nmo*nfrag,nspin))
      elec = 0.0d+00
      nant(1:nspin) = 0
      do igroup=1,nfrag
        natoms = atoms_in_frag(igroup)
        inside(1:natoms) = ifrag(1:natoms,igroup)
c
c-------Ordering the atoms of the fragment.
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
          allocate (aomab(1:ncent,1:nab,1:nab))
          allocate (sg(1:nab,1:nab))
          allocate (SArhoSA(1:nab,1:nab))
          allocate (sgo(1:nab,1:nab))
          allocate (wvecx(1:nab,1:nab))
          allocate (onerdm(1:nab,1:nab))
          allocate (epsix(1:nab))
          allocate (iaom(1:nab))
          allocate (ocup(nab))
          allocate (dcoef(1:nab,1:nprims))
c
c---------Extract from the current AOM matrix the part associated
c         to electrons with spin ALPHA or BETA
c
          aomab(:,1:nab,1:nab)=aominp(:,isigma(1:nab),isigma(1:nab))
c
c---------Extract from the current 1RDM matrix the part associated
c         to electrons with spin ALPHA or BETA
c
          onerdm(1:nab,1:nab)=rho(isigma(1:nab),isigma(1:nab))
c
c---------Group overlap matrix (S^A) in the fragment
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
c---------Write the diagonal values of S^A ordered in reversed order
c
          forall (ip=1:nab) iaom(ip) = ip
          forall (ip=1:nab) epsix(ip) = sg(ip,ip)
          call qcksort (epsix, iaom, 1, nab)
          write (lw,*) '# Ordered diagonal S^A(i,i) values'
          write (lw,'(1x,a)') '# Original index & value'
          do i=1,nab
            j=iaom(nab-i+1)
            write (lw,'(I4,1x,F17.10)') j,epsix(j)
          enddo
c
c---------Built the triple product S^A * 1RDM * S^A and diagonalize it
c
          SArhoSA = matmul(matmul(sg,onerdm),sg)
          allocate (dswork(max(1,3*nab-1)))
          call dsyev ('V','U',nab,SArhoSA,nab,epsix,dswork,
     &                max(1,3*nab-1),info)
          if (info.ne.0) write (lerr,*)
     &      '# neweos.f: dsyev(SArhoSA) failed, info=',info
          deallocate (dswork)
          wvecx(1:nab,1:nab) = SArhoSA(1:nab,1:nab)
          newnab = 0
          do j=1,nab
            if (epsix(j).ge.epseig) newnab = newnab +1
          enddo
          allocate (epsi(1:nab))
          allocate (wvec(1:nab,1:nab))
          forall (ip=1:nab) iaom(ip) = ip
          call qcksort (epsix, iaom, 1, nab)
          do i=1,nab
            j = iaom(nab-i+1)
            epsi(i) = epsix(j)
            wvec(:,i) = wvecx(:,j)
          enddo
          deallocate (iaom)
          write (lw,988) (epsi(i),i=1,nab)
          write (lw,*) '# The total sum is ',sum(epsi(1:nab))
          if (verbose) then
            write (lw,989)
            do i=1,nab
              write (lw,990) i,(wvec(j,i),j=1,nab)
            enddo
          endif
c
c---------Compute AOM between normalized in R^3 FNOs in all the atoms
c
          allocate (aom(1:ncent,1:nab,1:nab))
          do i=1,ncent
            aom(i,:,:)=matmul(transpose(wvec),matmul(aomab(i,:,:),wvec))
          enddo
c
c---------Ordering FNOs by decreasing value of the diagonal overlap in the
c         fragment. 
c
          allocate (diagaom(1:nab))
          allocate (diagord(1:nab))
          allocate (iaom(1:nab))
          diagaom=0.0d+00
          do ip=1,nab
            iaom(ip)=ip
            do i=1,natoms
              diagaom(ip) = diagaom(ip) + aom(inside(i),ip,ip)
            enddo
          enddo
          call qcksort (diagaom,iaom,1,nab)
          forall (i=1:nab) diagord(i)=diagaom(iaom(nab-i+1))
          write (lw,332) (diagord(i),i=1,nab)
          allocate (xloci(1:nab))
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
                xloci(ip) = xloci(ip) + aom(k,i,i) * aom(k,i,i)
              enddo
              write (lw,412) ip,(100.0d+00*aom(k,i,i),k=inic,ifin)
            enddo
            inic=ifin+1
            ifin=min(ifin+10,ncent)
          enddo
          nn=nant(ispin)
          iatom(nn+1:nn+nab,ispin) = igroup
          do ip=1,nab
            elec(nn+ip,ispin) = epsi(iaom(nab-ip+1))
            sumepsi(nn+ip,ispin)  = epsi(iaom(nab-ip+1))
          enddo
c
c---------Store all FNOs in the fno() array
c
          fno(1:nab,nn+1:nn+nab,ispin) = wvec(1:nab,1:nab)
          nant(ispin) = nn + nab
          forall (ip=1:nab) ocup(ip)=epsi(iaom(nab-ip+1))
          write (lw,991) (ocup(ip),ip=1,nab)
          write (lw,4445) sum(ocup(1:nab))
          write (lw,520) (1d0/xloci(ip),ip=1,nab)
c
c---------Find FNOs in the primitive basis, cc(). 
c
          allocate (cmatord(1:nab,1:nab))
          do i=1,nab
            cmatord(:,i) = wvec(:,iaom(nab-i+1))
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
          allocate (cc(1:nab,nprims))
          cc = matmul(transpose(cmatord),dcoef)
c
c---------write a WFN file with FNOs
c
          inpr = nprims
          wfnloc = trim(wfnfile)//"-fno"//spin(1:1)
          do i=1,natoms
            wfnloc = trim(wfnloc)//"-"//fourchar(inside(i))
          enddo
          udatnw = udat+10
          open (unit=udatnw,file=trim(wfnloc),status='unknown')
          write (lw,*)
          write (lw,351) trim(wfnloc)
          call cdafhmos (udat,udatnw,0,ocup,cc,nab,inpr)

          if (ispin.eq.1) malmbe = mal
          if (ispin.eq.2) malmbe = mbe
          goto 1234






C
C---------Ordering the above MOs by decreasing value of the diagonal overlap in the
C         fragment. 
C
C         forall (ip=1:nab) iaom(ip) = ip
C         call qcksort (epsi, iaom, 1, nab)
C         write (lw,3320) (epsi(iaom(nab-i+1)),i=1,nab)

c
c---------Store temporarily in wvec() the product wvec * €^{1/2}
c         (Only eigenvalues €_i greater than 'epseig'
c
          do j=1,newnab
            wvec(:,j) = wvec(:,j) * sqrt(epsi(j))
          enddo
          allocate (rhoa(1:newnab,1:newnab))
          rhoa = matmul(transpose(wvec(1:nab,1:newnab)),
     &           matmul(onerdm,wvec(1:nab,1:newnab)))
          if (verbose) then
            write (lw,*) 
     &       '# Matrix \bar \rho_A = €^{1/2}*W^t * 1RDM * W*€^{1/2}'
            do i=1,newnab
              write (lw,4444) (rhoa(i,j),j=1,i)
            enddo
          endif
c
c---------Diagonalize to obtain the FNO occupations and eigenvectors
c
          allocate (uvec(1:newnab,1:newnab))
          allocate (lambda(1:newnab))
          allocate (dswork(max(1,3*newnab-1)))
          call dsyev ('V','U',newnab,rhoa,newnab,lambda,dswork,
     &                max(1,3*newnab-1),info)
          if (info.ne.0) write (lerr,*)
     &      '# neweos.f: dsyev(rhoa) failed, info=',info
          deallocate (dswork)
          uvec(1:newnab,1:newnab) = rhoa(1:newnab,1:newnab)
c
c
c---------Store now in wvec() the product wvec * €^{-1/2}
c         (Only eigenvalues €_i greater than 'epseig'
c         i.e : wvec <---- wvec * €^{-1} to undo the previos modification
c
          do j=1,newnab
            wvec(:,j) = wvec(:,j)/epsi(j)
          enddo
          allocate (cmat(1:nab,1:newnab))
          cmat = matmul(wvec(1:nab,1:newnab),uvec)
          if (verbose) then
            write (lw,*) '# Transformation matrix from CANMOs to FNOs'
            write (lw,*) '# C = W * €^{-1/2} * U'
            do k=1,newnab
              write (lw,*) '# FNO ',k
              write (lw,'(5(1x,E17.10))') (cmat(l,k),l=1,nab)
            enddo
          endif
c
c---------AOM in R^3 of FNOs
c
          if (verbose) then
            write (lw,*) '# U^t * €^{-1} * U = AOM of FNOs in R^3'
          endif
          allocate (ueu(1:newnab,1:newnab))
          do i=1,newnab
            do j=1,newnab
              value=0.0d+00
              do k=1,newnab
                value=value+uvec(k,i)*uvec(k,j)/max(1.0D-20,epsi(k))
              enddo
              ueu(i,j)=value
            enddo
            if (verbose) write (lw,'(5(1x,E17.10))') (ueu(i,j),j=1,i)
          enddo
c
c---------Normalized FNOs in R^3
c
          forall (j=1:newnab) cmat(:,j)=cmat(:,j)/sqrt(ueu(j,j))
c
c---------cmat() contains now the Normalized FNOs in R^3
c---------Overlaps in R^3 of normalized in R^3 FNOs
c
          do i=1,newnab
            do j=1,newnab
              ueu(i,j)=ueu(i,j)/sqrt(ueu(i,i)*ueu(j,j))
            enddo
          enddo
          if (verbose) then
            write (lw,37)
            do k=1,newnab
              write (lw,'(5(1x,F17.10))') (ueu(k,l),l=1,k)
            enddo
          endif
c
c---------Compute AOM between normalized in R^3 FNOs in all the atoms
c
          allocate (aom(1:ncent,1:newnab,1:newnab))
          do i=1,ncent
            aom(i,:,:)=matmul(transpose(cmat),matmul(aomab(i,:,:),cmat))
          enddo
c
c---------Ordering FNOs by decreasing value of the diagonal overlap in the
c         fragment. 
c
          allocate (diagaom(1:newnab))
          allocate (diagord(1:newnab))
          allocate (iaom(1:newnab))
          diagaom=0.0d+00
          do ip=1,newnab
            iaom(ip)=ip
            do i=1,natoms
              diagaom(ip) = diagaom(ip) + aom(inside(i),ip,ip)
            enddo
          enddo
          call qcksort (diagaom, iaom, 1, newnab)
          forall (i=1:newnab) diagord(i)=diagaom(iaom(newnab-i+1))
          write (lw,332) (diagord(i),i=1,newnab)

          allocate (xloci(1:newnab))
          write (lw,68) 'normalized in R^3 FNO'
          nrepcent = ncent/10
          if (mod(ncent,10).gt.0) nrepcent=nrepcent+1
          inic=1
          ifin=min(10,ncent)
          do m=1,nrepcent
            write (lw,161) (i,i=inic,ifin)
            do ip=1,newnab
              i=iaom(newnab-ip+1)
              xloci(ip) = 0.0d+00  ! Effective number of centers of each FNO
              do k = 1, ncent
                xloci(ip) = xloci(ip) + aom(k,i,i) * aom(k,i,i)
              enddo
              write (lw,412) ip,(100.0d+00*aom(k,i,i),k=inic,ifin)
            enddo
            inic=ifin+1
            ifin=min(ifin+10,ncent)
          enddo

          nn=nant(ispin)
          iatom(nn+1:nn+nab,ispin) = igroup
          do ip=1,newnab
            elec(nn+ip,ispin) = lambda(iaom(newnab-ip+1))
            sumepsi(nn+ip,ispin)  = epsi(iaom(newnab-ip+1))
          enddo
c
c---------Store all FNOs in the fno() array
c
          fno(1:nab,nn+1:nn+newnab,ispin) = cmat(1:nab,1:newnab)
          nant(ispin) = nn + newnab

          forall (ip=1:newnab) ocup(ip)=lambda(iaom(newnab-ip+1))
          write (lw,991) (ocup(ip),ip=1,newnab)
          write (lw,4445) sum(lambda(1:newnab))
          write (lw,520) (1d0/xloci(ip),ip=1,newnab)
c
c---------Find FNOs in the primitive basis, cc(). 
c
          allocate (cmatord(1:nab,1:newnab))
          do i=1,newnab
            cmatord(:,i) = cmat(:,iaom(newnab-i+1))
          enddo
          if (verbose) then
            write (lw,*) '# Normalized in R^3 FNOs from canonical MOs'
            do i=1,newnab
              write (lw,*) '  Normalized in R^3 FNO number ',i
              write (lw,'(5(1x,F17.10))') (cmatord(j,i),j=1,nab)
            enddo
          endif
          do i=1,nab
            dcoef(i,1:nprims) = coef(nmo+isigma(i),1:nprims)
          enddo
          allocate (cc(1:newnab,nprims))
          cc = matmul(transpose(cmatord),dcoef)
c
c---------write a WFN file with FNOs
c
          inpr = nprims
          wfnloc = trim(wfnfile)//"-fno"//spin(1:1)
          do i=1,natoms
            wfnloc = trim(wfnloc)//"-"//fourchar(inside(i))
          enddo
          udatnw = udat+10
          open (unit=udatnw,file=trim(wfnloc),status='unknown')
          write (lw,*)
          write (lw,351) trim(wfnloc)
          call cdafhmos (udat,udatnw,0,ocup,cc,newnab,inpr)
          deallocate (ocup,stat=ier)
          if (ier.ne.0) stop '# neweos.f: Cannot deallocate ocup()'

          if (ispin.eq.1) malmbe = mal
          if (ispin.eq.2) malmbe = mbe

 1234     continue
          deallocate (aom)
          deallocate (aomab)
          deallocate (sg)
          deallocate (sgo)
*         deallocate (cmat)
          deallocate (cmatord)
          deallocate (wvec)
          deallocate (wvecx)
*         deallocate (uvec)
          deallocate (onerdm)
*         deallocate (rhoa)
*         deallocate (ueu)
          deallocate (dcoef)
          deallocate (cc)
          deallocate (epsi)
          deallocate (epsix)
*         deallocate (lambda)
          deallocate (xloci)
          deallocate (diagaom)
          deallocate (diagord)
          deallocate (iaom)
          deallocate (sarhosa)
          deallocate (ocup)
        enddo
      enddo
c
c-----Ordering the AOM eigenvalues by decreasing value
c
      allocate (epord(nmo*nfrag))
      allocate (ipord(nmo*nfrag))
      write (lw,555) 
      do is=1,nspin
        nval = nant(is)
        epord(1:nval) = sumepsi(1:nval,is)
        forall (j=1:nval) ipord(j) = j
        call qcksort (epord, ipord, 1, nval)
        if (is.eq.1) write (lw,556) 'ALPHA'
        if (is.eq.2) write (lw,556) 'BETA '
        write (lw,772) (epord(ipord(j)),iatom(ipord(j),is),j=nval,1,-1)
      enddo
      deallocate (sumepsi)
      deallocate (epord)
      deallocate (ipord)
 555  format (' #',/,
     &  ' # Ordered Eigenvalues of the full set of Occupations',/,
     &  ' # In parenthesis the fragment to which it belongs')
 556  format (' # BLOCK ',a)
      allocate (nelec(nfrag,2))
      allocate (zcharge(nfrag))
      allocate (elecsave(nmo*nfrag,nspin))
      nelec(1:nfrag,1:2) = 0
      zcharge(1:nfrag) = 0.0d+00
      do i=1,nfrag
        do j=1,atoms_in_frag(i)
          zcharge(i) = zcharge(i) + charge(ifrag(j,i))
        enddo
      enddo
      write (lw,111)
      do ispin=1,nspin
        if (ispin.eq.1) then
          malmbe = mal
          write (lw,'(a,/,a)') ' #',' # ALPHA set'
        else
          malmbe = mbe
          write (lw,'(a,/,a)') ' #',' # BETA  set'
        endif
        if (verbose) then
          write (lw,300)
          write (lw,770) (elec(j,ispin),iatom(j,ispin),j=1,nant(ispin))
        endif
c
c-------Ordering occupations
c
        allocate (elecord(nant(ispin)))
        allocate (ieos(nant(ispin)))
        forall (i=1:nant(ispin)) elecord(i)=elec(i,ispin)
        forall (i=1:nant(ispin)) ieos(i)=i
        call qcksort (elecord,ieos,1,nant(ispin))
        write (lw,'(a,a)') 
     &   ' # Full set of ordered electron populations,',
     &   ' associated fragments and original order'       
        do j=1,nant(ispin)
          elecsave(j,ispin) = elecord(ieos(nant(ispin)-j+1))
        enddo
        write (lw,771) (elecord(ieos(nant(ispin)-j+1)),
     &      iatom(ieos(nant(ispin)-j+1),ispin),
     &     ieos(nant(ispin)-j+1),j=1,nant(ispin))
        total = sum(elecord(1:nant(ispin)))
        write (lw,'(a,F17.10)') ' # Total= ',total
c
c-------Analyze degenerations
c
        write (lw,*) '#'
        write (lw,*) '# +++++++++++++++++++++++++++++++++++++++++++++++'
        write (lw,*) '# Analyze degeneration of FNOs'
        write (lw,*) '# +++++++++++++++++++++++++++++++++++++++++++++++'
        allocate (ocupdeg(nant(ispin)))
        allocate (iset(nant(ispin)))
        allocate (nindex(nant(ispin),10))
        nsets = 1
        iset(nsets) = 1       
        ocupdeg(nsets) = elecsave(1,ispin)
        nindex(nsets,iset(nsets)) = 1
        cyclei: do i=2,nant(ispin)
          cyclej: do j=1,nsets
            dif1 = abs(elecsave(i,ispin)-ocupdeg(j))
            dif2 = abs(elecsave(i,ispin))
            if (dif1.lt.epspop.and.dif2.gt.1.0d+00) then
              iset(j) = iset(j) + 1
              if (iset(j).gt.10) then
                stop '# neweos.f: Increase the second dim of nindex()'
              endif
              nindex(j,iset(j)) = i
              cycle cyclei
            endif
          enddo cyclej
          nsets = nsets + 1
          iset(nsets) = 1
          ocupdeg(nsets) = elecsave(i,ispin)
          nindex(nsets,iset(nsets)) = i
        enddo cyclei
        do i=1,nsets
          write (lw,1617) ocupdeg(i),iset(i),(nindex(i,j),j=1,iset(i)) 
        enddo
        deallocate (ocupdeg)
        deallocate (iset)
        deallocate (nindex)
c
c-------Compute electrons that are associated to each fragment.
c       All the eigenvalues are analyzed by decreasing order.
c
        allocate (MOinG(nfrag))
        allocate (iMOinG(nfrag,40))
        allocate (listG(20))
        k = 0
        nelec = 0
        write (lw,240)
        MOinG(1:nfrag) = 0
        do i=1,nant(ispin)
          ii=ieos(nant(ispin)-i+1)
          jj=ieos(nant(ispin)-i)
c
c---------If there are two similar fragments, as determined by their
c         respective total nuclear charges, and the diffrence between two
c         successive occupations differ by less than 'epspop', the ALPHA
c         electron is assigned to fragment 1 and the BETA electron to 
c         fragment 2.
c
          difel = elecord(ii)-elecord(jj)
          samez = anint(zcharge(1)).eq.anint(zcharge(2))
          samepopandz = abs(difel).lt.epspop.and.nfrag.eq.2.and.samez
c
          if (samepopandz) then
            if (ispin.eq.1) then
              igr=iatom(ii,ispin)
              iatomo=iatom(ii,ispin)
              write (lw,24) i,elecord(ii),iatomo,ii
              MOinG(iatomo) = MOinG(iatomo) + 1
              if (MOinG(iatomo).gt.40) then
                stop '# neweos.f: Increase second dimension of iMOinG'
              endif
              iMOinG(iatomo,MOinG(iatomo)) = i
              k=k+1
              if (k.le.malmbe) nelec(igr,ispin)=nelec(igr,ispin)+1
              if (k.eq.malmbe) then
                rsigma(ispin) = 100.0*min(1d0,max(0.0d+00,0.5d0))
                exit
              endif
            else
              igr=iatom(jj,ispin)
              iatomo=iatom(jj,ispin)
              write (lw,24) i,elecord(jj),iatomo,jj
              MOinG(iatomo) = MOinG(iatomo) + 1
              if (MOinG(iatomo).gt.40) then
                stop '# neweos.f: Increase second dimension of iMOinG'
              endif
              iMOinG(iatomo,MOinG(iatomo)) = i
              k=k+1
              if (k.le.malmbe) nelec(igr,ispin)=nelec(igr,ispin)+1
              if (k.eq.malmbe) then
                rsigma(ispin) = 100.0*min(1d0,max(0.0d+00,0.5d0))
                exit
              endif
            endif
          else
            igr=iatom(ii,ispin)
            iatomo=iatom(ii,ispin)
            write (lw,24) i,elecord(ii),iatomo,ii
            MOinG(iatomo) = MOinG(iatomo) + 1
            if (MOinG(iatomo).gt.40) then
              stop '# neweos.f: Increase second dimension of iMOinG'
            endif
            iMOinG(iatomo,MOinG(iatomo)) = i
            k=k+1
            if (k.le.malmbe) nelec(igr,ispin)=nelec(igr,ispin)+1
            if (k.eq.malmbe) then
              rsigma(ispin) = 100.0*min(1d0,max(0.0d+00,0.5d0+difel))
              exit
            endif
          endif
        enddo
c
c-------Reliability EOS value
c
        do igr=1,nfrag
          if (ispin.eq.1) then
            eleca(igr)=dble(nelec(igr,ispin))
          else
            elecb(igr)=dble(nelec(igr,ispin))
          endif
        enddo

        ninlist = 0
        do i1=1,nfrag-1
          do i2=i1+1,nfrag
            if (MOinG(i1).ne.MOinG(i2) .and. 
     &        (abs(zcharge(i1)-zcharge(i2)).lt.1.0d-04)) then
              write (lw,98) 
              write (lw,99) i1,(iMOinG(i1,k),k=1,MOinG(i1))
              write (lw,99) i2,(iMOinG(i2,k),k=1,MOinG(i2))
              if (ninlist.eq.0) then
                listG(1)=i1
                listG(2)=i2
                ninlist = 2
              else
                do j=1,ninlist
                  if (i1.eq.listG(j)) goto 343
                enddo
                ninlist=ninlist+1
                listG(ninlist)=i1
 343            continue
                do j=1,ninlist
                  if (i2.eq.listG(j)) goto 344
                enddo
                ninlist=ninlist+1
                listG(ninlist)=i2
 344            continue
              endif
            endif
          enddo
        enddo
        if (ninlist.gt.0) then
          write (lw,1004) (listG(j),j=1,ninlist)
          maxFNOs = -100
          do j=1,ninlist
            if (MOinG(listG(j)).gt.maxFNOs) then
              maxFNOs = MOinG(listG(j))
              imaxFNO = listG(j)
            endif
          enddo
          write (lw,1005) imaxFNO,maxFNOs
c
c---------Recompute electrons of these conflicting fragments
c
          ntotFNOs = 0
          do j=1,ninlist
            ntotFNOs = ntotFNOs + MOinG(listG(j))
          enddo
          write (lw,1006) ntotFNOs
          addelec = dble(ntotFNOs)/dble(ninlist)
          if (ispin.eq.1) then
*           do j=1,ninlist
*             igr = listG(j)
*             eleca(igr) = 0.0d+00
*           enddo
*           addelec = dble(ntotFNOs)/dble(ninlist)
            do j=1,ninlist
              igr = listG(j)
              eleca(igr) = addelec
            enddo
          else
*           do j=1,ninlist
*             igr = listG(j)
*             elecb(igr) = 0.0d+00
*           enddo
            do j=1,ninlist
              igr = listG(j)
              elecb(igr) = addelec
            enddo
          endif
          naddfno=ninlist*maxFNOs - ntotFNOs
          write (lw,1007) naddfno
          write (lw,'(a,1000I4)') ' #',(j,j=malmbe+1,malmbe+naddfno)
          write (lw,240)
          do j=malmbe+1,malmbe+naddfno
            ii=ieos(nant(ispin)-j+1)
            iatomo=iatom(ii,ispin)
            write (lw,24) j,elecord(ii),iatomo,ii
          enddo
c
c---------Recompute the Reliability index
c
          ii=ieos(nant(ispin)-(malmbe+naddfno)+1)
          jj=ieos(nant(ispin)-(malmbe+naddfno)+0)
          difel = elecord(ii)-elecord(jj)
          rsigma(ispin) = 100.0*min(1d0,max(0.0d+00,0.5d0+difel))
          write (lw,*) '# +++++++++++++++++++++++++++++++++++++++++++++'
        endif

        deallocate (elecord)
        deallocate (ieos)
        deallocate (MOinG)
        deallocate (iMOinG)
        deallocate (listG)
      enddo
c
c-----Effective oxidation state of the fragment
c
      rsigmin=min(rsigma(1),rsigma(2))
      do i=1,nfrag
        eos(i)=zcharge(i)-eleca(i)-elecb(i)
      enddo
      write (lw,*) '#'
      write (lw,*) '#  ++++++++++++++++++++++++++++++++++'
      write (lw,*) '#  OXIDATION STATES OF ALL THE GROUPS'
      write (lw,*) '#  ++++++++++++++++++++++++++++++++++'
      write (lw,*) '#'
      write (lw,51) (i,eos(i),i=1,nfrag)
      write (lw,52) rsigma(1),rsigma(2),rsigmin

      deallocate (elec)
      deallocate (nelec)
      deallocate (zcharge)

      allocate (savemo(nmo*nfrag,nspin))
      savemo(1:nmo*nfrag,1:2) = .true.
      allocate (similarFNO(1:nmo*nfrag))
      
c
c-----Analyze set of quasi-similar normalized in R^3 FNOs
c
      write (lw,335) critoverlap
      do ispin=1,nspin
        if (ispin.eq.1) nab = nalpha
        if (ispin.eq.2) nab = nbeta
        if (ispin.eq.1) write (lw,600) 'ALPHA block'
        if (ispin.eq.2) write (lw,600) 'BETA  block'
        ocupt(:,ispin) = 0.0d+00
        iblock = 0
        do i=1,nant(ispin)
          if (savemo(i,ispin)) then
            iblock = iblock + 1
            savemo(i,ispin) = .false.
            nsim = 1
            similarFNO(nsim) = i
            do j=i+1,nant(ispin)
              overlap = dot_product(fno(:,i,ispin),fno(:,j,ispin))
              if (overlap.gt.critoverlap) then 
                savemo(j,ispin) = .false.
                nsim = nsim + 1
                similarFNO(nsim) = j
              endif
            enddo
            if (verbose.and.nsim.gt.1) write (lw,334)
            write (lw,331) iblock,nsim,(similarFNO(j),j=1,nsim)
c
c-----------Overlaps between quasi-equal FNOs
c
            do k=1,nsim
              kk = similarFNO(k)
              ocupt(iblock,ispin) = ocupt(iblock,ispin) 
     &                            + elecsave(kk,ispin)
              do l=1,nsim
                ll = similarFNO(l)
                ov(k,l) = dot_product(fno(:,kk,ispin),fno(:,ll,ispin))
              enddo
            enddo
            if (verbose) then
              if (nsim.gt.1) then
                write (lw,*) '# --------------'
                write (lw,*) '# Overlap matrix'
                write (lw,*) '# --------------'
                do k=1,nsim
                  write (lw,'(6(1x,F17.10))') (ov(k,l),l=1,nsim)
                enddo
              endif
            endif
          endif
        enddo
*       write (lw,*) '#'
*       write (lw,*) '# Total population of blocks'
*       write (lw,'(5(1x,F17.10))') (ocupt(ll,ispin),ll=1,iblock)
*       write (lw,'(1x,F17.10,a)') sum(ocupt(1:iblock,ispin)),' = SUM'
      enddo

      deallocate (iatom)
      deallocate (elecsave)
      deallocate (savemo)
      deallocate (fno)
      deallocate (ocupt)
      deallocate (similarFNO)
      write (lw,102)
      call timer (4,ineweos,'_neweos   ',-1)
      return
c
c.....Formats
c
 410  format (' # ',I4)
 411  format (1000(17x,10(1x,F6.2,'%'),/))
 68   format (/' # Localization of each ',a,' in each atom')
 711  format (' # SUM  ',17x,'= ',F18.8)
 351  format (1x,80('+'),/,1x,"+ File '",a,"'",
     & " contains the Normalized in R^3 FNOs",/,1x,80('+'))
 3510  format (1x,80('+'),/,1x,"+ File '",a,"'",
     & " contains all different Normalized in R^3 FNOs",/,1x,80('+'))
 353  format (1x,"# Writing file '",a,"'",1x,"with AOM for FNOs")
 1124 format (' # Eigenvalue ',I4,' of sg() is negative, ',E17.10)
 28   format (' # The value of EPSNEG is ',E15.8)
 222  format (' # Negative eigenvalue ',E15.8,
     &  ' set to 0.0 for fragment ',100I4)
 101  format (1000(' # ',5(1x,F15.8),/))
 80   format (6(1x,e16.10))
 3340 format (/,' # Second & fourth moment orbital spread of ',a)
 520  format (' # Atoms expanded by each FNO',/,1000(5(1x,F17.10),/))
 987  format (//' # ',77('+'),/,
     & ' # FNO ANALYSIS - FNO ANALYSIS - FNO ANALYSIS - FNO ANALYSIS',/
     & ' # FNOs are not normalized in R^3',/,
     & ' # NFNOs ARE normalized in R^3',/,' # ',77('+'),/,
     & ' # CURRENT Fragment (A in the following) formed by atoms:',/,
     & 1000(' # ',20(1x,I3),/))
 988  format (' Eigenvalues of S^A * 1RDM * S^A by decreasing order',/,
     &   1000(5(1x,F17.10),/))
 994  format (' # Full set of FNOs TOTAL POPULATION = ',F15.8,/,
     & ' # ',77('-'))
 991  format (' # ',77('-'),/,
     & ' # Fragment Natural Orbitals (FNO) occupations',/,' # ',
     & 77('-'),/,1000(5(1x,F17.10),/))
 989  format (' # W = Eigenvectors of S^A * 1RDM * S^A')
 990  format (' # Eigenvector ',I3,/,1000(5(1x,F17.10),/))
 42   format (' # Eigenvalues    > ',1PE15.8,' are ',I4)
 43   format (' # EPSEIGEN value = ',1PE15.8,/,
     &  ' # Number of Eigenvalues of S^A > EPSEIGEN = ',I4)
 332  format (' # Diagonal Overlaps in the fragment'
     & ' of normalized in R^3 FNOs',/,1000(5(1x,F17.10),/))
 3320 format (' # Ordered eigenvalues of S^A:',/,1000(5(1x,F17.10),/))
 34   format (' # ',a,' FNOs are determined')
 37   format (' # AOM in R^3 of the normalized in R^3 FNOs (NFNOs)',/,
     &        ' # i.e. AOM in R^3 of NFNOs')
 4444 format (5(1x,F17.10))
 4445 format (1x,F17.10,
     &  '  <-- Total Electron Population',/,' # ',117('-'))
 77   format (1000(2x,5(1x,F17.10),/))
 78   format (1000(2x,10(1x,I4),/))
 79   format (10x,F17.10,4x,I4,3x,I4)
 82   format (' # Oxidation state of fragment ',I2,' is ',F8.0)
 820  format (' # Fragment ',I2,' EOS = ',I2,
     & ' Reliability index = ',F9.2,' %')
 83   format (' # Reliability index = ',F12.5,' %')
 161  format (' #  MO\ATOM   ',10I9)
 412  format (' # ',I4,10x,10(1x,F7.2,'%'))
 770  format (1000(1x,5(1x,F13.10,' ( ',I2,' ) '),/))
 771  format (1000(1x,3(1x,F13.10,' ( ',I4,1x,I6,' ) '),/))
 772  format (1000(1x,1(1x,F13.10,' ( ',I4,' ) '),/))
 24   format (I7,2x,F17.10,2x,I4,6x,I4)
 240  format (' #',/,' # Ordered largest electron populations'
     & ,/,' # Number     Eigenvalue   Group    Original order',/,
     &    ' # -----------------------------------------------')
 111  format (//,' # ',40('+'),/,' # Oxidation States Analysis',/,
     &        ' # ',40('+'),/)
 51   format (' #  Group ',I2,4x,'EOS = ',F15.8)
 52   format (1x,'#  Reliability indices (R_alpha, R_beta, R)% = ',
     &   3(1x,F9.4))
 100  format (//,
     & ' ',/,
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' +                                                         +',/,
     & ' +      B E G I N   (OQS)  & (EOS)  A N A L Y S I S        +',/,
     & ' +                                                         +',/,
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' ')
 102  format (//,
     & ' ',/,
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' +                                                         +',/,
     & ' +         E N D   (OQS)  & (EOS)  A N A L Y S I S         +',/,
     & ' +                                                         +',/,
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' ')
 331  format (' # Block ',I3,' formed by ',I4,' FNOs:',1000(1x,I4))
 334  format (' #',/,' # ',80('+'))
 335  format (' #',/,' # ',80('+'),/,
     & ' # Sets of quasi-equal FNOs (i.e. with overlap close to 1.0)',/,
     & ' # Critital Overlap is ',F12.5,/,
     &              ' # ',80('+'))
 1617 format (' # Occup ',F13.10,' Degeneration = ',I4,
     &    1x,'FNOs --> ',1000I3)
 98   format (' # ',84('-'),/,
     & ' # Conflicting EOS assigmment: Two equivalent fragments',
     &   ' have a different number of FNOs')
 99   format (' # Group ',I2,' has the MOs ',1000I4)
 300  format (' # Full set of electron populations',/,
     & ' # In parenthesis the fragment to which it belongs')
 1004 format (' # Final list of fragments involved in the conflict: ',
     &    1000(20I4,/))
 1005 format (' # The fragment with more FNOs is the fragment ',I4,
     &  ' with ',I4,' FNOs')  
 1006 format (' # Total Number of Conflicting FNOs is ',I4)
 1007 format (' # Number of new FNOs added to the list ',I4)
 600  format (' #',/,' # ',11('-'),/,' # ',a,/,' # ',11('-'))
c
      end
