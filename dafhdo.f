c
c-----------------------------------------------------------------------
c
      subroutine dafhdo (aom,dafh,vrot,nloc,nelhalf,udat,
     &   wfnfile,largwr,warndafh,lw)
c
c.....Performs the localization of a set of NMO molecular orbitals into 
c     NCENT atoms by diagonalizing the AOM of single atoms, pairs, trios, 
c     ... of atoms. After each diagonalization, the transformed MO
c     sufficiently localized are excluded from the total DAH, if desired.
c                            
c                            
c-----PARAMETERS--------------------------------------------------------
c
c-----aom()                            INPUT/OUTPUT
c
c     Atomic overlap matrix between the input or output MOs in all the
c     center of the molecule.
c
c.....dafh()                           INPUT
c
c     Domain Averaged Fermi Hole in all the centers of the molecule.
c
c.....vrot()                           OUTPUT 
c
c     Rotation matrix giving the transformed MO's from the input MO's. 
c
c.....nloc                             OUTPUT
c
c     Number of Localized MOs found.
c
c.....nelhalf                          INPUT
c
c     Half the number of electrons of the molecule.
c
c-----udat                             INPUT
c
c     Logical unit of the file containing the WFN of the molecule
c
c.....wfnfile                          INPUT
c
c     Name of the file with the wavefunction of the molecule
c
c.....largwr                           INPUT
c
c     .TRUE. if large output is requested, .FALSE. otherwise
c
c.....warndafh                         OUTPUT
c
c     This logical variable will be .TRUE. if something goes wrong
c
c.....lw                               INPUT
c
c     Output logical unit
c-----------------------------------------------------------------------
c                         
      USE      space_for_wfncoef
      USE      space_for_wfnbasis
      USE      space_for_dafh
      include 'implicit.inc'
      include 'constants.inc'
      include 'param.inc'
      include 'wfn.inc'
      include 'error.inc'
      include 'mline.inc'
      real(kind=8)  aom(ncent,nmo,nmo),dafh(ncent,nmo,nmo),vrot(nmo,nmo)
      real(kind=8), allocatable, dimension (:)     :: eigs,work,eigloc
      real(kind=8), allocatable, dimension (:)     :: sumcent,elec
      real(kind=8), allocatable, dimension (:,:)   :: v
      real(kind=8), allocatable, dimension (:,:)   :: aomt,dafht
      real(kind=8), allocatable, dimension (:,:)   :: eta
      real(kind=8), allocatable, dimension (:,:,:) :: aom0,sgx,dafh0
      real(kind=8), allocatable, dimension (:,:)   :: dafhloc,dafhloc0
      real(kind=8), allocatable, dimension (:,:)   :: ctmp
      real(kind=8), allocatable, dimension  (:,:)  :: cc
      real(kind=8), allocatable, dimension (:)     :: ocup
      logical         warndafh,warnofmo
      integer(kind=4) p,q,udat,udatnw
      logical   largwr
      character(len=*) wfnfile
*     parameter   (mline   = 200)
      character(len=mline) wfnloc,locform(nmo)
      character(len=1) dp
c
      call timer (2,idafhdo,'_dafhdo   ',-1)
      write (lw,10) ' B E G I N  '
      write (lw,9) ndafh
      do i=1,ndafh
        if (deplet(i)) then
          dp='T'
        else
          dp='F'
        endif
        write (lw,105) i,cutoff(i),dp,nbcen(i),(ibcen(i,j),j=1,nbcen(i))
      enddo
      write (lw,115)
c
c.....Allocate arrays
c
      n = nmo
      allocate (eigs(n))
      allocate (eigloc(n))
      allocate (work(3*n-1))
      allocate (v(n,n))
      allocate (ctmp(n,n))
      allocate (dafh0(ncent,n,n))
      allocate (aom0(ncent,n,n))
      allocate (sgx(ncent,n,n))
      allocate (dafhloc(n,n))
      allocate (dafhloc0(n,n))
      allocate (sumcent(ncent))
      allocate (elec(ncent))
c
c-----Print the input DAFH and AOM arrays
c
      if (largwr) then
        allocate (aomt(n,n))
        allocate (dafht(n,n))
        aomt=zero
        write (lw,*) '# Atomic Overlap matrix (AOM)'
        write (lw,*) '#'
        do i=1,ncent
          write (lw,*) '     ATOM ',i
          do j=1,n
            do k=1,n
              aomt(j,k)=aomt(j,k)+aom(i,j,k)
            enddo
            write (lw,'(6(1x,F15.8))') (aom(i,j,k),k=1,j)
          enddo
        enddo
        write (lw,*) '#'
        write (lw,*) '#    TOTAL'
        do j=1,n
          write (lw,'(6(1x,F15.8))') (aomt(j,k),k=1,j)
        enddo
        write (lw,*) '#'
        dafht=zero
        write (lw,*) '#'
        write (lw,*) '# Domain Averaged Fermi Hole (DAFH)'
        write (lw,*) '#'
        do i=1,ncent
          write (lw,*) '     ATOM ',i
          do j=1,n
            do k=1,n
              dafht(j,k)=dafht(j,k)+dafh(i,j,k)
            enddo
            write (lw,'(6(1x,F15.8))') (dafh(i,j,k),k=1,j)
          enddo
        enddo
        write (lw,*) '#'
        write (lw,*) '#    TOTAL'
        do j=1,n
          write (lw,'(6(1x,F15.8))') (dafht(j,k),k=1,j)
        enddo
        deallocate (aomt,dafht)
      endif
c
c.....Start the process.
c
      dafh0    = dafh 
      dafhloc  = zero
      dafhloc0 = zero
      ctmp     = zero
      n3       = 3*n-1
      nloc     = 0
      sumcent  = zero
      warndafh = .false.

      do nd=1,ndafh
        v=zero
        do j=1,nbcen(nd)
          ib=ibcen(nd,j)
          v=v+dafh0(ib,:,:)
        enddo
        v=v-dafhloc0
        eigs = 0d0
        call dsyev ('V','U',n,v,n,eigs(1:n),work,n3,info)
        call errdsyev (info,'dafhdo.f')
        if (largwr) then
          write (lw,*) '#'
          write (lw,432) (ibcen(nd,j),j=1,nbcen(nd))
          write (lw,*) '#'
          write (lw,'(6(1x,F15.8))') (eigs(k),k=1,n)
        endif
        do i=n,1,-1
          if (eigs(i).lt.cutoff(nd)) then
            exit
          else
            nloc = nloc + 1
            if (nloc.gt.nmo) then
              write (lw,*) '#'
              write (lw,*) '# Too many localized DAFH MOs'
              write (lw,*) '#'
              stop              
            endif
            eigloc(nloc)=eigs(i)
            ctmp(:,nloc)=v(:,i)
            sumcent(nbcen(nd))=sumcent(nbcen(nd))+eigs(i)
            do j=1,n
              do k=1,n
                prod=ctmp(j,nloc)*ctmp(k,nloc)*eigloc(nloc)
                dafhloc(j,k)=dafhloc(j,k)+prod
              enddo
            enddo
            if (nloc.le.nelhalf) then
              write (lw,116) nloc,eigs(i),
     &          (atnam(ibcen(nd,j))(3:4),ibcen(nd,j),j=1,nbcen(nd))
            else
              write (locform(nloc-nelhalf),116) nloc,eigs(i),
     &          (atnam(ibcen(nd,j))(3:4),ibcen(nd,j),j=1,nbcen(nd))
            endif
          endif
        enddo
        if (deplet(nd)) dafhloc0=dafhloc
      enddo

      if (nloc.lt.nelhalf) then
        write (lw,*)'# '
        write (lw,*)'# dafhdo.f: Unable to localize all the MOs'
        write (lw,*)'# '
        write (lw,302)
        warndafh = .true.
      elseif (nloc.eq.nelhalf) then
        write (lw,300)
      else
        write (lw,*)
        write (lw,*)'# dafhdo.f: ! WARNING: Too many localized MOs !'
        write (lw,*)'# '
        do i=nelhalf+1,nloc
          write (lw,'(a)') trim(locform(i-nelhalf))
        enddo
      endif
      vrot(1:nloc,:)=transpose(ctmp(:,1:nloc))
      if (largwr) then
        write (lw,20)
        do i=1,nloc
          write (lw,33) i
          write (lw,30) (vrot(i,j),j=1,n)
        enddo
      endif
c
c-----Write diagonal AOM elements in each atom or fragment
c
      sumeigen=sum(eigloc(1:nloc))
      do i=1,ncent
        aom0(i,1:nloc,1:nloc) = matmul(vrot(1:nloc,:),
     &  matmul(aom(i,:,:),transpose(vrot(1:nloc,:))))
      enddo
      write (lw,69) sumeigen
      do i=1,ncent
        if (abs(sumcent(i)).gt.zero) write (lw,70) i,sumcent(i)
      enddo
      write (lw,68)
      write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
      elec = zero
      do ip=1,nloc
        do j=1,ncent
          elec(j) = elec(j) + aom0(j,ip,ip)
        enddo
        write (lw,410,advance='no') ip
        if (ncent.gt.10) then
          write (lw,41 ) (100d0*aom0(k,ip,ip),k=1,10)
          write (lw,411) (100d0*aom0(k,ip,ip),k=11,ncent)
          write (lw,*) 
        else
          write (lw,41) (100d0*aom0(k,ip,ip),k=1,ncent)
        endif
      enddo
      aom(:,1:nloc,1:nloc)=aom0(:,1:nloc,1:nloc)
      write (lw,71) 
      do i=1,ncent
        write (lw,710) i,elec(i)
      enddo
      write (lw,711) sum(elec(1:ncent))
      write (lw,112) 'BEGINNING'
      goto 118
      if (.not.allocated(cc)) then
        allocate (cc(nloc+nloc,nprims),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate cc()'
      endif
      do i=1,nloc
        do ip=1,nprims
          tmp=zero
          do j=1,nmo
            tmp=tmp+coef(j+nmo,ip)*vrot(i,j)
          enddo
          cc(i,ip)=tmp
          cc(i+nloc,ip)=tmp
        enddo
      enddo
      write (lw,1000)
      if (.not.allocated(eta)) then
        allocate (eta(nloc,nloc),stat=ier)
        if (ier.ne.0) stop 'dafhdo.f: Cannot allocate eta()'
      endif
      call ofmortho (cc,eta,nloc,nprims,warnofmo,lw)
      if (warnofmo) stop ' # dafhdo.f: Singular matrix in ofmortho.f'
      do i=1,ncent
        sgx(i,1:nloc,1:nloc) = matmul(eta,
     &  matmul(aom(i,1:nloc,1:nloc),transpose(eta)))
      enddo
      write (lw,68)
      write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
      elec = zero
      do ip=1,nloc
        do j=1,ncent
          elec(j) = elec(j) + sgx(j,ip,ip)
        enddo
        write (lw,410,advance='no') ip
        if (ncent.gt.10) then
          write (lw,41) (100d0*sgx(k,ip,ip),k=1,10)
          write (lw,411) (100d0*sgx(k,ip,ip),k=11,ncent)
          write (lw,*) 
        else
          write (lw,41) (100d0*sgx(k,ip,ip),k=1,ncent)
        endif
      enddo
      write (lw,334) 'orthonormalized Localized MOs'
      call newdrvmult (nloc,nprims,ncent,cc(1:nloc,:),lw)
      inpr=nprims
      occ(1:nloc)=zero
      wfnloc=trim(wfnfile)//"dafh-ortho"    
      udatnw=udat+10
      open (unit=udatnw,file=trim(wfnloc),status='unknown')
      write (lw,351) trim(wfnloc)
      allocate (ocup(nloc),stat=ier)
      if (ier.ne.0) stop '# dafhdo.f: Cannot allocate ocup()'
      ocup(1:nloc)=occ(1:nloc)
      call cdafhmos (udat,udatnw,0,ocup,cc(1:nloc,:),nloc,inpr)
      write (lw,112) 'END'
      deallocate (ocup,stat=ier)
      if (ier.ne.0) stop '# dafhdo.f: Cannot deallocate ocup()'
      deallocate (eta)
      deallocate (aom0)
      deallocate (sgx)
      deallocate (dafh0)
      deallocate (eigs)
      deallocate (work)
      deallocate (v)
      deallocate (dafhloc)
      deallocate (dafhloc0)
      deallocate (ctmp)
      deallocate (eigloc)
      deallocate (sumcent)
      deallocate (elec)
 118  continue
      write (lw,10) '  E N D  '
      return
      call timer (4,idafhdo,'_dafhdo   ',-1)
c
c.....Formats
c
 9    format (//,' # ',20x,
     & 'Definition of the DODAFH order',/,' # ',80('-'),/,
     &  ' # Number of bonds to analyze = ',I3,/,' # ',80('-'))
 105   format (' # Bond ',I5,1x,'Cutoff = ',F7.5,1x,
     & '(deplet rho? = ',a,'): ',I2,' centers = ',1000I3)
 115   format (' # ',80('-'))
 410  format (' # ',I4)
 20   format (' # Rotation Matrix (READ BY ROWS)',
     &        '    DAFH MO_i =  SUM_j CANMO_j C_ji')
 30   format (6(1x,F12.6))
 41   format (1000(10x,10(1x,F5.1,'%')))
 411  format (1000(17x,10(1x,F5.1,'%'),/))
 300  format (' #',/,' # ',100('-'),/,
     &  ' # CONVERGED (1--N)-CENTER DAFH ANALYSIS',/,' #')
 302  format (' #',/,' # ',100('-'),/,
     &        ' # NON CONVERGED (1--N)-CENTER DAFH ANALYSIS',/,' #')
 66   format (/' # Diagonal overlaps of each MO in each atom')
 68   format (/' # Degree of localization of each MO in each atom')
 33   format (' # MO ',I3)
 10   format (/,' # ',100('-'),/,' # ',a,
     & '(1--N)-C E N T E R   D A F H   A N A L Y S I S',/,
     & ' # ',100('-'))
 11   format (' #',/,' # DAFH analysis on ',I2,' centers',/,' #')
 69   format (' #',/,' # SUM of eigenvalues = ',F18.8,/,' #')
 70   format (' # ',7x,I4,'-center = ',F18.8)
 71   format (' # ',/,
     &  ' # Electrons associated to successful DAFH diagonalizations')
 710  format (' # Atom ',I4,'   Electrons = ',F18.8)
 711  format (' # SUM  ',17x,'= ',F18.8)
 334  format (/,' # Second & fourth moment orbital spread of ',a,/,3x,
     &'I.M. Hoyvik & P. Jorgensen [Chem. Rev. 116, 3306-3326 (2016)]',
     &/,3x,
     & 'I.M. Hoyvik, B. Jansik & P. Jorgensen [JCP 137, 2224114 (2012)'
     & /,3x,61('-'),/,4x,'Orbital',1x,'sigma_2^p',5x,'sigma_4^p',
     & 5x,'beta_p',8x,'beta_p^4',/,3x,61('-'))
 351  format (/1x,"# Writing file ","'",a,"' with orthonormalized MOs")
 1000 format (' # DAFH MOs are orthonormalized using the Lowdin method')
 112  format (2(' #',/),' #',90('-'),/,' # ',a,
     & ' OF RESULTS WITH ORTHONORMALIZED OPTIMIZED MOs',/,' #',90('-'))
 51   format (' # In single atom analyses, H atoms are',a,'skipped')
 52   format (
     & " # In 'atom1-atom2' searches, maximum distance index is ",I2)
 432  format (' # Eigenvalues of the DAFH of atoms ',20(1x,I3))
 113  format (1x,'# ',I2,
     & '-center DAFH analysis, Critical overlap & atoms = ',F12.5,
     & 4x,20I3,/,' # ',78('-'))
 116  format (4x,'MO ',I4,' with eigenvalue = ',F15.8,
     & ') localized on atom(s)',20(1x,'(',A2,I2,')'))

      end
