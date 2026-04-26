c
c-----------------------------------------------------------------------
c
      subroutine otheroqs (aom,vrot,nloc,nelhalf,udat,epsneg,
     &   wfnfile,largwr,warn,lw,lerr)
c
c-----PARAMETERS--------------------------------------------------------
c
c     aom()                          INPUT/OUTPUT
c
c     Atomic Overlap Matrix (AOM) of Canonical MOs in all the centers.
c
c.....vrot()                          OUTPUT 
c
c     Rotation matrix giving the localized MOs in terms of the input MOs. 
c
c.....nloc                           OUTPUT
c
c     Number of localized MOs found. It should be equal to nelhalf
c     (see below) if the routine ends nicely
c
c.....nelhalf                        INPUT
c
c     Half the number of electrons of the molecule (assumed to have
c     a even number of electrons, i.e. nalpha = nbeta)
c
c.....udat                           INPUT
c
c     Unit number of original WFN file
c
c-----epsneg                         INPUT
c
c     Diagonalizations in this routine involve definite positive 
c     matrices that should have all of their eigvenvalues positive.
c     Numerical errors, however, make that this not necessarily happens 
c     and some of the eigenvalues can be marginally negative. The value 
c     of 'epsneg' is the largest value (in absolute value) allowed in order
c     that the method ends nicely. For instance, if epsneg=-1D-6
c     and all the eigenvalues are greater that this number the 
c     process continues up to the end. However, if an eigvenvalue
c     is smaller that -1D-6, for instance -1D-5, this routine stops.
c
c-----wfnfile                        INPUT
c
c     Name of the WFN file
c
c-----largwr                         INPUT
c
c     Logical variable: if largwr=.true. a large output is requested,
c     otherwise a short output is used.
c
c-----warn                           OUTPUT
c
c     Returned  .false. if the routine ends nicely. Otherwise warn is
c     returned .true.
c
c-----lw                             INPUT
c
c     Logical output unit
c
c-----lerr                           INPUT
c
c     Logical output unit of error messages
c
c-----------------------------------------------------------------------
c                         
      USE      space_for_wfncoef
      USE      space_for_wfnbasis
      USE      space_for_rdm1
      USE      space_for_dafh
      include 'implicit.inc'
      include 'constants.inc'
      include 'param.inc'
      include 'wfn.inc'
      include 'error.inc'
      include 'mline.inc'
      real(kind=8)  vrot(nmo,nmo)
      real(kind=8)  aom(ncent,nmo,nmo)
      real(kind=8)  sg(nmo,nmo),eigen(nmo),uvec(nmo,nmo)
      real(kind=8), allocatable, dimension (:,:)   :: rho
      real(kind=8), allocatable, dimension (:)     :: eigs,work
      real(kind=8), allocatable, dimension (:)     :: eigloc
      real(kind=8), allocatable, dimension (:)     :: uveceig
      real(kind=8), allocatable, dimension (:,:)   :: v
      real(kind=8), allocatable, dimension (:,:)   :: eta
      real(kind=8), allocatable, dimension (:,:,:) :: aom0,aomx
      real(kind=8), allocatable, dimension (:,:)   :: sloc,sloc0
      real(kind=8), allocatable, dimension (:,:)   :: ctmp
      integer(kind=4),allocatable, dimension (:)   :: nmoc
c
c-----Data used in the orthonormalization process
c
      real(kind=8), allocatable,dimension  (:,:)   :: cc
      real(kind=8)  epsneg
      logical       warn,warno
      integer(kind=4) p,q,udat,udatnw
      logical   largwr
      character(len=*) wfnfile
      character(len=mline) wfnloc,locform(nmo)
      character(len=1) dp

      character*1    digs(0:9)
      data digs / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
c
c-----------------------------------------------------------------------
c
      call timer (2,iiloc,'_itotherlo',-1)

      allocate (rho(nmo,nmo))
      if (cciqa) then
        rho(1:nmo,1:nmo)=c1et(1:nmo,1:nmo)*0.5d0
      else
c
c-------Read 1-RDM in the canonical MO basis
c
        lu1 = 55
        open (lu1,file=trim(wfnfile)//".1rdm",form='unformatted',
     &        iostat=ios,status='old')
        if (ios.ne.0) write (lerr,100) trim(wfnfile)//".1rdm"
        if (ios.ne.0) stop
        read (lu1) nmox
        if (nmox.ne.nmo) then
          write (lerr,1123) trim(wfnfile)//".1rdm",nmox,nmo
          stop
        endif
        rho = 0d0
        do
          read (lu1,iostat=ios) m,n,value
          if (ios.ne.0) exit
          rho(m,n) = value
          rho(n,m) = value
        enddo
        rho = rho * 0.5d0
      endif
      write (lw,10) ' B E G I N  '
      write (lw,9) ndafh
      do i=1,nrhopen
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
      mxcent = maxval(nbcen(1:ndafh))
      allocate (eigs(n)            )
      allocate (eigloc(n)          )
      allocate (work(3*n-1)        )
      allocate (v(n,n)             )
      allocate (ctmp(n,n)          )
      allocate (aom0(ncent,n,n)    )
      allocate (aomx(ncent,n,n)    )
      allocate (uveceig(n)         )
      allocate (nmoc(mxcent)       )
c
c-----Print the input AOM arrays
c
      if (largwr) then
        write (lw,*) '# Atomic Overlap matrix (AOM)'
        write (lw,*) '#'
        do i=1,ncent
          write (lw,*) '     ATOM ',i
          do j=1,n
            write (lw,'(6(1x,F15.8))') (aom(i,j,k),k=1,j)
          enddo
        enddo
        write (lw,*) '#'
        write (lw,*) '#    TOTAL'
        do j=1,n
          write (lw,'(6(1x,F15.8))') (sum(aom(:,j,k)),k=1,j)
        enddo
      endif
c
c.....Start the process.
c
      sloc  = zero
      sloc0 = zero
      ctmp     = zero
      n3       = 3*n-1
      nloc     = 0
      nmoc     = 0
      warn = .false.
c
c-----Group overlap integrals are obtained from the AOM
c
      write (lw,28) epsneg
      do nd=1,ndafh
        sg=zero
        do j=1,nbcen(nd)
          ib=ibcen(nd,j)
          sg=sg+aom(ib,:,:)
        enddo
        call dsyev ('V','U',nmo,sg,nmo,eigen(1:nmo),work,n3,info)
        call errdsyev (info,'otheroqs.f')
        uvec=sg
c
c-------Test that sg() is definite positive
c
        do m=1,nmo
          if (eigen(m).gt.0d0) then
          elseif (eigen(m).gt.epsneg) then
            write (lerr,222) eigen(m),(ibcen(nd,j),j=1,nbcen(nd))
            eigen(m)=0d0
          else
            write (lerr,1124) m,eigen(m)
            stop
          endif
        enddo
        sg = 0d0   ! This is really S^{1/2}
        do k=1,nmo
          do l=1,nmo
            forall (m=1:nmo) uveceig(m)=uvec(l,m)*sqrt(eigen(m))
            sg(k,l) = dot_product(uvec(k,:),uveceig(:))
          enddo
        enddo
c
c-------This is S^{1/2} * 1-RDM * S^{1/2}. The previous "localized"
c       density is substracted
c
        v = matmul(sg,matmul(rho,transpose(sg))) - sloc0
        eigs = 0d0
        call dsyev ('V','U',n,v,n,eigs(1:n),work,n3,info)
        call errdsyev (info,'otheroqs.f')
        if (largwr) then
          write (lw,*) '#'
          write (lw,432) (ibcen(nd,j),j=1,nbcen(nd))
          write (lw,*) '#'
          write (lw,'(6(1x,F15.8))') (eigs(k),k=1,n)
        endif
c
c-------Be carefull with the following loop. It only works if the 
c       used diagonalization routine (in this case, 'dsyev') 
c       returns the eigenvalues in ascencing order. Otherwise
c       before this loop eigenvalues and eigvectors should be
c       reordered such that eigs(1) is the smallest eigvenvalue
c       eigs(2) the following, and so on.
c
        do i=n,1,-1
          if (eigs(i).lt.cutoff(nd)) exit
          nloc = nloc + 1
          nmoc(nbcen(nd))=nmoc(nbcen(nd))+1
          if (nloc.gt.nmo) then
            write (lw,*) '#'
            write (lw,*) '# Too many localized MOs'
            write (lw,*) '#'
            stop              
          endif
          eigloc(nloc)=eigs(i)
          ctmp(:,nloc)=v(:,i)
*         sumcent(nbcen(nd))=sumcent(nbcen(nd))+eigs(i)
          do j=1,n
            do k=1,n
              prod=ctmp(j,nloc)*ctmp(k,nloc)*eigloc(nloc)
              sloc(j,k)=sloc(j,k)+prod
            enddo
          enddo
          if (nloc.le.nelhalf) then
            write (lw,116) nloc,eigs(i),
     &        (atnam(ibcen(nd,j))(3:4),ibcen(nd,j),j=1,nbcen(nd))
          else
            write (locform(nloc-nelhalf),116) nloc,eigs(i),
     &        (atnam(ibcen(nd,j))(3:4),ibcen(nd,j),j=1,nbcen(nd))
          endif
        enddo
        if (deplet(nd)) sloc0=sloc
      enddo

      if (nloc.lt.nelhalf) then
        write (lw,302)
        warn = .true.
      elseif (nloc.eq.nelhalf) then
        write (lw,300) (nmoc(i),i=1,mxcent)
      else
        write (lw,435)
        do i=nelhalf+1,nloc
          write (lw,'(a)') trim(locform(i-nelhalf))
        enddo
      endif
      vrot(1:nloc,:) = transpose(ctmp(:,1:nloc))
c
c-----Each OPENLOC localized MO (LMO) is normalized to 1.0 in the
c     fragment it was obtained, but is normalized to eigloc^{-1} in R^3. 
c     For this reason each LMO will be multiplied by sqrt(eigloc) such
c     that it will be normalized in R^3 instead of the fragment.
c
      forall (i=1:nloc) vrot(i,1:n)=vrot(i,1:n)*sqrt(eigloc(i))
      if (largwr) then
        write (lw,20)
        do i=1,nloc
          write (lw,33) i
          write (lw,30) (vrot(i,j),j=1,n)
        enddo
      endif
c
c-----write a WFN file for 1c, 2c, 3c, ... NON-orthonormal molecular orbitals
c
      allocate (cc(nloc,nprims))
      do i=1,nloc
        do ip=1,nprims
          cc(i,ip) = dot_product(coef(nmo+1:nmo+nmo,ip),vrot(i,:))
        enddo
      enddo
      wfnloc = trim(wfnfile)//"open"
      udatnw = udat+10
      inic=1
      nfin=0
      do i=1,mxcent
        if (nmoc(i).gt.0) then
          nfin=nfin+nmoc(i)
          open (unit=udatnw,file=trim(wfnloc)//'-'//digs(i)//'center')
          write (lw,353) trim(wfnloc)//'-'//digs(i),i
          call cdafhmos (udat,udatnw,0,cc(inic:nfin,:),nmoc(i),nprims)
          inic=inic+nmoc(i)
        endif
      enddo
      deallocate (cc)
c
c-----Write diagonal AOM elements in each atom or fragment
c
*     sumeigen = sum(eigloc(1:nloc))
      do i=1,ncent
        aom0(i,1:nloc,1:nloc) = matmul(vrot(1:nloc,:),
     &  matmul(aom(i,:,:),transpose(vrot(1:nloc,:))))
      enddo

      aom(:,1:nloc,1:nloc) = aom0(:,1:nloc,1:nloc)

      write (lw,68)
      write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
      do ip=1,nloc
        write (lw,410,advance='no') ip
        if (ncent.gt.10) then
          write (lw,41 ) (100d0*aom0(k,ip,ip),k=1,10)
          write (lw,411) (100d0*aom0(k,ip,ip),k=11,ncent)
          write (lw,*) 
        else
          write (lw,41) (100d0*aom0(k,ip,ip),k=1,ncent)
        endif
      enddo
      write (lw,112) 'BEGINNING'
*     goto 118
      allocate (cc(nloc+nloc,nprims))
      do i=1,nloc
        do ip=1,nprims
          cc(i,ip) = dot_product(coef(nmo+1:nmo+nmo,ip),vrot(i,:))
          cc(i+nloc,ip) = cc(i,ip)
        enddo
      enddo
      write (lw,1000) 'Open System Localized MOs '
      allocate (eta(nloc,nloc))
      call ofmortho (cc,eta,nloc,nprims,warno,lw)
      if (warno) stop ' # otheroqs.f: Singular matrix in ofmortho.f'
      do i=1,ncent
        aomx(i,1:nloc,1:nloc) = matmul(eta,
     &  matmul(aom(i,1:nloc,1:nloc),transpose(eta)))
      enddo
      write (lw,68)
      write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
      do ip=1,nloc
        write (lw,410,advance='no') ip
        if (ncent.gt.10) then
          write (lw,41) (100d0*aomx(k,ip,ip),k=1,10)
          write (lw,411) (100d0*aomx(k,ip,ip),k=11,ncent)
          write (lw,*) 
        else
          write (lw,41) (100d0*aomx(k,ip,ip),k=1,ncent)
        endif
      enddo
      write (lw,334) 'orthonormalized Localized MOs'
      call newdrvmult (nloc,nprims,ncent,cc(1:nloc,:),lw)
      inpr=nprims
      occ(1:nloc)=zero
      wfnloc=trim(wfnfile)//"open-ortho"    
      udatnw=udat+10
      open (unit=udatnw,file=trim(wfnloc),status='unknown')
      write (lw,351) trim(wfnloc)
      call cdafhmos (udat,udatnw,0,cc(1:nloc,:),nloc,inpr)
c
c-----write a WFN file for 1c, 2c, 3c, ... orthonormal molecular orbitals
c
      inic=1
      nfin=0
      do i=1,mxcent
        if (nmoc(i).gt.0) then
          nfin=nfin+nmoc(i)
          open (unit=udatnw,file=trim(wfnloc)//'-'//digs(i)//'center')
          write (lw,352) trim(wfnloc)//'-'//digs(i),i
          call cdafhmos (udat,udatnw,0,cc(inic:nfin,:),nmoc(i),nprims)
          inic=inic+nmoc(i)
        endif
      enddo
      deallocate (cc)

      write (lw,112) 'END'
c
      deallocate (rho        )
      deallocate (eta        )
      deallocate (aom0       )
      deallocate (aomx       )
      deallocate (eigs       )
      deallocate (work       )
      deallocate (v          )
      deallocate (ctmp       )
      deallocate (eigloc     )
      deallocate (cc         )
*118  continue
      write (lw,10) '  E N D  '
      return
      call timer (4,iiloc,'_itotherlo',-1)
c
c.....Formats
c
 9    format (//,' # ',20x,
     & 'Definition of the OTHEROPEN order',/,' # ',80('-'),/,
     &  ' # Number of bonds to analyze = ',I3,/,' # ',80('-'))
 105   format (' # Bond ',I5,1x,'Cutoff = ',F7.5,1x,
     & '(deplet rho? = ',a,'): ',I2,' centers = ',1000I3)
 115   format (' # ',80('-'))
 410  format (' # ',I4)
 20   format (' # Rotation Matrix (READ BY ROWS)',6x,
     & 'Open System MO_i =  SUM_j CANMO_j C_ji',/,' # ',
     & 'Open System MOs have been renormalized such that <i|i>_R^3=1')
 30   format (8(1x,F15.8))
 41   format (1000(10x,10(1x,F5.1,'%')))
 411  format (1000(17x,10(1x,F5.1,'%'),/))
 300  format (' #',/,' # ',100('-'),/,
     &  ' # CONVERGED (1--N)-CENTER OTHEROPEN ANALYSIS',/,' #',/,
     &  ' # Number of MOs involving 1, 2, 3, ... centers',/,
     &  ' # ',20(1x,I3))
 302  format (' #',/,' # otheroqs.f: Unable to localize all the MOs',
     & /,' #',/,' # ',100('-'),/,
     & ' # NON CONVERGED (1--N)-CENTER OPENLOC ANALYSIS',/,' #')
 66   format (/' # Diagonal overlaps of each MO in each atom')
 68   format (/' # Degree of localization of each MO in each atom')
 33   format (' # MO ',I3)
 10   format (/,' # ',100('-'),/,' # ',a,
     & '(1--N)-C E N T E R   D A F H   A N A L Y S I S',/,
     & ' # ',100('-'))
 11   format (' #',/,' # OPENLOC analysis on ',I2,' centers',/,' #')
 69   format (' #',/,' # SUM of eigenvalues = ',F18.8,/,' #')
 70   format (' # ',7x,I4,'-center = ',F18.8)
 71   format (' # ',/,
     & ' # Electrons associated to successful OPENLOC diagonalizations')
 710  format (' # Atom ',I4,'   Electrons = ',F18.8)
 711  format (' # SUM  ',17x,'= ',F18.8)
 334  format (/,' # Second & fourth moment orbital spread of ',a,/,3x,
     &'I.M. Hoyvik & P. Jorgensen [Chem. Rev. 116, 3306-3326 (2016)]',
     &/,3x,
     & 'I.M. Hoyvik, B. Jansik & P. Jorgensen [JCP 137, 2224114 (2012)'
     & /,3x,61('-'),/,4x,'Orbital',1x,'sigma_2^p',5x,'sigma_4^p',
     & 5x,'beta_p',8x,'beta_p^4',/,3x,61('-'))
 351  format (/1x,"# Writing file ","'",a,"' with orthonormalized MOs")
 1000 format (' # ',a,'are orthonormalized using the Lowdin method')
 112  format (2(' #',/),' #',90('-'),/,' # ',a,
     & ' OF RESULTS WITH ORTHONORMALIZED OPTIMIZED MOs',/,' #',90('-'))
 51   format (' # In single atom analyses, H atoms are',a,'skipped')
 52   format (
     & " # In 'atom1-atom2' searches, maximum distance index is ",I2)
 432  format (' # Eigenvalues of S^1/2*RHO*S^1/2 of atoms ',20(1x,I3))
 113  format (1x,'# ',I2,
     & '-center OPENLOC analysis, Critical overlap & atoms = ',F12.5,
     & 4x,20I3,/,' # ',78('-'))
 116  format (4x,'MO ',I4,' with eigenvalue = ',F15.8,
     & ') localized on atom(s)',20(1x,'(',A2,I2,')'))
 100  format (' # itotherrho.f: File ',a,' NOT FOUND')
 1123 format (' # NMO in file ',a,'(',I4,') # Current value = ',I4)
 1124 format (' # Eigenvalue ',I4,' of sg() is negative, ',E17.10)
 28   format (' # The value of EPSNEG is ',E15.8)
 222  format (' # Negative eigenvalue ',E15.8,
     &    ' set to 0.0 diagonalizing fragment ',100I4)
 435  format (' #',
     & /,' # otheroqs.f: ! WARNING: Too many localized MOs !',/,' #')
 352  format (" # Writing file '",a,"' with ",I2,
     &  "-center ortho-normalized MOs")
 353  format (" # Writing file '",a,"' with ",I2,
     &  "-center non-orthonormal MOs")
      end
