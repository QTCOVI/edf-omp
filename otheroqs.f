c
c-----------------------------------------------------------------------
c
      subroutine otheroqs (aom,aomo,rho,cret,coret,icode,nloc,
     &  nloctarget,nab,udat,epsneg,mal,mbe,largwr,warn,lw,lerr)
c
c-----PARAMETERS--------------------------------------------------------
c
c     aom()                          INPUT/OUTPUT
c
c     Atomic Overlap Matrix (AOM) of Canonical MOs in all the centers.
c
c     aomo()                         OUTPUT
c
c     AOM between orthonormalized Localized MOs in all the centers in the
c     output
c
c-----rho()                          INPUT
c
c     First order RDM of alpha (icode=1) or beta (icode=2) electrons
c
c-----CRET()                         OUTPUT 
c
c     Rotation matrix giving the nonorthonormal localized MOs in terms of 
c     the input MOs: Loc_i = SUM_j C_ij CANMO_j
c
c-----CORET()                        OUTPUT 
c
c     Rotation matrix giving the orthonormal localized MOs in terms of 
c     the input MOs: Loc_i = SUM_j C_ij CANMO_j
c
c-----ICODE                          INPUT 
c
c     = 1 ---> alpha block
c     = 2 ---> beat  block
c
c.....nloc                           OUTPUT
c
c     Number of localized MOs found. It should be equal to nloctarget
c     (see below) if the routine ends nicely
c
c.....nloctarget                     INPUT
c
c     Expected number of Localized MOs that should be found in case 
c     that the method works fine.

c-----nab                            input
c
c     Number of MOs involved in the localization  (equal to the number
c     of alpha or beta MOs in the full set of NMO orbitals of the WFN
c     file)
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
c-----mal,mbe                        INPUT
c
c     Number of alpha and beta electrons
c
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
      real   (kind=8), dimension(nab,nab)  :: vrot
      real   (kind=8), dimension(nab,nab)  :: vroto
      real   (kind=8), dimension(nab,nab)  :: rho
      real   (kind=8), dimension(nab,nab)  :: sg
      real   (kind=8), dimension(nab)      :: eigen
      real   (kind=8), dimension(nab,nab)  :: uvec
      real   (kind=8), dimension(nmo)      :: xloci
      real   (kind=8), dimension (nab)     :: eigs
      real   (kind=8), dimension (3*nab-1) :: work
      real   (kind=8), dimension (nab)     :: eigloc
      real   (kind=8), dimension (nab)     :: uveceig
      real   (kind=8), dimension (nab,nab) :: v
      real   (kind=8), dimension (nab,nab) :: ctmp
      real   (kind=8), dimension(nloctarget,nab)        :: cret
      real   (kind=8), dimension(nloctarget,nab)        :: coret
      real   (kind=8), dimension(ncent,nab,nab)         :: aom
      real   (kind=8), dimension(ncent,nab,nab)         :: aom0
      real   (kind=8), dimension(ncent,nab,nab)         :: aomo
      real   (kind=8), dimension(nloctarget,nloctarget) :: sumaom
      integer(kind=4), dimension(ncent,ncent)           :: enlaces
      integer(kind=4) iab(2,nmo*(nmo-1)/2)
      integer(kind=4) bond(nmo)

      real   (kind=8), allocatable, dimension (:,:)   :: sloc,sloc0
      integer(kind=4),allocatable, dimension (:)      :: nmoc
      real(kind=8)  epsneg
      logical       warn,warno,largwr
      integer(kind=4) p,q,udat,udatnw
      character(len=mline) wfnloc,locform(nmo)
      character (len=30) typebond
      character(len=1) dp

      character*1    digs(0:9)
      data digs / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
c
c-----------------------------------------------------------------------
c
      call timer (2,iiloc,'_otheroqs ',-1)

      if (icode.eq.1) write (lw,10) ' B E G I N  ','ALPHA block'
      if (icode.eq.2) write (lw,10) ' B E G I N  ','BETA  block'
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
      allocate (nmoc(mxcent)        )
c
c-----Print the input AOM arrays
c
      if (largwr) then
        write (lw,*) '# Atomic Overlap matrix (AOM)'
        write (lw,*) '#'
        do i=1,ncent
          write (lw,*) '     ATOM ',i
          do j=1,nab
            write (lw,'(6(1x,F15.8))') (aom(i,j,k),k=1,j)
          enddo
        enddo
        write (lw,*) '#'
        write (lw,*) '#    TOTAL'
        do j=1,nab
          write (lw,'(6(1x,F15.8))') (sum(aom(:,j,k)),k=1,j)
        enddo
      endif
c
c.....Start the process.
c
      sloc  = 0d0
      sloc0 = 0d0
      ctmp  = 0d0
      n3    = 3*nab-1
      nloc  = 0
      nmoc  = 0
      warn  = .false.
      nbab  = 0
      enlaces = 0
c
c-----Group overlap integrals are obtained from the AOM
c
      write (lw,28) epsneg
      do nd=1,ndafh
        sg=0d0
        do j=1,nbcen(nd)
          ib=ibcen(nd,j)
          sg=sg+aom(ib,:,:)
        enddo
        call dsyev ('V','U',nab,sg,nab,eigen(1:nab),work,n3,info)
        call errdsyev (info,'otheroqs.f')
        uvec = sg
c
c-------Test that sg() is definite positive
c
        do m=1,nab
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
        do k=1,nab
          do l=1,nab
            forall (m=1:nab) uveceig(m)=uvec(l,m)*sqrt(eigen(m))
            sg(k,l) = dot_product(uvec(k,:),uveceig(:))
          enddo
        enddo
c
c-------This is S^{1/2} * 1-RDM * S^{1/2}. The previous "localized"
c       density is substracted
c
        v = matmul(sg,matmul(rho,transpose(sg))) - sloc0
        eigs = 0d0
        call dsyev ('V','U',nab,v,nab,eigs(1:nab),work,n3,info)
        call errdsyev (info,'otheroqs.f')
        if (largwr) then
          write (lw,*) '#'
          write (lw,432) (ibcen(nd,j),j=1,nbcen(nd))
          write (lw,*) '#'
          write (lw,'(6(1x,F15.8))') (eigs(k),k=1,nab)
        endif
c
c-------Be carefull with the following loop. It only works if the 
c       used diagonalization routine (in this case, 'dsyev') 
c       returns the eigenvalues in ascencing order. Otherwise
c       before this loop eigenvalues and eigvectors should be
c       reordered such that eigs(1) is the smallest eigvenvalue
c       eigs(2) the following, and so on.
c
        do i=nab,1,-1
          if (eigs(i).lt.cutoff(nd)) exit
          nloc = nloc + 1
          if (nbcen(nd).eq.2) then
            nbab=nbab+1
            if (nbab.gt.nmo*(nmo-1)/2) then
              stop '# otheroqs.f.f:Increase dimension 1 of iab'
            endif
            iab(1:2,nbab)=ibcen(nd,1:2)
            bond(nbab)=nloc
          endif
          nmoc(nbcen(nd))=nmoc(nbcen(nd))+1
          if (nloc.gt.nab) write (lw,433)
          if (nloc.gt.nab) stop
          eigloc(nloc)=eigs(i)
          forall (j=1:nab) ctmp(j,nloc)=v(j,i)
          do j=1,nab
            do k=1,nab
              prod=ctmp(j,nloc)*ctmp(k,nloc)*eigloc(nloc)
              sloc(j,k)=sloc(j,k)+prod
            enddo
          enddo
          if (nloc.le.nloctarget) then
            write (lw,116) nloc,eigs(i),
     &        (atnam(ibcen(nd,j))(3:4),ibcen(nd,j),j=1,nbcen(nd))
            if (nbcen(nd).eq.2) then
              i1=min(ibcen(nd,1),ibcen(nd,2))
              i2=max(ibcen(nd,1),ibcen(nd,2))
              enlaces(i1,i2)=enlaces(i1,i2)+1
            endif
          else
            write (locform(nloc-nloctarget),116) nloc,eigs(i),
     &        (atnam(ibcen(nd,j))(3:4),ibcen(nd,j),j=1,nbcen(nd))
          endif
        enddo
        if (deplet(nd)) sloc0=sloc
      enddo

      if (nloc.lt.nloctarget) then
        write (lw,302)
        warn = .true.
        return
      elseif (nloc.eq.nloctarget) then
        write (lw,300) 
        write (lw,301) (nmoc(i),i=1,mxcent)
      else
        write (lw,435)
        do i=nloctarget+1,nloc
          write (lw,'(a)') trim(locform(i-nloctarget))
        enddo
      endif
c
c-----Print type of bond between each pair of atoms
c
      write (lw,320)
 320  format (' # Type of bond between pairs of atoms',/,' # ',70('-'))
      do i=1,ncent-1
        do j=i+1,ncent
          if (enlaces(i,j).eq.0) then
            typebond = 'NO BOND                       '
          elseif (enlaces(i,j).eq.1) then
            typebond = 'SINGLE BOND                   '
          elseif (enlaces(i,j).eq.2) then
            typebond = 'DOUBLE BOND                   '
          elseif (enlaces(i,j).eq.3) then
            typebond = 'TRIPLE BOND                   '
          elseif (enlaces(i,j).eq.4) then
            typebond = 'QUADRUPLE BOND                '
          elseif (enlaces(i,j).eq.5) then
            typebond = 'QUINTUPLE BOND                '
          else
            typebond = 'SEXTUPLE or higher order BOND '
          endif
          write (lw,321) i,j,typebond
        enddo
      enddo
 321  format (' # Type of bond betweeen atoms ',2I4,' : ',a30)

      if (largwr) write (lw,20)
      do i=1,nloc
        forall (j=1:nab) vrot(i,j)=ctmp(j,i)
        if (largwr) write (lw,33) i
        if (largwr) write (lw,30) (vrot(i,j),j=1,nab)
      enddo
      cret=vrot(1:nloc,1:nab)
c
c-----Compute AOM between the new nonorthonormal localized MOs
c
      aom0 = aom
      do i=1,ncent
        aom(i,1:nloc,1:nloc) = matmul(vrot(1:nloc,:),
     &  matmul(aom0(i,:,:),transpose(vrot(1:nloc,:))))
      enddo
c
c-----Write diagonal AOM elements in each atom or fragment. 
c
      write (lw,68) 'non-orthonormal'
      write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
      do ip=1,nloc
        write (lw,410,advance='no') ip
        if (ncent.gt.10) then
          write (lw,41 ) (hundred*aom(k,ip,ip),k=1,10)
          write (lw,411) (hundred*aom(k,ip,ip),k=11,ncent)
          write (lw,*)
        else
          write (lw,41) (hundred*aom(k,ip,ip),k=1,ncent)
        endif
      enddo
c
c-----------------------------------------------------------------------
c
c-----effective number of atoms expanded by each localized MO
c
      do ip=1,nloc
        xloci(ip)=0d0
        do k = 1, ncent
          spp = aom(k,ip,ip)
          xloci(ip) = xloci(ip) + spp * spp
        enddo
      enddo
      write (lw,520) 'non-orthonormal',(1d0/xloci(ip),ip=1,nloc)
c
c-----Approximate two-center delocalization indices
c
      write (lw,*)
      write (lw,*) '# Approximate two-center delocalization indices'
      write (lw,*) '# ---------------------------------------------'
      do i=1,nbab
        i1=iab(1,i)
        i2=iab(2,i)
        ib=bond(i)
        di = abs(occ(ib)) * 2d0 * aom(i1,ib,ib) * aom(i2,ib,ib)
        write (lw,65) iab(1:2,i),di
      enddo
 65   format (' # Atoms ',2I3,'   Delta = ',F15.8)
c
c-----Orthogonalize the nonorthonormal Localized MOs by Lowdin. This is
c     equivalent to the maximum overlap method. First compute the overlap
c     in R^3 of the nonorthonormal Localized MOs.
c
c-----First, obtain S, where S is the overlap matrix in R^3 betwen
c     the just found nonorthonormal Localized MOs.
c
      sumaom=0d0
      do k=1,ncent
        sumaom(:,:)=sumaom(:,:)+aom(k,1:nloc,1:nloc)
      enddo
      if (largwr) then
        write (lw,*) '#'
        write (lw,*) '# Overlaps in R^3 of OQS MOs'
        do i=1,nloc
          write (lw,133) (sumaom(i,j),j=1,i)
        enddo
      endif
      n3=3*nloc-1
      uvec(1:nloc,1:nloc) = sumaom
c
c-----Diagonalize S
c
      call dsyev ('V','U',nloc,uvec(1:nloc,1:nloc),
     &           nloc,eigen(1:nloc),work,n3,info)
      call errdsyev (info,'otheroqs.f')
c
c-----Obtain  S^{-1/2} (stored in sumaom())
c
      sumaom = 0d0
      do i=1,nloc
        do j=1,nloc
          value=0d0
          do k=1,nloc
            value=value+uvec(i,k)*uvec(j,k)/sqrt(eigen(k))
          enddo
          sumaom(i,j) = value
        enddo
      enddo
      vroto(1:nab,1:nloc)=matmul(transpose(vrot(1:nloc,1:nab)),sumaom)
      vroto(1:nloc,1:nab)=transpose(vroto(1:nab,1:nloc))
      do i=1,ncent
        aomo(i,1:nloc,1:nloc)=matmul(transpose(sumaom),
     &  matmul(aom(i,1:nloc,1:nloc),sumaom))
      enddo
c
c-----Write orthonormalized MOs in terms of the canonical MOs
c
      if (largwr) then
        write (lw,2001)
        do i=1,nloc
          write (lw,33) i
          write (lw,30) (vroto(i,j),j=1,nab)
        enddo
      endif
      coret=vroto(1:nloc,1:nab)
c
c-----------------------------------------------------------------------
c
c-----Write diagonal AOM elements in each atom or fragment between
c     orthonormalized localized MOs
c
      write (lw,68) 'orthonormal'
      write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
      do ip=1,nloc
        write (lw,410,advance='no') ip
        if (ncent.gt.10) then
          write (lw,41)  (hundred*aomo(k,ip,ip),k=1,10)
          write (lw,411) (hundred*aomo(k,ip,ip),k=11,ncent)
          write (lw,*)
        else
          write (lw,41) (hundred*aomo(k,ip,ip),k=1,ncent)
        endif
      enddo
c
c-----------------------------------------------------------------------
c
c-----effective number of atoms expanded by each orthonormal MO
c
      do ip=1,nloc
        xloci(ip)=0d0
        do k = 1, ncent
          spp = aomo(k,ip,ip)
          xloci(ip) = xloci(ip) + spp * spp
        enddo
      enddo
      write (lw,520) 'orthonormalized ',(1d0/xloci(ip),ip=1,nloc)
c
c-----------------------------------------------------------------------
c
c-----Overlaps in R^3 of the orthonormalized OQS MOs
c
      sumaom=0d0
      do k=1,ncent
        sumaom=sumaom+aomo(k,1:nloc,1:nloc)
      enddo
      if (largwr) then
        write (lw,*) '#'
        write (lw,*) '# Overlaps in R^3 of the orthonormalized OQS MOs'
        do i=1,nloc
          write (lw,133) (sumaom(i,j),j=1,i)
        enddo
      endif
      write (lw,112) 'END'
c
      deallocate (nmoc)
      write (lw,10) '  E N D  ',''
      return
      call timer (4,iiloc,'_otheroqs',-1)
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
 300  format (' #',/,' # ',70('-'),/,
     &  ' # CONVERGED (1--N)-CENTER OPEN SYSTEM ANALYSIS',/,' #')
 301  format (' # Number of MOs involving 1, 2, 3, ... centers',/,
     &        ' # ',20(1x,I3))
 302  format (' #',/,' # otheroqs.f: Unable to localize all the MOs',
     & /,' #',/,' # ',100('-'),/,
     & ' # NON CONVERGED (1--N)-CENTER OPENLOC ANALYSIS',/,' #')
 66   format (/' # Diagonal overlaps of each MO in each atom')
 68   format (/' # Degree of localization of each MO in each atom')
 33   format (' # MO ',I3)
 10   format (/,' # ',90('-'),/,' # ',a,
     & '(1--N)-CENTER OPEN SYSTEM LOCALIZATION ANALYSIS',3x,a,
     & /,' # ',90('-'))
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
 100  format (' # otheroqs.f: File ',a,' NOT FOUND')
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
 433  format (' #',/,' # Too many Open Systems Localized MOs')
 520  format (' # Atoms expanded by each ',a,' OQS MO',/,
     &    1000(5(3x,F12.5),/))
 133  format (5(1x,F15.8))
 2001 format (/' # Rotation Matrix (READ BY ROWS)',6x,
     & 'Orthonormal OQS MO_i =  SUM_j CANMO_j C_ij')
      end
