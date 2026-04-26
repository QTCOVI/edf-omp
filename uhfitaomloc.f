c
c-----------------------------------------------------------------------
c
      subroutine uhfitaomloc (sgab,mxcent,mxcrit,critov,largwr,lw,wfnf)
c
c.....Performs the localization of a set of NMO molecular orbitals into 
c     NCENT atoms by diagonalizing the AOM in single atoms and pairs of
c     atoms. After each diagonalization, the MOs which are localized in 
c     the atom or pair of atoms where the AOM was computed are excluded. 
c     The proccess goes on (progressively decreasing the CRITOV 
c     parameter) until there are NOT non localized MOs.
c-----------------------------------------------------------------------
c                            
c-----PARAMETERS--------------------------------------------------------
c
c.....sgab()                         INPUT
c
c     Group overlap matrix between input MOs.
c
c.....mxcent                         INPUT
c
c     Maximum number of N in the search of (Nc,2e) bonds.
c
c-----mxcrit                         INPUT
c
c     Maximum order 'n' in the search of (nc,2e) localized MOs
c
c.....critov()                       INPUT
c
c     If the overlap of an MO with itself in a domain is greater than
c     CRITOV(1), the routine considers that this MO is fully localized 
c     on that domain. If the overlap of an MO with itself in a pair of 
c     domains is greater than CRITOV(2), the routine considers that 
c     this MO is fully localized on that pair of domains.
c               
c.....lw                             INPUT
c
c     Output logical unit
c
c.....largwr                         INPUT
c
c     .true.  ====>  Large output
c     .false. ====>  Short output
c
c.....wfnf                           INPUT
c
c
c     Name of the WFN file
c
c-----------------------------------------------------------------------
c                         
      USE      space_for_wfnbasis
      USE      space_for_wfncoef
      include 'implicit.inc'
      include 'wfn.inc'
      include 'constants.inc'
      include 'error.inc'
      include 'corr.inc'
      real(kind=8)    sgab(ncent,nmo,nmo)
      real(kind=8),   allocatable, dimension (:)     :: eigs,work,diags
      real(kind=8),   allocatable, dimension (:,:)   :: v
      real(kind=8),   allocatable, dimension (:,:)   :: c
      real(kind=8),   allocatable, dimension (:,:)   :: cc
      real(kind=8),   allocatable, dimension (:,:,:) :: sg
      real(kind=8),   allocatable, dimension (:,:,:) :: sg0
      real(kind=8),   allocatable, dimension (:)     :: xloci
      real(kind=8),   allocatable, dimension (:,:,:) :: sgtmp
      real(kind=8),   allocatable, dimension (:,:)   :: ctmp
      integer(kind=4),allocatable, dimension (:)     :: iord
      integer(kind=4),allocatable, dimension (:)     :: inside
      real(kind=8) critov(mxcrit)
      integer(kind=8) ntuples
      integer(kind=4)  :: spin
      character*20 spinlbl
      character*(*) wfnf
      logical   largwr,dow,first
c
      call timer (2,ipid,'_itaomloc ',-1)
c
c.....Allocate arrays
c
      do spin=1,2
        if (spin.eq.1) spinlbl='ALPHA'
        if (spin.ne.1) spinlbl='BETA '
        if (spin.eq.1) n=nalpha
        if (spin.ne.1) n=nbeta
        allocate (sg(ncent,n,n))
        allocate (c(n,n))

        if (spin.eq.1) sg(:,:,:)=sgab(:,ialpha(1:n),ialpha(1:n))
        if (spin.ne.1) sg(:,:,:)=sgab(:,ibeta(1:n) ,ibeta(1:n) )

        allocate (eigs(n))
        allocate (diags(ncent))
        allocate (iord(ncent))
        allocate (inside(ncent))
        allocate (work(3*n-1))
        allocate (v(n,n))
        allocate (sg0(ncent,n,n))
        allocate (sgtmp(ncent,n,n))
        allocate (xloci(n))
        allocate (ctmp(n,n))
        sg0  =  sg
        write (lw,10) '         ',trim(spinlbl)
        ncenteff = min(mxcent,ncent)
        write (lw,21) ncenteff,(critov(j),j=1,ncenteff)
c
c.....Start the process.
c
        ctmp = zero
        forall (i=1:n) ctmp(i,i) = one
        nm     = n
        n3     = 3*nm-1
        nloc   = 0
        sgtmp  = sg
        forall (i=1:ncent) iord(i)=i
        natoms = 1
        do while (natoms.le.ncenteff.and.nloc.lt.n)
          write (lw,*) ' # Localization on ',natoms,' centers'
          xtuples=1d0
          do i=0,natoms-1
            xtuples=xtuples*dble(ncent-i)/dble(natoms-i)
          enddo
          ntuples=int(xtuples)
          first=.true.
          do i=1,ntuples
            call ngen (natoms,ncent,ncent,inside,first)
            v(1:nm,1:nm) = 0d0
            do j=1,natoms
              ji = inside(j)
              v(1:nm,1:nm)=v(1:nm,1:nm)+sgtmp(ji,1:nm,1:nm)
            enddo
            eigs = 0d0
            call dsyev ('V','U',nm,v,n,eigs(1:nm),work,n3,info)
            call errdsyev (info,'itaomloc.f')
            c(:,1:nm)=matmul(ctmp(:,1:nm),v(1:nm,1:nm))  ! Update C
            ctmp(:,1:nm)=c(:,1:nm)
            nfin=nm
            do k=nfin,1,-1
              if (eigs(k).lt.critov(natoms)) then
                exit
              else
                nm = nm - 1
                nloc = nloc+1
                if (dow) then
                  write (lw,113) natoms,critov(1)
                  dow = .false.
                endif
*               write (lw,116) k,eigs(k),(inside(j),j=1,natoms)
                write (lw,116) k,eigs(k),
     &            (atnam(inside(j))(3:4),inside(j),j=1,natoms)
              endif
            enddo
            n3=3*nm-1
            do k=1,ncent   ! Transform the AOM matrix
              sg(k,1:nm,1:nm)=matmul(transpose(v(1:nfin,1:nm)),
     &        matmul(sgtmp(k,1:nfin,1:nfin),v(1:nfin,1:nm)))
            enddo
            sgtmp(1:ncent,1:nm,1:nm)=sg(1:ncent,1:nm,1:nm)
          enddo
          natoms = natoms + 1 
        enddo
  
        if (nloc.lt.n) then
          write (lw,*)'# itaomloc.f: Unable to localize all the MOs'
          write (lw,302)
        else
          write (lw,300)
        endif
        c=transpose(ctmp)
        do k=1,ncent
          sg(k,:,:)=matmul(c,matmul(sg0(k,:,:),ctmp))
        enddo
c
c-------Write AOM between the transformed MOs
c
        if (largwr) call writeaom (sg,ncent,n,lw,5,.true.)

        do ip=1,n
          xloci(ip)=ddot(ncent,sg(:,ip,ip),1,sg(:,ip,ip),1)
        enddo
        write (lw,54) 
        write (lw,30) (xloci(i),i=1,n)
        write (lw,520)
        write (lw,30)  (one/xloci(i),i=1,n)
        if (largwr) then
          write (lw,20)
          do i=1,n
            write (lw,33) i
            write (lw,30) (c(i,j),j=1,n)
          enddo
        endif
cAQUI
c
c-------Test of the obtained AOM elements.
c
        call aomtest (sg(:,1:nloc,1:nloc),nloc,ncent,0.01d0,lw,largwr)
        allocate (cc(nloc,nprims))
        do i=1,nloc
          do ip=1,nprims
            cc(i,ip) = dot_product(coef(nmo+1:nmo+nmo,ip),c(i,:))
          enddo
        enddo
c
c-------Multipole moments of localized MOs are obtained.
c
        write (lw,334) 'Localized MOs'
        call newdrvmult (nloc,nprims,ncent,cc,lw)
        inmo=nloc
        inpr=nprims
        if (spin.eq.1) then
          open (unit=2,file=trim(wfnf)//'aom-alpha')
          call wrtwfn (1,2,0,57,cc,inmo,inpr)
          write (lw,351) trim(wfnf)//'aom-alpha'
        else
          open (unit=2,file=trim(wfnf)//'chem-beta')
          call wrtwfn (1,2,0,57,cc,inmo,inpr)
          write (lw,351) trim(wfnf)//'chem-beta'
        endif

        deallocate (eigs)
        deallocate (iord)
        deallocate (work)
        deallocate (v)
        deallocate (sg0)
        deallocate (sgtmp)
        deallocate (xloci)
        deallocate (ctmp)
c
c-------Write diagonal AOM elements in each atom or fragment
c
        write (lw,66)
        do k=1,ncent
          write (lw,*) '# Center ',k
          write (lw,30) (sg(k,ip,ip),ip=1,n)
        enddo
        write (lw,68)
        write (lw,*) '#    MO\ATOM --->   1      2      3      4  ...'
        do ip=1,n
          write (lw,410,advance='no') ip
          write (lw,41) (100d0*sg(k,ip,ip),k=1,ncent)
          if (ncent.gt.15) write (lw,*) 
        enddo
        write (lw,10) '  E N D  ',trim(spinlbl)
      enddo
      write (lw,'(a,/,a,/,a)') ' #',' # Done',' #'
      call timer (4,ipid,'_itaomloc ',-1)
      return
c
c.....Formats
c
 410  format (' # ',I4)
 20   format (' # Rotation Matrix (READ BY ROWS)',
     &        '    locMO_i =  SUM_j CANMO_j C_ji')
 200  format (' # ROTATION MATRIX (READ BY ROWS)',
     &        '    locMO_i =  SUM_j CANMO_j C_ji')
 30   format (6(1x,F12.6))
 41   format (1000(10x,15(1x,F5.1,'%')))
 301  format (5x,'Trace =  ',F12.6)
 31   format (2(/,2(5x,E14.7)))
 32   format (1x,'# Transformation matrix')
 300  format (' #',/,' # ',70('-'),/,' # CONVERGED LOCALIZATION',/,' #')
 302  format (' #',/,' # ',70('-'),/,
     &        ' # NON CONVERGED LOCALIZATION',/,' #')
 400  format (1x,'# CONVERGED LOCALIZATION')
 54   format (' # L_i (i=1..N)')
 520  format (' # Fragments expanded by MOs, L_i^(-1) (i=1..N)')
 66   format (/' # Diagonal overlaps of each MO in each atom')
 68   format (/' # Localization percentage of each MO in each atom')
 67   format (' # ATOM ',I2)
 33   format (' # MO ',I3)
 1026 format (//' # MOs localized more than',F7.2,'% on each fragment')
 1027 format (' # ATOM ',I4,' --> ',20(1x,20I4,/))
 1025 format (/1x,'# Diagonal overlaps in each atom')
 1029 format (' # Fragments ',10I4)
 1030 format (3x,'MO ',I4,5x,F5.2,'% + ',F5.2,'%')
 1033 format (3x,'MO ',I4,5x,F5.2,'% + ',F5.2,'%+ ',F5.2,'%')
 113  format (/,1x,70('-'),/,2x,I2,
     & '-center localizations, Critical overlap = ',F12.5,/,1x,70('-'))
 116  format (4x,'MO ',I4,' with eigenvalue = ',F12.5,
     & ') localized on atom(s)',10(1x,'(',A2,I2,')'))
 21   format (' # CRITOV (1 ..',I2,')  = ',10(1x,F5.3))
 10   format (/,' # ',80('-'),/,' # ',a,
     & '(1--N)-C E N T E R   A O M   A N A L Y S I S',2x,'----> ',a,
     & '   B L O C K'/,' # ',80('-'))
 334  format (/,' # Second & fourth moment orbital spread of ',a)
 351  format (1x,"#",/,1x,"# Writing file ","'",a,"' with",
     & " Localized MOs instead of Natural MOs")
      end
