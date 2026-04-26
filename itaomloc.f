c
c-----------------------------------------------------------------------
c
      subroutine itaomloc (sg,c,nmo,ncent,mxcent,maxcritov,critov,
     &    largwr,lw)
c
c.....Performs the localization of a set of NMO molecular orbitals into 
c     NCENT atoms by diagonalizing the AOM in single atoms and pairs of
c     atoms. After each diagonalization, the MOs which are localized in 
c     the atom or pair of atoms where the AOM was computed are excluded. 
c     The proccess goes on (progressively decreasing the CRITOV 
c     parameter) until there are NOT non localized MOs.
c                            
c                            
c-----PARAMETERS--------------------------------------------------------
c
c.....sg()                           INPUT/OUTPUT
c
c     Group overlap matrix between input MO's when entering the 
c     routine or between output MO's at the end of the routine.
c
c.....c()                            OUTPUT 
c
c     Rotation matrix giving the transformed MO's from the input MO's. 
c
c.....nmo                            INPUT
c
c     Number of MO's.
c
c.....ncent                          INPUT
c
c     Number of fragments 
c
c.....critov()                       INPUT
c
c     If the overlap of an MO with itself in a domain is greater than
c     CRITOV(1), the routine considers that this MO is fully localized 
c     on that domain. If the overlap of an MO with itself in a pair of 
c     domains is greater than CRITOV(2), the routine considers that 
c     this MO is fully localized on that pair of domains, and so on.
c
c.....largwr
c
c     If LARGWR=.TRUE. a slightly greater output is written
c
c.....lw                             INPUT
c
c     Output logical unit
c-----------------------------------------------------------------------
c                         
      USE      space_for_wfnbasis
      include 'implicit.inc'
      include 'constants.inc'
      include 'error.inc'
      real(kind=8)    c(nmo,nmo)
      real(kind=8)    sg(ncent,nmo,nmo)
      real(kind=8),   allocatable, dimension (:)     :: eigs,work,diags
      real(kind=8),   allocatable, dimension (:,:)   :: v
      real(kind=8),   allocatable, dimension (:,:,:) :: sg0
      real(kind=8),   allocatable, dimension (:)     :: xloci
      real(kind=8),   allocatable, dimension (:,:,:) :: sgtmp
      real(kind=8),   allocatable, dimension (:,:)   :: ctmp
      integer(kind=4),allocatable, dimension (:)     :: iord
      integer(kind=4),allocatable, dimension (:)     :: inside
      real(kind=8) critov(maxcritov)
      integer(kind=8) ntuples
      logical   largwr,dow,first
c
      call timer (2,ipid,'_itaomloc ',-1)
c
c.....Allocate arrays
c
      n = nmo
      if (.not.allocated(eigs)) then
        allocate (eigs(n),stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot allocate eigs()'
      endif
      if (.not.allocated(diags)) then
        allocate (diags(ncent),stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot allocate diags()'
      endif
      if (.not.allocated(iord)) then
        allocate (iord(ncent),stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot allocate iord()'
      endif
      if (.not.allocated(inside)) then
        allocate (inside(ncent),stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot allocate inside()'
      endif
      if (.not.allocated(work)) then
        allocate (work(3*n-1),stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot allocate work()'
      endif
      if (.not.allocated(v)) then
        allocate (v(n,n),stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot allocate v()'
      endif
      if (.not.allocated(sg0)) then
        allocate (sg0(ncent,n,n),stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot allocate sg0()'
      endif
      if (.not.allocated(sgtmp)) then
        allocate (sgtmp(ncent,n,n),stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot allocate sgtmp()'
      endif
      if (.not.allocated(xloci)) then
        allocate (xloci(n),stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot allocate xloci()'
      endif
      if (.not.allocated(ctmp)) then
        allocate (ctmp(n,n),stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot allocate ctmp()'
      endif
      sg0  =  sg
      write (lw,10) '         '
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
      do while (natoms.le.ncenteff.and.nloc.lt.nmo)
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
*             write (lw,116) k,eigs(k),(inside(j),j=1,natoms)
              write (lw,116) k,eigs(k),
     &          (atnam(inside(j))(3:4),inside(j),j=1,natoms)
            endif
          enddo
          n3=3*nm-1
          do k=1,ncent   ! Transform the AOM matrix
            sg(k,1:nm,1:nm)=matmul(transpose(v(1:nfin,1:nm)),
     &      matmul(sgtmp(k,1:nfin,1:nfin),v(1:nfin,1:nm)))
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
      if (allocated(eigs)) then
        deallocate (eigs,stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot deallocate eigs()'
      endif
      if (allocated(iord)) then
        deallocate (iord,stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot deallocate iord()'
      endif
      if (allocated(work)) then
        deallocate (work,stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot deallocate work()'
      endif
      if (allocated(v)) then
        deallocate (v,stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot deallocate v()'
      endif
      if (allocated(sg0)) then
        deallocate (sg0,stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot deallocate sg0()'
      endif
      if (allocated (sgtmp)) then
        deallocate (sgtmp,stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot deallocate sgtmp()'
      endif
      if (allocated (xloci)) then
        deallocate (xloci,stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot deallocate xloci()'
      endif
      if (allocated (ctmp)) then
        deallocate (ctmp,stat=ier)
        if (ier.ne.0) stop 'itaomloc.f: Cannot deallocate ctmp()'
      endif
c
c-----Write diagonal AOM elements in each atom or fragment
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
      write (lw,10) '  E N D  '
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
*116  format (3x,'MO ',I4,' ( eig = ',F12.5,
*    &         ' )     localized on atom(s) ',100I4)
 116  format (4x,'MO ',I4,' with eigenvalue = ',F12.5,
     & ') localized on atom(s)',10(1x,'(',A2,I2,')'))
 21   format (' # CRITOV (1 ..',I2,')  = ',10(1x,F5.3))
*10   format (/,' # ',70('-'),/,' # ',a,/,' # ',70('-'))
 10   format (/,' # ',70('-'),/,' # ',a,
     & '(1--N)-C E N T E R   A O M   A N A L Y S I S',/,
     & ' # ',70('-'))
      end
