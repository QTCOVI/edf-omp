c
c-----------------------------------------------------------------------
c
      subroutine chemloc (sg,c,nmo,ncent,maxcritov,critov,damps,
     & critics,covx,hesel,lw,largwr,okw)
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
c-----maxcritov                      INPUT
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
c-----damps
c     
c     In every cycle of localization on pairs of domains the value of
c     CRITOV(2) is multiplied by DAMPS (a number that should be smaller
c     than 1.0) with the aim of reducing the localization requirement.
c
c.....critics                        INPUT
c
c     This routine determines the set of localized MOs that are partially
c     localized in ecah atom. A MO fulfills this condition if 
c     S_ii^A > critics.
c
c.....covx                           INPUT
c
c     If the distance between two atoms A and B is smaller than
c     (RA+RB)*COVX, a bond is assumed to exist between A and B.
c     RA and RB are the covalent radius stored in the array covrad()
c     of the 'connect.f' routine.
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
c.....okw                            INPUT
c
c     .false. ====> The routine writes nothing.
c     .true.  ====> Output according to the value of the LARGWR value
c-----------------------------------------------------------------------
c                         
      USE      space_for_wfnbasis
      include 'implicit.inc'
      include 'constants.inc'
      include 'error.inc'
c
      parameter (maxcoord = 20)
      real(kind=8)    pi,pi4
      real(kind=8)    covx
      real(kind=8)    c(nmo,nmo)
      real(kind=8)    sg(ncent,nmo,nmo)
      real(kind=8),   allocatable, dimension (:)     :: eigs,work,diags
      real(kind=8),   allocatable, dimension (:,:)   :: v
      real(kind=8),   allocatable, dimension (:,:,:) :: sg0
      real(kind=8),   allocatable, dimension (:)     :: xloci
      real(kind=8),   allocatable, dimension (:,:)   :: over
      real(kind=8),   allocatable, dimension (:,:,:) :: sgtmp
      real(kind=8),   allocatable, dimension (:,:)   :: ctmp
      integer(kind=4),allocatable, dimension (:)     :: iords
      integer(kind=4),allocatable, dimension (:,:)   :: wh
      integer(kind=4),allocatable, dimension (:)     :: ixion
      integer(kind=4),allocatable, dimension (:)     :: coord
      integer(kind=4),allocatable, dimension (:,:)   :: madis
      integer(kind=4),allocatable, dimension (:)     :: madaux
      integer, allocatable,dimension (:,:)           :: blg
      integer, allocatable,dimension (:)             :: nbatom,ncom
      real(kind=8) critov(maxcritov)

      logical   warn,largwr,okw,condition,dow
      logical   there_are_MOs,not_hydrogen
      logical   hesel

      call timer (2,ipid,'_chemloc  ',-1)
c
c.....Allocate arrays
c
      n = nmo
      if (.not.allocated(eigs)) then
        allocate (eigs(n),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate eigs()'
      endif
      if (.not.allocated(diags)) then
        allocate (diags(ncent),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate diags()'
      endif
      if (.not.allocated(iords)) then
        allocate (iords(ncent),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate iords()'
      endif
      if (.not.allocated(work)) then
        allocate (work(3*n-1),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate work()'
      endif
      if (.not.allocated(v)) then
        allocate (v(n,n),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate v()'
      endif
      if (.not.allocated(sg0)) then
        allocate (sg0(ncent,n,n),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate sg0()'
      endif
      if (.not.allocated(wh)) then
        allocate (wh(ncent,maxcoord),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate wh()'
      endif
      if (.not.allocated(blg)) then
        allocate (blg(ncent,n),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate blg()'
      endif
      if (.not.allocated(nbatom)) then
        allocate (nbatom(ncent),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate nbatom()'
      endif
      if (.not.allocated(ncom)) then
        allocate (ncom(n),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate ncom()'
      endif
      if (.not.allocated(over)) then
        allocate (over(n,n),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate over()'
      endif
      if (.not.allocated(sgtmp)) then
        allocate (sgtmp(ncent,n,n),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate sgtmp()'
      endif
      if (.not.allocated(xloci)) then
        allocate (xloci(n),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate xloci()'
      endif
      if (.not.allocated(ixion)) then
        allocate (ixion(n),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate ixion()'
      endif
      if (.not.allocated(ctmp)) then
        allocate (ctmp(n,n),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate ctmp()'
      endif
      if (.not.allocated(madis)) then
        allocate (madis(ncent,ncent),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate madis()'
      endif
      if (.not.allocated(madaux)) then
        allocate (madaux(ncent),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate madaux()'
      endif
      if (.not.allocated(coord)) then
        allocate (coord(ncent),stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate coord()'
      endif

      pi4  =  atan(one)
      four =  two+two
      pi   =  four*pi4
      warn = .false.
      sg0  =  sg

      write (lw,10)
      write (lw,21) critov(1),critov(2),damps,critics,covx
c
c.....Print input group overlap integrals.
c
      if (okw.and.largwr) then
        write (lw,*) '# INPUT ATOMIC OVERLAP MATRIX BETWEEN MOs'
      endif
      over(1:n,1:n)=zero
      diags=zero
      do k=1,ncent
        over(1:n,1:n)=over(1:n,1:n)+sg(k,1:n,1:n)
        if (okw.and.largwr) then
          write (lw,'(a,I4)') ' # ATOM ',k
          do ip=1,n
            diags(k)=diags(k)+sg(k,ip,ip)
            write (lw,'(5(1x,F22.15))')(sg(k,ip,iq),iq=1,ip)
          enddo
        endif
      enddo
      if (okw.and.largwr) then
        write (lw,*) '# SUM OVER ATOMS'
        do ip=1,n
          write (lw,'(5(1x,F22.15))')(over(ip,iq),iq=1,ip)
        enddo
      endif
c
c-----Order the A centers by increasing value of SUM_i S_ii^A
c
      forall (i=1:ncent) iords(i)=i
      call qcksort (diags, iords, 1, ncent)
      if (allocated(diags)) then
        deallocate (diags,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate diags()'
      endif
      write (lw,321) (iords(i),i=1,ncent)
c
c-----Obtain coordinations, connectivities, and distance arrays
c
      call connect (lw,covx,ncent,maxcoord,largwr,nbonds,coord,wh,madis)
c
c-----Largest number of bonds between two 'connected' atoms
c
      longest_chain=maxval(madis) 
c
c.....Start process. Initial rotation matrix is the unit matrix.
c
      ctmp=zero
      forall (i=1:n) ctmp(i,i)=one
      there_are_MOs = .true.
      ichain = 0
      nm     = n
      n3     = 3*nm-1
      nloc   = 0
      sgtmp  = sg
      iter   = 1
      maxit  = 1
      write (lw,131) damps

      if (.not.hesel) then
        if (critov(1).gt.0.98D0) maxit=2
      else
        maxit = 1000
      endif
      do while (there_are_MOs.and.iter.le.maxit)
        icent=1
        dow=.true.
        do while (icent.le.ncent.and.there_are_MOs)
          i=iords(icent)
          not_hydrogen = .not.nint(charge(i)).eq.1
          if (not_hydrogen) then
            v(1:nm,1:nm) = sgtmp(i,1:nm,1:nm)
            eigs=0d0
            call dsyev ('V','U',nm,v,n,eigs(1:nm),work,n3,info)
            call errdsyev (info,'chemloc.f')
            c(:,1:nm)=matmul(ctmp(:,1:nm),v(1:nm,1:nm))  ! Update C
            ctmp(:,1:nm)=c(:,1:nm)
            nfin=nm
            do k=nfin,1,-1
              if (eigs(k).lt.critov(1)) then
                if (largwr) then
                  if (dow) then
                    write (lw,113) ' One ',critov(1),iter
                    dow=.false.
                  endif
                  write (lw,117) k,eigs(k)
                else
                  exit
                endif
              else
                nm = nm - 1
                nloc = nloc+1
                if (dow) then
                  write (lw,113) ' One ',critov(1),iter
                  dow=.false.
                endif
                write (lw,116) k,eigs(k),i
              endif
            enddo
            if (nloc.eq.n) there_are_MOs = .false.
            n3=3*nm-1
            do k=1,ncent   ! Transform the AOM matrix
              sg(k,1:nm,1:nm)=matmul(transpose(v(1:nfin,1:nm)),
     &        matmul(sgtmp(k,1:nfin,1:nfin),v(1:nfin,1:nm)))
            enddo
            sgtmp(1:ncent,1:nm,1:nm)=sg(1:ncent,1:nm,1:nm)
          endif
          icent=icent+1
        enddo
        critov(1)=critov(1)*damps
        iter=iter+1
      enddo

      iter=1
      do while (there_are_MOs)
        i=1
        dow=.true.
        do while (i.le.ncent-1.and.there_are_MOs)
          if (coord(i).gt.0) then
            j=i+1
            do while (j.le.ncent.and.there_are_MOs)
              if (madis(i,j).eq.1) then
                v(1:nm,1:nm)=sgtmp(i,1:nm,1:nm)+sgtmp(j,1:nm,1:nm)
                eigs=0d0
                call dsyev ('V','U',nm,v,n,eigs(1:nm),work,n3,info)
                call errdsyev (info,'chemloc.f')
                c(:,1:nm)=matmul(ctmp(:,1:nm),v(1:nm,1:nm))  ! Update C
                ctmp(:,1:nm)=c(:,1:nm)
                nfin=nm
                do k=nfin,1,-1
                  if (eigs(k).lt.critov(2)) then
                    if (largwr) then
                      if (dow) then
                        write (lw,113) ' Two ',critov(2),iter
                        dow=.false.
                      endif
                      write (lw,117) k,eigs(k)
                    else
                      exit
                    endif
                  else
                    nm = nm - 1
                    nloc = nloc+1
                    if (dow) then
                      write (lw,113) ' Two ',critov(2),iter
                      dow=.false.
                    endif
                    write (lw,116) k,eigs(k),i,j
                  endif
                enddo
                if (nloc.eq.n) there_are_MOs = .false.
                n3=3*nm-1
                do k=1,ncent   ! Transform the AOM matrix
                  sg(k,1:nm,1:nm)=matmul(transpose(v(1:nfin,1:nm)),
     &            matmul(sgtmp(k,1:nfin,1:nfin),v(1:nfin,1:nm)))
                enddo
                sgtmp(1:ncent,1:nm,1:nm)=sg(1:ncent,1:nm,1:nm)
              endif
              j=j+1
            enddo
          endif
          i=i+1
        enddo
        critov(2)=critov(2)*damps
        iter=iter+1
      enddo

      if (nloc.lt.n) then
        write (lw,*)'# chemloc.f: Unable to localize all the MOs'
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
c
      if (allocated(eigs)) then
        deallocate (eigs,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate eigs()'
      endif
      if (allocated(iords)) then
        deallocate (iords,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate iords()'
      endif
      if (allocated(work)) then
        deallocate (work,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate work()'
      endif
      if (allocated(v)) then
        deallocate (v,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate v()'
      endif
      if (allocated(sg0)) then
        deallocate (sg0,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate sg0()'
      endif
      if (allocated(wh)) then
        deallocate (wh,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot allocate wh()'
      endif
      if (allocated(blg)) then
        deallocate (blg,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate blg()'
      endif
      if (allocated(nbatom)) then
        deallocate (nbatom,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate nbatom()'
      endif
      if (allocated(ncom)) then
        deallocate (ncom,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate ncom()'
      endif
      if (allocated (over)) then
        deallocate (over,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate over()'
      endif
      if (allocated (sgtmp)) then
        deallocate (sgtmp,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate sgtmp()'
      endif
      if (allocated (xloci)) then
        deallocate (xloci,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate xloci()'
      endif
      if (allocated (ixion)) then
        deallocate (ixion,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate ixion()'
      endif
      if (allocated (ctmp)) then
        deallocate (ctmp,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate ctmp()'
      endif
      if (allocated(madis)) then
        deallocate (madis,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate madis()'
      endif
      if (allocated(madaux)) then
        deallocate (madaux,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate madaux()'
      endif
      if (allocated(coord)) then
        deallocate (coord,stat=ier)
        if (ier.ne.0) stop 'chemloc.f: Cannot deallocate coord()'
      endif
c
      if (okw) then
        write (lw,*) '#'
        write (lw,*) '# Done'
        write (lw,*) '#'
      endif
      call timer (4,ipid,'_chemloc  ',-1)
      return
c
c.....Formats
c
 20   format (' # Rotation Matrix (READ BY ROWS)',
     &        '    locMO_i =  SUM_j CANMO_j C_ji')
 200  format (' # ROTATION MATRIX (READ BY ROWS)',
     &        '    locMO_i =  SUM_j CANMO_j C_ji')
 30   format (10(1x,F12.6))
 301  format (5x,'Trace =  ',F12.6)
 31   format (2(/,2(5x,E14.7)))
 32   format (1x,'# Transformation matrix')
 300  format (' #',/,' # ',70('-'),/,' # CONVERGED LOCALIZATION',/,' #')
 302  format (' #',/,' # ',70('-'),/,
     &        ' # NON CONVERGED LOCALIZATION',/,' #')
 400  format (1x,'# CONVERGED LOCALIZATION')
 54   format (' # L_i (i=1..N)')
 520  format (' # Fragments expanded by MOs, L_i^(-1) (i=1..N)')
 66   format (/' # Diagonal overlaps of each LMO in each atom')
 67   format (' # ATOM ',I2)
 33   format (' # LMO ',I3)
 1022 format (' # Center ',I4,/,1000(3x,'MOs',4x,20I4,/))
 1023 format (1x,'# Total MOs = ',I4)
 1024 format (/1x,'# MOs partially localized on each atom',
     &        '  (S_ii^A > ',1PE15.7,' )')
 1026 format (//' # MOs localized more than',F7.2,'% on each fragment')
 1027 format (' # ATOM ',I4,' --> ',20(1x,20I4,/))
 1025 format (/1x,'# Diagonal overlaps in each atom')
 1028 format (/' # MOs localized more than',F7.2,
     & '% on a pair of fragments')
 1031 format (/' # MOs localized more than',F7.2,
     & '% on a trio of fragments')
 1037 format (/' # MOs localized more than',F7.2,
     & '% on a quartet of fragments')
 1029 format (' # Fragments ',10I4)
 1030 format (3x,'MO ',I4,5x,F5.2,'% + ',F5.2,'%')
 1033 format (3x,'MO ',I4,5x,F5.2,'% + ',F5.2,'%+ ',F5.2,'%')
 1036 format (3x,
     &  'MO ',I4,5x,F5.2,'% + ',F5.2,'%+ ',F5.2,'%+ ',F5.2,'%')
 121  format (3x,'ATOM ',I4,a,' ATOM ',I4,' --> MOs ',25I4,/,
     &             1000(38x,25I4,/))
 122  format (' # MOs localized in two atoms simultaneously')
 124  format (' # MOs localized in one atom or the other')
 111  format (' # MO ( eig = ',F12.5,' ) localized on atom  ',100I3)
 112  format (' # MO ( eig = ',F12.5,' ) localized on atoms ',100I3)
 113  format (/,1x,70('-'),/,
     & 2x,a,'center localizations  Critical overlap = ',F12.5,
     & 2x,'ITER = ',I2,/,1x,70('-'))
 114  format (/,1x,70('-'),/,3x,
     &        'Functional groups localizations up to 1-',I1,/,
     &           1x,70('-'))
 116  format (3x,'MO ',I4,' ( eig = ',F12.5,
     &         ' )     localized on atom(s) ',100I4)
 117  format (3x,'MO ',I4,' ( eig = ',F12.5,' ) non localized yet')
 321  format (' #'/,' # Atoms are processed in the order',/,
     &  1000(' # ',20I4,/),/,' #')
 131  format (1x,'Damping factor of Critical Overlap = ',F12.5)
 21   format (/,
     &   ' # CRITOV(1)   = ',F16.10,/,' # CRITOV(2)   = ',F16.10,/,
     &   ' # DAMPS       = ',F16.10,/,' # CRITICS     = ',F16.10,/,
     &   ' # COVX        = ',F16.10)
 10   format (/,' # -----------------------------------------',/,
     &          ' # C H E M I C A L   L O C A L I Z A T I O N',/,
     &          ' # -----------------------------------------')
      end
