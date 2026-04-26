c
c-----------------------------------------------------------------------
c
      subroutine uhfchemloc (sgab,maxcritov,critov,damps,critics,covx,
     &  hesel,lw,largwr,okw,wfnfile)
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
c.....sgab()                         INPUT
c
c     Group overlap matrix between input MOs.
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
c.....hesel                          INPUT
c
c     .true.    The weights defined by Andreas Heselmann are used.
c               (A. Heselmann, JCTC 12, 2720-2741 (2016).)
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
c
c.....largwr                         INPUT
c
c     .true.  ====> Large Output
c     .false. ====> Shrt  Output.
c
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
c
      parameter (maxcoord = 20)
      real(kind=8)    pi,pi4
      real(kind=8)    covx
      real(kind=8)    sgab(ncent,nmo,nmo)
      real(kind=8),   allocatable, dimension (:)     :: eigs,work,diags
      real(kind=8),   allocatable, dimension (:,:)   :: v
      real(kind=8),   allocatable, dimension (:,:,:) :: sg0
      real(kind=8),   allocatable, dimension (:,:,:) :: sg
      real(kind=8),   allocatable, dimension (:,:)   :: c
      real(kind=8),   allocatable, dimension (:,:)   :: cc
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

      integer(kind=4)  :: spin
      logical   warn,largwr,okw,condition,dow
      logical   there_are_MOs,not_hydrogen
      logical   hesel
      character*20 spinlbl
      character*(*) wfnfile

      call timer (2,ipid,'_uhfchemlo',-1)
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
        allocate (iords(ncent))
        allocate (work(3*n-1))
        allocate (v(n,n))
        allocate (sg0(ncent,n,n))
        allocate (wh(ncent,maxcoord))
        allocate (blg(ncent,n))
        allocate (nbatom(ncent))
        allocate (ncom(n))
        allocate (over(n,n))
        allocate (sgtmp(ncent,n,n))
        allocate (xloci(n))
        allocate (ixion(n))
        allocate (ctmp(n,n))
        allocate (madis(ncent,ncent))
        allocate (madaux(ncent))
        allocate (coord(ncent))

        pi4  =  atan(one)
        four =  two+two
        pi   =  four*pi4
        warn = .false.
        sg0  =  sg
  
        write (lw,10) trim(spinlbl)
        write (lw,'(a)') ' # '//trim(spinlbl)// ' ORBITALS ARE'
        if (spin.eq.1) write (lw,'(20(1x,I3))') ialpha(1:n)
        if (spin.ne.1) write (lw,'(20(1x,I3))') ibeta(1:n)
        write (lw,21) critov(1),critov(2),damps,critics,covx
c
c.......Print input group overlap integrals.
c
        if (okw.and.largwr) then
          write (lw,*) '# INPUT ATOMIC OVERLAP MATRIX BETWEEN MOs'
        endif
        over(1:n,1:n)=zero
        diags=zero
        do k=1,ncent
          over(1:n,1:n)=over(1:n,1:n)+sg(k,1:n,1:n)
          do ip=1,n
            diags(k)=diags(k)+sg(k,ip,ip)
          enddo
          if (okw.and.largwr) then
            write (lw,'(a,I4)') ' # ATOM ',k
            do ip=1,n
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
c-------Order the A centers by increasing value of SUM_i S_ii^A
c
        forall (i=1:ncent) iords(i)=i
        call qcksort (diags, iords, 1, ncent)
        deallocate (diags)
        write (lw,321) (iords(i),i=1,ncent)
c
c-------Obtain coordinations, connectivities, and distance arrays
c
        mxc = maxcoord
        call connect (lw,covx,ncent,mxc,largwr,nbonds,coord,wh,madis)
c
c-------Largest number of bonds between two 'connected' atoms
c
        longest_chain=maxval(madis) 
c
c.......Start process. Initial rotation matrix is the unit matrix.
c
        ctmp = zero
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
              call errdsyev (info,'uhfchemloc.f')
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
     &          matmul(sgtmp(k,1:nfin,1:nfin),v(1:nfin,1:nm)))
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
                  call errdsyev (info,'uhfchemloc.f')
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
     &              matmul(sgtmp(k,1:nfin,1:nfin),v(1:nfin,1:nm)))
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
          write (lw,*)'# uhfchemloc.f: Unable to localize all the MOs'
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

        do ip=1,nloc
          xloci(ip)=ddot(ncent,sg(:,ip,ip),1,sg(:,ip,ip),1)
        enddo
        write (lw,54) 
        write (lw,30) (xloci(i),i=1,nloc)
        write (lw,520)
        write (lw,30)  (one/xloci(i),i=1,nloc)
        if (largwr) then
          write (lw,20)
          do i=1,nloc
            write (lw,33) i
            write (lw,30) (c(i,j),j=1,n)
          enddo
        endif
c
c-------Test of the obtained AOM elements.
c
        call aomtest (sg(:,1:nloc,1:nloc),nloc,ncent,critics,lw,largwr)
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
          open (unit=2,file=trim(wfnfile)//'chem-alpha')
          call wrtwfn (1,2,0,57,cc,inmo,inpr)
          write (lw,351) trim(wfnfile)//'chem-alpha'
        else
          open (unit=2,file=trim(wfnfile)//'chem-beta')
          call wrtwfn (1,2,0,57,cc,inmo,inpr)
          write (lw,351) trim(wfnfile)//'chem-beta'
        endif
c
        deallocate (eigs)
        deallocate (iords)
        deallocate (work)
        deallocate (v)
        deallocate (c)
        deallocate (cc)
        deallocate (sg)
        deallocate (sg0)
        deallocate (wh)
        deallocate (blg)
        deallocate (nbatom)
        deallocate (ncom)
        deallocate (over)
        deallocate (sgtmp)
        deallocate (xloci)
        deallocate (ixion)
        deallocate (ctmp)
        deallocate (madis)
        deallocate (madaux)
        deallocate (coord)
      enddo
      if (okw) write (lw,'(a,/,a,/,a)') ' #',' # Done',' #'
      call timer (4,ipid,'_uhfchemlo',-1)
      return
c
c.....Formats
c
 20   format (' # Rotation Matrix (READ BY ROWS)',
     &        '    locMO_i =  SUM_j CANMO_j C_ji')
 200  format (' # ROTATION MATRIX (READ BY ROWS)',
     &        '    locMO_i =  SUM_j CANMO_j C_ji')
 30   format (5(1x,F15.8))
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
 10   format (/,' # ',80('-'),/,
     &  ' # C H E M I C A L   L O C A L I Z A T I O N',
     &   2x,'----> ',a,'   B L O C K',/,' # ',80('-'))
 334  format (/,' # Second & fourth moment orbital spread of ',a)
 351  format (1x,"#",/,1x,"# Writing file ","'",a,"' with",
     & " Localized MOs instead of Natural MOs")
      end
