c
c.......................................................................
c
      subroutine edfstd (lr,lw,lerr,lu18)
c
c.......................................................................
c
c     EEEEEEEE  DDDDDD      FFFFFFFFF   SSSSSS   TTTTTTTTT   DDDDDD     
c     EEEEEEEE  DDDDDDDD    FFFFFFFFF  SSSSSSSS  TTTTTTTTT   DDDDDDDD   
c     EE        DD   DDDD   FF         SS           TT       DD   DDDD  
c     EE        DD      DD  FF        SSS           TT       DD      DD
c     EE        DD      DD  FF         SS           TT       DD      DD
c     EEEEEEEE  DD      DD  FFFFFFFFF   SSS         TT       DD      DD
c     EEEEEEEE  DD      DD  FFFFFFFFF     SSSS      TT       DD      DD
c     EE        DD      DD  FF               SS     TT       DD      DD
c     EE        DD      DD  FF               SSS    TT       DD      DD
c     EE        DD   DDDD   FF               SS     TT       DD   DDDD  
c     EEEEEEEE  DDDDDDDD    FF         SSSSSSSS     TT       DDDDDDDD   
c     EEEEEEEE  DDDDDD      FF          SSSSSS      TT       DDDDDD     
c
c.......................................................................
c
c.....THIS ROUTINE DEALS EXCLUSIVELY WITH CLOSED-SHELL MOLECULES DESCRI-
c     BED BY SINGLE-DETERMINANT WAVE FUNCTIONS.
c
c     After reading the IOVERLAP value from the first record of the input 
c     file, and the second one with the name of the file containing the 
c     AOM, the main program EDF calls this routine in case IOVERLAP=-11. 
c     The rest of the input data is described below (CAPITAL lettters are 
c     fixed names).
c
c
c=====Record 3: NMO,NCENT,NGROUP
c
c     NMO    ( > 0 )  = Number of Molecular Orbitals (MO).
c     NCENT  ( > 0 )  = Total Number of Atoms.
c     NGROUP ( > 0 )  = Number of Groups of Atoms in which the molecule
c                       is divided
c
c=====Record 4.i (i=1,NGROUP-1)   ! Note that final value of i is NGROUP-1
c
c     nfugrp(i),(ifugrp(j,i),j=1,nfugrp(i)) [ minelec [  maxelec ] ]
c
c     Only groups from 1 to NGROUP-1 are defined below. The atoms in the 
c     last group are defined by difference. None atom can simultaneously
c     belong to two or more groups.
c
c     nfugrp(i)    = Number of atoms in group i.
c     ifugrp(j,i)  = Indices of atoms that belong to group i.
c     minelec      = Minimum number of ALPHA=BETA electrons in group i.
c     maxelec      = Maximum number of ALPHA=BETA electrons in group i.
c
c         Reading of MINELEC AND MAXELEC is optional. In case they are 
c         not read in, they will take the values 0 and NMO, respectively. 
c         This means that group i can have a minimum number of electrons
c         equal to 0, and a maximum value (for each spin) equal to NMO 
c         (the total number of electrons for each spin)
c   
c=====Record 5: END
c
c     This indicates the end of the input of the present calculation.
c
c=====Record 6: LARGE
c
c     Large output is requested. The default is short output.
c
c=====Record 7: DNO
c
c     After reading the AOM of all the atoms that belong to all the
c     groups except the last one, this order performs the following  
c     tasks:
c     
c     1) Computes, by difference, the Group Overlap Matrix of the 
c        last group: S_last.
c     2) Diagonalizes S_last.
c     3) Selects the eigenvalues of S_last smaller than a critical
c        overlap value (critover, see below), as they are partially
c        delocalized in all the fragments, i.e. they are not exclusi-
c        vely localized in the last fragment.
c     4) Construct the eigenvectors correponding to these eigenvalues. 
c        These eigenvectors will be the DNOs (linear combinations of the 
c        original MOs) that are partially localized in the fragment 
c        1 + 2 +... + NGROUP-1 
c     5) Reconstruct the AOM in terms of the DNOs.
c     6) Computes the EDF and the localization and delocalization 
c        indices between all the groups.
c
c=====Record 8: CRITOVERLAP critover
c
c     critover is the parameter that controls the above order. The 
c     default value is 0.99, and the minimum and maximum values that 
c     can be given in this order are 0.8 and 1.0, respectively. A value
c     very close to 1.0 means that step 3 in the above order will select
c     all the MOs. Contrariry, a small value will reduced the number
c     of selected MOs.
c
c=====Record 9: FULLOC
c
c     Localizes the Molecular Orbitals using the defined groups (not the 
c     individual atoms) to contruct the localization function which is 
c     maximized, and computes the EDF of the molecule. A series of simpli-
c     fications to avoid numerical instabilities, very similar to the 
c     ones used in the order DNO above are carried out before the EDF is 
c     actually obtained:
c
c     1) Performs the isopycnic localization.
c     2) Divides the isopycnic MOs in two sets: (a) MOs partially deloc-
c        alized in two or more fragments, and (b) MOs almost fully loca-
c        lized in an unique fragment. 
c     3) Excludes the set (b) from the computation of the EDF, adding
c        two electrons to the core of a fragment each time an isopycnic
c        MO is localized in this fragment.
c
c=====Record 10: EPSLOC epsloc
c
c     epsloc is the parameter that controls the FULLOC order. The 
c     default value is 1D-2, and the minimum and maximum values
c     1D-4 and 1D-2, respectively.
c
c=====Record 11: DAFH
c
c     Performs a DAFH analysis in each fragment and from the DAFH eigen-
c     values greater than epsdafh (see below) determines the maximum 
c     number of electrons in this fragment. 
c
c=====Record 12: EPSDAFH epsdafh
c
c     epsdafh is the parameter that controls the EPSDAFH order. If a
c     DAFH eigenvalue in a group satisfies value > epsdafh, the maximum
c     number of electrons of this group increased by 2.0.
c
c=====Record 13: PROBCUT pcut
c
c     In the output of probabilities, only those probabilities 
c     greater than pcut are written. The default is 1D-3
c
c     The records 6-12 can be given in any order.
c
c-----COMMENTS CONCERNING THE COMPLEX ATOMIC OVERLAP MATRIX (AOM)
c
c-----------------------------------------------------------------------
c
      USE          space_for_wfnbasis
      USE          space_for_wfncoef
      USE          space_for_cidet
      include     'implicit.inc'
      include     'wfn.inc'
      include     'constants.inc'
      include     'mline.inc'
      real   (kind=8), allocatable, dimension (:,:,:) :: aom
      real   (kind=8), allocatable, dimension (:,:,:) :: sg
      integer(kind=4), allocatable, dimension (:,:)   :: ifugrp,mimael
      integer(kind=4), allocatable, dimension (:)     :: nfugrp,mocore
      integer(kind=4), allocatable, dimension (:,:)   :: resnc
      real   (kind=8)  critover,epsloc,epsdafh
      integer(kind=4)  lw,lr
      logical     setint,ok,goon,dno,largwr,localize,dodafh
      logical     setword
      character*(mline) line,word,uppcase
c
c-----------------------------------------------------------------------
c
      call timer (2,iedfstd,'_edfstd   ',-1)
      write (lw,57) 
      read (lr,*) nmo,ncent,ngroup
      allocate (mocore(ngroup))
      allocate (nfugrp(ngroup))
      allocate (ifugrp(ncent,ngroup))
      allocate (mimael(2,ngroup))
      do i=1,ngroup-1
        read (lr,'(a)') line
        lp=1
        if (.not.setint(nfugrp(i),line,lp)) then
          write (lerr,*) 'edfstd.f: Error defining fragment',i
          stop
        endif
        do j=1,nfugrp(i)
          if (setint(ifugrp(j,i),line,lp)) then
            if (ifugrp(j,i).lt.0.or.ifugrp(j,i).gt.ncent) then
              stop '# edfstd.f: Not valid atom has been read in'
            endif
          else
            stop '# edfstd.f: Not valid definition of a group of atoms'
          endif
        enddo
        ok=setint(mimael(1,i),line,lp)
        if (ok) then
          if (mimael(1,i).lt.0) then
            mimael(1,i)=0
            ok=setint(mimael(2,i),line,lp)
            if (ok) then
              if (mimael(2,i).lt.0)   mimael(2,i)=0
              if (mimael(2,i).gt.nmo) mimael(2,i)=nmo
            else
              mimael(2,i)=nmo
            endif
          else
            mimael(1,i)=min(mimael(1,i),nmo)
            ok=setint(mimael(2,i),line,lp)
            if (ok) then
              if (mimael(2,i).lt.0)   mimael(2,i)=0
              if (mimael(2,i).gt.nmo) mimael(2,i)=nmo
            else
              mimael(2,i)=nmo
            endif
          endif
        else
          mimael(1,i)=0
          mimael(2,i)=nmo
        endif
        if (mimael(1,i).gt.mimael(2,i)) then
          write (lerr,*) 'edfstd.f: Error defining fragment',i
          stop
        endif
      enddo
c
      dno      = .false.
      localize = .false.
      dodafh   = .false.
      largwr   = .false.
      critover = 0.99d+00
      epsloc   = 1.0D-2
      epsdafh  = 1.0D-3
c
c.....Read control parameters used in the calculation.
c
      pcut=1.0D-3
      goon=.true.
      do while (goon) 
        read (lr,'(a)',end=110) line
        line=uppcase(line)
        lp=1
        ok=setword(word,line,lp)
        if (word(1:1).eq.'#') then
        elseif (trim(word).eq.'END') then
          goon=.false.
        elseif (trim(word).eq.'DNO') then
*         if (ngroup.eq.3) dno=.true.
          dno=.true.
        elseif (trim(word).eq.'LARGE') then
          largwr=.true.
        elseif (trim(word).eq.'FULLOC') then
          localize=.true.
        elseif (trim(word).eq.'DAFH') then
          dodafh=.true.
        elseif (trim(word).eq.'CRITOVERLAP') then
          line=line(lp:)
          read (line,*) critover
          critover = abs(critover)
          if (critover.gt.1.0d+00) critover = 1.0d+00
          if (critover.lt.0.8d+00) critover = 0.8d+00
        elseif (trim(word).eq.'EPSDAFH') then
          line=line(lp:)
          read (line,*) epsdafh
          epsdafh=max(min(abs(epsdafh),1.0D-2),1.0D-6)
        elseif (trim(word).eq.'EPSLOC') then
          line=line(lp:)
          read (line,*) epsloc
          epsloc=max(min(abs(epsloc),1.0D-2),1.0D-4)
        elseif (trim(word).eq.'PROBCUT') then
          line=line(lp:)
          read (line,*) pcut
        else
          write (lerr,'(a)') ' # '//trim(word)
          write (lerr,*) '# Key word in input file not understood'
        endif
      enddo
 110  continue
      if (dno) then
        localize=.false.
        dodafh=.false.
      elseif (localize) then
        dodafh=.false.
        dno=.false.
      elseif (dodafh) then
        localize=.false.
        dno=.false.
      endif
c
c.....Determine the atoms in the last group
c
      nfugrp(ngroup)=0
      cyclel: do l=1,ncent
        do i=1,ngroup-1
          do j=1,nfugrp(i)
            if (l.eq.ifugrp(j,i)) cycle cyclel
          enddo
        enddo
        nfugrp(ngroup)=nfugrp(ngroup)+1
        ifugrp(nfugrp(ngroup),ngroup)=l
      enddo cyclel
c
c.....Test that none group is empty of atoms.
c
      do i=1,ngroup
        if (nfugrp(i).eq.0) then
          write (lerr,*) 'edfstd.f:  GROUP ',i,' is empty of atoms.'
          stop
        endif
      enddo
c
c.....Test that each atom only belongs to an unique group.
c
      do i=2,ngroup
        do k=1,i-1
          do m=1,nfugrp(k)
            do l=1,nfugrp(i)
              if (ifugrp(l,i).eq.ifugrp(m,k)) then
                 write (lerr,*)
                 write (lerr,432) l,i,m,k
                 write (lerr,*)
                 stop
              endif
            enddo
          enddo
        enddo
      enddo
c
c-----Write groups of atoms
c
      write (lw,92) ngroup
      do i=1,ngroup
        write (lw,93) i,nfugrp(i),(ifugrp(j,i),j=1,nfugrp(i))
      enddo

      if (dno) then
        call dnostd (nmo,ncent,ngroup,nfugrp,ifugrp,pcut,
     &   critover,lw,lerr,lu18)
         call timer (4,iedfstd,'_edfstd   ',-1)
         return
      endif
      if (dodafh) then
         if (ngroup.eq.2) call dafhstd (nmo,ncent,ngroup,
     &     nfugrp,ifugrp,pcut,epsdafh,lr,lw,lerr,lu18)
         if (ngroup.gt.2) call gendafhstd (nmo,ncent,naom,
     &     ngroup,nfugrp,ifugrp,pcut,epsdafh,lr,lw,lerr,lu18)
         call timer (4,iedfstd,'_edfstd   ',-1)
         return
      endif
      write (lw,766) nmo,ncent,ngroup
c
c-----Compute the minimum and maximum population of the last group
c
      do k=1,ngroup-1
        mimael(1,k)=max(mimael(1,k),0)
        mimael(2,k)=min(mimael(2,k),nmo)
      enddo
      mimael(2,ngroup)=min(nmo-sum(mimael(1,1:ngroup-1)),nmo)
      mimael(1,ngroup)=max(nmo-sum(mimael(2,1:ngroup-1)),0)
      write (lw,*) '# Min/Max number of electrons in each group'
      do k=1,ngroup
        write (lw,171) k,mimael(:,k)
      enddo
 171  format (' # Group ',I3,3x,'MinElec,MaxElec = ',2I3)
c
      allocate (aom(ncent,nmo,nmo))
      allocate (sg(ngroup,nmo,nmo))
c
      npab=1
      do i=1,ngroup-1
        npab=npab*(mimael(2,i)-mimael(1,i)+1)
      enddo
      nprev=npab
c
      allocate (resnc(nprev,ngroup))
c
c.....Computation of resonance structures.
c
      call rnprobs (mimael,resnc,nmo,ngroup,npab,lw)
      if (nprev.lt.npab) then
        write (lerr,*) 'edfstd.f: NPREV < NPAB returned by rnprobs.f'
        stop
      endif
      write (lw,*) '# Number of RSRSs = ',npab
      do i=1,npab
        write (lw,'(1000(" # ",20I3))') resnc(i,1:ngroup)
      enddo
c
c.....Read AOM
c
      call readaomstd (aom,nfugrp,ifugrp,ncent,ngroup,nmo,naom,
     &   lu18,lerr)
c
c.....Compute Group Overlap integrals of all but the last one.
c
      sg(1:ngroup,1:nmo,1:nmo)=0.0d+00
      do i=1,ngroup-1
        do j=1,nfugrp(i)
          k=ifugrp(j,i)
          sg(i,:,:)=sg(i,:,:)+aom(k,:,:)
        enddo
      enddo
c
c-----Compute Group Overlap integrals for the last group
c
      do i=1,nmo
       do j=1,nmo
         if (i.eq.j) delta=1.0d+00
         if (i.ne.j) delta=0.0d+00
         sg(ngroup,i,j)=delta-sum(sg(1:ngroup-1,i,j))
       enddo
      enddo
c
c.....Compute EDF using an isopycnic localization
c
      if (localize) then
        call locastd (sg,nmo,ngroup,ncent,
     &      nfugrp,ifugrp,pcut,lw,epsloc,.true.)
         call timer (4,iedfstd,'_edfstd   ',-1)
         return
      endif

      mocore(1:ngroup)=0
      call compedfstd (pcut,npab,nprev,ngroup,nmo,lw,mimael,sg,
     &                 resnc,mocore)
      deallocate (nfugrp)
      deallocate (ifugrp)
      deallocate (aom)
      deallocate (sg)
      deallocate (mimael)
      deallocate (resnc)
      deallocate (mocore)
      call timer (4,iedfstd,'_edfstd   ',-1)
      return
 766  format (
     & ' # Number of Molecular Orbitals               = ',I4,/,
     & ' # Total number of atoms                      = ',I4,/,
     & ' # Number of groups that divide the molecule  = ',I4)
 57   format (1x,'# Atomic Overlap Matrix is read in')
 432  format (1x,'# Atom number ',I3,' of group ',I2,
     &        1x,'is equal to atom number ',I3,' of group ',I2)
 92   format (' # Molecule formed by ',I3,' framents')
 93   format(' # Fragmet ',I3,' has ',I3,' atoms :',1000(20I3,/,' # '))
      end
