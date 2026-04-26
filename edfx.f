c
c.......................................................................
c
      subroutine edfx (lr,lw,lu18)
c
c.......................................................................
c
c    EEEEEEEEEEEEEEE  DDDDDDDDDDD        FFFFFFFFFFFFFFF  XXX       XXX
c    EEEEEEEEEEEEEEE  DDDDDDDDDDDDD      FFFFFFFFFFFFFFF   XXX     XXX
c    EEE              DDD       DDDD     FFF                XX     XX
c    EEE              DDD          DDD   FFF                 XX   XX
c    EEE              DDD          DDD   FFF                  XX XX
c    EEEEEEEEEEEEEE   DDD          DDD   FFFFFFFFFFFFFF        XXX
c    EEEEEEEEEEEEEE   DDD          DDD   FFFFFFFFFFFFFF         X 
c    EEE              DDD          DDD   FFF                   XXX
c    EEE              DDD          DDD   FFF                  XX XX
c    EEE              DDD       DDDD     FFF                 XX   XX
c    EEEEEEEEEEEEEEE  DDDDDDDDDDDDD      FFF               XXX     XXX
c    EEEEEEEEEEEEEEE  DDDDDDDDDDD        FFF              XXX       XXX
c
c.......................................................................
c
c.....This routine deals with COMPLEX Atomic Overlap Matrices (AOMs).
c     Only closed shell RHF molecules are allowed.
c
c     After reading the first record with the value of IOVERLAP from the
c     (standard unit) input file, and the second one with the name of 
c     the file containing the AOM, this routine  is called in case that 
c     IOVERLAP=-1. The rest of the input data is described below 
c     (CAPITAL lettters are fixed names).
c
c
c=====Record 3: NMO,NCENT,NAOM,NGROUP
c
c     NMO    ( > 0 )  = Number of Complex Molecular Orbitals (MO).
c     NCENT  ( > 0 )  = Total number of Atoms.
c     NAOM   ( > 0 )  = Number of Atoms for which the AOM will be read in.
c     NGROUP ( > 0 )  = Number of groups that partition the molecule.
c
c=====Record 4.i (i=1,NGROUP-1)  
c
c     nfugrp(i),(ifugrp(j,i),j=1,nfugrp(i)) [ minelec [  maxelec ] ]
c
c     Only groups from 1 to NGROUP-1 are defined below. The atoms in the 
c     last group are defined by difference. Atoms in the first NGROUP-1 
c     groups must necessarily belong to the set 1 <= atom <= NAOM; i.e.
c     ifugrp(j,i) .le. NAOM necessarily.
c
c     nfugrp(i)    = Number of atoms in group i.
c     ifugrp(j,i)  = Indices of atoms that belong to group i.
c     minelec      = Minimum number of ALPHA=BETA electrons in group i.
c     maxelec      = Maximum number of ALPHA=BETA electrons in group i.
c   
c     MINELEC and MAXELEC values are optional. If not given they are 
c     set to 0 and NMO, respectively. 
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
c     This order applies only when NGROUP=3. In this case, after reading
c     the AOM, this order performs the following tasks:
c     
c                                    1
c     1) Computes the matrix D_12 =  - [ S_1 S_2 + S_2 S_1 ]
c                                    2
c
c     where 1 and 2 labels the groups 1 and 2, respectively.
c
c     2) Diagonalizes D_12
c     3) Selects the eigenvalues of D_12 greater than epsbond (see below)
c     4) Construct the eigenvectors correponding to the above eigenvalues. 
c        These eigenvectors will be the DNOs (linear combinations of the 
c        original MOs), which are fully or partially localized in the joint 
c        fragment 1+2.
c     5) Reconstruct the AOM in terms of the DNOs.
c     6) Computes the EDF and the localization and delocalization 
c        indices
c
c
c=====Record 8: EPSBOND epsbond
c
c     epsbond is the parameter that controls the above order. The 
c     default value is 1D-7, and the minimum and maximum values that 
c     can be given in this order are 1D-10 and 1D-3, respectively.
c
c=====Record 9: FULLOC
c
c     Localizes the Molecular Orbitals using the defined groups (not the 
c     individual atoms) to contruct the localization function which is 
c     maximized, and computes the EDF of the molecule. A series of simpli-
c     fications to avoid numerical instabilities, very similar to the 
c     ones used in the order DNO above, are carried out before the EDF is 
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
      include     'stderr.inc'
      include     'constants.inc'
      include     'mline.inc'
      complex*16, allocatable,dimension (:,:,:)  :: aomcmplx
      complex*16, allocatable,dimension (:,:,:)  :: sgcmplx
      integer,    allocatable,dimension (:,:)    :: ifugrp,eleca
      integer,    allocatable,dimension (:)      :: nfugrp,mocore
      integer,    allocatable,dimension (:,:)    :: resnca,resncb
      complex*16  aom1cmplx
      real(kind=8)      aom1
      integer     lw,lr
      logical     setint,ok,notingroup,goon,dno,largwr,localize,dodafh
      logical     setword
      character*(mline) line,word,uppcase
      real(kind=8)      epsbond,epsloc,epsdafh
c
c-----------------------------------------------------------------------
c
      call timer (2,iedfx,'_edfx     ',-1)
      write (lw,57) 
      read (lr,*) nmo,ncent,naom,ngroup
      if (.not.allocated(mocore)) then
        allocate (mocore(ngroup),stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot allocate mocore()'
      endif
      if (.not.allocated(nfugrp)) then
        allocate (nfugrp(ngroup),stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot allocate nfugrp()'
      endif
      if (.not.allocated(ifugrp)) then
        allocate (ifugrp(ncent,ngroup),stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot allocate ifugrp()'
      endif
      if (.not.allocated(eleca)) then
        allocate (eleca(2,ngroup),stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot allocate eleca()'
      endif
      do i=1,ngroup-1
        read (lr,'(a)') line
        lp=1
        ok=setint(nfugrp(i),line,lp)
        if (.not.ok) then
           write (stderr,*) 'edfx.f: Error defining fragment',i
           stop
        endif
        do j=1,nfugrp(i)
          ok=setint(ifugrp(j,i),line,lp)
          if (.not.ok) then
            write (stderr,*) 'edfx.f: Error defining fragment',i
            stop
          endif
          if (ifugrp(j,i).lt.izero) then
            stop 'edfx.f: Not valid read in atom'
          endif
          if (ifugrp(j,i).gt.naom) then
            write (stderr,*) '# Read in atom > naom = ',naom
            stop
          endif
        enddo
        ok=setint(eleca(1,i),line,lp)
        if (ok) then
          if (eleca(1,i).lt.izero) then
            eleca(1,i)=izero
            ok=setint(eleca(2,i),line,lp)
            if (ok) then
              if (eleca(2,i).lt.0)   eleca(2,i)=izero
              if (eleca(2,i).gt.nmo) eleca(2,i)=nmo
            else
              eleca(2,i)=nmo
            endif
          else
            eleca(1,i)=min(eleca(1,i),nmo)
            ok=setint(eleca(2,i),line,lp)
            if (ok) then
              if (eleca(2,i).lt.0)   eleca(2,i)=izero
              if (eleca(2,i).gt.nmo) eleca(2,i)=nmo
            else
              eleca(2,i)=nmo
            endif
          endif
        else
          eleca(1,i)=izero
          eleca(2,i)=nmo
        endif
        if (eleca(1,i).gt.eleca(2,i)) then
          write (stderr,*) 'edfx.f: Error defining fragment',i
          stop
        endif
      enddo
c
      dno      = .false.
      localize = .false.
      dodafh   = .false.
      largwr   = .false.
      epsbond  = 1D-7
      epsloc   = 1D-2
      epsdafh  = 1D-3
c
c.....Read control parameters used in the calculation.
c
      pcut=1D-3
      goon=.true.
      do while (goon) 
        read (lr,'(a)',end=110) line
        line=uppcase(line)
        lp=1
        ok=setword(word,line,lp)
        if (word(1:1).eq.'#') then
        elseif (word(1:leng(word)).eq.'END') then
          goon=.false.
        elseif (word(1:leng(word)).eq.'DNO') then
*         dno=.true.
          if (ngroup.eq.3) dno=.true.
        elseif (word(1:leng(word)).eq.'LARGE') then
          largwr=.true.
        elseif (word(1:leng(word)).eq.'FULLOC') then
          localize=.true.
        elseif (word(1:leng(word)).eq.'DAFH') then
*         if (ngroup.eq.2) dodafh=.true.
          dodafh=.true.
        elseif (word(1:leng(word)).eq.'EPSBOND') then
          line=line(lp:)
          read (line,*) epsbond
          epsbond=max(min(abs(epsbond),1D-3),1D-10)
        elseif (word(1:leng(word)).eq.'EPSDAFH') then
          line=line(lp:)
          read (line,*) epsdafh
          epsdafh=max(min(abs(epsdafh),1D-2),1D-6)
        elseif (word(1:leng(word)).eq.'EPSLOC') then
          line=line(lp:)
          read (line,*) epsloc
          epsloc=max(min(abs(epsloc),1D-2),1D-4)
        elseif (word(1:leng(word)).eq.'PROBCUT') then
          line=line(lp:)
          read (line,*) pcut
        else
          write (stderr,'(a)') ' # '//word(1:leng(word))
          write (stderr,*) '# Key word in input file not understood'
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
      nfulast=izero
      do l=1,ncent
        notingroup=.true.
        do i=1,ngroup-1
          do j=1,nfugrp(i)
            if (l.eq.ifugrp(j,i)) notingroup=.false.
          enddo
        enddo
        if (notingroup) then
          nfulast=nfulast+1
          ifugrp(nfulast,ngroup)=l
        endif
      enddo
      nfugrp(ngroup)=nfulast
c
      if (dno) then
        call dnoedf (nmo,ncent,naom,ngroup,nfugrp,ifugrp,pcut,
     &   epsbond,lr,lw,stderr,lu18)
         call timer (4,iedfx,'_edfx     ',-1)
         return
      endif
      if (dodafh) then
         if (ngroup.eq.2) 
     &   call dafh (nmo,ncent,naom,ngroup,nfugrp,ifugrp,pcut,
     &   epsdafh,lr,lw,stderr,lu18)
         if (ngroup.gt.2) 
     &   call gendafh (nmo,ncent,naom,ngroup,nfugrp,ifugrp,pcut,
     &   epsdafh,lr,lw,stderr,lu18)
         call timer (4,iedfx,'_edfx     ',-1)
         return
      endif
c
c.....Test that the group is not empty of atoms.
c
      do i=1,ngroup
        if (nfugrp(i).eq.0) then
          write (stderr,*) 'edfx.f:  GROUP ',i,' is empty of atoms.'
          stop
        endif
      enddo
c
c.....Test that each atom only belongs to a group.
c
      do i=2,ngroup
        do k=1,i-1
          do m=1,nfugrp(k)
            do l=1,nfugrp(i)
              if (ifugrp(l,i).eq.ifugrp(m,k)) then
                 write (stderr,*) 
                 write (stderr,432) l,i,m,k
                 write (stderr,*) 
                 stop
              endif
            enddo
          enddo
        enddo
      enddo
c
      minpop=izero
      maxpop=izero
      do i=1,ngroup-1
        minpop=minpop+eleca(1,i)
        maxpop=maxpop+eleca(2,i)
      enddo
      eleca(1,ngroup)=nmo-maxpop
      eleca(2,ngroup)=nmo-minpop
c
      if (.not.allocated(aomcmplx)) then
        allocate (aomcmplx(naom,nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot allocate aomcmplx()'
      endif
      if (.not.allocated(sgcmplx)) then
        allocate (sgcmplx(ngroup,nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot allocate sgcmplx()'
      endif
      if (.not.allocated(eleca)) then
        allocate (eleca(2,ngroup),stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot allocate eleca()'
      endif
c
      npab=1
      do i=1,ngroup-1
        npab=npab*(eleca(2,i)-eleca(1,i)+1)
      enddo
      nprev=npab
c
      if (.not.allocated(resncb)) then
        allocate (resncb(nprev,ngroup),stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot allocate resncb()'
      endif
c
c.....Computation of resonance structures.
c
      call rnprobs (eleca,resncb,nmo,ngroup,npab,lw)
      if (nprev.lt.npab) then
        write (stderr,*) 'edfx.f: NPREV < NPAB returned by rnprobs.f'
        stop
      endif

      if (.not.allocated(resnca)) then
        allocate (resnca(npab,ngroup),stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot allocate resnca()'
      endif
      resnca(1:npab,1:ngroup)=resncb(1:npab,1:ngroup)
c
c.....Read complex AOM
c
*     call readcmplx (aomcmplx,ncent,nmo,lu18,stderr)
      call readaom   (aomcmplx,naom,ncent,nmo,lu18,stderr)
c
c.....Compute complex Group Overlap integrals of all but the last one.
c
      sgcmplx(1:ngroup,1:nmo,1:nmo)=cmplx(zero,zero)
      do i=1,ngroup-1
        do j=1,nfugrp(i)
          k=ifugrp(j,i)
          sgcmplx(i,1:nmo,1:nmo)=
     &    sgcmplx(i,1:nmo,1:nmo)+aomcmplx(k,1:nmo,1:nmo)
        enddo
      enddo
c
c.....Compute complex Group Overlap integrals of the last group.
c
      do m=1,nmo
        aom1=0d0
        do i=1,ngroup-1
          aom1=aom1+dreal(sgcmplx(i,m,m))
        enddo
        sgcmplx(ngroup,m,m)=cmplx(one-aom1,zero)
      enddo
  
      do j=2,nmo
        do m=1,j-1
          aom1cmplx=cmplx(zero,zero)
          do i=1,ngroup-1
            aom1cmplx=aom1cmplx+sgcmplx(i,j,m)
          enddo
          sgcmplx(ngroup,j,m)=-aom1cmplx
          sgcmplx(ngroup,m,j)=conjg(sgcmplx(ngroup,j,m))
        enddo
      enddo
c
c.....Compute EDF using an isopycnic localization
c
      if (localize) then
        call locmplx (sgcmplx,nmo,ngroup,ncent,
     &      nfugrp,ifugrp,pcut,lw,epsloc,.true.)
         call timer (4,iedfx,'_edfx     ',-1)
         return
      endif
c
c.....Compute EDF without using the DNO or ISOPYCNIC modules.
c
      write (lw,766) nmo,ncent,naom,ngroup
      do i=1,ngroup
        write (lw,764) i,eleca(1,i),eleca(2,i)
        write (lw,767) i,nfugrp(i)
        write (lw,765) (ifugrp(j,i),j=1,nfugrp(i))
      enddo
 767  format (1x,'# Number of atoms in group ',I3,' = ',I3)
      mocore(1:ngroup)=izero
      call cmplxedf (pcut,npab,nprob,ngroup,nmo,lw,eleca,sgcmplx,
     &               resnca,mocore)
      if (allocated(nfugrp)) then
        deallocate (nfugrp,stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot deallocate nfugrp()'
      endif
      if (allocated(ifugrp)) then
        deallocate (ifugrp,stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot deallocate ifugrp()'
      endif
      if (allocated(aomcmplx)) then
        deallocate (aomcmplx,stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot deallocate aomcmplx()'
      endif
      if (allocated(sgcmplx)) then
        deallocate (sgcmplx,stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot deallocate sgcmplx()'
      endif
      if (allocated(eleca)) then
        deallocate (eleca,stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot deallocate eleca()'
      endif
      if (allocated(resnca)) then
        deallocate (resnca,stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot deallocate resnca()'
      endif
      if (allocated(resncb)) then
        deallocate (resncb,stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot deallocate resncb()'
      endif
      if (allocated(mocore)) then
        deallocate (mocore,stat=ier)
        if (ier.ne.0) stop 'edfx.f: Cannot deallocate mocore()'
      endif
      call timer (4,iedfx,'_edfx     ',-1)
      return
 766  format (
     & ' # Number of Complex Molecular Orbitals       = ',I4,/,
     & ' # Total number of atoms                      = ',I4,/,
     & ' # Number of Atoms for which the AOM is read  = ',I4,/,
     & ' # Number of groups that divide the molecule  = ',I4)
 764  format (' # GROUP ',I2,
     &  ' MinElec and MaxElec (alpha or beta) = ',2I6)
 765  format (' # ',20I5)
 57   format (1x,'# Complex Atomic Overlap Matrix is read in')
 432  format (1x,'# Atom number ',I3,' of group ',I2,
     &        1x,'is equal to atom number ',I3,' of group ',I2)
      end
