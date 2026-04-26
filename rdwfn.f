c
c-----------------------------------------------------------------------
c
      SUBROUTINE RDWFN (linp,uout,stderr,ifilc,mal,mbe,wfnfile,cicoef)
c
c.....read the wave function with a slightly modified AIMPAC routine.
c
      USE       space_for_wfnbasis
      USE       space_for_wfncoef
      USE       space_for_cidet
      USE       space_for_primgto
      USE       space_for_rdm1
      USE       space_for_wfxedf
      implicit  real(kind=8) (a-h,o-z)
      include  'param.inc'
      include  'wfn.inc'
      include  'corr.inc'
      include  'lengrec.inc'
      include  'fact.inc'
      include  'primgto.inc'
      include  'mline.inc'
      parameter (tolocc=1.0d-6)
c
      integer leng,stderr,uout
      character(len = *) wfnfile
      character(len = mline)  line
      character(len = mline)  cicoef
      character(len = mline)  uppcase
c
c.....local variables
c
      character (len = 80) :: wfnttl,jobttl,word
      character (len =  4) :: mode
      character (len =  8) :: check
      character (len = 17) :: label
      character (len =  5) :: l15

      logical setint, setword, ok, setdble, nonum, onlysingly
      character(len = 1) ch

      real(kind=8),    allocatable,dimension (:)   :: oexpa
      real(kind=8),    allocatable,dimension (:,:) :: coefa
      integer,         allocatable,dimension (:)   :: icena,itypa

      integer          ntfu(0:4)
      character(len=1) iotipe(0:4)
      data  iotipe(0) /'S'/, iotipe(1) /'P'/, iotipe(2) /'D'/,
     &      iotipe(3) /'F'/, iotipe(4) /'G'/
c
      iwfn=linp
      icorr=.false.
c
      read (iwfn,101) wfnttl
      read (iwfn,102) mode,nmo,nprims,ncent
c
      call allocate_space_for_wfncoef (nprims,nmo)
      call allocate_space_for_wfnbasis (nprims,nmo,ncent)
      call allocate_space_for_primgto (nprims,ncent,mgrp,ngtoG)
      call allocate_space_for_rdm1 (nmo)
      do i = 1,ncent
        read (iwfn,'(a)') line
        read (line,'(a8)') atnam(i)
        line = line(24:)
        lp = 1
        ok = .true.
        ok = setdble(xyz(i,1),line,lp) .and. ok
        ok = setdble(xyz(i,2),line,lp) .and. ok
        ok = setdble(xyz(i,3),line,lp) .and. ok
        if (.not.ok) then
          backspace (iwfn)
          read (iwfn,'(a)') line
          line = line(25:)
          read (line,'(3F12.8,10x,F5.1)') (xyz(i,k),k=1,3),charge(i)
        else
          ch=line(lp:lp)
          nonum = ch.eq.'+'.or.ch.eq.'-'.or.ch.eq.'.'.or.
     &           (ch.ge.'0'.and.ch.le.'9')
          do while ((.not.nonum).and.lp.lt.leng(line))
            lp = lp+1
            ch = line(lp:lp)
            nonum = ch.eq.'+' .or. ch.eq.'-' .or. ch.eq.'.'
     &        .or. (ch.ge.'0'.and.ch.le.'9')
          enddo
          ok = setdble(charge(i),line,lp)
          if (.not.ok) write (0,182) i+2
          if (.not.ok) stop
        endif
      enddo
C
C-----------------------------------------------------------------------
C
      read (iwfn,104) (icen(i),i=1,nprims)
      read (iwfn,104) (ityp(i),i=1,nprims)
      read (iwfn,105) (oexp(i),i=1,nprims)
c
c-----Here, we read the Electron Density Function (EDF) associated to the
c     core electron with the same format used in the AIMALL code
c
*     nrrec = 0
*     call readedf (iwfn,nrrec,thereisEDF)
      call readedf (iwfn,thereisEDF)
c
c-----Print the function defining the core electron density
c
      call writeedf (uout,thereisEDF)
c
      do i = 1,nprims
        if(ityp(i).GT.maxtype) then 
          stop '# rdwfn.f: EDF cannot work with L>=5 primitives'
        endif
      enddo
c
c.....Fill the nlm() array
c
      call nlmfill ()
C
C-----read occ(), eorb() and MO coefficients in the primitive basis
C
      rhf      = .false.
      rohf     = .false.
      uhf      = .false.
      cciqa    = .false. ! By default WFN is not a Coupled Cluster one
      onlynat  = .false. ! By default WFN contains canonical and natural MOs
      ndou     = 0
      nsin     = 0
      nalpha   = 0
      nbeta    = 0
      mal      = 0
      mbe      = 0
      totel    = 0.0d+00
      occmin   = 2.0d+00
      occmax   = 0.0d+00
      n        = 0
      onlysingly = .true.
      do i = 1,nmo
        read (iwfn,106) occ(i),eorb(i) 
        totel = totel + abs(occ(i))
        read (iwfn,107) (coef(i,j),j=1,nprims)
        forall (j=1:nprims) coef(i+nmo,j)=coef(i,j)
c
c-------Test if the WFN is made exclusively of singly occupied MOs
c
        if (abs(occ(i)-1.0d+00).gt.tolocc) onlysingly = .false.
      enddo

      occmin = minval(occ(1:nmo))
      occmax = maxval(occ(1:nmo))
      if (abs(occmax-2.0d+00).lt.tolocc .and.
     &    abs(occmin-2.0d+00).lt.tolocc) then
        nsin    = 0
        ndou    = nmo
        nalpha  = nmo
        nbeta   = nmo
        mal     = nmo
        mbe     = nmo
        forall (i=1:nmo) ialpha(i)=i
        forall (i=1:nmo) ibeta (i)=i
        rhf = .true.
      elseif (abs(occmax-2.0d+00).lt.tolocc .and.
     &        abs(occmin-1.0d+00).lt.tolocc) then
        ndou = 0
        nsin = 0
        nalpha = 0
        nbeta  = 0
        do i=1,nmo 
          if (abs(occ(i)-2.0d+00).lt.tolocc) then
            ndou = ndou + 1
            nalpha = nalpha + 1
            nbeta  = nbeta + 1
            ialpha(nalpha) = i
            ibeta(nbeta) = i
          endif         
          if (abs(occ(i)-1.0d+00).lt.tolocc) then
            nsin = nsin + 1
            nalpha = nalpha + 1
            ialpha(nalpha) = i
          endif
        enddo
        if (ndou+nsin.ne.nmo) then
          stop '# rdwfn.f: Bad number of alpha and beta MOs'
        endif
        mal = nalpha
        mbe = nbeta
        rohf = .true.
      elseif (abs(occmax-1.0d+00).lt.tolocc .and.
     &        abs(occmin+1.0d+00).lt.tolocc) then
        nalpha = 0
        nbeta  = 0
        nsin   = nmo
        ndou   = 0
        do i=1,nmo
          if (occ(i).gt.0.0d+00) then
            nalpha = nalpha+1
            ialpha(nalpha) = i
          else
            nbeta = nbeta+1
            ibeta(nbeta) = i
          endif
        enddo
        mal = nalpha
        mbe = nbeta
        uhf = .true.
      elseif (abs(occmax-2.0d+00).lt.tolocc .and.
     &        abs(occmin-0.0d+00).lt.tolocc) then
        cciqa = .true.
      else
        onlynat = .true.
        write (uout,*) '# WFN file with exclusively natural MOs read in'
      endif

      read (iwfn,108) check
      if (uppcase(check) .ne. 'END DATA') then
         stop '# rdwfn.f : End card 1 not found '
      endif
c
c    Read in total scf energy and -V/T
c
      read (iwfn,'(a)') line
      lp = 1
      if (.not.setword(label,line,lp)) then
        write (stderr,*) '#'
        write (stderr,*) '# Error reading "label,tote,gamma" from input'
        write (stderr,*) '#'
        stop
      endif
      l15 = uppcase(label(1:5))
      if (l15.eq.'MCSCF' .or. l15.eq.'ALDET' .or. 
     &    l15.eq.'ORMAS' .or. l15.eq.'GENCI') then
c
c.......correlated case. Equal number of NMO
c
c-------Since this is a multideterminant wavefunction,
c       it is not a wavefunction with only natural MOs
c
        onlynat = .false.
c
        icorr=.true.
        read (iwfn,'(a)',end=9999) label
        do i=nmo+1,nmo+nmo
           read (iwfn,*) label
           read (iwfn,107) (coef(i,j),j=1,nprims)
        enddo
        read (iwfn,108) check
        if (uppcase(check) .ne. 'END DATA') then
          stop '# rdwfn.f : end card 2 not found '
        endif
        read (iwfn,*) label
        read (iwfn,*) nelact, ndets, ncore, nact
c
c.......Allocate space for cidet()
c
        call allocate_space_for_cidet (nelact+2*ncore)
c
        nel = nelact + 2 * ncore
c
c.......Read coefficients and occupied spin-MO's of the determinants
c       and write them in a direct access file.
c
        cicoef=wfnfile(1:leng(wfnfile))//".CIcoef"
        open (ifilc,file=cicoef,access='direct',
     &      recl=RecLength*(nelact+1),form='unformatted')
        read (iwfn,*) label
        do i=1,ndets
           read (iwfn,*) (cidet(j),j=0,nelact)
           write (ifilc,rec=i) (cidet(j),j=0,nelact)
           if (i.eq.1) then
             mal = ncore
             mbe = ncore
             do m=1,nelact
               mm=nint(cidet(m))
               if (mm.gt.0) mal=mal+1
               if (mm.lt.0) mbe=mbe+1
             enddo
           endif
        enddo
c
c-------In Multi_Det wave functions the number of ALPHA and BETA
c       electrons are exactly MAL and MBE. However, it is appropriate 
c       DO NOT initialize NALPHA and NBETA to these values, but lo leave 
c       them both equal to NMO, the total number of MOs, even thoug some
c       MOs may be singly occupied in some of the determinants. This is 
c       the reason why the folloging two lines are commented.
c
*       nalpha = mal
*       nbeta = mbe
        nalpha = nmo
        nbeta  = nmo
        if (nel.ne.mal+mbe) then
          stop 'rdwfn:  !!!!! (e- alpha + e- beta).ne.nelact !!!!! '
        endif
        forall (i=1:nalpha) ialpha(i) = i
        forall (i=1:nbeta) ibeta(i) = i
        close (ifilc)
      elseif (uppcase(label(1:5)).eq.'CCIQA' .or. 
     &        uppcase(label(1:3)).eq.'MP2') then
        nel = int(anint(totel))
c
c-------Since this is a multideterminant wavefunction,
c       it is not a wavefunction with only natural MOs
c
        onlynat = .false.
        icorr = .false.
        nelact = nel
        ncore = 0
        nact = nmo
        c1et = 0.0d+00
        read (iwfn,'(a)') label        
 125    read (iwfn,*,err=123) i,j,c1et(i,j)
        goto 125
 123    continue
        trzrdm=0.0d+00
        do i=1,nmo
          do j=1,i
            value=(c1et(i,j)+c1et(j,i))*0.5d0
            c1et(i,j)=value
            c1et(j,i)=value
          enddo
          trzrdm=trzrdm+c1et(i,i)
        enddo
        ntrz=int(anint(trzrdm))
        if (nel.ne.ntrz) then
          stop '# rdwfn.f: Inconsistent number of e- in the WFN file'
        endif
        if (mod(ntrz,2).ne.0) then
          stop '# rdwfn.f: Open-shell coupled cluster WFNs not allowed'
        else
          c1ea = c1et*0.5d0
          c1eb = c1et*0.5d0
          nalpha = nmo
          nbeta = nmo
          forall (i=1:nalpha) ialpha(i)=i
          forall (i=1:nbeta) ibeta(i)=i
          mal = ntrz/2
          mbe = ntrz/2
        endif
      elseif (onlynat) then
c
c-------We have now an RHF, ROHF or UHF wave function. We transform
c       them to the Multi-Determinant case with NDETS = 1
c
      else
        if (onlysingly) then
          write (stderr,400) tolocc
          stop
        endif
        icorr = .false.
        ndets = 1
c
c.......Test that the number of electrons is integer
c
        if (abs(totel-nint(totel)).gt.tolocc) then
          write (stderr,*) 
          write (stderr,*) 'rdwdn: ! Number of electrons is not integer'
          write (stderr,*) 
        endif
        nel = nint(totel)
        nelact = nel
        ncore = 0
        nact = nmo
c
c.......Allocate space for cidet()
c
        call allocate_space_for_cidet (nelact)
        cidet(0)=1.0d+00
c
        if (.not.(rhf.or.rohf.or.uhf)) then
          write (stderr,*) 
          write (stderr,'(a)') ' # rdwdn: !!! Monodeterminantal '//
     &     'WFN which is not of RHF, ROHF, or UHF type'
          write (stderr,*) 
        endif
c
        cidet(1:nalpha) = ialpha(1:nalpha)
        forall (i=1:nbeta) cidet(i+nalpha) = -ibeta(i)
        if (nalpha+nbeta.ne.nelact) then
          stop 'rdwfn: Improper number of electrons '
        endif
c
        if (ndou+ndou+nsin.ne.nelact) then
           stop 'rdwfn: Improper number of electrons '
        endif
        cicoef=wfnfile(1:leng(wfnfile))//".CIcoef"
        open (ifilc,file=cicoef,access='direct',
     &      recl=RecLength*(nelact+1),form='unformatted')
        write (ifilc,rec=1) (cidet(j),j=0,nelact)
        if (nel.ne.mal+mbe) then
          stop 'rdwfn:  !!!!! (e- alpha + e- beta).ne.nel !!!!! '
        endif
        close (ifilc)
      endif
 9999 continue
c
      rewind(iwfn)
c
c-----MOs in the WFN file are expressed in terms of unnormalized carte-
c     sian gaussians (CG). It could happens that a given CG is repeated.
c     In that case, we condensate all the coefficients in a single one.
c     
      if (.not.allocated(oexpa)) then
        allocate (oexpa(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate oexpa()'
      endif
      if (.not.allocated(coefa)) then
        allocate (coefa(nmo+nmo,nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate coefa()'
      endif
      if (.not.allocated(icena)) then
        allocate (icena(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate icena()'
      endif
      if (.not.allocated(itypa)) then
        allocate (itypa(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate itypa()'
      endif
      npa=0
      do 10 j=1,nprims
        if (j.eq.1) then
          npa=npa+1
          icena(npa)=icen(j)
          oexpa(npa)=oexp(j)
          itypa(npa)=ityp(j)
          do i=1,nmo+nmo
            coefa(i,npa)=coef(i,j)
          enddo
        else
          do m=1,npa
            if (icen(j).eq.icena(m) .and.
     &          ityp(j).eq.itypa(m) .and.
     &          abs(oexp(j)-oexpa(m)).le.1d-10)  then
                do i=1,nmo+nmo
                  coefa(i,m)=coefa(i,m)+coef(i,j)
                enddo
                goto 10
            endif
          enddo
          npa=npa+1
          icena(npa)=icen(j)
          oexpa(npa)=oexp(j)
          itypa(npa)=ityp(j)
          do i=1,nmo+nmo
            coefa(i,npa)=coef(i,j)
          enddo
        endif
 10   continue
c
c-----Recompute the original variables
c
      write (uout,*)'# '
      write (uout,311) nprims,npa
 311  format (1x,'# Input number of Primitives ',I8,' reduced to ',I8)
      nprims=npa
c
c-----Deallocate and allocate again those arrays depending on 'nprims',
c     whose dimensions might have changed with respect to their initial 
c     values due to the above line 'nprims=npa'
c
      if (nprims.ne.npa) then
        deallocate (sprim,oexp,icen,ityp)
        allocate (sprim(nprims,nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate sprim() array ?'
        allocate (oexp(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate oexp() array ?'
        allocate (icen(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate icen() array ?'
        allocate (ityp(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate ityp() array ?'
        deallocate (coef)
        allocate (coef(nmo+nmo,nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate coef() array ?'
        deallocate (xnorm,icenat)
        allocate (xnorm(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate xnorm() array ?'
        allocate (icenat(nprims,ncent),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Can not allocate icenat() array ?'
      endif
      do j=1,nprims
        icen(j)=icena(j)
        oexp(j)=oexpa(j)
        ityp(j)=itypa(j)
        do i=1,nmo+nmo
          coef(i,j)=coefa(i,j)
        enddo
      enddo
c
c.....Determine primitives corresponding to each center.
c
      do ic=1,ncent
        npc(ic)=0
      enddo
      do j=1,nprims
        ic=icen(j)
        npc(ic)=npc(ic)+1
        inda=npc(ic)
        icenat(inda,ic)=j
      enddo
c
c.....Classify primitives in each center by types and similar exponents.
c
      ngroup(1:ncent)=0
      do ic=1,ncent
        if (npc(ic).eq.0) cycle
        do 33 j=1,npc(ic)
          k=icenat(j,ic)
          itk=ityp(k)
          isuk=nlm(itk,1)+nlm(itk,2)+nlm(itk,3)
          if (j.eq.1) then
            ngroup(ic)=1
            nzexp(ic,1)=1
            nuexp(ic,1,1)=k
          else
            do m=1,ngroup(ic)
              inda=nuexp(ic,m,1)
              itd=ityp(inda)
              isud=nlm(itd,1)+nlm(itd,2)+nlm(itd,3)
              if (abs(oexp(k)-oexp(inda)).lt.1d-8) then
                if (itk.eq.1.and.itd.eq.1) then
                  write (uout,*)
     &            'Warning TWO S primitives with equal exponents'
                   stop
                else
                  if (isuk.eq.isud) then
                    nzexp(ic,m)=nzexp(ic,m)+1
                    nuexp(ic,m,nzexp(ic,m))=k
                    goto 33
                  endif
                endif
              endif
            enddo
            ngroup(ic)=ngroup(ic)+1
            if (ngroup(ic).gt.mgrp) then
              write (0,*) 'Faterr: Increase MGRP in primgto.inc file'
              stop
            endif
            nzexp(ic,ngroup(ic))=1
            nuexp(ic,ngroup(ic),1)=k
          endif
 33     continue
      enddo
      i = 0
      do ic=1,ncent
        do m=1,ngroup(ic)
          do k=1,nzexp(ic,m)
            j = nuexp(ic,m,k)
            alph = oexp(j)
            itip = ityp(j)
            icentro = icen(j)
            i = i + 1
            itypa (i) = itip
            oexpa (i) = alph
            icena (i) = icentro
            do n=1,nmo+nmo
              coefa(n,i) = coef(n,j)
            enddo
          enddo
        enddo
      enddo
c 
      do i = 1, nprims
        ityp(i) = itypa (i)
        oexp(i) = oexpa (i)
        icen(i) = icena (i)
        do n=1,nmo+nmo
           coef(n,i) = coefa(n,i)
        enddo
c
c.......Normalization constants of Cartesian Gaussians.
c
        it1=nlm(ityp(i),1)
        it2=nlm(ityp(i),2)
        it3=nlm(ityp(i),3)
        it1d=it1+it1-1
        it2d=it2+it2-1
        it3d=it3+it3-1
        isu=it1+it2+it3
        dblisu2=(isu+isu+3)*0.25d0
        xnorm(i)=(2.0d+00/pi)**0.75d0*2.0d+00**isu*oexp(i)**dblisu2
        xnorm(i)=xnorm(i)/sqrt(facd(it1d)*facd(it2d)*facd(it3d))
      enddo
c
      npcant = 0
      do ic=1,ncent
        if (npc(ic).eq.0) cycle
        do k=1,npc(ic)
          icenat(k,ic) = k + npcant
        enddo
        npcant = npcant + npc(ic)
      enddo
c
c.....Reconstruct the values of nuexp(). Now, they are ordered.
c.....Determine also which is the maximum value of ngroup(ic)
c
      i = 0
      maxgrp=0
      do ic=1,ncent
        if (ngroup(ic).eq.0) cycle
        if (ngroup(ic).gt.maxgrp) maxgrp=ngroup(ic)
        do m=1,ngroup(ic)
           do k=1,nzexp(ic,m)
             i = i + 1
             nuexp(ic,m,k) = i
           enddo
        enddo
      enddo
c
c-----nprimc(ic) is the number of primitives of center ic.
c
      write (uout,*)'# '
      write (uout,*)'# Description of the Primitive Basis Set'
      write (uout,604) nprims
      do ic=1,ncent
        nprimc(ic)=0
        if (ncent.lt.20) write (uout,210) ic
        ntfu(0:4) = 0
        do m=1,ngroup(ic)
          nprimc(ic)=nprimc(ic)+nzexp(ic,m)
          zz=oexp(nuexp(ic,m,1))
          ii=nuexp(ic,m,1)
          itip=ityp(ii)
          it1=nlm(itip,1)
          it2=nlm(itip,2)
          it3=nlm(itip,3)
          isu=it1+it2+it3
          ntfu(isu)=ntfu(isu)+1
          nzicm=nzexp(ic,m)
          if (ncent.le.20) then
            write (uout,211) iotipe(isu),zz,(nuexp(ic,m,k),k=1,nzicm)
          endif
        enddo
        if (ncent.le.20) then
          write (uout,212) (ntfu(nn),nn=0,4),nprimc(ic)
        else
          write (uout,2120) ic,(ntfu(nn),nn=0,4),nprimc(ic)
        endif
      enddo
c
c.....Deallocate arrays.
c
      if (allocated(oexpa)) then
        deallocate (oexpa,stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot deallocate oexpa()'
      endif
      if (allocated(coefa)) then
        deallocate (coefa,stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot deallocate coefa()'
      endif
      if (allocated(icena)) then
        deallocate (icena,stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot deallocate icena()'
      endif
      if (allocated(itypa)) then
        deallocate (itypa,stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot deallocate itypa()'
      endif

      return
c
 101  format (A80)
 102  format (4X,A4,10X,3(I5,15X))
 103  format (A8,11X,I3,2X,3F12.8,10X,F5.1)
 104  format (20X,20I3)
 105  format (10X,5E14.7)
 106  format (35X,F12.8,15X,F12.8)
 107  format (5E16.8)
 108  format (A8)
 109  format (a17,F20.12,18X,F13.8)
 2000 format (' rdwfn: !!! WFN file is not of RHF, MCSCF, ALDET 
     &  or GENCI type type   !!!')
 604  format (' # Total number of Primitive Gaussians: ',I6)
 210  format (1x,'# CENTER ',I3)
 211  format (1x,'# ',a,' Shell (Z=',e16.10,') : ',30I4)
 212  format (1x,'# This seems to be a [ ',I2,'s | ',I2,'p |',I2,
     &  'd | ',I2,'f | ',I2,'g ] basis, Total primitives = ', I4)
 2120 format (1x,'# CENTER ',I3,' Basis = [ ',I2,'s | ',I2,'p |',I2,
     &   'd | ',I2,'f | ',I2,'g ], Total primitives = ', I4)
 2121 format (1x,'# Primitives of this center = ',I4)
 182  format (1x,'# rdwfn.f: Error reading charge() of line ',
     &  I2,' of WFN file')
 400  format (/,' # ',56('+'),/,
     & ' # This seems to be an UHF wave function, even though all',/,
     & ' # the molecular orbital (MO) occupations differ from 1.000',/,
     & ' # by less than ',F13.6,//,
     & ' # It is very recomendable to put "by hand" the occupations',/,
     & ' # of BETA MOS equal to -1.000000, or transform the WFN',/,
     & ' # in such a way that it looks like a Multi-Det WFN',/,
     & ' # ',56('+'))
 
      end
