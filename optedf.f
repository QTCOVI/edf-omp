c
c-----Driver routine of the optimization process used to approximate 
c     EDF as a direct product of (2c,2ne) chemical bond EDFs.
c
c     Within the input file of the 'edf' program the following keywords
c     can be given that are exclusively relevant of this routine:
c
c        1) ENDBONDING
c        2) TYPE
c        3) TYPE3C
c        4) PAIR
c        5) TRIO
c        6) EPSBOND
c        7) INICVAR
c        8) PRIN
c
c         ) WEDF
c         ) WPOP
c         ) WDIS
c
c         ) PLUS
c         ) MINUS
c
c     The meaning of all of these orders is explained below.
c
c     1) 'ENDBONDING'
c
c     Sintax: ENDBONDING
c
c     Meaing: This order means that all the input that the routine 
c     'optedf' needs to work has been given, so that the calculation 
c     can start.
c
c     2) 'TYPE'
c
c     Sintax: TYPE type qq(type) ifixq(type) ffty(type) ifixf(type)
c
c     Meaning: This order defines a type of (2c,2e) bond. The maximum
c     value of 'type' is MAXBOND, a parameter whose current value is
c     maxbond=18. To increase the maximum possible number of types of
c     (2c,2e) bonds just increase the MVAR, MVARH, and MAXBOND values 
c     in the sentence 'parameter (mvar=36,mvarh=mvar/2,maxbond=18)' 
c     below and recompile.
c
c     'qq(type)' is 'the charge value of the bond', which is a measure
c     of the polarity of the bond. Its value is optimized in case that
c      the read in value of 'ifixq(type)' is non zero. Otherwise, the
c     value of 'qq(type)' is fixed. 
c
c     'ffty(type)' if the 'correlation factor of the bond'. Its value 
c     is optimized in case that the read in value of 'ifixf(type)' is
c     non zero. Otherwise, the value of 'ffty(type)' is fixed.
c
c     NO default values exist for 'qq(type)' and 'ffty(type)' and the
c     five values after the 'TYPE' keyword are mandatory. Tipical input
c     values for 'qq(type)' and 'ffty(type)' are 0.5 and 0.0,
c     respectively.
c
c     3) 'TYPE3C'
c
c     Sintax: TYPE3C type p200 p020 p002 p110 p101 f200 f020 f002 f110 f101
c
c     Meaning: This order defines a type of (3c,2e) bond. The maximum
c     value of 'type' is MAXBOND3C, a parameter whose current value is
c     maxbond3c=6. To increase the maximum possible number of types of
c     (3c,2e) bonds just increase the MAXBOND3C and MVAR3C values in the 
c     sentence 'parameter (mvar3c=30,maxbond3c=6)' below and recompile. 
c
c     'p200 p020 p002' are the input values for the probability that 
c     both electrons of the 3c bond are located in the first center (p200),
c     in the second center (p020), or in the third center (p002).
c
c     'p110' is the input value for the probability that one of
c     the electrons of the 3c bond is located in the first center and 
c     the second electron in the second center.
c
c     'p101' is the input value for the probability that one of
c     the electrons of the 3c bond is located in the first center and 
c     the second electron in the third center.
c
c     'f200 f020 f002 f110 f101' are five integer numbers that determine
c     if 'p200 p020 p002 p110 p101' are fixed or optimizable parameters.
c     For instance, f200=1 or 0 means that 'p200' is fixed or optimizable
c     respectively.
c
c      4) 'PAIR'
c
c      Sintax: PAIR fragment1 fragment2 type [ rtype ]
c
c      Meaning: This order indicates that a (2c,2e) bond of type 'type'
c      purportedly exists between the groups of atoms 'fragment1' and 
c      'fragment2'. 

c      'rtype' indicates the resonance structure used. The default is
c      rtype=1. It is recommended at this moment NOT TO ENTER THE
c      rtype VALUE IN THIS ORDER.
c
c     5) 'TRIO'
c
c     Sintax: TRIO fragment1 fragment2 fragment3 type [ rtype ]
c
c     Meaing: This order indicates that a (3c,2e) bond of type 'type'
c     purportedly exists between the groups of atoms 'fragment1',
c     'fragment2', and 'fragment3'.

c     'rtype' indicates the resonance structure used. The default is
c     rtype=1. It is recommended at this moment NOT TO ENTER THE
c     rtype VALUE IN THIS ORDER.
c
c     6) 'EPSBOND'
c
c     Sintax: EPSBOND epsbond
c
c     Convergence threshold in the EDF fitting. This order is not nece-
c     ssary, since a default value for 'epsbond' exists [10^(-4)].
c     
c
c     7) 'INICVAR'
c
c     Sintax: INICVAR h0
c
c     Meaning: 'h0' is a parameter used in the EDF fitting. Not necessary 
c     since a default value exists (0.1).
c
c     8) PRIN prin
c
c     Sintax: PRIN prin
c
c     Meaning: 'prin' is a printing level in the EDF fitting. Not necessary
c     since a default value exists (0).
c
c
c-----------------------------------------------------------------------
c
      subroutine optedf (
     &  epsdet,probcut,ndets,nelact,nact,ncent,nmo,ncore,ngroup,mal,
     &  mbe,mocore,nel,moval,ifilc,lw,lr,aom,sg,nfugrp,ifugrp,mocogrp,
     &  ival,short,wfnf)
      USE space_for_rsrs
      USE space_for_bonds
      USE space_for_wfnbasis

      include 'implicit.inc'
      include 'param.inc'
      include 'constants.inc'
      include 'mline.inc'
      parameter (mvar=36,mvarh=mvar/2,maxbond=18)
      parameter (mvar3c=30,maxbond3c=6)
      parameter (mvartot=mvar+mvar3c)

      real(kind=8)  aom(ncent,nmo,nmo)
      real(kind=8)  sg(ngroup,nmo,nmo)
      integer nfugrp(ngroup)
      integer ifugrp(ncent,ngroup)
      integer mocogrp(ngroup)
      integer ival(nmo)
*     parameter  (mline=200)
      character*(mline) wfnf
      character*(3) choq,chof,chot
      integer, allocatable,dimension (:,:)   :: eleca,elecb
      logical inlist,ok,ok1,ok2,ok3,ok4
      logical setword,setint,setdble,ditwocent
      logical plus,minus
      character*(mline) line,word,uppcase
      integer f200,f020,f002,f110,f101

      logical deftype(mvarh),defpair,thistyp,short,prevtyp
      logical deftype3c(maxbond3c)
      real(kind=8) t0,machep,h0,xvar(mvartot),fmin,charval,fffac
      real(kind=8) praxis,gedfit,edfit,edfdif
      integer prin,g1,g2,g3,typ,typ1,leng,rtype

      external gedfit
c
c-----------------------------------------------------------------------
c
      call timer (2,iopedf,'_optedf   ',-1)
c
c-----Read bonding pattern data from the input file.
c
      rewind (lr)
      largwr = .not.short
      plus   = .false.
      minus  = .false.
      nfrag  = ngroup
      lwrite = lw
c
c-----Default weights in the optimization function
c
      weigedf = 1D0
      weigpop = 0D0
      weigdis = 0D0
      call aspace_for_bonds (maxbond,mvarh,mvar,mvar3c,maxbond3c)
      mres    = maxres
      qrest(1:nfrag,1:maxres) = 0
      nbonds(1:maxres)        = 0
      nbonds3c(1:maxres)      = 0
      deftype(1:mvarh)        = .false.
      ditwocent               = .false.
      t0                      = 1D-4
      h0                      = 0.1D0
      prin                    = 0

 20   read (lr,'(a)',end=110) line
      line=uppcase(line)
      lp=1
      ok=setword(word,line,lp)
      if (word(1:1).eq.'#') then
      elseif (word(1:8).eq.'ENDBONDING') then
        goto 110
      elseif (word(1:leng(word)).eq.'EPSBOND') then
        if (setdble(t0,line,lp)) t0=abs(t0)
      elseif (word(1:leng(word)).eq.'WEDF') then
        if (setdble(weigedf,line,lp)) weigedf=abs(weigedf)
      elseif (word(1:leng(word)).eq.'WPOP') then
        if (setdble(weigpop,line,lp)) weigpop=abs(weigpop)
      elseif (word(1:leng(word)).eq.'WDIS') then
        if (setdble(weigdis,line,lp)) weigdis=abs(weigdis)
      elseif (word(1:leng(word)).eq.'INICVAR') then
        if (setdble(h0,line,lp)) h0=abs(h0)
      elseif (word(1:leng(word)).eq.'PRIN') then
        if (setint(prin,line,lp)) prin=abs(prin)
      elseif (word(1:leng(word)).eq.'TYPE') then
        line=line(lp:)
        lp=1
        ok=setint(typ,line,lp)
        if (ok) ok=ok.and.typ.gt.0.and.typ.le.maxbond
        if (ok) then
          deftype(typ)=.true.
          line=line(lp:)
          lp=1
          read (line,*,iostat=ios) qqty(typ),ifixq,ffty(typ),ifixf
          if (ios.ne.0) then
            stop ' # optedf.f: Incomplete or erroneus TYPE order'
          endif
          if (ifixq.eq.0) then
           fixeq(typ)=.true.
          else
           fixeq(typ)=.false.
          endif
          if (ifixf.eq.0) then
           fixef(typ)=.true.
          else
           fixef(typ)=.false.
          endif
        endif
      elseif (word(1:leng(word)).eq.'TYPE3C') then
        line=line(lp:)
        lp=1
        ok=setint(typ,line,lp)
        if (ok) ok=ok.and.typ.gt.0.and.typ.le.maxbond3c
        if (ok) then
          deftype3c(typ)=.true.
          line=line(lp:)
          lp=1
          read (line,*,iostat=ios) 
     &      p200,p020,p002,p110,p101,f200,f020,f002,f110,f101
          if (ios.ne.0) then
            stop ' # optedf.f: Incomplete or erroneus TYPE3C order'
          endif
          p200=abs(p200)
          p020=abs(p020)
          p002=abs(p002)
          p110=abs(p110)
          p101=abs(p101)
          prob3c(typ,1)=p200
          prob3c(typ,2)=p020
          prob3c(typ,3)=p002
          prob3c(typ,4)=p110
          prob3c(typ,5)=p101
          do k=1,5
            if (prob3c(typ,k).gt.1d0) then
              stop ' # optedf.f: TYPE3C order, probability > 1'
            endif
          enddo
          if (f200.eq.0) then
           fixet(typ,1)=.true.
          else
           fixet(typ,1)=.false.
          endif
          if (f020.eq.0) then
           fixet(typ,2)=.true.
          else
           fixet(typ,2)=.false.
          endif
          if (f002.eq.0) then
           fixet(typ,3)=.true.
          else
           fixet(typ,3)=.false.
          endif
          if (f110.eq.0) then
           fixet(typ,4)=.true.
          else
           fixet(typ,4)=.false.
          endif
          if (f101.eq.0) then
           fixet(typ,5)=.true.
          else
           fixet(typ,5)=.false.
          endif
        endif
      elseif (word(1:leng(word)).eq.'PAIR') then
        line=line(lp:)
        lp=1
        ok1=setint(g1,line,lp)
        ok2=setint(g2,line,lp)
        ok3=setint(typ,line,lp)
        if (ok1) ok1=ok1.and.g1.gt.0.and.g1.le.ngroup
        if (ok2) ok2=ok2.and.g2.gt.0.and.g2.le.ngroup
        if (ok3) ok3=ok3.and.typ.ge.1.and.typ.le.mvarh
        if (ok3) ok3=ok3.and.deftype(typ)
        if (ok1.and.ok2.and.ok3) then
          if (g1.eq.g2) then
            stop 'optedf.f: PAIR order: Both groups can not be equal'
          elseif (g1.gt.g2) then
            iming=min(g1,g2)
            imaxg=max(g1,g2)
            g1=iming
            g2=imaxg
          endif
          ok=setint(rtype,line,lp)
          if (ok) then
            ok=ok.and.rtype.gt.0.and.rtype.le.maxres
            if (.not.ok) stop 'optedf.f: Improper resonance structure'
          else
            rtype=1
          endif
          if (nbonds(rtype).lt.maxbond) then
              nbonds(rtype)=nbonds(rtype)+1
              ipair(nbonds(rtype),1,rtype)=g1
              ipair(nbonds(rtype),2,rtype)=g2
              ipair(nbonds(rtype),3,rtype)=typ
          else
            stop 'optedf.f: Too many bonds. Increase maxbond'
          endif
        endif
      elseif (word(1:leng(word)).eq.'TRIO') then
        line=line(lp:)
        lp=1
        ok1=setint(g1,line,lp)
        ok2=setint(g2,line,lp)
        ok3=setint(g3,line,lp)
        ok4=setint(typ,line,lp)
        if (ok1) ok1=ok1.and.g1.gt.0.and.g1.le.ngroup
        if (ok2) ok2=ok2.and.g2.gt.0.and.g2.le.ngroup
        if (ok3) ok3=ok3.and.g3.gt.0.and.g3.le.ngroup
        if (ok4) ok4=ok4.and.typ.ge.1.and.typ.le.maxbond3c
        if (ok4) ok4=ok4.and.deftype3c(typ)
        if (ok1.and.ok2.and.ok3.and.ok4.and.
     &      g1.lt.g2.and.g1.lt.g3.and.g2.lt.g3) then
          ok=setint(rtype,line,lp)
          if (ok) then
            ok=ok.and.rtype.gt.0.and.rtype.le.maxres
            if (.not.ok) stop 'optedf.f: Improper resonance structure'
          else
            rtype=1
          endif
          if (nbonds3c(rtype).lt.maxbond3c) then
              nbonds3c(rtype)=nbonds3c(rtype)+1
              itrio(nbonds3c(rtype),1,rtype)=g1
              itrio(nbonds3c(rtype),2,rtype)=g2
              itrio(nbonds3c(rtype),3,rtype)=g3
              itrio(nbonds3c(rtype),4,rtype)=typ
          else
            stop 'optedf.f: Too many 3C bonds. Increase maxbond3c'
          endif
        else
            stop 'optedf.f: Wrong TRIO order'
        endif
      elseif (word(1:4).eq.'PLUS'.or.word(1:4).eq.'MINU') then
        if (word(1:4).eq.'PLUS') then
          nn=+1
          plus=.true.
        endif
        if (word(1:4).eq.'MINU') then 
          nn=-1
          minus=.true.
        endif
        line=line(lp:)
        lp=1
        ok=setint(g1,line,lp)
        if (ok) then
          ok=ok.and.g1.gt.0.and.g1.le.ngroup
          if (.not.ok) then
            stop 'optedf.f: Bad value of GROUP in PLUS or MINUS order'
          else
            ok1=setint(isthisres,line,lp)
            if (ok1) then
              ok1=ok1.and.isthisres.gt.0.and.isthisres.le.maxres
              if (ok1) then
                qrest(g1,isthisres)=qrest(g1,isthisres)+nn
              else
                stop 'optedf.f: Bad value of RES in PLUS or MINUS order'
              endif
            else
              qrest(g1,1)=qrest(g1,1)+nn
            endif
          endif
        else
          stop 'optedf.f: Erroneus format in PLUS or MINUS order'
        endif
      elseif (word(1:leng(word)).eq.'TWOCENDI') then
        ditwocent = .true.
        epsbond=1D0
        line=line(lp:)
        read (line,*,err=339) epsbond
 339    continue
      else
      endif
      goto 20
 110  continue
c
c-----Number of resonance structures actually used.
c
      nonzres   = 0
      nbondgt   = 0
      nbondgt3c = 0
      do rtype=1,maxres
        if (nbonds(rtype).gt.0.or.nbonds3c(rtype).gt.0) then
           nonzres=nonzres+1
           resact(nonzres)=rtype
        endif
        if (nbonds(rtype).gt.nbondgt) nbondgt=nbonds(rtype)
        if (nbonds3c(rtype).gt.nbondgt3c) nbondgt3c=nbonds3c(rtype)
      enddo
c
c-----Determine the possible bonds from the intergroup DIs
c
      if (ditwocent) then
      endif
c
c-----Alpha, Beta, and total number of valence electrons.
c
      nval =nel-mocore*2
      naval=mal-mocore
      nbval=mbe-mocore
c
c-----Obtain number of alpha, beta, and global probabilities.
c
      nproba=1
      nprobb=1
      mmm=1
      do i=1,ngroup-1
        nproba=nproba*(naval+i)
        nprobb=nprobb*(nbval+i)
        mmm=mmm*i
      enddo
      nproba=nproba/mmm
      nprobb=nprobb/mmm
      combi=one
      do i=0,nval-1
        combi=combi*dble(nval+ngroup-1-i)/dble(nval-i)
      enddo
      nprobt=int(combi+1D-10)
c
c-----Allocate space for probabilities, and alpha, beta, and global 
c     resonance structures. Remember that resnca(), resncb(), and 
c     resnc() do not include CORE electrons. To include the latter the 
c     following changes are mandatory:
c
c     resnca(n,ig) <--- resnca(n,ig) +     mocogrp(ig)
c     resncb(n,ig) <--- resncb(n,ig) +     mocogrp(ig)
c
      call allocate_space_for_rsrs (nproba,nprobb,ngroup,nbonds)
c
c-----Obtain alpha and beta resonance structures.
c
      if (.not.allocated(eleca)) then
         allocate (eleca(2,ngroup),stat=ier)
         if (ier.ne.0) stop 'optedf.f: Cannot allocate eleca()'
      endif
      if (.not.allocated(elecb)) then
        allocate (elecb(2,ngroup),stat=ier)
        if (ier.ne.0) stop 'optedf.f: Cannot allocate elecb()'
      endif
      do i=1,ngroup
        eleca(1,i)=izero
        elecb(1,i)=izero
        eleca(2,i)=naval
        elecb(2,i)=nbval
      enddo

      call rnprobs (eleca,resnca,naval,ngroup,nproba,lw)
      call rnprobs (elecb,resncb,nbval,ngroup,nprobb,lw)

      if (allocated(eleca)) then
         deallocate (eleca,stat=ier)
         if (ier.ne.0) stop 'optedf.f: Cannot deallocate eleca()'
      endif
      if (allocated(elecb)) then
        deallocate (elecb,stat=ier)
        if (ier.ne.0) stop 'optedf.f: Cannot deallocate elecb()'
      endif
c
c-----Obtain global resonance structures.
c
      n=0
      do ia=1,nproba
        do ib=1,nprobb
          n=n+1
          if (n.eq.1) then
            np=1
            abtot(n)=np
            do i=1,ngroup
              resnc(np,i)=resnca(ia,i)+resncb(ib,i)
            enddo
          else
            j=1
            inlist=.false.
            do while (j.le.np.and..not.inlist)
              inlist=.true.
              do k=1,ngroup
                iab=resnca(ia,k)+resncb(ib,k)
                inlist=inlist.and.(iab.eq.resnc(j,k))
              enddo
              if (inlist) abtot(n)=j
              j=j+1
            enddo
            if (.not.inlist) then
              np=np+1
              abtot(n)=np
              do i=1,ngroup
                resnc(np,i)=resnca(ia,i)+resncb(ib,i)
              enddo
            endif
          endif
        enddo
      enddo
c
c-----Prepare the input of the optimization routine.
c
      nvar=nonzres-1
      if (nvar.gt.0) then
        do i=1,nvar
          xvar(i)=1d0/dble(nonzres)
        enddo
      endif
c
c-----2-center bonds
c
      do m=1,nonzres
        rtype=resact(m)
        do i=1,nbonds(rtype)
          typ=ipair(i,3,rtype)
          thistyp=.false.
c
c---------Test whether this type of bond is present or not in the
c         current resonance structure.
c
          j=1
          do while (j.lt.i.and..not.thistyp)
            if (ipair(j,3,rtype).eq.typ) thistyp=.true.
            j=j+1
          enddo
          if (.not.thistyp) then
c
c.......... ¿ and in previous resonance structures ?
c
            prevtyp=.false.
            do l=1,m-1
              typ1=resact(l)
              do j=1,nbonds(typ1)
                if (ipair(j,3,typ1).eq.rtype) prevtyp=.true.
              enddo
            enddo
            if (.not.prevtyp) then
              if (.not.fixeq(typ)) then
                nvar=nvar+1
                if (nvar.gt.mvartot) then
                  stop 'optedf.f: Too many variables to optimize'
                endif
                xvar(nvar)=qqty(typ)
              endif
              if (.not.fixef(typ)) then
                nvar=nvar+1
                if (nvar.gt.mvartot) then
                  stop 'optedf.f: Too many variables to optimize'
                endif
                xvar(nvar)=ffty(typ)
              endif
            endif
          endif
        enddo
      enddo
c
c-----3-center bonds
c
      do m=1,nonzres
        rtype=resact(m)
        do i=1,nbonds3c(rtype)
          typ=itrio(i,4,rtype)
          thistyp=.false.
c
c---------Test whether this type of bond is present or not in the
c         current resonance structure.
c
          j=1
          do while (j.lt.i.and..not.thistyp)
            if (itrio(j,4,rtype).eq.typ) thistyp=.true.
            j=j+1
          enddo
          if (.not.thistyp) then
c
c.......... ¿ and in previous resonance structures ?
c
            prevtyp=.false.
            do l=1,m-1
              typ1=resact(l)
              do j=1,nbonds3c(typ1)
                if (itrio(j,4,typ1).eq.rtype) prevtyp=.true.
              enddo
            enddo
            if (.not.prevtyp) then
              do k=1,5
                if (.not.fixet(typ,k)) then
                  nvar=nvar+1
                  if (nvar.gt.mvartot) then
                    stop 'optedf.f: Too many variables to optimize'
                  endif
                  xvar(nvar)=prob3c(typ,k)
                endif
              enddo
            endif
          endif
        enddo
      enddo
c
c-----Number of protons in the molecule
c
      ntotel=0
      do i=1,ncent
        ntotel=ntotel+anint(charge(i))
      enddo
c
c-----NTOTREST must be = 0 for a neutral molecule.
c-----NTOTREST must be > 0 for a negatively charged molecule
c-----NTOTREST must be < 0 for a negatively charged molecule
c
      do m=1,nonzres
        rtype=resact(m)
        ntotrest=0
        do i=1,ngroup
          ntotrest=ntotrest+qrest(i,rtype)
        enddo
        if (nel.eq.ntotel) then  ! neutral molecule
          write (lw,199) 
          if (ntotrest.gt.0) write (lw,198) iabs(ntotrest)
          if (ntotrest.lt.0) write (lw,197) iabs(ntotrest)
          if (ntotrest.ne.0) stop
        elseif (nel.gt.ntotel) then  ! Negatively charged molecule
          write (lw,200) dble(nel-ntotel)
          if (abs(ntotrest).gt.nel-ntotel) then
            write (lw,198) nel-ntotel
            stop
          elseif (abs(ntotrest).lt.nel-ntotel) then
            write (lw,197) nel-ntotel
            stop
          endif
        elseif (nel.lt.ntotel) then  ! Positively charged molecule
          write (lw,201) dble(nel-ntotel)
          if (abs(ntotrest).lt.ntotel-nel) then
            write (lw,198) ntotel-nel
            stop
          elseif (abs(ntotrest).gt.ntotel-nel) then
            write (lw,197) ntotel-nel
            stop
          endif
        endif
      enddo
c
c-----Minimal number of electrons in each atom
c
      do m=1,nonzres
        rtype=resact(m)
        qmin(1:ngroup,rtype)=qrest(1:ngroup,rtype)
        do i=1,ngroup
          do j=1,nfugrp(i)
            qmin(i,rtype)=qmin(i,rtype)+anint(charge(ifugrp(j,i)))
          enddo
        enddo
      enddo
      do m=1,nonzres
        rtype=resact(m)
        do i=1,nbonds(rtype)
          i1=ipair(i,1,rtype)
          i2=ipair(i,2,rtype)
          qmin(i1,rtype)=qmin(i1,rtype)-1
          qmin(i2,rtype)=qmin(i2,rtype)-1
        enddo
      enddo
      do m=1,nonzres
        rtype=resact(m)
        do i=1,ngroup
          if (qmin(i,rtype).lt.0) then
            write (0,33) rtype,i,qmin(i,rtype)
            stop
          endif
        enddo
      enddo
 33   format (1x,'# In resonance structure ',I1,' the minimal',/,
     & 1x,'# number of electrons in fragment ',I2,' is negative (',
     & I2,'). Bad definition of bonding pattern.')
c
c-----Obtain exact probabilities.
c
      if (icorr) then
        call corredf (epsdet,nproba,nprobb,ngroup,ndets,nmo,
     &    ncore,nelact,nel,moval,ival,naval,nbval,ifilc,lw,wfnf,sg)
      else
        call hfedf (nproba,nprobb,ngroup,nmo,moval,ival,nel,
     &    naval,nbval,ifilc,sg)
      endif
c
c-----Write information corresponding to this RUN
c
      write (lw,10) weigedf,weigpop,weigdis,
     &    (nbonds(resact(m)),m=1,nonzres)
      write (lw,11)
 10   format (//1x,96('-'),/1x,
     & '# Fitting the EDF as a direct product of BONDING (2c,2ne) EDFs'
     & /1x,96('-'),/,1x,'Objective function WEIGHTS:',
     & ' (EDF,ElecPopul,DelocInd) = ',3F12.6,/,1x,
     & 'ALLEGED CHEMICAL BONDS IN ALL RESONANCE STRUCTURES   = ',20I4)
 11   format (1x,'INPUT PARAMETERS OF THE FITTING PROCEDURE.',
     & 1x,'An asterisk "(*)" means a fixed parameter')
      do m=1,nonzres
        rtype=resact(m)
        write (lw,29) rtype
        do i=1,nbonds(rtype)
          typ=ipair(i,3,rtype)
          xpi=(qqty(typ)+1D0)*half
          del=4D0*(1-ffty(typ))*xpi*(one-xpi)
          ip1=ipair(i,1,rtype)
          ip2=ipair(i,2,rtype)
          choq='   '
          if (fixeq(typ)) choq='(*)'
          chof='   '
          if (fixef(typ)) chof='(*)'
          chot='   '
          if (fixeq(typ).and.fixef(typ)) chot='(*)'
          ffis=ffty(typ)
          qqis=qqty(typ)
          write (lw,30) ip1,ip2,qqis,choq,xpi,choq,ffis,chof,del,chot
        enddo
      enddo
 29   format (1x,96('-'),/,1x,'RESONANCE STRUCTURE ',I1,/,1x,96('-'),/,
     &  9x,'group 1',1x,'group 2',5x,'Q',15x,'PI',16x,'f',15x,'Delta',
     & /1x,96('-'))
 30   format (1x,'BOND',I7,1x,I7,6(1x,F12.6,1x,a))
c
c-----Call to the optimization routine (praxis)
c
      machep=2.220446049250313D-16
      compedf=.true.
      if (nvar.gt.0) then
        edfdif = praxis (t0,machep,h0,nvar,prin,xvar,gedfit,fmin)
      endif
      largwr=.true.
      write (lw,40)
 40   format (/1x,96('-'),/,29x,'FINAL OPTIMIZED EDF',/1x,96('-'))
      edfitlast = gedfit (xvar,nvar)

      call deallocate_space_for_rsrs ()
      call dspace_for_bonds ()
      call timer (4,iopedf,'_optedf   ',-1)
      return

 199  format (1x,'# NEUTRAL MOLECULE')
 200  format (1x,'# NEGATIVELY CHARGED MOLECULE, CHARGE = ',F12.6)
 201  format (1x,'# POSITIVELY CHARGED MOLECULE, CHARGE = ',F12.6)
 197  format (1x,'# ADD ',I2,' PLUS  orders in the input file')
 198  format (1x,'# ADD ',I2,' MINUS orders in the input file')
      end
