c
c-----Obtain approximate EDF as a direct product of (2c,2ne) bond EDFs.
c
      real(kind=8) function gedfit (xvar,nvar)
      USE space_for_rsrs
      USE space_for_bonds
      include 'implicit.inc'
      include 'corr.inc'
      include 'param.inc'
      include 'constants.inc'
      integer,allocatable,dimension (:)     :: jloop,iloop,floop,sloop
      integer,allocatable,dimension (:)     :: mloop
      integer,allocatable,dimension (:)     :: jsum
      integer,allocatable,dimension (:,:)   :: eleca,elecb
      real(kind=8), allocatable,dimension (:)     :: qq,xp,ff
      real(kind=8), allocatable,dimension (:,:,:) :: p
      real(kind=8), allocatable,dimension (:,:)   :: prob
      real(kind=8), allocatable,dimension (:)     :: sumprob
      real(kind=8), allocatable,dimension (:)     :: popul,xliex,p1ex
      real(kind=8), allocatable,dimension (:)     :: xliap,p1ap,xliexa
      real(kind=8), allocatable,dimension (:)     :: p1exa
      real(kind=8), allocatable,dimension (:,:)   :: xliapf,p1apf
      real(kind=8), allocatable,dimension (:,:)   :: delind,sumdi
      real(kind=8) xvar(nvar)
      real(kind=8) resw(maxres),npap(maxres),npex(maxres)
      logical isinc,nodoedf,somex,thistyp,prevtyp
      integer rtype,typ,typ1
c
c-----Allocate memory for arrays used in the calculation.
c
      if (.not.allocated(jloop)) then
        allocate (jloop(nbondgt),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate jloop()'
      endif
      if (.not.allocated(iloop)) then
        allocate (iloop(nbondgt),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate iloop()'
      endif
      if (.not.allocated(floop)) then
        allocate (floop(nbondgt),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate floop()'
      endif
      if (.not.allocated(sloop)) then
        allocate (sloop(nbondgt),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate sloop()'
      endif
      if (.not.allocated(mloop)) then
        allocate (mloop(nbondgt),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate mloop()'
      endif
      if (.not.allocated(jsum)) then
        allocate (jsum(nfrag),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate jsum()'
      endif
      if (.not.allocated(eleca)) then
        allocate (eleca(2,nfrag),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate eleca()'
      endif
      if (.not.allocated(elecb)) then
        allocate (elecb(2,nfrag),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate elecb()'
      endif
      if (.not.allocated(qq)) then
        allocate (qq(nbondgt),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate qq()'
      endif
      if (.not.allocated(xp)) then
        allocate (xp(nbondgt),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate xp()'
      endif
      if (.not.allocated(ff)) then
        allocate (ff(nbondgt),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate ff()'
      endif
      if (.not.allocated(p)) then
        allocate (p(nbondgt,0:2,0:2),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate p()'
      endif
      if (.not.allocated(prob)) then
        allocate (prob(nprobt,nonzres),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate prob()'
      endif
      if (.not.allocated(sumprob)) then
        allocate (sumprob(nprobt),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate sumprob()'
      endif
      if (.not.allocated(xliex)) then
        allocate (xliex(nfrag),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate xliex()'
      endif
      if (.not.allocated(p1ex)) then
        allocate (p1ex(nfrag),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate p1ex()'
      endif
      if (.not.allocated(xliexa)) then
        allocate (xliexa(nfrag),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate xliexa()'
      endif
      if (.not.allocated(p1exa)) then
        allocate (p1exa(nfrag),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate p1exa()'
      endif
      if (.not.allocated(xliap)) then
        allocate (xliap(nfrag),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate xliap()'
      endif
      if (.not.allocated(p1ap)) then
        allocate (p1ap(nfrag),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate p1ap()'
      endif
      if (.not.allocated(xliapf)) then
        allocate (xliapf(nfrag,nonzres),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate xliapf()'
      endif
      if (.not.allocated(p1apf)) then
        allocate (p1apf(nfrag,nonzres),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate p1apf()'
      endif
      if (.not.allocated(delind)) then
        allocate (delind(nfrag,nfrag),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate delind()'
      endif
      if (.not.allocated(sumdi)) then
        allocate (sumdi(nfrag,nfrag),stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot allocate sumdi()'
      endif

      lw=lwrite
c
c-----Recover chemical bonding parameters from the current values
c     of the optimized variables.
c
      sumres=0d0
      if (nonzres.gt.1) then
        do i=1,nonzres-1
          resw(i)=xvar(i)
          sumres=sumres+resw(i)
        enddo
      endif
      resw(nonzres)=1d0-sumres
c
c-----Recover bonding parameters from the optimizable variables.
c
      ncur=nonzres-1
      do m=1,nonzres
        rtype=resact(m)
        do i=1,nbonds(rtype)
          typ=ipair(i,3,rtype)
          thistyp=.false.
c
c---------Test whether thys type of bond is present or not in the
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
                ncur=ncur+1
                qqty(typ)=xvar(ncur)
              endif
              if (.not.fixef(typ)) then
                ncur=ncur+1
                ffty(typ)=xvar(ncur)
              endif
            endif
          endif
        enddo
      enddo
      if (ncur.ne.nvar) then
        stop 'gedfit.f: Bad number of optimized variables'
      endif
c
c-----LOOP over all the resonance structures.
c
      idxp(1:nprobt,1:nonzres)=.false.
      do m=1,nonzres
        rtype=resact(m)
        do i=1,nbonds(rtype)
          typ=ipair(i,3,rtype)
          qq(i)=qqty(typ)
          ff(i)=ffty(typ)
          xp(i)=(qq(i)+1D0)/2D0
          oxp=1D0-xp(i)
          ffp=ff(i)+1D0
          p(i,1,1)=ffp*2D0*xp(i)*oxp
          p(i,2,0)=xp(i)*(1D0-ffp*oxp)
          p(i,0,2)=oxp*(1D0-ffp*xp(i))
          iloop(i)=0
          floop(i)=2
          sloop(i)=1
          jloop(i)=iloop(i)
        enddo
c
c-------Performs the core of the algorithm using the Nested Summation
c       Symbol (NSS) by Emili Besalu.
c
        k=nbonds(rtype)
        prob(1:nprobt,m)=0d0
        do while (k.gt.0)
          if ((jloop(k)-floop(k))*sloop(k).gt.0) then
             jloop(k)=iloop(k)
             k=k-1
          else
            jsum(1:nfrag)=qmin(1:nfrag,rtype)
            do i=1,nbonds(rtype)
              mloop(i)=2-jloop(i)
              i1=ipair(i,1,rtype)
              i2=ipair(i,2,rtype)
              jsum(i1)=jsum(i1)+jloop(i)
              jsum(i2)=jsum(i2)+mloop(i)
            enddo
            probal=1D0
            do i=1,nbonds(rtype)
              probal=probal*p(i,jloop(i),mloop(i))
            enddo
            isinc=.false.
            i=0
            do while (i.lt.nprobt.and.(.not.isinc))
              i=i+1
              isinc=.true.
              do j=1,nfrag
                isinc=isinc.and.(jsum(j).eq.resnc(i,j))
              enddo
              if (isinc) prob(i,m)=prob(i,m)+probal
              if (isinc) idxp(i,m)=.true.
            enddo
            k=nbonds(rtype)
          end if
          if (k.gt.0) jloop(k)=jloop(k)+sloop(k)
        end do
      enddo
c
c-----Objective function (Part of the probabilities)
c
      fun1=0d0
      do i=1,nprobt
        pthisp=pexact(i)
        sumprob(i)=0d0
        do m=1,nonzres
          probweig=prob(i,m)*resw(m)
          pthisp=pthisp-probweig
          sumprob(i)=sumprob(i)+probweig
        enddo
        fun1=fun1+pthisp**2
      enddo
c
c-----Average electron populations
c
      p1apf (1:nfrag,1:nonzres)=0d0
      do i=1,nprobt
        somex=.false.
        do m=1,nonzres
          somex=somex.or.idxp(i,m)
        enddo
        if (somex) then
          do j=1,nfrag
            irsrs=resnc(i,j)
            do m=1,nonzres
              p1apf(j,m)=p1apf(j,m)+prob(i,m)*irsrs
            enddo
          enddo
        endif
      enddo
      p1ap(1:nfrag)=0d0
      do j=1,nfrag
        do m=1,nonzres
          p1ap(j)=p1ap(j)+resw(m)*p1apf(j,m)
        enddo
      enddo
c
c-----Delocalization indices
c
      sumdi=0d0
      ncurrent=nonzres-1
      do m=1,nonzres
        rtype=resact(m)
        delind(1:nfrag,1:nfrag)=0d0
        do i=1,nbonds(rtype)
          typ=ipair(i,3,rtype)
          i1=ipair(i,1,rtype)
          i2=ipair(i,2,rtype)
          imin=min(i1,i2)
          imax=max(i1,i2)
          qq(i)=qqty(typ)
          ff(i)=ffty(typ)
          xp(i)=(qq(i)+1D0)/2D0
          oxp=1D0-xp(i)
          ffp=ff(i)+1D0
          del=4D0*(1-ff(i))*xp(i)*(one-xp(i))
          delind(imin,imax)=delind(imin,imax)+del
          delind(imax,imin)=delind(imin,imax)
        enddo
        sumdi=sumdi+delind*resw(m)
      enddo
c
c-----Objective function (part of electron populations + DIs)
c
      fun2=0d0
      do j=1,nfrag
        fun2=fun2+(p1ap(j)-p1(j))**2
      enddo    
      fun3=0d0
      do j=2,nfrag
        do k=1,j-1
          fun3 = fun3 + (sumdi(j,k)-diexact(j,k))**2
        enddo
      enddo
c
c-----Final value of the objective function
c
      gedfit = fun1 * weigedf + fun2 * weigpop + fun3 * weigdis
c
c-----After convergence ... DO THIS IF
c
      if (largwr) then
        write (lw,30) 
        sumex=0d0
        xliex (1:nfrag)=0d0
        p1ex  (1:nfrag)=0d0
        xliexa(1:nfrag)=0d0
        p1exa (1:nfrag)=0d0
        xliapf(1:nfrag,1:nonzres)=0d0
        p1apf (1:nfrag,1:nonzres)=0d0
        do i=1,nprobt
          pex=pexact(i)
          do j=1,nfrag
            irsrs=resnc(i,j)
            p1exa(j) =p1exa(j) +pex*irsrs
            xliexa(j)=xliexa(j)+pex*irsrs**2
          enddo
          somex=.false.
          do m=1,nonzres
            somex=somex.or.idxp(i,m)
          enddo
          if (somex) then
            sumex=sumex+pexact(i)
            do j=1,nfrag
              irsrs=resnc(i,j)
              p1ex(j)=p1ex(j)+pex*irsrs
              xliex(j)=xliex(j)+pex*irsrs**2
              do m=1,nonzres
                p1apf(j,m)=p1apf(j,m)+prob(i,m)*irsrs
                xliapf(j,m)=xliapf(j,m)+prob(i,m)*irsrs**2
              enddo
            enddo
            write (lw,42,advance='no') pex,sumprob(i)
            write (lw,41) (resnc(i,j),j=1,nfrag)
          endif
        enddo
        do j=1,nfrag
          xliex(j)=-(xliex(j)-p1ex(j)*(p1ex(j)+1d0))
          xliexa(j)=-(xliexa(j)-p1exa(j)*(p1exa(j)+1d0))
        enddo
        write (lw,20) sumex
        write (lw,21) gedfit
        write (lw,50)
        sumdi=0d0
        ncurrent=nonzres-1
        do m=1,nonzres
          rtype=resact(m)
          write (lw,71) rtype,resw(m)*100d0
          rtype=resact(m)
          delind(1:nfrag,1:nfrag)=0d0
          do i=1,nbonds(rtype)
            typ   = ipair(i,3,rtype)
            i1    = ipair(i,1,rtype)
            i2    = ipair(i,2,rtype)
            imin  = min(i1,i2)
            imax  = max(i1,i2)
            qq(i) = qqty(typ)
            ff(i) = ffty(typ)
            xp(i) = (qq(i)+1D0)/2D0
            oxp   = 1D0-xp(i)
            ffp   = ff(i)+1D0
            del   = 4D0*(1-ff(i))*xp(i)*(one-xp(i))
            delind(imin,imax) = delind(imin,imax)+del
            delind(imax,imin) = delind(imin,imax)
            p11   = ffp*2D0*xp(i)*oxp
            p20   = xp(i)*(1D0-ffp*oxp)
            p02   = oxp*(1D0-ffp*xp(i))
            ip1   = ipair(i,1,rtype)
            ip2   = ipair(i,2,rtype)
            write (lw,31) ip1,ip2,qq(i),xp(i),ff(i),del,p02,p11,p20
          enddo
          write (lw,*)
          write (lw,51)
*         write (lw,51,advance='no')
*         write (lw,510)
          do j=1,nfrag
            xliapf(j,m)=-(xliapf(j,m)-p1apf(j,m)*(p1apf(j,m)+1d0))
            perex=100d0*xliex(j)/p1ex(j)
            perap=100d0*xliapf(j,m)/p1apf(j,m)
            xap=xliapf(j,m)
            xex=xliex(j)
*           write (lw,54) j,p1ex(j),p1apf(j,m),xex,xap,perex,perap
            write (lw,54) j,p1apf(j,m),xap,perap
          enddo

          sumdi=sumdi+delind*resw(m)
          write (lw,*)
          write (lw,52)
          do i=2,nfrag
            do j=1,i-1
              write (lw,53) i,j,delind(i,j)
            enddo
          enddo
        enddo
c
c------Total approximate averaged populations 
c
        if (nonzres.ne.1) then
          write (lw,63)
          p1ap(1:nfrag)=0d0
          xliap(1:nfrag)=0d0
          do j=1,nfrag
            do m=1,nonzres
              p1ap(j)=p1ap(j)+resw(m)*p1apf(j,m)
              xliap(j)=xliap(j)+resw(m)*xliapf(j,m)
            enddo
            papav=100d0*xliap(j)/p1ap(j)
            write (lw,54) j,p1ap(j),xliap(j),papav
          enddo
          write (lw,*)
          write (lw,52)
          do i=2,nfrag
            do j=1,i-1
              write (lw,53) i,j,sumdi(i,j)
            enddo
          enddo
        endif
c
c-------Exact populations and LIs including all RSRSs
c
        write (lw,*)
        write (lw,610)
        do j=1,nfrag
          perex=100d0*xliexa(j)/p1exa(j)
          write (lw,54) j,p1exa(j),xliexa(j),perex
        enddo
      endif

      if (allocated(jloop)) then
        deallocate (jloop,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate jloop()'
      endif
      if (allocated(iloop)) then
        deallocate (iloop,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate iloop()'
      endif
      if (allocated(floop)) then
        deallocate (floop,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate floop()'
      endif
      if (allocated(sloop)) then
        deallocate (sloop,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate sloop()'
      endif
      if (allocated(mloop)) then
        deallocate (mloop,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate mloop()'
      endif
      if (allocated(jsum)) then
        deallocate (jsum,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate jsum()'
      endif
      if (allocated(eleca)) then
        deallocate (eleca,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate eleca()'
      endif
      if (allocated(elecb)) then
        deallocate (elecb,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate elecb()'
      endif
      if (allocated(qq)) then
        deallocate (qq,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate qq()'
      endif
      if (allocated(xp)) then
        deallocate (xp,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate xp()'
      endif
      if (allocated(ff)) then
        deallocate (ff,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate ff()'
      endif
      if (allocated(p)) then
        deallocate (p,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate p()'
      endif
      if (allocated(prob)) then
        deallocate (prob,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate prob()'
      endif
      if (allocated(sumprob)) then
        deallocate (sumprob,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate sumprob()'
      endif
      if (allocated(xliex)) then
        deallocate (xliex,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate xliex()'
      endif
      if (allocated(p1ex)) then
        deallocate (p1ex,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate p1ex()'
      endif
      if (allocated(xliexa)) then
        deallocate (xliexa,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate xliexa()'
      endif
      if (allocated(p1exa)) then
        deallocate (p1exa,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate p1exa()'
      endif
      if (allocated(xliap)) then
        deallocate (xliap,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate xliap()'
      endif
      if (allocated(p1ap)) then
        deallocate (p1ap,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate p1ap()'
      endif
      if (allocated(xliapf)) then
        deallocate (xliapf,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate xliapf()'
      endif
      if (allocated(p1apf)) then
        deallocate (p1apf,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate p1apf()'
      endif
      if (allocated(delind)) then
        deallocate (delind,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate delind()'
      endif
      if (allocated(sumdi)) then
        deallocate (sumdi,stat=ier)
        if (ier.ne.0) stop ' # gedfit.f: Cannot deallocate sumdi()'
      endif

      return

 30   format (13x,'EXACT AND APPROXIMATE EDFs',/1x,96('-'),/,
     & 12x,'EXACT',12x,'APPROXIMATE',6x,'RSRS')
 41   format (5x,20I3)
 42   format (7x,20(1x,F16.10))
 20   format (1x,'SUM  = ',F16.10)
 21   format (1x,'OBJETIVE FUNCTION = ',E16.10)
 50   format (/1x,96('-'),/21x,
     & 'PARAMETERS OF THE FITTING PROCEDURE',/,1x,96('-'))
 71   format (//1x,96('-'),/,1x,'RESONANCE STRUCTURE ',I1,
     & 6x,'(WEIGHT = ',F5.1,' % )',/,1x,'Q(1) = 2 PI - 1',/,
     & 1x,'Delta = 4 (1 - f) PI (1 - PI) = (1-f) (1-Q(1)^2)',/,
     &  1x,96('-'),/7x,'group 1',1x,'group 2',2x,'Q(1)',9x,'PI',
     & 10x,'f',8x,'Delta',6x,'p(0,2)',5x,'p(1,1)',5x,'p(2,0)',/,
     & 1x,96('-'))
*51   format (1x,'AVERAGE POPULATIONS AND LIs. "pexac" ARE ')
*510  format ('EXACT VALUES INCLUDING ONLY THE BONDING RSRSs',/,
*    & 1x,86('-'),/22x,'<N>_pexac',3x,'<N>_appr',3x,'LI_pexac',4x,
*    & 'LI_appr',2x,'%LOC pexac',2x,'%LOC appr',/,22x,65('-'))
 51   format (' AVERAGE POPULATIONS AND LIs',/,
     & ' ----------------------------------------------------',/,
     & 23x,'<N>_appr',3x,'LI_appr',2x,'%LOC appr',/,1x,52('-'))





 54   format (1x,'Fragment ',I2,8x,8(1x,F10.6))
 610  format (//
     & ' EXACT AVERAGE POPULATIONS AND LIs INCLUDING ALL RSRSs',/,
     & ' -----------------------------------------------------',/,
     & 22x,'<N>_exact',3x,'LI_exact',2x,'%LOC exact',/22x,32('-'))

 53   format (6x,'delta_(',2I3,' )     =   ',F12.6)
 52   format (1x,'APPROXIMATE INTERGROUP DELOCALIZATION INDICES',/,
     &        1x,'---------------------------------------------')
 31   format (1x,'BOND',I5,1x,I7,2x,7(1x,F10.6))
 63   format (//
     &  ' AVERAGED RESULTS COUNTING ALL THE RESONANCE STRUCTURE',/,
     &  ' -----------------------------------------------------',/,
     &  23x,'<N>_appr  LI_appr    %LOC appr',/,23x,30('-'))
      end
