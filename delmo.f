c
c-----------------------------------------------------------------------
c
      subroutine delmo (aom,mal,mbe,line,lw)
c 
c.......................................................................
c
      USE      space_for_wfnbasis
      USE      space_for_wfncoef
      include  'implicit.inc'
      include  'param.inc'
      include  'constants.inc'
      include  'wfn.inc'
      include  'corr.inc'
      include  'stderr.inc'

      real(kind=8), allocatable,dimension (:,:,:)  :: aomx
      real(kind=8), allocatable,dimension (:,:)    :: cc
      integer,      allocatable,dimension (:)      :: isup,isupx,ordsup

      character*(*) line
      real(kind=8)       aom(ncent,nmo,nmo)
      logical       setint,ok,todelete,notinlist
      integer       lp,nsup,isupress,nosup
c
      if (ncore.le.0) return
c
c.....Determine the core MOs to be elliminated.
c
      if (.not.allocated(isup)) then
        allocate (isup(nmo),stat=ier)
        if (ier.ne.0) stop 'delmo.f: Cannot allocate isup()'
      endif
      if (.not.allocated(ordsup)) then
        allocate (ordsup(nmo),stat=ier)
        if (ier.ne.0) stop 'delmo.f: Cannot allocate ordsup()'
      endif
      if (.not.allocated(isupx)) then
        allocate (isupx(nmo),stat=ier)
        if (ier.ne.0) stop 'delmo.f: Cannot allocate isupx()'
      endif

      lp=1
      nsup=0
      do while (setint(n,line,lp))
        if (n.le.ncore.and.n.gt.0) then
          if (nsup.eq.0) then
            nsup = 1
            isup(1) = n
          else
            i=1
            notinlist=.true.
            do while (i.le.nsup.and.notinlist) 
              if (n.eq.isup(i)) notinlist=.false.
              i=i+1
            enddo
            if (notinlist) then
              nsup=nsup+1
              isup(nsup)=n
            endif
          endif
        endif
      enddo
c
c.....Determine the core MOs that ARE NOT elliminated.
c
      if (nsup.gt.0) then
        nosup=0
        do i=1,nmo
          todelete=.false.
          do k=1,nsup
            if (i.eq.isup(k)) todelete=.true.
          enddo
          if (.not.todelete) then
            nosup=nosup+1
            isupx(nosup)=i
          endif
        enddo
      else
        nosup=nmo
        do i=1,nmo
          isupx(i)=i
        enddo
      endif
c
c.....Reconstruct arrays of the wavefunction.
c
      if (.not.allocated(cc)) then
        allocate (cc(nmo+nmo,nprims),stat=ier)
        if (ier.ne.0) stop 'delmo.f: Cannot allocate cc()'
      endif
      cc(1:nmo+nmo,1:nprims)=coef(1:nmo+nmo,1:nprims)
      if (allocated(coef)) then
        deallocate (coef,stat=ier)
        if (ier.ne.0) stop 'delmo.f: Cannot deallocate coef()'
      endif
      nmold=nmo
      nmo=nmo-nsup
      ncore=ncore-nsup
      nel=nelact+2*ncore
      mal=mal-nsup
      mbe=mbe-nsup
      nalpha=nalpha-nsup
      nbeta=nbeta-nsup
      if (.not.allocated(coef)) then
        allocate (coef(1:nmo+nmo,nprims),stat=ier)
        if (ier.ne.0) stop 'delmo.f: Cannot allocate coef()'
      endif
      do i=1,nmo
        coef(i,    1:nprims)=cc(isupx(i),      1:nprims)
        coef(i+nmo,1:nprims)=cc(isupx(i)+nmold,1:nprims)
      enddo

      if (.not.allocated(aomx)) then
        allocate (aomx(1:ncent,nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'delmo.f: Cannot allocate aomx()'
      endif

      do i=1,nmo
        do j=1,nmo
          aomx(1:ncent,i,j)=aom(1:ncent,isupx(i),isupx(j))
        enddo
      enddo
      aom(1:ncent,1:nmo,1:nmo) = aomx(1:ncent,1:nmo,1:nmo)
c
      if (allocated(cc)) then
        deallocate (cc,stat=ier)
        if (ier.ne.0) stop 'delmo.f: Cannot allocate cc()'
      endif
      if (allocated(aomx)) then
        deallocate (aomx,stat=ier)
        if (ier.ne.0) stop 'delmo.f: Cannot allocate aomx()'
      endif

      write (lw,313) (isup(i),i=1,nsup)

      if (allocated(isup)) then
        deallocate (isup,stat=ier)
        if (ier.ne.0) stop 'delmo.f: Cannot deallocate isup()'
      endif
      if (allocated(isupx)) then
        deallocate (isupx,stat=ier)
        if (ier.ne.0) stop 'delmo.f: Cannot deallocate isupx()'
      endif
      if (allocated(ordsup)) then
        deallocate (ordsup,stat=ier)
        if (ier.ne.0) stop 'delmo.f: Cannot deallocate ordsup()'
      endif
      return
 313  format (1x,'# CORE MOs elliminated from the WFN: ',100I4)
      end
