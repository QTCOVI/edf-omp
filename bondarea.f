c
c-----------------------------------------------------------------------
c
      module space_for_bonds
      save
      integer, allocatable,dimension (:,:,:)   :: ipair,itrio
      integer,allocatable,dimension  (:,:)     :: qmin
      integer, allocatable,dimension (:,:)     :: itype
      integer,allocatable,dimension  (:)       :: nbonds,nbonds3c,resact
      integer,allocatable,dimension  (:,:)     :: qrest
      integer,allocatable,dimension  (:)       :: numprob
      real(kind=8), allocatable,dimension  (:)       :: qqty,ffty
      logical,allocatable,dimension  (:)       :: fixeq,fixef
      logical,allocatable,dimension  (:,:)     :: fixet
      real(kind=8),allocatable,dimension (:,:) :: prob3c
      integer  nfrag,lwrite,nbondgt,nbondgt3c,maxres,nonzres
      real(kind=8) weigedf,weigpop,weigdis
      logical  compedf,largwr
      end module space_for_bonds
c
c-----------------------------------------------------------------------
c
      subroutine aspace_for_bonds (maxbond,mvarh,mvar,mvar3c,maxbond3c)
      USE space_for_bonds
c
c-----Maximun number of resonance structures
c
      maxres = 4
      if (.not.allocated(qqty)) then
        allocate (qqty(mvarh),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate qqty()'
      endif
      if (.not.allocated(ffty)) then
        allocate (ffty(mvarh),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate ffty()'
      endif
      if (.not.allocated(fixeq)) then
        allocate (fixeq(mvarh),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate fixeq()'
      endif
      if (.not.allocated(fixef)) then
        allocate (fixef(mvarh),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate fixef()'
      endif
      if (.not.allocated(ipair)) then
        allocate (ipair(maxbond,3,maxres),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate ipair()'
      endif
      if (.not.allocated(itrio)) then
        allocate (itrio(maxbond3c,4,maxres),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate ipair()'
      endif
      if (.not.allocated(itype)) then
        allocate (itype(mvarh,2),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate itype()'
      endif
      if (.not.allocated(qmin)) then
        allocate (qmin(nfrag,maxres),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate qmin()'
      endif
      if (.not.allocated(qrest)) then
        allocate (qrest(nfrag,maxres),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate qrest()'
      endif
      if (.not.allocated(nbonds)) then
        allocate (nbonds(maxres),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate nbonds()'
      endif
      if (.not.allocated(nbonds3c)) then
        allocate (nbonds3c(maxres),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate nbonds3c()'
      endif
      if (.not.allocated(resact)) then
        allocate (resact(maxres),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate resact()'
      endif
      if (.not.allocated(numprob)) then
        allocate (numprob(maxres),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate numprob()'
      endif
      if (.not.allocated(fixet)) then
        allocate (fixet(maxbond3c,5),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate fixet()'
      endif
      if (.not.allocated(prob3c)) then
        allocate (prob3c(maxbond3c,5),stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate prob3c()'
      endif
      end subroutine aspace_for_bonds
c
c-----------------------------------------------------------------------
c
      subroutine dspace_for_bonds ()
      USE space_for_bonds
      if (allocated(qqty)) then
        deallocate (qqty,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate qqty()'
      endif
      if (allocated(ffty)) then
        deallocate (ffty,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate ffty()'
      endif
      if (allocated(fixeq)) then
        deallocate (fixeq,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate fixeq()'
      endif
      if (allocated(fixef)) then
        deallocate (fixef,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate fixef()'
      endif
      if (allocated(ipair)) then
        deallocate (ipair,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate ipair()'
      endif
      if (allocated(itrio)) then
        deallocate (itrio,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot allocate itrio()'
      endif
      if (allocated(itype)) then
        deallocate (itype,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate itype()'
      endif
      if (allocated(qmin)) then
        deallocate (qmin,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate qmin()'
      endif
      if (allocated(qrest)) then
        deallocate (qrest,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate qrest()'
      endif
      if (allocated(nbonds)) then
        deallocate (nbonds,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate nbonds()'
      endif
      if (allocated(nbonds3c)) then
        deallocate (nbonds3c,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate nbonds3c()'
      endif
      if (allocated(resact)) then
        deallocate (resact,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate resact()'
      endif
      if (allocated(numprob)) then
        deallocate (numprob,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate numprob()'
      endif
      if (allocated(fixet)) then
        deallocate (fixet,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate fixet()'
      endif
      if (allocated(prob3c)) then
        deallocate (prob3c,stat=ier)
        if (ier.ne.0) stop 'bondarea.f: Cannot deallocate prob3c()'
      endif
      end subroutine dspace_for_bonds
c
c-----------------------------------------------------------------------
c
