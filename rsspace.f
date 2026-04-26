c
c-----------------------------------------------------------------------
c
      module space_for_rsrs
      save
      real(kind=8), allocatable,dimension (:)     :: pexact,p1
      real(kind=8), allocatable,dimension (:,:)   :: diexact
      integer, allocatable,dimension (:,:)   :: resnc,resnca,resncb
      integer, allocatable,dimension (:,:,:) :: resap
      integer, allocatable,dimension (:)     :: abtot
      integer, allocatable,dimension (:,:)   :: edfinc,apinc
      logical, allocatable,dimension (:,:)   :: idxp
      integer nprobt,mres
      end module space_for_rsrs
c
c-----------------------------------------------------------------------
c
      subroutine allocate_space_for_rsrs (nproba,nprobb,ngroup,nbonds)
      USE space_for_rsrs
      if (.not.allocated(pexact)) then
        allocate (pexact(nprobt),stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot allocate pexact()'
      endif
      if (.not.allocated(p1)) then
        allocate (p1(ngroup),stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot allocate p1()'
      endif
      if (.not.allocated(diexact)) then
        allocate (diexact(ngroup,ngroup),stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot allocate diexact()'
      endif
      if (.not.allocated(edfinc)) then
        allocate (edfinc(nprobt,mres),stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot allocate edfinc()'
      endif
      if (.not.allocated(apinc)) then
        allocate (apinc(3**nbonds,mres),stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot allocate apinc()'
      endif
      if (.not.allocated(resnc)) then
        allocate (resnc(nprobt,ngroup),stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot allocate resnc()'
      endif
      if (.not.allocated(resap)) then
        allocate (resap(nprobt,ngroup,mres),stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot allocate resap()'
      endif
      if (.not.allocated(resnca)) then
        allocate (resnca(nproba,ngroup),stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot allocate resnca()'
      endif
      if (.not.allocated(resncb)) then
        allocate (resncb(nprobb,ngroup),stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot allocate resncb()'
      endif
      if (.not.allocated(abtot)) then
        allocate (abtot(nproba*nprobb),stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot allocate abtot()'
      endif
      if (.not.allocated(idxp)) then
        allocate (idxp(nprobt,mres),stat=ier)
        if (ier.ne.0) stop ' # rsspace.f: Cannot allocate idxp()'
      endif
      end subroutine allocate_space_for_rsrs
c
c-----------------------------------------------------------------------
c
      subroutine deallocate_space_for_rsrs ()
      USE space_for_rsrs
      if (allocated(pexact)) then
        deallocate (pexact,stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot deallocate pexact()'
      endif
      if (allocated(p1)) then
        deallocate (p1,stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot deallocate p1()'
      endif
      if (allocated(diexact)) then
        deallocate (diexact,stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot deallocate diexact()'
      endif
      if (allocated(edfinc)) then
        deallocate (edfinc,stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot deallocate edfinc()'
      endif
      if (allocated(apinc)) then
        deallocate (apinc,stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot deallocate apinc()'
      endif
      if (allocated(resnc)) then
        deallocate (resnc,stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot deallocate resnc()'
      endif
      if (allocated(resap)) then
        deallocate (resap,stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot deallocate resap()'
      endif
      if (allocated(resnca)) then
        deallocate (resnca,stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot deallocate resnca()'
      endif
      if (allocated(resncb)) then
        deallocate (resncb,stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot deallocate resncb()'
      endif
      if (allocated(abtot)) then
        deallocate (abtot,stat=ier)
        if (ier.ne.0) stop 'rsspace.f: Cannot deallocate abtot()'
      endif
      if (allocated(idxp)) then
        deallocate (idxp,stat=ier)
        if (ier.ne.0) stop ' # rsspace.f: Cannot deallocate idxp()'
      endif
      end subroutine deallocate_space_for_rsrs
