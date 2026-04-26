      module space_for_sym
      save
      integer, allocatable, dimension (:)     :: ineq,mult
      integer, allocatable, dimension (:,:)   :: idx,indp,ipair,table
      real(kind=8), allocatable, dimension (:,:,:) :: chartab
      integer, allocatable, dimension (:,:)   :: irrep
      integer, allocatable, dimension (:,:)   :: idict,jpair
      end module space_for_sym
c-----------------------------------------------------------------------
      subroutine allocate_space_for_sym (nmo,ncent)
      USE        space_for_sym
c
c.....I think that this is true....:
c.....ineq(i): wfn index of non-eq atom i
c.....idx(i,1): non-eq index of wfn atom i
c.....idx(i,2): symop index taking the non-eq atom into wfn atom i
c.....idict(i,1): non-eq index of wfn atom i (up to nmono non-eq)
c.....idict(i,2): wfn index of non-eq atom i (up to nmono non-eq)
c
      allocate (chartab(nmo,6,48),stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot allocate chartab()'
      allocate (mult(ncent),stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot allocate mult()'
      allocate (idx(ncent,2),stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot allocate idx()'
      allocate (indp(ncent,ncent),stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot allocate indp()'
      allocate (ipair(ncent*ncent,2),stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot allocate ipair()'
      allocate (table(ncent,48),stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot allocate table()'
      allocate (irrep(nmo,0:6),stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot allocate irrep()'
      allocate (idict(ncent,2),stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot allocate idict()'
      allocate (jpair(ncent,0:ncent),stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot allocate jpair()'
      allocate (ineq(ncent),stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot allocate ineq()'
      end subroutine allocate_space_for_sym
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_sym
      USE        space_for_sym
      deallocate (chartab,stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot deallocate chartab()'
      deallocate (mult,stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot deallocate mult()'
      deallocate (idx,stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot deallocate idx()'
      deallocate (indp,stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot deallocate indp()'
      deallocate (ipair,stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot deallocate ipair()'
      deallocate (table,stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot deallocate table()'
      deallocate (irrep,stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot deallocate irrep()'
      deallocate (idict,stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot deallocate idict()'
      deallocate (jpair,stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot deallocate jpair()'
      deallocate (ineq,stat=ier) 
      if (ier.ne.0) stop 'symarea.f: Cannot deallocate ineq()'
      end subroutine deallocate_space_for_sym
