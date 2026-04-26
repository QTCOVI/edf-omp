      module space_for_primgto
      save
      real(kind=8), allocatable, dimension (:)     :: xnorm
      integer, allocatable, dimension (:)          :: npc,ngroup
      integer, allocatable, dimension (:,:)        :: icenat,nzexp
      integer, allocatable, dimension (:,:,:)      :: nuexp
      end module space_for_primgto
c-----------------------------------------------------------------------
      subroutine allocate_space_for_primgto (nprims,ncent,mgrp,ngtoG)
      USE        space_for_primgto
c
      if (.not.allocated(xnorm)) then
        allocate (xnorm(nprims),stat=ier) 
        if (ier.ne.0) then
          stop 'primtgto.f: Already allocated xnorm() array ?'
        endif
      endif
c
      if (.not.allocated(npc)) then
        allocate (npc(ncent),stat=ier) 
        if (ier.ne.0) then
          stop 'primtgto.f: Already allocated npc() array ?'
        endif
      endif
c
      if (.not.allocated(ngroup)) then
        allocate (ngroup(ncent),stat=ier) 
        if (ier.ne.0) then
          stop 'primtgto.f: Already allocated ngroup() array ?'
        endif
      endif
c
      if (.not.allocated(icenat)) then
        allocate (icenat(nprims,ncent),stat=ier) 
        if (ier.ne.0) then
          stop 'primtgto.f: Already allocated icenat() array ?'
        endif
      endif
c
      if (.not.allocated(nzexp)) then
        allocate (nzexp(ncent,mgrp),stat=ier) 
        if (ier.ne.0) then
          stop 'primtgto.f: Already allocated nzexp() array ?'
        endif
      endif
c
      if (.not.allocated(nuexp)) then
        allocate (nuexp(ncent,mgrp,ngtoG),stat=ier) 
        if (ier.ne.0) then
          stop 'primtgto.f: Already allocated nuexp() array ?'
        endif
      endif
      end subroutine allocate_space_for_primgto
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_primgto
      USE        space_for_primgto
c
      if (allocated(xnorm)) then
        deallocate (xnorm,stat=ier) 
        if (ier.ne.0) then
          stop 'primtgto.f: Already deallocated xnorm() array ?'
        endif
      endif
c
      if (allocated(npc)) then
        deallocate (npc,stat=ier) 
        if (ier.ne.0) then
          stop 'primtgto.f: Already deallocated npc() array ?'
        endif
      endif
c
      if (allocated(ngroup)) then
        deallocate (ngroup,stat=ier) 
        if (ier.ne.0) then
          stop 'primtgto.f: Already deallocated ngroup() array ?'
        endif
      endif
c
      if (allocated(icenat)) then
        deallocate (icenat,stat=ier) 
        if (ier.ne.0) then
          stop 'primtgto.f: Already deallocated icenat() array ?'
        endif
      endif
c
      if (allocated(nzexp)) then
        deallocate (nzexp,stat=ier) 
        if (ier.ne.0) then
          stop 'primtgto.f: Already deallocated nzexp() array ?'
        endif
      endif
c
      if (allocated(nuexp)) then
        deallocate (nuexp,stat=ier) 
        if (ier.ne.0) then
          stop 'primtgto.f: Already deallocated nuexp() array ?'
        endif
      endif
      end subroutine deallocate_space_for_primgto
