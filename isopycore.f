c
c-----------------------------------------------------------------------
c
      subroutine isopycore (aom,lw,largwr)
c 
c.......................................................................
c
      USE          space_for_wfnbasis
      USE          space_for_wfncoef
      include     'implicit.inc'
      include     'param.inc'
      include     'wfn.inc'
      include     'corr.inc'
      include     'stderr.inc'
      include     'constants.inc'
      real(kind=8), allocatable,dimension (:,:)    :: vrot
      real(kind=8), allocatable,dimension (:,:,:)  :: aomx
      real(kind=8), allocatable,dimension (:,:)    :: aomcv
      real(kind=8), allocatable,dimension (:,:)    :: cc
      real(kind=8), allocatable,dimension (:)      :: oc

      real(kind=8)     aom(ncent,nmo,nmo)
      logical     largwr,okw
c
      okw = .false.
      if (.not.allocated(vrot)) then
        allocate (vrot(ncore,ncore),stat=ier)
        if (ier.ne.0) stop 'isopycore.f: Cannot allocate vrot()'
      endif
      if (.not.allocated(aomx)) then
        allocate (aomx(1:ncent,ncore,ncore),stat=ier)
        if (ier.ne.0) stop 'isopycore.f: Cannot allocate aomx()'
      endif
      if (.not.allocated(aomcv)) then
        allocate (aomcv(ncore,nmo-ncore),stat=ier)
        if (ier.ne.0) stop 'isopycore.f: Cannot allocate aomcv()'
      endif
      if (.not.allocated(oc)) then
        allocate (oc(ncore),stat=ier)
        if (ier.ne.0) stop 'isopycore.f: Cannot allocate oc()'
      endif
      aomx(1:ncent,1:ncore,1:ncore)=aom(1:ncent,1:ncore,1:ncore)
      oc(1:ncore)=occ(1:ncore)
      write (lw,642,advance='no')
      call ruedmis (aomx,vrot,oc,ncore,ncent,lw,largwr,okw)
      aom(1:ncent,1:ncore,1:ncore)=aomx(1:ncent,1:ncore,1:ncore)
      occ(1:ncore)=oc(1:ncore)
c
c.....Obtain the AOM between OLD valence electrons and NEW CORE 
c     localized MOs
c
      do ic=1,ncent
        do i=1,ncore
          do j=ncore+1,nmo
            tmp = zero
            do k=1,ncore
              tmp = tmp + vrot(i,k) * aom (ic,k,j)
            enddo
            aomcv(i,j-ncore)=tmp
          enddo
        enddo
        do i=1,ncore
          do j=ncore+1,nmo
            aom(ic,i,j)=aomcv(i,j-ncore)
            aom(ic,j,i)=aom(ic,i,j)
          enddo
        enddo
      enddo
c
c.....Expand the localized CORE MOs in terms of primitive Gaussians.
c      
      if (.not.allocated(cc)) then
        allocate (cc(ncore,nprims),stat=ier)
        if (ier.ne.0) stop 'isopycore.f: Cannot allocate cc()'
      endif
      do i=1,ncore
        do j=1,nprims
          temp = zero
          do k=1,ncore
            temp = temp + vrot(i,k) * coef(k+nmo,j)
          enddo
          cc(i,j) = temp
        enddo
      enddo
      coef(1:ncore,1:nprims)=cc(1:ncore,1:nprims)
      coef(nmo+1:nmo+ncore,1:nprims)=cc(1:ncore,1:nprims)

      if (allocated(oc)) then
        deallocate (oc,stat=ier)
        if (ier.ne.0) stop 'isopycore.f: Cannot deallocate oc()'
      endif
      if (allocated(aomx)) then
        deallocate (aomx,stat=ier)
        if (ier.ne.0) stop 'isopycore.f: Cannot deallocate aomx()'
      endif
      if (allocated(aomcv)) then
        deallocate (aomcv,stat=ier)
        if (ier.ne.0) stop 'isopycore.f: Cannot deallocate aomcv()'
      endif
      if (allocated(cc)) then
        deallocate (cc,stat=ier)
        if (ier.ne.0) stop 'isopycore.f: Cannot deallocate cc()'
      endif
      if (allocated(vrot)) then
        deallocate (vrot,stat=ier)
        if (ier.ne.0) stop 'isopycore.f: Cannot deallocate vrot()'
      endif
 642  format (1x,'# Isopycnic Localization of CORE MOs: ')
      return
      end
