c
c-----------------------------------------------------------------------
c
      SUBROUTINE ofmortho (coef,eta,nmo,nprims,warnofmo,lw)
      USE  space_for_wfnbasis
      implicit none
      integer(kind=4) nmo,nprims,j,info,nat,lw
      integer(kind=4) nmo1(2),nmo2(2)
      integer(kind=4) iofmo
      real   (kind=8) coef(nmo+nmo,nprims)
      real   (kind=8) eta(nmo,nmo)
      real   (kind=8), allocatable,dimension (:,:) :: v,u,ux,cx,vx
      real   (kind=8), allocatable,dimension (:)   :: work,eig
      real   (kind=8)                              :: xeig
      logical warnofmo

      call timer (2,iofmo,'_ofmortho ',-1)
c
c.....allocate arrays.
c
      warnofmo = .false.
      allocate (v(nmo,nprims))
      allocate (vx(nmo,nprims))
      allocate (cx(nmo,nprims))
      allocate (u(nmo,nmo))
      allocate (ux(nmo,nmo))
      allocate (work(3*nmo-1))
      allocate (eig(nmo))
c
c     computes overlaps between primitive cartesian GTOs
c
      if (.not.sprimok) then
        call gtogto ()
        sprimok=.true.
      endif

      nmo1(1)  = nmo + 1
      nmo2(1)  = nmo + nmo
      nmo1(2)  = 1
      nmo2(2)  = nmo
      do nat=1,2
        cx(1:nmo,1:nprims)=coef(nmo1(nat):nmo2(nat),1:nprims)
        v(1:nmo,1:nprims)=cx(1:nmo,1:nprims)
        vx=matmul(v,sprim(1:nprims,1:nprims))
        u=matmul(vx,transpose(v))
*       u=matmul(v,matmul(sprim,transpose(v)))
        call dsyev ('V','L',nmo,u,nmo,eig,work,3*nmo-1,info)
        do j=1,nmo
*         if (abs(eig(j)).lt.1D-10) then
*           warnofmo = .true.
*           ux(:,j)=0.0d+00
*         else
*           ux(:,j)=u(:,j)/sqrt(eig(j))
*         endif
          xeig = max(1.0D-20,eig(j))
          ux(:,j)=u(:,j)/xeig
        enddo
        eta=matmul(ux,transpose(u))
        cx=matmul(eta,v)
        coef(nmo1(nat):nmo2(nat),1:nprims)=cx
      enddo
c
      deallocate (v)
      deallocate (u)
      deallocate (ux)
      deallocate (work)
      deallocate (eig)
      call timer (4,iofmo,'_ofmortho ',-1)
      return
      end
