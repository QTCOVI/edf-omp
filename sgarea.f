c
c-----------------------------------------------------------------------
c
      module space_for_sgarea
      save
      real(kind=8), allocatable,dimension (:,:,:) :: sg,sgnat,sgloc
      real(kind=8), allocatable,dimension (:,:,:) :: sgalpha,sgbeta
      real(kind=8), allocatable,dimension (:,:,:) :: sgnatalpha
      real(kind=8), allocatable,dimension (:,:,:) :: sgnatbeta
      end module space_for_sgarea
c
c-----------------------------------------------------------------------
c
      subroutine allocatesg (ngroup,nmo,nalpha,nbeta)
      USE space_for_sgarea
      integer (kind=4) ngroup,nmo,nalpha,nbeta
      allocate (sg(ngroup,nmo,nmo))
      allocate (sgnat(ngroup,nmo,nmo))
      allocate (sgloc(ngroup,nmo,nmo))
      allocate (sgalpha(ngroup,nalpha,nalpha))
      allocate (sgbeta(ngroup,nbeta,nbeta))
      allocate (sgnatalpha(ngroup,nalpha,nalpha))
      allocate (sgnatbeta(ngroup,nbeta,nbeta))
      end subroutine allocatesg
c
c-----------------------------------------------------------------------
c
      subroutine computesg 
     &   (ncent,ngroup,nmo,nalpha,nbeta,aom,nfugrp,ifugrp)
      USE space_for_wfnbasis
      USE space_for_sgarea
      real   (kind=8) aom(ncent,nmo,nmo)
      integer(kind=4) ncent,ngroup,nmo,nalpha,nbeta
      integer(kind=4) nfugrp(ngroup)
      integer(kind=4) ifugrp(ncent,ngroup)
c
c.....Obtain group overlap integrals.
c
      sg=0.0d+00
      if (ngroup.eq.1) then
        forall (m=1:nmo) sg(1,m,m)=1.0d+00
      else
        do k=1,ncent
          do ngr=1,ngroup
            do i=1,nfugrp(ngr)
              if (ifugrp(i,ngr).eq.k) then
                 sg(ngr,:,:)=sg(ngr,:,:)+aom(k,:,:)
              endif
            enddo
          enddo
        enddo
      endif
      do i=1,ngroup
        sgalpha(i,:,:) = sg(i,ialpha(1:nalpha),ialpha(1:nalpha))
        sgbeta (i,:,:) = sg(i,ibeta (1:nbeta ),ibeta (1:nbeta ))
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine computesgnat
     &  (ncent,ngroup,nmo,nalpha,nbeta,nfugrp,ifugrp,v1mata) 
      USE space_for_wfnbasis
      USE space_for_sgarea
      integer(kind=4) :: ncent,ngroup,nmo,nalpha,nbeta
      integer(kind=4) :: nfugrp(ngroup)
      integer(kind=4) :: ifugrp(ncent,ngroup)
      real   (kind=8) :: v1mata(nmo,nmo)

      real(kind=8), allocatable,dimension (:,:)   :: sgtmp
      real(kind=8), allocatable,dimension (:)     :: vtmp

      allocate (sgtmp(nmo,nmo))
      allocate (vtmp(nmo))
      sgnat=0.0d+00
      do igr=1,ngroup
        sgtmp=sg(igr,:,:)
        do k=1,nmo
          vtmp(:)=matmul(sgtmp,v1mata(:,k))
          do i=1,nmo
            sgnat(igr,i,k)=dot_product(v1mata(:,i),vtmp(:))
          enddo
        enddo
      enddo
      deallocate (sgtmp)
      deallocate (vtmp)
      do i=1,ngroup
        sgnatalpha(i,:,:) = sgnat(i,ialpha(1:nalpha),ialpha(1:nalpha))
        sgnatbeta (i,:,:) = sgnat(i,ibeta (1:nbeta ),ibeta (1:nbeta ))
      enddo
      return
      end
