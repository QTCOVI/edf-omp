c
c-----------------------------------------------------------------------
c
      subroutine allrhoc (wfnfile,iwfn,ngroup,stdout,aom,nfugrp,ifugrp)
      include     'implicit.inc'
      include     'wfn.inc'
      include     'corr.inc'
      integer, allocatable,dimension (:,:)   :: eleca
      integer, allocatable,dimension (:,:)   :: resnca
      real(kind=8),  allocatable,dimension (:)     :: pnew
      character*(80) res
      real(kind=8) aom(ncent,nmo,nmo)
      integer nfugrp(ngroup)
      integer ifugrp(ncent,ngroup)
      character*(*) wfnfile
      integer iwfn,ngroup
      integer stdout,leng
c
c-----------------------------------------------------------------------
c
      if (.not.allocated(eleca)) allocate (eleca(2,ngroup))
      do i=1,ngroup
        eleca(1,i)=0
        eleca(2,i)=nel-1
      enddo
      mmm=1
      nproba=1
      do i=1,ngroup-1
        nproba=nproba*(nel+i)
        mmm=mmm*i
      enddo
      nproba=nproba/mmm
      if (.not.allocated(resnca)) allocate (resnca(nproba,ngroup))
      call rnprobs (eleca,resnca,nel-1,ngroup,nproba,lw)

      if (.not.allocated(resnca)) allocate (resnca(nproba,ngroup))
      do i=1,nproba
        do j=1,ngroup
          write (res(4*(j-1)+1:4*j),'(I4)') resnca(i,j)
        enddo
        forall(k=4*ngroup+1:80) res(k:k) = ' '
        call rhocond (trim(res),wfnfile,iwfn,ngroup,stdout,aom,
     &       nfugrp,ifugrp)
      enddo
      if (allocated(resnca)) deallocate (resnca)
      if (allocated(eleca))  deallocate (eleca)
      return
      end

