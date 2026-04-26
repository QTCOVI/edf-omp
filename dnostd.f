c
c-----------------------------------------------------------------------
c
      subroutine dnostd (nmo,ncent,ngroup,nfugrp,ifugrp,pcut,
     &                   critover,lw,ler,lu18)
c
c.......................................................................
c
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
      real   (kind=8), allocatable,dimension (:,:)   :: sdiag,bdiag
      real   (kind=8), allocatable,dimension (:,:,:) :: sg,sgn
      real   (kind=8), allocatable,dimension (:,:,:) :: aom
      real   (kind=8), allocatable,dimension (:)     :: seigen
      integer(kind=4), allocatable,dimension (:)     :: eignz
      integer(kind=4), allocatable,dimension (:,:)   :: mimael
      integer(kind=4), allocatable,dimension (:,:)   :: resnc
      integer(kind=4), allocatable,dimension (:)     :: mocore
      real   (kind=8)  val1
      integer(kind=4)  nzeig,zeig
      integer(kind=4)  nfugrp(ngroup),ifugrp(ncent,ngroup)
c
      call timer (2,idno,'_dnostd   ',-1)
c
c-----Read AOM
c
      rewind (lu18)
      allocate (aom(ncent,nmo,nmo),stat=ier)
      if (ier.ne.0) stop 'dnostd.f: Cannot allocate aom()'
      call readaomstd (aom,nfugrp,ifugrp,ncent,ngroup,nmo,naom,lu18,ler)
      rewind (lu18)
c
c.....Compute Group Overlap integrals for all groups but the last one.
c
      allocate (sg(ngroup,nmo,nmo),stat=ier)
      if (ier.ne.0) stop 'dnostd.f: Cannot allocate sg()'
      sg(1:ngroup,1:nmo,1:nmo)=0.0d+00
      do i=1,ngroup-1
        do j=1,nfugrp(i)
          k=ifugrp(j,i)
          sg(i,1:nmo,1:nmo)=sg(i,1:nmo,1:nmo)+aom(k,1:nmo,1:nmo)
        enddo
      enddo
c
c-----Compute Group Overlap integrals for the last group
c
      do i=1,nmo
       do j=1,nmo
         if (i.eq.j) delta=1.0d+00
         if (i.ne.j) delta=0.0d+00
         sg(ngroup,i,j)=delta-sum(sg(1:ngroup-1,i,j))
       enddo
      enddo
      allocate (sdiag(nmo,nmo))
      allocate (seigen(nmo))
      allocate (bdiag(nmo,nmo))
      sdiag (:,:)=sg(ngroup,:,:)
      call jacobi (sdiag,nmo,nmo,seigen,bdiag,nrot)
      write (lw,1120) critover,ngroup
      write (lw,222) (seigen(i),i=1,nmo)
c
c-----Determine the non-zero elements of seigen().
c
      allocate (eignz(nmo))
      nzeig=0
      do i=1,nmo
        if (seigen(i).lt.critover) then
          nzeig=nzeig+1
          eignz(nzeig)=i
        endif
      enddo
      zeig=nmo-nzeig
      write (lw,1111) ngroup-1,nzeig,(eignz(k),k=1,nzeig)
c
c-----Reconstruct the Group Overlap Matrix.
c
      allocate (sgn(ngroup,nzeig,nzeig))
      sgn(ngroup,nzeig,nzeig) = 0.0d+00
      do k=1,ngroup-1
        do i=1,nzeig
          i1=eignz(i)
          do j=1,nzeig
            j1=eignz(j)
            val1=0.0d+00
            do l=1,nmo
              do m=1,nmo
                val1=val1+bdiag(l,i1)*bdiag(m,j1)*sg(k,l,m)
              enddo
            enddo
            sgn(k,i,j)=val1
          enddo
        enddo
      enddo
      do i=1,nzeig
       do j=1,nzeig
         if (i.eq.j) delta=1.0d+00
         if (i.ne.j) delta=0.0d+00
         sgn(ngroup,i,j)=delta-sum(sgn(1:ngroup-1,i,j))
       enddo
      enddo

      allocate (mimael(2,ngroup))
      mimael(1,1:ngroup) = 0
      mimael(2,1:ngroup) = nzeig
      allocate (mocore(ngroup))
      mocore(1:ngroup-1)=0
      mocore(ngroup)=zeig
      write (lw,766) ngroup
      do i=1,ngroup
        write (lw,800) i,nfugrp(i),mimael(1,i),mimael(2,i)
      enddo
c
c-----Compute number of probabilities.
c
      npab=1
      do i=1,ngroup-1
        npab=npab*(mimael(2,i)-mimael(1,i)+1)
      enddo
      nprev=npab
      allocate (resnc(nprev,ngroup))
c
c.....Computation of resonance structures.
c
      call rnprobs (mimael,resnc,nzeig,ngroup,npab,lw)
      if (nprev.lt.npab) then
        write (ler,*) 'dnostd.f: NPREV < NPAB returned by rnprobs.f'
        stop
      endif
      write (lw,801) npab
*     do i=1,npab
*       write (lw,802) resnc(i,1:ngroup)
*     enddo
c
c-----Compute EDF.
c
      call compedfstd 
     & (pcut,npab,nprev,ngroup,nzeig,lw,mimael,sgn,resnc,mocore)
c
c-----Deallocate arrays.
c
      deallocate (aom)
      deallocate (sg)
      deallocate (sdiag)
      deallocate (bdiag)
      deallocate (seigen)
      deallocate (eignz)
      deallocate (sgn)
      deallocate (mimael)
      deallocate (resnc)
      deallocate (mocore)
c
      call timer (4,idno,'_dnostd   ',-1)
      return
c
 222  format (5(1x,F17.10))
 1111 format (' # MOs partially localized in the first ',
     &   I4, ' groups = ',I4,1000(/,' # ',20I4))
 766  format (' # NUMBER OF GROUPS = ',I2)
 1120 format (/,' # ',70('+'),/,' #',25x,'DNOEDF MODULE',/,
     & ' # ',70('+'),/,' # Critical Overlap Value = ',E15.8,/,
     & ' # Eigenvalues of AOM in group ',I3)
 800  format (' # GROUP ',i4,' has ',I4,' atoms',4x,
     &    'MinElec,MaxElec (ALPHA or BETA) = ',2I4)
 801  format (' # Number of (alpha or beta) RSRSs is ',I8)
      end
