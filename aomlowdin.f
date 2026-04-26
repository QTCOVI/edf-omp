c
c-----------------------------------------------------------------------
c
      subroutine aomlowdin (aom,nat,iwrite,lw,wfnfile)
c
c-----------------------------------------------------------------------
c
c     Determine the atomic overlap matrix for Lowdin atoms.
c
c-----------------------------------------------------------------------
c
      USE        space_for_wfnbasis
      USE        space_for_wfncoef
      include   'implicit.inc'      
      include   'wfn.inc'      
      include   'param.inc'      
      include   'constants.inc'      
      include   'mline.inc'
*     integer,      parameter :: mline = 200
      real(kind=8), parameter :: epsdiag=1d-14
      character (len=*)     wfnfile
      character (len=mline) lowfile
      real(kind=8)  aom(ncent,nmo,nmo)
      real(kind=8)  term
      real(kind=8), allocatable,dimension (:,:)  :: sprimn,uvec,pij
      real(kind=8), allocatable,dimension (:,:)  :: sps,coefn
      real(kind=8), allocatable,dimension (:)    :: diag,xnorm
      real(kind=8), allocatable,dimension (:,:)  :: aomt
      logical iwrite
c
c-----Test if the AOM file already exists
c
      lowfile=trim(wfnfile)//'.lowdinAOM'
      open (unit=18,file=lowfile,iostat=ierr,status='new')
c
c-----The AOM file exists, so that we read it and continue.
c
      if (ierr.ne.0) then
        open (unit=18,file=lowfile,iostat=ierr,status='old')
        do i=1,ncent
          read (18,'(I5,a)') ic
          read (18,80) ((aom(ic,m,j),m=1,j),j=1,nmo)
          do j=1,nmo
            do m=1,j
              aom(ic,j,m)=aom(ic,m,j)
            enddo
          enddo
        enddo
        close (18)
        goto 100
      endif
c
c.....Overlap matrix between Cartesian Gaussian primitives.
c
      if (.not.sprimok) then
        call gtogto ()
        sprimok=.true.
      endif
c
c.....Normalize the primitive Cartesian Gaussians.
c
      if (.not.allocated(sprimn)) then
        allocate (sprimn(nprims,nprims),stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot allocate sprimn()'
      endif
      if (.not.allocated(uvec)) then
        allocate (uvec(nprims,nprims),stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot allocate uvec()'
      endif
      if (.not.allocated(xnorm)) then
        allocate (xnorm(nprims),stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot allocate xnorm()'
      endif
      forall (i=1:nprims) xnorm(i)=1d0/sqrt(sprim(i,i)) 
      forall (i=1:nprims,j=1:nprims) 
     &       sprimn(i,j)=sprim(i,j)*xnorm(i)*xnorm(j)
c
c.....diagonalizes overlap matrix
c
      if (.not.allocated(diag)) then
        allocate (diag(nprims),stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot allocate diag()'
      endif
      call jacobi (sprimn,nprims,nprims,diag,uvec,nrot)
      do i=1,nprims
        if (diag(i).lt.0d0) then
          if (abs(diag(i)).lt.epsdiag) then
            diag(i)=0d0
          else
            write (0,*) '# i,diag(i)=',i,diag(i)
            write (0,*) '# aomlowdin.f: ! Negative eigenvalue of S()'
            stop
          endif
        else
          diag(i)=sqrt(diag(i))
        endif
      enddo
c
      do i=1,nprims
        do j=1,i
          term=0d0
          do k=1,nprims
            term=term+uvec(i,k)*uvec(j,k)*diag(k)
          enddo
          sprimn(i,j)=term
          sprimn(j,i)=term
        enddo
      enddo
c
c.....Determine AOM
c
      if (.not.allocated(pij)) then
        allocate (pij(nprims,nprims),stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot allocate pij()'
      endif
      if (.not.allocated(sps)) then
        allocate (sps(nprims,nprims),stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot allocate sps()'
      endif
      if (.not.allocated(coefn)) then
        allocate (coefn(nmo,nprims),stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot allocate coefn()'
      endif
c
c.....Initialize aom()
c
      aom(1:ncent,1:nmo,1:nmo)=zero
      if (nat .eq. 1) then
        nmois = 0
      else
        nmois = nmo
      endif
      forall (i=1:nmo,j=1:nprims) coefn(i,j)=coef(i+nmois,j)/xnorm(j)
      do k=1,nmo
        do l=1,nmo
          forall (i=1:nprims,j=1:nprims) pij(i,j)=coefn(k,i)*coefn(l,j)
          sps=matmul(matmul(sprimn,pij),sprimn)
          do i=1,nprims
            ic=icen(i)
            aom(ic,k,l)=aom(ic,k,l)+sps(i,i)
          enddo
        enddo
      enddo
c
c.....Write results to an AOM file called wfnfile'.lowdinAOM'
c
      lowfile=wfnfile(1:leng(wfnfile))//'.lowdinAOM'
      open (18,file=lowfile,iostat=ierr)
      if (ierr.ne.0) then
        write (0,*) 'aomlowdin.f: Error opening *.lowdinAOM file'
        return
      endif
      do i=1,ncent
        write (18,'(I5,a)') i,' <--- AOM within this center '
        write (18,80) ((aom(i,m,j),m=1,j),j=1,nmo)
      enddo
      close (18)
c
c.....Deallocate arrays.
c
      if (allocated(sprimn)) then
        deallocate (sprimn,stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot deallocate sprimn()'
      endif
      if (allocated(uvec)) then
        deallocate (uvec,stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot deallocate uvec()'
      endif
      if (allocated(xnorm)) then
        deallocate (xnorm,stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot deallocate xnorm()'
      endif
      if (allocated(diag)) then
        deallocate (diag,stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot deallocate diag()'
      endif
      if (allocated(pij)) then
        deallocate (pij,stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot deallocate pij()'
      endif
      if (allocated(sps)) then
        deallocate (sps,stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot deallocate sps()'
      endif
      if (allocated(coefn)) then
        deallocate (coefn,stat=ier)
        if (ier.ne.0) stop 'aomlowdin.f: Cannot deallocate coefn()'
      endif
c
c-----Print AOM
c
 100  continue
      if (iwrite) then
        if (.not.allocated(aomt)) then
          allocate (aomt(nmo,nmo),stat=ier)
          if (ier.ne.0) stop 'aomindef.f: Cannot allocate aomt()'
        endif
        do m=1,nmo
          do n=1,nmo
            aomt (n,m) = sum(aom(1:ncent,m,n))
          enddo
        enddo
        write (lw,*) '# Results of Lowdin integrations'
        write (lw,*) '#'
        write (lw,*) '# AOM'
        do i=1,ncent
          write (lw,*) '# Center ',i
          do j=1,nmo
            write (lw,10) (aom(i,j,k),k=1,j)
          enddo
        enddo
        write (lw,*) '#'
        write (lw,*) '# AOMT'
        do j=1,nmo
          write (lw,10) (aomt(j,k),k=1,j)
        enddo
      endif
      return
 80   format (6(1x,e16.10))
 10   format (6(1x,F15.8))
      end
