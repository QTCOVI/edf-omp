c
c-----------------------------------------------------------------------
c
      subroutine openrho 
     &   (c,aom,nfugrp,ifugrp,nprims,ncent,nmo,ngroup,epsneg,lw,lerr,
     &   ifilc,udat,wfnfile,largwr)
      USE          space_for_wfncoef
      USE          space_for_wfnbasis
      include     'implicit.inc'
      include     'constants.inc'
      include     'mline.inc'
      real(kind=8) c(nmo,nmo)
      real(kind=8) aom(ncent,nmo,nmo)
      real(kind=8) epsneg
      integer(kind=4) nfugrp(ngroup)
      integer(kind=4) ifugrp(ncent,ngroup)
      logical largwr
      real(kind=8), allocatable,dimension (:,:)  :: sg,shalf,rhop,dno
      real(kind=8), allocatable,dimension (:,:)  :: vecprop,sminushalf
      real(kind=8), allocatable,dimension (:,:)  :: vecan
      real(kind=8), allocatable,dimension (:)    :: valprop,sumeigen
      real(kind=8), allocatable,dimension (:)    :: teigen,percent
      integer(kind=4), allocatable,dimension (:) :: ilabeig,iord
      real(kind=8), parameter ::   eps   = 1D-14
*     integer(kind=4), parameter ::   mline   =  200
      integer(kind=4) udat,udatnw
      character*(mline) wfnloc
      character*(*) wfnfile
      character*(4) fourchar
c
      if (.not.allocated(sg)) then
        allocate (sg(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate sg()'
      endif
      if (.not.allocated(shalf)) then
        allocate (shalf(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate shalf()'
      endif
      if (.not.allocated(sminushalf)) then
        allocate (sminushalf(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate sminushalf()'
      endif
      if (.not.allocated(rhop)) then
        allocate (rhop(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate rhop()'
      endif
      if (.not.allocated(vecprop)) then
        allocate (vecprop(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate vecprop()'
      endif
      if (.not.allocated(vecan)) then
        allocate (vecan(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate vecan()'
      endif
      if (.not.allocated(dno)) then
        allocate (dno(nmo+nmo,nprims),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate dno()'
      endif
      if (.not.allocated(valprop)) then
        allocate (valprop(nmo),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate valprop()'
      endif
      if (.not.allocated(sumeigen)) then
        allocate (sumeigen(nmo),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate sumeigen()'
      endif
      if (.not.allocated(teigen)) then
        allocate (teigen(nmo*ngroup),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate teigen()'
      endif
      if (.not.allocated(percent)) then
        allocate (percent(nmo),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate percent()'
      endif
      if (.not.allocated(ilabeig)) then
        allocate (ilabeig(nmo*ngroup),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate ilabeig()'
      endif
      if (.not.allocated(iord)) then
        allocate (iord(nmo*ngroup),stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot allocate iord()'
      endif
c
      write (lw,100)
c
c.....Obtain group overlap integrals.
c
      inmo   = nmo+nmo
      inpr   = nprims
      udatnw = 2
c
c-----Run over the different groups
c
      sumeigen = zero
      init=1
      ifin=nmo
      do ig=1,ngroup
        sg = zero
        do j=1,nfugrp(ig)
          sg(:,:)=sg(:,:)+aom(ifugrp(j,ig),:,:)
        enddo
c
c-------Diagonalize group overlap matrix
c
        call jacobi (sg,nmo,nmo,valprop,vecprop,nrot)
        shalf = zero
        sminushalf = zero
        do k=1,nmo
          if (valprop(k).gt.zero) then
          elseif (valprop(k).gt.epsneg) then
            write (lerr,222) k,valprop(k),eps+eps
            valprop(k) = eps+eps
          elseif (abs(valprop(k)).lt.eps) then
            write (lerr,223) k,valprop(k),eps+eps
            valprop(k) = eps+eps
          else
            write (lerr,*) '# openrho.f: Negative eigenvalue of sg()'
            stop
          endif
          do i=1,nmo
            do j=1,i
              pvec=vecprop(i,k)*vecprop(j,k)
              shalf(i,j)=shalf(i,j)+pvec*sqrt(valprop(k))
              sminushalf(i,j)=sminushalf(i,j)+pvec/sqrt(valprop(k))
            enddo
          enddo
        enddo
        do i=1,nmo
          do j=1,i
            shalf(j,i)=shalf(i,j)
            sminushalf(j,i)=sminushalf(i,j)
          enddo
        enddo
c
c-------Built S^1/2 * 1-RDM * S^1/2
c
        rhop = matmul (shalf,matmul(c,transpose(shalf)))
c
c-------Diagonalize the proyected density matrix S^1/2 * 1-RDM * S^1/2
c       and obtain the eigenvalues lambda.
c
        call jacobi (rhop,nmo,nmo,valprop,vecprop,nrot)
c
c-------Accumulate eigenvalues
c
        teigen(init:ifin)  = valprop(1:nmo)
        ilabeig(init:ifin) = ig
        init = init+nmo
        ifin = ifin+nmo
c
c-------Add eigenvalues of this group to the total
c
        sumeigen = sumeigen+valprop
        write (lw,500) ig,nfugrp(ig),(ifugrp(j,ig),j=1,nfugrp(ig))
        write (lw,*) '#'
        write (lw,*) '# FRAGMENT First-order density matrix Eigenvalues'
        write (lw,445) (valprop(i),i=1,nmo)
        write (lw,446) sum(valprop(1:nmo))
 446    format (' # ELECTRONS =   ',F15.8)
c
c-------Matrix to express the DNOs in the canonical basis, S^{-1/2} \tilde{U}.
c
        vecan = matmul(sminushalf,vecprop)
c
c-------Compute S^{-1}
c
        sminushalf=matmul(sminushalf,sminushalf)
c
        if (largwr) then
          write (lw,'(a)') ' #'
          write (lw,'(a)') ' # Eigenvectors in the canonical MO basis'
          write (lw,'(a)') ' #'
        endif
        do i=1,nmo
          if (largwr) then
            write (lw,'(a,I3)') ' # Eigenvector ',i
            write (lw,444) (vecan(j,i),j=1,nmo)
          endif
          do j=1,nprims
            tmp = zero
            do k=1,nmo
              tmp = tmp + vecan(k,i) * coef(k+nmo,j)
            enddo
            dno(i,j) = tmp
          enddo
        enddo
c
c-------Overlaps of DNOs in R^3
c
        vecan = matmul(transpose(vecan),vecan)   ! vecan is now <dno|dno>_R^3
        do i=1,nmo
          percent(i)=100d0/vecan(i,i)
        enddo
        if (largwr) then
          write (lw,'(a)') ' #'
          write (lw,'(a)') ' # Overlaps of DNOs in R^3'
          write (lw,'(a)') ' #'
          do i=1,nmo
            percent(i)=100d0/vecan(i,i)
            write (lw,444) (vecan(i,j),j=1,i)
          enddo
          write (lw,'(a)') ' #'
          write (lw,'(a)') ' #'
          do i=1,nmo
            write (lw,'(a,I3,a)') ' # DNO ',i,' in the primitive basis'
            write (lw,444) (dno(i,j),j=1,nprims)
          enddo
        endif
        write (lw,26) (percent(i),i=1,nmo)
c
c-------Write a WFN file with DNOs. The i^th DNO is multiplied by 
c       vecan(i,i)^{-1/2} such that it is normalized to 1.0 in R^3
c
        occ(1:nmo)=valprop(1:nmo)
        do i=1,nmo
          xnorm=1.0D0/sqrt(vecan(i,i))
          forall (j=1:nprims) dno(i,j)=dno(i,j)*xnorm
          forall (j=1:nprims) dno(i+nmo,j)=coef(i+nmo,j)
        enddo
        wfnloc = trim(wfnfile)//"-frag"//fourchar(ig)
        write (lw,'(a)') ' # '
        write (lw,22) trim(wfnloc)
        write (lw,'(a)') ' # '
        open (unit = udatnw,file = trim(wfnloc),status='unknown')
        call wrtwfn (udat,udatnw,lerr,ifilc,dno,inmo,inpr)
        close (unit = udatnw)
      enddo
      write (lw,23) sum(sumeigen(1:nmo))
c
c-----Order the accumulated eigenvalues
c
      forall (i=1:ngroup*nmo) iord(i)=i
      call qcksort (teigen, iord, 1,ngroup*nmo)
      write (lw,240) 
      do i=ngroup*nmo,1,-1
        ii=iord(i)
        write (lw,24) i,teigen(ii),ilabeig(ii)
      enddo
      write (lw,101)
c

      if (allocated(sg)) then
        deallocate (sg,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate sg()'
      endif
      if (allocated(shalf)) then
        deallocate (shalf,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate shalf()'
      endif
      if (allocated(sminushalf)) then
        deallocate (sminushalf,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate sminushalf()'
      endif
      if (allocated(rhop)) then
        deallocate (rhop,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate rhop()'
      endif
      if (allocated(vecprop)) then
        deallocate (vecprop,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate vecprop()'
      endif
      if (allocated(vecan)) then
        deallocate (vecan,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate vecan()'
      endif
      if (allocated(dno)) then
        deallocate (dno,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate dno()'
      endif
      if (allocated(valprop)) then
        deallocate (valprop,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate valprop()'
      endif
      if (allocated(sumeigen)) then
        deallocate (sumeigen,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate sumeigen()'
      endif
      if (allocated(teigen)) then
        deallocate (teigen,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate teigen()'
      endif
      if (allocated(percent)) then
        deallocate (percent,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate percent()'
      endif
      if (allocated(ilabeig)) then
        deallocate (ilabeig,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate ilabeig()'
      endif
      if (allocated(iord)) then
        deallocate (iord,stat=ier)
        if (ier.ne.0) stop 'openrho.f: Cannot deallocate iord()'
      endif
      return
 100  format (//,
     & ' ',/,
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' +                                                         +',/,
     & ' + B E G I N   O P E N   S Y S T E M S   A N A L Y S I S   +',/,
     & ' +                                                         +',/,
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' ')
 101  format (//,
     & ' ',/,
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' +                                                         +',/,
     & ' +    E N D   O P E N   S Y S T E M S   A N A L Y S I S    +',/,
     & ' +                                                         +',/,
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' ',4/)
 444  format (5(1x,E15.8))
 445  format (5(1x,F15.8))
 20   format (///' # ',77('-'),/,' # FRAGMENT ',I2)
 21   format (' # ',77('-'))
 22   format (' # ',77('-'),/,
     & ' # Writing file "',a,'". Each DNO is divided by <DNO|DNO>^{1/2}'
     & /' # such that it is normalized to 1.0 in R^3',/,' # ',77('-'))
 23   format (' # Total number of electrons, N = ',F25.18,/,' #')
 240  format (//' # Ordered Eigenvalues of S^1/2 \rho S^1/2',/,
     &          ' # Number     Eigenvalue   Group',/,
     &          ' # -----------------------------')
 24   format (I7,2x,F15.8,2x,I4)
 500  format (///' # FRAGMENT ',I2,' FORMED BY ',I3,' BASINS:',
     & 15I4,/,1000(' # ',33x,15I4,/))
 26   format (/' # Percentage of location of the DNOs in the group',/,
     &   1000(' # ',10(1x,F6.2),/))
 222  format (' # Negative eigenvalue ',I4,' : ',E15.8,
     &  ' is set to ',E17.10)
 223  format (' # Too low eigenvalue ',I4,' of sg(): ',
     &  E15.8,', is set to',E15.8)
 1124 format (' # Eigenvalue ',I4,' of sg() is negative, ',E17.10)
      end
