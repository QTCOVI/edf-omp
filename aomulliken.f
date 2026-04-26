c
c-----------------------------------------------------------------------
c
      subroutine aomulliken (aom,nat,iwrite,lw,wfnfile)
c
c-----------------------------------------------------------------------
c
c     Determine the atomic overlap matrix for Mulliken atoms
c     between natural MOs (nat .eq. 1) or canonical MOs (nat .ne. 1)
c
c-----------------------------------------------------------------------
c
      USE        space_for_wfncoef
      USE        space_for_wfnbasis
      include   'implicit.inc'      
      include   'wfn.inc'      
      include   'param.inc'      
      include   'constants.inc'      
      include   'mline.inc'
*     parameter  (mline = 200)  
      character*(*) wfnfile
      character*(mline) mullfile
c
      real(kind=8)    aom(ncent,nmo,nmo)
      integer    it(3), jt(3)
      real(kind=8)     ax(3), bx(3), ori, orj, xyzint
      real(kind=8)     sabxyz
      real(kind=8), allocatable,dimension (:,:)  :: aomt
      real(kind=8), allocatable,dimension (:)  :: csi,csj
      logical iwrite
c
c-----Test if the AOM file already exists
c
      mullfile=trim(wfnfile)//'.mullikenAOM'
      open (unit=18,file=mullfile,iostat=ierr,status='new')
c
c-----The AOM file exists, so that we read it and continue.
c
      if (ierr.ne.0) then
        open (unit=18,file=mullfile,iostat=ierr,status='old')
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
        write (lw,*)'# SPRIM overlaps computation'
        call gtogto ()
        write (lw,*)'# done'
        sprimok=.true.
      else
        write (lw,*)'# SPRIM overlaps previously available'
      endif
c
c.....Initialize aom()
c
      aom(1:ncent,1:nmo,1:nmo)=zero
      nmois = nmo
      if (nat.eq.1) nmois=0
      mxnp = 0
      do ic=1,ncent
        if (nprimc(ic).gt.mxnp) mxnp=nprimc(ic)
      enddo
      if (.not.allocated(csi)) then
        allocate (csi(mxnp),stat=ier)
        if (ier.ne.0) stop '# aomulliken.f: Cannot allocate csi()'
      endif
      if (.not.allocated(csj)) then
        allocate (csj(mxnp),stat=ier)
        if (ier.ne.0) stop '# aomulliken.f: Cannot allocate csj()'
      endif

      do i=1,nmo
        ii=i+nmois
        do j=1,i
          jj=j+nmois
          npi=1
          npf=nprimc(1)
          do ic=1,ncent
            npt=nprimc(ic)
            do imu=1,npt
              csi(imu)=dot_product(coef(ii,:),sprim(:,npi+imu-1))
              csj(imu)=dot_product(coef(jj,:),sprim(:,npi+imu-1))
            enddo
            aom(ic,i,j) = (dot_product(coef(jj,npi:npf),csi(1:npt)) 
     &                  +  dot_product(coef(ii,npi:npf),csj(1:npt)))/2d0
            npi=npf+1
            if (ic.lt.ncent) npf=npf+nprimc(ic+1)
          enddo
          aom(1:ncent,j,i)=aom(1:ncent,i,j)
        enddo
      enddo
      if (allocated(csi)) then
        deallocate (csi,stat=ier)
        if (ier.ne.0) stop '# aomulliken.f: Cannot deallocate csi()'
      endif
      if (allocated(csj)) then
        deallocate (csj,stat=ier)
        if (ier.ne.0) stop '# aomulliken.f: Cannot deallocate csj()'
      endif
c
c.....Write results to an AOM file called wfnfile'.mullikenAOM'
c
      mullfile=wfnfile(1:leng(wfnfile))//'.mullikenAOM'
      open (18,file=mullfile,iostat=ierr)
      if (ierr.ne.0) then
        write (0,*) 'aomulliken.f: Error opening *.mullikenAOM file'
        return
      endif
      do i=1,ncent
        write (18,'(I5,a)') i,' <--- AOM within this center '
        write (18,80) ((aom(i,m,j),m=1,j),j=1,nmo)
      enddo
      close (18)
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
        write (lw,*) '# Results of Mulliken integrations'
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
 80   format (6(1x,e16.10))
 10   format (6(1x,F15.8))
      return
      end
