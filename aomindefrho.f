c
c-----------------------------------------------------------------------
c
c-----Computes the atomic overlap matrix (AOM) on every center (A) of 
c     the molecule (S^A)_ij using a fuzzy Becke partition of 3D space
c     but using as weighting factors 'w_A=rho_A/rho', where 'rho' is the
c     total electron density and 'rho_A' is the MinDef density of atom A.
c
c     (AOM^A)_ij = Int  w_A MO_i x MO_j dq
c
c     It also computes the spatial AOM matrix (SAOM) defined by:
c
c       (SAOM^A)_ij = Int  w_A |MO_i| x |MO_j| dq
c
c     The SAOM is defined as the AOM but changing the integrand  
c     'MO_i x MO_j' by '|MO_i| x |MO_j|', where the vertical bars mean 
c     absolute value.
c
c     The sum of AOMs and SAOMs over all the centers are also obtained.
c
c     AOMT  = SUM_A AOM^A
c     SAOMT = SUM_A SAOM^A
c
c-----------------------------------------------------------------------
c
      subroutine aomindefrho (aom,alphahesel,betahesel,
     &           irho,nrad,nang,irmesh,rhopow,wri,c,lw,wfnfile)
 
      USE space_for_wfnbasis

      implicit none
      include 'wfn.inc'
      include 'datatm.inc'
      include 'mline.inc'
*     real(kind=8) wgatm,rbragg

      real(kind=8), parameter :: pi = 3.141592653589793d0
      real(kind=8), parameter :: fourpi = 4d0*pi
*     integer(kind=4), parameter :: mline = 200
      character*(*) wfnfile
      character*(mline) aomfile

      real(kind=8) c(nmo+nmo,nprims)
      real(kind=8) xmo(nmo)
      integer(kind=4) :: nrad,nang,irmesh,nprims,nmo,ncent,lw
      integer(kind=4) :: i,nangx,ierr
      integer(kind=4), parameter :: numquad = 32
      integer(kind=4), dimension(numquad) :: index
      data (index(i),i=1,numquad) 
     &     /6,  14,  26,  38,  50,  74,  86, 110, 146,
     &    170, 194, 230, 266, 302, 350, 434, 590, 770,
     &    974,1202,1454,1730,2030,2354,2702,3074,3470,
     &   3890,4334,4802,5294,5810/
  
      real(kind=8), allocatable, dimension(:,:,:) :: saom
      real(kind=8), allocatable, dimension(:,:)   :: aomt,saomt
      real(kind=8) :: aom(ncent,nmo,nmo)
      real(kind=8) :: alphahesel,betahesel
  
      real(kind=8) :: smat(ncent,ncent)
      real(kind=8) :: rmid,r,wmindef,rhopow
      real(kind=8) :: pmos,pabs,weight
      integer(kind=4) :: j,k,m,n,ii
      real(kind=8), allocatable :: rads(:), wrads(:)
      real(kind=8), allocatable :: xang(:), yang(:), zang(:), wang(:)
      integer(kind=4) :: ir,il,imindef,istat,irho
      real(kind=8)    :: x(3),xr(3)
      logical         :: wri
c
c-----------------------------------------------------------------------
c
      call timer (2,imindef,'_aomindef_',-1)
  
      nangx=nang
      if (nangx.ge.index(numquad)) then
        nang=index(numquad)
      elseif  (nangx.le.index(1)) then
        nang=index(1)
      else
        do i=numquad,2,-1
          if (nangx.ge.index(i-1).and.nangx.lt.index(i)) then
             nang=index(i-1)
          endif
        enddo
      endif
c
c-----Generate the mesh and compute the AOM integrals
c
      allocate(saom(ncent,nmo,nmo),stat=istat)
      if (istat /= 0) then
        stop 'aomindefrho.f: Cannot allocate saom() array'
      endif
      allocate(saomt(nmo,nmo),stat=istat)
      if (istat /= 0) then
        stop 'aomindefrho.f: Cannot allocate saomt() array'
      endif
      allocate(aomt(nmo,nmo),stat=istat)
      if (istat /= 0) then
        stop 'aomindefrho.f: Cannot allocate aomt() array'
      endif

      aom   = 0d0
      saom  = 0d0
      saomt = 0d0
      aomt  = 0d0
c
c-------Radial part
c
      allocate(rads(nrad),wrads(nrad),stat=istat)
      if (istat /= 0) then
        stop 'aomindefrho.f: Cannot allocate radial meshes'
      endif
c
c-----Angular part
c
      allocate (xang(nang),yang(nang),zang(nang),stat=istat)
      if (istat /= 0) then
        stop 'aomindefrho.f: Cannot allocate angular meshes'
      endif
      if (.not.allocated(wang)) then
        allocate (wang(nang),stat=istat)
        if (istat /= 0) then
          stop 'aomindefrho.f: Cannot allocate angular weights'
        endif
      endif
c
      do i = 1, ncent
        write (0,*)'# AOM IN CENTER ',i
        rmid = rbragg(int(charge(i),4))
        if (irmesh .eq. 1) then
          call rmeshb (nrad,rmid,rads,wrads)
        else if (irmesh .eq. 2) then
          call rmeshh1 (nrad,rmid,rads,wrads)
        else if (irmesh .eq. 3) then
          call rmeshh2 (nrad,rmid,rads,wrads)
        else
          call rmeshb (nrad,rmid,rads,wrads)
        end if
c
        if (nang == 6) then
          call ld0006 (xang,yang,zang,wang,nang)
        else if (nang == 14) then
          call ld0014 (xang,yang,zang,wang,nang)
        else if (nang == 26) then
          call ld0026 (xang,yang,zang,wang,nang)
        else if (nang == 38) then
          call ld0038 (xang,yang,zang,wang,nang)
        else if (nang == 50) then
          call ld0050 (xang,yang,zang,wang,nang)
        else if (nang == 74) then
          call ld0074 (xang,yang,zang,wang,nang)
        else if (nang == 86) then
          call ld0086 (xang,yang,zang,wang,nang)
        else if (nang == 110) then
          call ld0110 (xang,yang,zang,wang,nang)
        else if (nang == 146) then
          call ld0146 (xang,yang,zang,wang,nang)
        else if (nang == 170) then
          call ld0170 (xang,yang,zang,wang,nang)
        else if (nang == 194) then
          call ld0194 (xang,yang,zang,wang,nang)
        else if (nang == 230) then
          call ld0230 (xang,yang,zang,wang,nang)
        else if (nang == 266) then
          call ld0266 (xang,yang,zang,wang,nang)
        else if (nang == 302) then
          call ld0302 (xang,yang,zang,wang,nang)
        else if (nang == 350) then
          call ld0350 (xang,yang,zang,wang,nang)
        else if (nang == 434) then
          call ld0434 (xang,yang,zang,wang,nang)
        else if (nang == 590) then
          call ld0590 (xang,yang,zang,wang,nang)
        else if (nang == 770) then
          call ld0770 (xang,yang,zang,wang,nang)
        else if (nang == 974) then
          call ld0974 (xang,yang,zang,wang,nang)
        else if (nang == 1202) then
          call ld1202 (xang,yang,zang,wang,nang)
        else if (nang == 1454) then
          call ld1454 (xang,yang,zang,wang,nang)
        else if (nang == 1730) then
          call ld1730 (xang,yang,zang,wang,nang)
        else if (nang == 2030) then
          call ld2030 (xang,yang,zang,wang,nang)
        else if (nang == 2354) then
          call ld2354 (xang,yang,zang,wang,nang)
        else if (nang == 2702) then
          call ld2702 (xang,yang,zang,wang,nang)
        else if (nang == 3074) then
          call ld3074 (xang,yang,zang,wang,nang)
        else if (nang == 3470) then
          call ld3470 (xang,yang,zang,wang,nang)
        else if (nang == 3890) then
          call ld3890 (xang,yang,zang,wang,nang)
        else if (nang == 4334) then
          call ld4334 (xang,yang,zang,wang,nang)
        else if (nang == 4802) then
          call ld4802 (xang,yang,zang,wang,nang)
        else if (nang == 5294) then
          call ld5294 (xang,yang,zang,wang,nang)
        else if (nang == 5810) then
          call ld5810 (xang,yang,zang,wang,nang)
        else
          stop 'aomindefrho.f: Unknown value of nang'
        end if
c
c-------Set origin and mix angular and radial part
c
        do ir = 1,nrad
          r = rads(ir)
          if (irho.eq.4) wmindef = 1d0-erf(alphahesel*r**betahesel)
          do il = 1, nang
            xr = r * (/xang(il),yang(il),zang(il)/)
            x  = xyz(i,:) + xr
            if (irho.eq.1.or.irho.eq.4) then
              call mindefw  (x,i,wmindef,irho,c,xmo)
            elseif (irho.eq.2) then
              call netrhow  (x,i,wmindef,c,xmo)
            elseif (irho.eq.3) then
              call promrhow (x,i,rhopow,wmindef,c,xmo)
            endif
            weight = wmindef * wrads(ir) * wang(il)
            do m=1,nmo
              do n=1,m
                pmos=xmo(m)*xmo(n)
                pabs=abs(pmos)
                aom(i,m,n)  = aom(i,m,n)  + pmos * weight
                saom(i,m,n) = saom(i,m,n) + pabs * weight
              enddo
            enddo
          end do
        end do
      end do
c
      if (allocated(rads)) then
        deallocate (rads,stat=istat)
        if (istat /= 0) then
          stop 'aomindefrho.f: Cannot deallocate rads() array'
        endif
      endif
      if (allocated(wrads)) then
        deallocate (wrads,stat=istat)
        if (istat /= 0) then
          stop 'aomindefrho.f: Cannot deallocate wrads() array'
        endif
      endif
      if (allocated(xang)) then
        deallocate (xang,stat=istat)
        if (istat /= 0) then
          stop 'aomindefrho.f: Cannot deallocate xang() array'
        endif
      endif
      if (allocated(yang)) then
        deallocate (yang,stat=istat)
        if (istat /= 0) then
          stop 'aomindefrho.f: Cannot deallocate yang() array'
        endif
      endif
      if (allocated(zang)) then
        deallocate (zang,stat=istat)
        if (istat /= 0) then
          stop 'aomindefrho.f: Cannot deallocate zang() array'
        endif
      endif
      if (allocated(wang)) then
        deallocate (wang,stat=istat)
        if (istat /= 0) then
          stop 'aomindefrho.f: Cannot deallocate wang() array'
        endif
      endif

      do m=1,nmo
        do n=1,m
          aom (1:ncent,n,m)  = aom(1:ncent,m,n)
          saom(1:ncent,n,m) = saom(1:ncent,m,n)
        enddo
      enddo
      do m=1,nmo
        do n=1,nmo
          aomt (n,m) = sum( aom(1:ncent,m,n))
          saomt(n,m) = sum(saom(1:ncent,m,n))
        enddo
      enddo
c
c.....Write results to an AOM file called wfnfile'.beckeAOM'
c
      if (irho.eq.1) aomfile=trim(wfnfile)//'.mindefrhoAOM'
      if (irho.eq.2) aomfile=trim(wfnfile)//'.netrhoAOM'
      if (irho.eq.3) aomfile=trim(wfnfile)//'.promrhoAOM'
      if (irho.eq.4) aomfile=trim(wfnfile)//'.heselAOM'
      open (18,file=trim(aomfile),iostat=ierr)
      if (ierr.ne.0) then
        write (0,*) 'aomindefrho.f: Error opening AOM file'
        return
      endif
      do i=1,ncent
        write (18,'(I5,a)') i,' <--- AOM within this center '
        write (18,80) ((aom(i,m,j),m=1,j),j=1,nmo)
      enddo
      close (18)
c
c-----write 
c
      if (wri) then
        write (lw,501)
        write (lw,*) '#'
        if (irho.eq.1) then
          write (lw,*) '# Results of MinDef-RHO-like integrations'
        endif
        if (irho.eq.2) then
          write (lw,*) '# Results of NETRHO-like integrations'
        endif
        if (irho.eq.3) then
          write (lw,*) '# Results of PROMRHO-like integrations'
        endif
        if (irho.eq.4) then
          write (lw,*) '# Results of Heselmann-like integrations'
        endif
        write (lw,*) '# Radial  points              = ',nrad
        write (lw,*) '# Angular points              = ',nang
        write (lw,*) '# Radial mapping              = ',irmesh
        if (irho.eq.3) then
          write (lw,*) '# RHO_A power                 = ',rhopow
        endif
        if (irho.eq.4) then
          write (lw,*) '# ALPHA_Heselmann             = ',alphahesel
          write (lw,*) '#  BETA_Heselmann             = ',betahesel
        endif
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
        write (lw,*) '#'
        write (lw,*) '# Spatial AOM'
        write (lw,*) '# [SAOM_ij]_A = <abs(Phi_i)|abs(Phi_j)>_A'
        do i=1,ncent
          write (lw,*) '# Center ',i
          do j=1,nmo
            write (lw,10) (saom(i,j,k),k=1,j)
          enddo
        enddo
        write (lw,*) '#'
        write (lw,*) '# Spatial AOMT'
        do j=1,nmo
          write (lw,10) (saomt(j,k),k=1,j)
        enddo
      endif

      deallocate(saom,stat=istat)
      if (istat /= 0) then
        stop 'aomindefrho.f: Cannot deallocate saom() array'
      endif
      deallocate(saomt,stat=istat)
      if (istat /= 0) then
        stop 'aomindefrho.f: Cannot deallocate saomt() array'
      endif
      deallocate(aomt,stat=istat)
      if (istat /= 0) then
        stop 'aomindefrho.f: Cannot deallocate aomt() array'
      endif

 80   format (6(1x,e16.10))
 10   format (6(1x,F15.8))
 501  format (1x,69('-'))

      call timer (4,imindef,'_aomindef_',-1)
      end 
