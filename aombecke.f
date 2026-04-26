c
c-----------------------------------------------------------------------
c
c-----Computes the atomic overlap matrix (AOM) on every center (A) of 
c     the molecule (S^A)_ij using a fuzzy Becke partition of 3D space;
c
c     (AOM^A)_ij = Int_A  MO_i x MO_j dq
c
c     It also computes the spatial AOM matrix (SAOM) defined by:
c
c       (SAOM^A)_ij = Int_A  |MO_i| x |MO_j| dq
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
      subroutine aombecke (aom,nrad,nang,pow,irmesh,verbose,c,lw,wfnf)
 
      USE space_for_wfnbasis

      implicit none
      include 'wfn.inc'
      include 'datatm.inc'
      include 'mline.inc'

      real(kind=8), parameter :: pi = 3.141592653589793d0
      real(kind=8), parameter :: cutoffa = 0.7d0
      real(kind=8), parameter :: fourpi = 4d0*pi
      character*(*) wfnf
      character*(mline) beckefile

      real(kind=8) c(nmo,nprims)
      real(kind=8) xmo(nmo)
      integer(kind=4) :: nrad,nang,pow,irmesh,nprims,nmo,ncent,lw
      integer(kind=4) :: i,atom,nangx,ierr
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
  
      real(kind=8) :: rr(ncent,ncent),dxyz(3)
      real(kind=8) :: smat(ncent,ncent)
      real(kind=8) :: pvec(ncent)
      real(kind=8) :: rmid,r,r1,r2,rmiu,rmiu2,tmps,beckew,chi,uij,aij
      real(kind=8) :: pmos,pabs,weight,wthiscenter,rmiu2n,xionemu
      integer(kind=4) :: j,k,m,n,kk,ii
      real(kind=8), allocatable :: rads(:), wrads(:)
      real(kind=8), allocatable :: xang(:), yang(:), zang(:), wang(:)
      integer(kind=4) :: ir,il,imeshb,istat,ib
      real(kind=8)    :: x(3)
      logical         :: verbose
c
c-----------------------------------------------------------------------
c
      call timer (2,imeshb,'_aombecke_',-1)
  
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
c-----write
c
      if (verbose) then
        write (lw,501)
        write (lw,*) '# Results of Becke-like integrations'
        write (lw,*) '# Radial  points              = ',nrad
        write (lw,*) '# Angular points              = ',nang
        write (lw,*) '# Radial mapping              = ',irmesh
        write (lw,*) '# Becke K iterative parameter = ',pow
        write (lw,*) '#'
      endif
c
c-----Compute interatomic distances.
c
      rr = 0.0d+00
      do i = 1,ncent
        do j = i+1,ncent
          dxyz(:) = xyz(i,:)-xyz(j,:)
          rr(i,j) = sqrt(dot_product(dxyz,dxyz))
          rr(j,i) = rr(i,j)
        end do
      end do
c
c-----Generate the mesh and compute the AOM integrals
c
      allocate(saom(ncent,nmo,nmo),stat=istat)
      if (istat /= 0) then
        stop 'aombecke.f: Cannot allocate saom() array'
      endif
      allocate(saomt(nmo,nmo),stat=istat)
      if (istat /= 0) then
        stop 'aombecke.f: Cannot allocate saomt() array'
      endif
      allocate(aomt(nmo,nmo),stat=istat)
      if (istat /= 0) then
        stop 'aombecke.f: Cannot allocate aomt() array'
      endif
      aom   = 0.0d+00
      saom  = 0.0d+00
      saomt = 0.0d+00
      aomt  = 0.0d+00
      kk    = 0
c
c-------Radial part
c
      allocate(rads(nrad),wrads(nrad),stat=istat)
      if (istat /= 0) then
        stop 'aombecke.f: Cannot allocate radial meshes'
      endif
c
c-----Angular part
c
      allocate (xang(nang),yang(nang),zang(nang),stat=istat)
      if (istat /= 0) then
        stop 'aombecke.f: Cannot allocate angular meshes'
      endif
      if (.not.allocated(wang)) then
        allocate (wang(nang),stat=istat)
        if (istat /= 0) then
          stop 'aombecke.f: Cannot allocate angular weights'
        endif
      endif
c
      do atom = 1, ncent
        write (0,'(a,I3)') ' # AOM in center ',atom
        rmid = rbragg(int(charge(atom),4))
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
          stop 'aombecke.f: Unknown value of nang'
        end if
c
c-------Set origin and mix angular and radial part
c
        do ir = 1,nrad
          r = rads(ir)
          do il = 1, nang
            x = xyz(atom,:) + r * (/xang(il),yang(il),zang(il)/)
CCCC        call wbec (atom,pow,x,weight,.true.,.false.,.false.)
            smat = 1.0d+00
            do j = 1,ncent
              r1 = sqrt(dot_product(x(:)-xyz(j,:),x(:)-xyz(j,:)))
              do k = 1,ncent
                if (j==k) cycle
                r2 = sqrt(dot_product(x(:)-xyz(k,:),x(:)-xyz(k,:)))
                rmiu = (r1-r2)/rr(j,k)
                chi = rbragg(int(charge(j),4))/rbragg(int(charge(k),4))
                uij = (chi-1)/(chi+1)
                aij = uij/(uij**2-1)
                if (aij >  0.5d0) aij =  0.5d0
                if (aij < -0.5d0) aij = -0.5d0
                rmiu2 = rmiu+aij*(1-rmiu**2)
                if (rmiu2.le.cutoffa) then
                  tmps=-1.0d+00
                elseif (rmiu2.ge.cutoffa) then
                  tmps=+1.0d+00
                else
                  rmiu2n=rmiu2/cutoffa
                  tmps = 1.5d0*rmiu2n-0.5d0*rmiu2n**3
                  do ib = 2,pow
                    tmps = 1.5d0*(tmps)-0.5d0*(tmps)**3
                  end do
                endif
                smat(j,k) = 0.5d0*(1.0d+00-tmps)
              end do
            end do
            pvec = 1.0d+00
            do ii = 1,ncent
              pvec = pvec*smat(:,ii)
            end do
            beckew = pvec(atom)/sum(pvec)
            kk = kk + 1
            wthiscenter = wrads(ir) * wang(il)
            weight = beckew * wthiscenter
c
c-----------Compute MOs a the mesh point x.
c
            call calcmos (x,c,xmo)
c
c-----------AOM integrals
c
            do m=1,nmo
              do n=1,m
                pmos=xmo(m)*xmo(n)
                pabs=abs(pmos)
                aom(atom,m,n)  = aom(atom,m,n)  + pmos * weight
                saom(atom,m,n) = saom(atom,m,n) + pabs * weight
              enddo
            enddo
          end do
        end do
      end do
c
      if (allocated(rads)) then
        deallocate (rads,stat=istat)
        if (istat /= 0) then
          stop 'aombecke.f: Cannot deallocate rads() array'
        endif
      endif
      if (allocated(wrads)) then
        deallocate (wrads,stat=istat)
        if (istat /= 0) then
          stop 'aombecke.f: Cannot deallocate wrads() array'
        endif
      endif
      if (allocated(xang)) then
        deallocate (xang,stat=istat)
        if (istat /= 0) then
          stop 'aombecke.f: Cannot deallocate xang() array'
        endif
      endif
      if (allocated(yang)) then
        deallocate (yang,stat=istat)
        if (istat /= 0) then
          stop 'aombecke.f: Cannot deallocate yang() array'
        endif
      endif
      if (allocated(zang)) then
        deallocate (zang,stat=istat)
        if (istat /= 0) then
          stop 'aombecke.f: Cannot deallocate zang() array'
        endif
      endif
      if (allocated(wang)) then
        deallocate (wang,stat=istat)
        if (istat /= 0) then
          stop 'aombecke.f: Cannot deallocate wang() array'
        endif
      endif
c
c-----Symmetrize AOM
c
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
c.....Write results to an AOM file called wfnf'.beckeAOM'
c
      beckefile=trim(wfnf)//'.beckeAOM'
      open (18,file=beckefile,iostat=ierr)
      if (ierr.ne.0) then
        write (0,*) 'aombecke.f: Error opening *.beckeAOM file'
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
      if (verbose) then
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
 80   format (6(1x,e16.10))
 10   format (6(1x,F15.8))
 501  format (1x,69('-'))

      call timer (4,imeshb,'_aombecke_',-1)
      end 
c
c-----------------------------------------------------------------------
c
c
c.....Malla Handy, Eulec-Maclaurin n=2
c
      subroutine rmeshh2(n,rmid,r,wintr)
      implicit none
!
!     radial mesh and integration weights,
!
!     the q-mesh is uniform on the interval (0,+1). transformation is
!
!                    q**2
!      r  =  rmid ---------
!                 ( 1 - q )**2
!
      real(kind=8), parameter :: pi = 3.141592653589793d0
      real(kind=8), parameter :: fourpi = 4d0*pi
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(in) :: rmid
      real(kind=8), intent(out), dimension(n) :: r, wintr
     
      real(kind=8) :: h, q
      integer(kind=4) :: i
     
      h = 1.0d+00/real(n+1,8)
      do i = 1,n
        q = h*i
        r(i) = q**2/(1-q)**2*rmid
        wintr(i) = 2*q**5*h/(1-q)**7*rmid**3*fourpi
      end do
      return
      end
c
c-----------------------------------------------------------------------
c
c
c.....Malla radial Euler-McLaurin with n = 1
c
      subroutine rmeshh1(n,rmid,r,wintr)
    
      implicit none
!
!     radial mesh and integration weights,
!
!     the q-mesh is uniform on the interval (0,+1). transformation is
!
!                   q
!     r  =  rmid ---------
!                ( 1 - q )
!
      real(kind=8), parameter :: pi = 3.141592653589793d0
      real(kind=8), parameter :: fourpi = 4d0*pi
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(in) :: rmid
      real(kind=8), intent(out), dimension(n) :: r, wintr
    
      real(kind=8) :: h, q
      integer(kind=4) :: i
    
      h = 1.0d+00/real(n+1,8)
      do i = 1,n
        q = h*i
        r(i) = rmid * q / (1.d0-q)
        wintr(i) = fourpi * h * r(i)**2 * rmid / (1.d0-q)**2
      end do
      return
      end
c
c-----------------------------------------------------------------------
c
c
c-----Malla Becke radial
c
      subroutine rmeshb(n,rmid,r,wintr)
    
      implicit none
!
!     radial mesh and integration weights,
!
!     the q-mesh is uniform on the interval (-1,+1). transformation is
!
!                 ( 1 + q )
!      r  =  rmid ---------
!                 ( 1 - q )
!
      real(kind=8), parameter :: pi = 3.141592653589793d0
      real(kind=8), parameter :: fourpi = 4d0*pi
      integer(kind=4), intent(in) :: n
      real(kind=8), intent(in) :: rmid
      real(kind=8), intent(out), dimension(n) :: r, wintr
    
      real(kind=8) :: h, q
      integer(kind=4) :: i
    
      h = 1.0d+00/real(n+1,8)
      do i = 1,n
        q = cos(h*i*pi)
        r(i) = ((1+q)/(1-q))*rmid
        wintr(i) = 2*pi*h*rmid**3*(1+q)**2.5d0/(1-q)**3.5d0*fourpi
      end do
      return
      end
