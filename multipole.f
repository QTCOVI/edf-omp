c
c----------------------------------------------------------------------
c
      subroutine multipole 
     &  (pnt,nmo,nprims,ncent,maxtype,nlm,nexp,ipow,coef,rexp)
c
c.....Computes different orbital multipoles with respect to the point 
c     pnt(), i.e.
c
c                    /
c                   |      
c     multip(i) = - | Phi_i Phi_i x_C^nxj y_C^nyj z_C^nzj dx dy dz
c                   |
c                   /
c
c     where i=1,...,NMO, 
c     [nxj,nyj,nzj] = [ipow(j,1),ipow(j,w),ipow(j,3)], j=1,...,nexp
c
c     and
c
c     x_C = x - pnt(1)
c     y_C = y - pnt(2)
c     z_C = z - pnt(3)
c
c     The integrals are computed using the expressions of Section 9.5.3
c     of "Molecular Electronic-Structure Theory" by Trygve Helgaker,
c     Poul Jorgensen, and Jeppe Olsen (HJJ).
c
c----------------------------------------------------------------------
c
      USE        space_for_wfnbasis
      USE        space_for_primgto
      include   'implicit.inc'
      include   'primgto.inc'
      real(kind=8), parameter :: epsden = 1.0D-10
      real(kind=8), parameter :: zero   = 0D0
      real(kind=8), parameter :: pi     = 3.14159265358979324D+00
c
c.....etijx() are the coefficients, except for factor EXP(-XMU*R_AB^2), 
c     where XMU=a*b/(a+b), that result from the expansion of the product 
c     of two primitive cartesian Gaussian
c
      real(kind=8), allocatable,dimension (:,:,:,:)  :: etijx
      real(kind=8), allocatable,dimension   (:,:,:)  :: mcoef
      real(kind=8), allocatable,dimension     (:,:)  :: sije

      real(kind=8)    coef(nmo,nprims)
      real(kind=8)    rexp(nmo,nexp)
      integer(kind=4) ipow(nexp,3)
c
      real(kind=8) rpc(3),pnt(3),ax(3),bx(3)
      integer(kind=4) nlm(maxtype,3)
      integer(kind=4) t,u,v,powe
c
c-----------------------------------------------------------------------
c
      mx=0
      do n=1,nexp
        do j=1,3
          if (ipow(n,j).gt.mx) mx=ipow(n,j)
          if (ipow(n,j).lt.0)  stop 'multipole.f: Illegal power of x'
        enddo
      enddo
      if (.not.allocated(etijx)) then
        allocate (etijx(-1:2*mx,-1:mx,-1:mx,3),stat=ier)
        if (ier.ne.0) stop 'multipole.f: Cannot allocate etijx()'
      endif
      if (.not.allocated(mcoef)) then
        allocate (mcoef(-1:2*mx+1,0:2*mx,3),stat=ier)
        if (ier.ne.0) stop 'multipole.f: Cannot allocate mcoef()'
      endif
      if (.not.allocated(sije)) then
        allocate (sije(0:mx,3),stat=ier)
        if (ier.ne.0) stop 'multipole.f: Cannot allocate sije()'
      endif

      mcoef(-1:2*mx+1,0:2*mx,1:3) = 0.0d+00
      rexp(1:nmo,1:nexp)          = 0.0d+00

      do ica=1,ncent
        ax(:)=xyz(ica,:)
        do ma=1,ngroup(ica)   
          nua=nuexp(ica,ma,1)
          za=oexp(nua)
          itipa=ityp(nua)
          la=nlm(itipa,1)+nlm(itipa,2)+nlm(itipa,3)
          do icb=1,ica
            bx(:)=xyz(icb,:)
            ab2 = 0.0d+00   ! Square of the intercenter A-B distance
            do j=1,3
              ab2 = ab2+(ax(j)-bx(j))**2
            enddo
            if (icb.lt.ica) then
              ngrpdo=ngroup(icb)
            else
              ngrpdo=ma
            endif
            do 3000 mb=1,ngrpdo
              nub=nuexp(icb,mb,1)
              zb=oexp(nub)
              itipb=ityp(nub)
              lb=nlm(itipb,1)+nlm(itipb,2)+nlm(itipb,3)
c
c.............p=exponent of this product of Gaussian Shells.
c
              p=za+zb
              xmu=za*zb/p
              prefactor=exp(-xmu*ab2)
              do j=1,3
                rpc(j) = ((za*ax(j)+zb*bx(j))/p)-pnt(j)
              enddo
              etijx=0.0d+00
              do j=1,3
                call etijcalc (j,mx,la,lb,etijx,za,zb,ax(j),bx(j))
              enddo
c
c.............Recursive evaluation of M_t^(e+1) coeffs, Eq. 9.5.36 HJJ.
c
              mcoef(0,0,1:3)=sqrt(pi/p)             
              do m=1,3
                xpc=rpc(m)
                do ix=0,mx
                  powe=ix
                  if (powe.gt.0) then
                    mcoef(0,1,m)=xpc*mcoef(0,0,m)
                    mcoef(1,1,m)=mcoef(0,0,m)                   
                  endif
                  if (powe.gt.1) then
                    do ie=1,powe-1
                      do it=0,ie+1
                        c1=mcoef(it-1,ie,m)
                        c2=mcoef(it,ie,m)
                        c3=mcoef(it+1,ie,m)
                        mcoef(it,ie+1,m)=dble(it)*c1+xpc*c2+c3/(2*p)
                      enddo
                    enddo
                  endif
                enddo
              enddo
c
c.............Compute S_ij^e integrals for all pairs of primitives,
c             Eq. 9.5.29 HJJ.
c
              do ka=1,nzexp(ica,ma)
                nua=nuexp(ica,ma,ka)
                itipa=ityp(nua)
                do 1000 kb=1,nzexp(icb,mb)
                  nub=nuexp(icb,mb,kb)
                  if (nua.lt.nub) goto 1000
                  itipb=ityp(nub)

                  do m=1,3
                    do ie=0,mx
                      sx=0.0d+00
                      do t=0,min(nlm(itipa,m)+nlm(itipb,m),ie)
                        nlma=nlm(itipa,m)
                        nlmb=nlm(itipb,m)
                        sx=sx+etijx(t,nlma,nlmb,m)*mcoef(t,ie,m)
                      enddo
                      sije(ie,m)=sx
                    enddo
                  enddo

                  do n=1,nexp
                    i=ipow(n,1)
                    j=ipow(n,2)
                    k=ipow(n,3)
                    prods=sije(i,1)*sije(j,2)*sije(k,3)
                    do l=1,nmo
                      cprod=coef(l,nua)*coef(l,nub)
                      if (nua.ne.nub) cprod=cprod+cprod
                      premat=cprod*prefactor
                      rexp(l,n)=rexp(l,n)+premat*prods
                    enddo
                  enddo
 1000           continue
              enddo
 3000       continue
          enddo
        enddo
      enddo
      if (allocated(etijx)) then
        deallocate (etijx,stat=ier)
        if (ier.ne.0) stop 'multipole.f: Cannot deallocate etijx()'
      endif
      if (allocated(mcoef)) then
        deallocate (mcoef,stat=ier)
        if (ier.ne.0) stop 'multipole.f: Cannot deallocate mcoef()'
      endif
      if (allocated(sije)) then
        deallocate (sije,stat=ier)
        if (ier.ne.0) stop 'multipole.f: Cannot deallocate sije()'
      endif
      return
      end
