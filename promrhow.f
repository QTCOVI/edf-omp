      subroutine promrhow (p,icent,rhopow,wnetrho,coef,xmo)
c
c----------------------------------------------------------------------
c     Computes at the point p() the canonical MOs defined in terms of 
c     CGFs by the array coef(). The results are returned in xmo().
c     It also computes the 'w_A=rho_A/rho' weight of the atom number 
c     'icent', which is returned in wnetrho
c----------------------------------------------------------------------
c
      USE          space_for_primgto
      USE          space_for_wfnbasis
      include     'implicit.inc'
      include     'param.inc'    
      include     'wfn.inc'    
      include     'primgto.inc'    

      real(kind=8) coef(nmo+nmo,nprims)
      real(kind=8) prim(nprims),rhoat(ncent)
      real(kind=8) xmo(nmo)
      dimension    it(3),p(3),xcoor(3),fun(3)
c
      xmo(1:nmo)=0d0
      do ic=1,ncent
        xcoor(:)=p(:)-xyz(ic,:)
        dis2=xcoor(1)*xcoor(1)+xcoor(2)*xcoor(2)+xcoor(3)*xcoor(3)
        do m=1,ngroup(ic)
          k=nuexp(ic,m,1)
          ori=-oexp(k)
          aexp=exp(ori*dis2)
          do jj=1,nzexp(ic,m)
            i=nuexp(ic,m,jj)
            itip=ityp(i)
            it(:)=nlm(itip,:)
            do j=1,3
               n =it(j)
               x=xcoor(j)
               if (n.ge.4) then
                  x2=x*x
                  pow=x**(n-2)
                  fun(j)=pow*x2
               else if  (n.eq.3) then
                  x2=x*x
                  fun (j)=x*x2
               else if  (n.eq.2) then
                  x2=x*x
                  fun(j)=x2
               else if  (n.eq.1) then
                  fun (j)=x
               else if  (n.eq.0) then
                  fun(j)=1d0
               endif
            enddo
            f123=fun(1)*fun(2)*fun(3)*aexp
            prim(i)=f123
            do j=1,nmo
              xmo(j)=xmo(j)+coef(j+nmo,i)*f123
            enddo
          enddo
        enddo
      enddo
c
c-----Weight factor of atom 'icent'
c
      rho=0d0
      rhoat=0d0
      do ic=1,ncent
        xcoor(:)=p(:)-xyz(ic,:)
        dis2=xcoor(1)*xcoor(1)+xcoor(2)*xcoor(2)+xcoor(3)*xcoor(3)
        rhoat(ic)=atomicrho(dis2,int(charge(ic)))
        rhoat(ic)=rhoat(ic)**rhopow
        rho=rho+rhoat(ic)
      enddo
      if (abs(rho).gt.0d0.and.rhoat(icent).gt.0d0) then
        wnetrho=rhoat(icent)/rho
      else
        wnetrho=0d0
      endif
      return
      end
