      subroutine mindefw (p,icent,wmindef,irho,coef,xmo)
c
c----------------------------------------------------------------------
c     Computes at the point p() the canonical MOs defined in terms of 
c     CGFs by the array coef(). The results are returned in xmo().
c     It also computes the 'w_A=rho_A/rho' weight of the MinDef atom
c     number 'icent', which is returned in wmindef
c----------------------------------------------------------------------
c
      USE          space_for_primgto
      USE          space_for_wfnbasis
      include     'implicit.inc'
      include     'param.inc'    
      include     'wfn.inc'    
      include     'primgto.inc'    

      real(kind=8) coef(nmo+nmo,nprims)
      real(kind=8) prim(nprims),dena(ncent)
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
            prim(i)=fun(1)*fun(2)*fun(3)*aexp
            do j=1,nmo
              xmo(j)=xmo(j)+coef(j+nmo,i)*f123
            enddo
          enddo
        enddo
      enddo
      if (irho.eq.4) return
c
c-----Weight factor of atom 'icent'
c
      dena=0d0
      do i=1,nprims
        ic=icen(i)
        ori=oexp(i)
        do j=1,i
          jc=icen(j)
          orj=oexp(j)
          pmat=0d0
          do k=1,nmo
            pmat=pmat+coef(k,i)*coef(k,j)*occ(k)
          enddo
          if (i.eq.j) pmat=pmat+pmat
          term=prim(i)*prim(j)*pmat
          if (ic.eq.jc) then
            dena(ic)=dena(ic)+term        ! Net term
          else
            if (ori.gt.orj) then
              dena(ic)=dena(ic)+term
            else if (ori.lt.orj) then
              dena(jc)=dena(jc)+term
            else
              dena(ic)=dena(ic)+term*0.5D0
              dena(jc)=dena(jc)+term*0.5D0
            endif
          endif
        enddo
      enddo
      densidad=0d0
      do ic=1,ncent
        densidad=densidad+dena(ic)
      enddo
      if (abs(densidad).gt.0d0.and.dena(icent).gt.0d0) then
        wmindef=dena(icent)/densidad
      else
        wmindef=0d0
      endif
      return
      end
