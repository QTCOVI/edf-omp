      subroutine calcmos (p,coef,xmo)
c
c----------------------------------------------------------------------
c     Computes at the point p() the MOs defined in terms of CGFs by the
c     array coef(). The results are returned in xmo().
c----------------------------------------------------------------------
c
      USE          space_for_primgto
      USE          space_for_wfnbasis
      include     'implicit.inc'
      include     'param.inc'    
      include     'wfn.inc'    
      include     'primgto.inc'    

      real(kind=8) coef(nmo,nprims)
      real(kind=8) xmo(nmo)
      dimension    it(3),p(3),xcoor(3),fun(3)
c
      xmo(1:nmo)=0d0
      do ic=1,ncent
        xcoor(:)=p(:)-xyz(ic,:)
        dis2=dot_product(xcoor(:),xcoor(:))
        do m=1,ngroup(ic)
          k=nuexp(ic,m,1)
          ori=-oexp(k)
          dp2=ori+ori
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
            do j=1,nmo
              xmo(j)=xmo(j)+coef(j,i)*f123
            enddo
          enddo
        enddo
      enddo
      return
      end
