C
c----------------------------------------------------------------------
c
      subroutine molorb1 (p,xmo,xmo1,canmos)
c
c----------------------------------------------------------------------
c     Compute the molecular orbitals (xmo) and their gradients (xmo1)
c     at the point p(). If canmo=.true. the MOs are the canonical ones. 
c     Otherwise, they are the natural MOs.
c----------------------------------------------------------------------
c
      USE        space_for_wfnbasis
      USE        space_for_wfncoef
      include   'implicit.inc'
      include   'param.inc'    
      include   'constants.inc'    
      include   'wfn.inc'    

      dimension  it(3)
      dimension  p(3),xcoor(3),fun(3),fun1(3)
      real(kind=8) xmo(nmo),xmo1(nmo,3)
      logical canmos
c
      call timer (2,imolorb1,'_molorb1  ',-1)
      xmo     =  zero
      xmo1    =  zero
      nmoplus = izero
      if (canmos) nmoplus = nmo
c
c     run over primitives
c
      do i=1,nprims
         ic=icen(i)
         itip=ityp(i)
         it(1:3)=nlm(itip,1:3)
         ori=-oexp(i)
         dp2=ori+ori
         xcoor(1:3)=p(1:3)-xyz(ic,1:3)
         dis2=xcoor(1)*xcoor(1)+xcoor(2)*xcoor(2)+xcoor(3)*xcoor(3)
         aexp=exp(ori*dis2)
         do j=1,3
            n=it(j)
            x=xcoor(j)
            x2=x*x
            dp2x2=dp2*x2
            if (n.ge.4) then
               pow=x**(n-3)
               nn1=n*(n-1)
               powx=pow*x
               fun1(j)=powx*x*(n+dp2x2)
               fun (j)=powx*x2
            else if (n.eq.3) then
               fun1(j)=x2*(3d0+dp2x2)
               fun (j)=x*x2
            else if (n.eq.2) then
               fun1(j)=x*(2d0+dp2x2)
               fun (j)=x2
            else if (n.eq.1) then
               fun1(j)=1d0+dp2x2
               fun (j)=x
            else if (n.eq.0) then
               dp2x=dp2*x
               fun1(j)=dp2x
               fun (j)=1d0
            endif
         enddo
         f23 = fun (2)*fun(3)*aexp  
         f13 = fun (1)*fun(3)*aexp
         f12 = fun (1)*fun(2)*aexp
c
c....... run over orbitals.
c
         do j=1,nmo
            cfj=coef(j+nmoplus,i)
            xmo(j)=xmo(j)+cfj*(fun(1)*f23)
            xmo1(j,1)=xmo1(j,1)+cfj*(fun1(1)*f23)  ! x
            xmo1(j,2)=xmo1(j,2)+cfj*(fun1(2)*f13)  ! y
            xmo1(j,3)=xmo1(j,3)+cfj*(fun1(3)*f12)  ! z
         enddo 
      enddo
      call timer (4,imolorb1,'_molorb1  ',-1)
      end
