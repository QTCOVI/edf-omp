C
c----------------------------------------------------------------------
c
      subroutine molorb0 (p,xmo,canmos)
c
c----------------------------------------------------------------------
c     Compute the molecular orbitals (xmo) at the point p(). 
c     If canmo=.true. the MOs are the canonical ones. Otherwise, they 
c     are the natural MOs.
c----------------------------------------------------------------------
c
      USE        space_for_wfnbasis
      USE        space_for_wfncoef
      include   'implicit.inc'
      include   'param.inc'    
      include   'constants.inc'    
      include   'wfn.inc'    

      dimension  it(3)
      dimension  p(3),xcoor(3),fun(3)
      real(kind=8) xmo(nmo)
      logical canmos
c
      call timer (2,imolorb0,'_molorb0  ',-1)
      xmo     =  zero
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
         xcoor(1:3)=p(1:3)-xyz(ic,1:3)
         dis2=xcoor(1)*xcoor(1)+xcoor(2)*xcoor(2)+xcoor(3)*xcoor(3)
         aexp=exp(ori*dis2)
         do j=1,3
            n=it(j)
            x=xcoor(j)
            x2=x*x
            if (n.ge.4) then
               fun (j)=x**n
            else if (n.eq.3) then
               fun (j)=x*x2
            else if (n.eq.2) then
               fun (j)=x2
            else if (n.eq.1) then
               fun (j)=x
            else if (n.eq.0) then
               fun (j)=1d0
            endif
         enddo
         f123 = fun(1)*fun(2)*fun(3)*aexp
c
c....... run over orbitals.
c
         do j=1,nmo
            cfj=coef(j+nmoplus,i)
            xmo(j)=xmo(j)+cfj*f123
         enddo 
      enddo
      call timer (4,imolorb0,'_molorb0  ',-1)
      end
