c----------------------------------------------------------------------
      subroutine       orbfun(p,gun,lcor)
c----------------------------------------------------------------------
c     computes the coulombic auxiliary functions.
c     
c
c        returns orbitals.
c........lcor : 0 standard (possibly natural)
c........lcor : 1 canonical 
c
c----------------------------------------------------------------------
c 
      USE               space_for_wfnbasis
      USE               space_for_wfncoef
      include          'implicit.inc'
      include          'param.inc'    
      include          'wfn.inc'    

c     local vars
      parameter        (half=0.5d0)
      dimension        it(3), jt(3)
c
      dimension        p(3), xcoor(3), fun (3),  gun(*)
c
c
c.....correlation? 
      if (lcor.eq.0) ifirst=0
      if (lcor.ne.0) ifirst=nmo
c
c
      do i=1,nmo
         gun(i)=0d0
      enddo 
    
         
c
c     run over primitives
      do i=1, nprims
         
c        compute common elements
         ic=icen(i)
         itip=ityp(i)
c
c........integer coefficients
c---------------------------------
         it(1)=nlm(itip,1)
         it(2)=nlm(itip,2)
         it(3)=nlm(itip,3)
c---------------------------------
c
         ori=-oexp(i)
         dp2=ori+ori
c
c........atomic coordinates
         xcoor(1)=p(1)-xyz(ic,1)
         xcoor(2)=p(2)-xyz(ic,2)
         xcoor(3)=p(3)-xyz(ic,3)
c
         dis2=xcoor(1)*xcoor(1)+xcoor(2)*xcoor(2)+xcoor(3)*xcoor(3)
         aexp=exp(ori*dis2)
c
         do j=1,3
            n =it(j)
            x=xcoor(j)
c
            if       (n.ge.4) then
               x2=x*x
               pow=x**(n-2)
               fun (j)=pow*x2
            else if  (n.eq.3) then
               x2=x*x
               fun (j)=x*x2
            else if  (n.eq.2) then
               x2=x*x
               fun (j)=x2
            else if  (n.eq.1) then
               fun (j)=x
            else if  (n.eq.0) then
               fun (j)=1d0
            endif
c
         enddo
c
         f12=fun (1)*fun (2)*aexp
         f123=f12*fun (3)
c
c....... run over orbitals.
         do j=1,nmo
            cfj=coef(j+ifirst,i)
            gun (j)=gun(j)+cfj*f123
         enddo      ! of orbitals
      enddo
c
c
      end
