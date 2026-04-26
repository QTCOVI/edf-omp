c
c-----------------------------------------------------------------------
c
      subroutine wscuseria (atom,p,weight)
c
c-----------------------------------------------------------------------
c
      USE         space_for_wfnbasis
*     USE         space_for_bicen
      implicit real(kind=8) (a-h,o-z)
      include    'param.inc'
      include    'wfn.inc'
      include    'datatm.inc'
*     include    'integ.inc'
      real(kind=8), parameter :: aparam = 0.64d0
      real(kind=8),     allocatable,dimension (:)  :: pb
      real(kind=8)      p(3)
      real(kind=8)      xio(3),xjo(3),xij(3)
      real(kind=8)      muij,nuij
      real(kind=8)      numod1,numod3,numod5,numod7
      integer           atom
c
c-----------------------------------------------------------------------
c
      call timer (2,iwscu,'_wscuseria',-1)
c
      allocate (pb(ncent),stat=ier)
      if (ier.ne.0) then
        write (0,*) 'wbec.f: Cannot allocate pb()'
        return
      endif
c
      sumpb=0.0d+00
c
c-----Compute P_i factor for all the atoms.
c
      do i=1,ncent
        pb(i)=1.0d+00
c
c-------Distance from p() to atom i.
c
        xio(:)=xyz(i,:)-p(:)
        ri=sqrt(dot_product(xio,xio))
        do j=1,ncent
          if (j.ne.i) then
c
c-----------Distance from p() to atom j.
c
            xjo(:)=xyz(j,:)-p(:)
            rj=sqrt(dot_product(xjo,xjo))
c
c-----------Rij distance.
c
            xij(:)=xyz(i,:)-xyz(j,:)
            rij=sqrt(dot_product(xij,xij))
c
c-----------nu_ij parameter.
c
            muij=(ri-rj)/rij
            xi = rbragg(int(charge(i)))/rbragg(int(charge(j)))
c
c-----------Original 'Nu' by Axel Becke
c
            uij=(xi-1.0d+00)/(xi+1.0d+00)
            aij=uij/(uij*uij-1.0d+00)
            if (aij.ge.(+0.5d0))    aij=+0.5d0
            if (aij.le.(-0.5d0))    aij=-0.5d0
            nuij=muij+aij*(1.0d+00-muij*muij)
c
c-----------Here the recipe by Stratmann, Scuseria, and Frisch
c
            if (nuij.le.-aparam) then
              h=-1.0d0
            elseif (nuij.ge.+aparam) then
              h=+1.0d0
            else
              numod1=nuij/aparam
              numod3=numod1**3
              numod5=numod1**5
              numod7=numod1**7
              h=(35.0d0*numod1-35.0d0*numod3+21.0d0*numod5-
     &            5.0d0*numod7)/16.0d0
            endif
c
c-----------S_k function.
c
            s=0.5d0*(1d0-h)
c
c-----------Finally, the P_i function is computed.
c
            pb(i)=pb(i)*s
          endif
        enddo
        sumpb=sumpb+pb(i)
      enddo
      if (abs(sumpb).eq.0.0d+00) then
         weight=0.0d+00
      else
         weight=pb(atom)/sumpb
      endif
      continue
      call timer (4,iwscu,'_wscuseria',-1)
      return
      end
