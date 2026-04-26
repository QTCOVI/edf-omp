c
c-----------------------------------------------------------------------
c
      subroutine wbec (atom,power,p,weight,exact,
     &   salvadornu,scuseriaw)
c
c-----------------------------------------------------------------------
c
      USE         space_for_wfnbasis
      implicit real(kind=8) (a-h,o-z)
      include    'param.inc'
      include    'wfn.inc'
      include    'datatm.inc'
      parameter  (EPS=1.0d-10)
      real(kind=8),     allocatable,dimension (:)  :: pb
      real(kind=8)      p(3)
      real(kind=8)      xio(3),xjo(3),xij(3)
      real(kind=8)      muij,nuij
      integer     power,atom
      logical     exact,salvadornu,scuseriaw
c
c-----------------------------------------------------------------------
c
c     Stratmann, Scuseria, Frisch weight in Becke-like integrations used
c
      if (scuseriaw) then
        call wscuseria (atom,p,weight)
        return
      endif
c
      call timer (2,ivect,'_wbec_____',-1)
c
      allocate (pb(ncent),stat=ier)
      if (ier.ne.0) then
        write (0,*) 'wbec.f: Cannot allocate pb()'
        return
      endif
c
      sumpb = 0.0d+00
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
          if (j.eq.i) cycle
c
c---------Distance from p() to atom j.
c
          xjo(:)=xyz(j,:)-p(:)
          rj=sqrt(dot_product(xjo,xjo))
c
c---------Rij distance.
c
          xij(:)=xyz(i,:)-xyz(j,:)
          rij=sqrt(dot_product(xij,xij))
c
c---------nu_ij parameter.
c
          muij=(ri-rj)/rij
          xi = rbragg(nint(charge(i)))/rbragg(nint(charge(j)))
c
c---------Modification of Salvador-and-Ramos_Cordoba.
c
          if (salvadornu) then
            xionemu = xi*(1.0d+00-muij)
            nuij = (1.0d+00+muij-xionemu)/(1.0d+00+muij+xionemu)
          else
c
c---------Original 'Nu' by Axel Becke
c
            uij=(xi-1.0d+00)/(xi+1.0d+00)
            aij=uij/(uij*uij-1.0d+00)
            if (aij.ge.(+0.5d0))    aij=+0.5d0
            if (aij.le.(-0.5d0))    aij=-0.5d0
            nuij=muij+aij*(1.0d+00-muij*muij)
          endif
c
c---------h{h[...h(nu_ij)]} function.
c
          h=1.5d0*nuij-0.5d0*nuij**3
          do m=2,power
            h=1.5d0*h-0.5d0*h**3
          enddo
c
c---------S_k function.
c
          s=0.5d0*(1.0d+00-h)
c
c---------Finally, the P_i function is computed.
c
          pb(i)=pb(i)*s
          if (.not.exact)  then
            if (pb(i).lt.EPS .and. i.eq.atom) then
              weight=0.0d+00
              goto 1
            endif
          endif
        enddo
        sumpb=sumpb+pb(i)
      enddo
      if (abs(sumpb).eq.0.0d+00) then
         weight=0.0d+00
      else
         weight=pb(atom)/sumpb
      endif
 1    continue
      call timer (4,ivect,'_wbec_____',-1)
      return
      end
