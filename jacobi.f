c
c-----------------------------------------------------------------------
c
      subroutine jacobi (a, n, np, d, v, nrot)
c
c.....jacobi - resolution of an eigenvalue equation.
c
      implicit real(kind=8) (a-h,o-z)
c
c.....Constants used in the calculation.
c
      parameter  (zero       = 0.0d+00)
      parameter  (one        = 1.0d+00)
      parameter  (half       = 0.5d+00)
      parameter  (afifth     = 0.2d+00)
c
      real(kind=8),  allocatable,dimension(:) :: b,z
c
c-----ax() is an auxiliary array that will be used to store the input
c     matrix a(). Before doing nothing we do ax(1:n,1:n) = a(1:n,1:n).
c     At the end of the routine we do a(1:n,1:n) = ax(1:n,1:n)
c
c-----In this way, the input matrix a() is not destroyed, but returned
c     as its original value.
c

      real(kind=8),  allocatable,dimension(:,:) :: ax
      real(kind=8)   a(np,np), d(np), v(np,np)
c
      allocate (b(n),stat=ier)
      if (ier.ne.0) stop 'jacobi.f: Cannot allocate b()'
      allocate (z(n),stat=ier)
      if (ier.ne.0) stop 'jacobi.f: Cannot allocate z()'
c
c-----Added on Nov 17, 2022 to avoid that the input matrix is destroyed
c
      allocate (ax(n,n),stat=ier)
      if (ier.ne.0) stop 'jacobi.f: Cannot allocate ax()'
      ax(1:n,1:n) = a(1:n,1:n)
c
      do 12 ip = 1, n
        do 11 iq = 1, n
          v(ip,iq) = zero
 11     continue
        v(ip,ip) = one
 12   continue
      do 13 ip = 1, n
        b(ip) = a(ip,ip)
        d(ip) = b(ip)
        z(ip) = zero
 13   continue
      nrot = 0
      do 24 i = 1, 500
        sm = zero
        do 15 ip = 1, n-1
          do 14 iq = ip+1, n
            sm = sm + abs(a(ip,iq))
 14       continue
 15     continue
        if (sm.eq.zero) then
          if (allocated(b)) then
            deallocate (b,stat=ier)
            if (ier.ne.0) stop 'jacobi.f: Cannot deallocate b()'
          endif
          if (allocated(z)) then
            deallocate (z,stat=ier)
            if (ier.ne.0) stop 'jacobi.f: Cannot deallocate z()'
          endif
c
c---------Recover the input matrix that was diagonalized
c
          a(1:n,1:n) = ax(1:n,1:n)
          return
        endif
        if (i.lt.4) then
          tresh = afifth*sm/n**2
        else
          tresh = zero
        endif
        do 22 ip = 1, n-1
          do 21 iq = ip+1, n
            g = 100*abs(a(ip,iq))
            if ((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))
     &         .and.(abs(d(iq))+g.eq.abs(d(iq)))) then
              a(ip,iq) = zero
            else if (abs(a(ip,iq)).gt.tresh) then
              h = d(iq)-d(ip)
              if (abs(h)+g.eq.abs(h)) then
                t = a(ip,iq)/h
              else
                theta = half*h/a(ip,iq)
                t = one/(abs(theta)+sqrt(one+theta**2))
                if (theta.lt.zero) t = -t
              endif
              c = one/sqrt(1+t**2)
              s = t*c
              tau = s/(one+c)
              h = t*a(ip,iq)
              z(ip) = z(ip)-h
              z(iq) = z(iq)+h
              d(ip) = d(ip)-h
              d(iq) = d(iq)+h
              a(ip,iq) = zero
              do 16 j = 1, ip-1
                g = a(j,ip)
                h = a(j,iq)
                a(j,ip) = g-s*(h+g*tau)
                a(j,iq) = h+s*(g-h*tau)
 16           continue
              do 17 j = ip+1, iq-1
                g = a(ip,j)
                h = a(j,iq)
                a(ip,j) = g-s*(h+g*tau)
                a(j,iq) = h+s*(g-h*tau)
 17           continue
              do 18 j = iq+1, n
                g = a(ip,j)
                h = a(iq,j)
                a(ip,j) = g-s*(h+g*tau)
                a(iq,j) = h+s*(g-h*tau)
 18           continue
              do 19 j = 1, n
                g = v(j,ip)
                h = v(j,iq)
                v(j,ip) = g-s*(h+g*tau)
                v(j,iq) = h+s*(g-h*tau)
 19           continue
              nrot = nrot+ 1
            endif
 21       continue
 22     continue
        do 23 ip = 1, n
          b(ip) = b(ip)+z(ip)
          d(ip) = b(ip)
          z(ip) = zero
 23     continue
 24   continue
      stop 'jacobi.f: 500 iterations should never happen'
      return
      end
