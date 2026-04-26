c
c-----------------------------------------------------------------------
c
      subroutine etijcalc (m,lamx,la,lb,ce,a,b,ax,bx)
c
c-----This routine computes recursively the E_t^ij coefficients
c     resulting from the expansion of a cartesian gaussian product in 
c     terms of hermite gaussians. The expressions used come from the 
c     book "Molecular Electronic-Structure Theory" by T. Helgaker, 
c     P.Jorgensen, and J. Olsen. All the coefficients arise from the 
c     first one [E_0^00], taking also into account that E_t^ij=0 for 
c     t<0 or t>i+j.
c
c-----Details
c     
c     G_i (x,a,A_x) = x_A^i exp(-a x_A^2)  ;  x_A = x -  A_x
c     G_j (x,b,B_x) = x_B^j exp(-b x_B^2)  ;  x_B = x -  B_x
c
c     where
c     x   = X coordinate of an electron in a general reference frame
c     a   = Exponente of the first  Cartesian Gaussian
c     b   = Exponente of the second Cartesian Gaussian
c     A_X = X coordinate of nucleus A in a general reference frame
c     B_X = X coordinate of nucleus B in a general reference frame
c
c     Omega_ij      = G_i (x,a,A_x) * G_j (x,b,B_x)
c
c     Omega_ij      = SUM [t=0,i+j] K_ab^x E_t^ij Lambda_t (x,p,P_x)
c
c     where
c
c     Lambda_t (x,p,P_x) = [d/d P_x]^t exp(-p x_P^2)  ;  x_P = x -  P_x
c     p                  = a+b
c     P_x                = (a A_x + b B_x)/p
c     mu                 = ab/(a+b)
c     X_AB               = A_x - B_x
c     K_ab^x             = exp(-mu X_AB^2)
c
c
c     The gives E_t^ij for  (0 <= i <= la; 0 <= j <= lb; 0 <= t <= i+j)
c
c-----------------------------------------------------------------------
c
      implicit none
      integer(kind=4) la,lb,lab,i,j,t,ij,i1,j1,t1,m,lamx
      real(kind=8) ce(-1:2*lamx,-1:lamx,-1:lamx,3)
      real(kind=8) a,b,ax,bx,p,ab,pa,pb,tp,zij,oij
c
      if (la.lt.0)    stop 'etijcalc.f: Fatal error, la < 0 '
      if (lb.lt.0)    stop 'etijcalc.f: Fatal error, lb < 0 '
      if (la.gt.lamx) stop 'etijcalc.f: Fatal error, la > lamx '
      if (lb.gt.lamx) stop 'etijcalc.f: Fatal error, lb > lamx '
      lab=la+lb
      ce(-1:lab,-1:la,-1:lb,m)=0D0
      ce(0,0,0,m)=1D0
      if (lab.eq.0) return
      p  = a + b
      ab = ax - bx
      pa = - b*ab/p
      pb = + a*ab/p
      tp = 1d0/(2d0*p)
      do i=0,la
        i1=i-1
        do j=0,lb
          j1=j-1
          do t=1,i+j
            t1=t-1
            ce(t,i,j,m)=tp*(i*ce(t1,i1,j,m)+j*ce(t1,i,j1,m))/dble(t)
          enddo
          if (i.lt.la) ce(0,i+1,j,m)=pa*ce(0,i,j,m)+ce(1,i,j,m)
          if (j.lt.lb) ce(0,i,j+1,m)=pb*ce(0,i,j,m)+ce(1,i,j,m)
        enddo
      enddo 
      return
      end
