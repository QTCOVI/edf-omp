c-----------------------------------------------------------------------
c
c-----Overlap integral between 2 un-normalized cartesian gaussians.
c     
c     xa^i ya^k za^m exp[-za ra^2] and xb^j yb^l zb^n exp[-zb rb^2]
c
c     where xa=x-ax, ya=y-ay, za=z-az, xb=x-bx, yb=y-by, zb=z-bz
c
c     and (ax,ay,az) and (ay,by,bz) are the coordinates of both centers.
c
c-----Input parameters are:
c
c     it(1..3) = i,k,m
c     jt(1..3) = j,l,n
c     ax(1..3) = ax,ay,az
c     bx(1..3) = bx,by,bz
c
c-----------------------------------------------------------------------
c
      real(kind=8) function  sabxyz (it,jt,ax,bx,za,zb)
      implicit real(kind=8) (a-h,o-z)
      parameter (pi=3.14159265358979323846D0)
      real(kind=8) za,zb,p,mu,res
      real(kind=8) xab(3),xpa(3),xpb(3),sabx(-1:9,-1:9,3),ax(3),bx(3)
      integer u,i,j,it(3),jt(3),ito,jto,m,q
c
      if (it(1).gt.9.or.it(2).gt.9.or.it(3).gt.9) then
         stop 'sabxyz: Center A: Too high L in the cartesian GTO'
      endif
      if (jt(1).gt.9.or.jt(2).gt.9.or.jt(3).gt.9) then
        stop 'sabxyz: Center B: Too high L in the cartesian GTO'
      endif
      p=za+zb
      sqpip = sqrt(pi/p)
      mu=za*zb/p
      res=1d0
      do u=1,3
        xab(u)=ax(u)-bx(u)
      enddo
      do u=1,3
        xpa(u)=-zb/p*xab(u)
        xpb(u)=+za/p*xab(u)
        sabx(-1,-1,u)=0d0
        sabx(-1, 0,u)=0d0
        sabx( 0,-1,u)=0d0
c
c.......S00 integral
c
        sabx(0,0,u) = sqpip 
        isum=it(u)+jt(u)
        do m=1,it(u)+jt(u)
          q=0
          do ito=m,0,-1
            q=q+1
            jto=m-ito
            if (q.le.m/2+1) then
              i=ito-1
              j=jto
              xp=xpa(u)
            else
              i=ito
              j=jto-1
              xp=xpb(u)
            endif
            sec=(i*sabx(i-1,j,u)+j*sabx(i,j-1,u))/2d0/p
            if (q.le.m/2+1) then
              sabx(i+1,j,u)=xp*sabx(i,j,u)+sec
            else
              sabx(i,j+1,u)=xp*sabx(i,j,u)+sec
            endif
          enddo
        enddo
        res=res*sabx(it(u),jt(u),u)
      enddo
      sabxyz = res*exp(-mu*(xab(1)**2+xab(2)**2+xab(3)**2))
      return
      end
