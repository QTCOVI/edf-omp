c
c-----------------------------------------------------------------------
c
c-----Driver routine for multipole.f
c
      subroutine drvmult (nmo,nprims,ncent,cc,stdout)
c
c.......................................................................
c
      include         'param.inc'
      real(kind=8),    allocatable,dimension (:,:)   :: rexp
      integer(kind=4), allocatable,dimension (:,:)   :: ipow
      real(kind=8)     pnt(3)
      real(kind=8)     cc(nmo,nprims)
      integer(kind=4)  stdout

      nexp=12
      allocate (rexp(nmo,nexp),stat=ier)
      if (ier.ne.0) stop 'drvmult.f: Cannot allocate rexp()'
      allocate (ipow(nexp,3),stat=ier)
      if (ier.ne.0) stop 'drvmult.f: Cannot allocate rexp()'
      rexp(:,:)  = 0d0
      pnt (:)    = 0d0
      ipow( 1,1:3) = (/1,0,0/)
      ipow( 2,1:3) = (/0,1,0/)
      ipow( 3,1:3) = (/0,0,1/)
      ipow( 4,1:3) = (/2,0,0/)
      ipow( 5,1:3) = (/0,2,0/)
      ipow( 6,1:3) = (/0,0,2/)
      ipow( 7,1:3) = (/3,0,0/)
      ipow( 8,1:3) = (/0,3,0/)
      ipow( 9,1:3) = (/0,0,3/)
      ipow(10,1:3) = (/4,0,0/)
      ipow(11,1:3) = (/0,4,0/)
      ipow(12,1:3) = (/0,0,4/)
      maxty      = maxtype
      call multipole (pnt,nmo,nprims,ncent,maxty,nlm,nexp,ipow,cc,rexp)
      do i=1,nmo
        xmu2p=rexp(i,4)-rexp(i,1)**2
        ymu2p=rexp(i,5)-rexp(i,2)**2
        zmu2p=rexp(i,6)-rexp(i,3)**2
        if (xmu2p.lt.0d0) then
          stop ' # drvmult.f xmu2p iz negative !!!'
        endif
        if (ymu2p.lt.0d0) then
          stop ' # drvmult.f ymu2p iz negative !!!'
        endif
        if (zmu2p.lt.0d0) then
          stop ' # drvmult.f zmu2p iz negative !!!'
        endif
        xyzmu2=sqrt(xmu2p+ymu2p+zmu2p)
        xmu4p=rexp(i,10)
     &       - 4d0 * rexp(i,7)*rexp(i,1)
     &       + 6d0 * rexp(i,4)*rexp(i,1)**2
     &       - 3d0 * rexp(i,1)**4
        ymu4p=rexp(i,11)
     &       - 4d0 * rexp(i,8)*rexp(i,2)
     &       + 6d0 * rexp(i,5)*rexp(i,2)**2
     &       - 3d0 * rexp(i,2)**4
        zmu4p=rexp(i,12)
     &       - 4d0 * rexp(i,9)*rexp(i,3)
     &       + 6d0 * rexp(i,6)*rexp(i,3)**2
     &       - 3d0 * rexp(i,3)**4
        if (xmu4p.lt.0d0) then
          stop ' # drvmult.f xmu4p iz negative !!!'
        endif
        if (ymu4p.lt.0d0) then
          stop ' # drvmult.f ymu4p iz negative !!!'
        endif
        if (zmu4p.lt.0d0) then
          stop ' # drvmult.f zmu4p iz negative !!!'
        endif
        xyzmu4=sqrt(sqrt(xmu4p+ymu4p+zmu4p))
        betap=xyzmu4/xyzmu2
        betap4=xyzmu4**4/xyzmu2**4
        write (stdout,333) i,xyzmu2,xyzmu4,betap,betap4
      enddo

*     write (stdout,332) 
*     do i=1,nmo
*       xmu2p=rexp(i,4)-rexp(i,1)**2  ! <x^2> - <x>^2
*       ymu2p=rexp(i,5)-rexp(i,2)**2  ! <y^2> - <y>^2
*       zmu2p=rexp(i,6)-rexp(i,3)**2  ! <z^2> - <z>^2
*       write (stdout,333) i,(rexp(i,j),j=1,6),xmu2p,ymu2p,zmu2p,
*    &   xmu2p+ymu2p+zmu2p
*     enddo
*332  format (



      deallocate (rexp,stat=ier)
      if (ier.ne.0) stop 'drvmult.f: Cannot allocate rexp()'
      deallocate (ipow,stat=ier)
      if (ier.ne.0) stop 'drvmult.f: Cannot allocate rexp()'
      return
 333  format (I8,10(1x,E13.6))
      end
