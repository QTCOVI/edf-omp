c
c-----------------------------------------------------------------------
c
c-----Driver routine for multipole.f
c
      subroutine newdrvmult (nmo,nprims,ncent,cc,stdout)
c
c.......................................................................
c
      include         'param.inc'
      real(kind=8),    allocatable,dimension (:,:)   :: rexp
      integer(kind=4), allocatable,dimension (:,:)   :: ipow
      real(kind=8)     pnt(3)
      real(kind=8)     cc(nmo,nprims)
      integer(kind=4)  stdout
      real(kind=8)     mu2px,mu2py,mu2pz
      real(kind=8)     mu4pxx,mu4pyy,mu4pzz,mu4pxy,mu4pxz,mu4pyz

      nexp=24
      allocate (rexp(nmo,nexp),stat=ier)
      if (ier.ne.0) stop 'newdrvmult.f: Cannot allocate rexp()'
      allocate (ipow(nexp,3),stat=ier)
      if (ier.ne.0) stop 'newdrvmult.f: Cannot allocate rexp()'
      rexp(:,:)  = 0.0d+00
      pnt (:)    = 0.0d+00
c
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
c
      ipow(13,1:3) = (/1,1,0/)
      ipow(14,1:3) = (/1,0,1/)
      ipow(15,1:3) = (/0,1,1/)
c
      ipow(16,1:3) = (/2,2,0/)
      ipow(17,1:3) = (/2,0,2/)
      ipow(18,1:3) = (/0,2,2/)
c
      ipow(19,1:3) = (/2,1,0/)
      ipow(20,1:3) = (/2,0,1/)
      ipow(21,1:3) = (/0,2,1/)
      ipow(22,1:3) = (/1,2,0/)
      ipow(23,1:3) = (/1,0,2/)
      ipow(24,1:3) = (/0,1,2/)
c
      maxty      = maxtype
      call multipole (pnt,nmo,nprims,ncent,maxty,nlm,nexp,ipow,cc,rexp)

      write (stdout,332) 
      do i=1,nmo
        xsig=sqrt(rexp(i,4)-rexp(i,1)**2)  ! SQRT [ <x^2> - <x>^2 ]
        ysig=sqrt(rexp(i,5)-rexp(i,2)**2)  ! SQRT [ <y^2> - <y>^2 ]
        zsig=sqrt(rexp(i,6)-rexp(i,3)**2)  ! SQRT [ <z^2> - <z>^2 ]
        write (stdout,333) i,(rexp(i,j),j=1,6),xsig,ysig,zsig
      enddo

      write (stdout,329)
      do i=1,nmo
        xsig=sqrt(rexp(i,4)-rexp(i,1)**2)
        ysig=sqrt(rexp(i,5)-rexp(i,2)**2)
        zsig=sqrt(rexp(i,6)-rexp(i,3)**2)
        mu4pxx=rexp(i,10) - 4d0 * rexp(i,7)*rexp(i,1)
     &  + 6d0 * rexp(i,4)*rexp(i,1)**2 - 3d0 * rexp(i,1)**4
        mu4pyy=rexp(i,11) - 4d0 * rexp(i,8)*rexp(i,2)
     &  + 6d0 * rexp(i,5)*rexp(i,2)**2 - 3d0 * rexp(i,2)**4
        mu4pzz=rexp(i,12) - 4d0 * rexp(i,9)*rexp(i,3)
     &  + 6d0 * rexp(i,6)*rexp(i,3)**2 - 3d0 * rexp(i,3)**4
        sig4xx=sqrt(sqrt(mu4pxx))
        sig4yy=sqrt(sqrt(mu4pyy))
        sig4zz=sqrt(sqrt(mu4pzz))
        bpxx=sig4xx/xsig
        bpyy=sig4yy/ysig
        bpzz=sig4zz/zsig
        write (stdout,333) i,mu4pxx,mu4pyy,mu4pzz,sig4xx,
     &       sig4yy,sig4zz,bpxx,bpyy,bpzz
      enddo
      deallocate (rexp,stat=ier)
      if (ier.ne.0) stop 'drvmult.f: Cannot allocate rexp()'
      deallocate (ipow,stat=ier)
      if (ier.ne.0) stop 'drvmult.f: Cannot allocate rexp()'
      return
 333  format (I8,10(1x,F13.6))
 330  format (3x,
     &'I.M. Hoyvik & P. Jorgensen [Chem. Rev. 116, 3306-3326 (2016)]',
     &/,3x,
     & 'I.M. Hoyvik, B. Jansik & P. Jorgensen [JCP 137, 2224114 (2012)'
     & /,3x,61('-'),/,4x,'Orbital',1x,'sigma_2^p',5x,'sigma_4^p',
     & 5x,'beta_p',8x,'beta_p^4',/,3x,61('-'))
 332  format (3x,
     &'I.M.Hoyvik & P.Jorgensen [Chem.Rev. 116, 3306-3326 (2016)]',3x,
     &'I.M.Hoyvik, B.Jansik & P.Jorgensen [JCP 137, 2224114 (2012)'
     & /,3x,'sigma(x) = SQRT [ <x^2> - <x>^2 ], ',
     & 'sigma(y) = SQRT [ <x^2> - <x>^2 ], '
     & 'sigma(z) = SQRT [ <z^2> - <z>^2 ]',
     & /,3x,131('-'),/,4x,'Orbital',6x,'<x>',11x,'<y>',11x,'<z>',10x,
     & '<x^2>',9x,'<y^2>',9x,'<z^2>',7x,'sigma(x)',6x,'sigma(y)',
     & 6x,'sigma(z)',/,3x,131('-'))
 329  format (3x,131('-'),/,4x,
     & 'mu_4^p(x)=<p| [x-<x>]^4 |p>  sigma_4(x)=mu_4^p(x)^(1/4)  ',
     & 'beta(x)=sigma_4(x)/sigma(x),...',/,3x,131('-'),
     & /,4x,'Orbital',3x,'mu_4^p(x)',5x,'mu_4^p(y)',5x,'mu_4^p(z)',5x,
     & 'sigma_4(x)',4x,'sigma_4(y)',4x,'sigma_4(z)',4x,
     & 'beta(x)',7x,'beta(y)',7x,'beta(z)',/,3x,131('-'))
      end
