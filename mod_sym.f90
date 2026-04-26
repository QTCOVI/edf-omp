module mod_sym

  implicit none

  private

  ! ejemplo de uso
  !  call symheader(uout) --> opcional
  !  call syminit(uout, .true.) --> importante que aqui se definen precisiones y demas
  !  call symrwcoord(uout, mol%x(1,:), mol%x(2,:), mol%x(3,:), mol%z(:), &
  !                  mol%mass(:), .true., .false., mol%n) --> coordenadas,numero atomico y masa
  !  call symdeter(uout, mol%x(1,:), mol%x(2,:), mol%x(3,:), mol%n) --> tipo de molecula
  !  if (mol_linear) then
  !    call sym1d(uout, mol%x(1,:), mol%x(2,:), mol%x(3,:), mol%z(:), mol%n)
  !  else if (mol_planar) then
  !    call sym2d(uout, mol%x(1,:), mol%x(2,:), mol%x(3,:), mol%z(:), mol%n)
  !  else
  !    call sym3d(uout, mol%x(1,:), mol%x(2,:), mol%x(3,:), mol%z(:), mol%n)
  !  endif 
  !  call symprint(uout,.false.)
  !  mol_linear y mol_planar son variables del modulo !!!
  
  !
  ! symmetry params
  !
  ! TOLdist ........ Check of a distance.
  ! TOLeqvm ........ Comparison of matrices.
  ! TOLsng ......... Singularity of a 3x3 matrix.
  ! TOLisint ....... Closeness to an integer.
  ! TOLnull ........ Null angles and projections.
  ! TOLeigen ....... Eigenvalues of the orthogonal matrices.
  ! TOLdirty ....... Activate the algorithm changes to deal with
  !                  noisy data.
  ! inf_order ...... Infinity order will be transformed into this.
  ! TOLtriplet...... Used to select an initial triplet of atoms to
  !                  determine the symmetry matrices. This value
  !                  should be large (0.1 is the default). Increase
  !                  this value to ask for even better triplets.
  !
  logical, public :: TOLdirty
  integer, public :: inf_order
  real(kind=8), public :: TOLdist, TOLeqvm, TOLsng, TOLisint, TOLnull
  real(kind=8), public :: TOLeigen, TOLtriplet
  !
  ! The largest errors on the different tests.
  ! 
  ! ERRortho ....... Error on the orthogonality of the matrices.
  ! ERRsng ......... Best nonsingular triplet found.
  ! ERReigen ....... Error on the eigenvalues.
  !
  real(kind=8), public :: ERRsng, ERReigen
  ! 
  ! nopsym ........ Number of symmetry operations.
  ! optype() ...... Symmetry operation types (see below).
  ! oporder() ..... Rotation order (the n in C_n^m or S_n^m).
  ! opm() ......... The m in C_n^m or S_n^m.
  ! opinv() ....... Inverse of an operation.
  ! opsym(i,,) .... The 3x3 matrix of the operation in the xyz rep.
  ! opaxis(i,) .... Rotation axis.
  ! opeuler(i,).... Euler angles for the rotation axis.
  ! opangle()...... Rotation angle (radians).
  ! opsymbol() .... Symbol for the operation.
  ! opproper() .... Proper/improper rotation (redundant).
  ! opisgener() ... Is this a generator for the group?
  ! optable(,)..... Multiplication (Cayley) table.
  ! linear_mol .... True if the molecule is linear.
  ! point_group ... Point group symbol.
  ! point_g2 ...... Point group symbol (short version).
  ! opmainord ..... Order of the main axis.
  ! 
  integer, parameter, public :: MOPSYM = 480
  logical, public  :: opproper(MOPSYM), opisgener(MOPSYM)
  integer, public  :: nopsym, optype(MOPSYM), oporder(MOPSYM)
  integer, public  :: opm(MOPSYM), opinv(MOPSYM)
  integer, public  :: optable(MOPSYM,MOPSYM), opmainord
  real(kind=8), public :: opsym(MOPSYM,3,3)
  real(kind=8), public :: opaxis(MOPSYM,3), opeuler(MOPSYM,3)
  real(kind=8), public :: opangle(MOPSYM)
  character(len=10), public :: opsymbol(MOPSYM)
  character(len=20), public :: point_group
  character(len=8), public  :: point_g2
  !
  ! Matrices used to convert the molecule to a symmetry orientation:
  !
  real(kind=8), public :: or_mat(3,3), or_imat(3,3)
  ! 
  ! Classes of symmetry operations:
  !
  integer, parameter, public :: MCLAS = 120, MINCLAS = 200
  integer, public :: nclas, opclas(MOPSYM)
  integer, public :: ninclas(MCLAS), iinclas(MCLAS,MINCLAS)
  !
  ! Sym operator types:
  !
  integer, parameter, public :: opt_identity     = 10
  integer, parameter, public :: opt_rotation     = 11
  integer, parameter, public :: opt_inversion    = 12
  integer, parameter, public :: opt_sigma        = 13
  integer, parameter, public :: opt_imp_rotation = 14
  integer, parameter, public :: opt_unknown      = 15
  !
  ! xcm .. zcm ...... center of mass coordinates.
  ! norbit .......... number of orbits.
  ! natorb() ........ number of atoms in an orbit.
  ! iatorb(,) ....... list of atoms in an orbit.
  ! orbZ() .......... atomic number of the atoms in this orbit.
  ! orbdis() ........ distance to the center for all
  !                   atoms in this orbit.
  ! orbmol() ........ orbit to which an atom belongs.
  ! molradius ....... maximum distance from the center to an atom.
  !
  integer, parameter, public :: MORBIT = 2000, MATORB = 300, MMOL1 = 300000
  integer, public  :: norbit, natorb(MORBIT), iatorb(MORBIT,MATORB)
  integer, public  :: orbZ(MORBIT), orbmol(MMOL1)
  real(kind=8), public :: orbdis(MORBIT), molradius, xcm, ycm, zcm

  ! public routines
  public :: sym1d, sym2d, sym3d, symgetgroup, syminit, symheader
  public :: symback, symcleancoord, symclosure, symdeter, symeigen 
  public :: symeuler, symimagic, syminv, symmatfill, symmatprod 
  public :: symnoise, symopadd, symorb, symprint, symproject, sympurify
  public :: symreplicate, symrwcoord, symtransform, symorient, symopchk

  ! aux logical data
  logical, public :: mol_linear, mol_planar, linear_mol

  ! internal variables and functions
  real(kind=8), parameter, public :: dpi = 3.14159265358979323846d0
  integer, parameter :: faterr = -1 !< fatal error flag
  integer, parameter :: warning = 1 !< warning flag
  integer, parameter :: noerr = 0 !< info flag
  integer :: nwarns = 0 !< Number of warnings
  integer :: ncomms = 0 !< Number of comments
  character(len=1), parameter :: null = char(0) !< null character (ascii 0)
  character(len=1), parameter :: blank = " " !< blank

contains

  !> Header of module
  subroutine symheader(uout)

    integer, intent(in) :: uout

!    write (uout,800)
!800 format (/                                                       &
!    1x, ' SSSSS YY  YY MM   MM MM   MM EEEEE TTTTTT RRRRRR YY  YY'/ &
!    1x, 'SS     YY  YY MMMMMMM MMMMMMM EEEEE TTTTTT RRRRRR YY  YY'/ &
!    1x, ' SS    YY  YY MM M MM MM M MM EE      TT   RR  RR YY  YY'/ &
!    1x, '  SS   YY  YY MM   MM MM   MM EE      TT   RR  RR YY  YY'/ &
!    1x, '   SS   YYYY  MM   MM MM   MM EEEEE   TT   RRRRRR  YYYY '/ &
!    1x, '    SS   YY   MM   MM MM   MM EE      TT   RR RR    YY  '/ &
!    1x, '    SS   YY   MM   MM MM   MM EEEEE   TT   RR  RR   YY  '/ &
!    1x, 'SSSSS    YY   MM   MM MM   MM EEEEE   TT   RR  RR   YY  '//)

      WRITE (uout,111)
 111  format (/                                            &
     1x,'# ---------------------------------'/ &
     1x,'# S Y M M E T R Y   A N A L Y S I S'/ &
     1x,'# ---------------------------------')

  end subroutine symheader

  !> Init data needed for symmetry
  subroutine syminit(uout, verbose)
  
    implicit none
    logical, intent(in) :: verbose
    integer, intent(in) :: uout
      
    ! default values  
    TOLdist  = 3E-5
    TOLeqvm  = 2.0d0 * TOLdist
    TOLsng   = 1E-5
    TOLisint = 3E-5
    TOLnull  = 3E-6
    TOLeigen = 3E-5
    TOLdirty = .false.
    TOLtriplet = 0.1d0
    ERRsng   = 0d0
    ERReigen = 0d0
    inf_order  = 16
    mol_linear = .false.
    mol_planar = .false.

    if (verbose) then
      write (uout,820) TOLdist, TOLeqvm, TOLsng, TOLisint, TOLnull, &
                       TOLeigen, TOLtriplet
    endif

820  format (/                                                  &
     1x, 'ALGORITHM(symmetry) Set of tolerances:'/              &
     1x, 'TOLdist (distance comparison) ....... ', 1p, e12.5 /  &
     1x, 'TOLeqvm (matrix equivalence) ........ ', 1p, e12.5 /  &
     1x, 'TOLsng (matrix singularity) ......... ', 1p, e12.5 /  &
     1x, 'TOLisint (closeness to an integer) .. ', 1p, e12.5 /  &
     1x, 'TOLnull (null angles & projections) . ', 1p, e12.5 /  &
     1x, 'TOLeigen (eigenvalues of sym matr) .. ', 1p, e12.5 /  &
     1x, 'TOLtriplet (quality of atom triplet). ', 1p, e12.5 )

  end subroutine syminit

  !> Sym1d - check symmetry of linear molecules.
  subroutine sym1d(uout, ax, ay, az, atZmol, nmol)

  implicit none

  integer :: uout, nmol
  integer, dimension(nmol) :: atZmol
  real(kind=8), dimension(nmol) :: ax, ay, az
 
  real(kind=8), dimension(3,3) :: xmat
  real(kind=8) :: xx, yy, zz, xnorm, alfa, theta, phi
  real(kind=8) ::  sa, ca, st, ct, sp, cp, s2t, c2t, s2p, c2p
  integer :: i, j, norder!, leng
  !logical :: symopchk
 
  ! Classify all atoms into orbits. All atoms in an orbit have the
  ! same atomic number and the same distance to the center of mass:
  call symorb (uout, ax, ay, az, atZmol, nmol)
 
  linear_mol = .true.
  nopsym = 0

  ! The identity is always a sym operator:
  call symmatfill(xmat, 1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0)
  call symopadd(xmat)
  !
  ! Is the inversion a sym op?
  !
  call symmatfill (xmat, -1d0,0d0,0d0, 0d0,-1d0,0d0, 0d0,0d0,-1d0)
  if (symopchk(xmat, ax, ay, az, atZmol, nmol)) then
    call symopadd (xmat)
    point_group = "D_inf_h"
  else
    point_group = "C_inf_v"
  endif
  point_g2 = point_group
  opmainord = 999
  !
  ! I don't know how to manage an infinity number of rotation
  ! operations. I will enforce a Dnh or Cnv group instead with n=4,
  ! for instance.
  !
  norder = inf_order
  !
  ! Determine the molecular axis and transform the direction into
  ! spherical polar angles:
  !
  xnorm = 0d0
  do i = 1, nmol
    if (ax(i)*ax(i) + ay(i)*ay(i) + az(i)*az(i) .gt. xnorm) then
      xx = ax(i)
      yy = ay(i)
      zz = az(i)
      xnorm = xx*xx + yy*yy + zz*zz
    endif
  enddo
  if(xnorm .le. 0d0) call error('sym1d', 'Molecular axis not found!', faterr)
  xnorm = sqrt(xnorm)
  theta = acos(zz/xnorm)
  phi = atan2(yy, xx)
  write(uout,600) xx/xnorm, yy/xnorm, zz/xnorm, theta*180d0/dpi, phi*180d0/dpi
  !write (uout,605) point_group(1:leng(point_group)), norder
  write (uout,605) point_group, norder
  !
  ! Compose the C_n^1 rotation matrix around the molecular axis:
  !
  alfa = (dpi+dpi) / norder
  ca = cos(alfa)
  sa = sin(alfa)
  ct = cos(theta)
  st = sin(theta)
  xmat(1,1) = ca*ct*ct + st*st
  xmat(1,2) = -sa*ct
  xmat(1,3) = (1-ca)*st*ct
  xmat(2,1) = sa*ct
  xmat(2,2) = ca
  xmat(2,3) = -sa*st
  xmat(3,1) = (1-ca)*st*ct
  xmat(3,2) = sa*st
  xmat(3,3) = ct*ct + ca*st*st
  if (symopchk(xmat, ax, ay, az, atZmol, nmol)) then
    call symopadd (xmat)
  else
    write (uout,610) ((xmat(i,j), j = 1, 3), i = 1, 3)
    call error ('sym1d', 'Algorithm error: C_n^1 not sym op!', faterr)
  endif
  !
  ! Compose a sigma_v operation and test it:
  !
  cp = cos(phi)
  sp = sin(phi)
  c2t = cos(theta+theta)
  s2t = sin(theta+theta)
  c2p = cos(phi+phi)
  s2p = sin(phi+phi)
  xmat(1,1) = ct*ct*c2p + st*st
  xmat(1,2) = ct*s2p
  xmat(1,3) = s2t*sp*sp
  xmat(2,1) = ct*s2p
  xmat(2,2) = -c2p
  xmat(2,3) = -2d0 * st*sp*cp
  xmat(3,1) = s2t*sp*sp
  xmat(3,2) = -2d0 * st*sp*cp
  xmat(3,3) = ct*ct + st*st*c2p
  if(symopchk(xmat, ax, ay, az, atZmol, nmol)) then
    call symopadd (xmat)
  endif
!
! Check the closure of the symmetry operators set:
!
  call symclosure (uout)
  call symgetgroup (uout)
 
600 format (/                                                    &
    1x, '++SYM1D: Symmetry of a linear molecule'/                &
    1x, 'Molecular axis direction:       ', 3f15.6/              &
    1x, 'Direction in polar coordinates: ', 2f15.6, ' degrees')

605 format (                                                     &
    1x, 'True symmetry group: ', a/                              &
    1x, 'The infinite order axis will be changed to order: ', i4)

610 format (/                                                    &
    1x, 'ALG(sym1d) Rotation around the molecular axis not',     &
    1x, 'recognized as a symmetry operation:'/                   &
    (1x, 3f15.6))
 
end subroutine sym1d

subroutine sym2d (uout, ax, ay, az, atZmol, nmol)
!
! Sym2d - check symmetry of planar molecules.
!
  implicit none
 
  integer :: uout, nmol 
  integer, dimension(nmol) :: atZmol
  real(kind=8), dimension(nmol) :: ax(nmol), ay(nmol), az(nmol)
 
  real(kind=8), dimension(3,3) :: v, vinv, xmat, xm, xv, xop
  real(kind=8) :: vdet, xx, yy, zz, xnorm, ztest, xdet
  integer :: i, j, ierror, ii, jj, i1, i2, j1, j2
  !logical :: symopchk

!
! Classify all atoms into orbits. All atoms in an orbit have the
! same atomic number and the same distance to the center of mass:
!
  call symorb (uout, ax, ay, az, atZmol, nmol)
!
!
!
!.....Get the new coordinates with Z being perpendicular to the
!     molecular plane:
!
  call symmatfill (v, 1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0)
  call symmatfill (vinv, 1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0)
  do i = 1, nmol
    do j = i+1, nmol
      xx = ay(i)*az(j) - az(i)*ay(j)
      yy = az(i)*ax(j) - ax(i)*az(j)
      zz = ax(i)*ay(j) - ay(i)*ax(j)
      xnorm = xx*xx + yy*yy + zz*zz
      if(xnorm .gt. 1d-1) goto 1100
    enddo
  enddo
  call error ('sym2d', 'Normal to the molecular plane not found!', faterr)
1100 continue

  xnorm = 1/sqrt(xnorm)
  v(1,3) = xx * xnorm
  v(2,3) = yy * xnorm
  v(3,3) = zz * xnorm

!
! Check if the molecule is already in the XY plane and no
! reorientation is needed:
!
  if (abs(v(1,3)).le.TOLdist .and. abs(v(2,3)).le.TOLdist .and. &
      abs(v(3,3)-1d0).le.TOLdist) goto 1103
!
! Otherwise get the two  directions of the molecular plane:
! (they will be converted into the new XY axes)
!
  do i = 1, nmol
    xx = ax(i)
    yy = ay(i)
    zz = az(i)
    xnorm = xx*xx + yy*yy + zz*zz
    if (xnorm .gt. 1d-2) goto 1101
  enddo
  call error ('sym2d', 'New x axis not found!', faterr)
1101 continue

  xnorm = 1/sqrt(xnorm)
  v(1,1) = xx * xnorm
  v(2,1) = yy * xnorm
  v(3,1) = zz * xnorm

!
! The vectorial product v3 x v1 = v2 should be normalized:
!
  xx = v(2,3)*v(3,1) - v(3,3)*v(2,1)
  yy = v(3,3)*v(1,1) - v(1,3)*v(3,1)
  zz = v(1,3)*v(2,1) - v(2,3)*v(1,1)
  xnorm = xx*xx + yy*yy + zz*zz
  xnorm = 1/sqrt(xnorm)
  v(1,2) = xx * xnorm
  v(2,2) = yy * xnorm
  v(3,2) = zz * xnorm
!
! The V and V^{-1} matrices produce the change of coordinates and
! transform the symmetry matrices from the old to the new
! coordinates:
!
  call syminv (v, vinv, vdet, ierror)

  if (ierror.ne.0) call error ('sym2d', &
  'Error in the inversion of the V matrix!', faterr)
  ztest = 0d0
  do i = 1, nmol
    xx = vinv(1,1)*ax(i) + vinv(1,2)*ay(i) + vinv(1,3)*az(i)
    yy = vinv(2,1)*ax(i) + vinv(2,2)*ay(i) + vinv(2,3)*az(i)
    zz = vinv(3,1)*ax(i) + vinv(3,2)*ay(i) + vinv(3,3)*az(i)
    ztest = ztest + abs(zz)
    ax(i) = xx
    ay(i) = yy
    az(i) = zz
  enddo
  if(ztest .gt. 1d-4) then
    call error ('sym2d','Transformed atoms do not fulfill z=0 test!', warning)
    write (0,500) ((v(i,j), j = 1, 3), i = 1, 3), ((vinv(i,j), j = 1, 3), i = 1, 3)
    write (0,505) (i, ax(i), ay(i), az(i), i = 1, nmol)
  endif

500 format (/                                                   &
    1x, 'DBG(sym2d) Planar molecule transformation problem: '/  &
    1x, 'Transformation matrix and inverse: '/                  &
    (1x, 3f16.9))

505 format (                                                   &
    1x, 'DBG(sym2d) Transformed molecular coordinates:'/       &
    1x, i5, 3f16.9)

!
!.....Start looking for symmetry operations:
!
1103 continue

!
!.....The identity is always a sym operator:
!
  linear_mol = .false.
  nopsym = 0
  call symmatfill (xmat, 1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0)
  call symopadd (xmat)
!
!.....Is the inversion a sym op?
!
  call symmatfill (xmat, -1d0,0d0,0d0, 0d0,-1d0,0d0, 0d0,0d0,-1d0)
  if (symopchk(xmat, ax, ay, az, atZmol, nmol)) then
    call symopadd (xmat)
  endif
!
!.....Reflection in the molecular plane should be a symmetry operation:
!
  call symmatfill (xmat, 1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,-1d0)
  if (symopchk(xmat, ax, ay, az, atZmol, nmol)) then
    call symopadd (xmat)
  endif
!
!.....Run over all doblet pairs. Each pair migth produce a new symmetry
!     operator. Only doblets having the same atoms and in the same
!     order are compatible.
!     The position matrix of the first doblet must be invertible for
!     the algorithm to work.
!
    do i1 = 1, nmol
      xmat(1,1) = ax(i1)
      xmat(2,1) = ay(i1)
      do i2 = i1+1, nmol
        xmat(1,2) = ax(i2)
        xmat(2,2) = ay(i2)
        xdet = xmat(1,1)*xmat(2,2) - xmat(1,2)*xmat(2,1)
        if(abs(xdet) .le. TOLsng) goto 1001
        xv(1,1) = xmat(2,2)/xdet
        xv(1,2) = -xmat(1,2)/xdet
        xv(2,1) = -xmat(2,1)/xdet
        xv(2,2) = xmat(1,1)/xdet
        ERRsng = max(ERRsng, abs(xdet))
        do j1 = 1, nmol
          if(orbmol(i1).ne.orbmol(j1)) goto 1002
          xm(1,1) = ax(j1)
          xm(2,1) = ay(j1)
          do j2 = 1, nmol
            if(orbmol(i2).ne.orbmol(j2)) goto 1003
            xm(1,2) = ax(j2)
            xm(2,2) = ay(j2)
            xop(1,1) = xm(1,1)*xv(1,1) + xm(1,2)*xv(2,1)
            xop(1,2) = xm(1,1)*xv(1,2) + xm(1,2)*xv(2,2)
            xop(2,1) = xm(2,1)*xv(1,1) + xm(2,2)*xv(2,1)
            xop(2,2) = xm(2,1)*xv(1,2) + xm(2,2)*xv(2,2)
            xop(1,3) = 0d0
            xop(2,3) = 0d0
            xop(3,1) = 0d0
            xop(3,2) = 0d0
            xop(3,3) = 1d0
!.................Check if this is a new sym operator:
            if(symopchk(xop,ax,ay,az,atZmol,nmol)) then
              call symopadd (xop)
              if(TOLdirty) call symclosure (uout)
            endif
1003        continue
          enddo
1002      continue
        enddo
1001    continue
      enddo
    enddo
!
!.....Transform back the coordinates and the symmetry operations to the
!     original orientation:
!
    do i = 1, nmol
    xx = v(1,1)*ax(i) + v(1,2)*ay(i) + v(1,3)*az(i)
    yy = v(2,1)*ax(i) + v(2,2)*ay(i) + v(2,3)*az(i)
    zz = v(3,1)*ax(i) + v(3,2)*ay(i) + v(3,3)*az(i)
    ax(i) = xx
    ay(i) = yy
    az(i) = zz
  enddo
  call symtransform (v, vinv)
!
!.....Check the closure of the symmetry operators set:
!
  call symclosure (uout)
  call symgetgroup (uout)

end subroutine sym2d

  subroutine sym3d (uout, ax, ay, az, atZmol, nmol)
!
!.....sym3d - check symmetry of nonlinear molecules.
!
    implicit none

    integer :: nmol, uout
    integer, dimension(nmol) :: atZmol
    real(kind=8), dimension(nmol) :: ax, ay, az

    integer :: i, j, k, ii, jj, i1, i2, i3, j1, j2, j3, ierror, ii1, ii2, ii3
    real(kind=8) :: dbest, xdet
    real(kind=8), dimension(3,3) :: xmat, xm, xv, xop
    logical :: symsingular
    integer :: ntest1, ntest2

!
!.....Classify all atoms into orbits. All atoms in an orbit have the
!     same atomic number and the same distance to the center of mass:
!
    call symorb (uout, ax, ay, az, atZmol, nmol)

!
!.....The identity is always a sym operator:
!
    linear_mol = .false.
    nopsym = 0
    call symmatfill (xmat, 1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0)
    call symopadd (xmat)

!
!.....Is the inversion a sym op?
!
    call symmatfill (xmat, -1d0,0d0,0d0, 0d0,-1d0,0d0, 0d0,0d0,-1d0)
    if (symopchk(xmat, ax, ay, az, atZmol, nmol)) call symopadd (xmat)

!
!.....Find a linear independent triplet of atoms.
!     Try first to find a very good triplet, or use the best available
!     one anyway.
!
    dbest = 0d0
    do i1 = 1, nmol
      xmat(1,1) = ax(i1)
      xmat(2,1) = ay(i1)
      xmat(3,1) = az(i1)
      do i2 = i1+1, nmol
        xmat(1,2) = ax(i2)
        xmat(2,2) = ay(i2)
        xmat(3,2) = az(i2)
        do i3 = i2+1, nmol
          xmat(1,3) = ax(i3)
          xmat(2,3) = ay(i3)
          xmat(3,3) = az(i3)
          call syminv (xmat, xv, xdet, ierror)
          if(ierror .eq. 0) then
            if(abs(xdet) .gt. TOLtriplet) goto 1001
            if(abs(xdet) .gt. dbest) then
              dbest = abs(xdet)
              ii1 = i1
              ii2 = i2
              ii3 = i3
            endif
          end if
        enddo ! i3
      enddo ! i2
    enddo ! i1
    xmat(1,1) = ax(ii1)
    xmat(2,1) = ay(ii1)
    xmat(3,1) = az(ii1)
    xmat(1,2) = ax(ii2)
    xmat(2,2) = ay(ii2)
    xmat(3,2) = az(ii2)
    xmat(1,3) = ax(ii3)
    xmat(2,3) = ay(ii3)
    xmat(3,3) = az(ii3)
    call syminv (xmat, xv, xdet, ierror)
    if(ierror .ne. 0) then
      call error('sym3d','singular triplet matrix',faterr)
    end if

!
!.....Run over all triplets that are compatible with the linear
!     independent triplet already found. Each compatible triplet migth
!     produce a new symmetry operator. To be compatible, both triplets
!     must be formed by the same atoms in the same order.
!
1001 continue

    ERRsng = abs(xdet)
    ntest1 = 0
    ntest2 = 0
    do j1 = 1, nmol
      if(orbmol(i1).ne.orbmol(j1)) goto 1002
      xm(1,1) = ax(j1)
      xm(2,1) = ay(j1)
      xm(3,1) = az(j1)
      do j2 = 1, nmol
        if(orbmol(i2).ne.orbmol(j2)) goto 1003
        if (j1.eq.j2) goto 1003
        xm(1,2) = ax(j2)
        xm(2,2) = ay(j2)
        xm(3,2) = az(j2)
        do j3 = 1, nmol
          if(orbmol(i3).ne.orbmol(j3)) goto 1004
          if(j1.eq.j3 .or. j2.eq.j3) goto 1004
          ntest1 = ntest1 + 1
          xm(1,3) = ax(j3)
          xm(2,3) = ay(j3)
          xm(3,3) = az(j3)
          do ii = 1, 3
            do jj = 1, 3
              xop(ii,jj) = 0d0
              do k = 1, 3
                xop(ii,jj) = xop(ii,jj) + xm(ii,k) * xv(k,jj)
              enddo
            enddo
          enddo
!..............Check if this is a new sym operator:
          if(symopchk(xop,ax,ay,az,atZmol,nmol)) then
            ntest2 = ntest2 + 1
            call symopadd (xop)
            if(TOLdirty) call symclosure (uout)
          endif
1004      continue
        enddo  !j3
1003    continue
      enddo  !j2
1002  continue
    enddo  !j1
 
    write (uout,600) ntest1, ntest2

600 format (                                           &
    1x, 'ALG(sym3d) Triplet pairs tested:     ', i12/  &
    1x, 'ALG(sym3d) Possible operators found: ', i12)

!
!.....Check the closure of the symmetry operators set:
!
    call symclosure (uout)
    call symgetgroup (uout)

  end subroutine sym3d

subroutine symback(uout, ax, ay, az, n)
!
! symback - put back the molecule to the original orientation.
! The routine uses the orientation matrices obtained by symorient()
! to reverse its action.
!
  implicit none

  integer  :: n
  integer  :: uout
  real(kind=8) :: ax(n), ay(n), az(n)

  integer  :: i, j, iat
  real(kind=8) :: xx, yy, zz, v(3,3), vinv(3,3)

  if (uout.gt.0) write (uout,600)
! 
! Get the transformation matrices from the symmetry database:
!
  do i = 1, 3
    do j = 1, 3
      v(i,j) = or_mat(i,j)
      vinv(i,j) = or_imat(i,j)
    enddo
  enddo
!
! Transform the local copy of the coordinates:
!
  do iat = 1, n
    xx = v(1,1)*ax(iat) + v(1,2)*ay(iat) + v(1,3)*az(iat)
    yy = v(2,1)*ax(iat) + v(2,2)*ay(iat) + v(2,3)*az(iat)
    zz = v(3,1)*ax(iat) + v(3,2)*ay(iat) + v(3,3)*az(iat)
    ax(iat) = xx
    ay(iat) = yy
    az(iat) = zz
  enddo
!
! Transform the symmetry matrices too:
!
      call symtransform (v, vinv)
!
! Print out the transformed coordinates?
!
  if (uout .lt. 0) return

  write (uout,615) (i, ax(i), ay(i), az(i), i = 1, n)
  write (uout,620)
 
600 format (/                                               &
    1x, '++SYMBACK:'/                                       &
    1x, 'Return the molecule to the original orientation.')

615 format (/                                               &
    1x, 'Transformed molecular coordinates:'/               &
    (1x, i5, 3f16.9))

620 format (/                                               &
    1x, 'END OF SYMBACK'/)
 
end subroutine symback
subroutine symcleancoord (uout, ax, ay, az, n)
!
!.....symcleancoord - Write down the coordinates of the perfectly
!     symmetric molecule with the same orientation and atomic numbering.
!
    implicit none
 
    integer :: uout, n
    real(kind=8), dimension(n) :: ax, ay, az

    integer :: iter, ip, iq, ir
    integer :: i, j, k, kmin, ncount, ntrouble, norb
    integer :: iorb, iop, iat, kop, kbest
    integer, dimension(MOPSYM) :: igat
    real(kind=8) :: alpha, ca, sa, tt, nx, ny, nz
    real(kind=8) :: ximg, yimg, zimg, dd, dd1
    real(kind=8) :: xxx, yyy, zzz, dbest, dj, TOL2
    real(kind=8), dimension(3,3) :: xmat
    real(kind=8), dimension(n) :: axx, ayy, azz
    logical :: doagain, found
    logical, dimension(MOPSYM) :: generated
    logical, dimension(n) :: indep

!
!.....Before the coordinates can be cleaned, the symmetry matrices
!     must be cleaned themselves. We start by cleaning the rotation
!     angles of the generating matrices:
!
    do i = 1, nopsym
      if(optype(i).eq.opt_identity .or.  &
         optype(i).eq.opt_inversion) then
        generated(i) = .true.
      else if(opisgener(i)) then
        alpha = opm(i) * 2d0 * dpi / oporder(i)
        alpha = sign(alpha,opangle(i))
        ca = cos(alpha)
        sa = sin(alpha)
        nx = opaxis(i,1)
        ny = opaxis(i,2)
        nz = opaxis(i,3)
        if (opproper(i)) then
           tt = 1d0 - ca
        else
           tt = -1d0 - ca
        endif
        xmat(1,1) = nx*nx*tt + ca
        xmat(1,2) = nx*ny*tt - nz*sa
        xmat(1,3) = nx*nz*tt + ny*sa
        xmat(2,1) = nx*ny*tt + nz*sa
        xmat(2,2) = ny*ny*tt + ca
        xmat(2,3) = ny*nz*tt - nx*sa
        xmat(3,1) = nx*nz*tt - ny*sa
        xmat(3,2) = ny*nz*tt + nx*sa
        xmat(3,3) = nz*nz*tt + ca
        do j = 1, 3
          do k = 1, 3
            opsym(i,j,k) = xmat(j,k)
          enddo
        enddo
        generated(i) = .true.
      else
        generated(i) = .false.
      endif
    enddo

!
!.....Regenerate the symmetry operations from the cleaned generators:
!
    doagain = .true.
    iter = 0d0
    do while(doagain .and. iter.lt.100)
      doagain = .false.
      iter = iter + 1
      do i = 1, nopsym
        if(generated(i)) then
          do j = 1, nopsym
            k = optable(i,j)
            if(generated(j) .and. (.not.generated(k))) then
!.............This is a non regenerated opsym:
              doagain = .true.
              do ip = 1, 3
                do iq = 1, 3
                  xmat(ip,iq) = 0d0
                  do ir = 1, 3
                    xmat(ip,iq) = xmat(ip,iq) + opsym(i,ip,ir) * opsym(j,ir,iq)
                  enddo
                enddo
              enddo
              do ip = 1, 3
                do iq = 1, 3
                  opsym(k,ip,iq) = xmat(ip,iq)
                enddo
              enddo
              generated(k) = .true.
            endif
          enddo
        endif
      enddo
    enddo
    if (doagain) call error ('symcleancoord',                                &
                 'The group was not regenerated in 100 iterations!', faterr)

!
!.....Regenerate the symmetry matrices properties:
!
    call symmatfill (xmat,1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0)
    call symtransform (xmat, xmat)

!
!.....Clean the atomic coordinates:
!
    goto 101

!
!.....Method 1: Use an atom from each orbit to regenerate its orbit.
!.....Make sure that each atom is counted only once. When an atom is
!     counted as a symmetric image of another one, it is labelled as
!     already counted.
!
101 continue
    do i = 1, n
      indep(i) = .true.
    enddo

!
!.....Run over the uncounted atoms:
!
    write (uout,700)
    ncount = 0
    ntrouble = 0
    norb = 0
    do i = 1, n
      if (indep(i)) then
!
!...........Apply all the symmetry operations to get all images of the
!           atom:
!
        norb = norb + 1
        do j = 1, nopsym
          ximg = opsym(j,1,1)*ax(i) + opsym(j,1,2)*ay(i) + opsym(j,1,3)*az(i)
          yimg = opsym(j,2,1)*ax(i) + opsym(j,2,2)*ay(i) + opsym(j,2,3)*az(i)
          zimg = opsym(j,3,1)*ax(i) + opsym(j,3,2)*ay(i) + opsym(j,3,3)*az(i)
          found = .false.
          dd = 1e30
          kmin = -1
          k = 1
          do while(.not.found .and. k.le.n)
            dd1 = (ximg-ax(k))**2 + (yimg-ay(k))**2 + (zimg-az(k))**2
              if(dd1 .le. dd) then
                dd = dd1
                kmin = k
              endif
              if(sqrt(dd1) .le. TOLdist) then
                found = .true.
              else
                k = k + 1
              endif
          enddo

!
!..............Store the image in lieu of the closest old atom:
!
          if(found .and. indep(k)) then
            indep(k) = .false.
            axx(k) = ximg
            ayy(k) = yimg
            azz(k) = zimg
            ncount = ncount + 1
            write(uout,705) k, norb, j, ximg, yimg, zimg
          else if(indep(kmin)) then
            k = kmin
            indep(kmin) = .false.
            axx(k) = ximg
            ayy(k) = yimg
            azz(k) = zimg
            ncount = ncount + 1
            ntrouble = ntrouble + 1
            write(uout,705) k, norb, j, ximg, yimg, zimg
          endif
        enddo
      endif
    enddo
 
!     goto 103

! 
!.....Method 2: Substitute each atom by the average of all its images.
!
102 continue
    TOL2 = TOLdist * TOLdist
    do iorb = 1, norbit
      do i = 1, natorb(iorb)
        iat = iatorb(iorb,i)
        axx(iat) = 0d0
        ayy(iat) = 0d0
        azz(iat) = 0d0
        do iop = 1, nopsym
!.........Action of the opsym over this atom
          xxx = opsym(iop,1,1)*ax(iat) + opsym(iop,1,2)*ay(iat) &
                + opsym(iop,1,3)*az(iat)
          yyy = opsym(iop,2,1)*ax(iat) + opsym(iop,2,2)*ay(iat) &
                + opsym(iop,2,3)*az(iat)
          zzz = opsym(iop,3,1)*ax(iat) + opsym(iop,3,2)*ay(iat) &
                + opsym(iop,3,3)*az(iat)
!.........Identify the image atom:
          kbest = 0
          dbest = 1d30
          j = 1
          found = .false.
          do while(.not.found .and. j.le.natorb(iorb))
            k = iatorb(iorb,j)
            dj = (xxx-ax(k))**2 + (yyy-ay(k))**2 + (zzz-az(k))**2
            if(dj.le.dbest) then
              kbest = k
              dbest = dj
            endif
            if(dj .le. TOL2) then
              found = .true.
            else
              j = j + 1
            endif
          enddo
          if(.not.found) then
            write(uout,650) iat, iop, kbest, dbest
            k = kbest
          endif
!.........The inverse operation will convert the image atom
!         into the current atom:
          kop = opinv(iop)
          xxx = opsym(kop,1,1)*ax(k) + opsym(kop,1,2)*ay(k) + opsym(kop,1,3)*az(k)
          yyy = opsym(kop,2,1)*ax(k) + opsym(kop,2,2)*ay(k) + opsym(kop,2,3)*az(k)
          zzz = opsym(kop,3,1)*ax(k) + opsym(kop,3,2)*ay(k) + opsym(kop,3,3)*az(k)
          axx(iat) = axx(iat) + xxx
          ayy(iat) = ayy(iat) + yyy
          azz(iat) = azz(iat) + zzz
        enddo
        axx(iat) = axx(iat) / nopsym
        ayy(iat) = ayy(iat) / nopsym
        azz(iat) = azz(iat) / nopsym
      enddo
    enddo

!
!.....write down the coordinates and actualize the improved positions:
!
103 continue

    write (uout,605) n, ncount, ntrouble
    do i = 1, n
      write (uout,610) i, axx(i), ayy(i), azz(i)                  &
                       , axx(i)-ax(i), ayy(i)-ay(i), azz(i)-az(i)
      ax(i) = axx(i)
      ay(i) = ayy(i)
      az(i) = azz(i)
    enddo
 
605 format (/                                                        &
    1x, '++SYMCLEANCOORD:'/                                          &
    1x, 'Write down the coordinates of the perfectly symmetric'/     &
    1x, 'molecule with the same orientation and atomic numbering.'/  &
    1x, 'Original/transformed/problematic atoms: ', 3i6/             &
    1x, 'New coordinates and (new-old) displacements:'/              &
    1x, '-i------new-x----------new-y----------new-z------',         &
        '---dif-x-----dif-y-----dif-z--')

610 format (i4, 3f15.9, 1p, 3e10.2)

650 format (                                                         &
    1x, 'DBG(symcleancoord): sym/at image not found', 3i5, 1p, e15.6)

700 format (/                                                        &
    1x, '++SYMCLEANCOORD/ORBIT:'/                                    &
    1x, 'Classification of the molecule into orbits of atoms',       &
    1x, 'equivalent by symmetry.'/                                   &
    1x, '--iat--orbit-symop',                                        &
        '----new-x----------new-y----------new-z------')

705 format (1x, 3i6, 6f15.9)
 
  end subroutine symcleancoord
subroutine symclosure (uout)
!
! Symclosure - check the closure of the symmetry operators set
! and get the properties of the symmetry operators.
!

  implicit none

  integer :: uout

  real(kind=8) :: xmat(3,3), xinv(3,3), xdet, xtrace, angle, ang1
  real(kind=8) :: xval(3), xvect(3,3), xnm, xnorm
  integer  :: i, j, ii, jj, kk, inum, old_nopsym, ini_nopsym
  logical  :: newop, found, generated(MOPSYM)
  integer  :: ierror, iop, jinv, iclas
  character*(40) :: fmt
                   
  if (uout.ge.0) write (uout,600)
!
! Initialize the Cayley table:
!
  do i = 1, MOPSYM
    do j = 1, MOPSYM
      optable(i,j) = 0
    enddo
  enddo
!
! Build up the Cayley table and look for new operators at the same
! time:
!
  ini_nopsym = nopsym
  newop = .true.
  do while (newop)
    newop = .false.
    old_nopsym = nopsym
    do i = 1, old_nopsym
      do j = 1, old_nopsym
        if (optable(i,j) .le. 0) then
!         Get the product of operators i x j and determine if
!         this is a new or an old symmetry operator:
          do ii = 1, 3
            do jj = 1, 3
              xmat(ii,jj) = 0d0
              do kk = 1, 3
                xmat(ii,jj) = xmat(ii,jj) + opsym(i,ii,kk) * opsym(j,kk,jj)
              enddo
            enddo
          enddo
          inum = symnumber (xmat)
          if (inum .le. 0) then
            newop = .true.
            optable(i,j) = nopsym + 1
            call symopadd (xmat)
          else
            optable(i,j) = inum
          endif
        endif
      enddo
    enddo
  enddo
!
! Determine the operation inverses:
!
  do i = 1, nopsym
    j = 1
    found = .false.
    do while (.not.found .and. j.le.nopsym)
      kk = optable(i,j)
      if (optype(kk) .eq. opt_identity) then
        opinv(i) = j
        found = .true.
      else
        j = j + 1
      endif
    enddo
    if (.not.found) write (uout,605) i
  enddo
!
! Check for a set of generators:
!
  do i = 1, nopsym
    generated(i) = .false.
    opisgener(i) = .false.
  enddo
  do i = 1, nopsym
    if (.not.generated(i) .and. optype(i).ne.opt_identity) then
      opisgener(i) = .true.
    endif
    do j = 1, i
      kk = optable(i,j)
      generated(kk) = .true.
      kk = optable(j,i)
      generated(kk) = .true.
    enddo
  enddo
!
! Determine the classes of the symmetry operations:
! (Non-essential)
! A and B belong in the same class if there exist any operation X in
! the group such that: X^{-1} A X = B.
!
  do i = 1, nopsym
    opclas(i) = -1
  enddo
  nclas = 0
  do i = 1, nopsym
    if (opclas(i) .le. 0) then
!     This op belongs to a new class:
      nclas = nclas + 1
      if (nclas.gt.MCLAS) call error ('symclosure',  &
      'Too many classes. Increase MCLAS!', faterr)
      opclas(i) = nclas
      iclas = nclas
      ninclas(nclas) = 1
      iinclas(nclas,1) = i
    else
!     The op belongs to a class already found:
      iclas = opclas(i)
    endif
!   Apply similarity transforms to get all ops in the same class:
    do j = 1, nopsym
      jinv = opinv(j)
      iop = optable(optable(jinv, i), j)
      if (opclas(iop) .le. 0) then
        opclas(iop) = iclas
        ninclas(nclas) = ninclas(nclas) + 1
        if (ninclas(nclas).gt.MINCLAS) call error ('symclosure', &
        'Too many ops in a class. Increase MINCLAS!', faterr)
        iinclas(nclas,ninclas(nclas)) = iop
      else if (opclas(iop) .ne. iclas) then
!       This is an error: iop and i should be in the same class,
!       but they are not.
        write (uout,620) i, iop
      endif
    enddo
  enddo
 
  if (nopsym.gt.ini_nopsym .and. uout.ge.0) then
    write (uout,600)
    write (uout,602) nopsym-ini_nopsym
  endif
 
600 format (/                                                     &
    1x, '++SYMCLOSURE: Complete the group by multiplying the',    &
    1x, 'matrices already known.')

602 format (                                                      &
    1x, 'ALG(symclosure) ', i3,                                   &
    ' new operations found because of the group closure')

605 format (                                                      &
    1x, 'DBG(symclosure): Inverse not found for opsym ', i3)

620 format (                                                      &
    1x, 'WARNING(symclosure): The ops', 2i5, ' should be in the', &
    1x, 'same class but they are not. We continue anyway!')
 
210 format (a)
 
end subroutine symclosure
subroutine symdeter(uout, ax, ay, az, n)
!
! Determine linear, planar, 3d
!
  implicit none

  integer ::  uout, n
  real(kind=8) :: ax(n), ay(n), az(n)

  real(kind=8) :: an(n)
  real(kind=8) :: nor
  integer  :: i, j, k
  real(kind=8) :: sum1, sum2, p_esc, p_mixto, contri

  if (uout.gt.0) write (uout,800)
  do i = 1, n
     nor = sqrt(ax(i)*ax(i) + ay(i)*ay(i) + az(i)*az(i))
     if (nor > 1d-9) then
        an(i) = 1d0 / nor
     else
        an(i) = 1d0
     end if
  enddo
!
! Test for linear molecules:
! The normalized scalar product of the position vectors of each pair
! of atoms must be +1 or -1.
!
  sum1 = 0d0
  do i = 1, n
    do j = i+1, n
      p_esc = (ax(i)*ax(j)+ay(i)*ay(j)+az(i)*az(j)) * an(i)*an(j)
      contri = abs(abs(p_esc) - 1d0)
      sum1 = sum1 + contri
    enddo
  enddo
  if (uout.gt.0) write (uout,820) sum1, TOLdist
  if (sum1 .le. TOLdist) then
    mol_linear = .true.
    if (uout.gt.0) write (uout,810)
  endif
!
! Test for planar molecules:
! The triple (volume) product of the position vectors of each group
! of three different atoms must be zero.
!
  sum2 = 0d0
  do i = 1, n
    do j = i+1, n
      do k = j+1, n
        p_mixto = ax(i) * (ay(j)*az(k) - az(j)*ay(k)) &
                + ay(i) * (az(j)*ax(k) - ax(j)*az(k)) &
                + az(i) * (ax(j)*ay(k) - ay(j)*ax(k))
        sum2 = sum2 + p_mixto
      enddo
    enddo
  enddo
  if (uout.gt.0) write (uout,821) sum2, TOLdist
  if (sum2 .le. TOLdist) then
    mol_planar = .true.
    if (uout.gt.0) write (uout,815)
  endif

 800 format (/                                             &
     1x, '++MOLSYM: Determining the point group symmetry') 
 810 format (                                              &
     1x, 'Linear molecule detected')                       
 815 format (                                              &
     1x, 'Planar molecule detected')
 820 format (                                              &
     1x, 'Linearity test value (TOL): ', 1p, e10.3, '   (',e10.3,')')
 821 format (                                              &
     1x, 'Planarity test value (TOL): ', 1p, e10.3, '   (',e10.3,')')

end subroutine symdeter
subroutine symeigen(xmat, xval, xvect, xdet, ierror)
!
! symeigen - obtain the eigenvalues and eigenvectors of the
! orthogonal matrix xmat(,) by calling to the dgeev() routine in
! the lapack package. If xdet is +/-1, the routine uses this
! value as the determinant of xmat(,). Otherwise the determinant
! is obtained.
!
! The eigenvalues are sorted with xval(3) being the equivalent
! to the matrix determinant and xvect(:,3) being the normal vector
! corresponding to the rotation axis.
! The two other eigenvalues are meaningless (the imaginary part is
! not returned) but their eigenvectors define the plane normal to
! the rotation axis.
! The three eigenvectors are real and orthonormal.
!
  implicit none

  integer :: ierror
  real(kind=8), dimension(3,3) :: xmat
  real(kind=8), dimension(3,3) :: xvect
  real(kind=8), dimension(3) :: xval
  real(kind=8) :: xdet

  integer, parameter :: lwork = 102
  integer :: i, j, ii, ij, ik
  real(kind=8) :: dif, difmax, xnorm
  real(kind=8), dimension(3) :: wr, wi
  real(kind=8), dimension(3,3) :: ar, vr, vl
  real(kind=8), dimension(lwork) ::  work
!
! Get the determinant if it is needed:
!
  if (abs(abs(xdet)-1d0) .ge. TOLeigen) then
    xdet = xmat(1,1)*(xmat(2,2)*xmat(3,3) - xmat(2,3)*xmat(3,2)) &
         + xmat(1,2)*(xmat(2,3)*xmat(3,1) - xmat(2,1)*xmat(3,3)) &
         + xmat(1,3)*(xmat(2,1)*xmat(3,2) - xmat(2,2)*xmat(3,1))
  endif

  do i = 1, 3
    do j = 1, 3
      ar(i,j) = xmat(i,j)
    enddo
  enddo
  
  call dgeev ("N", "V", 3, ar, 3, wr, wi, vl, 3, vr, 3, work, lwork, ierror)
!
! Get the eigenvalue identical to the determinant and return it as
! the third one:
!
  ik = 0
  difmax = 1d40
  do i = 1, 3
    dif = abs(wr(i)-xdet) + abs(wi(i))
    if (dif .lt. difmax) then
      ik = i
      difmax = dif
    endif
  enddo
  ii = mod(ik,3) + 1
  ij = mod(ii,3) + 1
  ERReigen = max(ERReigen, difmax)
  xval(3) = wr(ik)
  xvect(1,3) = vr(1,ik)
  xvect(2,3) = vr(2,ik)
  xvect(3,3) = vr(3,ik)
  xval(1) = wr(ii)
  xvect(1,1) = vr(1,ii)
  xvect(2,1) = vr(2,ii)
  xvect(3,1) = vr(3,ii)
  xval(2) = wr(ij)
!
!.....Get explicitely the v2 = v3 x v1 eigenvector:
!
  xvect(1,2) = xvect(2,3) * xvect(3,1) - xvect(3,3) * xvect(2,1)
  xvect(2,2) = xvect(3,3) * xvect(1,1) - xvect(1,3) * xvect(3,1)
  xvect(3,2) = xvect(1,3) * xvect(2,1) - xvect(2,3) * xvect(1,1)
!
!.....Enforce the normalization:
!
  do i = 1, 3
    xnorm = xvect(1,i)**2 + xvect(2,i)**2 + xvect(3,i)**2
    xnorm = 1d0 / sqrt(xnorm)
    xvect(1,i) = xvect(1,i) * xnorm
    xvect(2,i) = xvect(2,i) * xnorm
    xvect(3,i) = xvect(3,i) * xnorm
  enddo
 
end subroutine symeigen
subroutine symeuler (matrix, euler)
!
! symeuler - get the (z,x,z) counterclockwise Euler angles that
! correspond to the rotation input "matrix".
!
  implicit none

  real(kind=8), dimension(3,3) :: matrix
  real(kind=8), dimension(3) :: euler

  if (abs(matrix(3,3)) .gt. 1d0) then
    euler(2) = acos(sign(1d0,matrix(3,3)))
  else
    euler(2) = acos(matrix(3,3))
  endif

  if (abs(abs(matrix(3,3))-1d0) .le. 1d-6) then
    euler(3) = 0d0
    if (matrix(3,3) .gt. 0d0) then
      euler(1) = atan2(matrix(2,1), matrix(1,1))
    else
      euler(1) = atan2(-matrix(2,1), matrix(1,1))
    endif
  else
    euler(1) = atan2(matrix(3,1), matrix(3,2))
    euler(3) = atan2(matrix(1,3), -matrix(2,3))
  endif

  if (euler(1) .lt. 0d0) euler(1) = euler(1) + dpi + dpi
  if (euler(3) .lt. 0d0) euler(3) = euler(3) + dpi + dpi

end subroutine symeuler
subroutine symimagic (uout)
!
! Symimagic - Check the elements in the symmetry matrices against
! the icosahedral magic numbers internally known by the routine.
!

  implicit none
 
  integer :: uout
 
  real(kind=8) :: x, ax, xbest, xdif, TOLxxx
  real(kind=8) :: xmat(3,3)
  integer  :: iop, i, j, k, imin, imax, imed
  integer  :: kk, kbest, iter, ip, iq, ir
  logical  :: found, doagain, generated(MOPSYM)
 
  integer  :: nmagic
  real(kind=8) :: magic(24)

  data nmagic /24/
  data (magic(i), i = 1, 24) /0.00000000000000000000d0, &
  0.05278640450004206072d0, 0.13819660112501051518d0,   &
  0.16245984811645316308d0, 0.26286555605956680302d0,   &
  0.27639320225002103036d0, 0.30901699437494742410d0,   &
  0.36180339887498948482d0, 0.42532540417601996609d0,   &
  0.44721359549995793928d0, 0.50000000000000000000d0,   &
  0.52573111211913360603d0, 0.58778525229247312917d0,   &
  0.63819660112501051518d0, 0.67082039324993690892d0,   &
  0.68819096023558676910d0, 0.72360679774997896964d0,   &
  0.80901699437494742410d0, 0.85065080835203993218d0,   &
  0.86180339887498948482d0, 0.89442719099991587856d0,   &
  0.94721359549995793928d0, 0.95105651629515357212d0,   &
  1.00000000000000000000d0 /

!
! 0.003842921 is the smallest difference between two magic numbers.
! Use a larger tolerance for the comparisons.
!
  TOLxxx = min(4d-3, TOLeqvm)
!
! Apply the magic cleaning to the group generators. Identity and
! inversion matrices are clean by construction.
!
  do iop = 1, nopsym
    if (optype(iop).eq.opt_identity .or. &
        optype(iop).eq.opt_inversion) then
      generated(iop) = .true.
!
    else if (opisgener(iop)) then
      generated(iop) = .true.
      do i = 1, 3
        do j = 1, 3
          x = opsym(iop,i,j)
          ax = abs(x)
          imin = 1
          imax = nmagic
          k = 0
          found = .false.
          do while (.not.found)
            k = k + 1
            imed = (imin+imax)/2
            if (abs(ax-magic(imed)) .le. TOLxxx) then
              opsym(iop,i,j) = sign(magic(imed), x)
              found = .true.
            else if ((imax-imin).lt.1 .or. k.gt.20) then
              write (uout,700) iop, i, j, ax, imed,imin,imax,k
              kbest = -1
              xbest = 1d30
              do kk = max(1,imed-3), min(nmagic,imed+3)
                xdif = abs(ax-magic(kk))
                if (xdif .lt. xbest) then
                  xbest = xdif
                  kbest = kk
                endif
              enddo
              if (kbest .lt. 0) then
                write (uout,705) iop, i, j, x
                call error ('symimagic', 'Magic number not found!', faterr)
              endif
              opsym(iop,i,j) = sign(magic(kbest), x)
              write (uout,710) iop,i,j, x,opsym(iop,i,j),xbest
              found = .true.
            else if (magic(imed) .gt. ax) then
              imax = imed - 1
            else
              imin = imed + 1
            endif
          enddo
        enddo
      enddo
!
    else
      generated(iop) = .false.
    endif
  enddo
!
! Regenerate the symmetry operations from the cleaned generators:
!
  doagain = .true.
  iter = 0
  do while (doagain .and. iter.lt.20)
    doagain = .false.
    iter = iter + 1
    do i = 1, nopsym
      if (generated(i)) then
        do j = 1, nopsym
          k = optable(i,j)
          if (generated(j) .and. (.not.generated(k))) then
!           This is a non regenerated opsym:
            doagain = .true.
            do ip = 1, 3
              do iq = 1, 3
                xmat(ip,iq) = 0d0
                do ir = 1, 3
                  xmat(ip,iq) = xmat(ip,iq) + opsym(i,ip,ir) * opsym(j,ir,iq)
                enddo
              enddo
            enddo
            do ip = 1, 3
              do iq = 1, 3
                opsym(k,ip,iq) = xmat(ip,iq)
              enddo
            enddo
            generated(k) = .true.
          endif
        enddo
      endif
    enddo
  enddo
  if (doagain) call error ('symimagic', &
  'The group was not regenerated in 20 iterations!', faterr)
!
!.....Regenerate the symmetry matrices properties:
!
  call symmatfill (xmat,1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0)
  call symtransform (xmat, xmat)
 
700 format (/                                                           &
    1x, '**SYMIMAGIC: Problem finding an icosahedral magic number'/     &
    1x, 'Debug data: ', 3i5, 1p, e17.9, 0p, 4i5/                        &
    1x, 'TROUBLE! The binary search has failed. We will use the',       &
    1x, 'magic number closest to the dirty unknown!')

705 format (                                                            &
    1x, 'PANIC! The best approach search has failed too:',              &
    1x, 'op(', i3, ',', i1, ',', i1, ') dirty value is ', 1p, e17.9)

710 format (                                                            &
    1x, 'Best approach result: op(', i3, ',', i1, ',', i1, ')',         &
    1x, 'dirty, clean, and diff. values ', 1p, 2e17.9, e11.3)

end subroutine symimagic
subroutine syminv (xmat, xinv, xdet, ierror)
!
! Get the inverse of a 3x3 nonsingular matrix.
! The input matrix is not modified.
!
  implicit none

  real(kind=8) :: xmat(3,3), xinv(3,3), xdet
  integer :: ierror

  integer :: i, j
  real(kind=8) :: adj(3,3), det, ddet

  ierror = 0
  adj(1,1) = xmat(2,2)*xmat(3,3) - xmat(2,3)*xmat(3,2)
  adj(1,2) = xmat(2,3)*xmat(3,1) - xmat(2,1)*xmat(3,3)
  adj(1,3) = xmat(3,2)*xmat(2,1) - xmat(2,2)*xmat(3,1)
  det = xmat(1,1)*adj(1,1) + xmat(1,2)*adj(1,2) + xmat(1,3)*adj(1,3)
  if (abs(det) .le. TOLsng) then
     ierror = -1
     return
  endif
  adj(2,1) = xmat(1,3)*xmat(3,2) - xmat(1,2)*xmat(3,3)
  adj(2,2) = xmat(1,1)*xmat(3,3) - xmat(1,3)*xmat(3,1)
  adj(2,3) = xmat(1,2)*xmat(3,1) - xmat(1,1)*xmat(3,2)
  adj(3,1) = xmat(1,2)*xmat(2,3) - xmat(1,3)*xmat(2,2)
  adj(3,2) = xmat(1,3)*xmat(2,1) - xmat(1,1)*xmat(2,3)
  adj(3,3) = xmat(1,1)*xmat(2,2) - xmat(1,2)*xmat(2,1)
  ddet = 1d0/det
  do i = 1, 3
     do j = 1, 3
        xinv(i,j) = adj(j,i) * ddet
     enddo
  enddo
  xdet = det
  return
 
end subroutine syminv
subroutine symmatfill (xmat, a11, a12, a13, a21, a22, a23, a31, a32, a33)
!
! symmatfill - fill in the xmat matrix.
!
  implicit none

  real(kind=8), dimension(3,3) :: xmat
  real(kind=8) :: a11, a12, a13, a21, a22, a23, a31, a32, a33

  xmat(1,1) = a11
  xmat(1,2) = a12
  xmat(1,3) = a13
  xmat(2,1) = a21
  xmat(2,2) = a22
  xmat(2,3) = a23
  xmat(3,1) = a31
  xmat(3,2) = a32
  xmat(3,3) = a33

end subroutine symmatfill
subroutine symmatprod(a, b, c)
!
! symmatprod - perform the "C = A B" matrix multiplication between
! 3x3 matrices. The final matrix "C" can be equal to either "A" or
! "B", so the product will be stored on a secondary matrix before
! being moved to "C".
!
  implicit none

  real(kind=8), dimension(3,3) :: a, b, c, dummy

  integer :: i, j, k

  dummy = matmul(a,b)
  c(:,:) = dummy(:,:)
 
end subroutine symmatprod
subroutine symnoise (uout, ax, ay, az, n)
!
!.....symnoise - Estimate the average deviation of atoms with respect
!     to the perfectly symmetric structure.
!
  implicit none
 
  integer :: uout, n
  real(kind=8), dimension(n) :: ax, ay, az
 
  integer :: i, j, j1, k, iat, iorb, iop, kbest, kat
  real(kind=8) :: R1, R2, R3
  real(kind=8) :: dorb, dj, alpha, xxx, yyy, zzz, dbest, TOL2
  logical :: done

!
!.....All atoms in a orbit should have the same distance to the
!     center of mass. Get a simple noise measurement from checking
!     this assumption:
!
  R1 = 0d0
  R2 = 0d0
  R3 = 0d0
  k = 0
  do i = 1, norbit
    dorb = orbdis(i)
    do j = 1, natorb(i)
      k = k + 1
      j1 = iatorb(i,j)
      dj = sqrt(ax(j1)*ax(j1) + ay(j1)*ay(j1) + az(j1)*az(j1))
      R1 = R1 + (dj - dorb)
      R2 = R2 + (dj - dorb)**2
      R3 = R3 + abs(dj - dorb)
    enddo
  enddo
!  R1 = R1 / n
!  R2 = R2 / n
  R1 = R1 / k
  R2 = R2 / k
  R3 = R3 / k
  write (uout,605) R3, sqrt(R2 - R1*R1)

!
! Are the rotation angles of symmetry matrices exact?
!
  R1 = 0d0
  R2 = 0d0
  R3 = 0d0
  do i = 1, nopsym
    alpha = opm(i) * 2d0 * dpi / oporder(i)
    if (opangle(i) .lt. 0d0) alpha = -alpha
    R1 = R1 + (alpha - opangle(i))
    R2 = R2 + (alpha - opangle(i))**2
    R3 = R3 + abs(alpha - opangle(i))
  enddo
  R1 = R1 / nopsym
  R2 = R2 / nopsym
  R3 = R3 / nopsym
  write (uout,610) R3, sqrt(R2 - R1*R1)

!
! A desperate movement:
!   Produce all copies of all atoms operated by all symmetry
!   matrices. Get the average distance between every atom and all
!   of its copies:
!
  TOL2 = TOLdist * TOLdist
  R1 = 0d0
  kat = 0
  do iop = 1, nopsym
    if(optype(iop) .ne. opt_identity) then
      do iorb = 1, norbit
        do i = 1, natorb(iorb)
          iat = iatorb(iorb,i)
          xxx = opsym(iop,1,1)*ax(iat) + opsym(iop,1,2)*ay(iat) &
                + opsym(iop,1,3)*az(iat)
          yyy = opsym(iop,2,1)*ax(iat) + opsym(iop,2,2)*ay(iat) &
                + opsym(iop,2,3)*az(iat)
          zzz = opsym(iop,3,1)*ax(iat) + opsym(iop,3,2)*ay(iat) &
                + opsym(iop,3,3)*az(iat)
          kbest = 0
          dbest = 1d30
          j = 1
          done = .false.
          do while(.not.done .and. j.le.natorb(iorb))
            k = iatorb(iorb,j)
            dj = (xxx-ax(k))**2+(yyy-ay(k))**2+(zzz-az(k))**2
            if(dj.le.dbest) then
              kbest = k
              dbest = dj
            endif
            if(dj .le. TOL2) then
              done = .true.
            else
              j = j + 1
            endif
          enddo
          if(done) then
            kat = kat + 1
            R1 = R1 + sqrt(dj)
          else
            write (uout,650) iat, iop, kbest, dbest
          endif
        enddo
      enddo
    endif
  enddo
  R1 = R1 / kat
  write (uout,615) R1, n*(nopsym-1)-kat

605 format (/                &
    1x, '++SYMNOISE:'/       &
    1x, 'Error on the orbit radii. Average and sdev: ', 1p, 2e15.6)

610 format (1x, 'Error on the rot. angles. Average and sdev: ', 1p, 2e15.6)

615 format (1x, 'Distance between sym/at images. Average reject: ', 1p, e15.6,0p, i8)

650 format (1x, 'DBG(symnoise): sym/at image not found', 3i5, 1p, e15.6)

end subroutine symnoise
integer function symnumber (xmat)
!
!.....symnumber - get the number of xmat in the symmetry operator list.
!     A -1 value will be returned if xmat is not in the list.
!
  implicit none
 
  real(kind=8), dimension(3,3) :: xmat
 
  integer :: i, j, k
  real(kind=8) :: diff
 
  symnumber = -1
  do i = 1, nopsym
    diff = 0d0
    do j = 1, 3
      do k = 1, 3
        diff = diff + abs(xmat(j,k) - opsym(i,j,k))
      enddo
    enddo
    if(diff .le. TOLeqvm) then
      symnumber = i
      return
    endif
  enddo
 
end function symnumber
subroutine symopadd (xmat)
!
! Symopadd - if xmat is yet unknown add it to the list of sym ops.
!

  implicit none
 
  real(kind=8) :: xmat(3,3)
  integer, parameter :: MOP = 20
  integer  :: i, j, k
  real(kind=8) :: diff
 
  integer  :: ierror, ii, ij, ik, ibest
  logical  :: found
  real(kind=8) :: xtrace, xdet, xnorm, xnm, ang1, angle
  real(kind=8) :: err, besterr
  real(kind=8) :: xinv(3,3), xval(3), xvect(3,3), xeul(3)
  character*(40) :: fmt
 
  do i = 1, nopsym
    diff = 0d0
    do j = 1, 3
      do k = 1, 3
        diff = diff + abs(xmat(j,k) - opsym(i,j,k))
      enddo
    enddo
    if (diff .le. TOLeqvm) return  ! This is not a new operator
  enddo

  nopsym = nopsym + 1

  if (nopsym .gt. MOPSYM) call error ('symopadd', &
  'Too many sym operators. Increase MOPSYM', faterr)

  do j = 1, 3
    do k = 1, 3
      opsym(nopsym,j,k) = xmat(j,k)
    enddo
  enddo
 
!
! Get the properties of the symmetry matrix:
!
  xtrace = 0d0
  do i = 1, 3
    xtrace = xtrace + xmat(i,i)
  enddo
  call syminv (xmat, xinv, xdet, ierror)
  if (ierror .lt. 0) call error ('symopadd', &
  'Operator has a singular matrix!!!', faterr)
!
! Get the Euler angles for the rotation:
!
  call symeuler (xmat, xeul)
  opeuler(nopsym,1) = xeul(1)
  opeuler(nopsym,2) = xeul(2)
  opeuler(nopsym,3) = xeul(3)
!
! Is this a proper or improper rotation?
!
  if (abs(abs(xdet)-1d0) .gt. TOLisint) then
    call error ('symopadd', 'determinant is not +-1', warning)
  else if (xdet .gt. 0d0) then
    opproper(nopsym) = .true.
  else
    opproper(nopsym) = .false.
  endif
!
! Get the rotation axis except for E and i, for which it is
! fully degenerated:
!
  if (abs(abs(xtrace)-3d0) .le. TOLisint) then
    if (xtrace .gt. 0d0) then
      opsymbol(nopsym) = 'E'
      opaxis(nopsym,1) = 0d0
      opaxis(nopsym,2) = 0d0
      opaxis(nopsym,3) = 1d0
      oporder(nopsym) = 1
      opm(nopsym) = 1
      opangle(nopsym) = dpi + dpi
      optype(nopsym) = opt_identity
    else
      opsymbol(nopsym) = 'i'
      opaxis(nopsym,1) = 0d0
      opaxis(nopsym,2) = 0d0
      opaxis(nopsym,3) = 1d0
      oporder(nopsym) = 2
      opm(nopsym) = 1
      opangle(nopsym) = dpi
      optype(nopsym) = opt_inversion
    endif
  else
! Get the rotation angle and the n and m of the C_n^m or S_n^m
! symmetry operator:
    ang1 = 0.5d0*(xtrace-xdet)
    if (abs(ang1) .gt. 1d0) ang1 = sign(1d0, ang1)
    angle = acos(ang1)
    if (abs(angle) .le. TOLnull) then
      xnm = 1d0
      angle = dpi + dpi
    else
      xnm = 2d0 * dpi / angle
    endif
    opangle(nopsym) = angle
    ii = 0
    found = .false.
    do while (.not. found .and. ii.le.MOP)
      ii = ii + 1
      found = abs(xnm*ii - nint(xnm*ii)) .le. TOLisint
    enddo
    if (found) then
      oporder(nopsym) = nint(xnm*ii)
      opm(nopsym) = ii
    else
!     call error ('symopadd', 'Using best approx order', warning)
      ibest = 1
      besterr = abs(xnm - nint(xnm))
      do ii = 2, MOP
        err = abs(xnm*ii - nint(xnm*ii))
        if (err .le. besterr) then
          besterr = err
          ibest = ii
        endif
      enddo
      oporder(nopsym) = nint(xnm*ibest)
      opm(nopsym) = ibest
    endif
    write (fmt,200) 1+int(log10(1d0*oporder(nopsym))) &
                    , 1+int(log10(1d0*opm(nopsym)))
    call symeigen (xmat, xval, xvect, xdet, ierror)
    if (ierror .ne. 0) then
      call error ('symopadd', 'Trouble finding the rotation axis', warning)
    endif
    opaxis(nopsym,1) = xvect(1,3)
    opaxis(nopsym,2) = xvect(2,3)
    opaxis(nopsym,3) = xvect(3,3)
    if (xdet .gt. 0d0) then
      write (opsymbol(nopsym),fmt) 'C', oporder(nopsym), opm(nopsym)
      optype(nopsym) = opt_rotation
    else if (abs(xtrace - 1d0) .le. TOLisint) then
      write (opsymbol(nopsym),210) 'sigma'
      optype(nopsym) = opt_sigma
    else if (oporder(nopsym).eq.1 .and. opm(nopsym).eq.1) then
      write (opsymbol(nopsym),210) 'sigma'
      optype(nopsym) = opt_sigma
    else
      write (opsymbol(nopsym),fmt) 'S', oporder(nopsym), opm(nopsym)
      optype(nopsym) = opt_imp_rotation
    endif
  endif
 
200 format ('(a, "_", i', i2, ', "^", i', i2, ')')
210 format (a)
 
end subroutine symopadd

logical function symopchk(xmat, ax, ay, az, atgroup, nmol)
!
!.....symopchk - check if xmat is a sym operator.
!
  implicit none

  integer :: nmol
  integer, dimension(nmol) :: atgroup
  real(kind=8), dimension(nmol) :: ax, ay, az
  real(kind=8), dimension(3,3) :: xmat

  integer :: i, j, k
  real(kind=8) :: ximg, yimg, zimg, diff
  logical :: found
!
!.....First check: The matrix must be orthogonal:
!
  symopchk = .false.

!
!.....If the number of atoms is large, perhaps is better to check
!     for repeated operations of symmetry before testing if it
!     really transforms the molecule into itself:
!
  if(nmol .gt. nopsym+nopsym) then
    do i = 1, nopsym
      diff = 0d0
      do j = 1, 3
        do k = 1, 3
          diff = diff + abs(xmat(j,k) - opsym(i,j,k))
        enddo
      enddo
      if(diff .le. TOLeqvm) return  ! This is not a new operator
    enddo
  endif

!
!.....Transform every atom and check for the symmetric image:
!
  do i = 1, nmol
    ximg = xmat(1,1)*ax(i) + xmat(1,2)*ay(i) + xmat(1,3)*az(i)
    yimg = xmat(2,1)*ax(i) + xmat(2,2)*ay(i) + xmat(2,3)*az(i)
    zimg = xmat(3,1)*ax(i) + xmat(3,2)*ay(i) + xmat(3,3)*az(i)
    found = .false.
    j = 1
    do while(.not.found .and. j.le.nmol)
      if(atgroup(i).eq.atgroup(j)) then
        found = abs(ximg-ax(j))+abs(yimg-ay(j))+abs(zimg-az(j)) .le. TOLdist
      endif
      j = j + 1
    enddo
    if(.not.found) return
  enddo
  symopchk = .true.
 
end function symopchk 

subroutine symorb (uout, ax, ay, az, atZmol, nmol)
!
! symorb - classify all atoms into orbits. All atoms in an orbit
! have the same atomic number and the same distance to the center.
!
  implicit none
 
  integer :: uout
  integer :: nmol  
  integer, dimension(nmol) :: atZmol
  real(kind=8), dimension(nmol) :: ax, ay, az
 
  real(kind=8) :: dis2
  integer :: i, j, iorb

!
! Classify all atoms into orbits. All atoms in an orbit have the
! same atomic number and the same distance to the center of mass:
!
  norbit = 0
  molradius = 0d0
  do i = 1, nmol
    dis2 = sqrt(ax(i)*ax(i) + ay(i)*ay(i) + az(i)*az(i))
    if (dis2 .gt. molradius) molradius = dis2
    iorb = -1
    j = 1
    do while (j.le.norbit .and. iorb.lt.0)
      if (abs(dis2-orbdis(j)) .le. TOLdist .and. atZmol(i) .eq. orbZ(j)) then
        iorb = j
      endif
      j = j + 1
    enddo
    if (iorb .lt. 0) then
      norbit = norbit + 1
      if (norbit .gt. MORBIT) call error ('symorb', &
         'Too many orbits. Increase MORBIT!', faterr)
      orbdis(norbit) = dis2
      orbZ(norbit) = atZmol(i)
      natorb(norbit) = 0
      iorb = norbit
    endif
    orbmol(i) = iorb
    natorb(iorb) = natorb(iorb) + 1
    if (natorb(iorb) .gt. MATORB) call error ('symorb', &
       'Too many atoms in an orbit. Increase MATORB!', faterr)
     iatorb(iorb,natorb(iorb)) = i
  enddo

end subroutine symorb
subroutine symprint(uout,verbose)
!
! symprint - print the symmetry operators.
!
  implicit none

  integer :: uout
  logical :: verbose
  character(len=20) :: group

  integer :: i, j, k

  if (uout.lt.0) return
  write (uout,700) point_group, nopsym
  if (verbose) then
    do i = 1, nopsym
      write (uout,705) i, opsymbol(i), opangle(i) * 180d0 / dpi,       &  
                       (opaxis(i,k), k=1,3), opisgener(i), opinv(i),  &
                       opclas(i), (opeuler(i,k) * 180d0 / dpi, k=1,3), &
                       ((opsym(i,j,k), k=1,3), j=1,3)
    enddo
  endif
  write (uout,780) nclas
  do i = 1, nclas
     write (uout,785) i, ninclas(i), opsymbol(iinclas(i,1)), &
                      (iinclas(i,j), j = 1, ninclas(i))
  enddo
  if (verbose) then
    write (uout,750)
    do i = 1, nopsym
      write (uout,755) (optable(i,j), j=1,nopsym)
    enddo
  endif

700 format (/                                        &
    1x, '++SYMPRINT: List of symmetry operators:'/   &
    1x, 'Point group: ', a/                          &
    1x, 'Number of operators found: ', i5)

705 format (/                                        &
    1x, 'Operator #', i4, ' :'/                      &
    1x, 'Symbol: ', a/                               &
    1x, 'Rotation angle:   ', f15.9, ' degrees'/     &
    1x, 'Rotation axis:    ', 3f15.9/                &
    1x, 'Is a generator?   ', l6/                    &
    1x, 'Inverse operation:', i6/                    &
    1x, 'Operations class: ', i6/                    &
    1x, 'Euler zxz angles: ', 3f15.9, 3x, 'degrees'/ &
    1x, 'Operator 3x3 matrix (xyz representation):'/ &
    (1x, 3f15.9))

750 format (/                                        &
    1x, 'Multiplication (Cayley) table:')

755 format (120i4)

780 format (/                                        &
    1x, 'Number of symmetry operation classes:', i5/ &
    1x, 'Class...#op....Repres..Operations...')      

785 format (2i6, 2x, a10, (20i4))

end subroutine symprint
subroutine symproject (uout, ax, ay, az, n)
!
! Symproject - Clean the atomic coordinates by projecting the atoms
! quite close to the symmetry elements onto those symmetry elements.
! This routine is suppossed to work after sympurify() has cleaned
! the symmetry matrices, but the algorithm do not depends on that.
!

  implicit none
 
  integer  :: n
  integer  :: uout
  real(kind=8) :: ax(n), ay(n), az(n)
 
  integer  :: iop, iorb, i, iat, k, kmin, ntrouble
  real(kind=8) :: dd, xproj, aproj(3), dmin, dtrouble, dproj
  real(kind=8) :: ximg, yimg, zimg, xx, yy, zz
  real(kind=8) :: axx(n), ayy(n), azz(n)
  logical  :: cleaned(n), found, trouble
!
! Keep a copy of the original coordinates:
!
  do i = 1, n
    axx(i) = ax(i)
    ayy(i) = ay(i)
    azz(i) = az(i)
  enddo
!
! Project first the first atom in each atomic orbit.
!
  if (uout.ge.0) write (uout,600)
  do iop = 1, nopsym
!
    if (optype(iop).eq.opt_inversion) then
!     Move atoms quite close to the inversion center:
      do iorb = 1, norbit
        do i = 1, natorb(iorb)
          iat = iatorb(iorb,i)
          dd = sqrt(ax(iat)**2 + ay(iat)**2 + az(iat)**2)
          if (dd .le. TOLdist) then
            if (uout.ge.0) write (uout,605) iat, dd, ax(iat), ay(iat), az(iat)
            ax(iat) = 0d0
            ay(iat) = 0d0
            az(iat) = 0d0
          endif
        enddo
      enddo
!
    else if (optype(iop).eq.opt_rotation .or.       &
             optype(iop).eq.opt_imp_rotation) then
!     Atoms too close to a rotation axis will be projected
!     onto the axis:
      do iorb = 1, norbit
        do i = 1, natorb(iorb)
          iat = iatorb(iorb,i)
          xproj = ax(iat) * opaxis(iop,1)+ ay(iat) * opaxis(iop,2) + &
                  az(iat) * opaxis(iop,3)
          aproj(1) = xproj * opaxis(iop,1)
          aproj(2) = xproj * opaxis(iop,2)
          aproj(3) = xproj * opaxis(iop,3)
          xx = aproj(1) - ax(iat)
          yy = aproj(2) - ay(iat)
          zz = aproj(3) - az(iat)
          dd = sqrt(xx*xx + yy*yy + zz*zz)
          if (dd .le. TOLdist) then
            if (uout.ge.0) write (uout,610) iat, iop, dd, ax(iat), ay(iat), &
                                    az(iat), aproj(1), aproj(2), aproj(3)
            ax(iat) = aproj(1)
            ay(iat) = aproj(2)
            az(iat) = aproj(3)
          endif
        enddo
      enddo
!
    else if (optype(iop).eq.opt_sigma) then
!     Atoms too close to a symmetry plane will be projected
!     onto the plane:
      do iorb = 1, norbit
        do i = 1, natorb(iorb)
          iat = iatorb(iorb,i)
          xproj = ax(iat) * opaxis(iop,1) + ay(iat) * opaxis(iop,2) + az(iat) &
                  * opaxis(iop,3)
          if (abs(xproj) .le. TOLdist) then
            aproj(1) = ax(iat) - xproj * opaxis(iop,1)
            aproj(2) = ay(iat) - xproj * opaxis(iop,2)
            aproj(3) = az(iat) - xproj * opaxis(iop,3)
            if (uout.ge.0) write (uout,615) iat, iop, xproj, ax(iat), ay(iat), &
                                            az(iat), aproj(1), aproj(2), aproj(3)
            ax(iat) = aproj(1)
            ay(iat) = aproj(2)
            az(iat) = aproj(3)
          endif
        enddo
      enddo
!
    endif
!
  enddo
!
! Regenerate each orbit by applying the symmetry matrices to the
! first orbit atom.
!
  do i = 1, n
    cleaned(i) = .false.
  enddo
  ntrouble = 0
  dtrouble = 0d0
  trouble = .false.
 
  do iorb = 1, norbit
    do i = 1, natorb(iorb)
      iat = iatorb(iorb,i)
      cleaned(iat) = .true.
      do iop = 2, nopsym
!       Get the image of the atom by the symmetry operation:
        ximg = opsym(iop,1,1)*ax(iat) + opsym(iop,1,2)*ay(iat) &
               + opsym(iop,1,3)*az(iat)
        yimg = opsym(iop,2,1)*ax(iat) + opsym(iop,2,2)*ay(iat) &
               + opsym(iop,2,3)*az(iat)
        zimg = opsym(iop,3,1)*ax(iat) + opsym(iop,3,2)*ay(iat) &
               + opsym(iop,3,3)*az(iat)
!       Find this image atom or the closest approximation:
        found = .false.
        dmin = 1d30
        kmin = -1
        k = 1
        do while (.not.found .and. k.le.n)
          dd = (ximg-ax(k))**2+(yimg-ay(k))**2+(zimg-az(k))**2
          dd = sqrt(dd)
          if (dd .le. TOLdist) then
            found = .true.
            dmin = dd
            kmin = k
          else if (dd .lt. dmin) then
            dmin = dd
            kmin = k
          endif
          k = k + 1
        enddo
        if (.not.found) then
          ntrouble = ntrouble + 1
          dtrouble = max(dtrouble, dmin)
          if (.not.cleaned(kmin)) trouble = .true.
        else if (.not.cleaned(kmin)) then
          cleaned(kmin) = .true.
          dtrouble = max(dtrouble, dmin)
          ax(kmin) = ximg
          ay(kmin) = yimg
          az(kmin) = zimg
        endif
      enddo
    enddo
  enddo

  if (uout.ge.0 .or. ntrouble.gt.0) write (uout,630) ntrouble, dtrouble

  if (trouble) then
!   Check that all atoms have been cleaned:
    ntrouble = 0
    do i = 1, n
      if (.not.cleaned(i)) then
        ntrouble = ntrouble + 1
        write (uout,635) i
      endif
    enddo
    if (ntrouble.gt.0) then
      write (uout,636) ntrouble
      call error ('symproject', 'Error cleaning atoms', faterr)
    endif
  endif
!
! Write down the cleaned coordinates:
!
  if (uout.lt.0) return
  write (uout,640)
  dproj = 0d0
  do i = 1, n
    xx = ax(i)-axx(i)
    yy = ay(i)-ayy(i)
    zz = az(i)-azz(i)
    dproj = dproj + sqrt(xx*xx + yy*yy + zz*zz)
    write (uout,645) i, ax(i), ay(i), az(i), ax(i)-axx(i), ay(i)-ayy(i), &
                     az(i)-azz(i)
  enddo
  dproj = dproj / n
  write (uout,650) dproj
    
600 format (/                                                           &
    1x, '++ SYMPROJECT: Project onto the symmetry elements atoms',      &
    1x, 'too close to them.')

605 format (                                                            &
    1x, 'ALG(symproject) Move atom', i4, ' to the origin.',             &
    1x, 'Distance:', 1p, e12.4, 0p, 'Original coordinates: ', 3f15.9)

610 format (                                                            &
    1x, 'ALG(symproject) Move atom', i4, ' to the rotation axis', i4,   &
    1x, 'Distance:', 1p, e12.4, 0p,                                     &
    1x, 'Old and new coordinates: '/ (3f15.9))

615 format (                                                            &
    1x, 'ALG(symproject) Move atom', i4, ' to the mirror plane ', i4,   &
    1x, 'Distance:', 1p, e12.4, 0p,                                     &
    1x, 'Old and new coordinates: '/ (3f15.9))

630 format (                                                            &
    1x, 'ALG(symproject): Atomic images not found and max error ->',    &
    1x, i6, 1p, e12.4)

635 format (                                                            &
    1x, 'PROBLEM(symproject): Atom ', i4, ' was not cleaned!')

636 format (                                                            &
    1x, 'PANIC(symproject): Number of atoms not cleaned: ', i4)

640 format (/                                                           &
    1x, 'ALG(symproject): Cleaned molecular coordinates.'/              &
    1x, '-i------new-x----------new-y----------new-z------',            &
            '---dif-x-----dif-y-----dif-z--')

645 format (i4, 3f15.9, 1p, 3e10.2)

650 format (/                                                           &
    1x, 'ALG(symproject): Average displacement of atoms: ', 1p,e15.6/   &
    1x, 'ALG(symproject): end of cleaning ok!')
 
end subroutine symproject
subroutine sympurify (uout)
!
!.....sympurify - purify the symmetry operators. Produce the group from
!     a clean set of generators. This routine assumes a particular
!     orientation for the molecule and it should follow a successful
!     calling to symorient().
!
  implicit none
 
  integer :: uout
 
  real(kind=8), dimension(3,3) :: xmat
  real(kind=8) :: ze, on, c, s, theta, phi, iphi, aa, bb, xx, yy, zz
 
  if(uout.ge.0) write (uout,600)
!
!.....Restart the symmetry group. Keep the identity matrix:
!
  if(point_group(1:1) .ne. 'I') then
!    nopsym = 0
!    call symmatfill (xmat, 1d0,0d0,0d0, 0d0,1d0,0d0, 0d0,0d0,1d0)
!    call symopadd (xmat)
    nopsym = 1
    ze = 0d0
    on = 1d0
  endif

! 
!.....Special orientation for the different symmetry groups:
!
  if(point_group(1:1) .eq. 'I') then
!
!........Icosahedral groups.
!        The icosahedral groups (either I or Ih) are quite special
!        among the finite point groups. There is not a single way for
!        defining a special orientation and a corresponding set of
!        generators. For instance, if we use a C_5^1 axis to define z,
!        there are five equivalent but different to select a C_3^1 axis
!        to fix the xz plane. In the simplest scheme, we can use three
!        mutually perpendicular C_2^1 axes to define {x,y,z}, but we
!        still face two equivalent but different ways to orient the
!        group, i.e. two forms to select an additional C_5^1 or C_3^1
!        axis. As a consequence of this ambiguity, it is far from being
!        easy to ensure that a set of generators matches the orientation
!        produced in symorient(). Therefore, we turn to use a specially
!        designed method for the icosahedral groups. The method is based
!        on the property that all the elements in the icosahedral
!        symmetry matrices are members of a small set of "magic numbers"
!        when the group is appropriately oriented. We only have to
!        convert each element into the closest magic number. The natural
!        orientation of the group is conserved in this way.
!
    call symimagic (uout)
    return

  else if(point_group(1:1) .eq. 'O') then
!
!........Octahedral groups.
!        Orientation: The first two C_4^1 axes determine z and x.
!        Generators: C_4^1(z), C_3^1(111) for O; add i for Oh.
!
!........i:
    if(point_group(1:2) .eq. 'Oh') then
      call symmatfill (xmat, -on,ze,ze, ze,-on,ze, ze,ze,-on)
      call symopadd (xmat)
    endif
!........C_4^1(z):
    call symmatfill (xmat, ze,on,ze, -on,ze,ze, ze,ze,on)
    call symopadd (xmat)
!........C_3^1(111):
    call symmatfill (xmat, ze,on,ze, ze,ze,on, on,ze,ze)
    call symopadd (xmat)

  else if(point_group(1:1) .eq. 'T') then
!........Tetrahedral groups.
!        Orientation: The first two C_2^1 axes determine z and x
!        Generators (T): C_3^1(111) and C_2^1(z).
!        Generators (Td): C_3^1(111) and S_4^3(z).
!        Generators (Th): C_3^1(111), C_2^1(z), and i.
!
!........C_3^1(111):
    call symmatfill (xmat, ze,on,ze, ze,ze,on, on,ze,ze)
    call symopadd (xmat)
    if(point_group(1:2) .eq. 'Td') then
!...........S_4^3(z):
      call symmatfill (xmat, ze,-on,ze, on,ze,ze, ze,ze,-on)
      call symopadd (xmat)
    else
!...........C_2^1(z):
      call symmatfill (xmat, -on,ze,ze, ze,-on,ze, ze,ze,on)
      call symopadd (xmat)
      if(point_group(1:2) .eq. 'Th') then
!..............i:
        call symmatfill (xmat, -on,ze,ze, ze,-on,ze, ze,ze,-on)
        call symopadd (xmat)
      endif
    endif
 
  else if(point_group(1:2) .eq. 'D2') then
!........Dihedral D2, D2h, and D2d groups.
!        The main axis determines z.
!        The next non-collinear one determines y.
!        Generators (D2): C_2^1(z) and C_2^1(x).
!        Generators (D2h): gen. of D2 plus sigma_h.
!        Generators (D2d): gen. of D2 plus sigma_d.
!
!........C_2^1(z):
    call symmatfill (xmat, -on,ze,ze, ze,-on,ze, ze,ze,on)
    call symopadd (xmat)
!........C_2^1(x):
    call symmatfill (xmat, on,ze,ze, ze,-on,ze, ze,ze,-on)
    call symopadd (xmat)
    if(point_group(1:3) .eq. 'D2h') then
!...........sigma_h:
      call symmatfill (xmat, on,ze,ze, ze,on,ze, ze,ze,-on)
      call symopadd (xmat)
    else if(point_group(1:3) .eq. 'D2d') then
!...........sigma_d:
      call symmatfill (xmat, ze,on,ze, on,ze,ze, ze,ze,on)
      call symopadd (xmat)
    endif
      
  else if(point_group(1:1) .eq. 'D') then
!........Dihedral Dn, Dnh, and Dnd groups with n>2.
!        The main axis fixes z. Any of the two perpendicular C2 gives x.
!        Generators (Dn): C_n^1(z) and C_2^1(x).
!        Generators (Dnh): gen. of Dn plus sigma_h.
!        Generators (Dnd): gen. of Dn plus sigma_d.
!
!........C_n^1(z):
    theta = 2d0 * dpi / opmainord
    c = cos(theta)
    s = sin(theta)
    call symmatfill (xmat, c,-s,ze, s,c,ze, ze,ze,on)
    call symopadd (xmat)
!........C_2^1(x):
    call symmatfill (xmat, on,ze,ze, ze,-on,ze, ze,ze,-on)
    call symopadd (xmat)
    if(point_g2(1:3) .eq. 'Dnh') then
!...........sigma_h:
      call symmatfill (xmat, on,ze,ze, ze,on,ze, ze,ze,-on)
      call symopadd (xmat)
    else if(point_g2(1:3) .eq. 'Dnd') then
!...........sigma_d:
      call symmatfill (xmat, ze,on,ze, on,ze,ze, ze,ze,on)
      call symopadd (xmat)
    endif
 
  else if (point_g2(1:2) .eq. 'Cn') then
!........Use the single main axis as z.
!        Generators (Cn); C_n^1(z).
!        Generators (Cnh); C_n^1(z) and sigma_h.
!        Generators (Cnv); C_n^1(z) and sigma_v(xz).
!
!
!........C_n^1(z):
    theta = 2d0 * dpi / opmainord
    c = cos(theta)
    s = sin(theta)
    call symmatfill (xmat, c,-s,ze, s,c,ze, ze,ze,on)
    call symopadd (xmat)
    if(point_g2(1:3) .eq. 'Cnv') then
!...........sigma_v:
      call symmatfill (xmat, -on,ze,ze, ze,on,ze, ze,ze,on)
      call symopadd (xmat)
    else if (point_g2(1:3) .eq. 'Cnh') then
!...........sigma_h:
      call symmatfill (xmat, on,ze,ze, ze,on,ze, ze,ze,-on)
      call symopadd (xmat)
    endif
 
  else if (point_g2(1:3) .eq. 'S2n') then
!........Use the single main axis as z.
!        Generators (S2n): S_{2n}^1(z).
!
!........S_{2n}^1(z):
    theta = 2d0 * dpi / opmainord
    c = cos(theta)
    s = sin(theta)
    call symmatfill (xmat, c,-s,ze, s,c,ze, ze,ze,-on)
    call symopadd (xmat)
 
  else if (point_group(1:2) .eq. 'Cs') then
!........The normal to sigma_h defines z.
!        Generators: sigma_h.
!
!........sigma_h:
    call symmatfill (xmat, on,ze,ze, ze,on,ze, ze,ze,-on)
    call symopadd (xmat)
 
  else if (point_group(1:2) .eq. 'Ci') then
!.......................................................................
!........Generators: i.
!
!........i:
    call symmatfill (xmat, -on,ze,ze, ze,-on,ze, ze,ze,-on)
    call symopadd (xmat)
   
  else if (point_group(1:2) .eq. 'C1') then
!.......................................................................
!........Use the inertia eigenvectors to orient the molecule:

  else
!........C_inf_v or D_inf_h groups.
!........As a practical consideration, change the C_inf into a
!        C(inf_order)
!        Generators (C_inf_v): C_{inf_order}^1(z) and sigma_v(plane x).
!        Generators (D_inf_h): those of C_inf_v plus sigma_h.
!
!........C_{inf}^1(z):
    theta = 2d0 * dpi / inf_order
    c = cos(theta)
    s = sin(theta)
    call symmatfill (xmat, c,-s,ze, s,c,ze, ze,ze,on)
    call symopadd (xmat)
!........sigma_v (plane xz):
    call symmatfill (xmat, on,ze,ze, ze,-on,ze, ze,ze,on)
    call symopadd (xmat)
    if(point_group(1:1) .eq. 'D') then
!...........sigma_h:
      call symmatfill (xmat, on,ze,ze, ze,on,ze, ze,ze,-on)
      call symopadd (xmat)
    endif
 
  endif

!
!.....Check the closure of the symmetry operators set:
!
  call symclosure (uout)
  call symgetgroup (uout)
 
600 format (/                                               &
    1x, '++SYMPURIFY: Purify the symmetry operators.'/      &
    1x, 'Produce the group from a clean set of generators', &
    1x, 'with the chemical orientation.')

end subroutine sympurify
subroutine symreplicate (uout, xa1, ya1, za1)
!
! symreplicate - apply the molecular symmetry operations on the
! (xa1,ya1,za1) point and print all the images.
!
  implicit none

  integer :: uout
  real(kind=8) :: xa1, ya1, za1

  real(kind=8) :: xx, yy, zz, xcpy(MOPSYM,3)
  integer :: i, j, ncpy
  logical :: newcpy

  write (uout,600) xa1, ya1, za1
  ncpy = 0
  do i = 1, nopsym
    xx = opsym(i,1,1)*xa1 + opsym(i,1,2)*ya1 + opsym(i,1,3)*za1
    yy = opsym(i,2,1)*xa1 + opsym(i,2,2)*ya1 + opsym(i,2,3)*za1
    zz = opsym(i,3,1)*xa1 + opsym(i,3,2)*ya1 + opsym(i,3,3)*za1
    newcpy = .true.
    j = 1
    do while (newcpy .and. j.le.ncpy)
      if (abs(xcpy(j,1)-xx) .le. TOLdist .and. &
          abs(xcpy(j,2)-yy) .le. TOLdist .and. &
          abs(xcpy(j,3)-zz) .le. TOLdist) then
        newcpy = .false.
      else
        j = j + 1
      endif
    enddo
    if (newcpy) then
      ncpy = ncpy + 1
      xcpy(ncpy,1) = xx
      xcpy(ncpy,2) = yy
      xcpy(ncpy,3) = zz
      if (opproper(i)) then
        write (uout,605)  i, xx, yy, zz
      else
        write (uout,605) -i, xx, yy, zz
      endif
    endif
  enddo
 
 600 format (/                                                          &
     1x, '++SYMREPLICATE: Different symmetrical images of the point ',  &
     3f12.8)
 605 format (1x, 'Opsym: ', i5, '  Image: ', 3f18.12) 
 
end subroutine symreplicate
subroutine symrwcoord(uout, ax, ay, az, atZmol, mass, doread, dowrite, n)
!
! symrwcoord - read in the atomic coordinates from the molecular
! database or write them back to the database from the local
! matrices.
!
  implicit none

  integer :: n, uout
  integer, dimension(n) :: atZmol
  real(kind=8), dimension(n) :: ax, ay, az, mass
  logical :: doread, dowrite

  real(kind=8) :: wi, wt, zi, zt, xcq, ycq, zcq
  real(kind=8), dimension(n) :: xmol(n), ymol(n), zmol(n)
  integer :: i, qi, atmZi

  if (doread) then
    if (uout.ge.0) write (uout,600)
    write (uout,600)
    ! get the center of mass:
    xcm = 0d0
    ycm = 0d0
    zcm = 0d0
    xcq = 0d0
    ycq = 0d0
    zcq = 0d0 
    wt = 0d0
    zt = 0d0
    do i = 1, n
      zi = atZmol(i)
      wi = mass(i)
      xmol(i) = ax(i)
      ymol(i) = ay(i)
      zmol(i) = az(i)
      wt = wt + wi
      xcm = xcm + wi * xmol(i)
      ycm = ycm + wi * ymol(i)
      zcm = zcm + wi * zmol(i)
      zt = zt + zi
      xcq = xcq + zi * xmol(i)
      ycq = ycq + zi * ymol(i)
      zcq = zcq + zi * zmol(i)
    enddo
    xcm = xcm / wt
    ycm = ycm / wt
    zcm = zcm / wt
    xcq = xcq / zt
    ycq = ycq / zt
    zcq = zcq / zt
    if (uout.gt.0) write (uout,605) xcm, ycm, zcm, xcq, ycq, zcq
    ! move the origin to the center of mass and store locally the
    do i = 1, n
      ax(i) = xmol(i) - xcm
      ay(i) = ymol(i) - ycm
      az(i) = zmol(i) - zcm
    enddo
  endif

  if (dowrite) then
    if (uout.ge.0) write (uout,650)
    do i = 1, n
      xmol(i) = ax(i) + xcm
      ymol(i) = ay(i) + ycm
      zmol(i) = az(i) + zcm
    enddo
  endif

 600 format (/                                                        &
     1x, '++SYMRWCOORD (R):'/                                         &
     1x, 'Copy coordinates from the molecular database to the local ',&
     1x, 'symmetry arrays.')
 605 format (                                                         &
     1x, 'Center of mass:            ', 3f12.6/                       &
     1x, 'Center of nuclear charges: ', 3f12.6)
 650 format (/                                                        &
     1x, '++SYMRWCOORD (W):'/                                         &
     1x, 'Copy back coordinates to the molecular database')

end subroutine symrwcoord
logical function symsingular(xmat)
!
! symsingular - check if this 3x3 matrix is singular
!
  implicit none

  real(kind=8), dimension(3,3) :: xmat
   
  real(kind=8) :: det
  real(kind=8), dimension(3,3) :: adj

  adj(1,1) = xmat(2,2)*xmat(3,3) - xmat(2,3)*xmat(3,2)
  adj(1,2) = xmat(2,3)*xmat(3,1) - xmat(2,1)*xmat(3,3)
  adj(1,3) = xmat(3,2)*xmat(2,1) - xmat(2,2)*xmat(3,1)

  det = xmat(1,1)*adj(1,1) + xmat(1,2)*adj(1,2) + xmat(1,3)*adj(1,3)

  if (abs(det) .le. TOLsng) then
    symsingular = .true.
  else
    symsingular = .false.
  endif

  return

end function symsingular
subroutine symtransform (xleft, xright)
!
! Symtransform - multiply all symmetry operators by the xleft matrix
! from the left and by the xright matrix from the right. This is
! intended to produce a similarity transformation.
!
  implicit none

  real(kind=8), dimension(3,3) :: xleft, xright

  integer, parameter :: MOP = 20
  integer :: i, j, k, l, iop, ierror, ii, ij, ik, ibest
  real(kind=8), dimension(3,3) :: xmat, xinv, xvect
  real(kind=8), dimension(3) :: xeul, xval
  real(kind=8) :: xtrace, xdet, xnorm, xnm, ang1, angle, err, besterr
  character(len=40) :: fmt
  logical :: found

  do iop = 1, nopsym
    do i = 1, 3
      do j = 1, 3
        xmat(i,j) = 0d0
        do k = 1, 3
          do l = 1, 3
            xmat(i,j) = xmat(i,j) + xleft(i,k)*opsym(iop,k,l)*xright(l,j)
          enddo
        enddo
      enddo
    enddo
    do i = 1, 3
      do j = 1, 3
        opsym(iop,i,j) = xmat(i,j)
      enddo
    enddo
  enddo

!
! Recalculate the properties of the transformed matrices:
!
  do iop = 1, nopsym
    do i = 1, 3
      do j = 1, 3
        xmat(i,j) = opsym(iop,i,j)
      enddo
    enddo
    xtrace = 0d0
    do i = 1, 3
      xtrace = xtrace + xmat(i,i)
    enddo
    call syminv (xmat, xinv, xdet, ierror)
    if(ierror .lt. 0) call error ('symtransform', &
    'Operator has a singular matrix!!!', faterr)
!
!........Get the Euler angles for the rotation:
!
  call symeuler (xmat, xeul)
  opeuler(iop,1) = xeul(1)
  opeuler(iop,2) = xeul(2)
  opeuler(iop,3) = xeul(3)
!
!........Is this a proper or improper rotation?
!
  if(abs(abs(xdet)-1d0) .gt. TOLisint) then
    call error ('symtransform', 'determinant is not +-1', warning)
  else if(xdet .gt. 0d0) then
    opproper(iop) = .true.
  else
    opproper(iop) = .false.
  endif
!
!........Get the rotation axis except for E and i, for which it is
!        fully degenerated:
!
  if(abs(abs(xtrace)-3d0) .le. TOLisint) then
    if(xtrace .gt. 0d0) then
      opsymbol(iop) = 'E'
      opaxis(iop,1) = 0d0
      opaxis(iop,2) = 0d0
      opaxis(iop,3) = 1d0
      oporder(iop) = 1
      opm(iop) = 1
      opangle(iop) = dpi + dpi
      optype(iop) = opt_identity
    else
      opsymbol(iop) = 'i'
      opaxis(iop,1) = 0d0
      opaxis(iop,2) = 0d0
      opaxis(iop,3) = 1d0
      oporder(iop) = 2
      opm(iop) = 1
      opangle(iop) = dpi
      optype(iop) = opt_inversion
    endif
  else
!
!...........Get the rotation angle and the n and m of the C_n^m or S_n^m
!           symmetry operator:
!
        ang1 = 0.5d0*(xtrace-xdet)
        if(abs(ang1) .gt. 1d0) ang1 = sign(1d0, ang1)
        angle = acos(ang1)
        if(abs(angle) .le. TOLnull) then
          xnm = 1d0
          angle = dpi + dpi
        else
          xnm = 2d0 * dpi / angle
        endif
        opangle(iop) = angle
        ii = 0
        found = .false.
        do while(.not. found .and. ii.le.MOP)
          ii = ii + 1
          found = abs(xnm*ii - nint(xnm*ii)) .le. TOLisint
        enddo
        if(found) then
          oporder(iop) = nint(xnm*ii)
          opm(iop) = ii
        else
!          call error ('symtransform', 'Using best approx order', warning)
          ibest = 1
          besterr = abs(xnm - nint(xnm))
          do ii = 2, MOP
            err = abs(xnm*ii - nint(xnm*ii))
            if(err .le. besterr) then
              besterr = err
              ibest = ii
            endif
          enddo
          oporder(iop) = nint(xnm*ibest)
          opm(iop) = ibest
        endif
        write (fmt,200) 1+int(log10(1d0*oporder(iop))),1+int(log10(1d0*opm(iop)))

        call symeigen (xmat, xval, xvect, xdet, ierror)
        if(ierror .ne. 0) then
          call error ('symtransform', 'Trouble finding the rotation axis', warning)
        endif
        opaxis(iop,1) = xvect(1,3)
        opaxis(iop,2) = xvect(2,3)
        opaxis(iop,3) = xvect(3,3)
        if(xdet .gt. 0d0) then
          write (opsymbol(iop),fmt) 'C', oporder(iop), opm(iop)
          optype(iop) = opt_rotation
        else if(abs(xtrace - 1d0) .le. TOLisint) then
          write (opsymbol(iop),210) 'sigma'
          optype(iop) = opt_sigma
        else if (oporder(iop).eq.1 .and. opm(iop).eq.1) then
          write (opsymbol(iop),210) 'sigma'
          optype(iop) = opt_sigma
        else
          write (opsymbol(iop),fmt) 'S', oporder(iop), opm(iop)
          optype(iop) = opt_imp_rotation
        endif
      endif

    enddo
 
200 format ('(a, "_", i', i2, ', "^", i', i2, ')')
210 format (a)
 
  end subroutine symtransform

      subroutine symgetgroup (uout)
!
!.....symgetgroup - determine the point group from the set of symmetry
!     operators.
!
      implicit none
!
      integer           uout
!
      integer           MOP
      parameter         (MOP = 20)
      integer           mainorder, mainmult, imain(MOP), mainimproper
      integer           binmult, ibin(MOP), nsigma, isigma(MOP), isigmah
      integer           nbinperp, i, j, j1
      logical           inversion, sigmah
      real*8            xprod
      character*(4)     chmain, chimp
!
!      integer           leng
!
!.....Get the main axis, if any, and their multiplicity:
!     Only the C_n^1 operations will be considered.
!     Unfortunately, C_n^1 and C_n^(n-1) operations cannot be
!     distinguished.
!     Get also the list of binary axes, the sigma planes, and the
!     order of the main improper axis (n > 2).
!
      mainorder = 0
      mainmult = 0
      mainimproper = 0
      binmult = 0
      inversion = .false.
      nsigma = 0
      do i = 1, nopsym
         if (optype(i) .eq. opt_rotation .and. opm(i) .eq. 1) then
            if (oporder(i) .gt. mainorder) then
               mainorder = oporder(i)
               mainmult = 1
               imain(mainmult) = i
            else if (oporder(i) .eq. mainorder) then
               mainmult = mainmult + 1
               if (mainmult .gt. MOP) call error ('symgetgroup',   &
                 'Main axis mult too high. Increase MOP!', faterr)
               imain(mainmult) = i
            endif
            if (oporder(i) .eq. 2) then
               binmult = binmult + 1
               if (binmult .gt. MOP) call error ('symgetgroup', &
                 'Binary axis mult too high. Increase MOP!', faterr)
               ibin(binmult) = i
            endif
         else if (optype(i) .eq. opt_inversion) then
            inversion = .true.
         else if (optype(i) .eq. opt_sigma) then
            nsigma = nsigma + 1
            if (nsigma .gt. MOP) call error ('symgetgroup', &
              'Too many sigma planes. Increase MOP!', faterr)
            isigma(nsigma) = i
         else if (optype(i) .eq. opt_imp_rotation) then
            if (oporder(i) .gt. mainimproper) mainimproper = oporder(i)
         endif
      enddo
      if (mainorder .gt. 2) mainmult = mainmult/2
!
!.....If there is a single main order group, look for sigma_h planes:
!
      sigmah = .false.
      if (mainmult .eq. 1 .or. mainorder .eq. 2) then
         i = imain(1)
         j = 1
         do while (.not.sigmah .and. j .le. nsigma)
            j1 = isigma(j)
            xprod = opaxis(i,1)*opaxis(j1,1)  &
                  + opaxis(i,2)*opaxis(j1,2)  &
                  + opaxis(i,3)*opaxis(j1,3)
            if (abs(abs(xprod)-1d0) .le. TOLisint) then
               sigmah = .true.
               isigmah = j1
            endif
            j = j + 1
         enddo
      endif
      if (sigmah) then
         opsymbol(isigmah) = 'sigma_h'
      endif
!
!.....If there is a single main order group, look for n C2 perpendicular
!     axis:
!
      nbinperp = 0
      if (mainmult .eq. 1 .or. mainorder .eq. 2) then
         i = imain(1)
         do j = 1, binmult
            j1 = ibin(j)
            xprod = opaxis(i,1)*opaxis(j1,1) &
                  + opaxis(i,2)*opaxis(j1,2) &
                  + opaxis(i,3)*opaxis(j1,3)
            if (abs(xprod) .le. TOLnull) then
               nbinperp = nbinperp + 1
            endif
         enddo
      endif
!
!.....Store as a character the main proper and improper orders:
!
      if (mainorder .gt. 999) then
         write (chmain, '(i4)') mainorder
      else if (mainorder .gt. 99) then
         write (chmain, '(i3)') mainorder
      else if (mainorder .gt. 9) then
         write (chmain, '(i2)') mainorder
      else
         write (chmain, '(i1)') mainorder
      endif
      if (mainimproper .gt. 999) then
         write (chimp, '(i4)') mainimproper
      else if (mainimproper .gt. 99) then
         write (chimp, '(i3)') mainimproper
      else if (mainimproper .gt. 9) then
         write (chimp, '(i2)') mainimproper
      else
         write (chimp, '(i1)') mainimproper
      endif
      opmainord = mainorder
!
!.....Decision tree to get the point group. Start with the cubic-like
!     groups:
!
      if (mainmult .gt. 1 .and. mainorder .gt. 2) then
         if (mainorder .eq. 5 .and. mainmult .ge. 6) then
            if (inversion) then
               point_group = 'Ih'
            else
               point_group = 'I'
            endif
            point_g2 = point_group
            opmainord = 5
         else if (mainorder .eq. 4 .and. mainmult .ge. 3) then
            if (inversion) then
               point_group = 'Oh'
            else
               point_group = 'O'
            endif
            point_g2 = point_group
            opmainord = 4
         else if (mainorder .eq. 3 .and. mainmult .ge. 4) then
            if (inversion) then
               point_group = 'Th'
            else if (nsigma .eq. 6) then
               point_group = 'Td'
            else
               point_group = 'T'
            endif
            point_g2 = point_group
            opmainord = 3
         else
            call error ('symgetgroup', 'Unknown cubic-like group', &
                        warning)
            write (uout,700) mainorder, mainmult, inversion, binmult &
                           , nbinperp, nsigma, sigmah, mainimproper
            point_group = '??'
         endif
      else if (mainorder .ge. 2) then
         if (mainorder .eq. nbinperp) then
            if (sigmah) then
!AQUI
!              write (point_group,500) 'Dnh', mainorder
               point_group = "D" // chmain(1:leng(chmain)) // "h"
               !point_group = "D" // chmain // "h"
               point_g2 = 'Dnh'
            else if (mainorder .eq. nsigma) then
!              write (point_group,500) 'Dnd', mainorder
               point_group = "D" // chmain(1:leng(chmain)) // "d"
               !point_group = "D" // chmain // "d"
               point_g2 = 'Dnd'
            else
!              write (point_group,500) 'Dn', mainorder
               point_group = "D" // chmain(1:leng(chmain))
               !point_group = "D" // chmain
               point_g2 = 'Dn'
            endif
         else
            if (sigmah) then
!              write (point_group,500) 'Cnh', mainorder
               point_group = "C" // chmain(1:leng(chmain)) // "h"
               !point_group = "C" // chmain // "h"
               point_g2 = 'Cnh'
            else if (mainorder .eq. nsigma) then
!              write (point_group,500) 'Cnv', mainorder
               point_group = "C" // chmain(1:leng(chmain)) // "v"
               !point_group = "C" // chmain // "v"
               point_g2 = 'Cnv'
            else if (2*mainorder .eq. mainimproper) then
!              write (point_group,500) 'S2n', mainorder
               point_group = "S" // chimp(1:leng(chimp))
               !point_group = "S" // chmain 
               point_g2 = 'S2n'
               opmainord = mainorder
               opmainord = mainimproper
            else
!              write (point_group,500) 'Cn', mainorder
               point_group = "C" // chmain(1:leng(chmain))
               !point_group = "C" // chmain 
               point_g2 = 'Cn'
            endif
         endif
      else if (nopsym .eq. 2) then
         if (nsigma .eq. 1) then
            point_group = 'Cs'
            point_g2 = 'Cs'
         else if (inversion) then
            point_group = 'Ci'
            point_g2 = 'Ci'
         else
            call error ('symgetgroup', 'Unknown group', warning)
            write (uout,700) mainorder, mainmult, inversion, binmult &
                           , nbinperp, nsigma, sigmah, mainimproper
            point_group = '??'
         endif
      else if (nopsym .eq. 1) then
         point_group = 'C1'
         point_g2 = 'C1'
      else
         call error ('symgetgroup', 'Unknown group', warning)
         write (uout,700) mainorder, mainmult, inversion, binmult &
                        , nbinperp, nsigma, sigmah, mainimproper
         point_group = '??'
      endif
!
 500  format (a, ', n = ', i2)
!
 700  format (/                                                &
       1x, '++SYMGETGROUP:'/                                   &
       1x, 'DBG(symgetgroup) main order & mult.......:', 2i5/  &
       1x, 'DBG(symgetgroup) inversion...............:', l5/   &
       1x, 'DBG(symgetgroup) binary axes.............:', i5/   &
       1x, 'DBG(symgetgroup) binary axes perp to main:', i5/   &
       1x, 'DBG(symgetgroup) symmetry planes.........:', i5/   &
       1x, 'DBG(symgetgroup) sigma_h plane...........:', l5/   &
       1x, 'DBG(symgetgroup) main improper axis......:', i5)
!
      end subroutine

      subroutine symorient (uout, ax, ay, az, atZmol, mass, n)
!
!.....symorient - reorient the molecule according to the symmetry
!     elements.
!
      implicit none
!
      integer n
      integer           uout, atZmol(n)
      real*8            ax(n), ay(n), az(n), mass(n)
!
      integer           mainpro, imainpro, mainimp, imainimp
      integer           ifirst, isecond, isigmad
      integer           i, j, imin, ierror, nrot
      logical           found
      real*8            xprod, v(3,3), vinv(3,3), vdet
      real*8            imatrix(3,3), ilambda(3), irot(3,3)
      real*8            xmin, xortho, xnorm, coef, xang
      real*8            xx, yy, zz, wi, xxx, yyy, zzz, wt
!
      write (uout,600)
!
!.....Build the inertia matrix and get the rotation axes:
!
      do i = 1, 3
         do j = 1, 3
            imatrix(i,j) = 0d0
         enddo
      enddo
      do i = 1, n
         !wi = atmweight(atZmol(i))
         wi = mass(i)
         xx = ax(i)
         yy = ay(i)
         zz = az(i)
         imatrix(1,1) = imatrix(1,1) + wi * (yy*yy + zz*zz)
         imatrix(2,2) = imatrix(2,2) + wi * (zz*zz + xx*xx)
         imatrix(3,3) = imatrix(3,3) + wi * (xx*xx + yy*yy)
         imatrix(1,2) = imatrix(1,2) - wi * xx * yy
         imatrix(2,3) = imatrix(2,3) - wi * yy * zz
         imatrix(3,1) = imatrix(3,1) - wi * zz * xx
      enddo
      imatrix(2,1) = imatrix(1,2)
      imatrix(3,2) = imatrix(2,3)
      imatrix(1,3) = imatrix(3,1)
!
      write (uout,605) ((imatrix(i,j), j = 1, 3), i = 1, 3)
!
      call jacobi (imatrix, 3, 3, ilambda, irot, nrot)
      if (nrot .lt. 0) call error ('symorient',  &
        'Error on the diagonalization of the inertia matrix', faterr)
!
      write (uout,610) (ilambda(i), i = 1, 3), &
                       ((irot(i,j), j = 1, 3), i = 1, 3)
!
!.....Special orientation for the different symmetry groups:
!
      do i = 1, 3
         do j = 1, 3
            v(i,j) = 0d0
         enddo
      enddo
!
!.......................................................................
      if (point_group(1:1) .eq. 'O') then
!.......................................................................
!........Octahedral groups. The first two C_4^1 axes determine z and x
!
         ifirst = 1
         found = .false.
         do while (.not.found .and. ifirst.le.nopsym)
            if (optype(ifirst).eq.opt_rotation &
               .and. oporder(ifirst).eq.4      &
               .and. opm(ifirst).eq.1) then
               found = .true.
            else
               ifirst = ifirst + 1
            endif
         enddo
         if (.not.found) call error ('symorient', &
                              '(O,1) orientation not found', faterr)
         v(1,3) = opaxis(ifirst,1)
         v(2,3) = opaxis(ifirst,2)
         v(3,3) = opaxis(ifirst,3)
         isecond = 1
         found = .false.
         do while (.not.found .and. isecond.le.nopsym)
            if (isecond.ne.ifirst                     &
               .and. optype(isecond).eq.opt_rotation  &
               .and. oporder(isecond).eq.4            &
               .and. opm(isecond).eq.1) then
               xprod = abs(opaxis(ifirst,1)*opaxis(isecond,1) &
                         + opaxis(ifirst,2)*opaxis(isecond,2) &
                         + opaxis(ifirst,3)*opaxis(isecond,3))
               if (xprod.lt.0.05d0) then
                  found = .true.
               else
                  isecond = isecond + 1
               endif
            else
               isecond = isecond + 1
            endif
         enddo
         if (.not.found) call error ('symorient', &
                              '(O,2) orientation not found', faterr)
         v(1,1) = opaxis(isecond,1)
         v(2,1) = opaxis(isecond,2)
         v(3,1) = opaxis(isecond,3)
!
!.......................................................................
      else if (point_group(1:1) .eq. 'T') then
!.......................................................................
!........Tetrahedral groups. The first two C_2^1 axes determine z and x
!
         ifirst = 1
         found = .false.
         do while (.not.found .and. ifirst.le.nopsym)
            if (optype(ifirst).eq.opt_rotation &
               .and. oporder(ifirst).eq.2      &
               .and. opm(ifirst).eq.1) then
               found = .true.
            else
               ifirst = ifirst + 1
            endif
         enddo
         if (.not.found) call error ('symorient', &
                              '(T,1) orientation not found', faterr)
         v(1,3) = opaxis(ifirst,1)
         v(2,3) = opaxis(ifirst,2)
         v(3,3) = opaxis(ifirst,3)
         isecond = 1
         found = .false.
         do while (.not.found .and. isecond.le.nopsym)
            if (optype(isecond).eq.opt_rotation  &
               .and. oporder(isecond).eq.2       &
               .and. opm(isecond).eq.1) then
               xprod = abs(opaxis(ifirst,1)*opaxis(isecond,1) &
                         + opaxis(ifirst,2)*opaxis(isecond,2) &
                         + opaxis(ifirst,3)*opaxis(isecond,3))
               if (xprod.le.1d-5) then
                  found = .true.
               else
                  isecond = isecond + 1
               endif
            else
               isecond = isecond + 1
            endif
         enddo
         if (.not.found) call error ('symorient', &
                              '(T,2) orientation not found', faterr)
         v(1,1) = opaxis(isecond,1)
         v(2,1) = opaxis(isecond,2)
         v(3,1) = opaxis(isecond,3)
!
!.......................................................................
      else if (point_group(1:1) .eq. 'I') then
!.......................................................................
!........Icosahedral groups.
!
!........Dirty trick to select one of the possible orientations:
!
         goto 504
!
!      +------- Orientation 1 --------------+
!      | The first C_5^1 axis determines z. |
!      | A C_3^1 axis determines y.         |
!      +------------------------------------+
!
         continue
         ifirst = 1
         found = .false.
         do while (.not.found .and. ifirst.le.nopsym)
            if (optype(ifirst).eq.opt_rotation   &
               .and. oporder(ifirst).eq.5        &
               .and. opm(ifirst).eq.1) then
               found = .true.
            else
               ifirst = ifirst + 1
            endif
         enddo
         if (.not.found) call error ('symorient', &
                              '(I,1) orientation not found', faterr)
         v(1,3) = opaxis(ifirst,1)
         v(2,3) = opaxis(ifirst,2)
         v(3,3) = opaxis(ifirst,3)
         isecond = 1
         found = .false.
         do while (.not.found .and. isecond.le.nopsym)
            if (isecond.ne.ifirst                      &
               .and. optype(isecond).eq.opt_rotation   &
               .and. oporder(isecond).eq.3             &
               .and. opm(isecond).eq.1) then
               xprod = abs(opaxis(ifirst,1)*opaxis(isecond,1)  &
                         + opaxis(ifirst,2)*opaxis(isecond,2)  &
                         + opaxis(ifirst,3)*opaxis(isecond,3))
               xang = abs(acos(xprod)*180d0/dpi)
               if (abs(xang-37.377d0).lt.1d-1 .or.  &
                   abs(xang-142.623d0).lt.1d-1) then
                  found = .true.
               else
                  isecond = isecond + 1
               endif
            else
               isecond = isecond + 1
            endif
         enddo
         if (.not.found) call error ('symorient', &
                              '(I,2) orientation not found', faterr)
         v(1,2) = opaxis(isecond,1)
         v(2,2) = opaxis(isecond,2)
         v(3,2) = opaxis(isecond,3)
!
!........Get the x direction as a vectorial product of y and z:
!
         v(1,1) = v(2,2)*v(3,3) - v(3,2)*v(2,3)
         v(2,1) = v(3,2)*v(1,3) - v(1,2)*v(3,3)
         v(3,1) = v(1,2)*v(2,3) - v(2,2)*v(1,3)
         xnorm = v(1,1)*v(1,1) + v(2,1)*v(2,1) + v(3,1)*v(3,1)
         if (xnorm .le. 1d-5) call error ('symorient',  &
                   'Error normalizing orientation axis (I,x)', faterr)
         xnorm = 1d0 / sqrt(xnorm)
         v(1,1) = xnorm * v(1,1)
         v(2,1) = xnorm * v(2,1)
         v(3,1) = xnorm * v(3,1)
!
         goto 599
!
!      +------- Orientation 2 --------------+
!      | The first C_5^1 axis determines z. |
!      | A C_3^1 axis determines x.         |
!      +------------------------------------+
!
         continue
         ifirst = 1
         found = .false.
         do while (.not.found .and. ifirst.le.nopsym)
            if (optype(ifirst).eq.opt_rotation   &
               .and. oporder(ifirst).eq.5        &
               .and. opm(ifirst).eq.1) then
               found = .true.
            else
               ifirst = ifirst + 1
            endif
         enddo
         if (.not.found) call error ('symorient', &
                              '(I,1) orientation not found', faterr)
         v(1,3) = opaxis(ifirst,1)
         v(2,3) = opaxis(ifirst,2)
         v(3,3) = opaxis(ifirst,3)
         isecond = 1
         found = .false.
         do while (.not.found .and. isecond.le.nopsym)
            if (isecond.ne.ifirst                     &
               .and. optype(isecond).eq.opt_rotation  &
               .and. oporder(isecond).eq.3            &
               .and. opm(isecond).eq.1) then
               xprod = abs(opaxis(ifirst,1)*opaxis(isecond,1) &
                         + opaxis(ifirst,2)*opaxis(isecond,2) &
                         + opaxis(ifirst,3)*opaxis(isecond,3))
               xang = abs(acos(xprod)*180d0/dpi)
               if (abs(xang-37.377d0).lt.1d-1 .or.    &
                   abs(xang-142.623d0).lt.1d-1) then
                  found = .true.
               else
                  isecond = isecond + 1
               endif
            else
               isecond = isecond + 1
            endif
         enddo
         if (.not.found) call error ('symorient',   &
                              '(I,2) orientation not found', faterr)
         v(1,1) = opaxis(isecond,1)
         v(2,1) = opaxis(isecond,2)
         v(3,1) = opaxis(isecond,3)
         goto 599
!
!      +------- Orientation 3 --------------------+
!      | The first C_5^1 axis determines z.       |
!      | A perpendicular C_2^1 axis determines x. |
!      +------------------------------------------+
!
         continue
         ifirst = 1
         found = .false.
         do while (.not.found .and. ifirst.le.nopsym)
            if (optype(ifirst).eq.opt_rotation  &
               .and. oporder(ifirst).eq.5       &
               .and. opm(ifirst).eq.1) then
               found = .true.
            else
               ifirst = ifirst + 1
            endif
         enddo
         if (.not.found) call error ('symorient', &
                              '(I,1) orientation not found', faterr)
         v(1,3) = opaxis(ifirst,1)
         v(2,3) = opaxis(ifirst,2)
         v(3,3) = opaxis(ifirst,3)
         isecond = 1
         found = .false.
         do while (.not.found .and. isecond.le.nopsym)
            if (isecond.ne.ifirst                     &
               .and. optype(isecond).eq.opt_rotation  &
               .and. oporder(isecond).eq.2            &
               .and. opm(isecond).eq.1) then
               xprod = abs(opaxis(ifirst,1)*opaxis(isecond,1) &
                         + opaxis(ifirst,2)*opaxis(isecond,2) &
                         + opaxis(ifirst,3)*opaxis(isecond,3))
               if (xprod.lt.1d-4) then
                  found = .true.
               else
                  isecond = isecond + 1
               endif
            else
               isecond = isecond + 1
            endif
         enddo
         if (.not.found) call error ('symorient',  &
                              '(I,2) orientation not found', faterr)
         v(1,1) = opaxis(isecond,1)
         v(2,1) = opaxis(isecond,2)
         v(3,1) = opaxis(isecond,3)
         goto 599
!
!      +------- Orientation 4 ---------------------------+
!      | Two perpendicular C_2^1 axes determine z and x. |
!      | Actually, there are three mutually perp C_2^1   |
!      | axes forming {x,y,z}.                           |
!      +-------------------------------------------------+
!
 504     continue
         ifirst = 1
         found = .false.
         do while (.not.found .and. ifirst.le.nopsym)
            if (optype(ifirst).eq.opt_rotation  &
               .and. oporder(ifirst).eq.2       &
               .and. opm(ifirst).eq.1) then
               found = .true.
            else
               ifirst = ifirst + 1
            endif
         enddo
         if (.not.found) call error ('symorient', &
                              '(I,1) orientation not found', faterr)
         v(1,3) = opaxis(ifirst,1)
         v(2,3) = opaxis(ifirst,2)
         v(3,3) = opaxis(ifirst,3)
         isecond = 1
         found = .false.
         do while (.not.found .and. isecond.le.nopsym)
            if (isecond.ne.ifirst                      &
               .and. optype(isecond).eq.opt_rotation   &
               .and. oporder(isecond).eq.2             &
               .and. opm(isecond).eq.1) then
               xprod = abs(opaxis(ifirst,1)*opaxis(isecond,1) &
                         + opaxis(ifirst,2)*opaxis(isecond,2) &
                         + opaxis(ifirst,3)*opaxis(isecond,3))
               if (xprod.lt.1d-4) then
                  found = .true.
               else
                  isecond = isecond + 1
               endif
            else
               isecond = isecond + 1
            endif
         enddo
         if (.not.found) call error ('symorient',     &
                              '(I,2) orientation not found', faterr)
         v(1,1) = opaxis(isecond,1)
         v(2,1) = opaxis(isecond,2)
         v(3,1) = opaxis(isecond,3)
         goto 599
!
 599     continue
!
!.......................................................................
      else if (point_group(1:2) .eq. 'D2') then
!.......................................................................
!........Dihedral D2, D2h, and D2d groups.
!        The main axis determines z.
!        The next non-collinear one determines y.
!
         ifirst = -1
         isecond = -1
         if (point_group(1:3) .eq. 'D2h') then
!
!...........The normal to sigmah defines z. Any of the two C2
!           perpendicular to z defines x.
!
            i = 1
            do while (ifirst.le.0 .and. i.le.nopsym)
               if (optype(i).eq.opt_sigma) ifirst = i
               i = i + 1
            enddo
            if (ifirst.le.0) call error ('symorient', &
                                '(D2h,1) orientation not found', faterr)
            i = 1
            do while (isecond.le.0 .and. i.le.nopsym)
               if (optype(i).eq.opt_rotation) then
                  xprod = opaxis(ifirst,1)*opaxis(i,1) &
                        + opaxis(ifirst,2)*opaxis(i,2) &
                        + opaxis(ifirst,3)*opaxis(i,3)
                  if (abs(xprod) .le. 1d-2) isecond = i
               endif
               i = i + 1
            enddo
            if (isecond.le.0) call error ('symorient', &
                                '(D2h,2) orientation not found', faterr)
!
         else if (point_group(1:3) .eq. 'D2d') then
!
!...........We get any of the two sigmad first. The z direction is given
!           by the only C2 axis perpendicular to the normal of sigmad.
!           Any of the two C2 perpendicular to z defines x.
!
            isigmad = -1
            i = 1
            do while (isigmad.le.0 .and. i.le.nopsym)
               if (optype(i).eq.opt_sigma) isigmad = i
               i = i + 1
            enddo
            if (isigmad.le.0) call error ('symorient',  &
                                '(D2d,s) orientation not found', faterr)
            i = 1
            do while (ifirst.le.0 .and. i.le.nopsym)
               if (optype(i).eq.opt_rotation) then
                  xprod = opaxis(isigmad,1)*opaxis(i,1) &
                        + opaxis(isigmad,2)*opaxis(i,2) &
                        + opaxis(isigmad,3)*opaxis(i,3)
                  if (abs(xprod) .le. 1d-2) ifirst = i
               endif
               i = i + 1
            enddo
            if (ifirst.le.0) call error ('symorient',   &
                                '(D2d,1) orientation not found', faterr)
            i = 1
            do while (isecond.le.0 .and. i.le.nopsym)
               if (optype(i).eq.opt_rotation            &
                  .and. i.ne.ifirst) isecond = i
               i = i + 1
            enddo
            if (isecond.le.0) call error ('symorient',  &
                                '(D2d,2) orientation not found', faterr)
!
         else
!
!...........Any two of the three C2 axes define z and x.
!
            i = 1
            do while (ifirst.le.0 .and. i.le.nopsym)
               if (optype(i).eq.opt_rotation) ifirst = i
               i = i + 1
            enddo
            if (ifirst.le.0) call error ('symorient',  &
                                '(D2,1) orientation not found', faterr)
            i = 1
            do while (isecond.le.0 .and. i.le.nopsym)
               if (optype(i).eq.opt_rotation           &
                  .and. i.ne.ifirst) isecond = i
               i = i + 1
            enddo
            if (isecond.le.0) call error ('symorient', &
                                '(D2,2) orientation not found', faterr)
!     
         endif
!     
         v(1,3) = opaxis(ifirst,1)
         v(2,3) = opaxis(ifirst,2)
         v(3,3) = opaxis(ifirst,3)
         v(1,1) = opaxis(isecond,1)
         v(2,1) = opaxis(isecond,2)
         v(3,1) = opaxis(isecond,3)
!
!.......................................................................
      else if (point_group(1:1) .eq. 'D') then
!.......................................................................
!........Dihedral Dn, Dnh, and Dnd groups with n>2.
!        The main axis fixes z. Any of the two perpendicular C2 gives x.
!        The main axis can be a S2n improper axis in Dnd groups.
!
         ifirst = -1
         i = 1
         do while (ifirst.le.0 .and. i.le.nopsym)
            if ( (optype(i).eq.opt_rotation         &
                .or. optype(i).eq.opt_imp_rotation) &
                .and. (oporder(i).eq.opmainord)     &
                .and. (opm(i).eq.1) ) ifirst = i    
            i = i + 1
         enddo
         if (ifirst.le.0) call error ('symorient',  &
                              '(D,1) orientation not found', faterr)
         isecond = -1
         i = 1
         do while (isecond.le.0 .and. i.le.nopsym)
            if ( (optype(i).eq.opt_rotation)        &
                .and. (oporder(i).eq.2)             &
                .and. (opm(i).eq.1) ) then
               xprod = abs(opaxis(ifirst,1)*opaxis(i,1) &
                         + opaxis(ifirst,2)*opaxis(i,2) &
                         + opaxis(ifirst,3)*opaxis(i,3))
               if (abs(xprod) .le. 1d-2) isecond = i
            endif
            i = i + 1
         enddo
         if (isecond.le.0) call error ('symorient', &
                              '(D,2) orientation not found', faterr)
         v(1,3) = opaxis(ifirst,1)
         v(2,3) = opaxis(ifirst,2)
         v(3,3) = opaxis(ifirst,3)
         v(1,1) = opaxis(isecond,1)
         v(2,1) = opaxis(isecond,2)
         v(3,1) = opaxis(isecond,3)
!
!.......................................................................
      else if (point_g2(1:3) .eq. 'Cnv') then
!.......................................................................
!........Use the single main axis as z, and the normal axis to a
!        sigma_v plane as x:
!
         ifirst = 1
         found = .false.
         do while (.not.found .and. ifirst.le.nopsym)
            if (optype(ifirst).eq.opt_rotation     &
               .and. oporder(ifirst).eq.opmainord  &
               .and. opm(ifirst).eq.1) then
               found = .true.
            else
               ifirst = ifirst + 1
            endif
         enddo
         if (.not.found) call error ('symorient',  &
                              '(Cnv,1) orientation not found', faterr)
         v(1,3) = opaxis(ifirst,1)
         v(2,3) = opaxis(ifirst,2)
         v(3,3) = opaxis(ifirst,3)
!
         isecond = 1
         found = .false.
         do while (.not.found .and. isecond.le.nopsym)
            if (optype(isecond).eq.opt_sigma) then
               found = .true.
            else
               isecond = isecond + 1
            endif
         enddo
         if (.not.found) call error ('symorient', &
                              '(Cnv,2) orientation not found', faterr)
         v(1,1) = opaxis(isecond,1)
         v(2,1) = opaxis(isecond,2)
         v(3,1) = opaxis(isecond,3)
!
!.......................................................................
      else if (point_g2(1:2) .eq. 'Cn'            &
         .or. point_g2(1:3) .eq. 'S2n') then
!.......................................................................
!........Use the single main axis as z.
!        If n=2 the perpendicular inertia axis of lowest eigenvalue will
!        we used as x.
!        If n>2 the inertia eigenvalues are degenerated. The first atom
!        not lying in the symmetry axis will be used to determine x.
!
         ifirst = 1
         found = .false.
         do while (.not.found .and. ifirst.le.nopsym)
            if ((optype(ifirst).eq.opt_rotation          &
               .or. optype(ifirst).eq.opt_imp_rotation)  &
               .and. oporder(ifirst).eq.opmainord        &
               .and. opm(ifirst).eq.1) then
               found = .true.
            else
               ifirst = ifirst + 1
            endif
         enddo
         if (.not.found) call error ('symorient',       &
                              '(Cn,1) orientation not found', faterr)
         v(1,3) = opaxis(ifirst,1)
         v(2,3) = opaxis(ifirst,2)
         v(3,3) = opaxis(ifirst,3)
!
!        if (opmainord .eq. 2) then
            isecond = 1
            found = .false.
            do while (.not.found .and. isecond.le.3)
               xortho = 0d0
               do j = 1, 3
                  xortho = xortho + v(j,3) * irot(j,isecond)
               enddo
               if (abs(xortho) .lt. 1e-5) then
                  found = .true.
               else
                  isecond = isecond + 1
               endif
            enddo
            if (.not.found) call error ('symorient',    &
                                 '(Cn,2) orientation not found', faterr)
!        else
!---xxx---xxx---xxx
!        endif
         v(1,1) = irot(1,isecond)
         v(2,1) = irot(2,isecond)
         v(3,1) = irot(3,isecond)
!
!.......................................................................
      else if (point_group(1:2) .eq. 'Cs') then
!.......................................................................
!........The normal to sigma_h defines z.
!        The x axis is obtained from the inertia matrix: it is
!        perpendicular to z and it has the lowest rotation eigenvalue
!        otherwise:
!
         ifirst = 1
         found = .false.
         do while (.not.found .and. ifirst.le.nopsym)
            if (optype(ifirst).eq.opt_sigma) then
               found = .true.
            else
               ifirst = ifirst + 1
            endif
         enddo
         if (.not.found) call error ('symorient',   &
                              '(Cs,1) orientation not found', faterr)
         v(1,3) = opaxis(ifirst,1)
         v(2,3) = opaxis(ifirst,2)
         v(3,3) = opaxis(ifirst,3)
!
         isecond = 1
         found = .false.
         do while (.not.found .and. isecond.le.3)
            xortho = 0d0
            do j = 1, 3
               xortho = xortho + v(j,3) * irot(j,isecond)
            enddo
            if (abs(xortho) .le. 1e-1) then
               found = .true.
            else
               isecond = isecond + 1
            endif
         enddo
         if (.not.found) call error ('symorient',  &
                              '(Cs,2) orientation not found', faterr)
         v(1,1) = irot(1,isecond)
         v(2,1) = irot(2,isecond)
         v(3,1) = irot(3,isecond)
!
!.......................................................................
      else if (point_group(1:2) .eq. 'C1'         &
          .or. point_group(1:2) .eq. 'Ci') then
!.......................................................................
!........Use the inertia eigenvectors to orient the molecule:
!
         v(1,3) = irot(1,3)
         v(2,3) = irot(2,3)
         v(3,3) = irot(3,3)
         v(1,1) = irot(1,1)
         v(2,1) = irot(2,1)
         v(3,1) = irot(3,1)
!
!.......................................................................
      else
!.......................................................................
!........C_inf_v or D_inf_h groups.
!........Use the inertia eigenvectors to orient the molecule:
!
         v(1,3) = irot(1,3)
         v(2,3) = irot(2,3)
         v(3,3) = irot(3,3)
         v(1,1) = irot(1,1)
         v(2,1) = irot(2,1)
         v(3,1) = irot(3,1)
!
!.......................................................................
      endif
!
!.....Make sure that x is orthogonal to z and both are normalized:
!
      xnorm = v(1,3)*v(1,3) + v(2,3)*v(2,3) + v(3,3)*v(3,3)
      if (xnorm .le. 1d-5) call error ('symorient',  &
                'Error normalizing orientation axis (z)', faterr)
      xnorm = 1d0 / sqrt(xnorm)
      v(1,3) = xnorm * v(1,3)
      v(2,3) = xnorm * v(2,3)
      v(3,3) = xnorm * v(3,3)
!
      coef = v(1,1)*v(1,3) + v(2,1)*v(2,3) + v(3,1)*v(3,3)
      v(1,1) = v(1,1) - coef * v(1,3)
      v(2,1) = v(2,1) - coef * v(2,3)
      v(3,1) = v(3,1) - coef * v(3,3)
      xnorm = v(1,1)*v(1,1) + v(2,1)*v(2,1) + v(3,1)*v(3,1)
      if (xnorm .le. 1d-5) call error ('symorient', &
                'Error normalizing orientation axis (x)', faterr)
      xnorm = 1d0 / sqrt(xnorm)
      v(1,1) = xnorm * v(1,1)
      v(2,1) = xnorm * v(2,1)
      v(3,1) = xnorm * v(3,1)
!
!.....Get the y direction as a vectorial product of z and x:
!
      v(1,2) = v(2,3)*v(3,1) - v(3,3)*v(2,1)
      v(2,2) = v(3,3)*v(1,1) - v(1,3)*v(3,1)
      v(3,2) = v(1,3)*v(2,1) - v(2,3)*v(1,1)
!
!.....Get the inverse matrix of v(,):
!.....The V and V^{-1} matrices produce the change of coordinates and
!     transform the symmetry matrices from the old to the new
!     coordinates:
!
      call syminv (v, vinv, vdet, ierror)
      write (uout,625) ((v(i,j), j = 1, 3), i = 1, 3), &
                       ((vinv(i,j), j = 1, 3), i = 1, 3)
      if (ierror.ne.0) call error ('symorient',  &
          'Error in the inversion of the V matrix!', faterr)
      do i = 1, n
         xx = ax(i)
         yy = ay(i)
         zz = az(i)
         xxx = vinv(1,1)*xx + vinv(1,2)*yy + vinv(1,3)*zz
         yyy = vinv(2,1)*xx + vinv(2,2)*yy + vinv(2,3)*zz
         zzz = vinv(3,1)*xx + vinv(3,2)*yy + vinv(3,3)*zz
         ax(i) = xxx
         ay(i) = yyy
         az(i) = zzz
      enddo
      write (uout,615) (i, ax(i), ay(i), az(i), i = 1, n)
      write (uout,640)
!
!.....Transform the symmetry matrices too:
!
      call symtransform (vinv, v)
!
!.....Write the transformation matrices to the symmetry database:
!
      do i = 1, 3
         do j = 1, 3
            or_mat(i,j) = v(i,j)
            or_imat(i,j) = vinv(i,j)
         enddo
      enddo
!
 600  format (/                                               &
       1x, '++SYMORIENT:'/                                    &
       1x, 'Giving a chemical orientation to the molecule.')
 605  format (/                                               &
       1x, "Inertia matrix:"/                                 &
       (1x, 3f25.9))
 610  format (/                                               &
       1x, "Rotation eigenvalues and eigenvectors:"/          &
       (1x, 3f25.9))
 615  format (/                                               &
       1x, 'Newly oriented molecular coordinates:'/           &
       (1x, i5, 3f16.9))
 625  format (/                                               &
       1x, 'ALG(symorient) Molecular orientation matrices: '/ &
       1x, 'Transformation matrix: '/                         &
       3(1x, 3f16.9/)/                                        &
       1x, 'Inverse transformation matrix: '/                 &
       (1x, 3f16.9))
 640  format (/                                               &
       1x, 'END OF SYMORIENT'/)
!
      end subroutine

  !> Send an error message 'message' to stdout, coming from routine
  !> 'routine'. errortype is the error code (see mod_param.f90).
  subroutine error (routine,message,errortype)

    character(len=*), intent(in) :: routine !< Routine calling the error
    character(len=*), intent(in) :: message !< The message
    integer, intent(in) :: errortype !< Fatal, warning or info

    character(len=20) :: chtype

    integer, parameter :: uout = 6
    
    ! message styles shamelessly copied from abinit.
    if (errortype.eq.faterr) then
       chtype='ERROR'
    else if (errortype.eq.warning) then
       chtype='WARNING'
    else if (errortype.eq.noerr) then
       chtype='COMMENT'
    else
       chtype='UNKNOWN'
    endif
    write (uout,100) trim(adjustl(chtype)),trim(adjustl(routine)),&
       trim(adjustl(message))
    if (errortype.eq.faterr) then
       stop  
    else if(errortype.eq.warning) then
       nwarns = nwarns + 1
    else if(errortype.eq.noerr) then
       ncomms = ncomms + 1
    endif

100 format (A,"(",A,"): ",A)

  end subroutine error

  !> Obtains the leng of the string, not counting the blanks at the end.
  integer function leng (string)
        
    character(len=*), intent(in) :: string !< Input string

    integer :: i

    do i=len(string),1,-1
       if (string(i:i).ne.blank .and. string(i:i).ne.null) then
          leng=i
          return
       endif
    end do
    leng=0
    return
  end function leng

  !> Jacobi diagonalization of symmetric matrix.
  !> This algorithm should not be used for very large or nearly
  !> singular matrices.
  subroutine jacobi (a, n, np, d, v, nrot)
          
    implicit none
    integer :: n, np, nrot
    real(kind=8), dimension(np,np) :: a !<matrix to diagonalize
    real(kind=8), dimension(np) :: d !>eigen values vector output
    real(kind=8), dimension(np,np) :: v !>eigen vector matrix output
   
    integer, parameter :: nmax = 3
    integer :: ip, iq, i, j, idum
    real(kind=8), dimension(nmax) :: z, b
    real(kind=8) :: c, g, h, s, sm, t, tau, theta, tresh, xdum

    ! Initialize matrices:
    do ip = 1, n
      do iq = 1, n
        v(ip,iq) = 0d0
      enddo
      v(ip,ip) = 1d0
    enddo
    do ip = 1, n
      b(ip) = a(ip,ip)
      d(ip) = b(ip)
      z(ip) = 0d0
    enddo

    ! Loop over Jacobi rotations:
    nrot = 0
    do i = 1, 50
      ! Check for nullity the upper triangular matrix:
      sm = 0d0
      do ip = 1, n-1
        do iq = ip+1, n
          sm = sm + abs(a(ip,iq))
        enddo
      enddo
      ! If the current matrix is already triangular sort the eigenvalues
      ! and end:
      if (sm.eq.0d0) then
        do ip = 1, n-1
          do iq = ip+1, n
            if (d(iq) .lt. d(ip)) then
              xdum = d(ip)
              d(ip) = d(iq)
              d(iq) = xdum
              do j = 1, n
                xdum = v(j,ip)
                v(j,ip) = v(j,iq)
                v(j,iq) = xdum
              enddo
            endif
          enddo
        enddo
        return
      endif
      ! Find the pivot for a new Jacobi rotation:
      if (i.lt.4) then
        tresh = 0.2*sm/n**2
      else
        tresh = 0d0
      endif
      do ip = 1, n-1
        do iq = ip+1, n
          g = 100*abs(a(ip,iq))
          if ( (i.gt.4) .and. (abs(d(ip))+g.eq.abs(d(ip))) &
                        .and. (abs(d(iq))+g.eq.abs(d(iq))) ) then
            a(ip,iq) = 0d0
          else if (abs(a(ip,iq)) .gt. tresh) then
            h = d(iq)-d(ip)
            if (abs(h)+g .eq. abs(h)) then
              t = a(ip,iq)/h
            else
              theta = 0.5d0*h/a(ip,iq)
              t = 1d0/(abs(theta)+sqrt(1.+theta**2))
              if (theta .lt. 0d0) t = -t
            endif
            c = 1d0/sqrt(1+t**2)
            s = t*c
            tau = s/(1d0+c)
            h = t*a(ip,iq)
            z(ip) = z(ip)-h
            z(iq) = z(iq)+h
            d(ip) = d(ip)-h
            d(iq) = d(iq)+h
            a(ip,iq) = 0d0
            do j = 1, ip-1
              g = a(j,ip)
              h = a(j,iq)
              a(j,ip) = g-s*(h+g*tau)
              a(j,iq) = h+s*(g-h*tau)
            enddo
            do j = ip+1, iq-1
              g = a(ip,j)
              h = a(j,iq)
              a(ip,j) = g-s*(h+g*tau)
              a(j,iq) = h+s*(g-h*tau)
            enddo
            do j = iq+1, n
              g = a(ip,j)
              h = a(iq,j)
              a(ip,j) = g-s*(h+g*tau)
              a(iq,j) = h+s*(g-h*tau)
            enddo
            do j = 1, n
              g = v(j,ip)
              h = v(j,iq)
              v(j,ip) = g-s*(h+g*tau)
              v(j,iq) = h+s*(g-h*tau)
            enddo
            nrot = nrot + 1
          endif
        enddo
      enddo
      do ip = 1, n
        b(ip) = b(ip) + z(ip)
        d(ip) = b(ip)
        z(ip) = 0d0
      enddo
    enddo
    call error("mod_symm","jacobi: 50 iterations should never happen",faterr)

    return

  end subroutine jacobi

end module mod_sym
