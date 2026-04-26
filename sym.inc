c
c.....sym.inc - Symmetry variables include
c 
      integer(kind=4), parameter :: maxis=12
      integer(kind=4) :: nop,invop(48)
      real(kind=8) :: oper(3,3,48),oper2(3,3,48)
      real(kind=8) :: eulang(3,48),vaxis(3,48),vaxi
      integer(kind=4) :: ioper(48,3),masksym(-1:1,maxis)
      logical :: linear,planar,nosym
      integer(kind=4) :: nmono,mnnuc,neq,npair,iopcurr,neqsurf
      integer(kind=4) :: neqscalc
c-----------------------------------------------------------------
c....................... oper: symmetry matrices
c.......................ioper: indices : 1) (+1)proper or (-1)improper
c                                        2   order of axis (rot-refl)
c                                       3,4  the same (rot-inv)
c                       masksym: auxiliar. number of axes of each type.
c..................................................................
c     iopcurr: the inverse of this operation is the one from which
c              the current RSH rotation matrices are built
c..................................................................
      common /isyme/ ioper,masksym,nop,neq,npair,invop,nmono,mnnuc,
     &               iopcurr,neqsurf,neqscalc
      common /rsyme/ oper,oper2,eulang,vaxi,vaxis
      common /lsyme/ linear,planar,nosym
c
c.....Tolerance for the identity of point, matrices, ...
c
      real(kind=8) :: tolsym,epssym,epsorb
      common /symtol/ tolsym,epssym,epsorb
c
c.....Manually set point group
c
      logical :: lpointgroup
      character(len=3) :: cpointgroup
      common /symm1/ lpointgroup
      common /symm2/ cpointgroup
