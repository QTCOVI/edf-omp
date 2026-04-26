      subroutine pntgrp (ndim,  vector, ngroup, indx, tol, noperx,
     +                   noper, oper)
c
c ======================================================================
c
c purpose:  determine the point group operators
c
c input  :  ndim   - dimensionality of the system
c                    operators will be generated in ndim-dimensional
c                    space
c           vector - vectors (points) in ndim-dimensional space,
c                    describing the system.
c           ngroup - the vectors must be subdivided in ngroup groups.
c                    vectors in different groups are regarded
c                    inequivalent no operators will be generated which
c                    transform vectors from different groups into each
c                    other
c           indx   - index array (ngroup+1). the first vector of group
c                    ig, is no. (indx(ig)+1, and the total number of
c                    vectors is equal to indx(ngroup+1)
c                    furthermore, indx(1) must be zero
c           tol    - tolerance for identity (of points, matrices)
c           noperx - max. nr. of operators allowed
c                    (size of array 'oper').
c                    n.b.: noperx should be at least noper+1, where
c                    noper is the actual nr. of operators which will
c                    be found.
c
c output:   noper  - nr. of point group operators found
c           oper   - matrices of these operators (ndim,ndim,noper)
c
c *=====================================================================
c
*copy vartypes
c
      include        'implicit.inc'
      logical found
c
      parameter (zero  = 0d0,  one = 1d0,
     +           ndimx = 3,    eps = 1d-8)
c
      dimension vector(3,*), indx(ngroup+1)
c
      dimension basis(ndimx,ndimx), bimage(ndimx,ndimx),
     +          dummy(ndimx,ndimx), indgrp(ndimx),
     +          indbas(ndimx),      indimg(ndimx),
     +          result(ndimx),      umat(ndimx,ndimx),
     +          rwk(ndimx**2),    iwk(ndimx)
      dimension  work(ndimx),iwork(ndimx), xdet(2), rwork(ndimx)
      dimension  oper(3,3,noperx)
c
c ======================================================================
c
      if (ndim.lt.1 .or. ndim.gt.ndimx) stop 'ndim out of range'
c
      if (ngroup .le. 0) stop 'NGROUP LE 0. PNTGRP'
      if (indx(1) .ne. 0) stop 'FIRST INDEX.NE.0. PNTGRP'
c
c ======================================================================
c
c set unit matrix
c
      do 10 i = 1, ndimx
      do 10 j=1,ndimx
   10 umat(i,j) = zero
      umat(1,1)=one
      umat(2,2)=one
      umat(3,3)=one
c
c the first operator is the identity.
c
      noper = 1
      do 30 j = 1, ndim
        do 20 i = 1, ndim
          oper(i,j,noper) = umat(i,j)
   20   continue
   30 continue
c
c find an appropriate group combination (for basis and image)
c
c initialize the group indicators
c
      do 40 i = 1, ndim-1
   40 indgrp(i) = 1
      indgrp(ndim) = 0
c
c set the indicators for the next combination of groups.
c
   50 do 70 j = ndim, 1, -1
        if (indgrp(j) .lt. ngroup) then
          indgrp(j) = indgrp(j) + 1
          do 60 k = j+1, ndim
            indgrp(k) = indgrp(j)
   60     continue
          goto 80
        endif
   70 continue
      stop 'ERROR 1. PNTGRP'
c
c check if all indicated groups contain enough points.
c
   80 do 100 i = 1, ndim
        ig = indgrp(i)
        ntel = 1
        do 90 j = i+1, ndim
          if (indgrp(j).eq.ig) ntel = ntel+1
   90   continue
        if (ntel.gt.indx(ig+1)-indx(ig)) goto 50
  100 continue
c
c take a standard set (basis), to be transformed into the other
c sets (image).
c initialize the indicators: which vectors from the group will
c be in the basis.
c
      indbas(1) = 1
      do 110 i = 2, ndim
        indbas(i) = 1
        if (indgrp(i).eq.indgrp(i-1)) indbas(i) = indbas(i-1) + 1
  110 continue
c
c basis points.
c
  120 do 140 ivec = 1, ndim
        ig = indgrp(ivec)
        do 130 i = 1, ndim
          basis(i,ivec) = vector(i, indx(ig)+indbas(ivec))
  130   continue
  140 continue
c
c check linear independency of the basis.
c
      if (ndim.eq.1) then
        s = abs (basis(1,1))
      elseif (ndim.eq.2) then
        s = abs (basis(1,1) * basis(2,2) - basis(1,2) * basis(2,1))
      elseif (ndim.eq.3) then
        s = abs(r3volp(basis(1,1),basis(1,2),basis(1,3)))
      endif
c
      if (s .gt. tol) then
*           write (6,*) 'b', (indbas(ivec),ivec=1,3)
*        do i=1,3
*           write (6,*) (basis(i,j),j=1,3)
*        enddo
         goto 200
      endif
c
c the indicators for the next basis
c
  170 do 190 j = ndim, 1, -1
        ig = indgrp(j)
        if (indbas(j) .lt. indx(ig+1)-indx(ig)) then
          indbas(j) = indbas(j)+1
          do 180 k = j+1,ndim
            indbas(k) = 1
            if (indgrp(k).eq.indgrp(k-1)) then
              indbas(k) = indbas(k-1) + 1
              kg = indgrp(k)
              if (indbas(k).gt.indx(kg+1)-indx(kg)) goto 170
            endif
  180     continue
          goto 120
        endif
  190 continue
c
c all basis possibilities for this group combination have been
c exhausted apparently this combination is impossible.
c
      goto 50
c
c invert the basis
c a symmetry operator is defined by:
c        op * basis = set
c    :   op         = set * basis inverted
c
c
 200  continue
      call mydgeco (basis,ndimx,ndim,iwork,rcond,rwork)
      if (rcond.lt.tol)  stop 'ERROR2 FROM DGECO. PNTGRP'
      call mydgedi (basis,ndimx,ndim,iwork,xdet,rwork,1)
c
* 200 call minvr (ndimx, ndim, eps, det, basis, rwk, iwk, ier)
*     if (ier.ne.0) stop 'ERROR2 FROM MINVR. PNTGRP'

c
c a proper basis has been obtained now (and inverted).
c
c
c
c loop over the other sets, the basis might be transformed into:
c the images (or: projections).
c
c initialize the indicators.
c
      do 210 i = 1,ndim
  210 indimg(i) = indbas(i)
c
c the unit operator (identity) has already been obtained,
c take immediately the next projection.
c
      goto 430
c
c image points.
c
  220 do 240 ivec = 1,ndim
        ig = indgrp(ivec)
        do 230 i = 1,ndim
  230   bimage(i,ivec) = vector(i, indx(ig)+indimg(ivec))
  240 continue
c
c check linear independency of the projection.
c
      if (ndim.eq.1) then
        s = abs (bimage(1,1))
      elseif (ndim.eq.2) then
        s = abs (bimage(1,1)*bimage(2,2) - bimage(1,2)*bimage(2,1))
      elseif (ndim.eq.3) then
        s = abs (r3volp (bimage(1,1), bimage(1,2), bimage(1,3)))
      endif
c
      if (s .lt. tol) goto 430
c
c construct the operator which transforms this basis into
c this projection.
c
*           write (6,*) 'i', (indimg(ivec),ivec=1,3)
*        do i=1,3
*           write (6,*) (bimage(i,j),j=1,3)
*        enddo
c
      do 300 j = 1, ndim
        do 270 i = 1, ndim
  270   dummy(i,j) = bimage(i,1) * basis(1,j)
        do 290 k = 2, ndim
          do 280 i = 1, ndim
  280     dummy(i,j) = dummy(i,j) + bimage(i,k) * basis(k,j)
  290   continue
  300 continue
c
c
c check unitarity of the (real) operator.
c
      do 330 j = 1,ndim
        do 320 i = 1,ndim
          s = zero
          do 310 k = 1,ndim
  310     s = s  + dummy(k,i) * dummy(k,j)
          if (abs (s-umat(i,j)) .gt. tol) goto 430
  320   continue
  330 continue
c
c check if all vectors are transformed into each other by this
c operator.
c the vectors that are transformed into each other must belong
c to the same group.
c
      do 400 igroup = 1, ngroup
        do 390 ivec = indx(igroup)+1, indx(igroup+1)
c
        do 340 i = 1, ndim
  340   result(i) = dummy(i,1) * vector(1,ivec)
        do 360 k = 2, ndim
          do 350 i = 1, ndim
            result(i) = result(i) + dummy(i,k) * vector(k,ivec)
  350     continue
  360   continue
c
c check if the result vector occurs in the list (and belongs
c to the same group).
c
        do 380 jvec = indx(igroup)+1,indx(igroup+1)
          dev = zero
          do 370 i = 1, ndim
  370     dev = dev + abs (result(i)-vector(i,jvec))
          if (dev .lt. tol) goto 390
  380   continue
c
c the result vector has not been found in the list.
c
        goto 430
c
  390   continue
  400 continue
c
c add the new operator.
c
      do 420 j = 1, ndim
        do 410 i = 1, ndim
  410   rwk((j-1)*ndim+i) = dummy(i,j)
  420 continue
c
      call symadd (ndim, noperx, noper, oper, rwk)
c
c the indicators for the next image.
c
  430 continue
      do 450 j = ndim, 1, -1
        ig = indgrp(j)
        if (indimg(j).lt.indx(ig+1)-indx(ig)) then
          indimg(j) = indimg(j) + 1
          do 440 k = j+1, ndim
            indimg(k) = 1
            if (indgrp(k).eq.indgrp(k-1) .and.
     +          indimg(k).eq.indimg(k-1)) then
              indimg(k) = indimg(k-1) + 1
              kg = indgrp(k)
              if (indimg(k).gt.indx(kg+1)-indx(kg)) goto 430
            endif
  440     continue
          goto 220
        endif
  450 continue
c
c
c (not necessary:)
c if the inversion occurs in the set of symmetry operations,
c place it in the second postion
c
      do iop = 3, noper
c
        do j = 1, ndim
          do i = 1, ndim
            atmp=oper(i,j,iop)
            if (abs (atmp+umat(i,j)) .gt. tol) goto 500
*           if (abs (oper(i,j,iop)+umat(i,j)) .gt. tol) goto 500
          enddo
        enddo
c
c inversion found, interchange it with the second operator
c
        do j = 1, ndim
          do i = 1, ndim
            oper(i,j,iop) = oper(i,j,2)
            oper(i,j,2) = -umat(i,j)
          enddo
        enddo
        goto 510     !exit
  500   continue
      enddo
  510 continue
c
c check if created set of operators constitutes a group.
c
      call grpchk (ndim, noper, noperx, oper, tol, ier)
      if (ier.ne.0) stop 'GROUPCHECK ERROR. PNTGRP'
c
      return
      end
c-------------------------------------------------------
      subroutine symadd (ndim, noperx, noper, oper, operad)
c
c ======================================================================
c
c purpose:  add a symmetry operator to a group of symmetry operators,
c           and add also all product operators to make the output set
c           a group, supposing that the input set is one.
c
c input  :  ndim   - dimensionality of the operator matrices (ndim,ndim)
c           noperx - max. nr. of operators allowed (size of array oper)
c           operad - new operator, to be added to the group
c in-out :  noper  - in : initial nr. of symmetry operators
c                    out: final nr. of symmetry operators
c           oper   - operators
c
c *=====================================================================
c
*copy vartypes
c
      include        'implicit.inc'
      logical found
c
      character  rtname * (*)
      parameter (rtname = 'SYMADD (HFSLIB)')
c
      parameter (eps    = 1.0d-3,
     +           ndimx  = 4,
     +           zero   = 0.0d0)
c
c
      dimension  oper   (3,3,noperx),
     +           operad (3,3),operadn(3,3)
c
      dimension  prod   (ndimx,ndimx)
c
c ======================================================================
c
      if (ndim.eq.3) then
        do i=1,ndim
          do j=1,ndim
            operadn(i,j)=operad(i,j)
          enddo
        enddo
      elseif (ndim.eq.2) then
        operadn(1,1)=operad(1,1)
        operadn(2,1)=operad(2,1)
        operadn(1,2)=operad(3,1)
        operadn(2,2)=operad(1,2)
      else
        stop 'symadd: !! ndim is not 2 or 3 !!'
      endif
c
      if (noper.gt.noperx) then
        stop 'Error: symadd (noper,noperx)  '
      endif
c
      if (ndim.gt.ndimx) then
        stop 'ndim exceeded'
      endif
c
c search for the new operator in the old set
c
      do 20 ioper = 1,noper
c
        do 15 j = 1,ndim
          do 10 i = 1,ndim
            if (abs (oper(i,j,ioper) - operadn(i,j)).gt.eps) goto 20
   10     continue
   15   continue
c
c operator found, routine can be stopped
c
      return
c
   20 continue
c
c the new operator must be added
c
      noper = noper + 1
      if (noper.gt.noperx) then
        stop 'NOPERX TOO SMALL'
      endif
c
      do 30 j = 1,ndim
        do 25 i = 1,ndim
          oper(i,j,noper) = operadn(i,j)
   25   continue
   30 continue
c
c for each new operator (joper, starting at the current value of
c noper) consider the right and left products (as the operators
c not necessarily commute) with each (ioper) of the previous
c operators.
c
      joper = noper - 1
   40 joper = joper + 1
c
        ioper = 0
   50   ioper = ioper + 1
c
          do 70 j = 1,ndim
            do 65 i = 1,ndim
              s = zero
              do 60 k = 1,ndim
                s = s + oper(i,k,ioper) * oper(k,j,joper)
   60         continue
              prod(i,j) = s
   65       continue
   70     continue
c
c check if the computed product operator is already
c in the set.
c
          do 90 koper = 1,noper
            s = zero
            do 85 j = 1,ndim
              do 80 i = 1,ndim
                s = s + abs (oper(i,j,koper) - prod(i,j))
   80         continue
   85       continue
            if (s.lt.eps) goto 110
   90     continue
c
c product operator is not in the set, it will be added now
c
          noper = noper + 1
          if (noper.gt.noperx) then
            stop ' NOPERX TOO SMALL(2)'
          endif
c
          do 105 j = 1,ndim
            do 100 i = 1,ndim
              oper(i,j,noper) = prod(i,j)
  100       continue
  105     continue
c
c try also the commuted product
c
  110     do 130 j = 1,ndim
            do 125 i = 1,ndim
              s = zero
              do 120 k = 1,ndim
                s = s + oper(i,k,joper) * oper(k,j,ioper)
  120         continue
              prod(i,j) = s
  125       continue
  130     continue
c
c check if the computed product operator is already
c in the set.
c
          do 150 koper = 1,noper
            s = zero
            do 145 j = 1,ndim
              do 140 i = 1,ndim
                s = s + abs (oper(i,j,koper) - prod(i,j))
  140         continue
  145       continue
            if (s.lt.eps) goto 170
  150     continue
c
c product operator is not in the set, it will be added now
c
          noper = noper + 1
          if (noper.gt.noperx) then
            stop 'noperx too small'
          endif
c
          do 165 j = 1,ndim
            do 160 i = 1,ndim
              oper(i,j,noper) = prod(i,j)
  160       continue
  165     continue
c
  170     continue
c
        if (ioper.lt.noper) goto 50
c
      if (joper.lt.noper) goto 40
c
      return
      end
c-------------------------------------------------------
      subroutine grpchk (ndim, noper, noperx, oper, tol, ier)
c
c ======================================================================
c
c purpose:  check if a set of real ndim*ndim - matrices constitutes
c           a group. the matrices are supposed to be orthonormal.
c
c input  :  ndim  - dimension of the matrices
c           noper - nr. of matrices
c           oper  - matrices
c           tol   - tolerance parameter (identification of operators)
c
c output :  ier   - 0: group
c                   1: not a group
c
c remarks:  for each operator, two checks are made:
c           1. is its transpose in the set
c           2. is the product with each operator (incl. itself) in
c              the set
c
c *=====================================================================
c
*copy vartypes
c
      include        'implicit.inc'
      logical found
c
      parameter (zero = 0d0)
c
      dimension oper(3,3,noperx)
      dimension prod(3,3)
c
c ======================================================================
c
      if (ndim  .gt. 3) stop 'NDIM TOO LARGE. GRPCHK'
      if (noper .le. 0) stop 'NO OPERATORS. GRPCHK'
      ier = 0
c
      do 120 ioper = 1, noper
c
c check if its inverse is in the set
c for real, orthonormal matrices, the inverse is equal to
c the transpose
c
        do 30 joper = 1, noper
          s = zero
          do 20 j = 1, ndim
            do 10 i = 1, ndim
              s = s + abs (oper(i,j,ioper) - oper(j,i,joper))
   10       continue
   20     continue
          if (s .le. tol) goto 40
   30   continue
c
c the inverse is not found
c
        ier = 1
        return
c
c
c check the product with each operator
c
   40   do 110 joper = 1, noper
c
          do 70 j = 1, ndim
            do 60 i = 1, ndim
              s = zero
              do 50 k = 1, ndim
                s = s + oper(i,k,ioper) * oper(k,j,joper)
   50         continue
              prod(i,j) = s
   60       continue
   70     continue
c
c search among the operators for the product
c
          do 100 koper = 1, noper
c
            s = zero
            do 90 j = 1, ndim
              do 80 i = 1, ndim
                s = s + abs (oper(i,j,koper) - prod(i,j))
   80         continue
   90       continue
c
            if (s .le. tol) goto 110
  100     continue
c
c product not found
c
          ier = 1
          return
c
  110   continue
  120 continue
c
      return
      end
c--------------------------------------------------------
*compile v=0
      function r3volp (a, b, c)
c
c ======================================================================
c
c purpose: volume product in 3d space: (a*b).c
c
c remark * the result may be negative.
c
c *=====================================================================
c
*copy vartypes
      include        'implicit.inc'
c
      dimension  a(3), b(3), c(3)
c
c ======================================================================
c
c volume product (determinant)
c
      r3volp = (a(2)*b(3) - a(3)*b(2)) * c(1) +
     +         (a(3)*b(1) - a(1)*b(3)) * c(2) +
     +         (a(1)*b(2) - a(2)*b(1)) * c(3)
c
      return
      end
c.............................................................
      subroutine vecprod (a, b, c)
c
      include        'implicit.inc'
c
      dimension  a(3), b(3), c(3)
c
c ======================================================================
c
c     vector product
c
      c(1) = a(2)*b(3) - a(3)*b(2) 
      c(2) = a(3)*b(1) - a(1)*b(3) 
      c(3) = a(1)*b(2) - a(2)*b(1)
c
      return
      end
c
c....................................................
      subroutine rzero(a,n)
      include        'implicit.inc'
      dimension  a(n,n)
      do i=1,n
      do j=1,n
          a(i,j)=0d0
      enddo
      enddo
      end
c....................................................
      subroutine riden(a,n)
      include        'implicit.inc'
      dimension  a(n,n)
      do i=1,n
      do j=1,n
          if (i.ne.j) then
             a(i,j)=0d0
          else
             a(i,i)=1d0
          endif
      enddo
      enddo
      end
c....................................................
      subroutine rmul(a,b,c,n)
      include        'implicit.inc'
      dimension  a(n,n), b(n,n), c(n,n)
      do j=1,n
         do i=1,n
            c(i,j)=0d0
            do k=1,n
               c(i,j)=c(i,j)+a(i,k)*b(k,j)
            enddo
         enddo
      enddo
      end
c....................................................
      subroutine requal(a,b,n)
      include        'implicit.inc'
      dimension  a(n,n), b(n,n)
      do j=1,n
         do i=1,n
            a(i,j)=b(i,j)
         enddo
      enddo
      end
c....................................................
      subroutine triprod(a,b,c,n)
c     a^t b a = c
      include        'implicit.inc'
      dimension  a(n,n), b(n,n), c(n,n)
      do i=1,n 
         do j=1,n
            c(i,j)=0d0
            do k=1,n
               do l=1,n
                  c(i,j)=c(i,j)+
     &            a(k,i)*b(k,l)*a(l,j)
               enddo
            enddo
         enddo
      enddo
      end
