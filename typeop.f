      subroutine   typeop (rotm,tipo,order,order1,repet,euler,vec)
c
c.....Determines :
c            a) Order and directions of proper and improper axes.
c            b) Orientation of planes.
c     The matrix elements will be cleaned from numerical errors.
c
      include           'implicit.inc'
      include           'param.inc'
c
      parameter         (zero=0d0,one=1d0)
c
c.....Local variables
      integer           tipo, order, iord(2,2)
      integer           tipo1, order1, repet
      dimension         vec (3), rotm(3,3)
      real(kind=8)      mat,mati, aux(3)
      real(kind=8)      euler(3)
      dimension         mat(3,3), mati(3,3), eigen(3), eigeni(3), 
     &                  vprop (3,3), vpropi(3,3), fv1(3), fv2(3),
     &                  fv3(3)
      parameter         (epsilon=1d-6)
      logical           eqmat,notinteg
c
c.....We use eispack routine to diagonalize orthogonal but not 
c     necessarily symmetric matrices. Variables ending in "i" refer
c     to the imaginary parts of a complex result.
c
c
c.....Compute the trace and copy the rotation matrix to local var.
c-----------------------------------------------------------------
      trace=zero
      do i = 1,3
         trace = trace + rotm (i,i)
         fv1(i)=zero
         fv2(i)=zero
         fv3(i)=zero
         eigen(i)=zero
         eigeni(i)=zero
         vec(i)=zero
      enddo
      call requal(mat,rotm,3)
      call rzero(mati,3)
      call rzero(vprop,3)
      call rzero(vpropi,3)
      do i=1,2
         do j=1,2
            iord(i,j)=0
         enddo
      enddo
c
c.....Determine type of operation
c
      if (abs(trace+3) .lt. epsilon) then 
c
c........Inversion
c
         tipo = -1
         order= 2
         order1= 1
         repet=1
         call rzero(rotm,3)
         rotm(1,1)=-one
         rotm(2,2)=-one
         rotm(3,3)=-one
         vec(3)=one
         call geteuler(vec,euler,order1)
         return
      elseif (abs(trace-3) .lt.epsilon) then   
c
c........Identity
c
         tipo = 1
         order=1
         order1=1
         repet=1
         vec(3)=one
         call riden(rotm,3)
         call geteuler(vec,euler,order1)
         return
      else
c
         notinteg=.false.
         nm=3
c
c........It is an axis or a plane. Diagonalize
c
         call cg(nm,nm,mat,mati,eigen,eigeni,1,vprop,vpropi,
     &           fv1,fv2,fv3,ierr)
         if (ierr .ne. 0) then
            stop 'typeop.f: Error on diagonalizing symmetric matrix'
         endif
c
c........Determine actual operation.
c
         nones = 0
         nminusones = 0 
         do i = 1, 3
            if (abs(eigen (i)-one) .lt. epsilon .and.
     &          abs(eigeni(i))   .lt. epsilon ) then
               nones = nones + 1
               iord(1,nones) = i
            elseif (abs(eigen(i)+one) .lt. epsilon .and.
     &              abs(eigeni(i))  .lt. epsilon ) then
               nminusones = nminusones +1
               iord(2,nminusones) = i
            endif
         enddo
c
c........What is it?
c
         if (nones .eq. 1) then
c
c..........A proper axes, of course. Order?
c
           tipo=1
           temp=(trace-one)/2d0
           if (abs(temp).gt.1d0) temp=sign(1d0,temp)
           xord=abs(2d0*pi/(acos(temp)))
           order=nint(xord)
           iorder=int(xord+1d-6)
           if (abs(xord-iorder).gt.1d-3) notinteg=.true.
           order1= order
           repet=1
           do i = 1,3
              vec(i) = vprop (i,iord(1,1))
           enddo
            
         elseif (nones .eq. 2) then
c
c........A plane, isnt it? Find directions.
c
           tipo=-1
           order=1
           order1=2
           repet=1
           do i = 1,3
              vec(i) =vprop (i,iord(2,1))
           enddo
         elseif (nminusones .eq. 1) then
c
c...........An improper axes, of course. Order?
c
            tipo=-1
            temp=(trace+one)/2d0
            if (abs(temp).gt.1d0) temp=sign(1d0,temp)
            xord = abs(2d0*pi/(acos(temp)))
            order = nint(xord)
            iorder = int(xord+1d-6)
            if (abs(xord-iorder).gt.1d-3) notinteg=.true.
            order1= -order/(1+mod(abs(order/2),2))
            repet=1
            do i = 1,3
               vec(i) = vprop (i,iord(2,1))
            enddo
         else
            stop  'Symetry Element unknown'
         endif
c
c........cleaning strategies.
c        clean axis vector
c        reconstruct rotation
c        if it does not work, it is the inverse operation
c        clean rotation
c
         call cleanvec(vec)
         call cleanvec (rotm(1,1))
         call cleanvec (rotm(1,2))
         call cleanvec (rotm(1,3))
         call requal (mat,rotm,3)
 111     call rotz (rotm, vec)
         if (.not.notinteg) then
           call genoper  (rotm, tipo, order)
           repet=1
         else
           call genoper2 (rotm, tipo, xord, order, repet)
           if (nones .eq. 1) then
             order1= order
           elseif (nones .eq. 2) then
             order1=2
           elseif (nminusones .eq. 1) then
             order1= -order/(1+mod(abs(order/2),2))
           endif
         endif
         if (.not.eqmat (rotm,mat,3)) then
             order=-order
             xord=-xord
             goto 111
         endif
         call cleanvec (rotm(1,1))
         call cleanvec (rotm(1,2))
         call cleanvec (rotm(1,3))
c
c........pass the sign to rotation -inversion
c
         order1=order1*sign(1,order)
c
c........construct euler angles
c
         call geteuler(vec,euler,order1)
c         
         return 
      endif
c
c.....formats
      end
c---------------------------------------------------------
      subroutine writemat(mat,label)
      include           'implicit.inc'
      character*(*)     label
      real(kind=8)      mat(3,3)
      write (6,*) label
      do i=1,3
         write (6,'(3f12.6)') (mat(i,j),j=1,3)
      enddo
      write (6,*)
      end
 
c---------------------------------------------------------
      subroutine geteuler(vec,euler,order)
      include           'implicit.inc'
      parameter        (eps=1d-6)
c
      integer          order
      real(kind=8)     vec(3), euler(3)
c
c     (cos alp, sin alp, cos bet, sin bet, cos gam, sin gam)
c
      pi=acos(-1d0)
      rxy2=vec(1)*vec(1)+vec(2)*vec(2)
      rxy=sqrt(rxy2)
      sth=rxy
      cth=vec(3)
      if (rxy.lt.eps) then
         ph=0d0
      else
        if (sth.ge.0d0) then
           ph=acos(vec(1)/rxy)
        else
           ph=2d0*pi-acos(vec(1)/rxy)
        endif
      endif
      if (order.eq.1) then
         omg2=0d0
      else
         omg2=pi/dble(order)
      endif
      beta=2d0*asin(sth*sin(omg2))
      amg=2d0*ph-pi
      apg=2d0*atan(cth*tan(omg2))
      alpha=0.5d0*(apg+amg)
      gamma=0.5d0*(apg-amg)
      euler(1)=alpha
      euler(2)=beta
      euler(3)=gamma
      end
c
c---------------------------------------------------------
      logical function  eqmat(a,b,n)
      include           'implicit.inc'
      dimension         a(n,n),b(n,n)
      parameter         (eps=1d-4)
c
      eqmat=.true.
      do i=1,n
         do j=1,n
            if (abs(b(i,j)-a(i,j)).gt.eps) then
               eqmat=.false. 
               return
            endif
         enddo
      enddo
      end
c---------------------------------------------------------
      subroutine rotz(rotm,vec)
c
c.....find operation to turn vec axis into z axis
c
      include           'implicit.inc'
c
      real(kind=8)     rotm(3,3), vec(3)
      real(kind=8)     raux(3,3),raux2(3,3)
      parameter        (eps=1d-6)
c
      rxy=sqrt(vec(1)*vec(1)+vec(2)*vec(2))
      if (rxy.lt.eps) then
c        no orientation
         call riden(rotm,3)
      else if (abs(vec(1)).lt.eps) then
         cost=vec(3)
         sint=vec(2)
         call riden(rotm,3)
         rotm(2,2)=cost
         rotm(2,3)=-sint
         rotm(3,2)= sint
         rotm(3,3)=cost
      else
         cost=vec(2)/rxy
         sint=vec(1)/rxy
         call riden(rotm,3)
         rotm(1,1)=cost
         rotm(1,2)=-sint
         rotm(2,1)= sint
         rotm(2,2)=cost
         cost=vec(3)
         sint=rxy
         call riden(raux,3)
         raux(2,2)=cost
         raux(2,3)=-sint
         raux(3,2)= sint
         raux(3,3)=cost
         call rmul(raux,rotm,raux2,3)
         call requal(rotm,raux2,3)
      endif
      end
c
c---------------------------------------------------------
c
      subroutine genoper(rotm,tipo,order)
c
c.....generate clean operation
c
      include           'implicit.inc'
      include           'param.inc'
      real(kind=8)      rotm(3,3)
      real(kind=8)      raux(3,3),raux2(3,3)
      integer            tipo,order
c
       
      sig=dble(tipo)*1d0
      is=sign(1,order)
      aorder=abs(order)
      iaorder=int(aorder)
      if (iaorder.eq.1) then
         cost=1d0
         sint=0d0
      elseif(iaorder.eq.2) then
         cost=-1d0
         sint=0d0
      elseif(iaorder.eq.3) then
         cost=-0.5d0
         sint=sqrt(3d0)/2d0*is
      elseif(iaorder.eq.4) then
         cost=0d0
         sint=1d0*is
      elseif(iaorder.eq.5) then
         cost=(sqrt(5d0)-1d0)/4d0
         sint=sqrt(1d0-cost*cost)*is
      elseif(iaorder.eq.6) then
         cost=0.5d0
         sint=sqrt(3d0)/2d0*is
      elseif(iaorder.eq.8) then
         cost=1d0/sqrt(2d0)
         sint=cost*is
      else
         phi=2d0*pi/dble(order)
         cost=cos(phi)
         sint=sin(phi)
      endif
c
      call rzero(raux,3)
      raux(1,1)=cost
      raux(2,2)=cost
      raux(1,2)=-sint
      raux(2,1)= sint
      raux(3,3)=sig
c
      call triprod(rotm,raux,raux2,3)
      call requal(rotm,raux2,3)
c
      end


c
c---------------------------------------------------------
c
      subroutine genoper2 (rotm,tipo,xord,order,repet)
c
c.....generate clean operation
c
      include           'implicit.inc'
      include           'param.inc'
      real(kind=8)      rotm(3,3)
      real(kind=8)      raux(3,3),raux2(3,3)
      integer           tipo,order,repet
      logical           foundaxis
      real(kind=8)      xnm(22)
      integer            in(22),im(22)
      data xnm / 2.1428571428571432d0,2.1666666666666665d0,
     &           2.2000000000000000d0,2.2500000000000000d0,
     &           2.3333333333333330d0,2.4000000000000000d0,
     &           2.5000000000000000d0,2.6000000000000000d0,
     &           2.6666666666666666d0,2.7500000000000000d0,
     &           2.8000000000000000d0,3.2500000000000000d0,
     &           3.3333333333333333d0,3.5000000000000000d0,
     &           3.6666666666666667d0,3.7500000000000000d0,
     &           4.3333333333333333d0,4.5000000000000000d0,
     &           4.6666666666666667d0,5.5000000000000000d0,
     &           6.5000000000000000d0,7.5000000000000000d0/
      data in  / 15,13,11, 9, 7,12, 5,13, 8,11,14,13,10, 7,11,15,
     &           13, 9,14,11,13,15/
      data im  / 7,6,5,4,3,5,2,5,3,4,5,4,3,2,3,4,3,2,3,2,2,2/ 
c
      foundaxis=.false.
      do i=1,22
        if (abs(abs(xord)-xnm(i)).lt.1d-4) then
          foundaxis=.true.
          if (xord.gt.0d0) then
            order=+in(i)
            repet=+im(i)
            xdiv=+dble(in(i))/dble(im(i))
          else
            order=-in(i)
            repet=+im(i)
            xdiv=-dble(in(i))/dble(im(i))
          endif
        endif
        if (foundaxis) goto 10
      enddo
 10   continue
      if (.not.foundaxis) stop 'genoper:  !!! Axis not found !!!'
       
      sig=dble(tipo)*1d0
      cost=cos(2d0*pi/xdiv)
      sint=sin(2d0*pi/xdiv)
      call rzero(raux,3)
      raux(1,1)=cost
      raux(2,2)=cost
      raux(1,2)=-sint
      raux(2,1)= sint
      raux(3,3)=sig
c
      call triprod(rotm,raux,raux2,3)
      call requal(rotm,raux2,3)
c
      end

c---------------------------------------------------------
      subroutine cleanvec(vec)
c
      include           'implicit.inc'
      parameter         (eps=1d-6,nrat=20)
c
      dimension         vec(3)
      dimension         abv(3), siv(3), ifrac(3)
c.....cleans vector from inaccuracies
      nfrac=0
      ifrac(1)=0
      ifrac(2)=0
      ifrac(3)=0
      vnorm=sqrt(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3))
      do i=1,3
         vec(i)=vec(i)/vnorm
         abv(i)=abs(vec(i))
         siv(i)=sign(1d0,vec(i))
         if (abv(i).lt.eps) then
            vec(i)=0d0
            nfrac=nfrac+1
            ifrac(i)=1
         else
            do n=2,nrat
               temp=abv(i)*dble(n)
               itemp=nint(temp)
               if (abs(temp-itemp).lt.eps) then
                  vec(i)=siv(i)*dble(itemp)/dble(n)
                  nfrac=nfrac+1
                  ifrac(i)=1
                  goto 10
               endif
               temp=abv(i)**2*dble(n)
               itemp=nint(temp)
               if (abs(temp-itemp).lt.eps) then
                  vec(i)=siv(i)*sqrt(dble(itemp)/dble(n))
                  nfrac=nfrac+1
                  ifrac(i)=1
                  goto 10
               endif
            enddo
         endif
 10      continue
      enddo
      if (nfrac.eq.2) then
         do i=1,3
            if (ifrac(1).eq.0) then
               sum=vec(2)*vec(2)+vec(3)*vec(3)
               vec(1)=sqrt(1d0-sum)
            elseif (ifrac(2).eq.0) then
               sum=vec(1)*vec(1)+vec(3)*vec(3)
               vec(2)=sqrt(1d0-sum)
            elseif (ifrac(3).eq.0) then
               sum=vec(2)*vec(2)+vec(1)*vec(1)
               vec(3)=sqrt(1d0-sum)
            endif
         enddo
         nfrac=3
      endif
      do i=1,3
         abv(i)=abs(vec(i))
         siv(i)=sign(1d0,vec(i))
      enddo
      if (nfrac.ne.3) then 
c.....two equal ?
         if (abs(abv(1)-abv(2)).lt.eps) then
            vec(1)=siv(1)*abv(1)
            vec(2)=siv(2)*abv(1)
            vec(3)=siv(3)*sqrt(1d0-2d0*vec(1)*vec(1))
         elseif (abs(abv(1)-abv(3)).lt.eps) then
            vec(1)=siv(1)*abv(1)
            vec(3)=siv(3)*abv(1)
            vec(2)=siv(2)*sqrt(1d0-2d0*vec(1)*vec(1))
         elseif (abs(abv(2)-abv(3)).lt.eps) then
            vec(2)=siv(2)*abv(2)
            vec(3)=siv(3)*abv(2)
            vec(1)=siv(1)*sqrt(1d0-2d0*vec(2)*vec(2))
         endif
      endif
      end
