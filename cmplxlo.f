c
c-----------------------------------------------------------------------
c
      subroutine cmplxlo (sg,c,xloci,nmo,ngroup,lw,okw)
c
c.....ruedmis routine for complex sg() integrals. 
c-----------------------------------------------------------------------
c                         
      include 'implicit.inc'
      include 'constants.inc'
      include 'error.inc'
c
      parameter (maxit=200)

      real(kind=8)   c(nmo,nmo)
      real(kind=8)   xloci(nmo)
      real(kind=8)   pi,pi4
      complex*16  sg(ngroup,nmo,nmo)
      complex*16, allocatable,dimension (:,:,:)  :: sgtmp
      real(kind=8), allocatable,dimension (:)    :: xion
      integer, allocatable,dimension (:)         :: ixion
      real(kind=8), allocatable,dimension (:,:)  :: ctmp
      real(kind=8), allocatable,dimension (:)    :: delta,deltad,rotinc
      integer, allocatable,dimension (:)         :: indi
      real(kind=8)      spp,sqq,spq,temp
      complex*16  aom1,opp,oqq,opq,oqp,tmp1,tmp2

      logical   warn,okw,goon
      logical   moreiter
      integer   ipair(2)
c
      call timer (2,ixlo,'_cmplxlo  ',-1)
c
c.....Allocate arrays
c
      n    = nmo
      if (.not.allocated(sgtmp)) then
        allocate (sgtmp(ngroup,nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot allocate sgtmp()'
      endif

      if (.not.allocated(xion)) then
        allocate (xion(nmo),stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot allocate xion()'
      endif

      if (.not.allocated(ixion)) then
        allocate (ixion(nmo),stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot allocate ixion()'
      endif

      if (.not.allocated(ctmp)) then
        allocate (ctmp(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot allocate ctmp()'
      endif

      if (.not.allocated(delta)) then
        allocate (delta(ngroup*(ngroup-1)/2),stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot allocate delta()'
      endif

      if (.not.allocated(deltad)) then
        allocate (deltad(ngroup*(ngroup-1)/2),stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot allocate deltad()'
      endif

      if (.not.allocated(rotinc)) then
        allocate (rotinc(nmo*(nmo-1)/2),stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot allocate rotinc()'
      endif

      if (.not.allocated(indi)) then
        allocate (indi(nmo),stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot allocate indi()'
      endif

      pi4  = atan(one)
      four = two+two
      pi   = four*pi4
      warn =.false.
c
c.....Initial rotation matrix is the unit matrix.
c
      c(1:n,1:n)=zero
      do i = 1, n
        c(i,i)=one
      enddo
c
c.....Start process.
c
      moreiter = .true.
      maxitn=maxit*n*(n-1)/2
      nrot=1
      goon=.true.
      do while (nrot.le.maxitn.and.goon)
        ipq=0
c
c.......Compute and print (DE)Localization funtion before rotating.
c
        xloctot=zero
        do i=1,n
          xloci(i)=zero
          do k = 1, ngroup
            spp  = dreal(sg(k,i,i)) 
            xloci(i) = xloci(i) + spp * spp
          enddo
          xloctot=xloctot+xloci(i)
        enddo
        if (nrot.eq.1.and.okw) then
          write (lw,51) nrot,xloctot
          write (lw,54) 
          write (lw,30) (xloci(i),i=1,n)
          xloci(1:n)=one/(xloci(1:n))
          write (lw,520)
          write (lw,30)  (xloci(i),i=1,n)
        endif
c
c.......Find the pair of MOs increasing L the most.
c
        ipq=0
        xincmax = - 1d20
        do ip=1,n-1
          do iq=ip+1,n
            ipq=ipq+1
            apq = zero
            bpq = zero
            do k = 1, ngroup
              spp  = dreal(sg(k,ip,ip)) 
              sqq  = dreal(sg(k,iq,iq))
              spq  = dreal(sg(k,ip,iq)) 
              spp2 = spp * spp
              sqq2 = sqq * sqq
              spq2 = spq * spq
              apq  = apq + spq2 - (spp2+sqq2-two*spp*sqq)/four
              bpq  = bpq + spq * (spp-sqq)
            enddo
            raiz = sqrt (apq*apq + bpq*bpq)
            xinc = apq + raiz
            rotinc(ipq)=xinc
            if (xinc.gt.xincmax) then
              xincmax=xinc
              ipair(1)=ip
              ipair(2)=iq
              apqmax=apq
              bpqmax=bpq
              raizmax=raiz
            endif
          enddo
        enddo
c
        if (abs(xincmax).lt.epsq12) moreiter=.false.
        if (moreiter) then
c
c.........Sin(4*ALPHA) and Cos(4*ALPHA)
c
          sina=+bpqmax/raizmax
          cosa=-apqmax/raizmax
c
c.........Obtain the ALPHA angle in the range 0 le ALPHA le Pi/2
c
          if (bpqmax.ge.zero) then
            if (apqmax.lt.zero) then
              gam=asin(sina)/four
            else
              gam=(pi-asin(sina))/four
            endif
          else
            if (apqmax.lt.zero) then
              gam=(two*pi+asin(sina))/four
            else
              gam=(pi-asin(sina))/four
            endif
          endif
c
c.........Print rotation angle.
c
          ip=ipair(1)
          iq=ipair(2)
          gamdeg=gam*180d0/pi
          cc=cos(gam)
          ss=sin(gam)
          cc2=cc*cc
          ss2=ss*ss
          cs=ss*cc
c
c.........Update rotation matrix
c  
          do j=1,n
            ctmp(ip,j)=+cc*c(ip,j)+ss*c(iq,j)
            ctmp(iq,j)=-ss*c(ip,j)+cc*c(iq,j)
            c(ip,j)=ctmp(ip,j)
            c(iq,j)=ctmp(iq,j)
          enddo
c
c.........Update the group overlap matrix.
c
          m=0
          do i=1,n
            if (i.ne.ip.and.i.ne.iq) then
              m=m+1
              indi(m)=i
            endif
          enddo
          do k=1,ngroup
            sgtmp(k,1:n,1:n)=sg(k,1:n,1:n)
            do ii=1,n-2
              i=indi(ii)
              tmp1=sg(k,i,ip)
              tmp2=sg(k,i,iq)           
              sgtmp(k,i,ip)=+cc*tmp1+ss*tmp2
              sgtmp(k,i,iq)=-ss*tmp1+cc*tmp2
            enddo
            do jj=1,n-2
              j=indi(jj)
              tmp1=sg(k,ip,j)
              tmp2=sg(k,iq,j)           
              sgtmp(k,ip,j)=+cc*tmp1+ss*tmp2
              sgtmp(k,iq,j)=-ss*tmp1+cc*tmp2
            enddo
            opp=sg(k,ip,ip)
            oqq=sg(k,iq,iq)
            opq=sg(k,ip,iq)
            oqp=sg(k,iq,ip)
            tmp1=+cs*opq+cs*oqp
            tmp2=-cs*opp+cs*oqq
            sgtmp(k,ip,ip)=cc2*opp+ss2*oqq+tmp1
            sgtmp(k,iq,iq)=ss2*opp+cc2*oqq-tmp1
            sgtmp(k,ip,iq)=tmp2+cc2*opq-ss2*oqp
            sgtmp(k,iq,ip)=tmp2+cc2*oqp-ss2*opq
            sg(k,1:n,1:n)=sgtmp(k,1:n,1:n)
          enddo
        endif
c
c.......Compute the property that should be zero at convergence
c
        if (nrot.eq.maxitn) warn = .true.
        if ((.not.moreiter).or.warn) then
c
c.........Test convergence. Print final results if converged.
c
          if (warn) then
            call error 
     &      ('_cmplxlo','Too many iterations (Change maxit)',warning)
            write (lw,1120) 
          else
            if (okw) then
              write (lw,'(a)') ' #'
              write (lw,300)
            else
              write (lw,400)
            endif
          endif
          if (okw) write (lw,'(a)') ' #'
c
c---------If a diagonal element is > 1.0, this element is set to 1.0,
c         the diagonal elements of other groups are set to 0.0, and
c         the non-diagonal elements of all the groups are set to 0.0.
c
          do i=1,nmo
            do k=1,ngroup
              if (dreal(sg(k,i,i)).gt.one) then
                sg(1:ngroup,1:nmo,i)=cmplx(zero,zero)
                sg(1:ngroup,i,1:nmo)=cmplx(zero,zero)
                sg(k,i,i)=cmplx(one,zero)
              endif
            enddo
          enddo
c
c.........Compute and print (DE)Localization funtion after convergence
c
          xloctot=zero
          do ip=1,n
            xloci(ip)=zero
            do k = 1, ngroup
              spp = dreal(sg(k,ip,ip))
              xloci(ip) = xloci(ip) + spp * spp
            enddo
            xloctot=xloctot+xloci(ip)
          enddo
          if (okw) then
            write (lw,510) xloctot
            write (lw,'(a)') ' #'
            write (lw,54) 
            write (lw,30)  (xloci(ip),ip=1,n)
            xloci(1:n)=one/(xloci(1:n))
            write (lw,520) 
            write (lw,30)  (xloci(ip),ip=1,n)
            write (lw,'(a)') ' #'
          endif
          if (okw) then
            write (lw,200)
            do i=1,n
              write (lw,33) i
              write (lw,30) (c(i,j),j=1,n)
            enddo
          endif
          if (okw) then
            write (lw,'(a)') ' #'
            write (lw,'(a)') ' # Cioslowski-Mixon bond orders'
          endif
c
c.........Compute and write Cioslowski-Mixon bond orders and ionicity
c         indices
c
          npair=0
          do k=2,ngroup
            npair0=npair
            do l=1,k-1
              npair=npair+1
              del=zero
              delx=zero
              do ip=1,n
                del=del+4D0*dreal(sg(k,ip,ip))*dreal(sg(l,ip,ip))
                do iq=1,n
                  delx=delx+4D0*sg(k,ip,iq)*conjg(sg(l,ip,iq))
                enddo
              enddo
              deltad(npair)=del
              delta (npair)=delx
            enddo
            if (okw) write (lw,38) k,k,k-1,(deltad(i),i=npair0+1,npair)
          enddo
          if (okw) then
            write (lw,'(a)') ' #'
            write (lw,'(a)') ' # Corrected Cioslowski-Mixon bond orders'
          endif
          npair=0
          do k=2,ngroup
            npair0=npair
            npair=npair+k-1
            if (okw) write (lw,38) k,k,k-1,(delta(i),i=npair0+1,npair)
          enddo
          if (okw) then
            write (lw,'(a)') ' #'
            write (lw,'(a)') ' # Cioslowski-Mixon ionicity indices'
          endif
          do k=2,ngroup
            do l=1,k-1
              nb=0
              do ip=1,n
                deno=dreal(sg(k,ip,ip))+dreal(sg(l,ip,ip))
                if (abs(deno).gt.zero) then
                  nb=nb+1
                  ixion(nb)=ip
                  xion(nb)=(dreal(sg(k,ip,ip))-dreal(sg(l,ip,ip)))
                  xion(nb)=100d0*abs(xion(nb)/deno)
                endif
              enddo
              if (okw) write (lw,69) k,l,(ixion(ip),xion(ip),ip=1,nb)
            enddo
          enddo
          goon=.false.
        endif
        nrot=nrot+1
      enddo

      if (allocated (sgtmp)) then
        deallocate (sgtmp,stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot deallocate sgtmp()'
      endif
      if (allocated (xion)) then
        deallocate (xion,stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot deallocate xion()'
      endif
      if (allocated (ixion)) then
        deallocate (ixion,stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot deallocate ixion()'
      endif
      if (allocated (ctmp)) then
        deallocate (ctmp,stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot deallocate ctmp()'
      endif
      if (allocated (delta)) then
        deallocate (delta,stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot deallocate delta()'
      endif
      if (allocated (deltad)) then
        deallocate (deltad,stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot deallocate deltad ()'
      endif
      if (allocated (rotinc)) then
        deallocate (rotinc,stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot deallocate rotinc()'
      endif
      if (allocated(indi)) then
        deallocate (indi,stat=ier)
        if (ier.ne.0) stop 'cmplxlo.f: Cannot deallocate indi()'
      endif
      if (okw) then
        write (lw,*) '#'
        write (lw,*) '# Done _cmplxlo'
      endif
      call timer (4,ixlo,'_cmplxlo  ',-1)
      return
c
c.....Formats
c
 20   format (' # Rotation Matrix (READ BY ROWS)',
     &        '    localized MO_i =  SUM_j C_ij delocalized MO_j')
 200  format (' # ROTATION MATRIX (READ BY ROWS)',
     &        '    localized MO_i =  SUM_j C_ij delocalized MO_j')
 30   format (6(1x,F12.6))
 31   format (2(/,2(5x,E14.8)))
 32   format (1x,'# Transformation matrix')
 300  format (' #',/,' #',80('-'),/,' # CONVERGED LOCALIZATION',/,' #')
 400  format (1x,'# CONVERGED LOCALIZATION')
 51   format (1x,'# ',75('-'),/,' # L value before rotation ' ,I8,
     &        ' = ',E16.10)
 71   format (' # L value after  rotation ' ,I4,' = ',E16.10)
 53   format (' # L value after completing ',I4,' rotations = ',E16.10)
 510  format (' # L value             = ',E16.10)
 54   format (' # L_i (i=1..N)')
 520  format (' # Fragments expanded by MOs, L_i^(-1) (i=1..N)')
 56   format (1x,'# ROTATION ANGLE = ',F12.6,
     &  ' degree for pair ',2I3,/,' #')
 23   format (/1x,'# Rotated orbitals natural populations')
 233  format (/1x,'# INPUT orbitals natural populations')
 1120 format (1x,'#',/,1x,'# ',80('+'),/,
     &  ' # UNCONVERGED RESULTS',/,' #')
 66   format (' # Fragment occupation numbers of the localized MOs')
 67   format (' # FRAGMENT ',I2)
 38   format (' # delta(',I2,' 1) to delta(',2I2,') = ',12(1x,F8.5))
 69   format (' # Ionicity indices for pair ',2I2,/,
     &   ' # The associated MO is indicated previous to the index',/,
     &         6(2x,I5,1x,F5.1,'%'))
 610  format (' #',/,
     &   ' # Final potential increase in L for the rotation of all',
     &   ' pairs of MOs')
 612  format (' # ',a,a)
 611  format (10(1x,E12.6))
 33   format (' # Localized MO number ',I3)
      end
