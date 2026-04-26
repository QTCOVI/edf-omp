c
c-----------------------------------------------------------------------
c
      subroutine canlmo (sg,c,nmo,ngroup,critics,lw,largwr,okw)
c
c.....canlmo - Performs the localization of a set of molecular orbitals 
c     into ngroup by maximizing 
c                            
c     L = sum [i=1,NMO] nu_i^2 L_i
c                                   
c     L_i = SUM [klmn] A_klmn
c
c     A_klmn = U^*(ik) U(il) U^*(im) U(in) N(klmn) T(klmn)
c
c     N(klmn) = (n_k n_l n_m n_n)^(1/2)
c
c     T(klmn)=SUM_A <NMO_k|NMO_l>_A <NMO_m|NMO_n>_A
c
c.....where
c
c     NMO ....... Number of natural molecular orbitals (NMO)
c     U(rs) ....  The elements of an unitary matrix
c     A ......... Index that runs over all the fragments
c     NMO_r ..... The r^th molecular orbital (MO)
c     n_r ....... Natural population of r^th MO.
c     nu_i ...... Natural population of the i^th isopycnic MO
c                            
c-----PARAMETERS--------------------------------------------------------
c
c.....sg()                           INPUT/OUTPUT
c
c     Group overlap matrix between input MO's when entering the 
c     routine or between output MO's at the end of the routine.
c
c.....c()                            OUTPUT 
c
c     Rotation matrix giving the transformed MO's from the input MO's. 
c
c.....nmo                            INPUT
c
c     Number of MO's.
c
c.....ngroup                         INPUT
c
c     Number of fragments 
c
c     Maximum number of groups
c
c.....lw                             INPUT
c
c     Output logical unit
c
c.....largwr                         INPUT
c
c     .true.  ====>  Large output
c     .false. ====>  Short output
c
c.....okw                            INPUT
c
c     .false. ====> The routine writes nothing.
c     .true.  ====> Output according to the value of the LARGWR value
c-----------------------------------------------------------------------
c                         
      include 'implicit.inc'
      include 'constants.inc'
      include 'error.inc'
c
      parameter (maxit=200)
      parameter (xmoloc=0.05d0)

      real(kind=8)   c(nmo,nmo)
      real(kind=8)   pi,pi4
      real(kind=8)   sg(ngroup,nmo,nmo)
      real(kind=8)   apqsum(ngroup),bpqsum(ngroup)

      real(kind=8),allocatable,dimension (:)     :: xloci,xion,degloc
      real(kind=8),allocatable,dimension (:,:)   :: over
      real(kind=8),allocatable,dimension (:,:,:) :: sgtmp
      real(kind=8),allocatable,dimension (:,:)   :: ctmp
      real(kind=8),allocatable,dimension (:)     :: delta,deltad,rotinc
      integer, allocatable,dimension (:)         :: indi,ixion
      integer, allocatable,dimension (:,:)       :: blg
      integer, allocatable,dimension (:)         :: nbatom,ncom

      logical   warn,largwr,okw,goon,condition
      logical   moreiter
      integer   ipair(2)

      real(kind=8) crits(29)
      data crits /95D0,10D0,20D0,30D0,40D0,50D0,60D0,70D0,80D0,
     &  90D0,91D0,92D0,93D0,94D0,95D0,96D0,97D0,
     &  98D0,99D0,99.1D0,99.2D0,99.3D0,99.4D0,99.5D0, 
     &  99.6D0,99.7D0,99.8D0,99.9D0,99.99D0/
c
      call timer (2,ipid,'_canlmo   ',-1)
c
c.....Allocate arrays
c
      n = nmo
      allocate (over(nmo,nmo),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate over()'
      allocate (sgtmp(ngroup,nmo,nmo),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate sgtmp()'
      allocate (xloci(nmo),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate xloci()'
      allocate (xion(nmo),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate xion()'
      allocate (degloc(nmo),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate degloc()'
      allocate (ixion(nmo),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate ixion()'
      allocate (blg(ngroup,nmo),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate blg()'
      allocate (nbatom(ngroup),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate nbatom()'
      allocate (ncom(nmo),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate ncom()'
      allocate (ctmp(nmo,nmo),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate ctmp()'
      allocate (delta(ngroup*(ngroup-1)/2),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate delta()'
      allocate (deltad(ngroup*(ngroup-1)/2),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate deltad()'
      allocate (rotinc(nmo*(nmo-1)/2),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate rotinc()'
      allocate (indi(nmo),stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot allocate indi()'

      write (lw,500)
      pi4  = atan(one)
      four = two+two
      pi   = four*pi4
      warn =.false.
      over(1:n,1:n)=zero
      do k=1,ngroup
        over(1:n,1:n)=over(1:n,1:n)+sg(k,1:n,1:n)
      enddo
c
c.....Initial rotation matrix is the unit matrix.
c
      c=zero
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
        xloctot = zero
        do ip=1,n
          xloci(ip)=ddot(ngroup,sg(:,ip,ip),1,sg(:,ip,ip),1)
        enddo
        xloctot=sum(xloci)
c
c.......Find the pair of MOs increasing L the most.
c
        ipq=0
        xincmax = - 1q20
        do ip=1,n-1
          do iq=ip+1,n
            ipq=ipq+1
            do k = 1, ngroup
              spp  = sg(k,ip,ip)
              sqq  = sg(k,iq,iq)
              spq  = sg(k,ip,iq)
              sdif = spp-sqq
              apqsum(k) = spq * spq - sdif*sdif/four
              bpqsum(k) = spq * sdif
            enddo
            apq=sum(apqsum)
            bpq=sum(bpqsum)
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
            gam=asin(sina)/four
            if (apqmax.ge.zero) gam=pi4-gam
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
          cc     = cos(gam)
          ss     = sin(gam)
          cc2    = cc  * cc
          ss2    = ss  * ss
          cs     = cc  * ss
          cc2ss2 = cc2 + ss2
          cpp    =+cc
          cpq    =+ss
          cqp    =-ss
          cqq    =+cc
c
c.........Update rotation matrix
c  
          do j=1,n
            ctmp(ip,j)=+cpp*c(ip,j)+cpq*c(iq,j)
            ctmp(iq,j)=+cqp*c(ip,j)+cqq*c(iq,j)
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
            call error ('canlmo.f','Too many iterations',warning)
            write (lw,1120) 
          else
            write (lw,400)
          endif
c
c---------If a diagonal element is > 1.0, this element is set to 1.0,
c         the diagonal elements of other groups are set to 0.0, and
c         the non-diagonal elements of all the groups are set to 0.0.
c
          do i=1,nmo
            do k=1,ngroup
              if (sg(k,i,i).gt.one) then
                sg(1:ngroup,1:nmo,i)=zero
                sg(1:ngroup,i,1:nmo)=zero
                sg(k,i,i)=one
              endif
            enddo
          enddo
c
c.........Compute and print (DE)Localization funtion after convergence
c
          xloctot = zero
          do ip=1,n
            xloci(ip)=ddot(ngroup,sg(:,ip,ip),1,sg(:,ip,ip),1)
          enddo
          xloctot=sum(xloci)
          if (okw) then
            write (lw,510) xloctot
            write (lw,'(a)') ' #'
            write (lw,54) 
            write (lw,30)  (xloci(ip),ip=1,n)
            write (lw,520) 
            write (lw,30)  (one/xloci(ip),ip=1,n)
            write (lw,'(a)') ' #'
          endif
          over(1:n,1:n)=zero
          do k=1,ngroup
            over(1:n,1:n)=over(1:n,1:n)+sg(k,1:n,1:n)
          enddo
          goon=.false.
        endif
        nrot=nrot+1
      enddo
c
      deallocate (over,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate over()'
      deallocate (sgtmp,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate sgtmp()'
      deallocate (xloci,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate xloci()'
      deallocate (xion,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate xion()'
      deallocate (ixion,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate ixion()'
      deallocate (blg,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate blg()'
      deallocate (nbatom,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate nbatom()'
      deallocate (ncom,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate ncom()'
      deallocate (degloc,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate degloc()'
      deallocate (ctmp,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate ctmp()'
      deallocate (delta,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate delta()'
      deallocate (deltad,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate deltad ()'
      deallocate (rotinc,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate rotinc()'
      deallocate (indi,stat=ier)
      if (ier.ne.0) stop 'canlmo.f: Cannot deallocate indi()'
      call timer (4,ipid,'_canlmo   ',-1)
      return
c
c.....Formats
c
 20   format (' # Rotation Matrix (READ BY ROWS)',
     &        '    locMO_i =  SUM_j C_ij natMO_j')
 200  format (' # ROTATION MATRIX (READ BY ROWS)',
     &        '    locMO_i =  SUM_j C_ij natMO_j')
 30   format (10(1x,F12.6))
 301  format (5x,'Trace =  ',F12.6)
 31   format (2(/,2(5x,E14.7)))
 32   format (1x,'# Transformation matrix')
 300  format (' #',/,' #',80('-'),/,' # CONVERGED LOCALIZATION',/,' #')
 500  format (/1x,'# Performing Localization ------------->')
 400  format (1x,'# CONVERGED LOCALIZATION')
 51   format (' #',80('-'),/,' #',/,
     &        ' # L value before rotation ' ,I4,' = ',E16.9)
 71   format (' # L value after  rotation ' ,I4,' = ',E16.9)
 53   format (' # L value after completing ',I4,' rotations = ',E16.9)
 510  format (' # L value             = ',E16.9)
 54   format (' # L_i (i=1..N)')
 520  format (' # Fragments expanded by MOs, L_i^(-1) (i=1..N)')
 56   format (1x,'# ROTATION ANGLE = ',F12.6,
     &  ' degree for the pair of atoms ',2I3,/,' #')
 23   format (/1x,'# Rotated orbitals natural populations')
 1120 format (1x,'#',/,1x,'# ',80('+'),/,
     &  ' # UNCONVERGED RESULTS',/,' #')
 66   format (' # Fragment occupation numbers of the localized MOs')
 67   format (' # FRAGMENT ',I2)
 38   format (' # delta(',I4,' 1) to delta(',2I4,') = ',/,
     &        10(1x,F8.5))
 69   format (' # Pair of atoms ',2I4,/,
     &        5(1x,I3,1x,F5.1,'%',' (',F6.2'%)'))
 610  format (' #',/,
     &   ' # Final potential increase in L for the rotation of all',
     &   ' pairs of MOs')
 612  format (' # ',a,a)
 611  format (10(1x,E13.6))
 33   format (' # Localized MO number ',I3)
 1011 format (' #',/,
     &' # Cioslowski-Mixon ionicity indices: ',/,
     &' # ---------------------------------- ',/,
     &' # Only MOs localized more than ',F5.2,'% between the two',1x,
     &'atoms of the pair are considered',/,
     &' # The associated MO is indicated previous to the index',/,
     &' # and the degree of localization between the two atoms in',1x,
     &'parenthesis')
*1022  format (' # Center ',I4,/,1000(1x,20I4,/))
 1022  format (4x,I4,' --> ',25I4)
 1023  format (1x,'# Total MOs = ',I4)
 1024  format (/1x,'# MOs partially localized on each fragment',
     &        '  (S_ii^A > ',1PE15.7,' )')
 1026  format (//' # MOs localized more than',F7.2,'% on each fragment')
 1027  format (' # Fragment ',I4,' --> ',20(1x,20I4,/))
 1025  format (/1x,'# Diagonal overlaps in each fragment')
 1028  format (/' # MOs localized more than',F7.2,
     & '% on a pair of fragments')
 1031  format (/' # MOs localized more than',F7.2,
     & '% on a trio of fragments')
 1037  format (/' # MOs localized more than',F7.2,
     & '% on a quartet of fragments')
 1029  format (' # Fragments ',10I4)
 1030  format (3x,'MO ',I4,5x,F5.2,'% + ',F5.2,'%')
 1033  format (3x,'MO ',I4,5x,F5.2,'% + ',F5.2,'%+ ',F5.2,'%')
 1036  format (3x,
     &  'MO ',I4,5x,F5.2,'% + ',F5.2,'%+ ',F5.2,'%+ ',F5.2,'%')
 121   format (I4,a,I4,2x,'(',I4,' MOs) --> ',25I4)
 122   format (' # MOs partially localized in A & B (Intersection)') 
 124   format (' # MOs partially localized in A .or. B (Union)')
      end
