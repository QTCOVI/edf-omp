      subroutine computerdm (det,Ndet,coef,epsdet,nmo,Nint,c1ea,c1eb)
      implicit none
      integer, intent(in)         :: Ndet, Nint, nmo
      integer(kind=8), intent(in) :: det(Nint,2,Ndet)
      real(kind=8), intent(in)    :: coef(Ndet)
      real(kind=8), intent(in)    :: epsdet
      integer(kind=8)             :: det1(Nint,2),det2(Nint,2)
      real(kind=8), intent(out)   :: c1ea(nmo,nmo),c1eb(nmo,nmo)
      real(kind=8)                :: cd1,cd2,cd,cdt
      integer(kind=8)             :: buffer1,buffer2
      real(kind=8)                :: phase,c,ct
      integer                     :: i1,i2,k,l,n,il,ll,lm
      integer                     :: in1,in2,im,jm,jn
      integer                     :: inm
      integer                     :: in,ji,jl,km,kn,mn
      integer                     :: ishift1,ishift2
      integer                     :: deg,ispin
      integer                     :: exc(0:2,2,2)
      integer                     :: nexcitations

      real(kind=8)                :: indet,shundred,xhundred
      integer                     :: iper


!     exc(i,j,k)
!     k = 1,2   ==> alpha,beta
!     i = 1,2   ==> 1st excitation, 2nd excitation
!     exc(1,2,k) = 1st created particle in det2 with spin=k
!     exc(1,1,k) = 1st created hole     in det1 with spin=k
!     exc(2,2,k) = 2nd created particle in det2 with spin=k
!     exc(2,1,k) = 2nd created hole     in det1 with spin=k
!     exc(0,j,k) = Number of created particles in det2 (j=2) or
!                  number of created holes     in det1 (j=1) with spin=k
      c1ea = 0d0
      c1eb = 0d0

      iper     = 0
      shundred = (dble(Ndet)*(dble(Ndet)-1))/2d0/100d0
      xhundred = shundred
      indet    = 0d0

      do k=1,Ndet
        cd1 = coef(k)
        cd  = cd1 * cd1
        cdt = cd + cd
        if (abs(cd) <= abs(epsdet)) goto 1000
!
!       Diagonal terms
!
        det1(:,:) = det(:,:,k)
!
!       Alpha RDMs
!
        ishift1 = 0
        do in1=1,Nint
          buffer1 = det1(in1,1)
          do while (buffer1 /= 0_8)
            jn = trailz(buffer1) + ishift1 + 1
            c1ea(jn,jn) = c1ea(jn,jn) + cd
            buffer1 = iand(buffer1,buffer1-1_8)
          end do
          ishift1 = ishift1 + 63
        end do
!
!       Beta RDMs
!
        ishift1 = 0
        do in1=1,Nint
          buffer1 = det1(in1,2)
          do while (buffer1 /= 0_8)
            jn = trailz(buffer1) + ishift1 + 1
            c1eb(jn,jn) = c1eb(jn,jn) + cd
            buffer1 = iand(buffer1,buffer1-1_8)
          end do
          ishift1 = ishift1 + 63
        end do
 1000   continue
!
!       Non-Diagonal terms
!
        do l=1,k-1

          indet=indet+1d0
          if (indet.ge.xhundred.and.Ndet.gt.1000) then
            xhundred=xhundred+shundred
            iper=iper+1
            write (0,1212) iper
          endif

          cd2 = coef(l)
          cd  = cd1 * cd2
          cdt = cd + cd
          if (abs(cd) <= abs(epsdet)) cycle
          det2(:,:) = det(:,:,l)
          deg = nexcitations (det1,det2,Nint)
          select case (deg)
          case (0)
!
!           This cannot happen since we are in a Non-Diagonal term
!
            cycle
          case (1)
!
!           Single excitation
!
            call getexcitation (det1,det2,exc,deg,phase,Nint)
            c  = phase * cd
            ct = phase * cdt
!
!           The excitation corresponds to a ALPHA elecron ...
!
            if (exc(0,1,1) == 1) then
              ji = exc(1,1,1)
              km = exc(1,2,1)
              c1ea(km,ji) = c1ea(km,ji) + c
              c1ea(ji,km) = c1ea(ji,km) + c
            else
              ji = exc(1,1,2)
              km = exc(1,2,2)
              c1eb(km,ji) = c1eb(km,ji) + c
              c1eb(ji,km) = c1eb(ji,km) + c
            endif
          case(2:)
!
!           Double excitations do not contribute to the RDMs
!
            cycle
          end select
        end do
      end do
 1212 format (1x,'# Computing 1-RDM :',i6,'% done')
      end



      subroutine getexcitation (det1,det2,exc,degree,phase,Nint)
      implicit none
      integer, intent(in)         :: Nint
      integer(kind=8), intent(in) :: det1(Nint,2), det2(Nint,2)
      integer, intent(out)        :: exc(0:2,2,2)
      integer, intent(out)        :: degree
      real(kind=8), intent(out)   :: phase
      integer                     :: nexcitations

      degree = nexcitations (det1,det2,Nint)
      select case (degree)
      case (3:)
        degree = -1
        return
      case (2)
        call doublexcitation (det1,det2,exc,phase,Nint)
        return
      case (1)
        call singlexcitation (det1,det2,exc,phase,Nint)
        return
      case(0)
        return
      end select
      end



      integer function nexcitations (det1,det2,Nint)
      implicit none
      integer(kind=8), intent(in) :: det1(Nint,2), det2(Nint,2)
      integer, intent(in) :: Nint
      integer(kind=8)     :: d1a,d2a,d1b,d2b
      integer :: l

      d1a    = det1(1,1)
      d2a    = det2(1,1)
      d1b    = det1(1,2)
      d2b    = det2(1,2)
      nexcitations = popcnt(xor(d1a,d2a)) + popcnt(xor(d1b,d2b))
      do l=2,Nint
        d1a  = det1(l,1)
        d2a  = det2(l,1)
        d1b  = det1(l,2)
        d2b  = det2(l,2)
        nexcitations = nexcitations + popcnt(xor(d1a,d2a)) &
                                    + popcnt(xor(d1b,d2b))
      end do
      nexcitations = ishft(nexcitations,-1)
      end




      subroutine singlexcitation (det1,det2,exc,phase,Nint)
      implicit none
      integer, intent(in)           :: Nint
      integer(kind=8), intent(in)   :: det1(Nint,2)
      integer(kind=8), intent(in)   :: det2(Nint,2)
      integer, intent(out)          :: exc(0:2,2,2)
      double precision, intent(out) :: phase
      integer                       :: tz,l,ispin,ishift,nperm
      integer                       :: i,j,k,m,n,high,low
      integer(kind=8)               :: hole, particle, tmp
      real(kind=8), parameter       :: phase_dble(0:1) = (/1.d0,-1.d0/)

      exc(0,1,1) = 0
      exc(0,2,1) = 0
      exc(0,1,2) = 0
      exc(0,2,2) = 0
      do ispin = 1,2
        ishift = 0
        do l=1,Nint
          if (det1(l,ispin) == det2(l,ispin)) cycle
          tmp = xor( det1(l,ispin), det2(l,ispin) )
          particle = iand(tmp, det2(l,ispin))
          hole = iand(tmp, det1(l,ispin))
          if (particle /= 0_8) then
            tz = trailz(particle)
            exc(0,2,ispin) = 1
            exc(1,2,ispin) = tz+ishift+1
          end if
          if (hole /= 0_8) then
            tz = trailz(hole)
            exc(0,1,ispin) = 1
            exc(1,1,ispin) = tz+ishift+1
          end if
          if ( iand(exc(0,1,ispin),exc(0,2,ispin)) == 1 ) then
            low = min(exc(1,1,ispin),exc(1,2,ispin))
            high = max(exc(1,1,ispin),exc(1,2,ispin))
            j = ishft(low-1,-6)+1
            n = iand(low,63)
            k = ishft(high-1,-6)+1
            m = iand(high,63)
            if (j==k) then
              nperm = popcnt(iand(det1(j,ispin), &
              iand( ibset(0_8,m-1)-1_8, ibclr(-1_8,n)+1_8 ) ))
            else
              nperm = popcnt(iand(det1(k,ispin), ibset(0_8,m-1)-1_8)) &
                    + popcnt(iand(det1(j,ispin), ibclr(-1_8,n) +1_8))
              do i=j+1,k-1
                nperm = nperm + popcnt(det1(i,ispin))
              end do
            end if
            phase = phase_dble(iand(nperm,1))
            return
          end if
          ishift = ishift + 63
        end do
      end do
      end




      subroutine doublexcitation (det1,det2,exc,phase,Nint)
      implicit none
      integer, intent(in)           :: Nint
      integer(kind=8), intent(in)   :: det1(Nint,2), det2(Nint,2)
      integer, intent(out)          :: exc(0:2,2,2)
      double precision, intent(out) :: phase
      integer                       :: l,ispin,idx_hole
      integer                       :: idx_particle,ishift
      integer                       :: i,j,k,m,n,high,low,a,b,c,d
      integer                       :: nperm,tz,nexc
      integer(kind=8)               :: hole,particle,tmp
      real(kind=8), parameter       :: phase_dble(0:1) = (/1.d0,-1.d0/)

      exc(0,1,1) = 0
      exc(0,2,1) = 0
      exc(0,1,2) = 0
      exc(0,2,2) = 0
      nexc       = 0
      nperm      = 0
      do ispin = 1,2
        idx_particle = 0
        idx_hole     = 0
        ishift       = 0
        do l=1,Nint
          if (det1(l,ispin) == det2(l,ispin)) then
            cycle
          end if
          tmp      = xor( det1(l,ispin), det2(l,ispin) )
          particle = iand(tmp, det2(l,ispin))
          hole     = iand(tmp, det1(l,ispin))
          do while (particle /= 0_8)
            tz             = trailz(particle)
            nexc           = nexc+1
            idx_particle   = idx_particle + 1
            exc(0,2,ispin) = exc(0,2,ispin) + 1
            exc(idx_particle,2,ispin) = tz+ishift+1
            particle       = iand(particle,particle-1_8)
          end do
          do while (hole /= 0_8)
            tz                    = trailz(hole)
            nexc                  = nexc+1
            idx_hole              = idx_hole + 1
            exc(0,1,ispin)        = exc(0,1,ispin) + 1
            exc(idx_hole,1,ispin) = tz+ishift+1
            hole                  = iand(hole,hole-1_8)
          end do
          if (nexc == 4) exit
          ishift = ishift + 63
        end do
        do i=1,exc(0,1,ispin)
          low  = min(exc(i,1,ispin),exc(i,2,ispin))
          high = max(exc(i,1,ispin),exc(i,2,ispin))
          j    = ishft(low-1,-6)+1
          n    = iand(low,63)
          k    = ishft(high-1,-6)+1
          m    = iand(high,63)
          if (j==k) then
            nperm = nperm + popcnt(iand(det1(j,ispin), &
                    iand( ibset(0_8,m-1)-1_8, ibclr(-1_8,n)+1_8 ) ))
          else
            nperm = nperm + popcnt(iand(det1(k,ispin), &
                    ibset(0_8,m-1)-1_8))               &
            + popcnt(iand(det1(j,ispin),               &
            ibclr(-1_8,n) +1_8))
            do l=j+1,k-1
              nperm = nperm + popcnt(det1(l,ispin))
            end do
          end if
        end do
        if (exc(0,1,ispin) == 2) then
          a = min(exc(1,1,ispin), exc(1,2,ispin))
          b = max(exc(1,1,ispin), exc(1,2,ispin))
          c = min(exc(2,1,ispin), exc(2,2,ispin))
          d = max(exc(2,1,ispin), exc(2,2,ispin))
          if (c>a .and. c<b .and. d>b) nperm = nperm + 1
          exit
        end if
      end do
      phase = phase_dble(iand(nperm,1))
      end




      subroutine get_particles (det1,det2,particles,Nint,nmo)
      implicit none
      integer(kind=8), intent(in)  :: det1(Nint,2), det2(Nint,2)
      integer(kind=8), intent(out) :: particles(nmo,2)
      integer(kind=8)              :: p,position
      integer, intent(in)          :: Nint,nmo
      integer                      :: l,k,ispin

      do ispin = 1,2
        k = 1
        do l = 1, Nint
          p = and ( xor(det1(l,1),det2(l,1)), det2(l,1))
          do while ( p /= 0)
            position = trailz (p)
            particles (k,ispin) = 1 + 64 * (l-1) + position
            p = ibclr(p,position)
            k = k + 1
          enddo
        enddo
      enddo
      end





      subroutine get_holes (det1,det2,holes,Nint,nmo)
      implicit none
      integer(kind=8), intent(in)  :: det1(Nint,2), det2(Nint,2)
      integer(kind=8), intent(out) :: holes(nmo,2)
      integer(kind=8)              :: h,position
      integer, intent(in)          :: Nint,nmo
      integer                      :: l,k,ispin

      do ispin = 1,2
        k = 1
        do l = 1, Nint
          h = and ( xor(det1(l,1),det2(l,1)), det1(l,1))
          do while ( h /= 0)
            position = trailz (h)
            holes (k,ispin) = 1 + 64 * (l-1) + position
            h = ibclr(h,position)
            k = k + 1
          enddo
        enddo
      enddo
      end
