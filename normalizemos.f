c
c-----------------------------------------------------------------------
c
      SUBROUTINE normalizemos (cc,nmo,nprims)
      USE space_for_wfnbasis
      implicit none
      integer(kind=4) nmo,nprims,ma,mb,i
      real   (kind=8) cc(nmo,nprims)
      real   (kind=8) solap,prodc
c
c     computes overlaps between primitive cartesian GTOs
c
      if (.not.sprimok) then
        call gtogto ()
        sprimok=.true.
      endif
      do i=1,nmo
        solap=0d0
        do ma=1,nprims
          do mb=1,ma
            prodc=cc(i,ma)*cc(i,mb)
            if (ma.ne.mb) prodc=prodc+prodc
            solap=solap+prodc*sprim(ma,mb)
          enddo
        enddo
        cc(i,1:nprims)=cc(i,1:nprims)/sqrt(solap)
      enddo
      return
      end
