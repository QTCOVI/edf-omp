c
c-----------------------------------------------------------------------
c
      subroutine cutwfn (epswfn,ifilc,ndets,nelact)
c
c.....Cut WFN elliminating Slater dets with a coefficient < epswfn.
c.......................................................................
c
      USE       space_for_cidet
      include  'implicit.inc'
      include  'constants.inc'
      include  'param.inc'
      integer   ifilc,ndets,ndetsnew,nelact
c
c.....Read the coeffs for the first time just to compute the new value 
c     of the normalization constant.
c
      xnormnew=zero
      do i=1,ndets
        read (ifilc,rec=i) (cidet(k),k=0,nelact)
        cd1=cidet(0)
        if (abs(cd1).ge.epswfn) then
           xnormnew=xnormnew+cd1*cd1
        endif
      enddo
      ynormnew=one/sqrt(xnormnew)
c
c.....Read the coeffs for the second time and write the modified coeffs
c     that satisfy coefficient < epswfn.

      ndetsnew=0
      do i=1,ndets
        read (ifilc,rec=i) (cidet(k),k=0,nelact)
        cd1=cidet(0)
        if (abs(cd1).ge.epswfn) then
           ndetsnew=ndetsnew+1
c
c..........Modify coeff by the normalization constant.
c
           cidet(0)=cd1*ynormnew
           write (ifilc,rec=ndetsnew) (cidet(k),k=0,nelact)
        endif
      enddo
c
c.....New number of detS in terms of which the WFN is expressed.
c
      ndets=ndetsnew
      return
      end
