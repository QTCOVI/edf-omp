      subroutine uhfnorm (iren,ncent,nmo,nsigma,tol,lw,lerr,isigma,aom)
      implicit real(kind=8) (a-h,o-z)
      integer (kind=4) ncent,nsigma,nmo,i,j,m,ii,iren,lw,lerr
      integer (kind=4) isigma(nsigma) 
      real    (kind=8) aom(ncent,nmo,nmo),sumaom,valueaom,diffaom
*     character (len=*) line
      character (len=50) word
      logical setint,setword,ok
c
c.....Test consistency of AOM elements
c
      do m=1,nsigma
        do j=1,m
          sumaom=0d0
          do ii=1,ncent
            sumaom=sumaom+aom(ii,isigma(j),isigma(m))
          enddo
          valueaom=1d0
          if (j.ne.m) valueaom=0d0
          diffaom=abs(sumaom-valueaom)
          if (diffaom.gt.tol) then
            write (lerr,544) isigma(m),isigma(j),sumaom
          endif
        enddo
      enddo
c     
c.....Modify the AOM elements of center IREN such that 
c     SUM (i=1,ncent) aom(i,m,j) will be exactly 0 or 1.
c
      write (lw,335) iren
      do j=1,nsigma
        do m=1,nsigma
          aom2=aom(iren,isigma(m),isigma(j))
          aom1=0d0
          do i=1,ncent
            if (i.ne.iren) aom1=aom1+aom(i,isigma(m),isigma(j))
          enddo
          aomtot=aom1+aom2
          if (m.eq.j) then
            if (aomtot.lt.0d0) then
              aom(iren,isigma(m),isigma(j))=aom2-(1d0-abs(aomtot))
            endif
            if (aomtot.ge.0d0) then
              aom(iren,isigma(m),isigma(j))=aom2+(1d0-abs(aomtot))
            endif
          else
            aom(iren,isigma(m),isigma(j))=-aom1
          endif
        enddo
      enddo
      return
 544  format (' # !! AOM is not Ok !! : SUM_n AOM(n,',2I3,') = ',E15.8)
 335  format (1x,'# AOM is renormalized, the reference atom is ',I5)
      end
