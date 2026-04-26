c
c.....Obtain spinless multiple-fragment electron population covariances.
c
      subroutine ndelta (pt,pa,pb,ea,eb,ngr,npa,npb,lw)
c
      include   'implicit.inc'
      include   'param.inc'
      include   'constants.inc'
c
      integer    ea(npa,ngr),eb(npb,ngr)
      real(kind=8)    pa(ngr),pb(ngr)
      real(kind=8)    pt(npa*npb)
      integer    ig(20)
c
      if (ngr.gt.12) then
        write (0,*) ' # ndelta.f: Too many groups, Max value = ',ngr
        stop 
      endif
      if (ngr.gt.1) write (lw,112) 
      do i1=1,ngr
      ig(1)=i1
      if (ngr.gt.1) then 
      do i2=1,i1-1
      ig(2)=i2
      call covnce (pt,pa,pb,ea,eb,ngr,npa,npb,ig,2,lw)
      if (ngr.gt.2) then
      do i3=1,i2-1
      ig(3)=i3
      call covnce (pt,pa,pb,ea,eb,ngr,npa,npb,ig,3,lw)
      if (ngr.gt.3) then
      do i4=1,i3-1
      ig(4)=i4
      call covnce (pt,pa,pb,ea,eb,ngr,npa,npb,ig,4,lw)
      if (ngr.gt.4) then
      do i5=1,i4-1
      ig(5)=i5
      call covnce (pt,pa,pb,ea,eb,ngr,npa,npb,ig,5,lw)
      if (ngr.gt.5) then
      do i6=1,i5-1
      ig(6)=i6
      call covnce (pt,pa,pb,ea,eb,ngr,npa,npb,ig,6,lw)
      if (ngr.gt.6) then
      do i7=1,i6-1
      ig(7)=i7
      call covnce (pt,pa,pb,ea,eb,ngr,npa,npb,ig,7,lw)
      if (ngr.gt.7) then
      do i8=1,i7-1
      ig(8)=i8
      call covnce (pt,pa,pb,ea,eb,ngr,npa,npb,ig,8,lw)
      if (ngr.gt.8) then
      do i9=1,i8-1
      ig(9)=i9
      call covnce (pt,pa,pb,ea,eb,ngr,npa,npb,ig,9,lw)
      if (ngr.gt.9) then
      do i10=1,i9-1
      ig(10)=i10
      call covnce (pt,pa,pb,ea,eb,ngr,npa,npb,ig,10,lw)
      if (ngr.gt.10) then
      do i11=1,i10-1
      ig(11)=i11
      call covnce (pt,pa,pb,ea,eb,ngr,npa,npb,ig,11,lw)
      if (ngr.gt.11) then
      do i12=1,i11-1
      ig(12)=i12
      call covnce (pt,pa,pb,ea,eb,ngr,npa,npb,ig,12,lw)
      enddo
      endif
      enddo
      endif
      enddo
      endif
      enddo
      endif
      enddo
      endif
      enddo
      endif
      enddo
      endif
      enddo
      endif
      enddo
      endif
      enddo
      endif
      enddo 
      endif
      enddo
      return
 112  format (1x,'#',/,1x,'# Multiple-Fragment Covariances:',
     & ' COV (i,j..k) = < PROD (x=1..k) (n_x - <n_x>) >',/,' #')
      end
c
c-----------------------------------------------------------------------
c
      subroutine covnce (pt,pa,pb,ea,eb,ngr,npa,npb,ig,nig,lw)
      include 'implicit.inc'
      include 'param.inc'
      include 'constants.inc'
      integer  ig(20)
      integer  ea(npa,ngr),eb(npb,ngr)
      real(kind=8)  pa(ngr),pb(ngr),pt(npa*npb)
      real(kind=8)  delt
      integer  ij
c
      delt=zero
      ij=0
      do i=1,npa
        do j=1,npb
          ij=ij+1
          prod=one
          do ix=1,nig
            igx=ig(ix)
            prod=prod*(ea(i,igx)-pa(igx)+eb(j,igx)-pb(igx))
          enddo
          delt=delt+prod*pt(ij)
        enddo
      enddo
      write (lw,10) delt,(ig(k),k=1,nig)
      return
 10   format (1x,'# COVARIANCE = ',1x,F16.10,4x,'GROUPS --> ',20(I3))
      end
