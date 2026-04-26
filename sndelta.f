c
c-----------------------------------------------------------------------
c
c
c.....Obtain spinless multiple-fragment electron population covariances.
c
      subroutine sndelta (probs,p1a,p1b,ela,elb,ngr,npa,npb,npx,lw)
c
      include   'implicit.inc'
      include   'param.inc'
      include   'constants.inc'
c
      integer    ela(npa,ngr),elb(npb,ngr)
      real(kind=8)    p1a(ngr),p1b(ngr),probs(npx,2)
      integer    ig(20)
c
      call timer (2,isndel,'_sndelta  ',-1)
      if (ngr.gt.12) then
        write (0,*) ' # sndelta.f: Too many groups, Max value = ',ngr
        stop 
      endif
      if (ngr.gt.1) write (lw,112) 
      do i1=1,ngr
      ig(1)=i1
      if (ngr.gt.1) then 
      do i2=1,i1-1
      ig(2)=i2
      call scovnce (npa,npb,ngr,npx,ela,elb,p1a,p1b,probs,ig,2,lw)
      if (ngr.gt.2) then
      do i3=1,i2-1
      ig(3)=i3
      call scovnce (npa,npb,ngr,npx,ela,elb,p1a,p1b,probs,ig,3,lw)
      if (ngr.gt.3) then
      do i4=1,i3-1
      ig(4)=i4
      call scovnce (npa,npb,ngr,npx,ela,elb,p1a,p1b,probs,ig,4,lw)
      if (ngr.gt.4) then
      do i5=1,i4-1
      ig(5)=i5
      call scovnce (npa,npb,ngr,npx,ela,elb,p1a,p1b,probs,ig,5,lw)
      if (ngr.gt.5) then
      do i6=1,i5-1
      ig(6)=i6
      call scovnce (npa,npb,ngr,npx,ela,elb,p1a,p1b,probs,ig,6,lw)
      if (ngr.gt.6) then
      do i7=1,i6-1
      ig(7)=i7
      call scovnce (npa,npb,ngr,npx,ela,elb,p1a,p1b,probs,ig,7,lw)
      if (ngr.gt.7) then
      do i8=1,i7-1
      ig(8)=i8
      call scovnce (npa,npb,ngr,npx,ela,elb,p1a,p1b,probs,ig,8,lw)
      if (ngr.gt.8) then
      do i9=1,i8-1
      ig(9)=i9
      call scovnce (npa,npb,ngr,npx,ela,elb,p1a,p1b,probs,ig,9,lw)
      if (ngr.gt.9) then
      do i10=1,i9-1
      ig(10)=i10
      call scovnce (npa,npb,ngr,npx,ela,elb,p1a,p1b,probs,ig,10,lw)
      if (ngr.gt.10) then
      do i11=1,i10-1
      ig(11)=i11
      call scovnce (npa,npb,ngr,npx,ela,elb,p1a,p1b,probs,ig,11,lw)
      if (ngr.gt.11) then
      do i12=1,i11-1
      ig(12)=i12
      call scovnce (npa,npb,ngr,npx,ela,elb,p1a,p1b,probs,ig,12,lw)
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
      call timer (4,isndel,'_sndelta  ',-1)
      return
 112  format (1x,'#',/,1x,'# Multiple-Fragment Covariances:',
     & ' COV (i,j..k) = < PROD (x=1..k) (n_x - <n_x>) >',/,' #')
      end
c
c-----------------------------------------------------------------------
c
      subroutine scovnce 
     &   (npa,npb,ngr,npx,ela,elb,p1a,p1b,probs,ig,nig,lw)
c
      include   'implicit.inc'
      include   'param.inc'
      include   'constants.inc'
c
      integer    ela(npa,ngr),elb(npb,ngr)
      real(kind=8)    p1a(ngr),p1b(ngr),probs(npx,2)
      integer    ig(20)
c
      dela=0.0d+00
      do i=1,npa
        prod=1.0d+00
        do ix=1,nig
          prod=prod*(ela(i,ig(ix))-p1a(ig(ix)))
        enddo
        dela=dela+prod*probs(i,1)
      enddo
      delb=0.0d+00
      do i=1,npb
        prod=1.0d+00
        do ix=1,nig
          prod=prod*(elb(i,ig(ix))-p1b(ig(ix)))
        enddo
        delb=delb+prod*probs(i,2)
      enddo
      delt=dela+delb
      write (lw,10) dela,delb,delt,(ig(k),k=1,nig)
      return
 10   format (1x,'COV (alpha,beta,sum) = ',3(1x,F16.10),
     &  4x,'GROUPS --> ',20(I3))
      end
