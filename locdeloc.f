c
c-----------------------------------------------------------------------
c
      subroutine locdeloc (pnew,resnc,npa,npb,np,ngroup,lw)
      implicit none
      integer(kind=4) :: np,npa,npb,lw,ngroup,nprob
      real   (kind=8) :: pnew(npa*npb)
      integer(kind=4) :: resnc(npa*npb,ngroup)
c
      real   (kind=8), allocatable,dimension(:) :: p1,p2,p3,xli
      real   (kind=8), allocatable,dimension(:) :: d1,d2
      integer(kind=4), allocatable,dimension(:) :: npop
      real   (kind=8) :: pis
      integer(kind=4) :: i,j,k,iter,ipair

!
!.....Compute spinless localization and delocalization indices
!
      allocate (p1(ngroup))
      allocate (xli(ngroup))
      p1 = 0.0d+00
      xli = 0.0d+00
      allocate (npop(ngroup))
      if (ngroup.gt.1) then
        allocate (p2(ngroup*(ngroup-1)/2))
        allocate (d1(ngroup*(ngroup-1)/2))
        p2=0.0d+00
        d1=0.0d+00
      endif
      if (ngroup.gt.2) then
        allocate (p3(ngroup*(ngroup-1)*(ngroup-2)/6))
        allocate (d2(ngroup*(ngroup-1)*(ngroup-2)/6))
        p3=0.0d+00
        d2=0.0d+00
      endif
 
      nprob=np
      do np=1,nprob
        pis=pnew(np)
        npop(:)=resnc(np,:)
        do i=1,ngroup
          p1(i)=p1(i)+npop(i)*pis
          xli(i)=xli(i)+pis*npop(i)*npop(i)
        enddo
        if (ngroup.gt.1) then
          ipair=0
          do i=1,ngroup
            do j=1,i-1
              ipair=ipair+1
              p2(ipair)=p2(ipair)+npop(i)*npop(j)*pis
            enddo
          enddo
        endif
        if (ngroup.gt.2) then
          iter=0
          do i=1,ngroup
            do j=1,i-1
              do k=1,j-1
              iter=iter+1
              p3(iter)=p3(iter)+npop(i)*npop(j)*npop(k)*pis
              enddo
            enddo
          enddo
        endif
      enddo
 
      write (lw,11)
      do i=1,ngroup
        write (lw,15) i,p1(i)
      enddo

      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (lw,2) i,j,p2(ipair)
          enddo
        enddo
      endif
 
      if (ngroup.gt.2) then
        iter=0
        do i=1,ngroup
          do j=1,i-1
            do k=1,j-1
            iter=iter+1
            write (lw,3) i,j,k,p3(iter)
            enddo
          enddo
        enddo
      endif
 
      do np=1,nprob
        pis=pnew(np)
c
c-------Add core electrons to each resonant structure
c
        npop(:)=resnc(np,:)
        if (ngroup.gt.1) then
          ipair=0
          do i=1,ngroup
            do j=1,i-1
              ipair=ipair+1
              d1(ipair)=d1(ipair)-2*pis*(p1(i)-npop(i))*(p1(j)-npop(j))
            enddo
          enddo
        endif
        if (ngroup.gt.2) then
          iter=0
          do i=1,ngroup
            do j=1,i-1
              do k=1,j-1
                iter=iter+1
                d2(iter)=d2(iter)-
     .          2*pis*(p1(i)-npop(i))*(p1(j)-npop(j))*(p1(k)-npop(k))
              enddo
            enddo
          enddo
        endif
      enddo  
 
      do i=1,ngroup
        xli(i)=-(xli(i)-p1(i)*(p1(i)+1.0d+00))
        write (lw,6) i,i,xli(i),100.0d+00*xli(i)/p1(i)
      enddo
      write (lw,33)
      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (lw,4) i,j,d1(ipair)
          enddo
        enddo
      endif
      if (ngroup.gt.2) then
        iter=0
        do i=1,ngroup
          do j=1,i-1
            do k=1,j-1
            iter=iter+1
               write (lw,5) i,j,k,d2(iter)
            enddo
          enddo
        enddo
      endif
      return
 11   format (//1x,'Average populations and localization indices')
 15   format (1x,'# <n(',I3,')>               = ',F16.9)
 2    format (1x,'# <n(',I3,') n(',I3,')>        = ',F16.9)
 3    format (1x,'# <n(',I3,') n(',I3,') n(',I3,')> = ',F16.9)
 6    format (1x,'# delta_(',I3,I3,')         = ',F16.9,
     .        2x,'% Localization = ',F8.4)
 33   format (//1x,'Delocalization indices,',
     .             ' Eq. (28) J. Chem. Phys.  126, 094102 (2007)')
 4    format (1x,'# delta_(',I3,I3,')         = ',F16.9)
 5    format (1x,'# delta_(',I3,I3,I3,')      = ',F16.9)
      end
