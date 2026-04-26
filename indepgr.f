c
c-----------------------------------------------------------------------
c
      subroutine indepgr (
     &  epsdet,probcut,ndets,nelact,nact,ncent,nmo,ncore,ngroup,mal,
     &  mbe,mocore,nel,moval,ifilc,lw,lr,aom,nfugrp,ifugrp,mocogrp,
     &  ival,doentropy,orderp,short,wfnf)
      include     'implicit.inc'
      include     'param.inc'
      include     'mline.inc'
      real(kind=8)  aom(ncent,nmo,nmo)
      integer nfugrp(ngroup)
      integer ifugrp(ncent,ngroup)
      integer mocogrp(ngroup)
      integer ival(nmo)
      logical doentropy,orderp,short
*     parameter  (mline=200)
      character*(mline) wfnf

      real(kind=8),  allocatable,dimension (:,:,:) :: sgi 
      real(kind=8),  allocatable,dimension (:,:)   :: pnew
      real(kind=8),  allocatable,dimension (:)     :: pord
      real(kind=8),  allocatable,dimension (:)     :: p1,p2,d1
      integer, allocatable,dimension (:)     :: ioprob
      integer, allocatable,dimension (:,:)   :: resnca,resncb,resnc
      integer, allocatable,dimension (:)     :: jloop,iloop,floop,sloop
      integer coregrp(2)
c
c-----------------------------------------------------------------------
c
      nval=nel-2*mocore
      naval=mal-mocore
      nbval=mbe-mocore

      write (lw,'(1x,"# ",80("+"))')
      write (lw,23) 
      do i=1,ngroup
        write (lw,232) i,(ifugrp(k,i),k=1,nfugrp(i))
      enddo
      write (lw,'(1x,"# ",80("+"),/,1x)')
      nproba=mal-mocore+1
      nprobb=mbe-mocore+1
      if (.not.allocated(resnca)) then
        allocate (resnca(nproba,2),stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot allocate resnca()'
      endif
      if (.not.allocated(resncb)) then
        allocate (resncb(nprobb,2),stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot allocate resncb()'
      endif
      if (.not.allocated(resnc)) then
        allocate (resnc(nstack,ngroup),stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot allocate resnc()'
      endif
      if (.not.allocated(sgi)) then
        allocate (sgi(2,nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot allocate sgi()'
      endif

      do i=0,nproba-1
        resnca(i+1,1)=i
        resnca(i+1,2)=nproba-i-1
      enddo
      do i=0,nprobb-1
        resncb(i+1,1)=i
        resncb(i+1,2)=nprobb-i-1
      enddo

      combi=1D0
      do i=0,nval-1
        combi=combi*dble(nval+1-i)/dble(nval-i)
      enddo
      nprobt=int(combi+1D-10)
      if (.not.allocated(pnew)) then
        allocate (pnew(nprobt,ngroup),stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot allocate pnew()'
      endif
      if (.not.allocated(pord)) then
        allocate (pord(nstack),stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot allocate pord()'
      endif
      if (.not.allocated(ioprob)) then
        allocate (ioprob(nstack),stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot allocate ioprob()'
      endif
      if (.not.allocated(p1)) then
        allocate (p1(ngroup),stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot allocate p1()'
      endif
      if (.not.allocated(p2)) then
        allocate (p2(ngroup*(ngroup-1)/2),stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot allocate p2()'
      endif
      if (.not.allocated(d1)) then
        allocate (d1(ngroup*(ngroup-1)/2),stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot allocate d1()'
      endif
c
      do ig=1,ngroup-1
        sgi=0D0
        do m=1,nmo
          do j=1,nmo
            delij=0D0
            if (m.eq.j) delij=1D0
            sgi(1,m,j)=0D0
            do i=1,nfugrp(ig)
              sgi(1,m,j)=sgi(1,m,j)+aom(ifugrp(i,ig),m,j)
            enddo
            sgi(2,m,j)=delij-sgi(1,m,j)
          enddo
        enddo
        coregrp(1)=mocogrp(1)
        coregrp(2)=0
        do i=2,ngroup
          coregrp(2)=coregrp(2)+mocogrp(i)
        enddo
        write (lw,24) ig
        if (icorr) then
          call cindgr (epsdet,probcut,nproba,nprobb,2,ndets,nmo,
     &      ncore,nelact,nel,moval,ival,mocogrp,naval,nbval,
     &      ifilc,lw,wfnf,sgi,resnca,resncb,orderp,.true.,
     &      pnew(1,ig),nprobt)
        else
          call nbasabw (probcut,nproba,nprobb,2,nmo,short,moval,ival,
     &      coregrp,nel,nval,naval,nbval,ifilc,lw,sgi,resnca,resncb,
     &      doentropy,orderp,pnew(1,ig),nprobt)
        endif
      enddo

      allocate (jloop(ngroup  ))
      allocate (iloop(ngroup-1))
      allocate (floop(ngroup-1))
      allocate (sloop(ngroup-1))

      do k=1,ngroup-1
        iloop(k)=1
        floop(k)=nprobt
        sloop(k)=1
        jloop(k)=iloop(k)
      end do
      k=ngroup-1
      p1(1:ngroup)=0D0
      p2(1:ngroup*(ngroup-1)/2)=0D0
      d1(1:ngroup*(ngroup-1)/2)=0D0

      isto=0
      istx=0
      do while (k.gt.0)
        if ((jloop(k)-floop(k))*sloop(k).gt.0) then 
           jloop(k)=iloop(k)
           k=k-1
        else
           p=1D0
           jsum=0
           do i=1,ngroup-1
             jsum=jsum+jloop(i)-1
             p=p*pnew(jloop(i),i)
           enddo
           if (jsum.le.nel) then
             jloop(ngroup)=nel-jsum+1
             if (p.gt.probcut) then
               if (isto.lt.nstack) then
                 isto=isto+1
                 ioprob(isto)=isto
                 pord(isto)=p
                 resnc(isto,1:ngroup)=jloop(1:ngroup)-1
               else
                 istx=istx+1
                 if (istx.eq.1) write (lw,231) 
                 write (lw,22) p,(jloop(i)-1,i=1,ngroup)
               endif
             endif
             do i=1,ngroup
               p1(i)=p1(i)+p*(jloop(i)-1)
             enddo
             if (ngroup.gt.1) then
               ipair=0
               do i=1,ngroup
                 do j=1,i-1
                   ipair=ipair+1
                   p2(ipair)=p2(ipair)+p*(jloop(i)-1)*(jloop(j)-1)
                 enddo
               enddo
             endif
           endif
           k=ngroup-1
        end if
        if (k.gt.0) jloop(k)=jloop(k)+sloop(k)
      end do
      write (lw,'(1x,"# ",80("+"),/,1x)')

      write (lw,233)
      call qqsort (pord,ioprob,1,isto,nstack)
      spord=0D0
      spordtot=0D0
      ngtprob=0
      do n=isto,1,-1
        m=ioprob(n)
        spordtot=spordtot+pord(m)
        if (abs(pord(m)).gt.probcut) then
          spord=spord+pord(m)
          ngtprob=ngtprob+1
          write (lw,22) pord(m),(resnc(m,i),i=1,ngroup) 
        endif
      enddo
      write (lw,'(1x,"# ",80("+"),/,1x)')
      write (lw,234) spord,ngtprob,probcut,spordtot

      do k=1,ngroup-1
        iloop(k)=1
        floop(k)=nprobt
        sloop(k)=1
        jloop(k)=iloop(k)
      end do
      k=ngroup-1
      do while (k.gt.0)
        if ((jloop(k)-floop(k))*sloop(k).gt.0) then 
           jloop(k)=iloop(k)
           k=k-1
        else
          if (ngroup.gt.1) then
            p=1D0
            jsum=0
            do i=1,ngroup-1
              jsum=jsum+jloop(i)-1
              p=p*pnew(jloop(i),i)
            enddo
            if (jsum.le.nel) then
              jloop(ngroup)=nel-jsum+1
              ij=0
              do i=1,ngroup
                jli=jloop(i)-1
                do j=1,i-1
                  jlj=jloop(j)-1
                  ij=ij+1
                  d1(ij)=d1(ij)-2D0*p*(p1(i)-jli)*(p1(j)-jlj)
                enddo
              enddo
            endif
          endif
          k=ngroup-1
        end if
        if (k.gt.0) jloop(k)=jloop(k)+sloop(k)
      end do

      write (lw,11)
      do i=1,ngroup
        write (lw,15) i,p1(i)
      enddo
      if (ngroup.gt.1) then
        ij=0
        do i=1,ngroup
          do j=1,i-1
            ij=ij+1
            write (lw,2) i,j,p2(ij)
          enddo
        enddo
        write (lw,33)
        ij=0
        do i=1,ngroup
          do j=1,i-1
            ij=ij+1
            write (lw,4) i,j,d1(ij)
          enddo
        enddo
      endif
      write (lw,*) '#'
      write (lw,*) '#'
      write (lw,*) '#'

      deallocate (jloop)
      deallocate (iloop)
      deallocate (floop)
      deallocate (sloop)

      if (allocated(resnca)) then
        deallocate (resnca,stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot deallocate resnca()'
      endif
      if (allocated(resncb)) then
        deallocate (resncb,stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot deallocate resncb()'
      endif
      if (allocated(resnc)) then
        deallocate (resnc,stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot deallocate resnc()'
      endif
      if (allocated(sgi)) then
        deallocate (sgi,stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot deallocate sgi()'
      endif
      if (allocated(p1)) then
        deallocate (p1,stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot deallocate p1()'
      endif
      if (allocated(p2)) then
        deallocate (p2,stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot deallocate p2()'
      endif
      if (allocated(d1)) then
        deallocate (d1,stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot deallocate d1()'
      endif
      if (allocated(pnew)) then
        deallocate (pnew,stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot deallocate pnew()'
      endif
      if (allocated(pord)) then
        deallocate (pord,stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot deallocate pord()'
      endif
      if (allocated(ioprob)) then
        deallocate (ioprob,stat=ier)
        if (ier.ne.0) stop 'indepgr.f: Cannot deallocate ioprob()'
      endif

      return

 2    format (1x,'<n(',I3,') n(',I3,')>        = ',F16.10)
 4    format (1x,'delta_(',I3,I3,')         = ',F16.10)
 11   format (1x,'Average populations and pair populations')
 15   format (1x,'<n(',I3,')>               = ',F16.10)
 22   format (' # Probability = ',F22.16,4x,'for RSRS ',20I6)
 23   format (' # EDF assuming that electron populations '
     &          'groups 1,2,...,NGROUP-1 are not coupled')
 24   format (' # Computing EDF (Group',I2,',Rest of the molecule)')
 33   format (1x,'Delocalization indices,',
     &           ' Eq. (28) J. Chem. Phys.  126, 094102 (2007)')
 231  format (' # ',80('-'),/,' # ',15x,
     & ' PROBABILITIES NOT INCLUDED IN THE ORDERED SET BELOW ',
     &  /,' # ',80('-')) 
 232  format (' # GROUP ',I3,'  with atoms ',10I4,/,
     &  1000(' # ',22x,10I4,/))   
 233  format (' # ',25('-'),' ORDERED SET OF PROBABILITIES ',25('-')) 
 234  format (1x,F26.16,' <-- ','SUM,  ',I8,' PROBABILITIES > ',
     & F16.10,/,1x,F26.16,' <-- TOTAL SUM',/)
      end
