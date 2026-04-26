c
c-----------------------------------------------------------------------
c
      subroutine aomtest (sg,nmo,nfrag,critics,lw,verbose)
c
c-----------------------------------------------------------------------
c                         
      include 'implicit.inc'
      include 'constants.inc'
      include 'error.inc'
c
      real(kind=8)   sg(nfrag,nmo,nmo)

      integer(kind=4), allocatable,dimension (:)   :: ixion,nbatom,ncom
      integer(kind=4), allocatable,dimension (:,:) :: blg
      logical   condition, first, verbose
      real(kind=8) crits(29)
      data crits /95D0,10D0,20D0,30D0,40D0,50D0,60D0,70D0,80D0,
     &  90D0,91D0,92D0,93D0,94D0,95D0,96D0,97D0,
     &  98D0,99D0,99.1D0,99.2D0,99.3D0,99.4D0,99.5D0, 
     &  99.6D0,99.7D0,99.8D0,99.9D0,99.99D0/
c
      call timer (2,iaomtest,'_aomtest  ',-1)
c
c.....Allocate arrays
c
      n = nmo
      allocate (ixion(nmo),stat=ier)
      if (ier.ne.0) stop 'aomtest.f: Cannot allocate ixion()'
      allocate (blg(nfrag,nmo),stat=ier)
      if (ier.ne.0) stop 'aomtest.f: Cannot allocate blg()'
      allocate (nbatom(nfrag),stat=ier)
      if (ier.ne.0) stop 'aomtest.f: Cannot allocate nbatom()'
      allocate (ncom(nmo),stat=ier)
      if (ier.ne.0) stop 'aomtest.f: Cannot allocate ncom()'
c
c-------------------------------------------------------------------
c
c-----Write diagonal AOM elements in each atom or fragment
c
      write (lw,1025)
      do k=1,nfrag
        write (lw,*) '# Center ',k
        write (lw,30) (sg(k,ip,ip),ip=1,n)
      enddo
      if (.not.verbose) return
c
c-----Write the MOs partially localized on each atom or fragment
c
      write (lw,1024) critics
      do k=1,nfrag
        nbatom(k)=0
        do ip=1,n
          if (abs(sg(k,ip,ip)).gt.critics) then
            nbatom(k)=nbatom(k)+1
            blg(k,nbatom(k))=ip
          endif
        enddo
        write (lw,1022) k,nbatom(k),(blg(k,j),j=1,nbatom(k))
      enddo
c
c-----MOs localized simultaneously in two atoms or fragments
c
      write (lw,122) critics,critics
      do k1=1,nfrag-1
        do k2=k1+1,nfrag
          na=nbatom(k1)
          nb=nbatom(k2)
          nab=max(na,nb)
          call sharedmo (blg(k1,:),blg(k2,:),ncom,na,nb,nab,nc)
          if (nc.gt.0) then
            write (lw,121) k1,' and ',k2,nc,(ncom(k),k=1,nc)
          endif
        enddo
      enddo
c
c-----MOs localized in one atom or the other of a given pair
c
      write (lw,124) critics,critics
      do k1=1,nfrag-1
        do k2=k1+1,nfrag
          na=nbatom(k1)
          nb=nbatom(k2)
          nab=nmo
          call unionmo (blg(k1,:),blg(k2,:),ncom,na,nb,nab,nc)
          write (lw,121) k1,' or  ',k2,nc,(ncom(k),k=1,nc)
        enddo
      enddo
c
      critov=0.95D0
      if (nfrag.ge.2) write (lw,1028) critov*100D0,'pair   '
      do k1=2,nfrag
        do k2=1,k1-1
          nb=0
          do ip=1,n
            condition = sg(k1,ip,ip)+sg(k2,ip,ip).gt.critov
            condition = condition .and. sg(k1,ip,ip).lt.critov
            condition = condition .and. sg(k2,ip,ip).lt.critov
            if (condition) then
              if (nb.eq.0) write (lw,1029) k1,k2
              nb=nb+1
              ixion(nb)=ip
              write (lw,1030) ixion(nb),
     &             100*sg(k1,ip,ip),100*sg(k2,ip,ip)
            endif
          enddo
        enddo
      enddo

      critov=0.95D0
      if (nfrag.ge.3) write (lw,1028) critov*100D0,'trio   '
      do k1=3,nfrag
        do k2=1,k1-1
          do k3=1,k2-1
            nb=0
            do ip=1,n
              s1=sg(k1,ip,ip)
              s2=sg(k2,ip,ip)
              s3=sg(k3,ip,ip)
              condition = s1+s2+s3.gt.critov
              condition = condition .and. s1.lt.critov
              condition = condition .and. s2.lt.critov
              condition = condition .and. s3.lt.critov
              condition = condition .and. (s1+s2).lt.critov
              condition = condition .and. (s1+s3).lt.critov
              condition = condition .and. (s2+s3).lt.critov
              if (condition) then
                if (nb.eq.0) write (lw,1029) k1,k2,k3
                nb=nb+1
                ixion(nb)=ip
                write (lw,1033) ixion(nb),100*s1,100*s2,100*s3
              endif
            enddo
          enddo
        enddo
      enddo

      critov=0.95D0
      if (nfrag.ge.4) write (lw,1028) critov*100D0,'quartet'
      do k1=4,nfrag
        do k2=1,k1-1
          do k3=1,k2-1
            do k4=1,k3-1
              nb=0
              do ip=1,n
                s1=sg(k1,ip,ip)
                s2=sg(k2,ip,ip)
                s3=sg(k3,ip,ip)
                s4=sg(k4,ip,ip)
                condition = s1+s2+s3+s4.gt.critov
                condition = condition .and. s1.lt.critov
                condition = condition .and. s2.lt.critov
                condition = condition .and. s3.lt.critov
                condition = condition .and. s4.lt.critov
                condition = condition .and. (s1+s2).lt.critov
                condition = condition .and. (s1+s3).lt.critov
                condition = condition .and. (s1+s4).lt.critov
                condition = condition .and. (s2+s3).lt.critov
                condition = condition .and. (s2+s4).lt.critov
                condition = condition .and. (s3+s4).lt.critov
                condition = condition .and. (s1+s2+s3).lt.critov
                condition = condition .and. (s1+s2+s4).lt.critov
                condition = condition .and. (s1+s3+s4).lt.critov
                condition = condition .and. (s2+s3+s4).lt.critov
                if (condition) then
                  if (nb.eq.0) write (lw,1029) k1,k2,k3,k4
                  nb=nb+1
                  ixion(nb)=ip
                  write (lw,1036) ixion(nb),
     &              100*s1,100*s2,100*s3,100*s4
                endif
              enddo
            enddo
          enddo
        enddo
      enddo

 
      do is=1,29
        critov = crits(is)*0.01D0
        first = .true.
        nbtot=0
        do k=1,nfrag
          nb=0
          do ip=1,n
            if (abs(sg(k,ip,ip)).gt.critov) then
              nb=nb+1
              ixion(nb)=ip
            endif
          enddo
          nbtot=nbtot+nb
          if (nb.gt.0) then 
            if (first) then
              write (lw,1026) critov*100D0
              first=.false.
            endif
            if (nb.le.10) then
              write (lw,1027) k,(ixion(j),j=1,nb)
            else
              write (lw,1027) k,(ixion(j),j=1,10)
              write (lw,1043)   (ixion(j),j=11,nb)
            endif
          endif
        enddo
      enddo
c
      deallocate (ixion,stat=ier)
      if (ier.ne.0) stop 'aomtest.f: Cannot deallocate ixion()'
      deallocate (blg,stat=ier)
      if (ier.ne.0) stop 'aomtest.f: Cannot deallocate blg()'
      deallocate (nbatom,stat=ier)
      if (ier.ne.0) stop 'aomtest.f: Cannot deallocate nbatom()'
      deallocate (ncom,stat=ier)
      if (ier.ne.0) stop 'aomtest.f: Cannot deallocate ncom()'

      call timer (4,iaomtest,'_aomtest  ',-1)
      return
c
c.....Formats
c
 30   format (5(1x,F15.8))
 31   format (2(/,2(5x,E14.7)))
 1022 format (3x,'Frag ',I2,' has ',I3,' MOs: ',20I3,/,
     &        1000(24x,20I3,/))
 1032 format (1000(23x,10I4,/))
 1024 format (/1x,'# MOs partially localized on each fragment',
     &        '  (S_ii^A > ',1PE15.7,' )')
 1026 format (//' # MOs localized more than',F7.2,'% on each fragment')
 1027 format (' # Fragment ',I4,' --> ',10I4)
 1043 format (1000(21x,10I4))
 1025 format (/1x,'# Diagonal overlaps in each fragment')
 1028 format (/' # MOs localized more than',F7.2,
     &'% on a ',a,' of fragments')
 1029 format (' # Fragments ',10I4)
 1030 format (3x,'MO ',I4,5x,F5.2,'% + ',F5.2,'%')
 1033 format (3x,'MO ',I4,5x,F5.2,'% + ',F5.2,'%+ ',F5.2,'%')
 1036 format (3x,
     & 'MO ',I4,5x,F5.2,'% + ',F5.2,'%+ ',F5.2,'%+ ',F5.2,'%')
 121  format ('   Frags ',I2,a,I2,' have ',I3,' MOs:',
     &  20I3,/,1000(32x,20I3,/))
 131  format (1000(42x,10I4,/))
 122  format (/,' # MOs partially localized in A & B (Intersection)',/,
     & ' # S_ii^A > ',1PE15.7,' and S_ii^B > ',1PE15.7) 
 124  format (/,' # MOs partially localized in A .or. B (Union)',/,
     & ' #  (S_ii^A > ',1PE15.7,' or  S_ii^B > ',1PE15.7,' )')
      end
c
c-----Computes the intersection of the ordered arrays of integer 
c     numbers a[] and b[]. The nc common elements are returned in c[]
c
      subroutine sharedmo (a,b,c,na,nb,nab,nc)
      implicit none
      integer a(na),b(nb),c(nab),na,nb,nab,nc,i,j,k

      i=1
      j=1
      k=0
      do while (i.le.na .and. j.le.nb)
        if (a(i) < b(j)) then
          i=i+1
        elseif (b(j) < a(i)) then
          j=j+1
        elseif (a(i) == b(j)) then
          k=k+1
          c(k)=a(i)
          i=i+1
          j=j+1
        endif
      enddo
      nc=k
      return
      end
c
c-----Computes the union of the ordered arrays of integer  numbers a[] 
c     and b[]. The nc that are in a[] or in b[] are returned in c[]
c
c
      subroutine unionmo (a,b,c,na,nb,nab,nc)
      implicit none
      integer a(na),b(nb),c(nab),na,nb,nab,nc,i,j,k

      i=1
      j=1
      k=0
      do while (i.le.na .and. j.le.nb)
        if (a(i) < b(j)) then
          k=k+1
          c(k)=a(i)
          i=i+1
        elseif (b(j) < a(i)) then
          k=k+1
          c(k)=b(j)
          j=j+1
        else
          k=k+1
          c(k)=a(i)
          i=i+1
          j=j+1
        endif
      enddo
      do while (i.le.na) 
        k=k+1
        c(k)=a(i)
        i=i+1
      enddo
      do while (j.le.nb) 
        k=k+1
        c(k)=b(j)
        j=j+1
      enddo
      nc=k
      return
      end
