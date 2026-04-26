c
c-----This routine takes the ALPHA and BETA probabilities, stored
c     in 'filea' and 'fileb' respectively, of a molecule divided into
c     'ngroup' fragments and computes the SPINLESS probabilities.
c
c-----------------------------------------------------------------------
c
      subroutine pfromab (pcut,ngroup,lw,minp,maxp,filea,fileb)
      implicit none

      integer(kind=4), parameter :: nstack = 5000
      real   (kind=8) :: pcut,sumpcut
      real   (kind=8), allocatable,dimension (:)   :: pord
      integer(kind=4), allocatable,dimension (:)   :: iord
      integer(kind=4), allocatable,dimension (:,:) :: resncord
      integer(kind=4)    :: maxp(ngroup)
      integer(kind=4)    :: minp(ngroup)
      character(len=*)   :: filea,fileb
      integer(kind=4)    :: ngroup,ngroupx
      integer(kind=4)    :: lw
      character(len=200) :: line,word

      real   (kind=8), allocatable,dimension (:)  :: pnew
      real   (kind=8), allocatable,dimension (:)  :: probt
      integer(kind=4), allocatable,dimension(:,:) :: resnc
      integer(kind=4), allocatable,dimension(:,:) :: na,nb
      integer(kind=4), allocatable,dimension(:)   :: ind
      real   (kind=8), allocatable,dimension(:)   :: pa
      real   (kind=8), allocatable,dimension(:)   :: pb
      real   (kind=8), allocatable,dimension(:)   :: avga,avgb
      real   (kind=8), allocatable,dimension(:)   :: avgasq,avgbsq
      real   (kind=8), allocatable,dimension(:)   :: p2aa,p2ab
      real   (kind=8), allocatable,dimension(:)   :: p2ba,p2bb
      real   (kind=8), allocatable,dimension(:)   :: d1aa,d1ab
      real   (kind=8), allocatable,dimension(:)   :: d1ba,d1bb
      real   (kind=8), allocatable,dimension(:,:) :: xlis

      real   (kind=8)    :: pis,twop,terma,termb
      integer(kind=4)    :: i,npa,npb,np,m,j,ia,ib,k,n,nprob,ls,ios
      integer(kind=4)    :: ipair,iaib
      integer(kind=4)    :: napnb,ngreater
      logical            :: inlist,exfil

      ls  =  20
      exfil = .true.
      do while (exfil)
         inquire (unit=ls,opened=exfil)
         if (.not.exfil) then
           open (unit=ls,file=trim(filea),iostat=ios,status='old')
           if (ios.ne.0) then
             write (0,*) ' # pfromab.f: Error openning '//trim(filea)
             stop
           endif
         else
           ls=ls+1
         endif
      enddo
c
c-----Read ALPHA probabilities
c
      read (ls,*) ngroup,npa
      allocate (pa(npa))
      allocate (na(npa,ngroup))
      do j=1,npa
        read (ls,*) pa(j),(na(j,k),k=1,ngroup)
      enddo
      close (ls)

      exfil = .true.
      do while (exfil)
         inquire (unit=ls,opened=exfil)
         if (.not.exfil) then
           open (unit=ls,file=trim(fileb),iostat=ios,status='old')
           if (ios.ne.0) then
             write (0,*) ' # pfromab.f: Error openning '//trim(fileb)
             stop
           endif
         else
           ls=ls+1
         endif
      enddo
c
c-----Read BETA  probabilities
c
      read (ls,*) ngroupx,npb
      if (ngroup.ne.ngroupx) then
        stop '# pfromab.f: ALPHA/BETA probs different NGROUP values'
      endif
      allocate (pb(npb))
      allocate (nb(npb,ngroup))
      do j=1,npb
        read (ls,*) pb(j),(nb(j,k),k=1,ngroup)
      enddo
      close (ls)

      allocate (pnew(npa*npb))
      allocate (probt(npa*npb))
      allocate (resnc(npa*npb,ngroup))
      allocate (ind(ngroup))
c
c-----Write Spin-resolved probabilities
c
      allocate (avga(ngroup))
      allocate (avgb(ngroup))
      allocate (avgasq(ngroup))
      allocate (avgbsq(ngroup))
      allocate (p2aa(ngroup*(ngroup-1)/2))
      allocate (p2ab(ngroup*(ngroup-1)/2))
      allocate (p2ba(ngroup*(ngroup-1)/2))
      allocate (p2bb(ngroup*(ngroup-1)/2))
      allocate (d1aa(ngroup*(ngroup-1)/2))
      allocate (d1ab(ngroup*(ngroup-1)/2))
      allocate (d1ba(ngroup*(ngroup-1)/2))
      allocate (d1bb(ngroup*(ngroup-1)/2))
      allocate (xlis(ngroup,2))
      avga = 0.0d+00
      avgb = 0.0d+00
      avgasq = 0.0d+00
      avgbsq = 0.0d+00
      p2aa = 0.0d+00
      p2ab = 0.0d+00
      p2ba = 0.0d+00
      p2bb = 0.0d+00
      d1aa = 0.0d+00
      d1ab = 0.0d+00
      d1ba = 0.0d+00
      d1bb = 0.0d+00
      xlis = 0.0d+00
      probt = 0.0d+00
      write (lw,11) 
      iaib = 0
      do ia=1,npa
        do ib=1,npb
          iaib = iaib + 1
          probt(iaib) = pa(ia)*pb(ib)
          pis = pa(ia)*pb(ib)
          do j=1,ngroup
            avga(j)=avga(j)+pis*na(ia,j)
            avgb(j)=avgb(j)+pis*nb(ib,j)
            avgasq(j)=avgasq(j)+pis*na(ia,j)*na(ia,j)
            avgbsq(j)=avgbsq(j)+pis*nb(ib,j)*nb(ib,j)
            xlis(j,1)=xlis(j,1)+pis*na(ia,j)*na(ia,j)
            xlis(j,2)=xlis(j,2)+pis*nb(ib,j)*nb(ib,j)
          enddo
          write (lw,10) pis,(na(ia,j),nb(ib,j),j=1,ngroup)
          if (ngroup.gt.1) then
            ipair=0
            do i=1,ngroup
              do j=1,i-1
                ipair=ipair+1
                p2aa(ipair)=p2aa(ipair)+na(ia,i)*na(ia,j)*pis
                p2ab(ipair)=p2ab(ipair)+na(ia,i)*nb(ib,j)*pis
                p2ba(ipair)=p2ba(ipair)+nb(ib,i)*na(ia,j)*pis
                p2bb(ipair)=p2bb(ipair)+nb(ib,i)*nb(ib,j)*pis
              enddo
            enddo
          endif
        enddo
      enddo
      write (lw,119)
      do i=1,ngroup
        write (lw,19) i,avga(i),i,avgb(i)
      enddo
      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (lw,29)
     &       i,j,p2aa(ipair),i,j,p2ab(ipair),
     &       i,j,p2ba(ipair),i,j,p2bb(ipair)
          enddo
        enddo
      endif
      np=0
      do ia=1,npa
        do ib=1,npb
          np=np+1
          pis = pa(ia)*pb(ib)
          twop = pis+pis
          if (ngroup.gt.1) then
            ipair=0
            do i=1,ngroup
              do j=1,i-1
                ipair=ipair+1
                d1aa(ipair)=d1aa(ipair)-
     &                   twop*(avga(i)-na(ia,i))*(avga(j)-na(ia,j))
                d1ab(ipair)=d1ab(ipair)-
     &                   twop*(avga(i)-na(ia,i))*(avgb(j)-nb(ib,j))
                d1ba(ipair)=d1ba(ipair)-
     &                   twop*(avgb(i)-nb(ib,i))*(avga(j)-na(ia,j))
                d1bb(ipair)=d1bb(ipair)-
     &                   twop*(avgb(i)-nb(ib,i))*(avgb(j)-nb(ib,j))
              enddo
            enddo
          endif
        enddo
      enddo

      do i=1,ngroup
        xlis(i,1)=-(xlis(i,1)-avga(i)*(avga(i)+1.0d+00))
        xlis(i,2)=-(xlis(i,2)-avgb(i)*(avgb(i)+1.0d+00))
        terma=0.0d+00
        termb=0.0d+00
        if (avga(i).gt.1.0D-08) terma=100.0d+00*xlis(i,1)/avga(i)
        if (avgb(i).gt.1.0D-08) termb=100.0d+00*xlis(i,2)/avgb(i)
        write (lw,69) i,i,xlis(i,1),terma,i,i,xlis(i,2),termb
      enddo
      if (ngroup.gt.1) then
        ipair=0
        do i=1,ngroup
          do j=1,i-1
            ipair=ipair+1
            write (lw,49) i,j,d1aa(ipair),i,j,d1ab(ipair),
     &         i,j,d1ba(ipair),i,j,d1bb(ipair)
          enddo
        enddo
      endif
c
      pnew=0d0
      n=0
      do ia=1,npa
        cycleib: do ib=1,npb
c
c---------Skip this probability if minp(i) > elec(i) or maxp(i) < elec(i)
c
          do k=1,ngroup
            napnb=na(ia,k)+nb(ib,k)
            if (napnb.gt.maxp(k).or.napnb.lt.minp(k)) cycle cycleib
          enddo
          n=n+1
          if (n.eq.1) then
            np=1
            resnc(np,1:ngroup)=na(ia,1:ngroup)+nb(ib,1:ngroup)
            pnew(np)=pa(ia)*pb(ib)
          else
            ind(1:ngroup)=na(ia,1:ngroup)+nb(ib,1:ngroup)
            do j=1,np
              inlist=.true.
              do k=1,ngroup
                inlist=inlist.and.(ind(k).eq.resnc(j,k))
              enddo
              if (inlist) then
                pnew(j)=pnew(j)+pa(ia)*pb(ib)
                cycle cycleib
              endif
            enddo
            np=np+1
            pnew(np)=pa(ia)*pb(ib)
            resnc(np,1:ngroup)=na(ia,1:ngroup)+nb(ib,1:ngroup)
          endif
        enddo cycleib
      enddo

      allocate (pord(np))
      allocate (iord(np))
      allocate (resncord(np,ngroup))
      pord(1:np)=pnew(1:np)
      forall (i=1:np) iord(i)=i
      nprob=min(np,nstack)
      call qqsort (pord,iord,1,nprob,nprob)
      write (lw,20) nprob,pcut
      sumpcut=0d0
      ngreater=0
      do n=nprob,1,-1
        m=iord(n)
        if (pord(m).ge.pcut) then
          sumpcut=sumpcut+pord(m)
          ngreater=ngreater+1
          write (lw,10) pord(m),(resnc(m,j),j=1,ngroup)
        endif
      enddo
      write (lw,3) sum(pord(1:nprob)),sumpcut,pcut
      write (lw,4) ngreater

cAQUI
      call ndelta (probt,avga,avgb,na,nb,ngroup,npa,npb,lw)
      call locdeloc (pnew,resnc,npa,npb,np,ngroup,lw)
      deallocate (avga)
      deallocate (avgb)
      deallocate (avgasq)
      deallocate (avgbsq)
      deallocate (p2aa)
      deallocate (p2ab)
      deallocate (p2ba)
      deallocate (p2bb)
      deallocate (d1aa)
      deallocate (d1ab)
      deallocate (d1ba)
      deallocate (d1bb)
      deallocate (xlis)
      deallocate (probt)
      deallocate (pnew)
      deallocate (pord)
cAQUI



 10   format (1x,'# ',F22.15,3x,20I6)
 20   format (//,' #',6x,'THERE ARE ',I8,' Spinless probabilities',/,
     & ' #',6x,'Printing probabilities greater than ',F15.8,/,
     & ' #',6x,'probability             RSRS --->',/,' #',6x,80('-'))
 3    format (1x,'# ',F22.15,2x,'<--- SUM',/,
     &        1x,'# ',F22.15,2x,'<--- SUM ( p >',1x,F15.8,' )')
 4    format (1x,'# ',I10,' probabilities printed')
 11   format (/,' # Spin-Resolved probabilities',/,' # ',27('-'),/,
     &  ' #',12x,'Prob',11x,'G1(a) G1(b) ...',/,' #',6x,80('-'))
 119  format (//1x,'Average populations and delocalization indices')
 19   format (1x,'<n(',I3,')_alpha>',17x,'= ',F16.9,/,
     .   1x,'<n(',I3,')_beta>',18x,'= ',F16.9)
 29   format (/1x,'<n(',I3,')_alpha n(',I3,')_alpha>     = ',F16.9,/,
     .         1x,'<n(',I3,')_alpha n(',I3,')_beta>      = ',F16.9,/,
     .         1x,'<n(',I3,')_beta  n(',I3,')_alpha>     = ',F16.9,/,
     .         1x,'<n(',I3,')_beta  n(',I3,')_beta>      = ',F16.9)
 69   format (/1x,'delta_(',I3,I3,')_alpha            = ',F16.9,
     .         2x,'% Localization = ',F8.4,/,
     .         1x,'delta_(',I3,I3,')_beta             = ',F16.9,
     .         2x,'% Localization = ',F8.4)
 49   format (/1x,'delta_(',I3,I3,')_{alpha,alpha)    = ',F16.9,/,
     .         1x,'delta_(',I3,I3,')_{alpha,beta )    = ',F16.9,/
     .         1x,'delta_(',I3,I3,')_{beta ,alpha)    = ',F16.9,/
     .         1x,'delta_(',I3,I3,')_{beta ,beta )    = ',F16.9)
      end
