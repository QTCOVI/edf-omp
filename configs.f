c
c-----------------------------------------------------------------------
c
      subroutine configs (ndets,nmo,nelact,nact,ifilc,stderr,wfnfile)
c
      USE        space_for_cidet
      USE        space_for_conf
      include   'implicit.inc'
      include   'constants.inc'
      include   'mline.inc'
      integer, allocatable,dimension (:)     :: iorda,iordb
      integer, allocatable,dimension (:)     :: ordia,ordib,ordja,ordjb
      integer, allocatable,dimension (:)     :: nordia,nordib
      integer, allocatable,dimension (:)     :: idifa,idifb
      integer    stderr
      logical    inlist,exfil
*     parameter  (mline=200)
      character*(mline) wfnfile
c
      call timer (2,iconfigs,'_configs  ',-1)
c
c.....Analyze the wave function and determine the different ALPHA/BETA 
c     configurations which appear on it. Store these configurations in
c     two arrays. Besides this, the arrays kwa() and kwb() store the 
c     order number of the ALPHA and BETA configuration appearing in 
c     each Slater determinant.
c
      nelal=0
      nelbe=0
      read (ifilc,rec=1) (cidet(m),m=0,nelact)
      do k=1,nelact
        mm=int(cidet(k))
        if (mm.gt.0) nelal=nelal+1
        if (mm.lt.0) nelbe=nelbe+1
      enddo
      quot=1d0
      do i=0,nelal-1
        quot=quot*dble(nact-i)/dble(i+1)
      enddo
      ndetsa=min(nint(quot),ndets)
      quot=1d0
      do i=0,nelbe-1
        quot=quot*dble(nact-i)/dble(i+1)
      enddo
      ndetsb=min(nint(quot),ndets)
      call allocate_space_for_conf (ndets,ndetsa,ndetsb,nmo)
c
      allocate (ordia(nmo),ordja(nmo),iorda(nmo))
      allocate (ordib(nmo),ordjb(nmo),iordb(nmo))
      allocate (nordia(nmo),nordib(nmo))
      allocate (idifa(nmo),idifb(nmo))
      xnorm=zero
      do i=1,ndets
         nelal=0
         nelbe=0
         read (ifilc,rec=i) cdet(i),(cidet(m),m=1,nelact)
         xnorm=xnorm+cdet(i)*cdet(i)
         do m=1,nelact
           mm=int(cidet(m))
           if (mm.gt.0) then
             nelal=nelal+1
             iorda(nelal)=nelal
             ordia(nelal)=abs(mm)
           else
             nelbe=nelbe+1
             iordb(nelbe)=nelbe
             ordib(nelbe)=abs(mm)
           endif
         enddo
c
c........Alpha electrons
c
         call iqcksort (ordia,iorda,nmo,1,nelal)
         nordia(1:nelal)=ordia(iorda(1:nelal))
         call colord (ordia,nordia,nelal,nmo,npera,nsiga(i),ndifa,idifa)
         if (i.eq.1) then
           npa=1
           kwa(i)=npa
           nconfa(npa,1:nelal)=nordia(1:nelal)
         else
           do k=1,npa
             inlist=.true.
             do j=1,nelal
               inlist=inlist.and.(nconfa(k,j).eq.nordia(j))
             enddo
             kwa(i)=k
             if (inlist) goto 67
           enddo
           npa=npa+1
           kwa(i)=npa
           nconfa(npa,1:nelal)=nordia(1:nelal)
         endif
c
c........Beta electrons
c
 67      call iqcksort (ordib,iordb,nmo,1,nelbe)
         nordib(1:nelbe)=ordib(iordb(1:nelbe))
         call colord (ordib,nordib,nelbe,nmo,nperb,nsigb(i),ndifb,idifb)
         if (i.eq.1) then
           npb=1
           kwb(i)=npb
           nconfb(npb,1:nelbe)=nordib(1:nelbe)
         else
           do k=1,npb
             inlist=.true.
             do j=1,nelbe
               inlist=inlist.and.(nconfb(k,j).eq.nordib(j))
             enddo
             kwb(i)=k
             if (inlist) goto 68
           enddo
           npb=npb+1
           kwb(i)=npb
           nconfb(npb,1:nelbe)=nordib(1:nelbe)
         endif
 68   enddo
      deallocate (ordia,ordja,ordib,ordjb,iorda,iordb,nordia,nordib)
      deallocate (idifa,idifb)
      call timer (4,iconfigs,'_configs  ',-1)
      return
      END
