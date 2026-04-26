
c-----------------------------------------------------------------------
c
      subroutine srhon (epsdet,ifilc,ndets,ncore,nelact,nmo)
c
c.....Obtain the first-order density matrix coefficients.
c.......................................................................
c
      USE              space_for_cidet
      USE              space_for_rdm1
      include         'implicit.inc'
      include         'constants.inc'
      include         'param.inc'
      include         'lengrec.inc'
      integer          ifilc,ndets,ncore,nelact,nmo
      real(kind=8)     epsdet,indet
c
      integer,         allocatable,dimension(:)      :: idif
      integer,         allocatable,dimension(:)      :: ordia,ordib
      real(kind=8),    allocatable,dimension(:)      :: cd
      integer,         allocatable,dimension(:)      :: sob,sod,iord
      integer,         allocatable,dimension(:)      :: sodaux
      integer(kind=8), allocatable, dimension(:,:,:) :: det
c
      call timer (2,ipid,'_srhon    ',-1)
      if (.not.allocated(ordia)) then
        allocate (ordia(nmo),stat=ier)
        if (ier.ne.0) stop 'srhon.f: Cannot allocate ordia()'
      endif
      if (.not.allocated(ordib)) then
        allocate (ordib(nmo),stat=ier)
        if (ier.ne.0) stop 'srhon.f: Cannot allocate ordib()'
      endif
c
c.....Determine the dimension which is necessary for some arrays.
c
      read (ifilc,rec=1) (cidet(k),k=0,nelact)
      mal=0
      mbe=0
      do k=1,nelact
        mm=int(cidet(k))
        if (mm.gt.0) mal=mal+1
        if (mm.lt.0) mbe=mbe+1
      enddo
      c1ea = zero
      c1eb = zero
!
!-----INTEGER*8 numbers necessary to codify an ALPHA or BETA config
!
      nmobin = (nmo/64 + 1)
      if (.not.allocated(det)) then
        allocate (det(nmobin,2,ndets),stat=ier)
        if (ier.ne.0) stop 'binrdm.f: Cannot allocate det()'
      endif
      if (.not.allocated(cd)) then
        allocate (cd(ndets),stat=ier)
        if (ier.ne.0) stop 'binrdm.f: Cannot allocate cd()'
      endif
      
*     luc0=81
*     open (luc0,file='C0coef',access='direct',
*    &      recl=RecLength*ncdw,form='unformatted')
      if (ncore.gt.0) then
        do mm=1,ncore
          ordia(mm)=mm
          ordib(mm)=mm
        enddo
      endif
      if (.not.allocated(sob)) then
        allocate ( sob(nelact),stat=ier)
        if (ier.ne.0) stop 'srhon.f: Cannot allocate sob()'
      endif
      if (.not.allocated(sod)) then
        allocate ( sod(nelact),stat=ier)
        if (ier.ne.0) stop 'srhon.f: Cannot allocate sod()'
      endif
      if (.not.allocated(sodaux)) then
        allocate ( sodaux(nelact),stat=ier)
        if (ier.ne.0) stop 'srhon.f: Cannot allocate sodaux()'
      endif
      if (.not.allocated(iord)) then
        allocate (iord(nelact),stat=ier)
        if (ier.ne.0) stop 'srhon.f: Cannot allocate iord()'
      endif
      if (.not.allocated(idif)) then
        allocate (idif(nelact),stat=ier)
        if (ier.ne.0) stop 'srhon.f: Cannot allocate idif()'
      endif
      iper     = 0
      shundred = dble(ndets)/10.0d0
      xhundred = shundred
      indet    = 0.0d+00
      det      = 0_8
      xnorm    = 0.0d+00
      do i=1,ndets
        indet=indet+1.0d+00
        if (indet.ge.xhundred.and.ndets.gt.1000) then
          xhundred=xhundred+shundred
          iper=iper+10
          write (0,1212) iper
          call flush (0)
        endif
        read (ifilc,rec=i) (cidet(k),k=0,nelact)
        cd(i) = cidet(0)
        xnorm   = xnorm + cd(i)**2
        do m=1,nelact
          mm=int(cidet(m))
          iord(m)=m
          if (mm.gt.0) then
            sob(m)=mm
          else
            sob(m)=mm+9999
          endif
        enddo
        call iqcksort (sob,iord,nelact,1,nelact)
        do m=1,mal
          sod(m)=sob(iord(m))
        enddo
        nn=0
        do m=mal+1,nelact
          iom=iord(m)
          sob(iom)=sob(iom)-9999
          sod(nelact-nn)=sob(iom)
          nn=nn+1
        enddo
        sodaux(1:nelact)=sod(1:nelact)
        call colord (sob,sodaux,nelact,nelact,nper,nsig,ndif,idif)
        if (ndif.gt.0) stop 'binrdm.f: NDIF should be 0 at this point'
        cd(i) = cd(i) * nsig
        ordia(ncore+1:ncore+mal)=ncore+sod(1:mal)
        ordib(ncore+1:ncore+mbe)=ncore-sod(mal+1:nelact)
        do ia=1,ncore+mal
          ibin=(ordia(ia)/64 + 1)
          det(ibin,1,i)=ibset(det(ibin,1,i),mod(ordia(ia)-1,63))
        enddo
        do ib=1,ncore+mbe
          ibin=(ordib(ib)/64 + 1)
          det(ibin,2,i)=ibset(det(ibin,2,i),mod(ordib(ib)-1,63))
        enddo
      enddo
!
!.....Renormalize the WFN
!
      xnorm       = 1.0d+00/sqrt(xnorm)
      cd(1:ndets) = cd(1:ndets) * xnorm
!
!.....Compute 1-RDM and 2-RDM
!
      call computerdm (det,ndets,cd,epsdet,nmo,nmobin,c1ea,c1eb)
      c1et=c1ea+c1eb
      call timer (4,ipid,'_srhon    ',-1)
      return
 1212 format (1x,'# Codifying determinants :',i4,'% done')
      end
