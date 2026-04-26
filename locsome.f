c
c-----------------------------------------------------------------------
c
      subroutine locsome (aom,line,lw,largwr)
c 
c.......................................................................
c
      USE      space_for_wfnbasis
      USE      space_for_wfncoef
      include  'implicit.inc'
      include  'param.inc'
      include  'constants.inc'
      include  'wfn.inc'
      include  'corr.inc'
      include  'stderr.inc'

      real(kind=8), allocatable,dimension (:)      :: oc
      real(kind=8), allocatable,dimension (:,:)    :: vrot
      real(kind=8), allocatable,dimension (:,:,:)  :: sgallc
      real(kind=8), allocatable,dimension (:,:)    :: aomt
      integer,      allocatable,dimension (:)      :: iloc,ilocx

      character*(*) line
      real(kind=8)       aom(ncent,nmo,nmo)
      logical       setint,ok,largwr,notinlist
      integer       lp,nloc,iaddloc
c
      if (.not.allocated(iloc)) then
        allocate (iloc(nmo),stat=ier)
        if (ier.ne.0) stop 'locsome.f: Cannot allocate iloc()'
      endif
      lp=1
      nloc=0
      do while (setint(n,line,lp))
        if (n.le.nmo.and.n.gt.0) then
          if (nloc.eq.0) then
            nloc = 1
            iloc(1) = n
          else
            i=1
            notinlist=.true.
            do while (i.le.nloc.and.notinlist)
              if (n.eq.iloc(i)) notinlist=.false.
              i=i+1
            enddo
            if (notinlist) then
              nloc=nloc+1
              iloc(nloc)=n
            endif
          endif
        endif
      enddo
c
c.....There are no MOs to be localized
c
      if (nloc.eq.0) then
        if (allocated(iloc)) then
          deallocate (iloc,stat=ier)
          if (ier.ne.0) stop 'locsome.f: Cannot deallocate iloc()'
        endif
        return
      endif

      if (.not.allocated(aomt)) then
        allocate (aomt(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'locsome.f: Cannot allocate aomt()'
      endif
      write (lw,99) (iloc(i),i=1,nloc)
      if (largwr) then
        write (lw,100) 'BEFORE'
        aomt(1:nmo,1:nmo)=zero
        do igr=1,ncent
          aomt(1:nmo,1:nmo)=aomt(1:nmo,1:nmo)+aom(igr,1:nmo,1:nmo)
          write (lw,101) igr
          do i=1,nmo
            write (lw,444) (aom(igr,i,j),j=1,i)
          enddo
        enddo
        write (lw,*) 'SUM OVER CENTERS'
        do i=1,nmo
          write (lw,444) (aomt(i,j),j=1,i)
        enddo
      endif

      if (.not.allocated(oc)) then
        allocate (oc(nloc),stat=ier)
        if (ier.ne.0) stop 'locsome.f: Cannot allocate oc()'
      endif
      if (.not.allocated(vrot)) then
        allocate (vrot(nloc,nloc),stat=ier)
        if (ier.ne.0) stop 'locsome.f: Cannot allocate vrot()'
      endif
      if (.not.allocated(sgallc)) then
        allocate (sgallc(ncent,nloc,nloc),stat=ier)
        if (ier.ne.0) stop 'locsome.f: Cannot allocate sgallc()'
      endif
      sgallc=zero
      do igr=1,ncent
        do i=1,nloc
          oc(i)=two
          do k=1,nloc
            sgallc(igr,i,k)=aom(igr,iloc(i),iloc(k))
          enddo
        enddo
      enddo
      call ruedmis (sgallc,vrot,oc,nloc,ncent,lw,largwr,.false.)
      do igr=1,ncent
        do i=1,nloc
          do k=1,nloc
            aom(igr,iloc(i),iloc(k))=sgallc(igr,i,k)
          enddo
        enddo
      enddo
c
      if (largwr) then
        write (lw,100) 'AFTER'
        aomt(1:nmo,1:nmo)=zero
        do igr=1,ncent
          aomt(1:nmo,1:nmo)=aomt(1:nmo,1:nmo)+aom(igr,1:nmo,1:nmo)
          write (lw,101) igr
          do i=1,nmo
            write (lw,444) (aom(igr,i,j),j=1,i)
          enddo
        enddo
        write (lw,*) 'SUM OVER CENTERS'
        do i=1,nmo
          write (lw,444) (aomt(i,j),j=1,i)
        enddo
        write (lw,*) 
        write (lw,*)
      endif

      if (allocated(iloc)) then
        deallocate (iloc,stat=ier)
        if (ier.ne.0) stop 'locsome.f: Cannot deallocate iloc()'
      endif
      if (allocated(aomt)) then
        deallocate (aomt,stat=ier)
        if (ier.ne.0) stop 'locsome.f: Cannot deallocate aomt()'
      endif
      if (allocated(oc)) then
        deallocate (oc,stat=ier)
        if (ier.ne.0) stop 'locsome.f: Cannot deallocate oc()'
      endif
      if (allocated(vrot)) then
        deallocate (vrot,stat=ier)
        if (ier.ne.0) stop 'locsome.f: Cannot deallocate vrot()'
      endif
      if (allocated(sgallc)) then
        deallocate (sgallc,stat=ier)
        if (ier.ne.0) stop 'locsome.f: Cannot deallocate sgallc()'
      endif
      return

 99   format (//,' # LOCALIZATION OF A SUBSET OF MOs:',20I3)
 100  format (1x,'# AOM ',A,' LOCALIZATION')
 101  format (1x,'# CENTER ',I2)
 444  format (6(2x,F12.6))

      end
