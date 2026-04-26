c
c-----------------------------------------------------------------------
c
      subroutine aomindef (aom,nat,iwrite,lw,wfnfile)
c
c-----------------------------------------------------------------------
c
c     Determine the atomic overlap matrix for MinDef atoms.
c
c-----------------------------------------------------------------------
c
      USE        space_for_wfnbasis
      USE        space_for_wfncoef
      include   'implicit.inc'      
      include   'wfn.inc'      
      include   'param.inc'      
      include   'constants.inc'      
      include   'mline.inc'
*     parameter  (mline = 200)  
      character*(*) wfnfile
      character*(mline) mindefile
      real(kind=8)    aom(ncent,nmo,nmo)
      integer    it(3), jt(3)
      real(kind=8)     ax(3), bx(3), ori, orj, xyzint
      real(kind=8)     sabxyz
      real(kind=8), allocatable,dimension (:,:)  :: aomt
      logical iwrite
c
c.....Overlap matrix between Cartesian Gaussian primitives.
c
      if (.not.sprimok) then
        call gtogto ()
        sprimok=.true.
      endif
c
c.....Initialize aom()
c
      aom(1:ncent,1:nmo,1:nmo)=zero
c
c.....first loop over primitives
c
      if (nat .eq. 1) then
        nmois = 0
      else
        nmois = nmo
      endif
      do i=1,nprims
        ic=icen(i)
        itip=ityp(i)
        it(1:3)=nlm(itip,1:3)
        ori=oexp(i)
c
c.......second loop over primitives
c
        do j=1,nprims
          jc=icen(j)
          jtip=ityp(j)
          jt(1:3)=nlm(jtip,1:3)
          orj=oexp(j)
          ax(1:3)=xyz(ic,1:3)
          bx(1:3)=xyz(jc,1:3)
c
c.........Overlap integral between Cartesian Gaussian primitives.
c
*         xyzint=sabxyz(it,jt,ax,bx,ori,orj)    ! Obara-Saika
          xyzint=sprim(i,j)
          do io=1,nmo
            do jo=1,nmo
              wab=coef(io+nmois,i)*coef(jo+nmois,j)
              wabxyz = wab * xyzint
              if (ic.eq.jc) then
                aom(ic,io,jo) = aom(ic,io,jo) + wabxyz
              else
                if (ori.gt.orj) then
                  aom(ic,io,jo) = aom(ic,io,jo) + wabxyz
                elseif (ori.lt.orj) then
                  aom(jc,io,jo) = aom(jc,io,jo) + wabxyz
                else
                  aom(ic,io,jo) = aom(ic,io,jo) + wabxyz * half
                  aom(jc,io,jo) = aom(jc,io,jo) + wabxyz * half
                endif
              endif
            enddo
          enddo
        enddo
      enddo      
c
c.....Write results to an AOM file called wfnfile'.mindefAOM'
c
      mindefile=wfnfile(1:leng(wfnfile))//'.mindefAOM'
      open (18,file=mindefile,iostat=ierr)
      if (ierr.ne.0) then
        write (0,*) 'aomindef.f: Error opening *.mindefAOM file'
        return
      endif
      do i=1,ncent
        write (18,'(I5,a)') i,' <--- AOM within this center '
        write (18,80) ((aom(i,m,j),m=1,j),j=1,nmo)
      enddo
      close (18)
c
c-----Print AOM
c
      if (iwrite) then
        if (.not.allocated(aomt)) then
          allocate (aomt(nmo,nmo),stat=ier)
          if (ier.ne.0) stop 'aomindef.f: Cannot allocate aomt()'
        endif
        do m=1,nmo
          do n=1,nmo
            aomt (n,m) = sum(aom(1:ncent,m,n))
          enddo
        enddo
        write (lw,501)
        write (lw,*) '# Results of MinDef integrations'
        write (lw,*) '#'
        write (lw,*) '# AOM'
        do i=1,ncent
          write (lw,*) '# Center ',i
          do j=1,nmo
            write (lw,10) (aom(i,j,k),k=1,j)
          enddo
        enddo
        write (lw,*) '#'
        write (lw,*) '# AOMT'
        do j=1,nmo
          write (lw,10) (aomt(j,k),k=1,j)
        enddo
      endif
      if (iwrite) then
        if (allocated(aomt)) then
          deallocate (aomt,stat=ier)
          if (ier.ne.0) stop 'aomindef.f: Cannot deallocate aomt()'
        endif
      endif
 10   format (6(1x,F15.8))
 80   format (6(1x,e16.10))
 501  format (1x,69('-'))
      return
      end
