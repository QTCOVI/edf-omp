c
c-----------------------------------------------------------------------
c
      subroutine rhoalls 
     &  (line,wfnfile,iwfn,ngroup,stdout,aom,nfugrp,ifugrp)
      include     'implicit.inc'
      include     'wfn.inc'
      include     'corr.inc'
      integer, allocatable,dimension (:)     :: rsrs
      character*(80) res
      real(kind=8)  aom(ncent,nmo,nmo)
      integer nfugrp(ngroup)
      integer ifugrp(ncent,ngroup)
      character*(*) line,wfnfile
      integer iwfn,ngroup
      integer stdout,leng
      logical ok,setint
c
c-----------------------------------------------------------------------
c
      if (.not.allocated(rsrs)) allocate (rsrs(ngroup))
      lp=1
      nelread=0
      do i=1,ngroup
        ok=setint(rsrs(i),line,lp)
        if (ok) then
          if (rsrs(i).gt.nel) then
            write (0,*)'# rhoalls.f: Too high population of fragment ',i
            stop
          elseif (rsrs(i).le.0) then
            write (0,*)'# rhoalls.f: Too low population of fragment  ',i
            stop
          else
            nelread=nelread+rsrs(i)
          endif
        else
          stop ' # rhoalls.f: Bad format in the RHOALLS keyword'
        endif
      enddo
      if (nelread.ne.nel) then
        write (0,*) '# rhoalls.f: Improper Number of electrons in [S]'
      endif
      do i=1,ngroup
        do j=1,ngroup
          icon=rsrs(i)
          if (j.eq.i) icon=icon-1
          write (res(4*(j-1)+1:4*j),'(I4)') icon
        enddo
        call rhocond (res,wfnfile,iwfn,ngroup,stdout,aom,nfugrp,ifugrp)
      enddo
      if (allocated(rsrs)   ) deallocate (rsrs)
      return
      end

