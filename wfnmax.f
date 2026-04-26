c
c-----Driver routine for the location of the maxima of Psi^2
c
      subroutine wfnmax (lr,lw,lerr)
c
      include 'implicit.inc'
      include 'param.inc'
      include 'constants.inc'
      include 'opt.inc'
      include 'corr.inc'
      include 'mline.inc'
*     parameter  (mline=200)
      parameter  (msave=7)
      character*(mline) wfnf
      logical setword,setint,setdble,ok
      character*(mline) line,word,uppcase
      real(kind=8),  allocatable,dimension (:)     :: xvar,gpsi
      logical  warn,goon
      external mxwfn
c
c-----------------------------------------------------------------------
c
      call timer (2,iwfnmax,'_wfnmax   ',-1)
c
c-----Read bonding pattern data from the input file.
c
      rewind (lr)
      lwrite = lw
c
c-----Default weights in the optimization function
c
      mup    = 5
      iqpr   = 2
      gtol   = 1D-6
      maxit  = 500
      goon   = .true.
 20   do while (goon) 
        read (lr,'(a)',end=110) line
        line=uppcase(line)
        lp=1
        ok=setword(word,line,lp)
        if (word(1:1).eq.'#') then
        elseif (word(1:3).eq.'END') then
          goon = .false.
        elseif (word(1:leng(word)).eq.'GTOL') then
          if (setdble(gtol,line,lp)) gtol=abs(gtol)
        elseif (word(1:leng(word)).eq.'MAXIT') then
          if (setint(maxit,line,lp)) maxit=abs(maxit)
        elseif (word(1:leng(word)).eq.'MUP') then
          if (setint(mup,line,lp)) mup=min(iabs(mup),msave)
        elseif (word(1:leng(word)).eq.'IPRINT') then
          if (setint(iprint,line,lp)) iqpr=iprint
        elseif (word(1:leng(word)).eq.'CONFIG') then
          n=0
 21       read (lr,'(a)') line
          line=uppcase(line)
          lp=1
          ok=setword(word,line,lp)
          if (word(1:3).eq.'END') then
            if (n.lt.nel) then
              stop 'wfnmax.f: Few electron coordinates are read in'
            endif
          else
            lp=1
          endif
          ok=.true.
          ok=ok.and.setdble(xx,line,lp)
          ok=ok.and.setdble(yy,line,lp)
          ok=ok.and.setdble(zz,line,lp)
          if (ok) then
            n=n+1
            xyzel(n,1)=xx
            xyzel(n,2)=yy
            xyzel(n,3)=zz
          else
            stop 'wfnmax.f: Wrong line format in the CONFIG environment'
          endif
          iopt(n,1:3)=1
          if (.not.setint(io,line,lp)) then
            iopt(n,1:3)=1
          else
            if (io.eq.0) then
              iopt(n,1)=0
              if (.not.setint(io,line,lp)) then
                iopt(n,2:3)=1
              else
                if (io.eq.0) iopt(n,2)=0
                if (io.ne.0) iopt(n,2)=1
                if (.not.setint(io,line,lp)) then
                  iopt(n,3)=1
                else
                  if (io.eq.0) iopt(n,3)=0
                  if (io.ne.0) iopt(n,3)=1
                endif
              endif
            else
              iopt(n,1)=1
              if (.not.setint(io,line,lp)) then
                iopt(n,2:3)=1
              else
                if (io.eq.0) iopt(n,2)=0
                if (io.ne.0) iopt(n,2)=1
                if (.not.setint(io,line,lp)) then
                  iopt(n,3)=1
                else
                  if (io.eq.0) iopt(n,3)=0
                  if (io.ne.0) iopt(n,3)=1
                endif
              endif
            endif
          endif
          if (n.lt.nel) goto 21
        endif
      enddo
 110  continue
c
c-----Optimizable variables are set up
c
      nvar=0
      do i=1,n
        do j=1,3
          if (iopt(i,j).eq.1) nvar=nvar+1
        enddo
      enddo
      allocate (xvar(nvar))
      allocate (gpsi(nvar))
      k=0
      do i=1,n
        do j=1,3
          if (iopt(i,j).eq.1) then
            k=k+1
            xvar(k)=xyzel(i,j)
          endif 
        enddo
      enddo

      call fdlbfgs (nvar,mup,xvar,psi,gpsi,mxwfn,iqpr,
     &    gtol,lw,lerr,maxit,iter,warn)
      call timer (4,iwfnmax,'_wfnmax   ',-1)
      return
      end
