c
c.......................................................................
c
c     This routine recovers from the input file the variables that are
c     relevant to the localization performed in the 'dafhdo.f' routine.
c
c-----------------------------------------------------------------------
c
      subroutine dafhdrv (lr,lw,mxdafh)
c
c.......................................................................
c
      USE          space_for_dafh
      include     'implicit.inc'
      include     'wfn.inc'
      include     'mline.inc'
*     parameter   (mline = 200)  
      integer(kind=4)  lw,lr,leng
      logical          setint,ok,goon,setword,setdble
      character*(mline) line,word,uppcase
      character*1 dp
c
c-----------------------------------------------------------------------
c
      ndafh=0
      goon=.true.
 20   do while (goon) 
        read (lr,'(a)',end=110) line
        line=uppcase(line)
        lp=1
        lpa=lp
        ok=setword(word,line,lp)
        lword=leng(word)
        if (word(1:1).eq.'#'.or.lword.eq.0) then
        elseif (word(1:lword).eq.'END'.or.word(1:3).eq.'END') then
          goon=.false.
        else
          lp=lpa
          ok=setint(n,line,lp)
          if (ok) then
            if (n.gt.0.and.n.le.ncent) then
              ndafh=ndafh+1
              if (ndafh.gt.mxdafh) then
                write (0,*) '# dafhdrv.f: Increase MXDAFH parameter'
                stop
              endif
              nbcen(ndafh)=n
            else
              cycle
            endif
            ok=setdble(cutoff(ndafh),line,lp)
            if (ok) then
              cutoff(ndafh)=abs(cutoff(ndafh))
              if (cutoff(ndafh).ge.1d0) cycle
            else
              ndafh=ndafh-1
              cycle
            endif
            do i=1,n
              ok=setint(ibcen(ndafh,i),line,lp)
              if (.not.ok) cycle
            enddo
            ok=setword(word,line,lp)
            if (ok) then
              if (trim(uppcase(word)).eq.'DEPLET') then
                deplet(ndafh)=.true.
              else
                deplet(ndafh)=.false.
              endif
            else
              deplet(ndafh)=.false.
            endif
          else
            cycle
          endif
        endif
      enddo
 110  continue
      return
      end
