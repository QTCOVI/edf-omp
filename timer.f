c
c-----------------------------------------------------------------------
c
      subroutine timer (key,pid,name,lw)
c
c.....timer - accumulates and prints out the elapsed times of a series
c             of up to MPROCESS processes.
c
c     Input parameters are:
c     key ........... 0 = reset all time tables.
c                     1 = reset pid entry and begin counting for it.
c                     2 = continue the cont for pid process. Do not
c                         reset previous entry times.
c                     3 = end of pid process. Free pid entry.
c                     4 = end of pid process. Do not free entry.
c                     5 = end of run. Print out all time tables.
c                     6 = close all processes and print out tables.
c     pid ........... process identification number (1..MPROCESS).
c     name .......... process name (used only in the print out).
c     lw ............ printer logical unit. Negative if print out is not
c                     desired.
c
c.....key controls what data are used in the run, as the following table
c     resumes:
c
c     key value    pid      name       lw
c     ---------  -------   -------   -------
c        0       ignored   ignored   ignored
c       1,2      output    input     ignored
c       3,4      input     ignored   input
c       5,6      ignored   ignored   input
c
      integer          key,pid,lw
      character*(*)    name
c
      include  'error.inc'
c
      integer          mprocess
      parameter        (mprocess=100)
      character*(15)   pname(mprocess)
      real*4           time(mprocess),cumtime(mprocess),timedum
      logical          popen(mprocess),pused(mprocess)
      logical          firsttime
      integer          pcalls(mprocess)
      save             time,cumtime,pname,popen,pused,pcalls
      save             firsttime
c
      data firsttime /.true./
c
      if (key.eq.0 .or. firsttime) then
c
c........initiallize all entries:
c
         firsttime=.false.
         do 10 i=1,mprocess
            time(i)=0.0
            cumtime(i)=0.0
            popen(i)=.false.
            pused(i)=.false.
            pname(i)=' '
            pcalls(i)=0
 10      continue
      else if (key.eq.1 .or. key.eq.2) then
c
c........begin pid count:
c
         call cpu_time (timedum)
         i=1
 15      if (pused(i)) then
            if (pname(i)(1:10).ne.name(1:10)) then
               i=i+1
               if (i.gt.mprocess) then
                  call error ('timer','pid out of bonds',warning)
                  return
               endif
               goto 15
            endif
         endif
         pid=i
         if (key.eq.1) then
            cumtime(pid)=0.0
            pcalls(pid)=1
         else
            pcalls(pid)=pcalls(pid)+1
         endif
         time(pid)=timedum
         popen(pid)=.true.
         pused(pid)=.true.
         pname(pid)=name
      else if (key.eq.3 .or. key.eq.4) then
c
c........end pid accounting:
c
         if (pid.le.0 .or. pid.gt.mprocess) then
            call error ('timer','pid out of bonds',warning)
         else if (.not.popen(pid)) then
            call error ('timer','pid unused or closed',warning)
         else
            call cpu_time (timedum)
            time(pid)=timedum-time(pid)
            cumtime(pid)=cumtime(pid)+time(pid)
            popen(pid)=.false.
            if (lw.gt.0) write (lw,100) pname(pid),time(pid)
            if (key.eq.3) then
               pused(pid)=.false.
               pcalls(pid)=0
               cumtime(pid)=0.0
               time(pid)=0.0
            endif
         endif
      else if (key.eq.5 .or. key.eq.6) then
c
c........print out the time tables:
c
         write (lw,105)
         call cpu_time (timedum)
         do 25 i=1,mprocess
            if (pused(i)) then
               if (popen(i)) then
                  time(i)=timedum-time(i)
                  cumtime(i)=cumtime(i)+time(i)
                  if (key.eq.6) popen(i)=.false.
               endif
               write (lw,110) i,pname(i)(1:10),
     &                  cumtime(i),pcalls(i),popen(i)
            endif
 25      continue
         write (lw,115)
      else
         call error ('timer','key value not recognized',warning)
      endif
      return
c
 100  format (/' #------------------------------------',/'#*** timer:'/
     & ' #    process name:',a10,5x,'elapsed time (sec):',f10.3/)
 105  format (/' #',53('-'),/,' #    timer:'/,' #   '/
     & ' #    -pid----name--------------cumtime-------pcalls--popen-')
 110  format (' #   ',i3,5x,a10,3x,f14.6,1x,i10,4x,l2)
 115  format (' #   '/)
      end
