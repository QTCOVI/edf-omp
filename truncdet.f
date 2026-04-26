      program truncdet
      implicit real*8 (a-h,o-z)
      character*(120) wfnfile
      real*8 c(10000000)
      integer*8 n,m,ndets
      integer maxj(10000000)
      logical noyet
*     write (6,*) 'Input file ?'
      read (5,'(a)') wfnfile
      wfnfile=trim(wfnfile)
*     write (6,*) 'EPSDET ?'
      read (5,*) epsdet
      open (1,file=wfnfile)
      open (2,file='out')
      i=0
      maxj(1)=1
      m=0
      moredets=.true.
      do while (moredets)
 1      i=i+1
        read (1,*,end=2) c(i)
*       if (abs(c(i)).lt.epsdet) exit
        do j=1,maxj(i)
          cc=abs(c(i)*c(j))
          m=m+1
          maxj(i+1)=maxj(i)+1
          if (cc.lt.epsdet) then
            maxj(i+1)=j
            exit
          endif
        enddo
        moredets = abs(c(i)).gt.epsdet
      enddo
 2    continue
      ndets=i-1
      do i=1,ndets
        write (2,6) i,maxj(i)
      enddo
*     write (2,*) sum(maxj(1:ndets)),dble(m)
      write (2,7) ndets,dble(m)/dble(ndets*(ndets+1)/2)
 6    format (1x,'Maximum value of j for i = ',I8,' is ',I8)
 7    format (1x,'Dets and Frac of pairs of dets = ',I8,6x,F17.10)
      end
