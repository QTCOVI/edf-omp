c
c-----------------------------------------------------------------------
c
      subroutine gengen (ig,n,m,maxg,nel,iord,inside,outside,first)
c
c.....Extracts a subset of N elements from a set of M>N natural numbers
c     iord(1), ..., iord(M), non necessarily consecutive storing them in 
c     the set inside[1..N]. The remaining M-N elements are stored in the 
c     set outside[1..M-N]. Each time the routine is called, it produces 
c     new subsets of 'N' and 'M-N' elements. Before 'gengen' is called
c     for the first time the logical variable FIRST(IG) must be set to 
c     .TRUE. This first time the returned inside[1..N] and outside[1..M-N] 
c     sets will contain the first N and the last M-N numbers of the set. 
c     This first time, 'gengen' makes FIRST(IG) equal to .FALSE. One 
c     must be careful that the routine is called as much (M over N) times.
c
      include 'implicit.inc'
      parameter (maxinout=100)
      integer combin(maxinout,maxg),combout(maxinout,maxg)
      integer n,m,i,j,isout,ig
      integer inside(nel),outside(nel)
      integer iord(nel)
      logical first(maxg)
c
      if (n.gt.maxinout.or.m-n.gt.maxinout) then
        stop 'gengen.f: Increase MAXINOUT parameter in gengen.f'
      endif
      if (first(ig)) then
        first(ig)=.false.
        do i=1,n
          combin(i,ig)=i
          inside(i)=iord(i)
        enddo
        do i=n+1,m
          in=i-n
          combout(in,ig)=i
          outside(in)=iord(i)
        enddo
      else
        i=n
        combin(1:n,ig) = inside(1:n)
        combout(1:m-n,ig) = outside(1:m-n)
 1      do while (i.ge.1)
          if (combin(i,ig).lt.m) then
            combin(i,ig)=combin(i,ig)+1
            do j=i+1,n
              combin(j,ig)=combin(j-1,ig)+1
              if (combin(j,ig).gt.m) then
                 i=i-1
                 goto 1
              endif
            enddo
            isout=0
            do j=1,m
              do i=1,n
                if (combin(i,ig).eq.j) goto 2
              enddo
              isout=isout+1
              combout(isout,ig)=j
 2            continue
            enddo
            inside (1:n)  =iord(combin(1:n,ig))
            outside(1:m-n)=iord(combout(1:m-n,ig))
            return
          else
            i=i-1
          endif
        enddo
      endif
      return
      end
