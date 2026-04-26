      subroutine        symgr(mask,maxis,linear,symbol)
c
c.....Finds the appropiate point group symbol by means of the
c     previous operation of pointgroup and typeop.
c
      include           'implicit.inc'
c
c.....Local variables
      character*6       symbol
      integer           mask(-1:1,maxis)
      integer           highest,highests
      logical           linear
c
c.....We work mainly with matrix mask.
c
      if (linear) then
         if (mask(-1,2).eq.0) then
             symbol='Cinfv'
         else
             symbol='Dinfh'
         endif
         return
      endif
c.....Is it a cubic point group?
      if (mask(1,3) .gt.2) then
c
c........cubic
         if (mask(1,4).ne.0) then
c...........O
            if (mask(-1,2).ne.0) then
c..............Oh
               symbol='Oh'
            else
c..............O
               symbol='O '
            endif
         elseif (mask(-1,4).ne.0) then
c...........Td
            symbol= 'Td'
         elseif (mask(-1,2).ne.0) then
c...........Th
            symbol = 'Th'
         else
c...........T
            symbol = 'T '
         endif
      else

c........Compute highest order proper axis.
         highest=0
         highests=0
         if (mask(1,2).ne.0) highest=2
         do i=3,12
            if (mask( 1,i).ne.0) highest=i
            if (mask(-1,i).ne.0) highests=i
         enddo
         if (highest.eq.0) then
            if (mask(-1,2).ne.0) then
c..............i
               symbol='i '
            elseif (mask(-1,1).ne.0) then
c..............Cs
               symbol='Cs'   
            else
c..............C1
               symbol='C1' 
            endif
         elseif (mask(1,2) .ge. highest) then
            if (mask(-1,1).eq.0) then 
c..............Dn
               symbol='D'//char(ichar('0')+highest)
            elseif (mask(-1,1) .eq. highest) then
c.................Dnd
                  symbol= 'D'//char(ichar('0')+highest)//'d'
            else
c.................Dnh
                  symbol= 'D'//char(ichar('0')+highest)//'h'
            endif
         elseif (mask(-1,1).eq.0) then
            if (highests .ne. 0) then
c..............Sn
               symbol= 'S'//char(ichar('0')+highests/2)
            else
c..............Cn
               symbol= 'C'//char(ichar('0')+highest)
            endif
         elseif (mask(-1,1) .lt. highest) then
c...........Cnh
               symbol= 'C'//char(ichar('0')+highest)//'h'
         else
c...........Cnv
               symbol= 'C'//char(ichar('0')+highest)//'v'
         endif
      endif
      end
