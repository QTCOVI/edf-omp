c
c-----------------------------------------------------------------------
c
      character (len=4) function fourchar (i)
      integer (kind=4) i,i1,i2,i3,i4,i5,i6
      character*1    digs(0:9),lbl
      data digs / '0','1','2','3','4','5','6','7','8','9'/

      i1 = i/1000
      i2 = mod(i,1000)
      i3 = i2/100
      i4 = mod(i2,100)
      i5 = i4/10
      i6 = mod(i4,10)
      fourchar = digs(i1)//digs(i3)//digs(i5)//digs(i6)
      return
      end
