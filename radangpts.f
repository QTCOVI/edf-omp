c
c-----Radial points
c
      integer(kind=4) function z2nr(z)
    
      implicit none
      integer(kind=4), intent(in) :: z
    
      z2nr = 40
      if (z > 2)  z2nr = 60
      if (z > 10) z2nr = 80
      if (z > 18) z2nr = 100
      if (z > 36) z2nr = 120
      if (z > 54) z2nr = 140
      if (z > 86) z2nr = 160
    
      endfunction
    
c
c-----Angular points
c
      integer(kind=4) function z2nang(z)
    
      implicit none
      integer(kind=4), intent(in) :: z
  
      z2nang = 194
    
      endfunction
