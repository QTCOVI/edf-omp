      character(len=2)   namatm(0:106)
      real(kind=8)       wgatm(0:106)
      real(kind=8)       rbragg(0:106)
      integer            atminstall
      common   /datatm/  wgatm
      common   /datatc/  namatm
      common   /datati/  atminstall
      common   /radii/   rbragg
