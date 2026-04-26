module mod_param

  use mod_prec, only: real

  implicit none

  private
  
  ! physical constants
  real*8, parameter, public :: pi      = 3.14159265358979323846
  real*8, parameter, public :: rad     = pi / 180                 
  real*8, parameter, public :: cte     = 2.71828182845904523536
  real*8, parameter, public :: ctsq2   = 1.41421356237309504880
  real*8, parameter, public :: ctsq3   = 1.73205080756887729352
  real*8, parameter, public :: cteuler = 0.57721566490153286061
  real*8, parameter, public :: ctgold  = 1.61803398874989484820
  real*8, parameter, public :: bohr2a  = 0.5291771
  real*8, parameter, public :: bohr2pm = 52.91771
  real*8, parameter, public :: har2cm1 = 219474.6
  real*8, parameter, public :: har2ev  = 27.212
  real*8, parameter, public :: kb      = 1.3806488d-23
  real*8, parameter, public :: na      = 6.02214129d23
  real*8, parameter, public :: r       = kb*na        
  real*8, parameter, public :: h       = 6.626d-34
  real*8, parameter, public :: c       = 299.792458d8
  real*8, parameter, public :: tokcalmol = 1.0/4184.0


  ! parameters
  real*8, parameter, public :: zero    = 0
  real*8, parameter, public :: one     = 1
  real*8, parameter, public :: two     = 2
  real*8, parameter, public :: three   = 3
  real*8, parameter, public :: half    = 0.5
  real*8, parameter, public :: third   = one / three
  real*8, parameter, public :: fourth  = 0.25
  real*8, parameter, public :: undef   = 9.72 + 21
  real*8, parameter, public :: EPS     = 1E-8
  real*8, parameter, public :: pusheps = 1E-14

end module mod_param
