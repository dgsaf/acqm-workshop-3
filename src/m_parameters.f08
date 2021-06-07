!>
module m_parameters

  implicit none

  ! numeric parameters
  double precision , parameter :: TOLERANCE = TINY(1.0d0)
  double precision , parameter :: INFINITY = HUGE(1.0d0)
  double precision , parameter :: PI = 4.0d0*DATAN(1.0d0)

  ! pretty printing parameters
  integer , parameter :: DP_MAX = 4

end module m_parameters
