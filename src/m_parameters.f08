!>
module m_parameters

#ifdef f2003
  use , intrinsic :: iso_fortran_env , only : input_unit, output_unit, &
      error_unit
#else
#define input_unit 5
#define output_unit 6
#define error_unit 0
#endif

  implicit none

  ! file units
  integer , parameter :: STDIN = input_unit
  integer , parameter :: STDOUT = output_unit
  integer , parameter :: STDERR = error_unit

  ! pretty printing parameters
  integer , parameter :: DP_MAX = 4

  ! numeric parameters
  double precision , parameter :: NEAR_ZERO = TINY(1.0d0)
  double precision , parameter :: INFINITY = HUGE(1.0d0)
  double precision , parameter :: PI = 4.0d0*DATAN(1.0d0)
  double precision , parameter :: TOLERANCE = 1.0d-3
  double precision , parameter :: SMALL = 1.0d-4

  ! computational parameters
  integer , parameter :: MAX_ITERATIONS = 1000

end module m_parameters
