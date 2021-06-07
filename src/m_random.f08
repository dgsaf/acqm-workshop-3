!>
module m_random

  use m_parameters , only : PI

  implicit none

  private
  public :: prepare_random, random_uniform, random_normal

contains

  ! prepare_random
  subroutine prepare_random (seed)
    integer , intent(in) :: seed
    integer :: n
    integer , allocatable :: s(:)

    call random_seed(size=n)

    allocate(s(n))
    s(:) = seed

    call random_seed(put=s)

  end subroutine prepare_random

  ! random_uniform
  subroutine random_uniform (a, b, x)
    double precision , intent(in) :: a, b
    double precision , intent(out) :: x
    double precision :: r

    call random_number(r)

    x = (1 - r)*(b - a) + a

  end subroutine random_uniform

  ! random_normal
  subroutine random_normal (mu, sigma, x1, x2)
    double precision , intent(in) :: mu, sigma
    double precision , intent(out) :: x1, x2
    double precision :: u1, u2, r

    u1 = random_uniform(0.0d0, 1.0d0)
    u2 = random_uniform(0.0d0, 1.0d0)

    r = sigma*sqrt(-2.0d0 * log(u1))

    x1 = mu + r*cos(2.0d0*PI*u2)
    x2 = mu + r*sin(2.0d0*PI*u2)

  end subroutine random_normal

end module m_random
