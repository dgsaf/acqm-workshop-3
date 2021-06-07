!>
module m_integrate

  implicit none

contains

  function integrate_trapezoid (n, x, f) result (integral)
    integer , intent(in) :: n
    double precision , intent(in) :: x(n), f(n)
    double precision :: w(n)
    double precision :: integral
    integer :: ii
    integer :: status

    ! check if arguments are valid
    status = 0

    if (n < 2) then
      status = 1
    end if

    ! terminate subroutine if arguments are invalid, otherwise, proceed
    if (status /= 0) then
      integral = 0.0d0
      return
    end if

    ! calculate weights
    w(1) = (x(2) - x(1)) / 2.0d0

    do ii = 2, n-1
      w(ii) = (x(ii+1) - x(ii-1)) / 2.0d0
    end do

    w(n) = (x(n) - x(n-1)) / 2.0d0

    ! calculate integral
    integral = 0.0d0

    do ii = 1, n
      integral = integral + (w(ii) * f(ii))
    end do

  end function integrate_trapezoid

end module m_integrate
