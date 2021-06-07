!>
module m_diffeq

  use m_parameters, only : INFINITY

  implicit none

  private
  public :: numerov_f, numerov_b

contains

  ! numerov_f
  !
  ! Brief:
  ! Calculates the values of $y(x)$, defined by the differential equation
  ! $$ y''(x) = -g(x)y(x) + s(x) $$ ,
  ! on a grid to fourth-order accuracy using the forward Numerov method.
  !
  ! Summary:
  ! Suppose, $X = \{x_{1}, \dotsc, x_{n_{x}}\}$, is a grid with $n_{x}$ points,
  ! such that $x_{i+1} - x_{i} = h$ for all $i = 1, \dotsc, n_{x}-1$,
  ! and that $g(x), s(x) : \mathbb{R} \to \mathbb{R}$ have been evaluated on
  ! this grid.
  ! Suppose further that $y : \mathbb{R} \to \mathbb{R}$ is defined by the
  ! differential equation,
  ! $$ y''(x) = -g(x)y(x) + s(x) $$ ,
  ! and the values of $y(x_{1}), y(x_{2})$ are known.
  ! The values of $y(x_{i})$ for $i = 3, \dotsc, n_{x}$ are calculated to
  ! fourth-order accuracy using the forward Numerov method.
  !
  ! Input:
  ! - `n_x` is the number of grid points.
  ! - `step_size` is the distance between consecutive points on the grid $X$.
  ! - `s_grid` is the evaluation of $s(x)$ on the grid $X$.
  ! - `g_grid` is the evaluation of $g(x)$ on the grid $X$.
  !
  ! Output:
  ! - `y_grid` is the evaluation of $y(x)$ on the grid $X$, calculated via the
  !   forward Numerov method.
  ! - `status` integer status code which takes the following values:
  !   - `status == 0` indicates successful execution;
  !   - `status > 0` indicates that the arguments were invalid;
  !   - `status == -1` indicates that a numerical error (NaN or infinity)
  !     occured during execution.
  !
  subroutine numerov_f (n_x, step_size, s_grid, g_grid, y_grid, status)
    integer , intent(in) :: n_x
    double precision , intent(in) :: step_size
    double precision , intent(in) :: s_grid(n_x)
    double precision , intent(in) :: g_grid(n_x)
    double precision , intent(out) :: y_grid(n_x)
    integer , intent(out) :: status
    double precision :: step_s_grid(n_x)
    double precision :: step_g_grid(n_x)
    integer :: ii

    ! check if arguments are valid
    status = 0

    if (n_x < 1) then
      ! no grid points
      status = 1
    end if

    if (step_size < 0.0d0) then
      ! non-positive step_size
      status = 2
    end if

    ! terminate subroutine if arguments are invalid, otherwise proceed
    if (status /= 0) then
      y_grid(:) = 0.0d0
      return
    end if

    ! perform forward-numerov method
    step_s_grid(:) = s_grid(:) * (step_size ** 2) / 12.0d0
    step_g_grid(:) = g_grid(:) * (step_size ** 2) / 12.0d0

    ! assume that y_grid(1), y_grid(2) already set
    ! y_grid(1) = y_1

    ! if (n_x >= 2) then
    !   y_grid(2) = y_2
    ! end if

    if (n_x >= 3) then
      do ii = 2, n_x - 1
        y_grid(ii+1) = &
            ((2*y_grid(ii)*(1.0d0 - 5.0d0*step_g_grid(ii))) &
            - (y_grid(ii-1)*(1.0d0 + step_g_grid(ii-1))) &
            + (step_s_grid(ii+1) + 10.0d0*step_s_grid(ii) + step_s_grid(ii-1)) &
            ) / (1.0d0 + step_g_grid(ii+1))

        ! check for and handle numerical error (NaN or Infinity)
        if ((y_grid(ii+1) /= y_grid(ii+1)) .or. (y_grid(ii+1) > INFINITY)) then
          status = -1

          if (ii + 2 <= n_x) then
            y_grid(ii+2:n_x) = 0.0d0
          end if

          return
        end if
      end do
    end if

  end subroutine numerov_f

  ! numerov_b
  !
  ! Brief:
  ! Calculates the values of $y(x)$, defined by the differential equation
  ! $$ y''(x) = -g(x)y(x) + s(x) $$ ,
  ! on a grid to fourth-order accuracy using the backward Numerov method.
  !
  ! Summary:
  ! Suppose, $X = \{x_{1}, \dotsc, x_{n_{x}}\}$, is a grid with $n_{x}$ points,
  ! such that $x_{i+1} - x_{i} = h$ for all $i = 1, \dotsc, n_{x}-1$,
  ! and that $g(x), s(x) : \mathbb{R} \to \mathbb{R}$ have been evaluated on
  ! this grid.
  ! Suppose further that $y : \mathbb{R} \to \mathbb{R}$ is defined by the
  ! differential equation,
  ! $$ y''(x) = -g(x)y(x) + s(x) $$ ,
  ! and the values of $y(x_{n_{x}-1}), y(x_{n_{x}})$ are known.
  ! The values of $y(x_{i})$ for $i = 1, \dotsc, n_{x} - 2$ are calculated to
  ! fourth-order accuracy using the backward Numerov method.
  !
  ! Input:
  ! - `n_x` is the number of grid points.
  ! - `step_size` is the distance between consecutive points on the grid $X$.
  ! - `s_grid` is the evaluation of $s(x)$ on the grid $X$.
  ! - `g_grid` is the evaluation of $g(x)$ on the grid $X$.
  !
  ! Output:
  ! - `y_grid` is the evaluation of $y(x)$ on the grid $X$, calculated via the
  !   backward Numerov method.
  ! - `status` integer status code which takes the following values:
  !   - `status == 0` indicates successful execution;
  !   - `status > 0` indicates that the arguments were invalid;
  !   - `status == -1` indicates that a numerical error (NaN or infinity)
  !     occured during execution.
  !
  subroutine numerov_b (n_x, step_size, s_grid, g_grid, y_grid, status)
    integer , intent(in) :: n_x
    double precision , intent(in) :: step_size
    double precision , intent(in) :: s_grid(n_x)
    double precision , intent(in) :: g_grid(n_x)
    double precision , intent(out) :: y_grid(n_x)
    integer , intent(out) :: status
    double precision :: step_s_grid(n_x)
    double precision :: step_g_grid(n_x)
    integer :: ii

    ! check if arguments are valid
    status = 0

    if (n_x < 1) then
      ! no grid points
      status = 1
    end if

    if (step_size < 0.0d0) then
      ! non-positive step_size
      status = 2
    end if

    ! terminate subroutine if arguments are invalid, otherwise, proceed
    if (status /= 0) then
      y_grid(:) = 0.0d0
      return
    end if

    ! perform forward-numerov method
    step_s_grid(:) = s_grid(:) * (step_size ** 2) / 12.0d0
    step_g_grid(:) = g_grid(:) * (step_size ** 2) / 12.0d0

    ! assume that y_grid(n_x), y_grid(n_x-1) already set
    ! y_grid(n_x) = y_1

    ! if (n_x >= 2) then
    !   y_grid(n_x-1) = y_2
    ! end if

    if (n_x >= 3) then
      do ii = n_x - 1, 2, -1
        y_grid(ii-1) = &
            ((2*y_grid(ii)*(1.0d0 - 5.0d0*step_g_grid(ii))) &
            - (y_grid(ii+1)*(1.0d0 + step_g_grid(ii+1))) &
            + (step_s_grid(ii+1) + 10.0d0*step_s_grid(ii) + step_s_grid(ii-1)) &
            ) / (1.0d0 + step_g_grid(ii-1))

        ! check for and handle numerical error (NaN or Infinity)
        if ((y_grid(ii-1) /= y_grid(ii-1)) .or. (y_grid(ii-1) > INFINITY)) then
          status = -1

          if (1 <= ii - 2) then
            y_grid(1:ii-2) = 0.0d0
          end if

          return
        end if
      end do
    end if

  end subroutine numerov_b

end module m_diffeq
