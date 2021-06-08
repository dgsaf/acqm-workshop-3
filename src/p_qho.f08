!>
program p_qho

  use m_io
  use m_parameters , only : PI
  use m_schrodinger , only : shooting_bisection, numerov_cooley

  implicit none

  double precision , parameter :: omega = 1.0d0
  integer , parameter :: n_wf = 4
  integer , parameter :: n_delta = 8

  integer :: n_x
  double precision :: step_size
  double precision , allocatable :: wf(:, :, :)
  double precision :: energies(n_wf, n_delta+1, 2)
  integer :: iterations(n_wf, n_delta, 2)

  double precision , allocatable :: x_grid(:)
  double precision , allocatable :: v_grid(:)
  integer :: status = 0
  integer :: ii, jj, kk, nn

  write (*, *) "qho"

  do jj = 1, n_delta
    ! set n_x and step_size
    n_x = 2 ** (jj + 9)
    step_size = 10.0d0/dble(n_x)

    ! allocate grids and wavefunctions
    allocate(x_grid(n_x))
    allocate(v_grid(n_x))
    allocate(wf(n_x, n_wf, 3))

    ! initialise x_grid
    do ii = 1, n_x
      x_grid(ii) = step_size*(dble(ii) - dble(n_x + 1)/2.0d0)
    end do

    ! initialise v_grid
    do ii = 1, n_x
      v_grid(ii) = 0.5d0*(omega ** 2)*(x_grid(ii) ** 2)
    end do

    ! initialise default wf
    wf(:, :, :) = 0.0d0

    ! calculate wavefunctions and energies using shooting_bisection
    do nn = 1, n_wf
      call shooting_bisection(n_x, step_size, x_grid, omega, v_grid, nn-1, &
          wf(:, nn, 1), energies(nn, jj, 1), iterations(nn, jj, 1), status)

      ! terminate subroutine if shooting_bisection failed
      if (status /= 0) then
        write (*, *) "shooting_bisection() failed"
        write (*, *) "exiting qho"
        call exit(status)
      end if
    end do

    ! calculate wavefunctions and energies using numerov_cooley
    do nn = 1, n_wf
      call numerov_cooley(n_x, step_size, x_grid, omega, v_grid, nn-1, &
          wf(:, nn, 2), energies(nn, jj, 2), iterations(nn, jj, 2), status)

      ! terminate subroutine if numerov_cooley failed
      if (status /= 0) then
        write (*, *) "numerov_cooley() failed"
        write (*, *) "exiting qho"
        call exit(status)
      end if
    end do

    ! calculate analytic wavefunctions of first 4 states
    do nn = 1, min(n_wf, 4)
      do kk = 1, n_x
        wf(kk, nn, 3) = hermite(nn, sqrt(omega)*x_grid(kk))
      end do
    end do

    ! ! display
    ! write (*, *) "wavefunctions:"
    ! do nn = 1, min(n_wf, 4)
    !   write (*, *) "n = ", int_trim(nn)
    !   write (*, *) "shooting_bisection:"
    !   call display_graph(n_x, x_grid, wf(:, nn, 1))
    !   write (*, *) "numerov_cooley:"
    !   call display_graph(n_x, x_grid, wf(:, nn, 2))
    !   write (*, *) "analytic:"
    !   call display_graph(n_x, x_grid, wf(:, nn, 3))
    !   write (*, *)
    ! end do

    ! write functions to file
    call write_wf(n_x, x_grid, n_wf, wf)

    ! deallocate grids
    deallocate(x_grid)
    deallocate(v_grid)
    deallocate(wf)
  end do

  ! calculate analytic energies of first 4 states
  do nn = 1, min(n_wf, 4)
    energies(nn, n_delta+1, 1) = omega*(dble(nn)-0.5d0)
    energies(nn, n_delta+1, 2) = omega*(dble(nn)-0.5d0)
  end do

  ! write energies to file
  call write_energies(n_wf, n_delta, energies)

  ! write iterations to file
  call write_iterations(n_wf, n_delta, iterations)


  ! ! display
  ! write (*, *) "energies:"
  ! write (*, *) "shooting_bisection:"
  ! call display_matrix(n_wf, n_delta, energies(:, :, 1))
  ! write (*, *) "numerov_cooley:"
  ! call display_matrix(n_wf, n_delta, energies(:, :, 2))
  ! write (*, *) "analytic:"
  ! call display_matrix(n_wf, n_delta, energies(:, :, 3))

  ! write (*, *) "iterations:"
  ! write (*, *) "shooting_bisection:"
  ! call display_matrix(n_wf, n_delta, 1.0d0*iterations(:, :, 1))
  ! write (*, *) "numerov_cooley:"
  ! call display_matrix(n_wf, n_delta, 1.0d0*iterations(:, :, 2))

  write (*, *) "end qho"

contains

  ! hermite functions up to n=4
  function hermite (n, x) result (f)
    integer , intent(in) :: n
    double precision , intent(in) :: x
    double precision :: f
    double precision :: a

    f = 0.0d0
    a = (1.0d0/sqrt(sqrt(PI)))*exp(-((x ** 2)/2.0d0))

    select case (n)
    case (1)
      f = a
    case (2)
      f = a*sqrt(2.0d0)*x
    case (3)
      f = a*(((2.0d0*x*x) - 1.0d0)/(sqrt(2.0d0)))
    case (4)
      f = a*(((2.0d0*x*x*x) - (3.0d0*x))/(sqrt(3.0d0)))
    end select

  end function hermite

  ! write wavefunctions to file
  subroutine write_wf (n_x, x_grid, n_wf, wf)
    integer , intent(in) :: n_x
    double precision , intent(in) :: x_grid(n_x)
    integer , intent(in) :: n_wf
    double precision , intent(in) :: wf(n_x, n_wf, 3)
    character(len=100) :: basename

    ! construct basename
    write (basename, *) "output/qho/wf.n_x-", int_trim(n_x), "."

    ! write shooting_bisection
    call write_functions(n_x, x_grid, n_wf, wf(:, :, 1), &
        trim(adjustl(basename))//"sb.txt")

    ! write numerov_cooley
    call write_functions(n_x, x_grid, n_wf, wf(:, :, 2), &
        trim(adjustl(basename))//"nc.txt")

    ! write analytic
    call write_functions(n_x, x_grid, n_wf, wf(:, :, 3), &
        trim(adjustl(basename))//"an.txt")

  end subroutine write_wf

  ! write energies to file
  subroutine write_energies (n_wf, n_delta, energies)
    integer , intent(in) :: n_wf
    integer , intent(in) :: n_delta
    double precision , intent(in) :: energies(n_wf, n_delta+1, 3)
    character(len=100) :: basename

    ! construct basename
    write (basename, *) "output/qho/en."

    ! write shooting_bisection
    call write_matrix(n_wf, n_delta+1, energies(:, :, 1), &
        trim(adjustl(basename))//"sb.txt")

    ! write numerov_cooley
    call write_matrix(n_wf, n_delta+1, energies(:, :, 2), &
        trim(adjustl(basename))//"nc.txt")

  end subroutine write_energies

  ! write iterations to file
  subroutine write_iterations (n_wf, n_delta, iterations)
    integer , intent(in) :: n_wf
    integer , intent(in) :: n_delta
    integer , intent(in) :: iterations(n_wf, n_delta, 2)
    character(len=100) :: basename

    ! construct basename
    write (basename, *) "output/qho/it."

    ! write shooting_bisection
    call write_matrix_int(n_wf, n_delta, iterations(:, :, 1), &
        trim(adjustl(basename))//"sb.txt")

    ! write numerov_cooley
    call write_matrix_int(n_wf, n_delta, iterations(:, :, 2), &
        trim(adjustl(basename))//"nc.txt")

  end subroutine write_iterations

end program p_qho
