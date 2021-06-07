!>
program p_test_io

  use m_random
  use m_io

  implicit none

  integer , parameter :: seed = 102131

  call test_display_vector(seed, 10, 0.0d0, 1.0d3)
  call test_display_vector(seed, 10, 0.0d0, 1.0d0)
  call test_display_vector(seed, 10, 0.0d0, 1.0d-4)

  call test_display_matrix(seed, 5, 5, 0.0d0, 1.0d3)
  call test_display_matrix(seed, 5, 5, 0.0d0, 1.0d0)
  call test_display_matrix(seed, 5, 5, 0.0d0, 1.0d-4)

contains

  ! test_display_vector
  subroutine test_display_vector (seed, n, mu, sigma)
    integer , intent(in) :: seed
    integer , intent(in) :: n
    double precision , intent(in) :: mu, sigma
    double precision :: x(n)
    double precision :: r1, r2
    integer :: ii

    write (*, *) "test_display_vector"
    write (*, *) "  `seed`  = ", seed
    write (*, *) "  `n`     = ", n
    write (*, *) "  `mu`    = ", mu
    write (*, *) "  `sigma` = ", sigma

    call prepare_random(seed)

    do ii = 1, n, 2
      call random_normal(mu, sigma, r1, r2)
      x(ii) = r1
      if (ii+1 <= n) then
        x(ii+1) = r2
      end if
    end do

    write (*, *) "display `x`"
    call display_vector(n, x)

    write (*, *) "display `x` with `dp`=7"
    call display_vector(n, x, dp=7)

    write (*, *) "end test_display_vector"

  end subroutine test_display_vector

  ! test_display_matrix
  subroutine test_display_matrix (seed, n_rows, n_cols, mu, sigma)
    integer , intent(in) :: seed
    integer , intent(in) :: n_rows, n_cols
    double precision , intent(in) :: mu, sigma
    double precision :: x(n_rows*n_cols)
    double precision :: A(n_rows, n_cols)
    double precision :: r1, r2
    integer :: ii

    write (*, *) "test_display_matrix"
    write (*, *) "  `seed`  = ", seed
    write (*, *) "  `n_rows = ", n_rows
    write (*, *) "  `n_cols = ", n_cols
    write (*, *) "  `mu`    = ", mu
    write (*, *) "  `sigma` = ", sigma

    call prepare_random(seed)

    do ii = 1, n_rows*n_cols, 2
      call random_normal(mu, sigma, r1, r2)
      x(ii) = r1
      if (ii+1 <= n_rows*n_cols) then
        x(ii+1) = r2
      end if
    end do

    A = reshape(x, (/n_rows, n_cols/))

    write (*, *) "display `A`"
    call display_matrix(n_rows, n_cols, A)

    write (*, *) "display `A` with `dp`=7"
    call display_matrix(n_rows, n_cols, A, dp=7)

    write (*, *) "end test_display_matrix"

  end subroutine test_display_matrix

end program p_test_io
