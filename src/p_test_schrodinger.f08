!>
program p_test_schrodinger

  use m_io
  use m_schrodinger

  ! debug
  use m_diffeq

  implicit none

  integer , parameter :: n_x = 10001
  double precision , parameter :: omega = 1.0d0
  double precision , parameter :: step_size = 20.0d0/dble(n_x)
  integer , parameter :: n_wf = 10
  double precision :: wf(n_x, n_wf)
  double precision :: energies(n_wf)
  integer :: status = 0

  double precision :: x_grid(n_x)
  integer :: ii

  ! initialise x_grid
  do ii = 1, n_x
    x_grid(ii) = step_size*(dble(ii) - dble(n_x + 1)/2.0d0)
  end do

  !
  call harmonic_oscillator(n_x, step_size, omega, n_wf, wf, energies, status)

  ! terminate program if harmonic_oscillator failed
  if (status /= 0) then
    write (*, *) "harmonic_oscillatory() failed"
    write (*, *) "exiting p_test_schrodinger"
    call exit(status)
  end if

  ! display energies
  write (*, *) "energies"
  call display_vector(n_wf, energies)

  ! visualise
  write (*, *) "wavefunctions"
  do ii = 1, n_wf
    write (*, *) "n = ", int_trim(ii)
    call display_graph(n_x, x_grid, wf(:, ii), width=80, height=30)
  end do

contains

end program p_test_schrodinger
