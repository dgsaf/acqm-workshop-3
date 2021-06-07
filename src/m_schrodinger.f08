!>
module m_schrodinger

  use m_parameters, only : TOLERANCE, MAX_ITERATIONS
  use m_diffeq, only : numerov_f, numerov_b
  use m_integrate, only : integrate_trapezoid

  implicit none

  private
  public :: harmonic_oscillator

contains

  ! harmonic_oscillator
  ! Brief:
  !
  ! Input:
  !
  ! Output:
  !
  subroutine harmonic_oscillator (n_x, step_size, omega, n_wf, &
      wf, energies, status)
    integer , intent(in) :: n_x
    double precision , intent(in) :: step_size
    double precision , intent(in) :: omega
    integer , intent(in) :: n_wf
    double precision , intent(out) :: wf(n_x, n_wf)
    double precision , intent(out) :: energies(n_wf)
    integer , intent(out) :: status
    double precision :: x_grid(n_x)
    double precision :: v_grid(n_x)
    integer :: ii, nn

    ! initialise x_grid
    if (mod(n_x, 2) == 0) then
      x_grid(n_x/2)     = - step_size
      x_grid(n_x/2 + 1) = step_size
    else
      x_grid(n_x/2 + 1) = 0.0d0
    end if

    do ii = 1, n_x/2
      x_grid(n_x/2 + ii) = x_grid(n_x/2 + ii - 1) + step_size
      x_grid(n_x/2 - ii) = x_grid(n_x/2 - ii + 1) + step_size
    end do

    ! initialise v_grid
    do ii = 1, n_x
      v_grid(ii) = 0.5d0*(omega ** 2)*(x_grid(ii) ** 2)
    end do

    do nn = 1, n_wf
      call numerov_cooley(n_x, step_size, x_grid, v_grid, nn, &
          wf(:, nn), energies(nn), status)

      ! terminate subroutine if numerov_cooley failed
      if (status /= 0) then
        wf(:, nn:n_wf) = 0.0d0
        energies(nn:n_wf) = 0.0d0
        return
      end if
    end do

  end subroutine harmonic_oscillator


  ! numerov_cooley
  !
  ! Brief:
  ! For a 1-D Schrodginer equation, this scheme determines the wavefunction and
  ! energy of the n-th eigenstate.
  !
  ! Input:
  ! - `n_x` is the number of grid points.
  ! - `step_size` is the distance between consecutive points on the grid $X$.
  ! - `x_grid` is the grid $X$.
  ! - `v_grid` is the potential plotted on the grid $X$.
  ! - `n` is the index of the eigenstate being searched for.
  !
  ! Output:
  ! - `psi_grid` is the numerically calculated $\psi_{n}(x)$ on the grid $X$.
  ! - `energy` is the numerically calculated energy of $\psi_{n}(x)$.
  ! - `status` integer status code which takes the following values:
  !   - `status == 0` indicates successful execution;
  !   - `status > 0` indicates that the arguments were invalid;
  !   - `status == -1` indicates that a numerical error (NaN or infinity)
  !     occured during execution.
  !
  subroutine numerov_cooley (n_x, step_size, x_grid, v_grid, n, &
      psi_grid, energy, status)
    integer , intent(in) :: n_x
    double precision , intent(in) :: step_size
    double precision , intent(in) :: x_grid(n_x)
    double precision , intent(in) :: v_grid(n_x)
    integer , intent(in) :: n
    double precision , intent(out) :: psi_grid(n_x)
    double precision , intent(out) :: energy
    integer , intent(out) :: status
    double precision :: energy_min, energy_max
    double precision :: correction
    double precision :: norm
    integer :: nodes
    logical :: found
    integer :: ii

    ! estimate energy_min, energy_max
    energy_min = 0.0d0
    energy_max = 100.0d0

    ! loop
    found = .false.
    ii = 0
    do while ((ii <= MAX_ITERATIONS) .and. (.not. found))
      ii = ii + 1

      energy = (energy_min + energy_max) / 2.0d0

      ! solve_symmetric
      call solve_symmetric(n_x, step_size, v_grid, energy, n, psi_grid, &
          correction, status)

      ! terminate subroutine if solve_symmetric failed
      if (status /= 0) then
        psi_grid(:) = 0.0d0
        return
      end if

      ! count nodes
      nodes = count_nodes(n_x, psi_grid)

      ! if wrong number of nodes, recalibrate energy_min, energy_max
      if (nodes < n) then
        energy_min = energy
      else if (nodes > n) then
        energy_max = energy
      else
        found = .true.
      end if
    end do

    ! terminate subroutine if correct nodes not found in time
    if (.not. found) then
      status = -1
      psi_grid(:) = 0.0d0
      return
    end if

    ! with right number of nodes, determine correct energy using cooley
    found = .false.
    ii = 0
    do while ((ii <= MAX_ITERATIONS) .and. (.not. found))
      ii = ii + 1

      energy = energy + correction

      call solve_symmetric(n_x, step_size, v_grid, energy, n, &
          psi_grid, correction, status)

      ! terminate subroutine if solve_symmetric failed
      if (status /= 0) then
        psi_grid(:) = 0.0d0
        return
      end if

      found = (abs(correction) < TOLERANCE)
    end do

    ! terminate subroutine if convergent energy not found in time
    if (.not. found) then
      status = -1
      ! todo
      return
    end if

    ! with right energy, normalise psi_grid
    norm = integrate_trapezoid(n_x, x_grid, (abs(psi_grid(:)) ** 2))
    psi_grid(:) = psi_grid(:) / norm

  end subroutine numerov_cooley


  ! solve_symmetric
  !
  ! Brief:
  ! The 1-D Schrodinger equation, with the Hamiltonian
  ! $$ - 0.5*\psi''(x) + V(x)*\psi(x) = E*\psi(x) $$ ,
  ! is solved for bound wavefunctions using the forwards and backwards Numerov
  ! method.
  ! The potential must be evaluated on a grid $X$, consisting of points
  ! consecutively equidistant.
  ! The predicted energy and number of nodes are required.
  ! The calculated wavefunction is not normalised.
  !
  ! Input:
  ! - `n_x` is the number of grid points.
  ! - `step_size` is the distance between consecutive points on the grid $X$.
  ! - `v_grid` is the potential plotted on the grid $X$.
  ! - `energy` is the predicted energy of $\psi_{n}(x)$.
  ! - `nodes` is the predicted number of nodes of $\psi_{n}(x)$.
  !
  ! Output:
  ! - `psi_grid` is the numerically calculated $\psi_{n}(x)$ on the grid $X$.
  ! - `correction` is the numerically calculated Cooley's energy correction.
  ! - `status` integer status code which takes the following values:
  !   - `status == 0` indicates successful execution;
  !   - `status > 0` indicates that the arguments were invalid;
  !   - `status == -1` indicates that a numerical error (NaN or infinity)
  !     occured during execution.
  !
  subroutine solve_symmetric (n_x, step_size, v_grid, energy, nodes, &
      psi_grid, correction, status)
    integer , intent(in) :: n_x
    double precision , intent(in) :: step_size
    double precision , intent(in) :: v_grid(n_x)
    double precision , intent(in) :: energy
    integer , intent(in) :: nodes
    double precision , intent(out) :: psi_grid(n_x)
    double precision , intent(out) :: correction
    integer , intent(out) :: status
    double precision :: g_grid(n_x), s_grid(n_x)
    double precision :: psi_l_grid(n_x), psi_r_grid(n_x)
    double precision , parameter :: small = 1.0d-4
    double precision :: sum, numerov_term
    logical :: matched
    integer :: i_m
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
      psi_grid(:) = 0.0d0
      return
    end if

    ! set up g(x), s(x) for calls to numerov
    g_grid(:) = 2.0d0*(v_grid(:) - energy)
    s_grid(:) = 0.0d0

    ! set up boundary conditions
    psi_grid(1) = 0.0d0
    psi_grid(2) = small * ((-1) ** nodes)
    psi_grid(3:n_x-2) = 0.0d0
    psi_grid(n_x-1) = small
    psi_grid(n_x) = 0.0d0

    ! solve forwards numerov
    psi_l_grid(:) = psi_grid(:)
    call numerov_f(n_x, step_size, s_grid, g_grid, psi_l_grid, status)

    ! terminate subroutine if numerov_f failed
    if (status /= 0) then
      psi_grid(:) = 0.0d0
      return
    end if

    psi_r_grid(:) = psi_grid(:)
    call numerov_b(n_x, step_size, s_grid, g_grid, psi_r_grid, status)

    ! terminate subroutine if numerov_b failed
    if (status /= 0) then
      psi_grid(:) = 0.0d0
      return
    end if

    ! locate where psi_l and psi_r match up
    matched = .false.
    i_m = 0
    ii = 1
    do while ((ii <= n_x) .and. (.not. matched))
      if (abs(psi_l_grid(ii) - psi_r_grid(ii)) < TOLERANCE) then
        matched = .true.
        i_m = ii
      else
        ii = ii + 1
      end if
    end do

    ! terminate subroutine if psi_l, psi_r could not be matched
    if (.not. matched) then
      status = -2
      psi_grid(:) = 0.0d0
      return
    end if

    ! stich psi together from psi_l, psi_r
    do ii = 1, i_m
      psi_grid(ii) = psi_l_grid(ii) / psi_l_grid(i_m)
    end do

    do ii = i_m+1, n_x
      psi_grid(ii) = psi_r_grid(ii) / psi_r_grid(i_m)
    end do

    ! calculate cooley's energy correction
    sum = 0.0d0
    do ii = 1, n_x
      sum = sum + (abs(psi_grid(ii)) ** 2)
    end do

    numerov_term = &
        ((1.0d0 - (g_grid(i_m+1)*(step_size ** 2)/12.0d0)) &
        +(1.0d0 - (g_grid(i_m-1)*(step_size ** 2)/12.0d0)) &
        -2.0d0*(1.0d0 - (g_grid(i_m)*(step_size ** 2)/12.0d0)))

    correction = &
        ((v_grid(i_m) - energy)*psi_grid(i_m) &
        -(0.5d0*numerov_term/(step_size ** 2)) &
        )*psi_grid(i_m)/sum

  end subroutine solve_symmetric

  ! count_nodes
  !
  ! Brief:
  ! Counts the number of nodes in an evaluated function.
  function count_nodes (n_x, psi_grid) result (nodes)
    integer , intent(in) :: n_x
    double precision , intent(in) :: psi_grid(n_x)
    integer :: nodes
    integer :: ii

    nodes = 0
    do ii = 2, n_x - 1
      if (((psi_grid(ii-1) < psi_grid(ii)) &
          .and. (psi_grid(ii+1) < psi_grid(ii))) &
          .or. &
          ((psi_grid(ii-1) > psi_grid(ii)) &
          .and. (psi_grid(ii+1) > psi_grid(ii)))) then
        nodes = nodes + 1
      end if
    end do

  end function count_nodes

end module m_schrodinger
