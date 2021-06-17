!>
module m_schrodinger

  use m_parameters , only : TOLERANCE, MAX_ITERATIONS, SMALL, INFINITY
  use m_diffeq , only : numerov_f, numerov_b
  use m_integrate , only : integrate_trapezoid

  ! debugging
  use m_io

  implicit none

  ! private
  ! public :: shooting_bisection, numerov_cooley

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
    integer :: ii, nn, iterations

    ! debug
    write (*, *) repeat('=', 80)
    write (*, *) "harmonic_oscillator()"

    ! initialise x_grid
    do ii = 1, n_x
      x_grid(ii) = step_size*(dble(ii) - dble(n_x + 1)/2.0d0)
    end do

    ! initialise v_grid
    do ii = 1, n_x
      v_grid(ii) = 0.5d0*(omega ** 2)*(x_grid(ii) ** 2)
    end do

    do nn = 1, n_wf
      ! call shooting_bisection(n_x, step_size, x_grid, omega, v_grid, nn-1, &
      !     wf(:, nn), energies(nn), iterations, status)
      call numerov_cooley(n_x, step_size, x_grid, omega, v_grid, nn-1, &
          wf(:, nn), energies(nn), iterations, status)

      ! terminate subroutine if numerov_cooley failed
      if (status /= 0) then
        wf(:, nn:n_wf) = 0.0d0
        energies(nn:n_wf) = 0.0d0
        write (*, *) "shooting_bisection() failed"
        write (*, *) "exiting harmonic_oscillator()"
        return
      end if
    end do

    ! debug
    write (*, *) "end harmonic_oscillator()"
    write (*, *) repeat('=', 80)

  end subroutine harmonic_oscillator

  ! shooting_bisection
  subroutine shooting_bisection (n_x, step_size, x_grid, omega, v_grid, n, &
      psi_grid, energy, iterations, status)
    integer , intent(in) :: n_x
    double precision , intent(in) :: step_size
    double precision , intent(in) :: x_grid(n_x)
    double precision , intent(in) :: omega
    double precision , intent(in) :: v_grid(n_x)
    integer , intent(in) :: n
    double precision , intent(out) :: psi_grid(n_x)
    double precision , intent(out) :: energy
    integer , intent(out) :: iterations
    integer , intent(out) :: status
    double precision :: energy_min, energy_max
    double precision :: norm
    integer :: nodes
    logical :: found
    integer :: ii

    ! debug
    write (*, *) repeat('-', 80)
    write (*, *) "shooting_bisection()"
    write (*, *) "<n_x> = ", int_trim(n_x)
    write (*, *) "<step_size> = ", dp_trim(step_size)
    write (*, *) "<n> = ", int_trim(n)

    ! estimate energy_min, energy_max
    energy_min = 0.0d0
    energy_max = 5.0d0*dble(n+1)*omega

    ! loop
    found = .false.
    ii = 0
    do while ((ii <= MAX_ITERATIONS) .and. (.not. found))
      ii = ii + 1

      ! debug
      write (*, *)
      write (*, *) "iteration: ", int_trim(ii)
      write (*, *) "<energy_min> = ", dp_trim(energy_min, dp=6)
      write (*, *) "<energy_max> = ", dp_trim(energy_max, dp=6)

      energy = (energy_min + energy_max) / 2.0d0

      ! debug
      write (*, *) "<energy> = ", dp_trim(energy)

      ! solve schrodinger differential equation for psi_grid
      call solve_finite_diff(n_x, step_size, v_grid, energy, n, &
          psi_grid, status)

      ! terminate subroutine if solve_finite_diff failed
      if (status /= 0) then
        psi_grid(:) = 0.0d0
        write (*, *) "solve_finite_diff() failed"
        write (*, *) "exiting shooting_bisection"
        return
      end if

      ! count nodes
      nodes = count_nodes(n_x, psi_grid(:))

      ! debug
      write (*, *) "<nodes> = ", int_trim(nodes)

      ! if wrong number of nodes, recalibrate energy_min, energy_max and cycle
      ! otherwise proceed
      if (nodes < n) then
        energy_min = energy

        ! debug
        write (*, *) "too few nodes"
        write (*, *) "setting <energy_min> = ", dp_trim(energy_min)
        write (*, *) "cycling"
        cycle
      else if (nodes > n) then
        energy_max = energy

        ! debug
        write (*, *) "too many nodes"
        write (*, *) "setting <energy_max> = ", dp_trim(energy_max)
        write (*, *) "cycling"
        cycle
      end if

      ! debug
      write (*, *) "correct number of nodes found"

      ! check right-boundary condition
      if (psi_grid(n_x) > 0.0d0) then
        energy_min = energy

        ! debug
        write (*, *) "psi -> +ve infinity"
        write (*, *) "setting <energy_min> = ", dp_trim(energy_min)
      else if (psi_grid(n_x) < 0.0d0) then
        energy_max = energy

        ! debug
        write (*, *) "psi -> -ve infinity"
        write (*, *) "setting <energy_max> = ", dp_trim(energy_max)
      end if

      ! check if convergent
      ! debug
      write (*, *) "<energy_diff> = ", dp_trim(abs(energy_max - energy_min))
      write (*, *) "<psi_grid(n_x)> = ", dp_trim(abs(psi_grid(n_x)))

      if ((abs(energy_max - energy_min) <= TOLERANCE) &
          .or. (abs(psi_grid(n_x)) < TOLERANCE)) then
        found = .true.
      else
        ! debug
        write (*, *) "cycling"
        cycle
      end if
    end do

    ! terminate subroutine if correct nodes not found in time
    ! otherwise proceed
    if (.not. found) then
      ! debug
      call display_graph(n_x, x_grid, psi_grid, width=80, height=30)

      status = -1
      psi_grid(:) = 0.0d0
      write (*, *) "convergent <psi_grid> not found in time"
      write (*, *) "exiting shooting_bisection"
      return
    end if

    ! record iterations
    iterations = ii

    ! with right energy, normalise psi_grid
    norm = sqrt(integrate_trapezoid(n_x, x_grid, (abs(psi_grid(:)) ** 2)))
    psi_grid(:) = psi_grid(:) / norm

    ! debug
    write (*, *) "end shooting_bisection()"
    write (*, *) repeat('-', 80)

  end subroutine shooting_bisection

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
  subroutine numerov_cooley (n_x, step_size, x_grid, omega, v_grid, n, &
      psi_grid, energy, iterations, status)
    integer , intent(in) :: n_x
    double precision , intent(in) :: step_size
    double precision , intent(in) :: x_grid(n_x)
    double precision , intent(in) :: omega
    double precision , intent(in) :: v_grid(n_x)
    integer , intent(in) :: n
    double precision , intent(out) :: psi_grid(n_x)
    double precision , intent(out) :: energy
    integer , intent(out) :: iterations
    integer , intent(out) :: status
    double precision :: energy_min, energy_max
    double precision :: correction
    double precision :: norm
    integer :: nodes
    logical :: found
    integer :: ii

    ! debug
    write (*, *) repeat('-', 80)
    write (*, *) "numerov_cooley()"
    write (*, *) "<n_x> = ", int_trim(n_x)
    write (*, *) "<step_size> = ", dp_trim(step_size)
    write (*, *) "<n> = ", int_trim(n)

    ! estimate energy_min, energy_max
    energy_min = 0.0d0
    energy_max = 2.0d0*dble(n+1)*omega
    ! energy_min = dble(n)*omega
    ! energy_max = dble(n+1)*omega

    ! loop
    found = .false.
    ii = 0
    do while ((ii <= MAX_ITERATIONS) .and. (.not. found))
      ii = ii + 1

      ! debug
      write (*, *)
      write (*, *) "node-finding iteration: ", int_trim(ii)
      write (*, *) "<energy_min> = ", dp_trim(energy_min, dp=6)
      write (*, *) "<energy_max> = ", dp_trim(energy_max, dp=6)

      energy = (energy_min + energy_max) / 2.0d0

      ! debug
      write (*, *) "<energy> = ", dp_trim(energy)

      ! solve schrodinger differential equation for psi_grid
      ! call solve_finite_diff(n_x, step_size, v_grid, energy, n, &
      !     psi_grid, status)
      call solve_numerov(n_x, step_size, x_grid, v_grid, energy, n, &
          psi_grid, correction, status)

      ! terminate subroutine if solve_finite_diff failed
      if (status /= 0) then
        psi_grid(:) = 0.0d0
        write (*, *) "solve_finite_diff() failed"
        write (*, *) "exiting numerov_cooley()"
        return
      end if

      ! count nodes
      nodes = count_nodes(n_x, psi_grid(:))

      ! debug
      write (*, *) "<nodes> = ", int_trim(nodes)

      ! if wrong number of nodes, recalibrate energy_min, energy_max and cycle
      ! otherwise proceed
      if (nodes < n) then
        energy_min = energy

        ! debug
        write (*, *) "too few nodes"
        write (*, *) "setting <energy_min> = ", dp_trim(energy_min)
        write (*, *) "cycling"
        cycle
      else if (nodes > n) then
        energy_max = energy

        ! debug
        write (*, *) "too many nodes"
        write (*, *) "setting <energy_max> = ", dp_trim(energy_max)
        write (*, *) "cycling"
        cycle
      else
        found = .true.
      end if
    end do

    ! terminate subroutine if correct nodes not found in time
    if (.not. found) then
      status = -1
      psi_grid(:) = 0.0d0
      write (*, *) "correct nodes not found in time"
      write (*, *) "exiting numerov_cooley"
      return
    end if

    ! ! debug
    ! call display_graph(n_x, x_grid, psi_grid, width=80, height=30)

    ! debug
    write (*, *) "correct number of nodes found"
    write (*, *) "<energy> = ", dp_trim(energy, dp=6)

    ! record node-finding iterations
    iterations = ii

    ! with right number of nodes, determine correct energy using cooley
    found = .false.
    ii = 0
    do while ((ii <= MAX_ITERATIONS) .and. (.not. found))
      ii = ii + 1

      ! debug
      write (*, *)
      write (*, *) "energy-finding iteration: ", int_trim(ii)

      call solve_numerov(n_x, step_size, x_grid, v_grid, energy, n, &
          psi_grid, correction, status)

      ! debug
      write (*, *) "<energy> = ", dp_trim(energy, dp=6)
      write (*, *) "<correction> = ", dp_trim(correction, dp=6)

      ! ! debug
      ! call display_graph(n_x, x_grid, psi_grid, width=80, height=30)

      ! terminate subroutine if solve_numerov failed
      if (status /= 0) then
        psi_grid(:) = 0.0d0
        write (*, *) "solve_numerov() failed"
        write (*, *) "exiting numerov_cooley()"
        return
      end if

      energy = energy + correction

      found = (abs(correction) < TOLERANCE)
    end do

    ! terminate subroutine if convergent energy not found in time
    if (.not. found) then
      status = -1
      write (*, *) "convergent energy not found in time"
      write (*, *) "exiting numerov_cooley()"
      return
    end if

    ! record energy-finding iterations
    iterations = iterations + ii

    ! with right energy, normalise psi_grid
    norm = sqrt(integrate_trapezoid(n_x, x_grid, (abs(psi_grid(:)) ** 2)))
    psi_grid(:) = psi_grid(:) / norm

    ! ! debug
    ! call display_graph(n_x, x_grid, psi_grid)

    ! debug
    write (*, *) "end numerov_cooley()"
    write (*, *) repeat('-', 80)

  end subroutine numerov_cooley

  ! solve_finite_diff
  subroutine solve_finite_diff (n_x, step_size, v_grid, energy, n, &
      psi_grid, status)
    integer , intent(in) :: n_x
    double precision , intent(in) :: step_size
    double precision , intent(in) :: v_grid(n_x)
    double precision , intent(in) :: energy
    integer , intent(in) :: n
    double precision , intent(out) :: psi_grid(n_x)
    integer , intent(out) :: status
    integer :: ii

    ! debug
    write (*, *) "solve_finite_diff()"
    ! write (*, *) "<n_x> = ", int_trim(n_x)
    ! write (*, *) "<step_size> = ", dp_trim(step_size)
    ! write (*, *) "<energy> = ", dp_trim(energy)
    ! write (*, *) "<n> = ", int_trim(n)

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

    ! set up left-boundary conditions
    psi_grid(1) = 0.0d0
    psi_grid(2) = SMALL*((-1) ** (n))

    ! solve forwards using finite difference method
    do ii = 2, n_x-1
      psi_grid(ii+1) = &
          2.0d0*((step_size ** 2)*(v_grid(ii)-energy) + 1.0d0)*psi_grid(ii) &
          - psi_grid(ii-1)

      ! check for and handle numerical error (NaN or Infinity) and terminates
      if ((psi_grid(ii+1) /= psi_grid(ii+1)) &
          .or. (psi_grid(ii+1) > INFINITY)) then
        status = -1
        if (ii + 2 <= n_x) then
          psi_grid(ii+2:n_x) = 0.0d0
        end if
        return
      end if
    end do

    ! debug
    write (*, *) "end solve_finite_diff()"

  end subroutine solve_finite_diff

  ! solve_numerov
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
  ! - `n` is the predicted number of nodes of $\psi_{n}(x)$.
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
  subroutine solve_numerov (n_x, step_size, x_grid, v_grid, energy, n, &
      psi_grid, correction, status)
    integer , intent(in) :: n_x
    double precision , intent(in) :: step_size
    double precision , intent(in) :: x_grid(n_x)
    double precision , intent(in) :: v_grid(n_x)
    double precision , intent(in) :: energy
    integer , intent(in) :: n
    double precision , intent(out) :: psi_grid(n_x)
    double precision , intent(out) :: correction
    integer , intent(out) :: status
    double precision :: g_grid(n_x), s_grid(n_x)
    double precision :: psi_l_grid(n_x), psi_r_grid(n_x), temp(n_x)
    logical :: matched
    integer :: i_m
    integer :: ii

    ! debug
    write (*, *) "solve_numerov()"
    ! write (*, *) "<n_x> = ", int_trim(n_x)
    ! write (*, *) "<step_size> = ", dp_trim(step_size)
    ! write (*, *) "<energy> = ", dp_trim(energy)
    ! write (*, *) "<n> = ", int_trim(n)

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
    g_grid(:) = -2.0d0*(v_grid(:) - energy)
    s_grid(:) = 0.0d0

    ! set up boundary conditions
    psi_grid(1) = 0.0d0
    psi_grid(2) = SMALL * dble((-1) ** n)
    psi_grid(3:n_x-2) = 0.0d0
    psi_grid(n_x-1) = SMALL
    psi_grid(n_x) = 0.0d0

    ! solve forwards numerov
    psi_l_grid(:) = psi_grid(:)
    call numerov_f(n_x, step_size, s_grid, g_grid, psi_l_grid, status)

    ! terminate subroutine if numerov_f failed
    if (status /= 0) then
      psi_grid(:) = 0.0d0
      write (*, *) "numerov_f() failed"
      write (*, *) "exiting solve_numerov()"
      return
    end if

    ! solve backwards numerov
    psi_r_grid(:) = psi_grid(:)
    call numerov_b(n_x, step_size, s_grid, g_grid, psi_r_grid, status)

    ! terminate subroutine if numerov_b failed
    if (status /= 0) then
      psi_grid(:) = 0.0d0
      write (*, *) "numerov_b() failed"
      write (*, *) "exiting solve_numerov()"
      return
    end if

    ! locate where psi_l and psi_r match up
    i_m = minloc(abs(psi_l_grid(:) - psi_r_grid(:)), 1)

    matched = .false.
    do while ((i_m <= n_x) .and. (.not. matched))
      matched = (abs(psi_l_grid(i_m)) > SMALL) &
          .and. (abs(psi_r_grid(i_m)) > SMALL)

      if (.not. matched) i_m = i_m + 1
    end do

    ! terminate subroutine if psi_l, psi_r could not be matched
    if (.not. matched) then
      status = -2
      psi_grid(:) = 0.0d0
      write (*, *) "failed to match"
      write (*, *) "<psi_l_grid>"
      ! call display_vector(n_x, psi_l_grid)
      call display_graph(n_x, x_grid, psi_l_grid)
      write (*, *) "<psi_r_grid>"
      ! call display_vector(n_x, psi_r_grid)
      call display_graph(n_x, x_grid, psi_r_grid)
      write (*, *) "exiting solve_numerov()"
      return
    end if

    ! stich psi together from psi_l, psi_r
    do ii = 1, i_m
      psi_grid(ii) = psi_l_grid(ii) / psi_l_grid(i_m)
    end do

    do ii = i_m+1, n_x
      psi_grid(ii) = psi_r_grid(ii) / psi_r_grid(i_m)
    end do

    ! ! debug
    ! write (*, *) "<psi_grid>"
    ! call display_graph(n_x, x_grid, psi_grid)

    ! calculate cooley's energy correction
    correction = cooley_correction(n_x, step_size, v_grid, energy, &
        psi_grid, i_m)

    write (*, *) "<correction> = ", dp_trim(correction, dp=6)

    ! debug
    write (*, *) "end solve_numerov()"

  end subroutine solve_numerov

  function cooley_correction (n_x, step_size, v_grid, energy, &
      psi_grid, i_m) result (correction)
    integer , intent(in) :: n_x
    double precision , intent(in) :: step_size
    double precision , intent(in) :: v_grid(n_x)
    double precision , intent(in) :: energy
    double precision , intent(in) :: psi_grid(n_x)
    integer , intent(in) :: i_m
    double precision :: correction
    double precision :: g_grid(n_x), y_grid(0:n_x+1), temp_grid(n_x)
    double precision :: overlap, energy_diff, numerov_term
    integer :: ii

    ! debug
    ! write (*, *) "cooley_correction()"

    ! grid variables
    g_grid(:) = -2.0d0*(v_grid(:) - energy)

    y_grid(0) = 0.0d0
    do ii = 1, n_x
      y_grid(ii) = (1.0d0 + (((step_size ** 2)*g_grid(ii))/12.0d0))*psi_grid(ii)
    end do
    y_grid(n_x+1) = 0.0d0

    ! calculate overlap
    overlap = 0.0d0
    do ii = 1, n_x
      overlap = overlap + (abs(psi_grid(ii)) ** 2)
    end do
    overlap = overlap*step_size

    !general method
    temp_grid(1) = 0.0d0
    temp_grid(n_x) = 0.0d0
    energy_diff = 0.0d0
    do ii = 2, n_x-1
      numerov_term = y_grid(ii+1) - 2.0d0*y_grid(ii) + y_grid(ii-1)

      temp_grid(ii) = ((v_grid(ii) - energy)*psi_grid(ii)) &
          - ((0.5d0*numerov_term)/(step_size ** 2))
      energy_diff = energy_diff + (psi_grid(ii) * temp_grid(ii))
    end do
    energy_diff=energy_diff*step_size

    ! ! specific method
    ! numerov_term = y_grid(i_m+1) - 2.0d0*y_grid(i_m) + y_grid(i_m-1)

    ! energy_diff = psi_grid(i_m)*( &
    !     (v_grid(i_m) - energy)*psi_grid(i_m) &
    !     - ((0.5d0*numerov_term)/(step_size ** 2)))
    ! energy_diff=energy_diff*step_size

    correction = energy_diff / overlap

    ! debug
    ! write (*, *) "<overlap> = ", dp_trim(overlap, dp=6)
    ! write (*, *) "<energy_diff> = ", dp_trim(energy_diff, dp=6)
    ! write (*, *) "<correction> = ", dp_trim(correction, dp=6)
    ! write (*, *) "end cooley_correction()"

  end function cooley_correction

  ! count_nodes
  !
  ! Brief:
  ! Counts the number of nodes in an evaluated function.
  function count_nodes (n_x, psi_grid) result (nodes)
    integer , intent(in) :: n_x
    double precision , intent(in) :: psi_grid(n_x)
    integer :: nodes
    integer :: definite_nodes, possible_nodes
    integer :: sgn(n_x)
    integer :: ii, jj

    ! debug
    write (*, *) "count_nodes()"

    ! method: sign change advanced
    where (abs(psi_grid(:)) <= SMALL*1d1)
      sgn(:) = 0
    elsewhere (psi_grid(:) > SMALL*1d1)
      sgn(:) = 1
    elsewhere (psi_grid(:) < SMALL*1d1)
      sgn(:) = -1
    endwhere

    ! ! debug
    ! do ii = 1, n_x
    !   write (*, *) ii, psi_grid(ii), sgn(ii)
    ! end do

    ! initialise node counting
    nodes = 0
    possible_nodes = 0
    definite_nodes = 0

    ! loop through sgn(:) to find sign changes
    ii = 1
    do while (ii < n_x)
      ! write (*, *) "ii, sgn, psi = ", ii, sgn(ii), psi_grid(ii)

      ! search for first definitively non-zero point
      if (sgn(ii) == 0) then
        ! write (*, *) "cycling for non-zero point"
        ii = ii + 1
        cycle
      end if

      ! write (*, *) "is non-zero"

      ! check if next step has change in sgn(:)
      if (sgn(ii) /= sgn(ii+1)) then
        ! write (*, *) "change at ", int_trim(ii), &
        !     " (", dp_trim(psi_grid(ii)), " -> ", dp_trim(psi_grid(ii+1)), ")"

        possible_nodes = possible_nodes + 1

        ! step past a run of zero or more 0s
        jj = ii + 1
        do while (jj <= n_x)
          ! write (*, *) "zero run: ", jj, sgn(jj)
          if (sgn(jj) == 0) then
            ! write (*, *) "continues"
            jj = jj + 1
          else
            ! write (*, *) "finishes"
            exit
          end if
        end do
        ! write (*, *) "run of ", int_trim(jj - ii - 1), " zeroes"

        ! if sign has changed, is definitely a nodal region
        if (jj <= n_x) then
          if (sgn(jj) == -sgn(ii)) then
            definite_nodes = definite_nodes + 1
          end if
        end if

        ii = jj
      else
        ii = ii + 1
      end if
    end do

    nodes = definite_nodes
    ! nodes = possible_nodes

    ! debug
    ! write (*, *) "definite_nodes = ", int_trim(definite_nodes)
    ! write (*, *) "possible_nodes = ", int_trim(possible_nodes)
    ! write (*, *) "nodes = ", int_trim(nodes)
    write (*, *) "nodes = ", int_trim(nodes), &
        " (", int_trim(definite_nodes), " def / ", &
        int_trim(possible_nodes), " pos)"

    ! debug
    write (*, *) "end count_nodes()"

  end function count_nodes

end module m_schrodinger
