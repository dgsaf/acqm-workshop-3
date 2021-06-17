!>
program p_h2_dissociate

  use m_io
  use m_parameters , only : PI
  use m_schrodinger , only : count_nodes

  implicit none

  double precision , parameter :: mu = 918.07635
  ! integer , parameter :: n_mu = 6
  ! double precision , parameter :: mu(n_mu) = &
  !     (/918.07635, 1223.89925, 1376.39236, 1835.24151, 2200.87999, 2748.46079/)

  integer :: n_r
  double precision :: r_max
  double precision :: d_r

  integer , parameter :: n_wf = 20
  double precision , allocatable :: wf(:, :, :)
  double precision :: ek_grid(n_wf)
  double precision :: D(2)
  double precision :: asym_amp, k

  double precision , allocatable :: r_grid(:)
  double precision , allocatable :: v_grid(:, :)

  double precision :: dcs(n_wf, 10)

  integer :: status = 0
  integer :: ii, jj, nn

  write (*, *) "h2_dissociate"

  ! set n_r and d_r
  r_max = 100.0d0
  d_r = 0.01d0
  n_r = ceiling(r_max / d_r) + 1

  ! allocate grids and wavefunctions
  allocate(r_grid(n_r))
  allocate(v_grid(n_r, 2))
  allocate(wf(n_r, n_wf, 2))

  ! initialise r_grid
  do ii = 1, n_r
    r_grid(ii) = d_r * dble(ii-1)
  end do

  ! initialise v_grid (for 1ssg/2psu)
  call interpolate_potential(n_r, r_grid, "PECS/PEC.1ssg", &
      v_grid(:, 1), status)

  call interpolate_potential(n_r, r_grid, "PECS/PEC.2psu", &
      v_grid(:, 2), status)

  ! initialise ek_grid (0 to 30ev)
  do nn = 1, n_wf
    ek_grid(nn) = 0.1d0 + (30.0d0 * (dble(nn - 1) / dble(n_wf - 1)))
  end do

  ! scale to Ha
  ek_grid(:) = ek_grid(:) / 27.21136d0

  ! initialise dissociative energy
  D(:) = v_grid(n_r, :)

  ! initialise default wf
  wf(:, :, :) = 0.0d0

  ! calculate wavefunctions using numerov method
  do nn = 1, n_wf
    k = sqrt(2.0d0 * mu * ek_grid(nn))
    asym_amp = sqrt((2.0d0 * mu) / (k * PI))

    ! 1ssg
    call solve_h2_numerov(n_r, d_r, r_grid, v_grid(:, 1), &
        ek_grid(nn) + D(1), mu, wf(:, nn, 1), status)

    ! normalise (scale by asymptotic amplitude)
    wf(:, nn, 1) = (wf(:, nn, 1) / maxval(abs(wf(n_r/2:, nn, 1)))) * asym_amp

    ! terminate subroutine if numerov_cooley failed
    if (status /= 0) then
      write (*, *) "numerov_cooley() failed"
      write (*, *) "exiting h2_dissociate"
      call exit(status)
    end if

    ! 2psu
    call solve_h2_numerov(n_r, d_r, r_grid, v_grid(:, 2), &
        ek_grid(nn) + D(2),  mu, wf(:, nn, 2), status)

    ! normalise (scale by asymptotic amplitude)
    wf(:, nn, 2) = (wf(:, nn, 2) / maxval(abs(wf(n_r/2:, nn, 2)))) * asym_amp

    ! terminate subroutine if numerov_cooley failed
    if (status /= 0) then
      write (*, *) "numerov_cooley() failed"
      write (*, *) "exiting h2_dissociate"
      call exit(status)
    end if
  end do

  ! ! display
  ! write (*, *) "wavefunctions:"

  ! write (*, *) "1ssg:"
  ! call display_graph(n_r, r_grid, v_grid(:, 1))
  ! do nn = 1, min(n_wf, 10)
  !   write (*, *) "n  = ", int_trim(nn)
  !   write (*, *) "en = ", dp_trim(ek_grid(nn) + D(1))
  !   call display_graph(n_r, r_grid, wf(:, nn, 1))
  !   write (*, *)
  ! end do

  ! write (*, *) "2psu:"
  ! call display_graph(n_r, r_grid, v_grid(:, 2))
  ! do nn = 1, min(n_wf, 10)
  !   write (*, *) "n = ", int_trim(nn)
  !   write (*, *) "en = ", dp_trim(ek_grid(nn) + D(2))
  !   call display_graph(n_r, r_grid, wf(:, nn, 2))
  !   write (*, *)
  ! end do

  ! franck cordon-approx (using values from assignment 2)
  call franck_condon(n_r, r_grid, n_wf, wf, dcs, status)

  ! write functions to file

  ! write energies to file

  ! write iterations to file

  ! deallocate grids
  deallocate(r_grid)
  deallocate(v_grid)
  deallocate(wf)

  write (*, *) "end h2_dissociate"

end program p_h2_dissociate


! interpolate_potential
!
! For given <n_r>, <r_grid>, <filename>, read the potential from
! "<filename>" and interpolate the potential onto <r_grid>.
subroutine interpolate_potential (n_r, r_grid, filename, v_grid, status)
  implicit none

  integer , intent(in) :: n_r
  double precision , intent(in) :: r_grid(n_r)
  character(len=*) , intent(in) :: filename
  double precision , intent(out) :: v_grid(n_r)
  integer , intent(out) :: status
  double precision , allocatable :: r_grid_raw(:), v_grid_raw(:)
  double precision :: r, v
  integer :: n_r_raw
  integer :: fileunit
  integer :: ii
  integer :: io_status

  ! open file
  fileunit = 10

  open (unit=fileunit, file=trim(adjustl(filename)), action="read")

  n_r_raw = 0
  io_status = 0
  do while (io_status == 0)
    read (fileunit, *, iostat=io_status) r, v
    if (io_status == 0) then
      n_r_raw = n_r_raw + 1
    end if
  end do

  ! handle invalid file read
  if (n_r_raw == 0) then
    status = 1
  end if

  ! handle invalid file read
  if (status /= 0) then
    status = 1
    return
  end if

  ! re-open file at beginning
  close (fileunit)
  open (unit=fileunit, file=trim(adjustl(filename)), action="read")

  ! allocate <r_grid_raw>, <v_grid_raw>
  allocate(r_grid_raw(n_r_raw))
  allocate(v_grid_raw(n_r_raw))

  ! read potential file
  ii = 1
  io_status = 0
  do while ((ii <= n_r_raw) .and. (io_status == 0))
    read (fileunit, *, iostat=io_status) r_grid_raw(ii), v_grid_raw(ii)
    if (io_status == 0) then
      ii = ii + 1
    end if
  end do

  ! handle invalid file read
  if (ii <= n_r_raw) then
    status = 1
  end if

  ! handle invalid file read
  if (status /= 0) then
    status = 1
    return
  end if

  ! close file
  close (fileunit)

  ! interpolate <v_grid> on <r_grid> from <v_grid_raw> on <r_grid_raw>
  call INTRPL(n_r_raw-1, r_grid_raw, v_grid_raw, n_r, r_grid, v_grid)

end subroutine interpolate_potential

! solve_h2_numerov
subroutine solve_h2_numerov (n_r, step_size, r_grid, v_grid, energy, mu, &
    wf, status)
  use m_parameters
  use m_diffeq
  use m_integrate

  implicit none

  integer , intent(in) :: n_r
  double precision , intent(in) :: step_size
  double precision , intent(in) :: r_grid(n_r)
  double precision , intent(in) :: v_grid(n_r)
  double precision , intent(in) :: energy
  double precision , intent(in) :: mu
  double precision , intent(out) :: wf(n_r)
  integer , intent(out) :: status
  double precision :: g_grid(n_r), s_grid(n_r)

  ! debug
  write (*, *) "solve_h2_numerov()"

  ! check if arguments are valid
  status = 0

  if (n_r < 1) then
    ! no grid points
    status = 1
  end if

  if (step_size < 0.0d0) then
    ! non-positive step_size
    status = 2
  end if

  ! terminate subroutine if arguments are invalid, otherwise proceed
  if (status /= 0) then
    wf(:) = 0.0d0
    return
  end if

  ! set up g(x), s(x) for calls to numerov
  g_grid(:) = 2.0d0*mu*(energy - v_grid(:))
  s_grid(:) = 0.0d0

  ! set up boundary conditions
  wf(1) = 0.0d0
  ! wf(2) = 0.0d0
  wf(2) = 1.0d-6

  ! solve forwards numerov
  call numerov_f(n_r, step_size, s_grid, g_grid, wf, status)

  ! terminate subroutine if numerov_f failed
  if (status /= 0) then
    wf(:) = 0.0d0
    write (*, *) "numerov_f() failed"
    write (*, *) "exiting solve_h2_numerov()"
    return
  end if

  ! debug
  write (*, *) "end solve_h2_numerov()"

end subroutine solve_h2_numerov

! franck_condon
! Read vibrational wavefunctions from file for 1ssg
subroutine franck_condon (n_r, r_grid, n_wf, wf, dcs, status)
  use m_integrate

  implicit none

  integer , intent(in) :: n_r
  double precision , intent(in) :: r_grid(n_r)
  integer , intent(in) :: n_wf
  double precision , intent(in) :: wf(n_r, n_wf, 2)
  double precision , intent(out) :: dcs(n_wf, 10)
  integer , intent(out) :: status
  character(len=*) , parameter :: filename = "VIB/vibrational_wf_128.txt"
  double precision :: vib_wf(n_r, 128)
  double precision :: r
  integer :: fileunit
  integer :: ii, jj, nn
  integer :: io_status

  ! open file
  fileunit = 10

  open (unit=fileunit, file=trim(adjustl(filename)), action="read")

  ! read vibrational wavefunctions file
  ii = 1
  io_status = 0
  do while ((ii <= n_r) .and. (io_status == 0))
    read (fileunit, *, iostat=io_status) r, vib_wf(ii, :)
    if (io_status == 0) then
      ii = ii + 1
    end if
  end do

  ! handle invalid file read
  if (status /= 0) then
    status = 1
    return
  end if

  ! close file
  close (fileunit)

  ! calculate franck_condon factors (1ssg)
  dcs(:, :) = 0.0d0
  do jj = 1, 10
    do nn = 1, n_wf
      dcs(nn, jj) = integrate_trapezoid(n_r, r_grid, &
          abs(vib_wf(:, jj) * wf(:, nn, 1)) ** 2)
    end do
  end do

end subroutine franck_condon
