!>
module m_io

  use m_parameters, only : DP_MAX

  implicit none

contains

  ! write_vector
  subroutine write_vector (n, x, filename)
    integer , intent(in) :: n
    double precision , intent(in) :: x(n)
    character(len=*) , intent(in) :: filename
    integer :: fileunit
    integer :: ii

    ! open file
    fileunit = 10

    open (unit=fileunit, file=trim(adjustl(filename)), action="write")

    ! write vector to file
    do ii = 1, n
      write (fileunit, *) x(ii)
    end do

    ! close file
    close (fileunit)

  end subroutine write_vector

  ! write_matrix
  subroutine write_matrix (n_rows, n_cols, A, filename)
    integer , intent(in) :: n_rows, n_cols
    double precision , intent(in) :: A(n_rows, n_cols)
    character(len=*) , intent(in) :: filename
    integer :: fileunit
    integer :: ii

    ! open file
    fileunit = 10

    open (unit=fileunit, file=trim(adjustl(filename)), action="write")

    ! write matrix to file
    do ii = 1, n_rows
      write (fileunit, *) A(ii, :)
    end do

    ! close file
    close (fileunit)

  end subroutine write_matrix

  ! write_basis
  subroutine write_basis (n_r, r_grid, n_basis, basis, filename)
    integer , intent(in) :: n_r, n_basis
    double precision , intent(in) :: r_grid(n_r)
    double precision , intent(in) :: basis(n_r, n_basis)
    character(len=*) , intent(in) :: filename
    integer :: fileunit
    integer :: ii

    ! open file
    fileunit = 10

    open (unit=fileunit, file=trim(adjustl(filename)), action="write")

    ! write r_grid, and basis functions to file
    do ii = 1, n_r
      write (fileunit, *) r_grid(ii), " ", basis(ii, :)
    end do

    ! close file
    close (fileunit)

  end subroutine write_basis

  ! display_vector
  subroutine display_vector (n, x, dp)
    integer , intent(in) :: n
    double precision , intent(in) :: x(n)
    integer , optional , intent(in) :: dp
    integer :: w, d
    character(len=50) :: fmt, str_w, str_d, str_zero
    integer :: ii

    ! determine double precision number formatting
    if (present(dp)) then
      d = dp
    else
      d = DP_MAX
    end if

    w = max(ceiling(log10(maxval(abs(x(:))))), 1) + d + 3

    write (str_w, *) w
    write (str_d, *) d

    write (fmt, *) "(f", trim(adjustl(str_w)), ".", trim(adjustl(str_d)), ")"

    str_zero = repeat(' ', w)

    ! write out vector elements
    do ii = 1, n
      ! if x(ii) will be written as "0.00..0", replace with " .     "
      if (abs(x(ii)) > (10.0**(-d))) then
        write (*, fmt) x(ii)
      else
        write (*, "(a, a, a)") str_zero(1:w-d-1), ".", str_zero(w-d+1:w)
      end if
    end do

  end subroutine display_vector

  ! display_matrix
  subroutine display_matrix (n_rows, n_cols, A, dp)
    integer , intent(in) :: n_rows, n_cols
    double precision , intent(in) :: A(n_rows, n_cols)
    integer , optional, intent(in) :: dp
    integer :: w, d
    character(len=50) :: fmt, str_w, str_d, str_zero
    integer :: ii, jj

    ! determine double precision number formatting
    if (present(dp)) then
      d = dp
    else
      d = DP_MAX
    end if

    w = max(ceiling(log10(maxval(abs(A(:, :))))), 1) + d + 3

    write (str_w, *) w
    write (str_d, *) d

    write (fmt, *) "(f", trim(adjustl(str_w)), ".", trim(adjustl(str_d)), ")"

    str_zero = repeat(' ', w)

    ! write out matrix elements
    do ii = 1, n_rows
      do jj = 1, n_cols
        ! if A(ii, jj) will be written as "0.00..0", replace with " .     "
        if (abs(A(ii, jj)) > (10.0**(-d))) then
          write (*, fmt, advance="no") A(ii, jj)
        else
          write (*, "(a, a, a)", advance="no") &
              str_zero(1:w-d-1), ".", str_zero(w-d+1:w)
        end if
      end do
      write (*, *)
    end do

  end subroutine display_matrix

  ! display_basis
  subroutine display_basis (n_r, r_grid, n_basis, basis, dp)
    integer , intent(in) :: n_r, n_basis
    double precision , intent(in) :: r_grid(n_r)
    double precision , intent(in) :: basis(n_r, n_basis)
    integer , optional , intent(in) :: dp
    integer :: w, d
    character(len=50) :: fmt, str_w, str_d, str_zero
    integer :: ii, jj

    ! determine double precision number formatting
    if (present(dp)) then
      d = dp
    else
      d = DP_MAX
    end if

    w = max(ceiling(log10(maxval(abs(basis(:, :))))), 1) + d + 3

    write (str_w, *) w
    write (str_d, *) d

    write (fmt, *) "(f", trim(adjustl(str_w)), ".", trim(adjustl(str_d)), ")"

    str_zero = repeat(' ', w)

    ! write out radial grid and radial basis values
    do ii = 1, n_r
      write (*, fmt, advance="no") r_grid(ii)
      do jj = 1, n_basis
        ! if basis(ii, jj) will be written as "0.00..0", replace with " .     "
        if (abs(basis(ii, jj)) > (10.0**(-d))) then
          write (*, fmt, advance="no") basis(ii, jj)
        else
          write (*, "(a, a, a)", advance="no") &
              str_zero(1:w-d-1), ".", str_zero(w-d+1:w)
        end if
      end do
      write (*, *)
    end do

  end subroutine display_basis


  ! dp_format
  !
  ! Brief:
  ! Determines the decimal precision and width edit descriptors required for
  ! pretty-printing the double precision values of an array.
  ! The decimal precision can be optionally provided as an argument, and is
  ! otherwise determined by the global parameter `DP_MAX`.
  ! The width is determined to be the minimum width required to print any
  ! element of the array, positive or negative, with the required level of
  ! decimal precision.
  subroutine dp_format (n, x, dp, w, d)
    integer , intent(in) :: n
    double precision , intent(in) :: x(n)
    integer , optional , intent(in) :: dp
    integer , intent(out) :: w, d
    ! ! alternative method variables
    ! character(len=100) :: temp, fmt_temp, str_d
    ! integer :: l
    ! integer :: ii


    ! determine decimal places, either from optional argument or predefined
    ! global parameter
    if (present(dp)) then
      d = dp
    else
      d = DP_MAX
    end if

    ! determine minimum width required
    w = max(ceiling(log10(maxval(abs(x(:))))), 1) + d + 3

    ! ! alternative method
    ! write (str_d, "(i)") d
    ! write (fmt_temp, "(i)") "(f0.", trim(adjustl(str_d)), ")"
    ! l = 0
    ! do ii = 1, n
    !   write (temp, fmt_temp) x(ii)
    !   l = max(l, len(trim(adjustl(temp))))
    ! end do
    ! w = l

  end subroutine dp_format

  ! dp_format_string
  !
  ! Brief:
  ! Given the decimal precision and width edit descriptors for pretty-printing a
  ! double precision array, construct a format string suitable for using in
  ! i/o statements.
  function dp_format_string (w, d) result (str)
    integer , intent(in) :: w, d
    character(len=:) , allocatable :: str
    character(len=100) :: str_w, str_d, str_temp
    integer :: l

    ! write width, decimal precision to string
    write (str_w, *) w
    write (str_d, *) d

    ! write format string
    write (str_temp, *) &
        "(f", trim(adjustl(str_w)), ".", trim(adjustl(str_d)), ")"

    ! determine length of (non-empty part of) string
    l = len(trim(adjustl(str_temp)))

    ! write compacted format string (uses minimal space)
    allocate(character(len=l) :: str)

    write (str, *) "(f", trim(adjustl(str_w)), ".", trim(adjustl(str_d)), ")"

  end function dp_format_string

  ! dp_zero_string
  !
  ! Brief:
  ! Given the decimal precision and width edit descriptors for pretty-printing a
  ! double precision array, construct a string used to represent values which
  ! would otherwise be displayed as " ... 0.0...0".
  ! Currently, this string is constructed as " . " with the decimal place
  ! aligned with non-zero elements.
  function dp_zero_string (w, d) result (str)
    integer , intent(in) :: w, d
    character(len=w) :: str

    ! write zero string
    str = repeat(' ', w)
    str(w-d:w-d) = '.'

  end function dp_zero_string

end module m_io
