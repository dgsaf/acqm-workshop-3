!>
module m_io

  use m_parameters, only : DP_MAX, STDOUT

  implicit none

  private
  public :: write_vector, write_matrix, write_functions, display_vector, &
      display_matrix, display_functions

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

  ! write_functions
  subroutine write_functions (n_x, x_grid, n_f, f, filename)
    integer , intent(in) :: n_x
    double precision , intent(in) :: x_grid(n_x)
    integer , intent(in) :: n_f
    double precision , intent(in) :: f(n_x, n_f)
    character(len=*) , intent(in) :: filename
    integer :: fileunit
    integer :: ii

    ! open file
    fileunit = 10

    open (unit=fileunit, file=trim(adjustl(filename)), action="write")

    ! write x_grid, and functions to file
    do ii = 1, n_x
      write (fileunit, *) x_grid(ii), " ", f(ii, :)
    end do

    ! close file
    close (fileunit)

  end subroutine write_functions

  ! display_vector
  subroutine display_vector (n, x, unit, dp)
    integer , intent(in) :: n
    double precision , intent(in) :: x(n)
    integer , optional , intent(in) :: unit
    integer , optional , intent(in) :: dp
    integer :: u
    integer :: w, d
    character(len=:) , allocatable :: fmt_str, zero_str
    integer :: ii

    ! write to stdout if optional file unit not provided
    if (present(unit)) then
      u = unit
    else
      u = STDOUT
    end if

    ! determine double precision number formatting
    call dp_format(n, x, w, d, decimals=dp)
    call dp_format_string(w, d, fmt_str)
    call dp_zero_string(w, d, zero_str)

    ! write out vector elements
    do ii = 1, n
      ! if x(ii) will be written as "0.00..0", replace with " .     "
      if (abs(x(ii)) > (10.0d0**(-d))) then
        write (u, fmt_str) x(ii)
      else
        write (u, "(a)") zero_str
      end if
    end do

  end subroutine display_vector

  ! display_matrix
  subroutine display_matrix (n_rows, n_cols, A, unit, dp)
    integer , intent(in) :: n_rows, n_cols
    double precision , intent(in) :: A(n_rows, n_cols)
    integer , optional, intent(in) :: unit
    integer , optional, intent(in) :: dp
    integer :: u
    integer :: w, d
    character(len=:) , allocatable :: fmt_str, zero_str
    integer :: ii, jj

    ! write to stdout if optional file unit not provided
    if (present(unit)) then
      u = unit
    else
      u = STDOUT
    end if

    ! determine double precision number formatting
    call dp_format(n_rows*n_cols, reshape(A, (/n_rows*n_cols/)), w, d, &
        decimals=dp)
    call dp_format_string(w, d, fmt_str)
    call dp_zero_string(w, d, zero_str)

    ! write out matrix elements
    do ii = 1, n_rows
      do jj = 1, n_cols
        ! if A(ii, jj) will be written as "0.00..0", replace with " .     "
        if (abs(A(ii, jj)) > (10.0d0**(-d))) then
          write (u, fmt_str, advance="no") A(ii, jj)
        else
          write (u, "(a)", advance="no") zero_str
        end if
      end do
      write (u, *)
    end do

  end subroutine display_matrix

  ! display_functions
  subroutine display_functions (n_x, x_grid, n_f, f, unit, dp)
    integer , intent(in) :: n_x
    double precision , intent(in) :: x_grid(n_x)
    integer , intent(in) :: n_f
    double precision , intent(in) :: f(n_x, n_f)
    integer , optional , intent(in) :: unit
    integer , optional , intent(in) :: dp
    integer :: u
    integer :: w, d
    character(len=:) , allocatable :: fmt_str, zero_str
    integer :: ii, jj

    ! write to stdout if optional file unit not provided
    if (present(unit)) then
      u = unit
    else
      u = STDOUT
    end if

    ! determine double precision number formatting
    call dp_format(n_x*n_f, reshape(f, (/n_x*n_f/)), w, d, decimals=dp)
    call dp_format_string(w, d, fmt_str)
    call dp_zero_string(w, d, zero_str)

    ! write out radial grid and radial basis values
    do ii = 1, n_x
      write (u, fmt_str, advance="no") x_grid(ii)
      do jj = 1, n_f
        ! if f(ii, jj) will be written as "0.00..0", replace with " .     "
        if (abs(f(ii, jj)) > (10.0d0**(-d))) then
          write (u, fmt_str, advance="no") f(ii, jj)
        else
          write (u, "(a)", advance="no") zero_str
        end if
      end do
      write (u, *)
    end do

  end subroutine display_functions

  ! ! display_graph
  ! subroutine display_function_graph (n_x, x_grid, f_grid, dp)

  ! end subroutine display_function_graph


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
  subroutine dp_format (n, x, w, d, decimals)
    integer , intent(in) :: n
    double precision , intent(in) :: x(n)
    integer , intent(out) :: w, d
    integer , optional , intent(in) :: decimals
    ! ! alternative method variables
    ! character(len=100) :: temp, fmt_temp, str_d
    ! integer :: l
    ! integer :: ii

    ! determine decimal places, either from optional argument or predefined
    ! global parameter
    if (present(decimals)) then
      d = decimals
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
  subroutine dp_format_string (w, d, str)
    integer , intent(in) :: w, d
    character(len=:) , allocatable , intent(out) :: str
    character(len=100) :: str_w, str_d, str_temp
    integer :: l

    ! write width, decimal precision to string
    write (str_w, *) w
    write (str_d, *) d

    ! write format string
    write (str_temp, "(a,a,a,a,a)") &
        "(f", trim(adjustl(str_w)), ".", trim(adjustl(str_d)), ")"

    ! determine length of (non-empty part of) string
    l = len(trim(adjustl(str_temp)))

    ! write compacted format string (uses minimal space)
    allocate(character(len=l) :: str)

    write (str, "(a)") trim(adjustl(str_temp))

  end subroutine dp_format_string

  ! dp_zero_string
  !
  ! Brief:
  ! Given the decimal precision and width edit descriptors for pretty-printing a
  ! double precision array, construct a string used to represent values which
  ! would otherwise be displayed as " ... 0.0...0".
  ! Currently, this string is constructed as " . " with the decimal place
  ! aligned with non-zero elements.
  subroutine dp_zero_string (w, d, str)
    integer , intent(in) :: w, d
    character(len=:) , allocatable , intent(out) :: str

    ! allocate
    allocate(character(len=w) :: str)

    ! write zero string
    str = repeat(' ', w)
    str(w-d:w-d) = '.'

  end subroutine dp_zero_string

end module m_io
