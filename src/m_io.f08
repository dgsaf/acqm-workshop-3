!>
module m_io

  use m_parameters , only : DP_MAX, STDOUT

  implicit none

  private
  public :: write_vector, write_matrix, write_functions, &
      display_vector, display_matrix, display_functions, display_graph, &
      dp_trim, int_trim

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

  ! display_graph
  !
  ! Brief:
  ! Display the graph of a function on an ascii-based grid.
  subroutine display_graph (n, x, f, unit, &
      left, right, bottom, top, width, height)
    integer , intent(in) :: n
    double precision , intent(in) :: x(n)
    double precision , intent(in) :: f(n)
    integer , optional , intent(in) :: unit
    double precision , optional , intent(in) :: left, right
    double precision , optional , intent(in) :: bottom, top
    integer , optional , intent(in) :: width, height
    integer :: u
    double precision :: l_x, u_x, l_y, u_y
    integer :: n_x, n_y
    double precision , allocatable :: x_seg(:), y_seg(:)
    double precision :: temp
    integer :: ii, jj, kk
    character :: symbol

    ! write to stdout if optional file unit not provided
    if (present(unit)) then
      u = unit
    else
      u = STDOUT
    end if

    ! x-axis geometry
    if (present(left)) then
      l_x = left
    else
      l_x = x(1)
    end if

    if (present(right)) then
      u_x = right
    else
      u_x = x(n)
    end if

    ! y-axis geometry
    if (present(bottom)) then
      l_y = bottom
    else
      l_y = minval(f(:))
    end if

    if (present(top)) then
      u_y = top
    else
      u_y = maxval(f(:))
    end if

    ! check ordering of left < right, bottom < top
    if (l_x > u_x) then
      temp = u_x
      u_x = l_x
      l_x = temp
    end if

    if (l_y > u_y) then
      temp = u_y
      u_y = l_y
      l_y = temp
    end if

    ! grid width and height (max 80x80)
    n_x = 80
    if (present(width) .and. (width >= 1)) then
      n_x = min(width, n_x)
    end if

    n_y = 45
    if (present(height) .and. (height >= 1)) then
      n_y = min(height, n_y)
    end if

    ! grid x, y segments
    allocate(x_seg(n_x))
    allocate(y_seg(n_y))

    do ii = 1, n_x
      x_seg(ii) = l_x + (u_x - l_x)*(ii - 1)/(n_x - 1)
    end do

    do jj = 1, n_y
      y_seg(jj) = l_y + (u_y - l_y)*(jj - 1)/(n_y - 1)
    end do

    ! write graph details

    write (u, *) &
        "(", dp_trim(l_x), " , ", dp_trim(u_x), ") x ", &
        "(", dp_trim(l_y), " , ", dp_trim(u_y), ") on ", &
        "(", int_trim(n_x), " x ", int_trim(n_y), ")"

    ! loop over grid boxes, with
    ! box(i,j) = (x_seg(i), x_seg(i+1)) x (y_seg(j), y_seg(j+1))
    do jj = n_y-1, 1, -1
      kk = 1
      do ii = 1, n_x-1
        ! if box(i,j) contains nothing, write ' '
        symbol = ' '

        ! if box(i,j) contains y=0, write '-'
        if ((y_seg(jj) <= 0.0d0) .and. (0.0d0 < y_seg(jj+1))) then
          symbol = '-'
        end if

        ! if box(i,j) contains x=0, write '|', or '+' if it contains origin
        if ((x_seg(ii) <= 0.0d0) .and. (0.0d0 < x_seg(ii+1))) then
          if (symbol == '-') then
            symbol = '+'
          else
            symbol = '|'
          end if
        end if

        ! if box(i,j) contains f(x(k)), then write '*'
        do kk = 1, n
          if ((x_seg(ii) <= x(kk)) .and. (x(kk) < x_seg(ii+1)) &
              .and. (y_seg(jj) <= f(kk)) .and. (f(kk) < y_seg(jj+1))) then
            symbol = '*'
          end if
        end do

        ! write symbol
        write (u, "(a)", advance="no") symbol
      end do

      ! newline
      write(u, *)
    end do

  end subroutine display_graph


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


  ! int_trim
  !
  ! Brief:
  ! Trim an integer number.
  function int_trim (x) result (str)
    integer , intent(in) :: x
    character(len=:) , allocatable :: str
    character(len=100) :: str_temp
    integer :: l

    ! write integer to string
    write (str_temp, "(i90)") x
    str_temp = adjustl(str_temp)

    l = len(trim(adjustl(str_temp)))

    ! write compacted format string (uses minimal space)
    allocate(character(len=l) :: str)

    write (str, "(a)") trim(adjustl(str_temp))

  end function int_trim

  ! dp_trim
  !
  ! Brief:
  ! Trim a double precision number
  function dp_trim (x) result (str)
    double precision , intent(in) :: x
    character(len=:) , allocatable :: str
    character(len=:) , allocatable :: fmt_str
    character(len=100) :: str_temp
    logical :: trimmed
    integer :: ii

    ! write format string
    call dp_format_string(90, DP_MAX, fmt_str)

    ! write dp to string
    write (str_temp, fmt_str) x
    str_temp = trim(adjustl(str_temp))

    ! search for last redundant zero/whitespace
    trimmed = .false.
    ii = len(str_temp)
    do while ((ii > 1) .and. (.not. trimmed))
      if ((str_temp(ii:ii) == " ") .or. (str_temp(ii:ii) == "0")) then
        ii = ii - 1
      else
        trimmed = .true.
      end if
    end do

    ! handle case of all decimals being zero
    if (str_temp(ii:ii) == ".") then
      ii = ii + 1
    end if

    ! write compacted format string (uses minimal space)
    allocate(character(len=ii) :: str)

    write (str, "(a)") trim(str_temp(1:ii))

  end function dp_trim

end module m_io
