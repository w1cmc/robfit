!
! nuclides_module.f90 - Fortran interface to nuclide emission library
!
! Provides access to nuclide gamma-ray emission data stored in nuclides.xml
!
module nuclides_module
  use, intrinsic :: iso_c_binding
  implicit none
  private

  ! Public interface
  public :: nuclides_load, nuclides_free, nuclides_count
  public :: nuclides_get_name, nuclides_find_index
  public :: nuclides_get_lines, nuclides_line_count

  ! C function interfaces
  interface

    function nuclide_load_c(filename) bind(C, name='nuclide_load_f')
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: filename(*)
      integer(c_int) :: nuclide_load_c
    end function

    subroutine nuclide_free_c() bind(C, name='nuclide_free_f')
    end subroutine

    function nuclide_count_c() bind(C, name='nuclide_count_f')
      import :: c_int
      integer(c_int) :: nuclide_count_c
    end function

    function nuclide_get_name_c(index, name, name_len) bind(C, name='nuclide_get_name_f')
      import :: c_int, c_char
      integer(c_int), value, intent(in) :: index
      character(kind=c_char), intent(out) :: name(*)
      integer(c_int), value, intent(in) :: name_len
      integer(c_int) :: nuclide_get_name_c
    end function

    function nuclide_find_index_c(name, name_len) bind(C, name='nuclide_find_index_f')
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), value, intent(in) :: name_len
      integer(c_int) :: nuclide_find_index_c
    end function

    function nuclide_get_lines_c(name, name_len, e_min, e_max, &
                                  energies, intensities, types, max_lines) &
                                  bind(C, name='nuclide_get_lines_f')
      import :: c_int, c_char, c_double
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), value, intent(in) :: name_len
      real(c_double), value, intent(in) :: e_min, e_max
      real(c_double), intent(out) :: energies(*)
      real(c_double), intent(out) :: intensities(*)
      character(kind=c_char), intent(out) :: types(*)
      integer(c_int), value, intent(in) :: max_lines
      integer(c_int) :: nuclide_get_lines_c
    end function

    function nuclide_line_count_c(name, name_len) bind(C, name='nuclide_line_count_f')
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), value, intent(in) :: name_len
      integer(c_int) :: nuclide_line_count_c
    end function

  end interface

contains

  !> Load nuclide library from XML file
  logical function nuclides_load(filename)
    character(len=*), intent(in) :: filename
    character(len=512) :: c_filename
    integer :: i

    ! Make null-terminated C string
    c_filename = trim(filename) // c_null_char
    
    i = nuclide_load_c(c_filename)
    nuclides_load = (i /= 0)
  end function

  !> Free nuclide library
  subroutine nuclides_free()
    call nuclide_free_c()
  end subroutine

  !> Get number of nuclides in library
  integer function nuclides_count()
    nuclides_count = nuclide_count_c()
  end function

  !> Get nuclide name by index (1-based)
  !> Returns number of emission lines, or -1 if not found
  integer function nuclides_get_name(index, name)
    integer, intent(in) :: index
    character(len=*), intent(out) :: name
    
    nuclides_get_name = nuclide_get_name_c(index, name, len(name))
  end function

  !> Find nuclide index by name (returns 1-based index)
  !> Returns 0 if not found
  integer function nuclides_find_index(name)
    character(len=*), intent(in) :: name
    
    nuclides_find_index = nuclide_find_index_c(name, len_trim(name))
  end function

  !> Get emission lines for a nuclide within energy range
  !> Returns number of lines found
  integer function nuclides_get_lines(name, e_min, e_max, &
                                       energies, intensities, types, max_lines)
    character(len=*), intent(in) :: name
    real(8), intent(in) :: e_min, e_max
    real(8), intent(out) :: energies(max_lines)
    real(8), intent(out) :: intensities(max_lines)
    character(len=1), intent(out) :: types(max_lines)
    integer, intent(in) :: max_lines
    
    nuclides_get_lines = nuclide_get_lines_c(name, len_trim(name), &
                                              e_min, e_max, &
                                              energies, intensities, types, max_lines)
  end function

  !> Get total number of lines for a nuclide
  integer function nuclides_line_count(name)
    character(len=*), intent(in) :: name
    
    nuclides_line_count = nuclide_line_count_c(name, len_trim(name))
  end function

end module nuclides_module
