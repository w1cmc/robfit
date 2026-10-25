!
! config_module.f90 - Fortran interface to INI configuration
!
! This module provides Fortran-friendly wrappers around the C config API,
! allowing Fortran code to read and write INI configuration files.
!
! Usage:
!   use config_module
!   
!   if (config_load('robfit.ini')) then
!       spectrum_file = config_get_string('spectrum', 'file', 'NONE')
!       energy_begin = config_get_double('display', 'energy_begin', 0.0d0)
!       view_data = config_get_bool('display', 'view_data', .true.)
!   end if
!
module config_module
    use iso_c_binding
    implicit none
    private
    
    ! Public procedures
    public :: config_load
    public :: config_new
    public :: config_free
    public :: config_save
    public :: config_get_string
    public :: config_get_int
    public :: config_get_double
    public :: config_get_bool
    public :: config_set_string
    public :: config_set_int
    public :: config_set_double
    public :: config_set_bool
    public :: config_has_section
    public :: config_has_key
    
    ! C function interfaces
    interface
        function config_load_c(filename) bind(C, name="config_load_f")
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: filename(*)
            integer(c_int) :: config_load_c
        end function
        
        function config_new_c() bind(C, name="config_new_f")
            import :: c_int
            integer(c_int) :: config_new_c
        end function
        
        subroutine config_free_c() bind(C, name="config_free_f")
        end subroutine
        
        function config_save_c(filename) bind(C, name="config_save_f")
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: filename(*)
            integer(c_int) :: config_save_c
        end function
        
        function config_get_string_c(section, key, default_val, result, result_len) &
                bind(C, name="config_get_string_f")
            import :: c_int, c_char, c_ptr
            character(kind=c_char), intent(in) :: section(*)
            character(kind=c_char), intent(in) :: key(*)
            character(kind=c_char), intent(in) :: default_val(*)
            character(kind=c_char), intent(out) :: result(*)
            integer(c_int), value :: result_len
            integer(c_int) :: config_get_string_c
        end function
        
        function config_get_int_c(section, key, default_val) bind(C, name="config_get_int_f")
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: section(*)
            character(kind=c_char), intent(in) :: key(*)
            integer(c_int), value :: default_val
            integer(c_int) :: config_get_int_c
        end function
        
        function config_get_double_c(section, key, default_val) bind(C, name="config_get_double_f")
            import :: c_int, c_char, c_double
            character(kind=c_char), intent(in) :: section(*)
            character(kind=c_char), intent(in) :: key(*)
            real(c_double), value :: default_val
            real(c_double) :: config_get_double_c
        end function
        
        function config_get_bool_c(section, key, default_val) bind(C, name="config_get_bool_f")
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: section(*)
            character(kind=c_char), intent(in) :: key(*)
            integer(c_int), value :: default_val
            integer(c_int) :: config_get_bool_c
        end function
        
        function config_set_string_c(section, key, value) bind(C, name="config_set_string_f")
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: section(*)
            character(kind=c_char), intent(in) :: key(*)
            character(kind=c_char), intent(in) :: value(*)
            integer(c_int) :: config_set_string_c
        end function
        
        function config_set_int_c(section, key, value) bind(C, name="config_set_int_f")
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: section(*)
            character(kind=c_char), intent(in) :: key(*)
            integer(c_int), value :: value
            integer(c_int) :: config_set_int_c
        end function
        
        function config_set_double_c(section, key, value) bind(C, name="config_set_double_f")
            import :: c_int, c_char, c_double
            character(kind=c_char), intent(in) :: section(*)
            character(kind=c_char), intent(in) :: key(*)
            real(c_double), value :: value
            integer(c_int) :: config_set_double_c
        end function
        
        function config_set_bool_c(section, key, value) bind(C, name="config_set_bool_f")
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: section(*)
            character(kind=c_char), intent(in) :: key(*)
            integer(c_int), value :: value
            integer(c_int) :: config_set_bool_c
        end function
        
        function config_has_section_c(section) bind(C, name="config_has_section_f")
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: section(*)
            integer(c_int) :: config_has_section_c
        end function
        
        function config_has_key_c(section, key) bind(C, name="config_has_key_f")
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: section(*)
            character(kind=c_char), intent(in) :: key(*)
            integer(c_int) :: config_has_key_c
        end function
    end interface
    
contains

    !
    ! Convert Fortran string to C string (null-terminated)
    !
    function to_c_string(f_string) result(c_string)
        character(len=*), intent(in) :: f_string
        character(len=len_trim(f_string)+1) :: c_string
        
        c_string = trim(f_string) // c_null_char
    end function
    
    !
    ! Load configuration from INI file
    !
    function config_load(filename) result(success)
        character(len=*), intent(in) :: filename
        logical :: success
        
        success = (config_load_c(to_c_string(filename)) /= 0)
    end function
    
    !
    ! Create new empty configuration
    !
    function config_new() result(success)
        logical :: success
        
        success = (config_new_c() /= 0)
    end function
    
    !
    ! Free configuration
    !
    subroutine config_free()
        call config_free_c()
    end subroutine
    
    !
    ! Save configuration to INI file
    !
    function config_save(filename) result(success)
        character(len=*), intent(in) :: filename
        logical :: success
        
        success = (config_save_c(to_c_string(filename)) /= 0)
    end function
    
    !
    ! Get string value
    !
    function config_get_string(section, key, default_val) result(value)
        character(len=*), intent(in) :: section, key, default_val
        character(len=256) :: value
        integer :: actual_len
        
        actual_len = config_get_string_c(to_c_string(section), &
                                         to_c_string(key), &
                                         to_c_string(default_val), &
                                         value, len(value))
        ! Value is already space-padded by C function
    end function
    
    !
    ! Get integer value
    !
    function config_get_int(section, key, default_val) result(value)
        character(len=*), intent(in) :: section, key
        integer, intent(in) :: default_val
        integer :: value
        
        value = config_get_int_c(to_c_string(section), &
                                 to_c_string(key), &
                                 default_val)
    end function
    
    !
    ! Get double precision value
    !
    function config_get_double(section, key, default_val) result(value)
        character(len=*), intent(in) :: section, key
        double precision, intent(in) :: default_val
        double precision :: value
        
        value = config_get_double_c(to_c_string(section), &
                                    to_c_string(key), &
                                    default_val)
    end function
    
    !
    ! Get boolean value
    !
    function config_get_bool(section, key, default_val) result(value)
        character(len=*), intent(in) :: section, key
        logical, intent(in) :: default_val
        logical :: value
        integer :: default_int, result_int
        
        if (default_val) then
            default_int = 1
        else
            default_int = 0
        end if
        
        result_int = config_get_bool_c(to_c_string(section), &
                                       to_c_string(key), &
                                       default_int)
        value = (result_int /= 0)
    end function
    
    !
    ! Set string value
    !
    function config_set_string(section, key, value) result(success)
        character(len=*), intent(in) :: section, key, value
        logical :: success
        
        success = (config_set_string_c(to_c_string(section), &
                                       to_c_string(key), &
                                       to_c_string(value)) /= 0)
    end function
    
    !
    ! Set integer value
    !
    function config_set_int(section, key, value) result(success)
        character(len=*), intent(in) :: section, key
        integer, intent(in) :: value
        logical :: success
        
        success = (config_set_int_c(to_c_string(section), &
                                    to_c_string(key), &
                                    value) /= 0)
    end function
    
    !
    ! Set double precision value
    !
    function config_set_double(section, key, value) result(success)
        character(len=*), intent(in) :: section, key
        double precision, intent(in) :: value
        logical :: success
        
        success = (config_set_double_c(to_c_string(section), &
                                       to_c_string(key), &
                                       value) /= 0)
    end function
    
    !
    ! Set boolean value
    !
    function config_set_bool(section, key, value) result(success)
        character(len=*), intent(in) :: section, key
        logical, intent(in) :: value
        logical :: success
        integer :: int_val
        
        if (value) then
            int_val = 1
        else
            int_val = 0
        end if
        
        success = (config_set_bool_c(to_c_string(section), &
                                     to_c_string(key), &
                                     int_val) /= 0)
    end function
    
    !
    ! Check if section exists
    !
    function config_has_section(section) result(exists)
        character(len=*), intent(in) :: section
        logical :: exists
        
        exists = (config_has_section_c(to_c_string(section)) /= 0)
    end function
    
    !
    ! Check if key exists in section
    !
    function config_has_key(section, key) result(exists)
        character(len=*), intent(in) :: section, key
        logical :: exists
        
        exists = (config_has_key_c(to_c_string(section), &
                                   to_c_string(key)) /= 0)
    end function

end module config_module
