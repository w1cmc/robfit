!
! test_config.f90 - Test program for config_module
!
program test_config
    use config_module
    implicit none
    
    character(len=256) :: spectrum_file
    double precision :: energy_begin, energy_end
    integer :: poly_order
    logical :: view_data
    logical :: ok
    
    print *, "Testing config_module..."
    print *
    
    ! Try to load config
    if (config_load('robfit.ini')) then
        print *, "Loaded robfit.ini successfully"
    else
        print *, "Could not load robfit.ini, creating new config"
        ok = config_new()
        ok = config_set_string('spectrum', 'file', 'TEST.+ML')
        ok = config_set_double('display', 'energy_begin', 0.0d0)
        ok = config_set_double('display', 'energy_end', 3000.0d0)
        ok = config_set_int('fitting', 'poly_order', 5)
        ok = config_set_bool('display', 'view_data', .true.)
        ok = config_save('test_config.ini')
        print *, "Created test_config.ini"
        
        ! Reload it
        ok = config_load('test_config.ini')
    end if
    
    print *
    
    ! Read values
    spectrum_file = config_get_string('spectrum', 'file', 'NONE')
    print *, "spectrum.file = ", trim(spectrum_file)
    
    energy_begin = config_get_double('display', 'energy_begin', -1.0d0)
    print *, "display.energy_begin = ", energy_begin
    
    energy_end = config_get_double('display', 'energy_end', -1.0d0)
    print *, "display.energy_end = ", energy_end
    
    poly_order = config_get_int('fitting', 'poly_order', -1)
    print *, "fitting.poly_order = ", poly_order
    
    view_data = config_get_bool('display', 'view_data', .false.)
    print *, "display.view_data = ", view_data
    
    ! Test has_section / has_key
    print *
    if (config_has_section('spectrum')) then
        print *, "[spectrum] section exists"
    else
        print *, "[spectrum] section NOT found"
    end if
    
    if (config_has_key('spectrum', 'file')) then
        print *, "spectrum.file key exists"
    else
        print *, "spectrum.file key NOT found"
    end if
    
    if (config_has_key('spectrum', 'nonexistent')) then
        print *, "spectrum.nonexistent key exists (unexpected!)"
    else
        print *, "spectrum.nonexistent key NOT found (expected)"
    end if
    
    ! Clean up
    call config_free()
    
    print *
    print *, "Test complete."
    
end program test_config
