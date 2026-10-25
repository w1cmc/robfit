!
! test_nuclides.f90 - Test program for nuclide library
!
program test_nuclides
  use nuclides_module
  implicit none

  integer :: i, n, nlines
  character(len=32) :: name
  real(8) :: energies(500), intensities(500)
  character(len=1) :: types(500)

  print *, '=== Nuclide Library Test ==='
  print *

  ! Load the library
  if (.not. nuclides_load('nuclides.xml')) then
    print *, 'ERROR: Failed to load nuclides.xml'
    stop 1
  end if

  ! Get count
  n = nuclides_count()
  print '(A,I3)', ' Loaded nuclides: ', n

  ! List all nuclides
  print *
  print *, 'Nuclides in library:'
  do i = 1, n
    nlines = nuclides_get_name(i, name)
    print '(A,I3,A,A,A,I4,A)', '  ', i, ': ', trim(name), ' (', nlines, ' lines)'
  end do

  ! Test finding a nuclide
  print *
  i = nuclides_find_index('Co60')
  if (i > 0) then
    print '(A,I3)', ' Found Co60 at index: ', i
  else
    print *, ' Co60 not found'
  end if

  ! Get emission lines for Co60
  print *
  print *, 'Co60 emission lines (100-2000 keV):'
  nlines = nuclides_get_lines('Co60', 100.0d0, 2000.0d0, &
                               energies, intensities, types, 500)
  print '(A,I3)', ' Found lines: ', nlines
  do i = 1, min(nlines, 10)
    print '(A,F10.3,A,E12.5,A,A1)', '   E=', energies(i), ' keV, I=', &
          intensities(i), ', Type=', types(i)
  end do

  ! Get lines for Cs137
  print *
  print *, 'Cs137 emission lines (0-1000 keV):'
  nlines = nuclides_get_lines('Cs137', 0.0d0, 1000.0d0, &
                               energies, intensities, types, 500)
  print '(A,I3)', ' Found lines: ', nlines
  do i = 1, min(nlines, 10)
    print '(A,F10.3,A,E12.5,A,A1)', '   E=', energies(i), ' keV, I=', &
          intensities(i), ', Type=', types(i)
  end do

  ! Clean up
  call nuclides_free()
  print *
  print *, 'Test completed successfully!'

end program test_nuclides
