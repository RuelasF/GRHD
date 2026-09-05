program test_gw_stress
  use output, only: finn_evans_cell_stress
  implicit none

  real*8, parameter :: tolerance = 1.0d-13
  real*8 :: mass, particle_mass, radius, phase, omega
  real*8 :: position(3), velocity(3)
  real*8 :: Ixx_cell, Iyy_cell, Ixy_cell
  real*8 :: Ixx_expected, Iyy_expected, Ixy_expected
  integer :: failures

  failures = 0
  mass = 1.0d0
  particle_mass = 0.2d0
  radius = 5.0d0
  phase = 0.37d0
  omega = sqrt(mass / radius**3)

  position = [radius*cos(phase), radius*sin(phase), 0.0d0]
  velocity = [-radius*omega*sin(phase), radius*omega*cos(phase), 0.0d0]

  call finn_evans_cell_stress(particle_mass, 0.0d0, mass/radius**3, &
                              position, velocity, 1.0d0, Ixx_cell, &
                              Iyy_cell, Ixy_cell)

  ! Segunda derivada analítica del cuadrupolo STF de una masa puntual en
  ! órbita circular. La traza no altera las polarizaciones face-on.
  Ixx_expected = -2.0d0 * particle_mass * radius**2 * omega**2 * cos(2.0d0*phase)
  Iyy_expected =  2.0d0 * particle_mass * radius**2 * omega**2 * cos(2.0d0*phase)
  Ixy_expected = -2.0d0 * particle_mass * radius**2 * omega**2 * sin(2.0d0*phase)

  call check_close('Ixx circular orbit', Ixx_cell, Ixx_expected, failures)
  call check_close('Iyy circular orbit', Iyy_cell, Iyy_expected, failures)
  call check_close('Ixy circular orbit', Ixy_cell, Ixy_expected, failures)

  if (failures /= 0) then
    write(*,*) 'Finn-Evans stress tests FAILED: ', failures
    error stop 1
  end if
  write(*,*) 'Finn-Evans stress tests PASSED.'

contains

  subroutine check_close(label, value, expected, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: value, expected
    integer, intent(inout) :: failure_count

    if (abs(value - expected) > tolerance) then
      write(*,*) 'FAIL: ', trim(label), ' value=', value, ' expected=', expected
      failure_count = failure_count + 1
    end if
  end subroutine check_close

end program test_gw_stress
