program test_gw_geometry
  use metrics, only: event_horizon_radius
  use output, only: gw_cartesian_kinematics, kerr_schild_spheroidal_position
  implicit none

  integer :: failures
  real*8 :: r, theta, phi, spin, alpha, eps_fd
  real*8 :: beta(3), beta_log(3), eulerian_velocity(3), eulerian_velocity_log(3)
  real*8 :: rates(3), zero_velocity(3)
  real*8 :: position(3), position_log(3), velocity(3), velocity_log(3), jacobian
  real*8 :: position_plus(3), position_minus(3), velocity_dummy(3), mapped_position(3)
  real*8 :: jacobian_dummy, finite_difference_velocity(3)

  failures = 0
  zero_velocity = 0.0d0

  call check_close('Schwarzschild horizon', event_horizon_radius(1.0d0, 0.0d0), &
                   2.0d0, 1.0d-14, failures)
  call check_close('Kerr a=0.9 horizon', event_horizon_radius(1.0d0, 0.9d0), &
                   1.0d0 + sqrt(0.19d0), 1.0d-14, failures)
  call check_close('Extremal Kerr horizon', event_horizon_radius(1.0d0, 1.0d0), &
                   1.0d0, 1.0d-14, failures)

  r = 7.0d0
  theta = 1.1d0
  phi = 0.7d0
  spin = 0.8d0
  alpha = 0.8d0
  beta = [0.04d0, -0.02d0, 0.01d0]
  eulerian_velocity = beta / alpha
  call gw_cartesian_kinematics(r, theta, phi, spin, alpha, beta, &
                               eulerian_velocity, .false., position, velocity, jacobian)
  call check_vector('stationary physical-grid transport velocity', velocity, &
                    [0.0d0, 0.0d0, 0.0d0], 1.0d-14, failures)

  eulerian_velocity = [0.03d0, -0.04d0, 0.05d0]
  call gw_cartesian_kinematics(r, theta, phi, spin, alpha, beta, &
                               eulerian_velocity, .false., position, velocity, jacobian)
  beta_log = beta
  beta_log(1) = beta(1) / r
  eulerian_velocity_log = eulerian_velocity
  eulerian_velocity_log(1) = eulerian_velocity(1) / r
  call gw_cartesian_kinematics(r, theta, phi, spin, alpha, beta_log, &
                               eulerian_velocity_log, .true., position_log, &
                               velocity_log, jacobian_dummy)
  call check_vector('log-grid Cartesian position', position_log, position, &
                    1.0d-14, failures)
  call check_vector('log-grid transport velocity', velocity_log, velocity, &
                    1.0d-14, failures)

  r = 4.3d0
  theta = 1.1d0
  phi = 0.7d0
  spin = 0.8d0
  rates = [0.12d0, -0.03d0, 0.21d0]
  alpha = 1.0d0
  beta = 0.0d0
  call gw_cartesian_kinematics(r, theta, phi, spin, alpha, beta, rates, &
                               .false., position, velocity, jacobian)
  call kerr_schild_spheroidal_position(r, theta, phi, spin, mapped_position)
  call check_vector('VTK/GW spheroidal position consistency', mapped_position, &
                    position, 1.0d-14, failures)
  call check_close('KS flat Jacobian', jacobian, &
                   (r**2 + spin**2 * cos(theta)**2) * sin(theta), &
                   1.0d-13, failures)

  eps_fd = 1.0d-6
  call gw_cartesian_kinematics(r + eps_fd * rates(1), &
                               theta + eps_fd * rates(2), &
                               phi + eps_fd * rates(3), spin, alpha, beta, &
                               zero_velocity, .false., position_plus, &
                               velocity_dummy, jacobian_dummy)
  call gw_cartesian_kinematics(r - eps_fd * rates(1), &
                               theta - eps_fd * rates(2), &
                               phi - eps_fd * rates(3), spin, alpha, beta, &
                               zero_velocity, .false., position_minus, &
                               velocity_dummy, jacobian_dummy)
  finite_difference_velocity = (position_plus - position_minus) / (2.0d0 * eps_fd)
  call check_vector('KS Cartesian velocity chain rule', velocity, &
                    finite_difference_velocity, 2.0d-9, failures)

  call gw_cartesian_kinematics(r, theta, phi, 0.0d0, alpha, beta, zero_velocity, &
                               .false., position, velocity_dummy, jacobian)
  call check_vector('a=0 spherical Cartesian position', position, &
                    [r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)], &
                    1.0d-13, failures)
  call check_close('a=0 spherical Jacobian', jacobian, r**2*sin(theta), &
                   1.0d-13, failures)

  call kerr_schild_spheroidal_position(r, 0.0d0, phi, spin, mapped_position)
  call check_vector('exact north-pole collapse', mapped_position, &
                    [0.0d0, 0.0d0, r], 0.0d0, failures)
  call kerr_schild_spheroidal_position(r, acos(-1.0d0), phi, spin, mapped_position)
  call check_vector('exact south-pole collapse', mapped_position, &
                    [0.0d0, 0.0d0, -r], 0.0d0, failures)

  if (failures /= 0) then
    write(*,*) 'GW geometry tests FAILED: ', failures
    error stop 1
  end if
  write(*,*) 'GW geometry tests PASSED.'

contains

  subroutine check_close(label, value, expected, tolerance, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: value, expected, tolerance
    integer, intent(inout) :: failure_count

    if (abs(value - expected) > tolerance) then
      write(*,*) 'FAIL: ', trim(label), ' value=', value, ' expected=', expected
      failure_count = failure_count + 1
    end if
  end subroutine check_close

  subroutine check_vector(label, value, expected, tolerance, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: value(3), expected(3), tolerance
    integer, intent(inout) :: failure_count

    if (maxval(abs(value - expected)) > tolerance) then
      write(*,*) 'FAIL: ', trim(label)
      write(*,*) ' value=', value
      write(*,*) ' expected=', expected
      failure_count = failure_count + 1
    end if
  end subroutine check_vector

end program test_gw_geometry
