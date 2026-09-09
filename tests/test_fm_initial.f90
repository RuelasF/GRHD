program test_fm_initial
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use variables
  use metrics
  use conditions
  implicit none

  integer :: failures

  failures = 0
  adb_idx = 4.0d0 / 3.0d0
  g1 = adb_idx / (adb_idx - 1.0d0)
  rho_floor = 1.0d-10
  p_floor = 1.0d-13
  bh_mass = 1.0d0
  geom_type = 'Spherical'
  use_log_r = .false.

  call check_model('Eddington-Finkelstein', 0.0d0, 6.0d0, 12.0d0, 0.015d0, failures)
  call check_model('Kerr-Schild', 0.9d0, 6.25d0, 9.25d0, 0.0015d0, failures)

  if (failures /= 0) then
    write(*,'(A,I0)') 'FISHBONE-MONCRIEF TESTS FAILED: ', failures
    error stop 1
  end if
  write(*,'(A)') 'FISHBONE-MONCRIEF TESTS PASSED'

contains

  subroutine check_model(metric_name, spin, r_in, r_center, K_poly, failures)
    character(len=*), intent(in) :: metric_name
    real*8, intent(in) :: spin, r_in, r_center, K_poly
    integer, intent(inout) :: failures
    real*8 :: l_ang, W_in, alpha, beta(3), g(3,3)
    real*8 :: state_center(neq), state_left(neq), state_right(neq)
    real*8 :: state_north(neq), state_south(neq), state_atmosphere(neq)
    real*8 :: radial_offset
    logical :: ok, inside

    metric_type = metric_name
    a_spin = spin
    call set_metric_type()
    call prepare_fishbone_moncrief(r_in, r_center, l_ang, W_in, ok)
    call assert_true(trim(metric_name)//' preparation', ok, failures)
    if (.not. ok) return

    call calculate_metric(r_center, pi/2.0d0, alpha=alpha, beta=beta, gamma=g)
    call evaluate_fishbone_moncrief_state(r_center, r_in, K_poly, l_ang, W_in, &
                                          alpha, beta, g, state_center, inside)
    call assert_true(trim(metric_name)//' center belongs to torus', inside, failures)
    call assert_true(trim(metric_name)//' center is finite', &
                     all(ieee_is_finite(state_center)), failures)
    call assert_true(trim(metric_name)//' positive thermodynamics', &
                     state_center(eq_de) > rho_floor .and. &
                     state_center(eq_pr) > p_floor, failures)
    call check_specific_angular_momentum(metric_name, state_center, alpha, beta, g, &
                                         l_ang, failures)

    radial_offset = 1.0d-3 * r_center
    call evaluate_at_point(r_center-radial_offset, pi/2.0d0, r_in, K_poly, &
                           l_ang, W_in, state_left, inside)
    call assert_true(trim(metric_name)//' inner neighbor belongs to torus', inside, failures)
    call evaluate_at_point(r_center+radial_offset, pi/2.0d0, r_in, K_poly, &
                           l_ang, W_in, state_right, inside)
    call assert_true(trim(metric_name)//' outer neighbor belongs to torus', inside, failures)
    call assert_true(trim(metric_name)//' pressure maximum radius', &
                     state_center(eq_pr) >= state_left(eq_pr) .and. &
                     state_center(eq_pr) >= state_right(eq_pr), failures)

    call evaluate_at_point(r_center, pi/2.0d0-0.2d0, r_in, K_poly, &
                           l_ang, W_in, state_north, inside)
    call assert_true(trim(metric_name)//' north point belongs to torus', inside, failures)
    call evaluate_at_point(r_center, pi/2.0d0+0.2d0, r_in, K_poly, &
                           l_ang, W_in, state_south, inside)
    call assert_true(trim(metric_name)//' south point belongs to torus', inside, failures)
    call assert_close(trim(metric_name)//' equatorial symmetry', state_north, state_south, &
                      2.0d-12, failures)

    call evaluate_at_point(0.9d0*r_in, pi/2.0d0, r_in, K_poly, l_ang, W_in, &
                           state_atmosphere, inside)
    call assert_true(trim(metric_name)//' inner atmosphere classification', .not. inside, failures)
    call assert_true(trim(metric_name)//' inner atmosphere floors', &
                     state_atmosphere(eq_de) == rho_floor .and. &
                     state_atmosphere(eq_pr) == p_floor, failures)
  end subroutine check_model


  subroutine evaluate_at_point(radius, theta, r_in, K_poly, l_ang, W_in, state, inside)
    real*8, intent(in) :: radius, theta, r_in, K_poly, l_ang, W_in
    real*8, intent(out) :: state(neq)
    logical, intent(out) :: inside
    real*8 :: alpha, beta(3), g(3,3)

    call calculate_metric(radius, theta, alpha=alpha, beta=beta, gamma=g)
    call evaluate_fishbone_moncrief_state(radius, r_in, K_poly, l_ang, W_in, &
                                          alpha, beta, g, state, inside)
  end subroutine evaluate_at_point


  subroutine check_specific_angular_momentum(label, state, alpha, beta, g, expected, failures)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: state(neq), alpha, beta(3), g(3,3), expected
    integer, intent(inout) :: failures
    real*8 :: velocity(3), shift_cov(3), u_con(0:3)
    real*8 :: v_sq, lorentz_factor, g_tt, u_t, u_phi, actual

    velocity = state(eq_vx:eq_vz)
    v_sq = dot_product(velocity, matmul(g, velocity))
    lorentz_factor = 1.0d0 / sqrt(1.0d0-v_sq)
    shift_cov = matmul(g, beta)
    g_tt = -alpha**2 + dot_product(beta, shift_cov)
    u_con(0) = lorentz_factor / alpha
    u_con(1:3) = lorentz_factor * (velocity-beta/alpha)
    u_t = g_tt*u_con(0) + dot_product(shift_cov, u_con(1:3))
    u_phi = shift_cov(3)*u_con(0) + dot_product(g(3,:), u_con(1:3))
    actual = -u_phi/u_t
    call assert_scalar(trim(label)//' constant angular momentum', actual, expected, &
                       2.0d-12, failures)
  end subroutine check_specific_angular_momentum


  subroutine assert_true(label, condition, failures)
    character(len=*), intent(in) :: label
    logical, intent(in) :: condition
    integer, intent(inout) :: failures

    if (.not. condition) then
      failures = failures + 1
      write(*,'(A,1X,A)') 'FAIL', trim(label)
    end if
  end subroutine assert_true


  subroutine assert_scalar(label, actual, expected, tolerance, failures)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: actual, expected, tolerance
    integer, intent(inout) :: failures

    if (.not. ieee_is_finite(actual) .or. &
        abs(actual-expected) > tolerance*max(1.0d0,abs(expected))) then
      failures = failures + 1
      write(*,'(A,1X,A,2(1X,ES14.6))') 'FAIL', trim(label), actual, expected
    end if
  end subroutine assert_scalar


  subroutine assert_close(label, actual, expected, tolerance, failures)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: actual(:), expected(:), tolerance
    integer, intent(inout) :: failures

    if (.not. all(ieee_is_finite(actual)) .or. &
        maxval(abs(actual-expected)) > tolerance*max(1.0d0,maxval(abs(expected)))) then
      failures = failures + 1
      write(*,'(A,1X,A,1X,ES14.6)') 'FAIL', trim(label), maxval(abs(actual-expected))
    end if
  end subroutine assert_close

end program test_fm_initial
