program test_hllc_fluxes
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use variables
  use equations
  use fluxes
  implicit none

  integer :: failures

  failures = 0
  adb_idx = 4.0d0 / 3.0d0
  g1 = adb_idx / (adb_idx - 1.0d0)
  rho_floor = 1.0d-10
  p_floor = 1.0d-13
  riemann_solver_id = RS_HLLC

  call test_uniform_cartesian(failures)
  call test_uniform_general_metric(failures)
  call test_stationary_contact(failures)
  call test_moving_contact(failures)
  call test_shifted_contact(failures)
  call test_tangential_contact(failures)
  call test_general_metric_contact(failures)
  call test_general_metric_moving_contact(failures)
  call test_atmosphere_scale_contact(failures)
  call test_scale_invariance(failures)
  call test_strong_relativistic_jump(failures)

  if (failures /= 0) then
    write(*,'(A,I0)') 'HLLC TESTS FAILED: ', failures
    error stop 1
  end if
  write(*,'(A)') 'HLLC TESTS PASSED'

contains

  subroutine configure_metric(metric_name, geometry_name)
    character(len=*), intent(in) :: metric_name, geometry_name
    metric_type = metric_name
    geom_type = geometry_name
    call init_wavespeed_solver()
  end subroutine configure_metric


  subroutine set_identity_metric(g)
    real*8, intent(out) :: g(3,3)
    g = 0.0d0
    g(1,1) = 1.0d0
    g(2,2) = 1.0d0
    g(3,3) = 1.0d0
  end subroutine set_identity_metric


  real*8 function determinant_3x3(g)
    real*8, intent(in) :: g(3,3)
    determinant_3x3 = g(1,1)*(g(2,2)*g(3,3) - g(2,3)*g(3,2)) &
                    - g(1,2)*(g(2,1)*g(3,3) - g(2,3)*g(3,1)) &
                    + g(1,3)*(g(2,1)*g(3,2) - g(2,2)*g(3,1))
  end function determinant_3x3


  subroutine assert_close(label, actual, expected, tolerance, failures)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: actual(:), expected(:), tolerance
    integer, intent(inout) :: failures
    real*8 :: error, scale

    scale = max(1.0d0, maxval(abs(expected)))
    error = maxval(abs(actual - expected)) / scale
    if (.not. all(ieee_is_finite(actual)) .or. error > tolerance) then
      failures = failures + 1
      write(*,'(A,1X,A,1X,ES12.4)') 'FAIL', trim(label), error
      write(*,'(A,5(1X,ES14.6))') ' actual :', actual
      write(*,'(A,5(1X,ES14.6))') ' expected:', expected
    else
      write(*,'(A,1X,A,1X,ES12.4)') 'PASS', trim(label), error
    end if
  end subroutine assert_close


  subroutine assert_finite(label, values, failures)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: values(:)
    integer, intent(inout) :: failures

    if (.not. all(ieee_is_finite(values))) then
      failures = failures + 1
      write(*,'(A,1X,A)') 'FAIL', trim(label)
    else
      write(*,'(A,1X,A)') 'PASS', trim(label)
    end if
  end subroutine assert_finite


  subroutine test_uniform_cartesian(failures)
    integer, intent(inout) :: failures
    integer :: direction
    real*8 :: q(neq), g(3,3), beta(3), numerical(neq), exact(neq)

    call configure_metric('Minkowski', 'Cartesian')
    call set_identity_metric(g)
    beta = 0.0d0
    q = (/ 1.2d0, 0.25d0, 0.12d0, -0.07d0, 0.05d0 /)

    do direction = DIR_X, DIR_Z
      call resolve_riemann_problem(q, q, direction, 1.0d0, beta, g, 1.0d0, numerical)
      call calc_fluxes(q, exact, direction, 1.0d0, beta, g, 1.0d0)
      call assert_close('uniform Cartesian direction', numerical, exact, 2.0d-12, failures)
    end do
  end subroutine test_uniform_cartesian


  subroutine test_uniform_general_metric(failures)
    integer, intent(inout) :: failures
    integer :: direction
    real*8 :: q(neq), g(3,3), beta(3), alpha, sqg
    real*8 :: numerical(neq), exact(neq)

    call configure_metric('Kerr-Schild', 'Spherical')
    g = reshape((/ 2.0d0, 0.08d0, -0.35d0, &
                   0.08d0, 1.5d0, 0.10d0, &
                  -0.35d0, 0.10d0, 3.0d0 /), shape(g))
    alpha = 0.73d0
    beta = (/ 0.18d0, -0.05d0, 0.09d0 /)
    sqg = sqrt(determinant_3x3(g))
    q = (/ 0.9d0, 0.12d0, 0.08d0, -0.04d0, 0.06d0 /)

    do direction = DIR_X, DIR_Z
      call resolve_riemann_problem(q, q, direction, alpha, beta, g, sqg, numerical)
      call calc_fluxes(q, exact, direction, alpha, beta, g, sqg)
      call assert_close('uniform general metric direction', numerical, exact, 5.0d-12, failures)
    end do
  end subroutine test_uniform_general_metric


  subroutine test_stationary_contact(failures)
    integer, intent(inout) :: failures
    real*8 :: q_L(neq), q_R(neq), g(3,3), beta(3), numerical(neq), exact(neq)

    call configure_metric('Minkowski', 'Cartesian')
    call set_identity_metric(g)
    beta = 0.0d0
    q_L = (/ 1.0d0, 0.4d0, 0.0d0, 0.0d0, 0.0d0 /)
    q_R = (/ 0.125d0, 0.4d0, 0.0d0, 0.0d0, 0.0d0 /)

    call resolve_riemann_problem(q_L, q_R, DIR_X, 1.0d0, beta, g, 1.0d0, numerical)
    call calc_fluxes(q_L, exact, DIR_X, 1.0d0, beta, g, 1.0d0)
    call assert_close('stationary density contact', numerical, exact, 2.0d-12, failures)
  end subroutine test_stationary_contact


  subroutine test_moving_contact(failures)
    integer, intent(inout) :: failures
    real*8 :: q_L(neq), q_R(neq), g(3,3), beta(3), numerical(neq), exact(neq)

    call configure_metric('Minkowski', 'Cartesian')
    call set_identity_metric(g)
    beta = 0.0d0
    q_L = (/ 1.0d0, 0.1d0, 0.2d0, 0.03d0, -0.02d0 /)
    q_R = (/ 0.2d0, 0.1d0, 0.2d0, 0.03d0, -0.02d0 /)

    call resolve_riemann_problem(q_L, q_R, DIR_X, 1.0d0, beta, g, 1.0d0, numerical)
    call calc_fluxes(q_L, exact, DIR_X, 1.0d0, beta, g, 1.0d0)
    call assert_close('right-moving contact selects left state', numerical, exact, 3.0d-12, failures)
  end subroutine test_moving_contact


  subroutine test_shifted_contact(failures)
    integer, intent(inout) :: failures
    real*8 :: q_L(neq), q_R(neq), g(3,3), beta(3), numerical(neq), exact(neq)

    call configure_metric('Kerr-Schild', 'Spherical')
    call set_identity_metric(g)
    beta = (/ 0.30d0, 0.0d0, 0.0d0 /)
    q_L = (/ 1.0d0, 0.1d0, 0.2d0, 0.03d0, -0.02d0 /)
    q_R = (/ 0.2d0, 0.1d0, 0.2d0, 0.03d0, -0.02d0 /)

    call resolve_riemann_problem(q_L, q_R, DIR_X, 0.8d0, beta, g, 1.0d0, numerical)
    call calc_fluxes(q_R, exact, DIR_X, 0.8d0, beta, g, 1.0d0)
    call assert_close('shifted contact selects right state', numerical, exact, 3.0d-12, failures)
  end subroutine test_shifted_contact


  subroutine test_tangential_contact(failures)
    integer, intent(inout) :: failures
    real*8 :: q_L(neq), q_R(neq), g(3,3), beta(3), numerical(neq), exact(neq)

    call configure_metric('Minkowski', 'Cartesian')
    call set_identity_metric(g)
    beta = 0.0d0
    q_L = (/ 1.0d0, 0.2d0, 0.0d0, 0.30d0, -0.10d0 /)
    q_R = (/ 0.4d0, 0.2d0, 0.0d0, -0.20d0, 0.15d0 /)

    call resolve_riemann_problem(q_L, q_R, DIR_X, 1.0d0, beta, g, 1.0d0, numerical)
    exact = (/ 0.0d0, 0.0d0, 0.2d0, 0.0d0, 0.0d0 /)
    call assert_close('stationary tangential contact', numerical, exact, 3.0d-12, failures)
  end subroutine test_tangential_contact


  subroutine test_general_metric_contact(failures)
    integer, intent(inout) :: failures
    real*8 :: q_L(neq), q_R(neq), g(3,3), beta(3), numerical(neq), exact(neq), sqg

    call configure_metric('Kerr-Schild', 'Spherical')
    g = 0.0d0
    g(1,1) = 2.2d0
    g(2,2) = 1.4d0
    g(3,3) = 3.1d0
    g(1,3) = -0.45d0
    g(3,1) = g(1,3)
    beta = 0.0d0
    sqg = sqrt(determinant_3x3(g))
    q_L = (/ 1.0d0, 0.15d0, 0.0d0, 0.04d0, 0.08d0 /)
    q_R = (/ 0.3d0, 0.15d0, 0.0d0, -0.03d0, -0.05d0 /)

    call resolve_riemann_problem(q_L, q_R, DIR_X, 0.75d0, beta, g, sqg, numerical)
    call calc_fluxes(q_L, exact, DIR_X, 0.75d0, beta, g, sqg)
    call assert_close('general-metric stationary contact', numerical, exact, 5.0d-12, failures)
  end subroutine test_general_metric_contact


  subroutine test_general_metric_moving_contact(failures)
    integer, intent(inout) :: failures
    real*8 :: q_L(neq), q_R(neq), g(3,3), beta(3), numerical(neq), exact(neq), sqg

    call configure_metric('Kerr-Schild', 'Spherical')
    g = reshape((/ 2.0d0, 0.08d0, -0.35d0, &
                   0.08d0, 1.5d0, 0.10d0, &
                  -0.35d0, 0.10d0, 3.0d0 /), shape(g))
    beta = (/ 0.02d0, -0.01d0, 0.03d0 /)
    sqg = sqrt(determinant_3x3(g))
    q_L = (/ 1.1d0, 0.12d0, 0.08d0, -0.04d0, 0.06d0 /)
    q_R = (/ 0.25d0, 0.12d0, 0.08d0, -0.04d0, 0.06d0 /)

    call resolve_riemann_problem(q_L, q_R, DIR_X, 0.8d0, beta, g, sqg, numerical)
    call calc_fluxes(q_L, exact, DIR_X, 0.8d0, beta, g, sqg)
    call assert_close('general-metric moving contact', numerical, exact, 8.0d-12, failures)
  end subroutine test_general_metric_moving_contact


  subroutine test_atmosphere_scale_contact(failures)
    integer, intent(inout) :: failures
    real*8 :: q_L(neq), q_R(neq), g(3,3), beta(3), numerical(neq), exact(neq)

    call configure_metric('Minkowski', 'Cartesian')
    call set_identity_metric(g)
    beta = 0.0d0
    q_L = (/ 4.0d-10, 2.0d-13, 0.12d0, 0.0d0, 0.0d0 /)
    q_R = (/ 1.2d-10, 2.0d-13, 0.12d0, 0.0d0, 0.0d0 /)

    call resolve_riemann_problem(q_L, q_R, DIR_X, 1.0d0, beta, g, 1.0d0, numerical)
    call calc_fluxes(q_L, exact, DIR_X, 1.0d0, beta, g, 1.0d0)
    call assert_close('atmosphere-scale moving contact', numerical, exact, 1.0d-18, failures)
  end subroutine test_atmosphere_scale_contact


  subroutine test_scale_invariance(failures)
    integer, intent(inout) :: failures
    real*8, parameter :: scale_factor = 1.0d-10
    real*8 :: q_L(neq), q_R(neq), q_L_small(neq), q_R_small(neq)
    real*8 :: g(3,3), beta(3), flux_reference(neq), flux_small(neq)

    call configure_metric('Minkowski', 'Cartesian')
    call set_identity_metric(g)
    beta = 0.0d0
    q_L = (/ 1.0d0, 1.0d0, 0.30d0, 0.10d0, 0.0d0 /)
    q_R = (/ 0.125d0, 0.10d0, -0.10d0, -0.05d0, 0.02d0 /)
    q_L_small = q_L
    q_R_small = q_R
    q_L_small(eq_de:eq_pr) = scale_factor * q_L_small(eq_de:eq_pr)
    q_R_small(eq_de:eq_pr) = scale_factor * q_R_small(eq_de:eq_pr)

    call resolve_riemann_problem(q_L, q_R, DIR_X, 1.0d0, beta, g, 1.0d0, flux_reference)
    call resolve_riemann_problem(q_L_small, q_R_small, DIR_X, 1.0d0, beta, g, 1.0d0, flux_small)
    call assert_close('scale-invariant HLLC flux', flux_small / scale_factor, &
                      flux_reference, 2.0d-11, failures)
  end subroutine test_scale_invariance


  subroutine test_strong_relativistic_jump(failures)
    integer, intent(inout) :: failures
    real*8 :: q_L(neq), q_R(neq), g(3,3), beta(3), numerical(neq)

    call configure_metric('Minkowski', 'Cartesian')
    call set_identity_metric(g)
    beta = 0.0d0
    q_L = (/ 1.0d0, 1.0d3, 0.0d0, 0.0d0, 0.0d0 /)
    q_R = (/ 1.0d0, 1.0d-2, 0.0d0, 0.0d0, 0.0d0 /)

    call resolve_riemann_problem(q_L, q_R, DIR_X, 1.0d0, beta, g, 1.0d0, numerical)
    call assert_finite('strong relativistic pressure jump', numerical, failures)
  end subroutine test_strong_relativistic_jump

end program test_hllc_fluxes
