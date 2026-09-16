program test_minkowski_spherical
  use variables
  use metrics
  implicit none

  integer :: failures
  real*8, parameter :: tolerance = 2.0d-13, curvature_tolerance = 2.0d-8
  real*8 :: r, theta, xlog
  real*8 :: alpha, determinant
  real*8 :: beta(3), gamma(3,3), gmunu(0:3,0:3), dlnalpha(0:3)
  real*8 :: christoffel(0:3,0:3,0:3), derivatives(0:3,0:3,1:3)
  real*8 :: alpha_reference, determinant_reference
  real*8 :: beta_reference(3), gamma_reference(3,3)
  real*8 :: gmunu_reference(0:3,0:3), dlnalpha_reference(0:3)
  real*8 :: christoffel_reference(0:3,0:3,0:3)
  real*8 :: derivatives_reference(0:3,0:3,1:3)

  failures = 0
  r = 3.7d0
  theta = 0.91d0
  xlog = log(r)
  ny = 8
  bh_mass = 0.0d0
  a_spin = 0.0d0
  geom_type = 'Spherical'

  metric_type = 'Minkowski'
  use_log_r = .false.
  call set_metric_type()
  call calculate_metric(r, theta, alpha, beta, gamma, gmunu, determinant, dlnalpha)
  call calculate_christoffel_symbols(r, theta, christoffel)
  call calculate_metric_derivatives(r, theta, derivatives)
  call check_physical_form(r, theta, alpha, beta, gamma, gmunu, determinant, &
                           dlnalpha, christoffel, derivatives, failures)
  call check_flat_riemann(r, theta, .false., curvature_tolerance, failures)

  alpha_reference = alpha
  beta_reference = beta
  gamma_reference = gamma
  gmunu_reference = gmunu
  determinant_reference = determinant
  dlnalpha_reference = dlnalpha
  christoffel_reference = christoffel
  derivatives_reference = derivatives

  metric_type = 'Eddington-Finkelstein'
  call set_metric_type()
  call calculate_metric(r, theta, alpha, beta, gamma, gmunu, determinant, dlnalpha)
  call calculate_christoffel_symbols(r, theta, christoffel)
  call calculate_metric_derivatives(r, theta, derivatives)
  call check_close('physical M=0 EF lapse', alpha, alpha_reference, tolerance, failures)
  call check_vector('physical M=0 EF shift', beta, beta_reference, tolerance, failures)
  call check_matrix3('physical M=0 EF spatial metric', gamma, gamma_reference, tolerance, failures)
  call check_matrix4('physical M=0 EF inverse metric', gmunu, gmunu_reference, tolerance, failures)
  call check_close('physical M=0 EF determinant', determinant, determinant_reference, tolerance, failures)
  call check_vector4('physical M=0 EF lapse derivative', dlnalpha, dlnalpha_reference, tolerance, failures)
  call check_tensor3('physical M=0 EF Christoffel', christoffel, christoffel_reference, tolerance, failures)
  call check_derivatives('physical M=0 EF metric derivatives', derivatives, derivatives_reference, tolerance, failures)

  metric_type = 'Minkowski'
  use_log_r = .true.
  call set_metric_type()
  call calculate_metric(xlog, theta, alpha, beta, gamma, gmunu, determinant, dlnalpha)
  call calculate_christoffel_symbols(xlog, theta, christoffel)
  call calculate_metric_derivatives(xlog, theta, derivatives)
  call check_log_form(r, theta, alpha, beta, gamma, gmunu, determinant, &
                      dlnalpha, christoffel, derivatives, failures)
  call check_flat_riemann(xlog, theta, .true., curvature_tolerance, failures)

  alpha_reference = alpha
  beta_reference = beta
  gamma_reference = gamma
  gmunu_reference = gmunu
  determinant_reference = determinant
  dlnalpha_reference = dlnalpha
  christoffel_reference = christoffel
  derivatives_reference = derivatives

  metric_type = 'Eddington-Finkelstein'
  call set_metric_type()
  call calculate_metric(xlog, theta, alpha, beta, gamma, gmunu, determinant, dlnalpha)
  call calculate_christoffel_symbols(xlog, theta, christoffel)
  call calculate_metric_derivatives(xlog, theta, derivatives)
  call check_close('log M=0 EF lapse', alpha, alpha_reference, tolerance, failures)
  call check_vector('log M=0 EF shift', beta, beta_reference, tolerance, failures)
  call check_matrix3('log M=0 EF spatial metric', gamma, gamma_reference, tolerance, failures)
  call check_matrix4('log M=0 EF inverse metric', gmunu, gmunu_reference, tolerance, failures)
  call check_close('log M=0 EF determinant', determinant, determinant_reference, tolerance, failures)
  call check_vector4('log M=0 EF lapse derivative', dlnalpha, dlnalpha_reference, tolerance, failures)
  call check_tensor3('log M=0 EF Christoffel', christoffel, christoffel_reference, tolerance, failures)
  call check_derivatives('log M=0 EF metric derivatives', derivatives, derivatives_reference, tolerance, failures)

  if (failures /= 0) then
    write(*,*) 'MINKOWSKI SPHERICAL TESTS FAILED: ', failures
    error stop 1
  end if
  write(*,*) 'MINKOWSKI SPHERICAL TESTS PASSED'

contains

  subroutine check_physical_form(radius, angle, lapse, shift, spatial_metric, inverse_metric, &
                                 det_gamma, lapse_derivative, connection, metric_derivative, failure_count)
    real*8, intent(in) :: radius, angle, lapse, shift(3), spatial_metric(3,3)
    real*8, intent(in) :: inverse_metric(0:3,0:3), det_gamma, lapse_derivative(0:3)
    real*8, intent(in) :: connection(0:3,0:3,0:3), metric_derivative(0:3,0:3,1:3)
    integer, intent(inout) :: failure_count
    real*8 :: expected3(3,3), expected4(0:3,0:3), expected_connection(0:3,0:3,0:3)
    real*8 :: expected_derivative(0:3,0:3,1:3), sin_angle, cos_angle

    sin_angle = sin(angle)
    cos_angle = cos(angle)
    expected3 = 0.0d0
    expected3(1,1) = 1.0d0
    expected3(2,2) = radius**2
    expected3(3,3) = radius**2*sin_angle**2
    expected4 = 0.0d0
    expected4(0,0) = -1.0d0
    expected4(1,1) = 1.0d0
    expected4(2,2) = 1.0d0/radius**2
    expected4(3,3) = 1.0d0/(radius**2*sin_angle**2)
    expected_connection = 0.0d0
    expected_connection(1,2,2) = -radius
    expected_connection(1,3,3) = -radius*sin_angle**2
    expected_connection(2,1,2) = 1.0d0/radius
    expected_connection(2,2,1) = 1.0d0/radius
    expected_connection(2,3,3) = -sin_angle*cos_angle
    expected_connection(3,1,3) = 1.0d0/radius
    expected_connection(3,3,1) = 1.0d0/radius
    expected_connection(3,2,3) = cos_angle/sin_angle
    expected_connection(3,3,2) = cos_angle/sin_angle
    expected_derivative = 0.0d0
    expected_derivative(2,2,1) = 2.0d0*radius
    expected_derivative(3,3,1) = 2.0d0*radius*sin_angle**2
    expected_derivative(3,3,2) = 2.0d0*radius**2*sin_angle*cos_angle

    call check_close('physical lapse', lapse, 1.0d0, tolerance, failure_count)
    call check_vector('physical shift', shift, [0.0d0,0.0d0,0.0d0], tolerance, failure_count)
    call check_matrix3('physical spatial metric', spatial_metric, expected3, tolerance, failure_count)
    call check_matrix4('physical inverse metric', inverse_metric, expected4, tolerance, failure_count)
    call check_close('physical determinant', det_gamma, radius**4*sin_angle**2, tolerance, failure_count)
    call check_vector4('physical lapse derivative', lapse_derivative, [0.0d0,0.0d0,0.0d0,0.0d0], &
                       tolerance, failure_count)
    call check_tensor3('physical Christoffel', connection, expected_connection, tolerance, failure_count)
    call check_derivatives('physical metric derivatives', metric_derivative, expected_derivative, tolerance, failure_count)
  end subroutine check_physical_form


  subroutine check_log_form(radius, angle, lapse, shift, spatial_metric, inverse_metric, &
                            det_gamma, lapse_derivative, connection, metric_derivative, failure_count)
    real*8, intent(in) :: radius, angle, lapse, shift(3), spatial_metric(3,3)
    real*8, intent(in) :: inverse_metric(0:3,0:3), det_gamma, lapse_derivative(0:3)
    real*8, intent(in) :: connection(0:3,0:3,0:3), metric_derivative(0:3,0:3,1:3)
    integer, intent(inout) :: failure_count
    real*8 :: expected3(3,3), expected4(0:3,0:3), expected_connection(0:3,0:3,0:3)
    real*8 :: expected_derivative(0:3,0:3,1:3), sin_angle, cos_angle

    sin_angle = sin(angle)
    cos_angle = cos(angle)
    expected3 = 0.0d0
    expected3(1,1) = radius**2
    expected3(2,2) = radius**2
    expected3(3,3) = radius**2*sin_angle**2
    expected4 = 0.0d0
    expected4(0,0) = -1.0d0
    expected4(1,1) = 1.0d0/radius**2
    expected4(2,2) = 1.0d0/radius**2
    expected4(3,3) = 1.0d0/(radius**2*sin_angle**2)
    expected_connection = 0.0d0
    expected_connection(1,1,1) = 1.0d0
    expected_connection(1,2,2) = -1.0d0
    expected_connection(1,3,3) = -sin_angle**2
    expected_connection(2,1,2) = 1.0d0
    expected_connection(2,2,1) = 1.0d0
    expected_connection(2,3,3) = -sin_angle*cos_angle
    expected_connection(3,1,3) = 1.0d0
    expected_connection(3,3,1) = 1.0d0
    expected_connection(3,2,3) = cos_angle/sin_angle
    expected_connection(3,3,2) = cos_angle/sin_angle
    expected_derivative = 0.0d0
    expected_derivative(1,1,1) = 2.0d0*radius**2
    expected_derivative(2,2,1) = 2.0d0*radius**2
    expected_derivative(3,3,1) = 2.0d0*radius**2*sin_angle**2
    expected_derivative(3,3,2) = 2.0d0*radius**2*sin_angle*cos_angle

    call check_close('log lapse', lapse, 1.0d0, tolerance, failure_count)
    call check_vector('log shift', shift, [0.0d0,0.0d0,0.0d0], tolerance, failure_count)
    call check_matrix3('log spatial metric', spatial_metric, expected3, tolerance, failure_count)
    call check_matrix4('log inverse metric', inverse_metric, expected4, tolerance, failure_count)
    call check_close('log determinant', det_gamma, radius**6*sin_angle**2, tolerance, failure_count)
    call check_vector4('log lapse derivative', lapse_derivative, [0.0d0,0.0d0,0.0d0,0.0d0], &
                       tolerance, failure_count)
    call check_tensor3('log Christoffel', connection, expected_connection, tolerance, failure_count)
    call check_derivatives('log metric derivatives', metric_derivative, expected_derivative, tolerance, failure_count)
  end subroutine check_log_form


  subroutine check_flat_riemann(coordinate1, coordinate2, logarithmic, allowed_error, failure_count)
    real*8, intent(in) :: coordinate1, coordinate2, allowed_error
    logical, intent(in) :: logarithmic
    integer, intent(inout) :: failure_count
    integer :: rho_index, sigma_index, mu_index, nu_index, lambda_index
    real*8 :: epsilon_fd, max_curvature, value
    real*8 :: connection0(0:3,0:3,0:3), connection_plus(0:3,0:3,0:3)
    real*8 :: connection_minus(0:3,0:3,0:3), derivative_connection(0:3,0:3,0:3,0:3)

    epsilon_fd = 1.0d-5
    derivative_connection = 0.0d0
    call calculate_christoffel_symbols(coordinate1, coordinate2, connection0)
    call calculate_christoffel_symbols(coordinate1 + epsilon_fd, coordinate2, connection_plus)
    call calculate_christoffel_symbols(coordinate1 - epsilon_fd, coordinate2, connection_minus)
    derivative_connection(:,:,:,1) = (connection_plus - connection_minus)/(2.0d0*epsilon_fd)
    call calculate_christoffel_symbols(coordinate1, coordinate2 + epsilon_fd, connection_plus)
    call calculate_christoffel_symbols(coordinate1, coordinate2 - epsilon_fd, connection_minus)
    derivative_connection(:,:,:,2) = (connection_plus - connection_minus)/(2.0d0*epsilon_fd)

    max_curvature = 0.0d0
    do rho_index = 0, 3
      do sigma_index = 0, 3
        do mu_index = 0, 3
          do nu_index = 0, 3
            value = derivative_connection(rho_index,nu_index,sigma_index,mu_index) - &
                    derivative_connection(rho_index,mu_index,sigma_index,nu_index)
            do lambda_index = 0, 3
              value = value + connection0(rho_index,mu_index,lambda_index)* &
                              connection0(lambda_index,nu_index,sigma_index) - &
                              connection0(rho_index,nu_index,lambda_index)* &
                              connection0(lambda_index,mu_index,sigma_index)
            end do
            max_curvature = max(max_curvature, abs(value))
          end do
        end do
      end do
    end do

    if (max_curvature > allowed_error) then
      if (logarithmic) then
        write(*,*) 'FAIL: log-coordinate Riemann tensor max=', max_curvature
      else
        write(*,*) 'FAIL: physical-coordinate Riemann tensor max=', max_curvature
      end if
      failure_count = failure_count + 1
    end if
  end subroutine check_flat_riemann


  subroutine check_close(label, value, expected, allowed_error, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: value, expected, allowed_error
    integer, intent(inout) :: failure_count
    if (abs(value - expected) > allowed_error*max(1.0d0,abs(expected))) then
      write(*,*) 'FAIL: ', trim(label), ' value=', value, ' expected=', expected
      failure_count = failure_count + 1
    end if
  end subroutine check_close


  subroutine check_vector(label, value, expected, allowed_error, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: value(3), expected(3), allowed_error
    integer, intent(inout) :: failure_count
    if (maxval(abs(value - expected)) > allowed_error*max(1.0d0,maxval(abs(expected)))) then
      write(*,*) 'FAIL: ', trim(label)
      failure_count = failure_count + 1
    end if
  end subroutine check_vector


  subroutine check_vector4(label, value, expected, allowed_error, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: value(0:3), expected(0:3), allowed_error
    integer, intent(inout) :: failure_count
    if (maxval(abs(value - expected)) > allowed_error*max(1.0d0,maxval(abs(expected)))) then
      write(*,*) 'FAIL: ', trim(label)
      failure_count = failure_count + 1
    end if
  end subroutine check_vector4


  subroutine check_matrix3(label, value, expected, allowed_error, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: value(3,3), expected(3,3), allowed_error
    integer, intent(inout) :: failure_count
    if (maxval(abs(value - expected)) > allowed_error*max(1.0d0,maxval(abs(expected)))) then
      write(*,*) 'FAIL: ', trim(label)
      failure_count = failure_count + 1
    end if
  end subroutine check_matrix3


  subroutine check_matrix4(label, value, expected, allowed_error, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: value(0:3,0:3), expected(0:3,0:3), allowed_error
    integer, intent(inout) :: failure_count
    if (maxval(abs(value - expected)) > allowed_error*max(1.0d0,maxval(abs(expected)))) then
      write(*,*) 'FAIL: ', trim(label)
      failure_count = failure_count + 1
    end if
  end subroutine check_matrix4


  subroutine check_tensor3(label, value, expected, allowed_error, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: value(0:3,0:3,0:3), expected(0:3,0:3,0:3), allowed_error
    integer, intent(inout) :: failure_count
    if (maxval(abs(value - expected)) > allowed_error*max(1.0d0,maxval(abs(expected)))) then
      write(*,*) 'FAIL: ', trim(label)
      failure_count = failure_count + 1
    end if
  end subroutine check_tensor3


  subroutine check_derivatives(label, value, expected, allowed_error, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: value(0:3,0:3,1:3), expected(0:3,0:3,1:3), allowed_error
    integer, intent(inout) :: failure_count
    if (maxval(abs(value - expected)) > allowed_error*max(1.0d0,maxval(abs(expected)))) then
      write(*,*) 'FAIL: ', trim(label)
      failure_count = failure_count + 1
    end if
  end subroutine check_derivatives

end program test_minkowski_spherical
