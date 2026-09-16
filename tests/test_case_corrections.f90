program test_case_corrections
  use variables
  use metrics
  use conditions
  use initialization, only: initialize_problem
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer :: failures

  failures = 0
  call check_axisymmetric_jet_preset(failures)
  call check_cylindrical_geometry(failures)
  call check_khi_perturbation(failures)
  call check_offaxis_logarithmic_radius(failures)
  call check_dust_solution_and_boundaries(failures)
  call check_shift_aware_outflow(failures)

  if (failures /= 0) then
    write(*,'(A,I0)') 'CASE CORRECTION TESTS FAILED: ', failures
    error stop 1
  end if
  write(*,'(A)') 'CASE CORRECTION TESTS PASSED'

contains

  subroutine check_axisymmetric_jet_preset(failure_count)
    integer, intent(inout) :: failure_count
    integer :: parameter_unit
    real*8, parameter :: tolerance = 2.0d-13
    character(len=*), parameter :: parameter_file = '/tmp/grhd_jet_preset_test.par'

    open(newunit=parameter_unit, file=parameter_file, status='replace', action='write')
    write(parameter_unit,'(A)') 'problem = axisymmetric_jet'
    close(parameter_unit)

    case_id = 5
    call initialize_problem(parameter_file, .true.)
    call assert_close('Jet azimuth lower face', y_min, 0.0d0, tolerance, failure_count)
    call assert_close('Jet azimuth upper face', y_max, pi, tolerance, failure_count)
    call assert_close('Jet meridional plane at phi=pi/2', 0.5d0*(y_min+y_max), &
      0.5d0*pi, tolerance, failure_count)

    open(newunit=parameter_unit, file=parameter_file, status='old', action='read')
    close(parameter_unit, status='delete')
  end subroutine check_axisymmetric_jet_preset

  subroutine check_cylindrical_geometry(failure_count)
    integer, intent(inout) :: failure_count
    real*8, parameter :: rho = 2.5d0, tolerance = 2.0d-13
    real*8 :: alpha, determinant, beta(3), gamma(3,3), gmunu(0:3,0:3)
    real*8 :: dlnalpha(0:3), connection(0:3,0:3,0:3)
    real*8 :: derivatives(0:3,0:3,1:3), expected_gamma(3,3)
    real*8 :: expected_inverse(0:3,0:3), expected_connection(0:3,0:3,0:3)
    real*8 :: expected_derivatives(0:3,0:3,1:3)

    metric_type = 'Minkowski'
    geom_type = 'Cylindrical'
    use_log_r = .false.
    call set_metric_type()
    call calculate_metric(rho, 0.0d0, alpha, beta, gamma, gmunu, determinant, dlnalpha)
    call calculate_christoffel_symbols(rho, 0.0d0, connection)
    call calculate_metric_derivatives(rho, 0.0d0, derivatives)

    expected_gamma = 0.0d0
    expected_gamma(1,1) = 1.0d0
    expected_gamma(2,2) = rho**2
    expected_gamma(3,3) = 1.0d0
    expected_inverse = 0.0d0
    expected_inverse(0,0) = -1.0d0
    expected_inverse(1,1) = 1.0d0
    expected_inverse(2,2) = 1.0d0/rho**2
    expected_inverse(3,3) = 1.0d0
    expected_connection = 0.0d0
    expected_connection(1,2,2) = -rho
    expected_connection(2,1,2) = 1.0d0/rho
    expected_connection(2,2,1) = 1.0d0/rho
    expected_derivatives = 0.0d0
    expected_derivatives(2,2,1) = 2.0d0*rho

    call assert_close('cylindrical lapse', alpha, 1.0d0, tolerance, failure_count)
    call assert_close('cylindrical determinant', determinant, rho**2, tolerance, failure_count)
    call assert_array('cylindrical shift', beta, [0.0d0,0.0d0,0.0d0], tolerance, failure_count)
    call assert_array('cylindrical metric order', gamma, expected_gamma, tolerance, failure_count)
    call assert_array('cylindrical inverse metric', gmunu, expected_inverse, tolerance, failure_count)
    call assert_array('cylindrical lapse derivative', dlnalpha, &
      [0.0d0,0.0d0,0.0d0,0.0d0], tolerance, failure_count)
    call assert_array('cylindrical connection', connection, expected_connection, tolerance, failure_count)
    call assert_array('cylindrical metric derivatives', derivatives, expected_derivatives, &
      tolerance, failure_count)
  end subroutine check_cylindrical_geometry


  subroutine check_khi_perturbation(failure_count)
    integer, intent(inout) :: failure_count
    real*8, parameter :: tolerance = 2.0d-13
    real*8 :: expected_transverse_velocity
    integer :: i

    nghost = 1
    nx = 2
    ny = 1
    nz = 3
    allocate(x(-nghost:nx+nghost), y(-nghost:ny+nghost), z(-nghost:nz+nghost))
    allocate(p(neq,-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost))
    x = 0.0d0
    y = 0.0d0
    z = 0.0d0
    x(1) = 0.05d0
    x(2) = 0.10d0
    z(1) = 0.20d0
    z(2) = 0.25d0
    z(3) = 0.40d0

    case_name = 'KHI'
    khi_half_width = 0.25d0
    khi_rho_inner = 2.0d0
    khi_rho_outer = 1.0d0
    khi_vx_inner = 0.5d0
    khi_vx_outer = -0.5d0
    khi_pressure = 2.5d0
    khi_perturbation_amplitude = 0.01d0
    khi_perturbation_wave_number = 10.0d0
    khi_perturbation_width = 0.05d0
    call set_initial_conditions()

    expected_transverse_velocity = khi_perturbation_amplitude*exp(-0.5d0)
    call assert_close('KHI inner base velocity unchanged', p(eq_vx,1,1,1), &
      khi_vx_inner, tolerance, failure_count)
    call assert_close('KHI outer base velocity unchanged', p(eq_vx,1,1,2), &
      khi_vx_outer, tolerance, failure_count)
    call assert_close('KHI perturbation is transverse', p(eq_vz,1,1,1), &
      expected_transverse_velocity, tolerance, failure_count)
    call assert_close('KHI perturbation peaks at interface', p(eq_vz,1,1,2), &
      khi_perturbation_amplitude, tolerance, failure_count)
    call assert_close('KHI longitudinal node remains exact', p(eq_vz,2,1,2), &
      0.0d0, tolerance, failure_count)
    do i = 1, nx
      call assert_close('KHI unused y velocity', p(eq_vy,i,1,1), 0.0d0, &
        tolerance, failure_count)
    end do

    deallocate(p, x, y, z)
  end subroutine check_khi_perturbation


  subroutine check_offaxis_logarithmic_radius(failure_count)
    integer, intent(inout) :: failure_count
    real*8, parameter :: tolerance = 2.0d-13
    real*8 :: expected_profile

    nghost = 1
    nx = 2
    ny = 1
    nz = 1
    allocate(x(-nghost:nx+nghost), y(-nghost:ny+nghost), z(-nghost:nz+nghost))
    allocate(p(neq,-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost))
    offaxis_radial_center = 15.0d0
    offaxis_phi_center = pi
    offaxis_width = 1.5d0
    offaxis_background_density = 0.1d0
    offaxis_background_pressure = 0.2d0
    offaxis_density_amplitude = 1.0d0
    offaxis_pressure_amplitude = 0.5d0
    use_log_r = .true.
    case_name = 'OffAxis'
    x = log(offaxis_radial_center)
    x(2) = log(offaxis_radial_center + offaxis_width)
    y = 0.5d0*pi
    z = offaxis_phi_center
    call set_initial_conditions()

    call assert_close('OffAxis log peak density', p(eq_de,1,1,1), &
      offaxis_background_density + offaxis_density_amplitude, tolerance, failure_count)
    call assert_close('OffAxis log peak pressure', p(eq_pr,1,1,1), &
      offaxis_background_pressure + offaxis_pressure_amplitude, tolerance, failure_count)
    expected_profile = exp(-1.0d0)
    call assert_close('OffAxis log physical radial width', p(eq_de,2,1,1), &
      offaxis_background_density + offaxis_density_amplitude*expected_profile, &
      tolerance, failure_count)

    deallocate(p, x, y, z)
  end subroutine check_offaxis_logarithmic_radius


  subroutine check_dust_solution_and_boundaries(failure_count)
    integer, intent(inout) :: failure_count
    real*8, parameter :: radius = 5.0d0, tolerance = 2.0d-12
    real*8 :: state_physical(neq), state_logarithmic(neq), expected_outer(neq)
    real*8 :: alpha, beta(3), gamma(3,3), lorentz_factor, u_radial, mass_flux
    real*8 :: u_time, u_covariant_time
    real*8, allocatable :: q(:,:,:,:)
    integer :: i

    bh_mass = 1.0d0
    a_spin = 0.0d0
    dust_accretion_constant = -0.5d0
    p_floor = 1.0d-12
    metric_type = 'Eddington-Finkelstein'
    geom_type = 'Spherical'
    ny = 1

    use_log_r = .false.
    call set_metric_type()
    call evaluate_dust_accretion_state(radius, state_physical)
    call calculate_metric(radius, 0.5d0*pi, alpha=alpha, beta=beta, gamma=gamma)
    lorentz_factor = 1.0d0/sqrt(1.0d0-gamma(1,1)*state_physical(eq_vx)**2)
    u_radial = lorentz_factor*(state_physical(eq_vx)-beta(1)/alpha)
    mass_flux = radius**2*state_physical(eq_de)*u_radial
    u_time = lorentz_factor/alpha
    u_covariant_time = -(1.0d0-2.0d0*bh_mass/radius)*u_time + &
      (2.0d0*bh_mass/radius)*u_radial
    call assert_close('Dust physical mass flux', mass_flux, dust_accretion_constant, &
      tolerance, failure_count)
    call assert_close('Dust free fall energy', u_covariant_time, -1.0d0, &
      tolerance, failure_count)
    call assert_close('Dust cold pressure floor', state_physical(eq_pr), p_floor, &
      tolerance, failure_count)

    use_log_r = .true.
    call set_metric_type()
    call evaluate_dust_accretion_state(log(radius), state_logarithmic)
    call calculate_metric(log(radius), 0.5d0*pi, alpha=alpha, beta=beta, gamma=gamma)
    lorentz_factor = 1.0d0/sqrt(1.0d0-gamma(1,1)*state_logarithmic(eq_vx)**2)
    u_radial = radius*lorentz_factor*(state_logarithmic(eq_vx)-beta(1)/alpha)
    mass_flux = radius**2*state_logarithmic(eq_de)*u_radial
    call assert_close('Dust log mass flux', mass_flux, dust_accretion_constant, &
      tolerance, failure_count)
    call assert_close('Dust log density', state_logarithmic(eq_de), &
      state_physical(eq_de), tolerance, failure_count)
    call assert_close('Dust log vector transformation', radius*state_logarithmic(eq_vx), &
      state_physical(eq_vx), tolerance, failure_count)

    use_log_r = .false.
    call set_metric_type()
    nghost = 2
    nx = 3
    ny = 2
    nz = 2
    allocate(x(-nghost:nx+nghost))
    allocate(q(neq,-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost))
    do i = -nghost, nx+nghost
      x(i) = 2.5d0 + dble(i)
    end do
    q = -huge(1.0d0)
    q(:,1:nx,1:ny,1:nz) = 0.0d0
    q(eq_de,1:nx,1:ny,1:nz) = 0.2d0
    q(eq_pr,1:nx,1:ny,1:nz) = 0.3d0
    q(eq_vx,1:nx,1:ny,1:nz) = 0.1d0
    case_name = 'Dust'
    call set_boundary_conditions(q)
    call assert_close('Dust inner outflow clips outward velocity', q(eq_vx,0,1,1), &
      0.0d0, tolerance, failure_count)
    call assert_close('Dust inner density copy', q(eq_de,0,1,1), &
      q(eq_de,1,1,1), tolerance, failure_count)
    call evaluate_dust_accretion_state(x(nx+1), expected_outer)
    call assert_array('Dust analytic outer injection', q(:,nx+1,1,1), &
      expected_outer, tolerance, failure_count)
    call assert_array('Dust transverse periodic boundary', q(:,2,0,1), &
      q(:,2,ny,1), tolerance, failure_count)

    deallocate(q, x)
  end subroutine check_dust_solution_and_boundaries


  subroutine check_shift_aware_outflow(failure_count)
    integer, intent(inout) :: failure_count
    real*8, parameter :: tolerance = 2.0d-13, expected_zero_flux_velocity = 0.25d0
    real*8, allocatable :: q(:,:,:,:)

    nghost = 1
    nx = 2
    ny = 1
    nz = 1
    allocate(q(neq,-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost))
    allocate(alpha_c(nx,ny,nz), beta_c(3,nx,ny,nz))
    q = 0.0d0
    q(eq_de,:,:,:) = 1.0d0
    q(eq_pr,:,:,:) = 1.0d0
    q(eq_vx,1,1,1) = 0.40d0
    q(eq_vx,2,1,1) = 0.00d0
    alpha_c = 0.80d0
    beta_c = 0.0d0
    beta_c(1,:,:,:) = 0.20d0
    case_name = 'OffAxis'

    call set_boundary_conditions(q)
    call assert_close('shift-aware inner zero-flux ceiling', q(eq_vx,0,1,1), &
      expected_zero_flux_velocity, tolerance, failure_count)
    call assert_close('shift-aware outer zero-flux floor', q(eq_vx,nx+1,1,1), &
      expected_zero_flux_velocity, tolerance, failure_count)

    deallocate(q, alpha_c, beta_c)
  end subroutine check_shift_aware_outflow


  subroutine assert_close(label, actual, expected, tolerance, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: actual, expected, tolerance
    integer, intent(inout) :: failure_count

    if (.not. ieee_is_finite(actual) .or. &
        abs(actual-expected) > tolerance*max(1.0d0,abs(expected))) then
      failure_count = failure_count + 1
      write(*,'(A,1X,A,2(1X,ES16.8))') 'FAIL', trim(label), actual, expected
    end if
  end subroutine assert_close


  subroutine assert_array(label, actual, expected, tolerance, failure_count)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: actual(..), expected(..), tolerance
    integer, intent(inout) :: failure_count

    select rank(actual)
    rank(1)
      select rank(expected)
      rank(1)
        if (any(.not. ieee_is_finite(actual)) .or. &
            maxval(abs(actual-expected)) > tolerance*max(1.0d0,maxval(abs(expected)))) then
          failure_count = failure_count + 1
          write(*,'(A,1X,A)') 'FAIL', trim(label)
        end if
      end select
    rank(2)
      select rank(expected)
      rank(2)
        if (any(.not. ieee_is_finite(actual)) .or. &
            maxval(abs(actual-expected)) > tolerance*max(1.0d0,maxval(abs(expected)))) then
          failure_count = failure_count + 1
          write(*,'(A,1X,A)') 'FAIL', trim(label)
        end if
      end select
    rank(3)
      select rank(expected)
      rank(3)
        if (any(.not. ieee_is_finite(actual)) .or. &
            maxval(abs(actual-expected)) > tolerance*max(1.0d0,maxval(abs(expected)))) then
          failure_count = failure_count + 1
          write(*,'(A,1X,A)') 'FAIL', trim(label)
        end if
      end select
    end select
  end subroutine assert_array

end program test_case_corrections
