! =======================================================================================
! Modulo: parameters
! ---------------------------------------------------------------------------------------
! Lector de archivos de parametros de GRHD2. El formato es deliberadamente sencillo:
!
!   nombre = valor       # comentario opcional
!
! Los nombres y los valores simbolicos no distinguen mayusculas de minusculas. Las rutas
! conservan su escritura original. Los limites theta_* y phi_* se expresan en unidades de
! pi en el archivo y se convierten a radianes al leerlos.
! =======================================================================================
module parameters
  use variables
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer, parameter :: MAX_PARAMETER_LINE = 1024
  integer, parameter :: MAX_PARAMETER_KEY = 64
  integer, parameter :: MAX_PARAMETER_KEYS = 128

  public :: read_parameter_header, read_parameter_file, validate_parameters, &
            print_parameter_summary

contains

  subroutine read_parameter_header(file_name)
    character(len=*), intent(in) :: file_name
    character(len=MAX_PARAMETER_LINE) :: line, key, value
    character(len=MAX_PARAMETER_KEY) :: seen_keys(MAX_PARAMETER_KEYS)
    integer :: unit_file, io_status, line_number, equals_position, comment_position
    integer :: seen_count
    logical :: problem_found, restart_file_found

    do_restart = .false.
    problem_found = .false.
    restart_file_found = .false.
    seen_count = 0

    open(newunit=unit_file, file=trim(file_name), status='old', action='read', &
         iostat=io_status)
    if (io_status /= 0) then
      write(*,'(A,1X,A)') 'ERROR: no se pudo abrir el archivo de parametros:', &
                           trim(file_name)
      error stop 1
    end if

    line_number = 0
    do
      read(unit_file, '(A)', iostat=io_status) line
      if (io_status < 0) exit
      if (io_status > 0) call parameter_error(file_name, line_number + 1, &
                                               'no se pudo leer la linea')
      line_number = line_number + 1
      comment_position = index(line, '#')
      if (comment_position > 0) line(comment_position:) = ' '
      if (len_trim(line) == 0) cycle

      equals_position = index(line, '=')
      if (equals_position <= 1 .or. equals_position == len_trim(line)) then
        call parameter_error(file_name, line_number, 'se esperaba nombre = valor')
      end if
      if (index(line(equals_position + 1:), '=') > 0) then
        call parameter_error(file_name, line_number, 'hay mas de un signo =')
      end if

      key = lowercase(trim(adjustl(line(:equals_position - 1))))
      value = trim(adjustl(line(equals_position + 1:)))
      call register_key(key, seen_keys, seen_count, file_name, line_number)

      select case(trim(key))
      case('problem')
        call parse_problem(value, case_id, file_name, line_number)
        problem_found = .true.
      case('restart')
        call parse_logical_value(value, do_restart, file_name, line_number)
      case('restart_file')
        call parse_string_value(value, restart_file, file_name, line_number)
        restart_file_found = .true.
      end select
    end do
    close(unit_file)

    if (.not. do_restart .and. .not. problem_found) then
      call parameter_error(file_name, 0, &
                           'falta el parametro obligatorio problem para un arranque en frio')
    end if
    if (do_restart .and. .not. restart_file_found) then
      call parameter_error(file_name, 0, &
                           'restart = true requiere el parametro restart_file')
    end if
  end subroutine read_parameter_header


  subroutine read_parameter_file(file_name, restart_mode)
    character(len=*), intent(in) :: file_name
    logical, intent(in) :: restart_mode
    character(len=MAX_PARAMETER_LINE) :: line, key, value
    character(len=MAX_PARAMETER_KEY) :: seen_keys(MAX_PARAMETER_KEYS)
    integer :: unit_file, io_status, line_number, equals_position, comment_position
    integer :: seen_count, parsed_case
    real*8 :: pi_coefficient
    logical :: parsed_restart

    seen_count = 0
    open(newunit=unit_file, file=trim(file_name), status='old', action='read', &
         iostat=io_status)
    if (io_status /= 0) then
      write(*,'(A,1X,A)') 'ERROR: no se pudo abrir el archivo de parametros:', &
                           trim(file_name)
      error stop 1
    end if

    line_number = 0
    do
      read(unit_file, '(A)', iostat=io_status) line
      if (io_status < 0) exit
      if (io_status > 0) call parameter_error(file_name, line_number + 1, &
                                               'no se pudo leer la linea')
      line_number = line_number + 1
      comment_position = index(line, '#')
      if (comment_position > 0) line(comment_position:) = ' '
      if (len_trim(line) == 0) cycle

      equals_position = index(line, '=')
      if (equals_position <= 1 .or. equals_position == len_trim(line)) then
        call parameter_error(file_name, line_number, 'se esperaba nombre = valor')
      end if
      if (index(line(equals_position + 1:), '=') > 0) then
        call parameter_error(file_name, line_number, 'hay mas de un signo =')
      end if

      key = lowercase(trim(adjustl(line(:equals_position - 1))))
      value = trim(adjustl(line(equals_position + 1:)))
      call register_key(key, seen_keys, seen_count, file_name, line_number)
      if (restart_mode .and. .not. is_restart_control(key)) then
        call parameter_error(file_name, line_number, &
          'este parametro pertenece al checkpoint y no puede cambiarse al reiniciar')
      end if
      if (.not. restart_mode .and. is_initial_condition_key(key) .and. &
          .not. initial_parameter_applies(key, case_id)) then
        call parameter_error(file_name, line_number, &
          'el parametro no corresponde al problema seleccionado: '//trim(key))
      end if

      select case(trim(key))
      case('problem')
        call parse_problem(value, parsed_case, file_name, line_number)
        case_id = parsed_case
      case('restart')
        call parse_logical_value(value, parsed_restart, file_name, line_number)
        if (parsed_restart .neqv. restart_mode) then
          call parameter_error(file_name, line_number, &
                               'restart no coincide con el modo de arranque')
        end if
        do_restart = parsed_restart
      case('restart_file')
        call parse_string_value(value, restart_file, file_name, line_number)
      case('metric')
        call parse_metric(value, metric_type, file_name, line_number)
      case('geometry')
        call parse_geometry(value, geom_type, file_name, line_number)
      case('bh_mass')
        call parse_real_value(value, bh_mass, file_name, line_number)
      case('spin')
        call parse_real_value(value, a_spin, file_name, line_number)
      case('adiabatic_index')
        call parse_real_value(value, adb_idx, file_name, line_number)
      case('nx')
        call parse_integer_value(value, nx, file_name, line_number)
      case('ny')
        call parse_integer_value(value, ny, file_name, line_number)
      case('nz')
        call parse_integer_value(value, nz, file_name, line_number)
      case('r_min')
        call parse_real_value(value, r_min, file_name, line_number)
      case('r_max')
        call parse_real_value(value, r_max, file_name, line_number)
      case('theta_min')
        call parse_pi_coefficient(value, pi_coefficient, file_name, line_number)
        y_min = pi_coefficient * pi
      case('theta_max')
        call parse_pi_coefficient(value, pi_coefficient, file_name, line_number)
        y_max = pi_coefficient * pi
      case('phi_min')
        call parse_pi_coefficient(value, pi_coefficient, file_name, line_number)
        z_min = pi_coefficient * pi
      case('phi_max')
        call parse_pi_coefficient(value, pi_coefficient, file_name, line_number)
        z_max = pi_coefficient * pi
      case('y_min')
        call parse_real_value(value, y_min, file_name, line_number)
      case('y_max')
        call parse_real_value(value, y_max, file_name, line_number)
      case('z_min')
        call parse_real_value(value, z_min, file_name, line_number)
      case('z_max')
        call parse_real_value(value, z_max, file_name, line_number)
      case('logarithmic_r')
        call parse_logical_value(value, use_log_r, file_name, line_number)
      case('sod_interface')
        call parse_real_value(value, sod_interface, file_name, line_number)
      case('sod_density_left')
        call parse_real_value(value, sod_rho_left, file_name, line_number)
      case('sod_density_right')
        call parse_real_value(value, sod_rho_right, file_name, line_number)
      case('sod_pressure_left')
        call parse_real_value(value, sod_pressure_left, file_name, line_number)
      case('sod_pressure_right')
        call parse_real_value(value, sod_pressure_right, file_name, line_number)
      case('strong_interface')
        call parse_real_value(value, strong_interface, file_name, line_number)
      case('strong_density_left')
        call parse_real_value(value, strong_rho_left, file_name, line_number)
      case('strong_density_right')
        call parse_real_value(value, strong_rho_right, file_name, line_number)
      case('strong_pressure_left')
        call parse_real_value(value, strong_pressure_left, file_name, line_number)
      case('strong_pressure_right')
        call parse_real_value(value, strong_pressure_right, file_name, line_number)
      case('shu_interface')
        call parse_real_value(value, shu_interface, file_name, line_number)
      case('shu_density_left')
        call parse_real_value(value, shu_rho_left, file_name, line_number)
      case('shu_pressure_left')
        call parse_real_value(value, shu_pressure_left, file_name, line_number)
      case('shu_density_right')
        call parse_real_value(value, shu_rho_right, file_name, line_number)
      case('shu_density_amplitude')
        call parse_real_value(value, shu_rho_amplitude, file_name, line_number)
      case('shu_wave_number')
        call parse_real_value(value, shu_wave_number, file_name, line_number)
      case('shu_pressure_right')
        call parse_real_value(value, shu_pressure_right, file_name, line_number)
      case('khi_half_width')
        call parse_real_value(value, khi_half_width, file_name, line_number)
      case('khi_density_inner')
        call parse_real_value(value, khi_rho_inner, file_name, line_number)
      case('khi_density_outer')
        call parse_real_value(value, khi_rho_outer, file_name, line_number)
      case('khi_velocity_inner')
        call parse_real_value(value, khi_vx_inner, file_name, line_number)
      case('khi_velocity_outer')
        call parse_real_value(value, khi_vx_outer, file_name, line_number)
      case('khi_pressure')
        call parse_real_value(value, khi_pressure, file_name, line_number)
      case('khi_perturbation_amplitude')
        call parse_real_value(value, khi_perturbation_amplitude, file_name, line_number)
      case('khi_perturbation_wave_number')
        call parse_real_value(value, khi_perturbation_wave_number, file_name, line_number)
      case('jet_ambient_density')
        call parse_real_value(value, jet_ambient_density, file_name, line_number)
      case('jet_ambient_pressure')
        call parse_real_value(value, jet_ambient_pressure, file_name, line_number)
      case('jet_ambient_velocity')
        call parse_real_value(value, jet_ambient_velocity, file_name, line_number)
      case('jet_nozzle_radius')
        call parse_real_value(value, jet_nozzle_radius, file_name, line_number)
      case('jet_nozzle_length')
        call parse_real_value(value, jet_nozzle_length, file_name, line_number)
      case('jet_density')
        call parse_real_value(value, jet_density, file_name, line_number)
      case('jet_pressure')
        call parse_real_value(value, jet_pressure, file_name, line_number)
      case('jet_velocity')
        call parse_real_value(value, jet_velocity, file_name, line_number)
      case('convergence_density')
        call parse_real_value(value, advected_wave_density, file_name, line_number)
      case('convergence_amplitude')
        call parse_real_value(value, advected_wave_amplitude, file_name, line_number)
      case('convergence_velocity')
        call parse_real_value(value, advected_wave_speed, file_name, line_number)
      case('convergence_pressure')
        call parse_real_value(value, advected_wave_pressure, file_name, line_number)
      case('convergence_wave_number')
        call parse_pi_coefficient(value, advected_wave_number, file_name, line_number)
      case('michel_critical_radius')
        call parse_real_value(value, michel_critical_radius, file_name, line_number)
      case('michel_critical_density')
        call parse_real_value(value, michel_critical_density, file_name, line_number)
      case('dust_accretion_constant')
        call parse_real_value(value, dust_accretion_constant, file_name, line_number)
      case('offaxis_radial_center')
        call parse_real_value(value, offaxis_radial_center, file_name, line_number)
      case('offaxis_phi_center')
        call parse_pi_coefficient(value, pi_coefficient, file_name, line_number)
        offaxis_phi_center = pi_coefficient * pi
      case('offaxis_width')
        call parse_real_value(value, offaxis_width, file_name, line_number)
      case('offaxis_background_density')
        call parse_real_value(value, offaxis_background_density, file_name, line_number)
      case('offaxis_background_pressure')
        call parse_real_value(value, offaxis_background_pressure, file_name, line_number)
      case('offaxis_density_amplitude')
        call parse_real_value(value, offaxis_density_amplitude, file_name, line_number)
      case('offaxis_pressure_amplitude')
        call parse_real_value(value, offaxis_pressure_amplitude, file_name, line_number)
      case('fm_inner_radius')
        call parse_real_value(value, fm_inner_radius, file_name, line_number)
      case('fm_pressure_max_radius')
        call parse_real_value(value, fm_pressure_max_radius, file_name, line_number)
      case('fm_polytropic_constant')
        call parse_real_value(value, fm_polytropic_constant, file_name, line_number)
      case('reconstruction')
        call parse_reconstruction(value, rec_method_id, file_name, line_number)
      case('riemann_solver')
        call parse_riemann_solver(value, riemann_solver_id, file_name, line_number)
      case('tvd_limiter')
        call parse_tvd_limiter(value, tvd_limiter_id, file_name, line_number)
      case('shock_sensor')
        call parse_logical_value(value, use_shock_sensor, file_name, line_number)
      case('final_time')
        call parse_real_value(value, final_time, file_name, line_number)
      case('cfl')
        call parse_real_value(value, CFL, file_name, line_number)
      case('save_interval')
        call parse_real_value(value, save_interval, file_name, line_number)
      case('checkpoint_interval')
        call parse_real_value(value, checkpoint_interval, file_name, line_number)
      case('output_prefix')
        call parse_string_value(value, output_prefix, file_name, line_number)
      case('output_folder')
        call parse_string_value(value, output_folder, file_name, line_number)
      case('vtk_mapping')
        call parse_vtk_mapping(value, vtk_mapping_id, file_name, line_number)
      case('extract_gw')
        call parse_logical_value(value, do_gw_extraction, file_name, line_number)
      case('extract_mdot')
        call parse_logical_value(value, do_mdot_extraction, file_name, line_number)
      case('ppi_diagnostics')
        call parse_logical_value(value, do_ppi_diagnostics, file_name, line_number)
      case('diagnostic_stride')
        call parse_integer_value(value, diagnostic_stride, file_name, line_number)
      case('perturbation')
        call parse_perturbation(value, perturbation_type, file_name, line_number)
        apply_perturbation = (perturbation_type /= PERT_NONE)
      case('perturbation_seed')
        call parse_integer_value(value, perturbation_seed, file_name, line_number)
      case('perturbation_time')
        call parse_real_value(value, perturbation_time, file_name, line_number)
      case('perturbation_amplitude')
        call parse_real_value(value, perturbation_amplitude, file_name, line_number)
      case('perturbation_mode')
        call parse_real_value(value, perturbation_mode, file_name, line_number)
      case default
        call parameter_error(file_name, line_number, 'parametro desconocido: '//trim(key))
      end select
    end do
    close(unit_file)
  end subroutine read_parameter_file


  subroutine validate_parameters(file_name)
    character(len=*), intent(in) :: file_name
    real*8 :: spin_tolerance

    if (.not. all([ieee_is_finite(r_min), ieee_is_finite(r_max), &
                   ieee_is_finite(y_min), ieee_is_finite(y_max), &
                   ieee_is_finite(z_min), ieee_is_finite(z_max), &
                   ieee_is_finite(adb_idx), ieee_is_finite(bh_mass), &
                   ieee_is_finite(a_spin), ieee_is_finite(final_time), &
                   ieee_is_finite(CFL), ieee_is_finite(save_interval), &
                   ieee_is_finite(checkpoint_interval), &
                   ieee_is_finite(perturbation_time), &
                   ieee_is_finite(perturbation_amplitude), &
                   ieee_is_finite(perturbation_mode)])) then
      call parameter_error(file_name, 0, 'los parametros reales deben ser finitos')
    end if

    if (nx < 1 .or. ny < 1 .or. nz < 1) then
      call parameter_error(file_name, 0, 'nx, ny y nz deben ser positivos')
    end if
    if (r_max <= r_min) then
      call parameter_error(file_name, 0, 'r_max debe ser mayor que r_min')
    end if
    if (use_log_r .and. r_min <= 0.0d0) then
      call parameter_error(file_name, 0, 'una malla radial logaritmica requiere r_min > 0')
    end if
    if (adb_idx <= 1.0d0) then
      call parameter_error(file_name, 0, 'adiabatic_index debe ser mayor que 1')
    end if
    if (bh_mass < 0.0d0) then
      call parameter_error(file_name, 0, 'bh_mass no puede ser negativa')
    end if
    if (final_time < 0.0d0 .or. CFL <= 0.0d0 .or. save_interval <= 0.0d0 .or. &
        checkpoint_interval <= 0.0d0) then
      call parameter_error(file_name, 0, &
        'se requiere final_time >= 0 y cfl/save_interval/checkpoint_interval > 0')
    end if
    if (len_trim(output_prefix) == 0 .or. len_trim(output_folder) == 0) then
      call parameter_error(file_name, 0, 'output_prefix y output_folder no pueden estar vacios')
    end if
    if (vtk_mapping_id /= VTK_MAP_PHYSICAL .and. &
        vtk_mapping_id /= VTK_MAP_UNTWISTED) then
      call parameter_error(file_name, 0, 'identificador de mapeo VTK invalido')
    end if
    if (diagnostic_stride < 1 .or. perturbation_seed < 0 .or. &
        perturbation_time < 0.0d0 .or. perturbation_amplitude < 0.0d0 .or. &
        perturbation_amplitude >= 1.0d0) then
      call parameter_error(file_name, 0, &
                           'controles de diagnostico o perturbacion invalidos')
    end if
    if (rec_method_id < REC_GODUNOV .or. rec_method_id > REC_WENO5) then
      call parameter_error(file_name, 0, 'identificador de reconstruccion invalido')
    end if
    if (riemann_solver_id /= RS_HLLE .and. riemann_solver_id /= RS_HLLC) then
      call parameter_error(file_name, 0, 'identificador de solucionador de Riemann invalido')
    end if
    if (rec_method_id == REC_TVD .and. &
        (tvd_limiter_id < LIM_MINMOD .or. tvd_limiter_id > LIM_MC)) then
      call parameter_error(file_name, 0, 'identificador de limitador TVD invalido')
    end if

    select case(trim(metric_type))
    case('Minkowski')
      if (abs(a_spin) > 0.0d0) then
        call parameter_error(file_name, 0, &
                             'Minkowski requiere spin = 0')
      end if
      if (trim(geom_type) /= 'Cartesian' .and. trim(geom_type) /= 'Cylindrical') then
        call parameter_error(file_name, 0, &
                             'Minkowski solo esta implementada en geometria cartesiana o cilindrica')
      end if
    case('Eddington-Finkelstein')
      if (abs(a_spin) > 0.0d0) then
        call parameter_error(file_name, 0, &
                             'Eddington-Finkelstein requiere spin = 0')
      end if
      if (bh_mass <= 0.0d0 .or. trim(geom_type) /= 'Spherical') then
        call parameter_error(file_name, 0, &
                             'Eddington-Finkelstein requiere bh_mass > 0 y geometria spherical')
      end if
    case('Kerr-Schild')
      if (bh_mass <= 0.0d0 .or. trim(geom_type) /= 'Spherical') then
        call parameter_error(file_name, 0, &
                             'Kerr-Schild requiere bh_mass > 0 y geometria spherical')
      end if
      spin_tolerance = 64.0d0 * epsilon(1.0d0) * max(1.0d0, bh_mass)
      if (abs(a_spin) > bh_mass + spin_tolerance) then
        call parameter_error(file_name, 0, 'Kerr-Schild requiere |spin| <= bh_mass')
      end if
    case default
      call parameter_error(file_name, 0, 'metrica interna invalida')
    end select

    if (trim(geom_type) == 'Spherical') then
      if (y_min < 0.0d0 .or. y_max > pi .or. y_max <= y_min) then
        call parameter_error(file_name, 0, &
                             'se requiere 0 <= theta_min < theta_max <= 1 (en unidades de pi)')
      end if
      if (z_max <= z_min .or. z_max - z_min > 2.0d0*pi + 64.0d0*epsilon(pi)) then
        call parameter_error(file_name, 0, &
                             'el intervalo azimutal debe tener longitud entre 0 y 2 pi')
      end if
    else
      if (y_max <= y_min .or. z_max <= z_min) then
        call parameter_error(file_name, 0, 'los limites de cada coordenada deben ser crecientes')
      end if
    end if
    call validate_initial_parameters(file_name)
  end subroutine validate_parameters


  subroutine validate_initial_parameters(file_name)
    character(len=*), intent(in) :: file_name
    real*8 :: angular_tolerance

    angular_tolerance = 128.0d0 * epsilon(pi)

    if (.not. initial_parameters_are_finite()) then
      call parameter_error(file_name, 0, &
                           'los parametros de condiciones iniciales deben ser finitos')
    end if

    select case(case_id)
    case(1)
      if (sod_interface <= r_min .or. sod_interface >= r_max .or. &
          min(sod_rho_left, sod_rho_right) <= 0.0d0 .or. &
          min(sod_pressure_left, sod_pressure_right) <= 0.0d0) then
        call parameter_error(file_name, 0, 'parametros invalidos para Sod')
      end if
    case(2)
      if (strong_interface <= r_min .or. strong_interface >= r_max .or. &
          min(strong_rho_left, strong_rho_right) <= 0.0d0 .or. &
          min(strong_pressure_left, strong_pressure_right) <= 0.0d0) then
        call parameter_error(file_name, 0, 'parametros invalidos para strong_shock')
      end if
    case(3)
      if (shu_interface <= r_min .or. shu_interface >= r_max .or. &
          shu_rho_left <= 0.0d0 .or. shu_rho_right <= abs(shu_rho_amplitude) .or. &
          min(shu_pressure_left, shu_pressure_right) <= 0.0d0 .or. &
          shu_wave_number <= 0.0d0) then
        call parameter_error(file_name, 0, 'parametros invalidos para Shu-Osher')
      end if
    case(4)
      if (khi_half_width <= 0.0d0 .or. min(khi_rho_inner, khi_rho_outer) <= 0.0d0 .or. &
          khi_pressure <= 0.0d0 .or. abs(khi_vx_inner) >= 1.0d0 .or. &
          abs(khi_vx_outer) >= 1.0d0 .or. khi_perturbation_amplitude < 0.0d0 .or. &
          khi_perturbation_amplitude >= 1.0d0 .or. &
          khi_perturbation_wave_number <= 0.0d0) then
        call parameter_error(file_name, 0, 'parametros invalidos para Kelvin-Helmholtz')
      end if
    case(5)
      if (min(jet_ambient_density, jet_density) <= 0.0d0 .or. &
          min(jet_ambient_pressure, jet_pressure) <= 0.0d0 .or. &
          jet_nozzle_radius <= 0.0d0 .or. jet_nozzle_length <= 0.0d0 .or. &
          abs(jet_ambient_velocity) >= 1.0d0 .or. abs(jet_velocity) >= 1.0d0) then
        call parameter_error(file_name, 0, 'parametros invalidos para el jet')
      end if
    case(6)
      if (advected_wave_density <= abs(advected_wave_amplitude) .or. &
          advected_wave_pressure <= 0.0d0 .or. abs(advected_wave_speed) >= 1.0d0 .or. &
          advected_wave_number <= 0.0d0) then
        call parameter_error(file_name, 0, 'parametros invalidos para convergencia')
      end if
    case(7)
      if (bh_mass <= 0.0d0 .or. michel_critical_radius <= 1.5d0*bh_mass .or. &
          michel_critical_density <= 0.0d0) then
        call parameter_error(file_name, 0, 'parametros invalidos para acrecion de Michel')
      end if
    case(8)
      if (bh_mass <= 0.0d0 .or. dust_accretion_constant >= 0.0d0) then
        call parameter_error(file_name, 0, 'parametros invalidos para acrecion de polvo')
      end if
    case(9)
      if (offaxis_radial_center <= 0.0d0 .or. offaxis_width <= 0.0d0 .or. &
          offaxis_background_density <= 0.0d0 .or. &
          offaxis_background_density + offaxis_density_amplitude <= 0.0d0 .or. &
          offaxis_background_pressure <= 0.0d0 .or. &
          offaxis_background_pressure + offaxis_pressure_amplitude <= 0.0d0) then
        call parameter_error(file_name, 0, 'parametros invalidos para offaxis_blast')
      end if
    case(10:12)
      if (bh_mass <= 0.0d0 .or. fm_inner_radius <= 0.0d0 .or. &
          fm_pressure_max_radius <= fm_inner_radius .or. &
          fm_polytropic_constant <= 0.0d0) then
        call parameter_error(file_name, 0, 'parametros invalidos para Fishbone-Moncrief')
      end if
      select case(case_id)
      case(10)
        if (ny /= 1 .or. abs(0.5d0*(y_min+y_max)-0.5d0*pi) > angular_tolerance) then
          call parameter_error(file_name, 0, &
            'fishbone_equatorial requiere ny=1 y la celda polar centrada en theta=pi/2')
        end if
      case(11)
        if (ny <= 1 .or. nz /= 1 .or. abs(y_min) > angular_tolerance .or. &
            abs(y_max-pi) > angular_tolerance) then
          call parameter_error(file_name, 0, &
            'fishbone_sagittal requiere ny>1, nz=1 y theta=[0,1] en unidades de pi')
        end if
      case(12)
        if (ny <= 1 .or. nz <= 1 .or. mod(nz,2) /= 0 .or. &
            abs(y_min) > angular_tolerance .or. abs(y_max-pi) > angular_tolerance .or. &
            abs((z_max-z_min)-2.0d0*pi) > angular_tolerance) then
          call parameter_error(file_name, 0, &
            'fishbone_3d requiere ny,nz>1, nz par, theta=[0,1] y longitud phi=2 [pi]')
        end if
      end select
    case default
      call parameter_error(file_name, 0, 'identificador de problema invalido')
    end select
  end subroutine validate_initial_parameters


  logical function initial_parameters_are_finite()
    initial_parameters_are_finite = all([ &
      ieee_is_finite(sod_interface), ieee_is_finite(sod_rho_left), &
      ieee_is_finite(sod_rho_right), ieee_is_finite(sod_pressure_left), &
      ieee_is_finite(sod_pressure_right), ieee_is_finite(strong_interface), &
      ieee_is_finite(strong_rho_left), ieee_is_finite(strong_rho_right), &
      ieee_is_finite(strong_pressure_left), ieee_is_finite(strong_pressure_right), &
      ieee_is_finite(shu_interface), ieee_is_finite(shu_rho_left), &
      ieee_is_finite(shu_pressure_left), ieee_is_finite(shu_rho_right), &
      ieee_is_finite(shu_rho_amplitude), ieee_is_finite(shu_wave_number), &
      ieee_is_finite(shu_pressure_right), ieee_is_finite(khi_half_width), &
      ieee_is_finite(khi_rho_inner), ieee_is_finite(khi_rho_outer), &
      ieee_is_finite(khi_vx_inner), ieee_is_finite(khi_vx_outer), &
      ieee_is_finite(khi_pressure), ieee_is_finite(khi_perturbation_amplitude), &
      ieee_is_finite(khi_perturbation_wave_number), &
      ieee_is_finite(jet_ambient_density), ieee_is_finite(jet_ambient_pressure), &
      ieee_is_finite(jet_ambient_velocity), ieee_is_finite(jet_nozzle_radius), &
      ieee_is_finite(jet_nozzle_length), ieee_is_finite(jet_density), &
      ieee_is_finite(jet_pressure), ieee_is_finite(jet_velocity), &
      ieee_is_finite(advected_wave_density), ieee_is_finite(advected_wave_amplitude), &
      ieee_is_finite(advected_wave_speed), ieee_is_finite(advected_wave_pressure), &
      ieee_is_finite(advected_wave_number), ieee_is_finite(michel_critical_radius), &
      ieee_is_finite(michel_critical_density), ieee_is_finite(dust_accretion_constant), &
      ieee_is_finite(offaxis_radial_center), ieee_is_finite(offaxis_phi_center), &
      ieee_is_finite(offaxis_width), ieee_is_finite(offaxis_background_density), &
      ieee_is_finite(offaxis_background_pressure), &
      ieee_is_finite(offaxis_density_amplitude), &
      ieee_is_finite(offaxis_pressure_amplitude), ieee_is_finite(fm_inner_radius), &
      ieee_is_finite(fm_pressure_max_radius), ieee_is_finite(fm_polytropic_constant)])
  end function initial_parameters_are_finite


  subroutine print_parameter_summary(file_name)
    character(len=*), intent(in) :: file_name

    write(*,'(/,A)') '=========================================='
    write(*,'(A,1X,A)') 'Parametros leidos de:', trim(file_name)
    write(*,'(A,I0,2A)') 'Problema: ', case_id, ' - ', trim(problem_label(case_id))
    write(*,'(A,L1)') 'Restart: ', do_restart
    if (do_restart) write(*,'(A,1X,A)') 'Checkpoint:', trim(restart_file)
    write(*,'(A,1X,A,2X,A,1X,A)') 'Metrica:', trim(metric_type), &
                                    'Geometria:', trim(geom_type)
    write(*,'(A,3(I0,1X))') 'Malla nx ny nz: ', nx, ny, nz
    write(*,'(A,2(ES13.5,1X))') 'r_min r_max: ', r_min, r_max
    if (trim(geom_type) == 'Spherical') then
      write(*,'(A,2(F10.6,1X))') 'theta_min theta_max [pi]: ', y_min/pi, y_max/pi
      write(*,'(A,2(F10.6,1X))') 'phi_min phi_max [pi]: ', z_min/pi, z_max/pi
    end if
    call print_initial_parameter_summary()
    write(*,'(A,1X,A)') 'Salida:', trim(output_folder)
    if (vtk_mapping_id == VTK_MAP_PHYSICAL) then
      write(*,'(A)') 'Mapeo VTK: physical'
    else
      write(*,'(A)') 'Mapeo VTK: untwisted'
    end if
    write(*,'(A)') '=========================================='
  end subroutine print_parameter_summary


  subroutine print_initial_parameter_summary()
    select case(case_id)
    case(1)
      write(*,'(A,5(ES12.4,1X))') 'IC Sod: x0 rhoL rhoR pL pR: ', sod_interface, &
        sod_rho_left, sod_rho_right, sod_pressure_left, sod_pressure_right
    case(2)
      write(*,'(A,5(ES12.4,1X))') 'IC Strong: x0 rhoL rhoR pL pR: ', strong_interface, &
        strong_rho_left, strong_rho_right, strong_pressure_left, strong_pressure_right
    case(3)
      write(*,'(A,7(ES12.4,1X))') 'IC Shu: x0 rhoL pL rhoR A k pR: ', shu_interface, &
        shu_rho_left, shu_pressure_left, shu_rho_right, shu_rho_amplitude, &
        shu_wave_number, shu_pressure_right
    case(4)
      write(*,'(A,8(ES12.4,1X))') 'IC KHI: width rhoi rhoo vxi vxo p A kpi: ', &
        khi_half_width, khi_rho_inner, khi_rho_outer, khi_vx_inner, khi_vx_outer, &
        khi_pressure, khi_perturbation_amplitude, khi_perturbation_wave_number
    case(5)
      write(*,'(A,8(ES12.4,1X))') 'IC Jet: rhoa pa va R L rhoj pj vj: ', &
        jet_ambient_density, jet_ambient_pressure, jet_ambient_velocity, &
        jet_nozzle_radius, jet_nozzle_length, jet_density, jet_pressure, jet_velocity
    case(6)
      write(*,'(A,5(ES12.4,1X))') 'IC Conv: rho A v p kpi: ', advected_wave_density, &
        advected_wave_amplitude, advected_wave_speed, advected_wave_pressure, &
        advected_wave_number
    case(7)
      write(*,'(A,2(ES12.4,1X))') 'IC Michel: rcrit rhocrit: ', &
        michel_critical_radius, michel_critical_density
    case(8)
      write(*,'(A,ES12.4)') 'IC Dust: constant: ', dust_accretion_constant
    case(9)
      write(*,'(A,7(ES12.4,1X))') 'IC Offaxis: r0 phi0[pi] width rhob pb Arho Ap: ', &
        offaxis_radial_center, offaxis_phi_center/pi, offaxis_width, &
        offaxis_background_density, offaxis_background_pressure, &
        offaxis_density_amplitude, offaxis_pressure_amplitude
    case(10:12)
      write(*,'(A,3(ES12.4,1X))') 'IC FM: r_in[M] r_pmax[M] K: ', &
        fm_inner_radius, fm_pressure_max_radius, fm_polytropic_constant
    end select
  end subroutine print_initial_parameter_summary


  subroutine register_key(key, seen_keys, seen_count, file_name, line_number)
    character(len=*), intent(in) :: key, file_name
    character(len=MAX_PARAMETER_KEY), intent(inout) :: seen_keys(:)
    integer, intent(inout) :: seen_count
    integer, intent(in) :: line_number
    integer :: i

    if (len_trim(key) == 0 .or. len_trim(key) > MAX_PARAMETER_KEY) then
      call parameter_error(file_name, line_number, 'nombre de parametro vacio o demasiado largo')
    end if
    do i = 1, seen_count
      if (trim(key) == trim(seen_keys(i))) then
        call parameter_error(file_name, line_number, 'parametro duplicado: '//trim(key))
      end if
    end do
    if (seen_count >= size(seen_keys)) then
      call parameter_error(file_name, line_number, 'demasiados parametros en el archivo')
    end if
    seen_count = seen_count + 1
    seen_keys(seen_count) = trim(key)
  end subroutine register_key


  logical function is_restart_control(key)
    character(len=*), intent(in) :: key

    select case(trim(key))
    case('restart', 'restart_file', 'final_time', 'cfl', 'save_interval', &
         'checkpoint_interval', 'output_prefix', 'output_folder', 'vtk_mapping', &
         'shock_sensor', &
         'extract_gw', 'extract_mdot', 'ppi_diagnostics', 'diagnostic_stride', &
         'perturbation', 'perturbation_seed', 'perturbation_time', &
         'perturbation_amplitude', 'perturbation_mode')
      is_restart_control = .true.
    case default
      is_restart_control = .false.
    end select
  end function is_restart_control


  logical function is_initial_condition_key(key)
    character(len=*), intent(in) :: key
    character(len=MAX_PARAMETER_KEY) :: normalized

    normalized = trim(key)
    is_initial_condition_key = &
      index(normalized, 'sod_') == 1 .or. &
      index(normalized, 'strong_') == 1 .or. &
      index(normalized, 'shu_') == 1 .or. &
      index(normalized, 'khi_') == 1 .or. &
      index(normalized, 'jet_') == 1 .or. &
      index(normalized, 'convergence_') == 1 .or. &
      index(normalized, 'michel_') == 1 .or. &
      index(normalized, 'dust_') == 1 .or. &
      index(normalized, 'offaxis_') == 1 .or. &
      index(normalized, 'fm_') == 1
  end function is_initial_condition_key


  logical function initial_parameter_applies(key, problem_id)
    character(len=*), intent(in) :: key
    integer, intent(in) :: problem_id

    select case(problem_id)
    case(1)
      initial_parameter_applies = index(trim(key), 'sod_') == 1
    case(2)
      initial_parameter_applies = index(trim(key), 'strong_') == 1
    case(3)
      initial_parameter_applies = index(trim(key), 'shu_') == 1
    case(4)
      initial_parameter_applies = index(trim(key), 'khi_') == 1
    case(5)
      initial_parameter_applies = index(trim(key), 'jet_') == 1
    case(6)
      initial_parameter_applies = index(trim(key), 'convergence_') == 1
    case(7)
      initial_parameter_applies = index(trim(key), 'michel_') == 1
    case(8)
      initial_parameter_applies = index(trim(key), 'dust_') == 1
    case(9)
      initial_parameter_applies = index(trim(key), 'offaxis_') == 1
    case(10:12)
      initial_parameter_applies = index(trim(key), 'fm_') == 1
    case default
      initial_parameter_applies = .false.
    end select
  end function initial_parameter_applies


  subroutine parse_integer_value(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    integer, intent(out) :: value
    integer, intent(in) :: line_number
    integer :: io_status

    read(text, *, iostat=io_status) value
    if (io_status /= 0) call parameter_error(file_name, line_number, &
                                              'se esperaba un entero: '//trim(text))
  end subroutine parse_integer_value


  subroutine parse_real_value(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    real*8, intent(out) :: value
    integer, intent(in) :: line_number
    integer :: io_status

    read(text, *, iostat=io_status) value
    if (io_status /= 0) call parameter_error(file_name, line_number, &
                                              'se esperaba un numero real: '//trim(text))
  end subroutine parse_real_value


  subroutine parse_pi_coefficient(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    real*8, intent(out) :: value
    integer, intent(in) :: line_number
    integer :: slash_position, io_status
    real*8 :: numerator, denominator

    slash_position = index(trim(text), '/')
    if (slash_position == 0) then
      call parse_real_value(text, value, file_name, line_number)
      return
    end if
    if (slash_position == 1 .or. slash_position == len_trim(text) .or. &
        index(text(slash_position + 1:), '/') > 0) then
      call parameter_error(file_name, line_number, &
                           'fraccion de pi invalida: '//trim(text))
    end if

    read(text(:slash_position - 1), *, iostat=io_status) numerator
    if (io_status /= 0) call parameter_error(file_name, line_number, &
                                              'numerador invalido: '//trim(text))
    read(text(slash_position + 1:), *, iostat=io_status) denominator
    if (io_status /= 0 .or. denominator == 0.0d0) then
      call parameter_error(file_name, line_number, &
                           'denominador invalido o cero: '//trim(text))
    end if
    value = numerator / denominator
  end subroutine parse_pi_coefficient


  subroutine parse_logical_value(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    logical, intent(out) :: value
    integer, intent(in) :: line_number
    character(len=MAX_PARAMETER_LINE) :: normalized

    normalized = lowercase(trim(adjustl(text)))
    select case(trim(normalized))
    case('1', 'true', '.true.', 'yes', 'si')
      value = .true.
    case('0', 'false', '.false.', 'no')
      value = .false.
    case default
      call parameter_error(file_name, line_number, &
                           'se esperaba true/false, yes/no o 1/0: '//trim(text))
    end select
  end subroutine parse_logical_value


  subroutine parse_string_value(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    character(len=*), intent(out) :: value
    integer, intent(in) :: line_number
    character(len=MAX_PARAMETER_LINE) :: parsed, unquoted
    integer :: text_length

    parsed = trim(adjustl(text))
    text_length = len_trim(parsed)
    if (text_length >= 2) then
      if ((parsed(1:1) == '"' .and. parsed(text_length:text_length) == '"') .or. &
          (parsed(1:1) == "'" .and. parsed(text_length:text_length) == "'")) then
        unquoted = ''
        if (text_length > 2) then
          unquoted(1:text_length - 2) = parsed(2:text_length - 1)
        end if
        parsed = unquoted
      end if
    end if
    if (len_trim(parsed) == 0) then
      call parameter_error(file_name, line_number, 'la cadena no puede estar vacia')
    end if
    if (len_trim(parsed) > len(value)) then
      call parameter_error(file_name, line_number, 'la cadena excede la longitud admitida')
    end if
    value = trim(parsed)
  end subroutine parse_string_value


  subroutine parse_problem(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    integer, intent(out) :: value
    integer, intent(in) :: line_number
    character(len=MAX_PARAMETER_LINE) :: normalized

    normalized = lowercase(trim(adjustl(text)))
    select case(trim(normalized))
    case('sod')
      value = 1
    case('strong_shock', 'blast_wave')
      value = 2
    case('shu_osher')
      value = 3
    case('kelvin_helmholtz', 'khi')
      value = 4
    case('axisymmetric_jet', 'jet')
      value = 5
    case('convergence', 'convergence_test')
      value = 6
    case('michel', 'michel_accretion')
      value = 7
    case('dust', 'dust_accretion')
      value = 8
    case('offaxis_blast', 'off_axis_blast')
      value = 9
    case('fishbone_equatorial', 'fishbone_moncrief_equatorial')
      value = 10
    case('fishbone_sagittal', 'fishbone_sagital', 'fishbone_moncrief_sagittal')
      value = 11
    case('fishbone_3d', 'fishbone_moncrief_3d')
      value = 12
    case default
      call parameter_error(file_name, line_number, 'problema desconocido: '//trim(text))
    end select
  end subroutine parse_problem


  subroutine parse_metric(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    character(len=*), intent(out) :: value
    integer, intent(in) :: line_number
    character(len=MAX_PARAMETER_LINE) :: normalized

    normalized = lowercase(trim(adjustl(text)))
    select case(trim(normalized))
    case('minkowski')
      value = 'Minkowski'
    case('eddington_finkelstein', 'ef')
      value = 'Eddington-Finkelstein'
    case('kerr_schild', 'ks')
      value = 'Kerr-Schild'
    case default
      call parameter_error(file_name, line_number, 'metrica desconocida: '//trim(text))
    end select
  end subroutine parse_metric


  subroutine parse_geometry(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    character(len=*), intent(out) :: value
    integer, intent(in) :: line_number
    character(len=MAX_PARAMETER_LINE) :: normalized

    normalized = lowercase(trim(adjustl(text)))
    select case(trim(normalized))
    case('cartesian')
      value = 'Cartesian'
    case('cylindrical')
      value = 'Cylindrical'
    case('spherical')
      value = 'Spherical'
    case default
      call parameter_error(file_name, line_number, 'geometria desconocida: '//trim(text))
    end select
  end subroutine parse_geometry


  subroutine parse_vtk_mapping(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    integer, intent(out) :: value
    integer, intent(in) :: line_number

    select case(trim(lowercase(trim(adjustl(text)))))
    case('physical', 'ks_cartesian', 'complete', 'completo')
      value = VTK_MAP_PHYSICAL
    case('untwisted', 'unrolled', 'desenrollado')
      value = VTK_MAP_UNTWISTED
    case default
      call parameter_error(file_name, line_number, &
                           'vtk_mapping debe ser physical o untwisted')
    end select
  end subroutine parse_vtk_mapping


  subroutine parse_reconstruction(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    integer, intent(out) :: value
    integer, intent(in) :: line_number

    select case(trim(lowercase(trim(adjustl(text)))))
    case('godunov'); value = REC_GODUNOV
    case('tvd');     value = REC_TVD
    case('weno3');   value = REC_WENO3
    case('mp5');     value = REC_MP5
    case('weno5');   value = REC_WENO5
    case default
      call parameter_error(file_name, line_number, &
                           'reconstruction debe ser godunov, tvd, weno3, mp5 o weno5')
    end select
  end subroutine parse_reconstruction


  subroutine parse_riemann_solver(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    integer, intent(out) :: value
    integer, intent(in) :: line_number

    select case(trim(lowercase(trim(adjustl(text)))))
    case('hlle'); value = RS_HLLE
    case('hllc'); value = RS_HLLC
    case default
      call parameter_error(file_name, line_number, 'riemann_solver debe ser hlle o hllc')
    end select
  end subroutine parse_riemann_solver


  subroutine parse_tvd_limiter(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    integer, intent(out) :: value
    integer, intent(in) :: line_number

    select case(trim(lowercase(trim(adjustl(text)))))
    case('minmod');   value = LIM_MINMOD
    case('superbee'); value = LIM_SUPERBEE
    case('mc');       value = LIM_MC
    case default
      call parameter_error(file_name, line_number, &
                           'tvd_limiter debe ser minmod, superbee o mc')
    end select
  end subroutine parse_tvd_limiter


  subroutine parse_perturbation(text, value, file_name, line_number)
    character(len=*), intent(in) :: text, file_name
    integer, intent(out) :: value
    integer, intent(in) :: line_number

    select case(trim(lowercase(trim(adjustl(text)))))
    case('none');           value = PERT_NONE
    case('pressure_noise'); value = PERT_PRESSURE_NOISE
    case('density_noise');  value = PERT_DENSITY_NOISE
    case('density_mode');   value = PERT_DENSITY_MODE
    case default
      call parameter_error(file_name, line_number, &
        'perturbation debe ser none, pressure_noise, density_noise o density_mode')
    end select
  end subroutine parse_perturbation


  function problem_label(problem_id) result(label)
    integer, intent(in) :: problem_id
    character(len=40) :: label

    select case(problem_id)
    case(1);  label = 'Sod'
    case(2);  label = 'Strong shock'
    case(3);  label = 'Shu-Osher'
    case(4);  label = 'Kelvin-Helmholtz'
    case(5);  label = 'Axisymmetric jet'
    case(6);  label = 'Convergence test'
    case(7);  label = 'Michel accretion'
    case(8);  label = 'Dust accretion'
    case(9);  label = 'Off-axis blast'
    case(10); label = 'Fishbone-Moncrief equatorial'
    case(11); label = 'Fishbone-Moncrief sagittal'
    case(12); label = 'Fishbone-Moncrief 3D'
    case default; label = 'unknown'
    end select
  end function problem_label


  pure function lowercase(text) result(lowered)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: lowered
    integer :: i, code

    lowered = text
    do i = 1, len(text)
      code = iachar(text(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        lowered(i:i) = achar(code + iachar('a') - iachar('A'))
      end if
    end do
  end function lowercase


  subroutine parameter_error(file_name, line_number, message)
    character(len=*), intent(in) :: file_name, message
    integer, intent(in) :: line_number

    if (line_number > 0) then
      write(*,'(A,1X,A,A,I0,A,1X,A)') 'ERROR en', trim(file_name), ':', &
                                      line_number, ':', trim(message)
    else
      write(*,'(A,1X,A,A,1X,A)') 'ERROR en', trim(file_name), ':', trim(message)
    end if
    error stop 1
  end subroutine parameter_error

end module parameters
