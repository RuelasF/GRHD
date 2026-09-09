! =======================================================================================
! Módulo: output (Visualización, Checkpoints y Extracción de Observables)
! ---------------------------------------------------------------------------------------
! Propósito: Gestiona toda la salida de datos del simulador hacia el disco duro.
! Características Clave:
! 1. Exportación VTK Binario: Máxima velocidad de escritura (I/O) minimizando el 
!    cuello de botella del disco duro. Genera archivos compatibles con VisIt y ParaView.
! 2. Tolerancia a Fallos: Sistema de Checkpoints que guarda un volcado binario exacto 
!    de la memoria RAM para pausar y reanudar simulaciones de semanas de duración.
! 3. Física Multimensajero: Extracción al vuelo de la señal de Ondas Gravitacionales 
!    (polarizaciones h+ y hx) y la tasa de acreción de masa al agujero negro.
! =======================================================================================
module output
  use variables
  use metrics
  use conditions
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  public :: save_vtk_at_time, save_checkpoint, read_checkpoint_metadata, &
            read_checkpoint_state, calc_gw_strain, calc_m_dot, calc_convergence_norms, &
            initialize_auxiliary_output, calc_ppi_modes, calc_global_diagnostics, &
            write_perturbation_event, finn_evans_cell_stress, gw_cartesian_kinematics, &
            kerr_schild_spheroidal_position

  ! La extracción Finn--Evans es un proxy de campo débil evaluado sobre el
  ! fluido GRHD. La distancia está expresada en las unidades geométricas del
  ! código y se registra en el manifiesto de cada corrida.
  real*8, parameter :: gw_observer_distance = 1000.0d0
  logical, parameter :: gw_exclude_inside_horizon = .true.
  logical, parameter :: gw_subtract_numerical_atmosphere = .true.

contains

  subroutine initialize_text_file(file_name, header, replace_existing)
    character(len=*), intent(in) :: file_name, header
    logical, intent(in) :: replace_existing
    logical :: file_exists
    integer :: unit_file

    inquire(file=trim(file_name), exist=file_exists)
    if (file_exists .and. .not. replace_existing) return

    open(newunit=unit_file, file=trim(file_name), status='replace', action='write')
    write(unit_file, '(A)') trim(header)
    close(unit_file)
  end subroutine initialize_text_file


  subroutine write_run_manifest()
    integer :: unit_file
    character(len=256) :: file_name
    real*8 :: manifest_spin

    file_name = trim(output_folder) // '/' // trim(scheme_name) // '/run_manifest.txt'
    open(newunit=unit_file, file=trim(file_name), status='unknown', &
         position='append', action='write')
    write(unit_file, '(A)') '[run]'
    write(unit_file, '(A,ES24.16E3)') 'start_time=', integration_time
    write(unit_file, '(A,A)') 'case=', trim(case_name)
    write(unit_file, '(A,A)') 'metric=', trim(metric_type)
    write(unit_file, '(A,A)') 'geometry=', trim(geom_type)
    if (vtk_mapping_id == VTK_MAP_PHYSICAL) then
      write(unit_file, '(A)') 'vtk_mapping=physical'
    else
      write(unit_file, '(A)') 'vtk_mapping=untwisted'
    end if
    write(unit_file, '(A,ES24.16E3)') 'bh_mass=', bh_mass
    write(unit_file, '(A,ES24.16E3)') 'a_spin=', a_spin
    if (trim(metric_type) /= 'Minkowski' .and. bh_mass > 0.0d0) then
      manifest_spin = 0.0d0
      if (trim(metric_type) == 'Kerr-Schild') manifest_spin = a_spin
      write(unit_file, '(A,ES24.16E3)') 'event_horizon_radius=', &
                                       event_horizon_radius(bh_mass, manifest_spin)
    end if
    write(unit_file, '(A,A)') 'scheme=', trim(scheme_name)
    write(unit_file, '(A,3(I0,1X))') 'grid_dimensions=', nx, ny, nz
    write(unit_file, '(A,2(ES24.16E3,1X))') 'radial_domain=', r_min, r_max
    write(unit_file, '(A,L1)') 'logarithmic_radial_grid=', use_log_r
    write(unit_file, '(A,ES24.16E3)') 'final_time=', final_time
    write(unit_file, '(A,ES24.16E3)') 'cfl=', CFL
    write(unit_file, '(A,L1)') 'ppi_diagnostics=', do_ppi_diagnostics
    write(unit_file, '(A,I0)') 'diagnostic_stride=', diagnostic_stride
    if (do_gw_extraction) then
      write(unit_file, '(A)') 'gw_method=Finn-Evans weak-field stress proxy'
      write(unit_file, '(A)') 'gw_observer_axis=positive-z'
      write(unit_file, '(A,ES24.16E3)') 'gw_observer_distance=', gw_observer_distance
      write(unit_file, '(A,L1)') 'gw_exclude_inside_horizon=', gw_exclude_inside_horizon
      write(unit_file, '(A,L1)') 'gw_subtract_numerical_atmosphere=', &
                                 gw_subtract_numerical_atmosphere
      write(unit_file, '(A)') 'gw_volume_measure=flat Cartesian Kerr-Schild Jacobian'
      write(unit_file, '(A)') 'gw_density_measure=Newtonian rest-mass density proxy'
      write(unit_file, '(A,I0)') 'gw_polar_cells=', ny
    end if
    write(unit_file, '(A,L1)') 'apply_perturbation=', apply_perturbation
    write(unit_file, '(A,L1)') 'perturbation_already_applied=', perturbation_applied
    write(unit_file, '(A,I0)') 'perturbation_type=', perturbation_type
    write(unit_file, '(A,I0)') 'perturbation_seed=', perturbation_seed
    write(unit_file, '(A,ES24.16E3)') 'perturbation_time=', perturbation_time
    write(unit_file, '(A,ES24.16E3)') 'perturbation_amplitude=', perturbation_amplitude
    write(unit_file, '(A,ES24.16E3)') 'perturbation_mode=', perturbation_mode
    write(unit_file, '(A)') '[/run]'
    close(unit_file)
  end subroutine write_run_manifest


  subroutine initialize_auxiliary_output(replace_existing)
    logical, intent(in) :: replace_existing
    integer :: exit_status
    character(len=256) :: directory, file_name

    directory = trim(output_folder) // '/' // trim(scheme_name)
    call execute_command_line('mkdir -p "' // trim(directory) // '"', &
                              exitstat=exit_status)
    if (exit_status /= 0) then
      write(*,*) 'CRITICAL ERROR: unable to create output directory: ', trim(directory)
      error stop 1
    end if

    if (do_gw_extraction) then
      file_name = trim(directory) // '/GW_signal.dat'
      call initialize_text_file(file_name, '# Time h_plus h_cross', replace_existing)
      if (ny == 1) then
        write(*,*) 'WARNING: GW amplitudes with ny=1 retain the single-cell polar quadrature.'
        write(*,*) '         Frequencies/phases are meaningful proxies; absolute 3D normalization is not fixed.'
      end if
    end if
    if (do_mdot_extraction) then
      file_name = trim(directory) // '/m_dot.dat'
      call initialize_text_file(file_name, '# Time M_dot', replace_existing)
    end if
    if (do_ppi_diagnostics) then
      file_name = trim(directory) // '/ppi_modes.dat'
      call initialize_text_file(file_name, &
        '# t C0 ReC1 ImC1 A1 phase1 ReC2 ImC2 A2 phase2 ReC3 ImC3 A3 phase3 ReC4 ImC4 A4 phase4', &
        replace_existing)
      file_name = trim(directory) // '/global_diagnostics.dat'
      call initialize_text_file(file_name, &
        '# t M_total M_disk J_phi rho_min p_min rho_max p_max v2_max N_atmosphere N_invalid', &
        replace_existing)
      file_name = trim(directory) // '/perturbation_events.dat'
      call initialize_text_file(file_name, &
        '# t type seed amplitude mode', replace_existing)
    end if

    call write_run_manifest()
  end subroutine initialize_auxiliary_output

  ! ===================================================================================
  ! 1. EXPORTADOR VTK (VisIt/ParaView) - FORMATO BINARIO
  ! ===================================================================================
  subroutine save_vtk_data(file_name)
    character(len=*), intent(in) :: file_name
    integer :: i, j, k, unit_file, io_status
    integer :: nx_points, ny_points, nz_points, n_points
    character(len=150) :: full_path
    character(len=100) :: header_line
    real*8, allocatable :: pts_buffer(:,:,:,:)
    real*8 :: r_phys, theta_point, phi_point
    real*8 :: geometry_spin, cartesian_position(3)
    real*8 :: topology_tolerance
    logical :: use_polar_unstructured_grid

    geometry_spin = 0.0d0
    if (trim(metric_type) == 'Kerr-Schild') geometry_spin = a_spin

    ! Una esfera completa no puede representarse con hexaedros estructurados
    ! regulares en los polos: allí distintos índices azimutales ocupan el mismo
    ! punto y VisIt puede producir nudos al cortar esas celdas degeneradas.
    topology_tolerance = 128.0d0 * epsilon(1.0d0)
    use_polar_unstructured_grid = trim(geom_type) == 'Spherical' .and. &
      ny > 1 .and. nz > 1 .and. &
      abs(y_face(0)) <= topology_tolerance .and. &
      abs(y_face(ny) - pi) <= topology_tolerance * max(1.0d0, abs(pi)) .and. &
      abs((z_face(nz) - z_face(0)) - 2.0d0*pi) <= &
        topology_tolerance * max(1.0d0, abs(2.0d0*pi))

    if (use_polar_unstructured_grid) then
      call save_vtk_spherical_unstructured(file_name, geometry_spin)
      return
    end if

    ! Las variables hidrodinámicas viven en celdas. La malla VTK se construye
    ! con sus caras verdaderas para que cada celda numérica corresponda a una
    ! celda geométrica y los polos no se formen con anillos artificiales.
    nx_points = nx + 1
    ny_points = 1
    nz_points = 1
    if (ny > 1) ny_points = ny + 1
    if (nz > 1) nz_points = nz + 1
    n_points = nx_points * ny_points * nz_points

    full_path = trim(output_folder) //  '/' // trim(scheme_name) // '/' // trim(file_name)
    print *, '>>> Guardando frame VTK Binario: ', trim(full_path)

    unit_file = 10
    open(unit_file, file=full_path, status='replace', access='stream', form='unformatted', convert='big_endian', iostat=io_status)
    if (io_status /= 0) then
      print *, 'CRITICAL ERROR: No se pudo crear el archivo VTK: ', trim(full_path)
      return
    end if

    write(unit_file) '# vtk DataFile Version 3.0' // char(10)
    write(unit_file) 'GRHD Simulation Data' // char(10)
    write(unit_file) 'BINARY' // char(10)   
    write(unit_file) 'DATASET STRUCTURED_GRID' // char(10)

    write(header_line, '(a, i0, 1x, i0, 1x, i0)') &
      'DIMENSIONS ', nx_points, ny_points, nz_points
    write(unit_file) trim(header_line) // char(10)

    write(header_line, '(a, i0, a)') 'POINTS ', n_points, ' double'
    write(unit_file) trim(header_line) // char(10)

    allocate(pts_buffer(3, 1:nx_points, 1:ny_points, 1:nz_points))
    
    ! -------------------------------------------------------------
    ! BUCLE DE MAPEO GEOMÉTRICO 3D (X, Y, Z)
    ! -------------------------------------------------------------
    do k = 1, nz_points
      if (nz > 1) then
        phi_point = z_face(k - 1)
      else
        phi_point = z(1)
      end if

      do j = 1, ny_points
        if (ny > 1) then
          theta_point = y_face(j - 1)
        else
          theta_point = y(1)
        end if

        do i = 1, nx_points
          if (use_log_r) then
            r_phys = exp(x_face(i - 1))
          else
            r_phys = x_face(i - 1)
          end if

          select case(trim(geom_type))
            case('Spherical')
              call vtk_spherical_position(r_phys, theta_point, phi_point, &
                                           geometry_spin, cartesian_position)
              pts_buffer(:, i, j, k) = cartesian_position

            case('Cylindrical')
              pts_buffer(1, i, j, k) = r_phys * cos(theta_point)
              pts_buffer(2, i, j, k) = r_phys * sin(theta_point)
              pts_buffer(3, i, j, k) = phi_point

            case default
              pts_buffer(1, i, j, k) = r_phys
              pts_buffer(2, i, j, k) = theta_point
              pts_buffer(3, i, j, k) = phi_point
          end select
        end do
      end do
    end do
    
    write(unit_file) pts_buffer
    write(unit_file) char(10) 
    deallocate(pts_buffer)

    call write_vtk_cell_data(unit_file)
    close(unit_file)
  end subroutine save_vtk_data


  subroutine save_vtk_spherical_unstructured(file_name, geometry_spin)
    character(len=*), intent(in) :: file_name
    real*8, intent(in) :: geometry_spin
    integer :: i, j, k, k_next, point_id, cell_id, offset
    integer :: unit_file, io_status, n_points, n_cells, cell_list_size
    integer(kind=4), allocatable :: ring_id(:,:,:), north_id(:), south_id(:)
    integer(kind=4), allocatable :: cell_buffer(:), cell_types(:)
    real*8, allocatable :: point_buffer(:,:)
    real*8 :: r_phys, position(3)
    character(len=150) :: full_path
    character(len=100) :: header_line

    ! Hay nz puntos por anillo: la conectividad periódica une el último con el
    ! primero. Cada cara radial tiene un único punto norte y uno sur.
    n_points = (nx + 1) * ((ny - 1) * nz + 2)
    n_cells = nx * ny * nz
    cell_list_size = nx * nz * (9 * ny - 4)

    allocate(point_buffer(3, n_points))
    allocate(ring_id(0:nx, 1:ny-1, 0:nz-1))
    allocate(north_id(0:nx), south_id(0:nx))

    point_id = 0
    do i = 0, nx
      if (use_log_r) then
        r_phys = exp(x_face(i))
      else
        r_phys = x_face(i)
      end if

      north_id(i) = point_id
      call vtk_spherical_position(r_phys, 0.0d0, z_face(0), geometry_spin, position)
      point_buffer(:, point_id + 1) = position
      point_id = point_id + 1

      do k = 0, nz - 1
        do j = 1, ny - 1
          ring_id(i, j, k) = point_id
          call vtk_spherical_position(r_phys, y_face(j), z_face(k), &
                                      geometry_spin, position)
          point_buffer(:, point_id + 1) = position
          point_id = point_id + 1
        end do
      end do

      south_id(i) = point_id
      call vtk_spherical_position(r_phys, pi, z_face(0), geometry_spin, position)
      point_buffer(:, point_id + 1) = position
      point_id = point_id + 1
    end do

    if (point_id /= n_points) error stop 'Internal VTK point-count mismatch.'

    allocate(cell_buffer(cell_list_size), cell_types(n_cells))
    offset = 1
    cell_id = 0
    do k = 0, nz - 1
      k_next = mod(k + 1, nz)
      do j = 1, ny
        do i = 0, nx - 1
          cell_id = cell_id + 1
          if (j == 1) then
            ! Prisma polar norte: caras triangulares a radio constante.
            cell_buffer(offset:offset+6) = [6, north_id(i), &
              ring_id(i, 1, k), ring_id(i, 1, k_next), north_id(i+1), &
              ring_id(i+1, 1, k), ring_id(i+1, 1, k_next)]
            cell_types(cell_id) = 13       ! VTK_WEDGE
            offset = offset + 7
          else if (j == ny) then
            ! En el polo sur se invierte el orden azimutal para conservar la
            ! orientación de la celda.
            cell_buffer(offset:offset+6) = [6, south_id(i), &
              ring_id(i, ny-1, k_next), ring_id(i, ny-1, k), south_id(i+1), &
              ring_id(i+1, ny-1, k_next), ring_id(i+1, ny-1, k)]
            cell_types(cell_id) = 13       ! VTK_WEDGE
            offset = offset + 7
          else
            cell_buffer(offset:offset+8) = [8, ring_id(i, j-1, k), &
              ring_id(i+1, j-1, k), ring_id(i+1, j, k), ring_id(i, j, k), &
              ring_id(i, j-1, k_next), ring_id(i+1, j-1, k_next), &
              ring_id(i+1, j, k_next), ring_id(i, j, k_next)]
            cell_types(cell_id) = 12       ! VTK_HEXAHEDRON
            offset = offset + 9
          end if
        end do
      end do
    end do

    if (cell_id /= n_cells .or. offset - 1 /= cell_list_size) then
      error stop 'Internal VTK cell-count mismatch.'
    end if

    full_path = trim(output_folder) // '/' // trim(scheme_name) // '/' // trim(file_name)
    print *, '>>> Guardando frame VTK Binario: ', trim(full_path)
    unit_file = 10
    open(unit_file, file=full_path, status='replace', access='stream', &
         form='unformatted', convert='big_endian', iostat=io_status)
    if (io_status /= 0) then
      print *, 'CRITICAL ERROR: No se pudo crear el archivo VTK: ', trim(full_path)
      return
    end if

    write(unit_file) '# vtk DataFile Version 3.0' // char(10)
    write(unit_file) 'GRHD Simulation Data' // char(10)
    write(unit_file) 'BINARY' // char(10)
    write(unit_file) 'DATASET UNSTRUCTURED_GRID' // char(10)
    write(header_line, '(a, i0, a)') 'POINTS ', n_points, ' double'
    write(unit_file) trim(header_line) // char(10)
    write(unit_file) point_buffer
    write(unit_file) char(10)

    write(header_line, '(a, i0, 1x, i0)') 'CELLS ', n_cells, cell_list_size
    write(unit_file) trim(header_line) // char(10)
    write(unit_file) cell_buffer
    write(unit_file) char(10)
    write(header_line, '(a, i0)') 'CELL_TYPES ', n_cells
    write(unit_file) trim(header_line) // char(10)
    write(unit_file) cell_types
    write(unit_file) char(10)

    call write_vtk_cell_data(unit_file)
    close(unit_file)
    deallocate(point_buffer, ring_id, north_id, south_id, cell_buffer, cell_types)
  end subroutine save_vtk_spherical_unstructured


  subroutine write_vtk_cell_data(unit_file)
    integer, intent(in) :: unit_file
    integer :: eq, i, j, k, n_cells
    character(len=100) :: header_line
    real*8, allocatable :: var_buffer(:,:,:)

    n_cells = nx * ny * nz
    write(header_line, '(a, i0)') 'CELL_DATA ', n_cells
    write(unit_file) trim(header_line) // char(10)
    allocate(var_buffer(1:nx, 1:ny, 1:nz))

    do eq = 1, neq
      write(header_line, '(a, 1x, a, a)') 'SCALARS', trim(var_names(eq)), ' double 1'
      write(unit_file) trim(header_line) // char(10)
      write(unit_file) 'LOOKUP_TABLE default' // char(10)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            var_buffer(i, j, k) = p(eq, i, j, k)
            if (use_log_r .and. eq == eq_vx) then
              var_buffer(i, j, k) = var_buffer(i, j, k) * exp(x(i))
            end if
          end do
        end do
      end do
      write(unit_file) var_buffer
      write(unit_file) char(10)
    end do

    write(unit_file) 'SCALARS Lapse_alpha double 1' // char(10)
    write(unit_file) 'LOOKUP_TABLE default' // char(10)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          var_buffer(i, j, k) = alpha_c(i, j, k)
        end do
      end do
    end do
    write(unit_file) var_buffer
    write(unit_file) char(10)

    write(unit_file) 'SCALARS Shift_beta_r double 1' // char(10)
    write(unit_file) 'LOOKUP_TABLE default' // char(10)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          var_buffer(i, j, k) = beta_c(1, i, j, k)
          if (use_log_r) var_buffer(i, j, k) = var_buffer(i, j, k) * exp(x(i))
        end do
      end do
    end do
    write(unit_file) var_buffer
    write(unit_file) char(10)

    deallocate(var_buffer)
  end subroutine write_vtk_cell_data


  subroutine save_vtk_at_time(current_time)
    real*8, intent(in) :: current_time
    character(len=50) :: file_name
    integer :: frame_id

    if (save_interval > 0.0d0) then
      frame_id = nint(current_time / save_interval)
    else
      frame_id = 0
    end if

    write(file_name, '(A, "_", A, "_", I4.4, ".vtk")') trim(output_prefix), trim(scheme_name), frame_id
    call save_vtk_data(file_name)
  end subroutine save_vtk_at_time

  ! ===================================================================================
  ! 2. SISTEMA DE CHECKPOINTS (Gestión de memoria binaria directa)
  ! ===================================================================================
  
  subroutine save_checkpoint(step_num, current_time)
    integer, intent(in) :: step_num
    real*8, intent(in) :: current_time
    integer :: unit_file, io_status
    character(len=100) :: file_name
    character(len=12), parameter :: checkpoint_magic = 'GRHDCP_V3'

    write(file_name, '(A, "/checkpoint_", A, "_", I5.5, ".rst")') trim(output_folder), trim(scheme_name), int(current_time)
    print *, '>>> Guardando Checkpoint de Respaldo: ', trim(file_name)

    unit_file = 20
    open(unit_file, file=file_name, status='replace', action='write', form='unformatted', iostat=io_status)
    if (io_status == 0) then
      write(unit_file) checkpoint_magic
      ! Topología de malla 3D y profundidad de fantasmas
      write(unit_file) nghost, nx, ny, nz 
      
      ! Fronteras lógicas (Actualizado a r_min/r_max)
      write(unit_file) r_min, r_max, y_min, y_max, z_min, z_max
      
      ! Física y Banderas Arquitectónicas (Actualizado con banderas de extracción)
      write(unit_file) bh_mass, adb_idx, a_spin
      write(unit_file) use_log_r, use_shock_sensor, do_gw_extraction, &
                       do_mdot_extraction, do_ppi_diagnostics
      write(unit_file) metric_type, geom_type
      
      ! Control de tiempo y salida
      write(unit_file) final_time, CFL, save_interval
      write(unit_file) output_prefix, output_folder
      
      ! Identificadores de esquemas
      write(unit_file) case_id, case_name, scheme_name
      write(unit_file) rec_method_id, tvd_limiter_id, riemann_solver_id

      ! Estado reproducible de la perturbación y frecuencia de diagnóstico
      write(unit_file) apply_perturbation, perturbation_applied
      write(unit_file) perturbation_type, perturbation_seed, diagnostic_stride
      write(unit_file) perturbation_time, perturbation_amplitude, perturbation_mode
      
      ! Estado temporal de la simulación
      write(unit_file) current_time
      write(unit_file) step_num
      
      ! Vuelco binario masivo del arreglo de estado conservado 4D
      write(unit_file) up 
      close(unit_file)
    else
      print *, 'CRITICAL ERROR: Fallo al crear el archivo de checkpoint.'
    end if
  end subroutine save_checkpoint

  subroutine read_checkpoint_metadata(unit_file)
    integer, intent(in) :: unit_file
    integer :: io_status
    logical :: is_v2, is_v3, is_versioned
    character(len=12) :: checkpoint_magic
    character(len=12), parameter :: checkpoint_magic_v2 = 'GRHDCP_V2'
    character(len=12), parameter :: checkpoint_magic_v3 = 'GRHDCP_V3'

    read(unit_file, iostat=io_status) checkpoint_magic
    is_v2 = (io_status == 0 .and. checkpoint_magic == checkpoint_magic_v2)
    is_v3 = (io_status == 0 .and. checkpoint_magic == checkpoint_magic_v3)
    is_versioned = is_v2 .or. is_v3
    if (.not. is_versioned) rewind(unit_file)

    ! Topología de malla
    read(unit_file) nghost, nx, ny, nz
    
    ! Fronteras lógicas
    read(unit_file) r_min, r_max, y_min, y_max, z_min, z_max
    
    ! Física y Banderas Arquitectónicas (Actualizado con banderas de extracción)
    if (is_versioned) then
      read(unit_file) bh_mass, adb_idx, a_spin
    else
      read(unit_file) bh_mass, adb_idx
      a_spin = 0.0d0
      print *, 'WARNING: legacy checkpoint has no Kerr spin; assuming a = 0.'
    end if
    if (is_v3) then
      read(unit_file) use_log_r, use_shock_sensor, do_gw_extraction, &
                      do_mdot_extraction, do_ppi_diagnostics
    else
      read(unit_file) use_log_r, use_shock_sensor, do_gw_extraction, do_mdot_extraction
      do_ppi_diagnostics = .false.
    end if
    read(unit_file) metric_type, geom_type
    
    ! Control
    read(unit_file) final_time, CFL, save_interval
    read(unit_file) output_prefix, output_folder
    
    ! Eschemas
    read(unit_file) case_id, case_name, scheme_name
    read(unit_file) rec_method_id, tvd_limiter_id, riemann_solver_id

    if (is_v3) then
      read(unit_file) apply_perturbation, perturbation_applied
      read(unit_file) perturbation_type, perturbation_seed, diagnostic_stride
      read(unit_file) perturbation_time, perturbation_amplitude, perturbation_mode
    else
      perturbation_applied = .false.
    end if
  end subroutine read_checkpoint_metadata

  subroutine read_checkpoint_state(unit_file, step_num, current_time)
    integer, intent(in) :: unit_file
    integer, intent(out) :: step_num
    real*8, intent(out) :: current_time

    read(unit_file) current_time
    read(unit_file) step_num
    
    ! El programa principal debe asegurar que 'up' ya esté alojado en memoria 
    ! dinámicamente con nx, ny, nz leídos en la metadata antes de llegar aquí.
    read(unit_file) up
  end subroutine read_checkpoint_state

  ! ===================================================================================
  ! 3. EXTRACCIÓN DE OBSERVABLES FÍSICOS (Física de Agujeros Negros)
  ! ===================================================================================

  pure subroutine kerr_schild_spheroidal_position(r_phys, theta, phi, spin, position)
    real*8, intent(in) :: r_phys, theta, phi, spin
    real*8, intent(out) :: position(3)
    real*8 :: sin_theta, cos_theta, sin_phi, cos_phi

    sin_theta = sin(theta)
    cos_theta = cos(theta)

    ! sin(pi) no es exactamente cero en punto flotante. En Kerr-Schild ese
    ! residuo genera un anillo diminuto cuya orientación depende de r y spin;
    ! VisIt puede mostrarlo como un nudo en el polo. Colapsamos ambos ejes de
    ! coordenadas de manera exacta.
    if (abs(sin_theta) <= 64.0d0 * epsilon(1.0d0)) then
      position(1) = 0.0d0
      position(2) = 0.0d0
      position(3) = sign(r_phys, cos_theta)
      return
    end if

    sin_phi = sin(phi)
    cos_phi = cos(phi)

    position(1) = (r_phys*cos_phi - spin*sin_phi) * sin_theta
    position(2) = (r_phys*sin_phi + spin*cos_phi) * sin_theta
    position(3) = r_phys*cos_theta
  end subroutine kerr_schild_spheroidal_position


  subroutine vtk_spherical_position(r_phys, theta, phi, spin, position)
    real*8, intent(in) :: r_phys, theta, phi, spin
    real*8, intent(out) :: position(3)
    real*8 :: sin_theta, cos_theta, cylindrical_radius

    if (vtk_mapping_id == VTK_MAP_PHYSICAL) then
      call kerr_schild_spheroidal_position(r_phys, theta, phi, spin, position)
      return
    end if

    sin_theta = sin(theta)
    cos_theta = cos(theta)
    if (abs(sin_theta) <= 64.0d0 * epsilon(1.0d0)) then
      position = [0.0d0, 0.0d0, sign(r_phys, cos_theta)]
      return
    end if

    ! Visualizacion oblata sin la rotacion radial atan(a/r). Conserva las
    ! superficies r=constante y colapsa exactamente ambos polos, pero no cambia
    ! las coordenadas ni los datos usados durante la evolucion.
    cylindrical_radius = sqrt(r_phys**2 + spin**2) * sin_theta
    position(1) = cylindrical_radius * cos(phi)
    position(2) = cylindrical_radius * sin(phi)
    position(3) = r_phys * cos_theta
  end subroutine vtk_spherical_position

  pure subroutine gw_cartesian_kinematics(r_phys, theta, phi, spin, alpha, beta, &
                                          eulerian_velocity, logarithmic_r, &
                                          position, velocity, flat_jacobian)
    real*8, intent(in) :: r_phys, theta, phi, spin, alpha
    real*8, intent(in) :: beta(3), eulerian_velocity(3)
    logical, intent(in) :: logarithmic_r
    real*8, intent(out) :: position(3), velocity(3), flat_jacobian
    real*8 :: coordinate_rates(3)
    real*8 :: sin_theta, cos_theta, sin_phi, cos_phi, sigma
    real*8 :: r_dot, theta_dot, phi_dot

    ! Las primitivas son velocidades eulerianas. Finn--Evans necesita las
    ! tasas coordenadas para construir velocidades cartesianas.
    coordinate_rates = alpha * eulerian_velocity - beta
    if (logarithmic_r) coordinate_rates(1) = r_phys * coordinate_rates(1)

    sin_theta = sin(theta)
    cos_theta = cos(theta)
    sin_phi = sin(phi)
    cos_phi = cos(phi)
    r_dot = coordinate_rates(1)
    theta_dot = coordinate_rates(2)
    phi_dot = coordinate_rates(3)

    ! Para spin=0 se recupera la transformación esférica usada por EF.
    call kerr_schild_spheroidal_position(r_phys, theta, phi, spin, position)

    velocity(1) = r_dot * cos_phi * sin_theta &
                + theta_dot * (r_phys * cos_phi - spin * sin_phi) * cos_theta &
                - phi_dot * (r_phys * sin_phi + spin * cos_phi) * sin_theta
    velocity(2) = r_dot * sin_phi * sin_theta &
                + theta_dot * (r_phys * sin_phi + spin * cos_phi) * cos_theta &
                + phi_dot * (r_phys * cos_phi - spin * sin_phi) * sin_theta
    velocity(3) = r_dot * cos_theta - r_phys * theta_dot * sin_theta

    sigma = r_phys**2 + spin**2 * cos_theta**2
    flat_jacobian = sigma * sin_theta
  end subroutine gw_cartesian_kinematics


  pure subroutine finn_evans_cell_stress(rho, pressure, central_mass_radius3, &
                                         position, velocity, cell_volume, &
                                         Ixx_cell, Iyy_cell, Ixy_cell)
    real*8, intent(in) :: rho, pressure, central_mass_radius3, cell_volume
    real*8, intent(in) :: position(3), velocity(3)
    real*8, intent(out) :: Ixx_cell, Iyy_cell, Ixy_cell

    Ixx_cell = (2.0d0 * rho * velocity(1)**2 + 2.0d0 * pressure &
               -2.0d0 * rho * central_mass_radius3 * position(1)**2) * cell_volume
    Iyy_cell = (2.0d0 * rho * velocity(2)**2 + 2.0d0 * pressure &
               -2.0d0 * rho * central_mass_radius3 * position(2)**2) * cell_volume
    Ixy_cell = (2.0d0 * rho * velocity(1) * velocity(2) &
               -2.0d0 * rho * central_mass_radius3 * position(1) * position(2)) * cell_volume
  end subroutine finn_evans_cell_stress

  subroutine calc_gw_strain(current_time)
    real*8, intent(in) :: current_time
    integer :: i, j, k
    real*8 :: Ixx_ddot, Iyy_ddot, Ixy_ddot
    real*8 :: r, th, phi, rho, pres, dV, mass_radius3
    real*8 :: h_plus, h_cross
    real*8 :: alpha, beta(3), eulerian_velocity(3)
    real*8 :: position(3), cartesian_velocity(3), flat_jacobian
    real*8 :: cartesian_radius2, radial_grid_jacobian
    real*8 :: horizon_radius, geometry_spin
    real*8 :: Ixx_cell, Iyy_cell, Ixy_cell

    if (trim(metric_type) == 'Minkowski' .or. bh_mass < 0.1d0) return

    ! EF es el caso Schwarzschild a=0. KS usa las coordenadas oblatas
    ! asociadas con el spin configurado en la métrica.
    geometry_spin = 0.0d0
    if (trim(metric_type) == 'Kerr-Schild') geometry_spin = a_spin
    horizon_radius = event_horizon_radius(bh_mass, geometry_spin)

    Ixx_ddot = 0.0d0
    Iyy_ddot = 0.0d0
    Ixy_ddot = 0.0d0

    !$OMP PARALLEL DO PRIVATE(i, j, k, r, th, phi, rho, pres, dV, mass_radius3, &
    !$OMP alpha, beta, eulerian_velocity, position, &
    !$OMP cartesian_velocity, flat_jacobian, cartesian_radius2, radial_grid_jacobian) &
    !$OMP PRIVATE(Ixx_cell, Iyy_cell, Ixy_cell) &
    !$OMP REDUCTION(+:Ixx_ddot, Iyy_ddot, Ixy_ddot)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (use_log_r) then
            r = exp(x(i))
          else
            r = x(i)
          end if

          ! El proxy representa únicamente la materia del disco exterior. La
          ! stress formula newtoniana no modela la desconexión causal del
          ! horizonte; incluir celdas interiores daría una fuente instantánea
          ! manifiestamente no física.
          if (gw_exclude_inside_horizon .and. r <= horizon_radius) cycle

          th = y(j)
          phi = z(k)

          rho = p(eq_de, i, j, k)
          pres = p(eq_pr, i, j, k)
          if (gw_subtract_numerical_atmosphere) then
            rho = max(0.0d0, rho - rho_floor)
            pres = max(0.0d0, pres - p_floor)
          end if
          if (rho <= 0.0d0 .and. pres <= 0.0d0) cycle

          alpha = alpha_c(i,j,k)
          beta = beta_c(:,i,j,k)
          eulerian_velocity = p(eq_vx:eq_vz, i, j, k)
          call gw_cartesian_kinematics(r, th, phi, geometry_spin, alpha, beta, &
                                       eulerian_velocity, use_log_r, position, &
                                       cartesian_velocity, flat_jacobian)

          ! d^3X usa el Jacobiano cartesiano plano de las coordenadas KS,
          ! tal como requiere la aproximación cuadrupolar de campo débil.
          radial_grid_jacobian = 1.0d0
          if (use_log_r) radial_grid_jacobian = r
          dV = flat_jacobian * radial_grid_jacobian * dx * dy * dz

          ! A orden newtoniano el potencial exterior sólo contiene el
          ! monopolo de masa. El spin afecta la dinámica GRHD y la geometría
          ! de coordenadas, pero no se inventa aquí un potencial escalar de
          ! frame dragging.
          cartesian_radius2 = sum(position**2)
          mass_radius3 = bh_mass / (cartesian_radius2 * sqrt(cartesian_radius2))

          call finn_evans_cell_stress(rho, pres, mass_radius3, position, &
                                      cartesian_velocity, dV, Ixx_cell, &
                                      Iyy_cell, Ixy_cell)
          Ixx_ddot = Ixx_ddot + Ixx_cell
          Iyy_ddot = Iyy_ddot + Iyy_cell
          Ixy_ddot = Ixy_ddot + Ixy_cell
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    h_plus  = (Ixx_ddot - Iyy_ddot) / gw_observer_distance
    h_cross = (2.0d0 * Ixy_ddot) / gw_observer_distance

    block
      character(len=150) :: gw_filename
      gw_filename = trim(output_folder) //  '/' // trim(scheme_name) // '/GW_signal.dat'
      open(30, file=gw_filename, position='append', status='unknown')
      write(30, '(3(E15.7, 2X))') current_time, h_plus, h_cross
      close(30)
    end block
  end subroutine calc_gw_strain


  subroutine calc_m_dot(current_time)
    implicit none
    real*8, intent(in) :: current_time
    real*8 :: m_dot
    real*8 :: alpha, beta(3), v_eff
    integer :: j, k

    ! Corrección: Actualización de bandera
    if (trim(metric_type) == 'Minkowski' .or. bh_mass < 0.1d0) return
    if (idx_probe < 1 .or. idx_probe > nx) then
      write(*,*) 'CRITICAL ERROR: invalid accretion-probe index: ', idx_probe
      error stop 1
    end if

    m_dot = 0.0d0
    ! Tasa de acreción integrada sobre una esfera a un radio fijo (idx_probe)
    ! Se propaga la coordenada polar j
    !$OMP PARALLEL DO REDUCTION(+:m_dot) PRIVATE(j, k, alpha, beta, v_eff)
    do k = 1, nz
      do j = 1, ny
        ! call calculate_metric(x(idx_probe), y(j), alpha=alpha, beta=beta)
        alpha = alpha_c(idx_probe, j, k)
        beta(:) = beta_c(:, idx_probe, j, k)

        v_eff = alpha * p(eq_vx, idx_probe, j, k) - beta(1)
        
        ! 'up' contiene el estado conservado densificado. up(eq_de) es D = rho * W * sqrt(gamma)
        m_dot = m_dot + up(eq_de, idx_probe, j, k) * v_eff * dy * dz
      end do
    end do
    !$OMP END PARALLEL DO
    
    m_dot = - m_dot 

    block
      character(len=150) :: m_dot_filename
      m_dot_filename = trim(output_folder) //  '/' // trim(scheme_name) // '/m_dot.dat'
      open(88, file=m_dot_filename, position='append', status='unknown')
      write(88, '(2(E15.7, 2X))') current_time, m_dot
      close(88)
    end block
  end subroutine calc_m_dot


  subroutine calc_ppi_modes(current_time)
    implicit none
    real*8, intent(in) :: current_time
    integer :: i, j, k, mode_m, base_index, unit_file
    real*8 :: cell_volume, density_cut, azimuthal_mass, mode_phase
    real*8 :: c0, c_real(PPI_MAX_MODE), c_imag(PPI_MAX_MODE)
    real*8 :: values(2 + 4 * PPI_MAX_MODE)
    character(len=256) :: file_name

    if (.not. do_ppi_diagnostics .or. nz < 2) return

    cell_volume = dx * dy * dz
    density_cut = 1.0d3 * rho_floor
    c0 = 0.0d0
    c_real = 0.0d0
    c_imag = 0.0d0

    ! Primero se integra la masa de cada sector azimutal. Esto evita evaluar
    ! funciones trigonométricas dentro del bucle tridimensional completo.
    do k = 1, nz
      azimuthal_mass = 0.0d0
      do j = 1, ny
        do i = 1, nx
          if (p(eq_de, i, j, k) > density_cut) then
            azimuthal_mass = azimuthal_mass + up(eq_de, i, j, k) * cell_volume
          end if
        end do
      end do
      c0 = c0 + azimuthal_mass
      do mode_m = 1, PPI_MAX_MODE
        c_real(mode_m) = c_real(mode_m) + azimuthal_mass * cos(dble(mode_m) * z(k))
        c_imag(mode_m) = c_imag(mode_m) - azimuthal_mass * sin(dble(mode_m) * z(k))
      end do
    end do

    values = 0.0d0
    values(1) = current_time
    values(2) = c0
    do mode_m = 1, PPI_MAX_MODE
      base_index = 3 + 4 * (mode_m - 1)
      values(base_index) = c_real(mode_m)
      values(base_index + 1) = c_imag(mode_m)
      if (c0 > tiny(1.0d0)) then
        values(base_index + 2) = sqrt(c_real(mode_m)**2 + c_imag(mode_m)**2) / c0
        mode_phase = -atan2(c_imag(mode_m), c_real(mode_m)) / dble(mode_m)
        values(base_index + 3) = mode_phase
      end if
    end do

    file_name = trim(output_folder) // '/' // trim(scheme_name) // '/ppi_modes.dat'
    open(newunit=unit_file, file=trim(file_name), status='unknown', &
         position='append', action='write')
    write(unit_file, '(*(ES24.16E3,1X))') values
    close(unit_file)
  end subroutine calc_ppi_modes


  subroutine calc_global_diagnostics(current_time)
    implicit none
    real*8, intent(in) :: current_time
    integer :: i, j, k, ii, jj, unit_file
    integer(kind=8) :: atmosphere_cells, invalid_cells
    real*8 :: cell_volume, density_cut, rho, pres, v_sq
    real*8 :: mass_total, disk_mass, angular_momentum
    real*8 :: rho_minimum, pressure_minimum, rho_maximum, pressure_maximum
    real*8 :: max_velocity_squared
    character(len=256) :: file_name

    if (.not. do_ppi_diagnostics) return

    cell_volume = dx * dy * dz
    density_cut = 1.0d3 * rho_floor
    mass_total = 0.0d0
    disk_mass = 0.0d0
    angular_momentum = 0.0d0
    rho_minimum = huge(1.0d0)
    pressure_minimum = huge(1.0d0)
    rho_maximum = -huge(1.0d0)
    pressure_maximum = -huge(1.0d0)
    max_velocity_squared = 0.0d0
    atmosphere_cells = 0_8
    invalid_cells = 0_8

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho = p(eq_de, i, j, k)
          pres = p(eq_pr, i, j, k)
          v_sq = 0.0d0
          do ii = 1, 3
            do jj = 1, 3
              v_sq = v_sq + gamma_c(ii, jj, i, j, k) * &
                     p(eq_vx + ii - 1, i, j, k) * p(eq_vx + jj - 1, i, j, k)
            end do
          end do

          if (.not. ieee_is_finite(rho) .or. .not. ieee_is_finite(pres) .or. &
              .not. ieee_is_finite(v_sq) .or. rho <= 0.0d0 .or. pres <= 0.0d0 .or. &
              v_sq < 0.0d0 .or. v_sq >= 1.0d0) then
            invalid_cells = invalid_cells + 1_8
          end if

          if (ieee_is_finite(rho)) then
            rho_minimum = min(rho_minimum, rho)
            rho_maximum = max(rho_maximum, rho)
          end if
          if (ieee_is_finite(pres)) then
            pressure_minimum = min(pressure_minimum, pres)
            pressure_maximum = max(pressure_maximum, pres)
          end if
          if (ieee_is_finite(v_sq)) max_velocity_squared = max(max_velocity_squared, v_sq)

          mass_total = mass_total + up(eq_de, i, j, k) * cell_volume
          angular_momentum = angular_momentum + up(eq_vz, i, j, k) * cell_volume
          if (rho > density_cut) then
            disk_mass = disk_mass + up(eq_de, i, j, k) * cell_volume
          else
            atmosphere_cells = atmosphere_cells + 1_8
          end if
        end do
      end do
    end do

    file_name = trim(output_folder) // '/' // trim(scheme_name) // '/global_diagnostics.dat'
    open(newunit=unit_file, file=trim(file_name), status='unknown', &
         position='append', action='write')
    write(unit_file, '(9(ES24.16E3,1X),2(I0,1X))') current_time, mass_total, disk_mass, &
      angular_momentum, rho_minimum, pressure_minimum, rho_maximum, pressure_maximum, &
      max_velocity_squared, atmosphere_cells, invalid_cells
    close(unit_file)
  end subroutine calc_global_diagnostics


  subroutine write_perturbation_event(current_time)
    implicit none
    real*8, intent(in) :: current_time
    integer :: unit_file
    character(len=256) :: file_name

    if (.not. do_ppi_diagnostics) return
    file_name = trim(output_folder) // '/' // trim(scheme_name) // '/perturbation_events.dat'
    open(newunit=unit_file, file=trim(file_name), status='unknown', &
         position='append', action='write')
    write(unit_file, '(ES24.16E3,1X,2(I0,1X),2(ES24.16E3,1X))') current_time, &
      perturbation_type, perturbation_seed, perturbation_amplitude, perturbation_mode
    close(unit_file)
  end subroutine write_perturbation_event


  subroutine calc_convergence_norms()
    implicit none
    integer :: i, i_start, i_end, j_idx, unit_out, idx_max_err
    real*8 :: rho_exact
    real*8 :: err_L1, err_L2, sum_L1, sum_L2
    real*8 :: err_Linf, diff, x_max_err
    real*8 :: sinc_factor, sinc_argument, celdas_validas, advected_coordinate

    real*8, allocatable :: p_conv(:,:,:,:)
    character(len=256) :: file_name

    if (trim(case_name) /= 'Conv' .and. trim(case_name) /= 'MichelA') return

    sum_L1 = 0.0d0
    sum_L2 = 0.0d0
    err_Linf = 0.0d0
    idx_max_err = 1
    j_idx = max(1, ny / 2)
    x_max_err = x(1)

    if (trim(case_name) == 'Conv') then
      ! Solución exacta en promedios celulares. Las fronteras son periódicas.
      sinc_argument = 0.5d0 * advected_wave_number * pi * dx
      if (abs(sinc_argument) > epsilon(1.0d0)) then
        sinc_factor = sin(sinc_argument) / sinc_argument
      else
        sinc_factor = 1.0d0
      end if
      celdas_validas = dble(nx)
      do i = 1, nx
        advected_coordinate = x(i) - advected_wave_speed * integration_time
        rho_exact = advected_wave_density + advected_wave_amplitude * sinc_factor * &
                    sin(advected_wave_number * pi * advected_coordinate)
        diff = abs(p(eq_de, i, j_idx, 1) - rho_exact)
        sum_L1 = sum_L1 + diff
        sum_L2 = sum_L2 + diff**2
        if (diff > err_Linf) then
          err_Linf = diff
          idx_max_err = i
          x_max_err = x(i)
        end if
      end do
      file_name = trim(output_folder) // '/convergence_' // trim(scheme_name) // '.dat'
    else
      ! Michel conserva la comparación histórica en la región interior segura.
      allocate(p_conv(neq, 1:nx, 1:ny, 1:nz))
      call setup_michel_accretion_initial(p_conv, michel_critical_radius, &
                                           michel_critical_density)
      i_start = max(1, int(0.10d0 * dble(nx)))
      i_end = min(nx, max(i_start, int(0.95d0 * dble(nx))))
      celdas_validas = dble(i_end - i_start + 1)
      do i = i_start, i_end
        rho_exact = p_conv(eq_de, i, j_idx, 1)
        diff = abs(p(eq_de, i, j_idx, 1) - rho_exact)
        sum_L1 = sum_L1 + diff
        sum_L2 = sum_L2 + diff**2
        if (diff > err_Linf) then
          err_Linf = diff
          idx_max_err = i
          if (use_log_r) then
            x_max_err = exp(x(i))
          else
            x_max_err = x(i)
          end if
        end if
      end do
      deallocate(p_conv)
      file_name = trim(output_folder) // '/michel_convergence_' // trim(scheme_name) // '.dat'
    end if

    err_L1 = sum_L1 / celdas_validas
    err_L2 = sqrt(sum_L2 / celdas_validas)

    print *, "N =", nx, "| L1 =", err_L1, "| L2 =", err_L2, &
             "| Linf =", err_Linf, " localizado en x =", x_max_err, &
             " indice = ", idx_max_err

    unit_out = 20

    ! ====================================================================
    ! LÓGICA DE ESCRITURA INTELIGENTE
    ! ====================================================================
    if (nx == 32) then
      open(unit=unit_out, file=trim(file_name), status='replace')
      write(unit_out, '(A)') '# N L1 L2 Linf'
    else
      open(unit=unit_out, file=trim(file_name), status='unknown', position='append')
    end if

    write(unit_out, '(I6, 3(2X, E20.12))') nx, err_L1, err_L2, err_Linf
    close(unit_out)

  end subroutine calc_convergence_norms

end module output
