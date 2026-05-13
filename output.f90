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
  implicit none
  public :: save_vtk_at_time, save_checkpoint, read_checkpoint_metadata, &
            read_checkpoint_state, calc_gw_strain, calc_m_dot, calc_convergence_norms

contains

  ! ===================================================================================
  ! 1. EXPORTADOR VTK (VisIt/ParaView) - FORMATO BINARIO ESTRUCTURADO
  ! ===================================================================================
  subroutine save_vtk_data(file_name)
    character(len=*), intent(in) :: file_name
    integer :: eq, i, j, k, j_src, k_src, unit_file, io_status
    integer :: ny_vtk, nz_vtk
    character(len=150) :: full_path
    character(len=100) :: header_line
    real*8 :: th_val  ! Variable para forzar el cierre topológico

    real*8, allocatable :: pts_buffer(:,:,:,:)
    real*8, allocatable :: var_buffer(:,:,:)
    real*8 :: r_phys

    ! -------------------------------------------------------------
    ! LÓGICA DE CIERRES TOPOLÓGICOS (Cierre Azimutal y Polar)
    ! -------------------------------------------------------------
    ny_vtk = ny
    nz_vtk = nz

    select case(trim(geom_type))
      case('Spherical')
        ! Cierre Azimutal (phi): Cierra el "hueco de pizza"
        if (nz > 1) nz_vtk = nz + 1
        ! Cierre Polar (theta): Cierra el "eje blanco" en los polos
        if (ny > 1) ny_vtk = ny + 2
      case('Cylindrical')
        ! Cierre Azimutal cilíndrico (aquí phi es el eje y)
        if (ny > 1) ny_vtk = ny + 1
    end select

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

    write(header_line, '(a, i0, 1x, i0, 1x, i0)') 'DIMENSIONS ', nx, ny_vtk, nz_vtk
    write(unit_file) trim(header_line) // char(10)

    write(header_line, '(a, i0, a)') 'POINTS ', nx * ny_vtk * nz_vtk, ' double'
    write(unit_file) trim(header_line) // char(10)

    allocate(pts_buffer(3, 1:nx, 1:ny_vtk, 1:nz_vtk))
    
    ! -------------------------------------------------------------
    ! BUCLE DE MAPEO GEOMÉTRICO 3D (X, Y, Z)
    ! -------------------------------------------------------------
    do k = 1, nz_vtk
      k_src = k
      if (k == nz_vtk .and. nz_vtk > nz) k_src = 1
      
      do j = 1, ny_vtk
        
        ! Inyección de los polos falsos para el renderizado continuo
        if (trim(geom_type) == 'Spherical' .and. ny > 1) then
          if (j == 1) then
            th_val = 0.0d0           ! Forzamos Polo Norte exacto
            j_src = 1                ! Copiamos gas de la primera celda
          else if (j == ny_vtk) then
            th_val = pi              ! Forzamos Polo Sur exacto
            j_src = ny               ! Copiamos gas de la última celda
          else
            th_val = y(j - 1)        ! Coordenadas reales desplazadas por el polo norte
            j_src = j - 1
          end if
        else
          ! Lógica Cilíndrica o Cartesiana
          j_src = j
          if (j == ny_vtk .and. ny_vtk > ny) j_src = 1
          th_val = y(j_src)
        end if
        
        do i = 1, nx
          ! Corrección: Actualización a use_log_r
          if (use_log_r) then
            r_phys = exp(x(i)) 
          else
            r_phys = x(i) 
          end if

          select case(trim(geom_type))
            case('Spherical')
              ! Mapeo Esférico usando th_val en lugar de y(j)
              pts_buffer(1, i, j, k) = r_phys * sin(th_val) * cos(z(k_src))
              pts_buffer(2, i, j, k) = r_phys * sin(th_val) * sin(z(k_src))
              pts_buffer(3, i, j, k) = r_phys * cos(th_val)
              
            case('Cylindrical')
              pts_buffer(1, i, j, k) = r_phys * cos(th_val)
              pts_buffer(2, i, j, k) = r_phys * sin(th_val)
              pts_buffer(3, i, j, k) = z(k_src)
              
            case default
              pts_buffer(1, i, j, k) = r_phys
              pts_buffer(2, i, j, k) = th_val
              pts_buffer(3, i, j, k) = z(k_src)
          end select
        end do
      end do
    end do
    
    write(unit_file) pts_buffer
    write(unit_file) char(10) 
    deallocate(pts_buffer)

    write(header_line, '(a, i0)') 'POINT_DATA ', nx * ny_vtk * nz_vtk
    write(unit_file) trim(header_line) // char(10)

    allocate(var_buffer(1:nx, 1:ny_vtk, 1:nz_vtk))
    
    ! -------------------------------------------------------------
    ! EXPORTACIÓN DE VARIABLES PRIMITIVAS
    ! -------------------------------------------------------------
    do eq = 1, neq
      write(header_line, '(a, 1x, a, a)') 'SCALARS', trim(var_names(eq)), ' double 1'
      write(unit_file) trim(header_line) // char(10)
      write(unit_file) 'LOOKUP_TABLE default' // char(10)
      
      do k = 1, nz_vtk
        k_src = k
        if (k == nz_vtk .and. nz_vtk > nz) k_src = 1
        
        do j = 1, ny_vtk
          ! Replicamos la misma lógica de j_src para la extracción de datos
          if (trim(geom_type) == 'Spherical' .and. ny > 1) then
            if (j == 1) then
              j_src = 1
            else if (j == ny_vtk) then
              j_src = ny
            else
              j_src = j - 1
            end if
          else
            j_src = j
            if (j == ny_vtk .and. ny_vtk > ny) j_src = 1
          end if
          
          do i = 1, nx
            ! Extracción del dato crudo lógico
            var_buffer(i, j, k) = p(eq, i, j_src, k_src)
            
            ! SANEAMIENTO VISUAL: Conversión de velocidad lógica a física (v^r = v^x * r)
            ! Corrección: Actualización a use_log_r
            if (use_log_r .and. eq == eq_vx) then
               var_buffer(i, j, k) = var_buffer(i, j, k) * exp(x(i))
            end if
          end do
        end do
      end do
      
      write(unit_file) var_buffer
      write(unit_file) char(10)
    end do

    ! -------------------------------------------------------------
    ! EXPORTACIÓN DE MÉTRICA PARA DIAGNÓSTICO
    ! -------------------------------------------------------------
    
    ! 1. Lapso Temporal (Alpha)
    write(header_line, '(a)') 'SCALARS Lapse_alpha double 1'
    write(unit_file) trim(header_line) // char(10)
    write(unit_file) 'LOOKUP_TABLE default' // char(10)

    do k = 1, nz_vtk
      k_src = k
      if (k == nz_vtk .and. nz_vtk > nz) k_src = 1
      do j = 1, ny_vtk
        if (trim(geom_type) == 'Spherical' .and. ny > 1) then
          if (j == 1) then
            j_src = 1
          else if (j == ny_vtk) then
            j_src = ny
          else
            j_src = j - 1
          end if
        else
          j_src = j
          if (j == ny_vtk .and. ny_vtk > ny) j_src = 1
        end if
        do i = 1, nx
          var_buffer(i, j, k) = alpha_c(i, j_src, k_src)
        end do
      end do
    end do
    write(unit_file) var_buffer
    write(unit_file) char(10)

    ! 2. Shift Radial (Beta^r)
    write(header_line, '(a)') 'SCALARS Shift_beta_r double 1'
    write(unit_file) trim(header_line) // char(10)
    write(unit_file) 'LOOKUP_TABLE default' // char(10)

    do k = 1, nz_vtk
      k_src = k
      if (k == nz_vtk .and. nz_vtk > nz) k_src = 1
      do j = 1, ny_vtk
        if (trim(geom_type) == 'Spherical' .and. ny > 1) then
          if (j == 1) then
            j_src = 1
          else if (j == ny_vtk) then
            j_src = ny
          else
            j_src = j - 1
          end if
        else
          j_src = j
          if (j == ny_vtk .and. ny_vtk > ny) j_src = 1
        end if
        do i = 1, nx
          var_buffer(i, j, k) = beta_c(1, i, j_src, k_src)
        end do
      end do
    end do
    write(unit_file) var_buffer
    write(unit_file) char(10)
    
    deallocate(var_buffer)
    close(unit_file)
  end subroutine save_vtk_data


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

    write(file_name, '(A, "/checkpoint_", A, "_", I5.5, ".rst")') trim(output_folder), trim(scheme_name), int(current_time)
    print *, '>>> Guardando Checkpoint de Respaldo: ', trim(file_name)

    unit_file = 20
    open(unit_file, file=file_name, status='replace', action='write', form='unformatted', iostat=io_status)
    if (io_status == 0) then
      ! Topología de malla 3D y profundidad de fantasmas
      write(unit_file) nghost, nx, ny, nz 
      
      ! Fronteras lógicas (Actualizado a r_min/r_max)
      write(unit_file) r_min, r_max, y_min, y_max, z_min, z_max
      
      ! Física y Banderas Arquitectónicas (Actualizado con banderas de extracción)
      write(unit_file) bh_mass, adb_idx
      write(unit_file) use_log_r, use_shock_sensor, do_gw_extraction, do_mdot_extraction
      write(unit_file) metric_type, geom_type
      
      ! Control de tiempo y salida
      write(unit_file) final_time, CFL, save_interval
      write(unit_file) output_prefix, output_folder
      
      ! Identificadores de esquemas
      write(unit_file) case_id, case_name, scheme_name
      write(unit_file) rec_method_id, tvd_limiter_id, riemann_solver_id
      
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

    ! Topología de malla
    read(unit_file) nghost, nx, ny, nz
    
    ! Fronteras lógicas
    read(unit_file) r_min, r_max, y_min, y_max, z_min, z_max
    
    ! Física y Banderas Arquitectónicas (Actualizado con banderas de extracción)
    read(unit_file) bh_mass, adb_idx
    read(unit_file) use_log_r, use_shock_sensor, do_gw_extraction, do_mdot_extraction
    read(unit_file) metric_type, geom_type
    
    ! Control
    read(unit_file) final_time, CFL, save_interval
    read(unit_file) output_prefix, output_folder
    
    ! Eschemas
    read(unit_file) case_id, case_name, scheme_name
    read(unit_file) rec_method_id, tvd_limiter_id, riemann_solver_id
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

  subroutine calc_gw_strain(current_time)
    real*8, intent(in) :: current_time
    integer :: i, j, k
    real*8 :: Ixx_ddot, Iyy_ddot, Ixy_ddot
    real*8 :: r, th, phi, vx, vy, x_coord, y_coord, sin_th
    real*8 :: rho, pres, vr, vth, vphi
    real*8 :: dV, mass_r3
    real*8 :: h_plus, h_cross
    real*8 :: observer_dist = 1000.0d0 

    ! Corrección: Actualización de bandera
    if (trim(metric_type) == 'Minkowski' .or. bh_mass < 0.1d0) return 

    Ixx_ddot = 0.0d0
    Iyy_ddot = 0.0d0
    Ixy_ddot = 0.0d0

    ! Se propaga j a todo el bloque y se calcula el volumen dV en coordenadas esféricas
    !$OMP PARALLEL DO PRIVATE(i, j, k, r, th, phi, sin_th, vr, vth, vphi, rho, pres, vx, vy, x_coord, y_coord, mass_r3, dV) &
    !$OMP REDUCTION(+:Ixx_ddot, Iyy_ddot, Ixy_ddot)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          
          ! Corrección vital para integrales físicas en mallas logarítmicas
          if (use_log_r) then
            r = exp(x(i))
          else
            r = x(i)
          end if
          
          th = y(j)
          phi = z(k)
          sin_th = sin(th)
          
          rho  = p(eq_de, i, j, k)
          pres = p(eq_pr, i, j, k)
          vr   = p(eq_vx, i, j, k)
          vth  = p(eq_vy, i, j, k)
          vphi = p(eq_vz, i, j, k) 

          ! Mapeo Esférico a Cartesiano
          x_coord = r * sin_th * cos(phi)
          y_coord = r * sin_th * sin(phi)
          
          ! Transformación tensorial básica de la velocidad a componentes cartesianos X e Y
          vx = vr * sin_th * cos(phi) + r * vth * cos(th) * cos(phi) - r * sin_th * vphi * sin(phi)
          vy = vr * sin_th * sin(phi) + r * vth * cos(th) * sin(phi) + r * sin_th * vphi * cos(phi)

          ! Elemento de volumen invariante adaptado a la topología
          if (use_log_r) then
            dV = (r**3) * sin_th * dx * dy * dz ! dr = r * dx
          else
            dV = (r**2) * sin_th * dx * dy * dz ! dr = dx
          end if
          
          mass_r3 = bh_mass / (r**3)

          ! Cuadripolo proyectado en el plano ecuatorial
          Ixx_ddot = Ixx_ddot + (2.0d0 * rho * vx * vx + 2.0d0 * pres - 2.0d0 * rho * mass_r3 * x_coord * x_coord) * dV
          Iyy_ddot = Iyy_ddot + (2.0d0 * rho * vy * vy + 2.0d0 * pres - 2.0d0 * rho * mass_r3 * y_coord * y_coord) * dV
          Ixy_ddot = Ixy_ddot + (2.0d0 * rho * vx * vy - 2.0d0 * rho * mass_r3 * x_coord * y_coord) * dV
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    h_plus  = (Ixx_ddot - Iyy_ddot) / observer_dist
    h_cross = (2.0d0 * Ixy_ddot) / observer_dist

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
    real*8 :: alpha, det, beta(3), g(3), v_eff
    integer :: j, k

    ! Corrección: Actualización de bandera
    if (trim(metric_type) == 'Minkowski' .or. bh_mass < 0.1d0) return 

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


  subroutine calc_convergence_norms()
    implicit none
    integer :: i, j_idx, unit_out, idx_max_err
    real*8 :: rho_exact
    real*8 :: err_L1, err_L2, sum_L1, sum_L2
    real*8 :: err_Linf, diff, x_max_err
    real*8 :: sinc_factor, celdas_validas
    
    ! Corrección: Sincronización de variables físicas con initial_conditions
    real*8 :: r_critical = 100.0d0
    real*8 :: rho_critical = 1.0d-1
    
    real*8, allocatable :: p_conv(:,:,:,:)
    character(len=100) :: file_name

    if (trim(case_name) /= 'Conv' .and. trim(case_name) /= 'MichelA') return

    sum_L1 = 0.0d0
    sum_L2 = 0.0d0
    err_Linf = 0.0d0
    sinc_factor = sin(pi * dx) / (pi * dx)

    allocate(p_conv(neq, 1:nx, 1:ny, 1:nz))
    ! Ojo: si en este setup se asume que MichelAccretion es unidimensional, la 
    ! evaluación de L1 asume simetría en el plano polar
    call setup_michel_accretion_initial(p_conv, r_critical, rho_critical)

    celdas_validas = dble(int(nx * 0.95) - int(nx * 0.10) + 1)
    
    ! ================================================================
    ! CORRECCIÓN: Blindaje contra división entera en 1D (ny=1 -> j=0)
    ! ================================================================
    j_idx = max(1, ny / 2)

    do i = int(nx * 0.10), int(nx * 0.95)
      ! Usamos los datos iniciales esféricos en el centro del eje polar seguro
      rho_exact = p_conv(eq_de, i, j_idx, 1)

      diff = abs(p(eq_de, i, j_idx, 1) - rho_exact)

      sum_L1 = sum_L1 + diff
      sum_L2 = sum_L2 + diff**2
    end do
    
    ! ================================================================
    ! CORRECCIÓN: Evitar sobre-división. Calculamos el promedio exacto
    ! ================================================================
    err_L1 = sum_L1 / celdas_validas
    err_L2 = sqrt(sum_L2 / celdas_validas)

    do i = int(nx * 0.10), int(nx * 0.95)
      rho_exact = p_conv(eq_de, i, j_idx, 1)
      diff = abs(p(eq_de, i, j_idx, 1) - rho_exact)
      
      if (diff > err_Linf) then
        err_Linf = diff
        idx_max_err = i
        
        ! Guardamos el radio físico real
        if (use_log_r) then
          x_max_err = exp(x(i))
        else
          x_max_err = x(i) 
        end if
      end if
    end do
    
    deallocate(p_conv)

    print *, "N =", nx, "| L1 =", err_L1, "| L2 =", err_L2, "| Linf =", err_Linf, " localizado en r =", x_max_err, " indice = ", idx_max_err

    file_name = 'MAconvergencia_' // trim(scheme_name) // '.dat'
    unit_out = 20

    ! ====================================================================
    ! LÓGICA DE ESCRITURA INTELIGENTE
    ! ====================================================================
    if (nx == 32) then
      ! Si es la primera malla del script, destruye el archivo de simulaciones 
      ! pasadas y lo crea en blanco.
      open(unit=unit_out, file=trim(file_name), status='replace')
    else
      ! Si N > 32, añade la fila de datos para armar la curva sin borrar lo anterior.
      open(unit=unit_out, file=trim(file_name), status='unknown', position='append')
    end if

    write(unit_out, '(I6, 3(2X, E20.12))') nx, err_L1, err_L2, err_Linf
    close(unit_out)

  end subroutine calc_convergence_norms

end module output