! =======================================================================================
! PROGRAMA PRINCIPAL: GRHD2 (General Relativistic Hydrodynamics 2D/3D)
! ---------------------------------------------------------------------------------------
! Orquestador central de la simulación. Controla el bucle temporal principal, la 
! integración espacial mediante el Método de Líneas (Method of Lines) delegada al 
! módulo 'evolution', la actualización temporal usando Runge-Kutta SSP de 3er orden, 
! y el paralelismo de memoria compartida masivo (OpenMP).
! =======================================================================================
program grhd2
  use variables
  use parameters
  use metrics
  use initialization
  use conditions
  use equations
  use fluxes
  use reconstruction
  use evolution
  use output
  use omp_lib
  implicit none

  real*8 :: start_time, end_time
  integer :: i, j, rk, k, argument_count, io_status
  character(len=1024) :: first_argument, parameter_file
  logical :: check_only

  ! ---  VARIABLES PARA EL CRONÓMETRO Y ETA ---
  real*8 :: session_start_time, session_progress
  real*8 :: current_wtime, elapsed_time, eta
  integer :: current_percent, last_percent
  integer :: eta_h, eta_m, eta_s  
  integer :: elap_h, elap_m, elap_s  

  real*8 :: total_mass

  argument_count = command_argument_count()
  check_only = .false.
  first_argument = ''
  parameter_file = ''

  if (argument_count == 1) then
    call get_command_argument(1, parameter_file)
  else if (argument_count == 2) then
    call get_command_argument(1, first_argument)
    call get_command_argument(2, parameter_file)
    if (trim(first_argument) /= '--check') then
      write(*,*) 'Uso: ./grhd2 [--check] archivo.par'
      error stop 1
    end if
    check_only = .true.
  else
    write(*,*) 'Uso: ./grhd2 [--check] archivo.par'
    error stop 1
  end if

  call read_parameter_header(trim(parameter_file))
  start_time = omp_get_wtime()

  ! =========================================================================
  ! 1. ARRANQUE DEL SISTEMA (Cold Start vs Checkpoint Restart)
  ! =========================================================================
  if (do_restart) then
    print *, '=========================================='
    print *, '>>> REANUDANDO SIMULACION DESDE BINARIO'
    print *, '>>> Archivo: ', trim(restart_file)
    print *, '=========================================='
    
    open(20, file=trim(restart_file), status='old', form='unformatted', iostat=io_status)
    if (io_status /= 0) then
      write(*,*) 'ERROR: no se pudo abrir el checkpoint: ', trim(restart_file)
      error stop 1
    end if
    call read_checkpoint_metadata(20)
    call read_parameter_file(trim(parameter_file), .true.)
    call validate_parameters(trim(parameter_file))
    call print_parameter_summary(trim(parameter_file))
    if (check_only) then
      close(20)
      print *, 'Archivo de parametros y checkpoint validos; no se inicio la simulacion.'
      stop
    end if

    ! Recalculamos parámetros termodinámicos auxiliares
    g1 = adb_idx / (adb_idx - 1.0d0)
    rho_floor = 1.0d-10
    p_floor   = rho_floor * 1.0d-3
    D_floor   = rho_floor
    tau_floor = p_floor * (g1 - 1.0d0)

    call set_metric_type()
    call allocate_and_grid()
    call precalculate_metric_cache()
    call init_wavespeed_solver()
    call read_checkpoint_state(20, n_steps, integration_time)
    close(20)
    
    ! Sincronizamos y recuperamos primitivas para el integrador RK3
    u = up

    ! Propagación de y(j) para la inversión termodinámica 3D
    !$OMP PARALLEL DO PRIVATE(i, j, k)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call cons_to_prim(u(:,i,j,k), p(:,i,j,k), i, j, k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call set_boundary_conditions(p)
    next_checkpoint = (int(integration_time / checkpoint_interval) + 1.0d0) * checkpoint_interval
  else
    ! --- ARRANQUE EN FRÍO NORMAL (T=0) ---
    integration_time = 0.0d0
    call initialize_problem(trim(parameter_file), check_only)
    if (check_only) stop
    next_checkpoint = checkpoint_interval

    ! Llenado inicial del vector de conservativas
    !$OMP PARALLEL DO PRIVATE(i, j, k)
    do k = 1, nz 
      do j = 1, ny
        do i = 1, nx 
          call prim_to_cons(p(:,i,j,k), up(:,i,j,k), i, j, k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    call save_vtk_at_time(integration_time)
  end if

  call initialize_auxiliary_output(.not. do_restart)
  if (do_gw_extraction) call calc_gw_strain(integration_time)
  if (do_mdot_extraction) call calc_m_dot(integration_time)
  if (do_ppi_diagnostics) then
    call calc_ppi_modes(integration_time)
    call calc_global_diagnostics(integration_time)
  end if

  ! -------------------------------------------
  next_save_time =  integration_time + save_interval
  call update_adaptive_dt() ! Inicialización del primer dt

  session_start_time = integration_time
  if (final_time > 0.0d0) then
    last_percent = int((integration_time / final_time) * 100.0d0)
  else
    last_percent = 100
  end if

    total_mass = 0.0d0
    
    ! Integral de volumen lógico 3D Coordenado
    ! (El Jacobiano dr/dx ya está matemáticamente embebido en el sqrt(gamma) dentro de u(eq_de))
    !$OMP PARALLEL DO REDUCTION(+:total_mass) PRIVATE(i, j, k)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          total_mass = total_mass + up(eq_de, i, j, k) * dx * dy * dz
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    
    write(*,*) 'Initial Mass in domain: ', total_mass

  ! =========================================================================
  ! 2. BUCLE PRINCIPAL DE INTEGRACIÓN TEMPORAL (SSP-RK3)
  ! =========================================================================
  do while (integration_time < final_time)
    call update_adaptive_dt()

    if (integration_time + dt > final_time) then
      dt = final_time - integration_time
    end if

    ! --- INYECCIÓN RETRASADA, CONFIGURABLE Y REPRODUCIBLE ---
    if (apply_perturbation .and. .not. perturbation_applied .and. &
        integration_time >= perturbation_time) then
      print *, ""
      print *, "=========================================================="
      print *, ">>> TIEMPO DE PERTURBACION ALCANZADO: ", integration_time
      print *, "=========================================================="

      select case(perturbation_type)
      case(PERT_PRESSURE_NOISE)
        call inject_pressure_noise(p, perturbation_amplitude, perturbation_seed)
      case(PERT_DENSITY_NOISE)
        call inject_density_noise(p, perturbation_amplitude, perturbation_seed)
      case(PERT_DENSITY_MODE)
        call inject_density_mode(p, perturbation_mode, perturbation_amplitude)
      case default
        write(*,*) 'Perturbation disabled: type = PERT_NONE.'
        apply_perturbation = .false.
      end select

      if (perturbation_type /= PERT_NONE) then
        ! Re-sincronizamos U y UP tras modificar P.
        !$OMP PARALLEL DO PRIVATE(i, j, k)
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              call prim_to_cons(p(:,i,j,k), u(:,i,j,k), i, j, k)
              up(:,i,j,k) = u(:,i,j,k)
            end do
          end do
        end do
        !$OMP END PARALLEL DO
        call set_boundary_conditions(p)

        perturbation_applied = .true.
        call write_perturbation_event(integration_time)
        if (do_ppi_diagnostics) then
          call calc_ppi_modes(integration_time)
          call calc_global_diagnostics(integration_time)
        end if
      end if
    end if

    n_steps = n_steps + 1
    integration_time = integration_time + dt

    u = up 
    call set_boundary_conditions(p)

    ! Integrador Runge-Kutta Fuerte de Preservación de Estabilidad (SSP)
    do rk = 1, 3
    
      ! A. Integración Espacial: Reconstrucción, Flujos y Términos Fuente
      call calc_rhs(p, rhs)

      ! B. Avance Temporal Low-Storage
      ! (Las operaciones matriciales vectorizadas funcionan sin bucles for)
      if (rk == 1) then
        up = u + rhs * dt
      else if (rk == 2) then
        up = 0.75d0 * u + 0.25d0 * (up + rhs * dt)
      else
        up = (u + 2.0d0*(up + rhs * dt)) / 3.0d0
      end if

      ! C. Recuperación Termodinámica (Conservativas -> Primitivas)
      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = 1, nz 
        do j = 1, ny
          do i = 1, nx
            call cons_to_prim(up(:,i,j,k), p(:,i,j,k), i, j, k)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! D. Actualización de Topología de Red
      call set_boundary_conditions(p)

    end do

    ! =====================================================================
    ! 3. DIAGNÓSTICOS Y SALIDA DE DATOS
    ! =====================================================================
    ! Corrección: Sustitución de 'is_flat' por la validación arquitectónica global
    if (trim(metric_type) /= 'Minkowski' .and. bh_mass > 0.0d0) then
      if (mod(n_steps, 10) == 0) then
        if (do_gw_extraction) call calc_gw_strain(integration_time)
        if (do_mdot_extraction) call calc_m_dot(integration_time)
      end if
    end if
    if (do_ppi_diagnostics .and. mod(n_steps, diagnostic_stride) == 0) then
      call calc_ppi_modes(integration_time)
      call calc_global_diagnostics(integration_time)
    end if

    ! --- CONSERVACIÓN DE MASA GLOBAL ---
    if (mod(n_steps, 1000) == 0) then
      total_mass = 0.0d0
      
      ! Integral de volumen lógico 3D Coordenado
      ! (El Jacobiano dr/dx ya está matemáticamente embebido en el sqrt(gamma) dentro de u(eq_de))
      !$OMP PARALLEL DO REDUCTION(+:total_mass) PRIVATE(i, j, k)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            total_mass = total_mass + u(eq_de, i, j, k) * dx * dy * dz
          end do
        end do
      end do
      !$OMP END PARALLEL DO
      
      write(*,*) 'Total Mass in domain: ', total_mass
    end if

    ! --- CRONÓMETRO INTELIGENTE (ETA) ---
    current_percent = int((integration_time / final_time) * 100.0d0)
    if (current_percent > last_percent) then
      current_wtime = omp_get_wtime()
      elapsed_time = current_wtime - start_time
      
      session_progress = (integration_time - session_start_time) / (final_time - session_start_time)
      
      if (session_progress > 0.0d0) then
        eta = (elapsed_time / session_progress) - elapsed_time
        
        eta_h = int(eta) / 3600
        eta_m = mod(int(eta), 3600) / 60
        eta_s = mod(int(eta), 60)
        
        elap_h = int(elapsed_time) / 3600
        elap_m = mod(int(elapsed_time), 3600) / 60
        elap_s = mod(int(elapsed_time), 60)
        
        write(*, '(A, I3, A, I3, A, I2.2, A, I2.2, A, I3, A, I2.2, A, I2.2, A, A, I8, A, F10.3)') &
          '>>> PROGRESO: ', current_percent, '% | Transcurrido: ', &
          elap_h, 'h ', elap_m, 'm ', elap_s, 's | Faltan aprox: ', &
          eta_h, 'h ', eta_m, 'm ', eta_s, 's  <<< | ', &
          'Paso ', n_steps, ' | Tiempo: ', integration_time
      end if
      
      last_percent = current_percent
    end if

    ! --- EJECUCIÓN DEL GUARDADO A DISCO ---
    if (integration_time >= next_save_time) then
      call save_vtk_at_time(next_save_time)
      next_save_time = next_save_time + save_interval
    end if

    if (integration_time >= next_checkpoint) then
      call save_checkpoint(n_steps, integration_time)
      next_checkpoint = next_checkpoint + checkpoint_interval
    end if
    
  end do

  ! Evaluación final de errores L1/L2 si corresponde al test
  call calc_convergence_norms() 

  ! =========================================================================
  ! 4. FINALIZACIÓN Y LIMPIEZA
  ! =========================================================================
  end_time = omp_get_wtime()
  elapsed_time = end_time - start_time
  
  elap_h = int(elapsed_time) / 3600
  elap_m = mod(int(elapsed_time), 3600) / 60
  elap_s = mod(int(elapsed_time), 60)

  write(*, '(/, A)') '========================================================='
  write(*, '(A, F10.3)') ' SIMULACIÓN FINALIZADA EXITOSAMENTE EN T = ', integration_time
  write(*, '(A, I8)')    ' TOTAL DE PASOS TEMPORALES: ', n_steps
  write(*, '(A, I3, A, I2.2, A, I2.2, A)') &
    ' TIEMPO TOTAL DE EJECUCIÓN: ', elap_h, 'h ', elap_m, 'm ', elap_s, 's'
  write(*, '(A, /)') '========================================================='

  deallocate(x,y,z,var_names)
  deallocate(u,up,s,rhs,p)
  deallocate(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c,sqrt_gamma_c,chris_c,dg_c,dlna_c)

  if (nx > 1) then
    deallocate(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
  end if

  if (ny > 1) then
    deallocate(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
  end if

  if (nz > 1) then
    deallocate(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
  end if

end program grhd2
