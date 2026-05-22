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
  integer :: i, j, rk, k

  ! ---  VARIABLES PARA EL CRONÓMETRO Y ETA ---
  real*8 :: session_start_time, session_progress
  real*8 :: current_wtime, elapsed_time, eta
  integer :: current_percent, last_percent
  integer :: eta_h, eta_m, eta_s  
  integer :: elap_h, elap_m, elap_s  

  real*8 :: total_mass 
  logical :: apply_perturbation = .false. 

  start_time = omp_get_wtime()

  ! =========================================================================
  ! 1. ARRANQUE DEL SISTEMA (Cold Start vs Checkpoint Restart)
  ! =========================================================================
  if (do_restart) then
    print *, '=========================================='
    print *, '>>> REANUDANDO SIMULACION DESDE BINARIO'
    print *, '>>> Archivo: ', trim(restart_file)
    print *, '=========================================='
    
    open(20, file=trim(restart_file), status='old', form='unformatted')
    call read_checkpoint_metadata(20)

    ! Recalculamos parámetros termodinámicos auxiliares
    g1 = adb_idx / (adb_idx - 1.0d0)
    rho_floor = 1.0d-10
    p_floor   = rho_floor * 1.0d-3
    D_floor   = rho_floor
    tau_floor = p_floor / (g1 - 1.0d0)

    call set_metric_type()
    call allocate_and_grid()
    call precalculate_metric_cache()
    call init_wavespeed_solver()
    call read_checkpoint_state(20, n_steps, integration_time)
    close(20)
    
    ! final_time = 1000.0d0 ! Activar si deseas extender el tiempo final al reanudar
    
    ! Si reanudamos después de T=1000, desactivamos el gatillo para
    ! no volver a inyectar perturbaciones.
    if (integration_time >= 1000.0d0) apply_perturbation = .false.

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
    next_checkpoint = checkpoint_interval
    call initialize_problem()

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
    if (do_gw_extraction) call calc_gw_strain(integration_time)
    if (do_mdot_extraction) call calc_m_dot(integration_time)
  end if

  ! -------------------------------------------
  next_save_time =  integration_time + save_interval
  call update_adaptive_dt() ! Inicialización del primer dt

  session_start_time = integration_time
  last_percent = int((integration_time / final_time) * 100.0d0)

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

    ! --- INYECCIÓN RETRASADA DE PERTURBACIÓN ---
    if (integration_time >= 1000.0d0 .and. apply_perturbation) then
      print *, ""
      print *, "=========================================================="
      print *, ">>> TIEMPO T=1000 ALCANZADO."
      print *, ">>> Inyectando superposición de modos (m=4) + Ruido Blanco"
      print *, "=========================================================="
      
      call save_vtk_at_time(integration_time) ! Guardar estado inmaculado justo antes
      ! call inject_density_mode(p, 2.0d0)      ! Inyectar modo
      call inject_pressure_noise(p)
      ! call inject_density_noise(p)

      ! Re-sincronizamos U y UP tras modificar P para evitar inestabilidades numéricas
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
      
      apply_perturbation = .false. ! Desarmamos el gatillo
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

end program grhd2