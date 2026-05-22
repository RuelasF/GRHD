! =======================================================================================
! Módulo: initialization (Configuración de la Malla y Parámetros del Problema)
! ---------------------------------------------------------------------------------------
! Este módulo orquesta el arranque de la simulación. Se encarga de alojar dinámicamente 
! la memoria RAM según la resolución deseada, construir la malla espacial (Grid) asegurando
! una topología Cell-Centered, y establecer el entorno físico base.
! =======================================================================================
module initialization
  use variables
  use metrics
  use conditions
  use fluxes
  implicit none
  private
  public :: initialize_problem, choose_numerical_architecture, allocate_and_grid, precalculate_metric_cache

contains

  ! Subrutina: choose_numerical_architecture
  ! Menú interactivo para definir la arquitectura espacial completa.
  ! Asigna los IDs numéricos del Reconstructor y del Solucionador de Riemann.
  subroutine choose_numerical_architecture()
    integer :: choice_rec, choice_limiter, choice_solver

    ! --- 1. SELECCIÓN DEL RECONSTRUCTOR ---
    print *, "=========================================="
    print *, "Choose reconstruction method:"
    print *, "1. Godunov (1st order) - Alta disipacion (Solo depuracion)"
    print *, "2. TVD    (2nd order) - Robustez estandar balanceada"
    print *, "3. WENO3  (3rd order) - Precision intermedia"
    print *, "4. MP5    (5th order) - Resolucion extrema de choques (Sharp)"
    print *, "5. WENO5  (5th order) - Maxima fidelidad en turbulencia (Smooth)"
    
    ! RECORDATORIO: Descomentar para producción interactiva
    ! read *, choice_rec
    choice_rec = 5

    select case(choice_rec)
    case(1)
      rec_method_id = REC_GODUNOV
      scheme_name = 'godunov'
      nghost = 1  
    case(2)
      rec_method_id = REC_TVD
      nghost = 2  
      
      print *, "Choose TVD limiter:"
      print *, "1. Minmod (Disipativo, estable)"
      print *, "2. Superbee (Compresivo, choques afilados)"
      print *, "3. MC (Estandar TVD, balanceado)"
      read *, choice_limiter

      select case(choice_limiter)
      case(1)
        tvd_limiter_id = LIM_MINMOD
        scheme_name = 'minmod'
      case(2)
        tvd_limiter_id = LIM_SUPERBEE
        scheme_name = 'superbee'
      case(3)
        tvd_limiter_id = LIM_MC
        scheme_name = 'mc'
      case default
        tvd_limiter_id = LIM_MC
        scheme_name = 'mc'
      end select
    case(3)
      rec_method_id = REC_WENO3
      scheme_name = 'weno3'
      nghost = 2
    case(4)
      rec_method_id = REC_MP5
      scheme_name = 'mp5'
      nghost = 3  
    case(5)
      rec_method_id = REC_WENO5
      scheme_name = 'weno5'
      nghost = 3  
    case default
      print *, "Invalid choice. Using Default: MP5."
      rec_method_id = REC_MP5
      scheme_name = 'mp5'
      nghost = 3
    end select

    ! --- 2. SELECCIÓN DEL SOLUCIONADOR DE RIEMANN ---
    print *, "=========================================="
    print *, "Choose Riemann Solver:"
    print *, "1. HLLE (Ultra-robusto, disipativo. Ideal para Agujeros Negros)"
    print *, "2. HLLC (Baja disipacion, resuelve contacto. Ideal para Jets/KHI)"
    
    ! RECORDATORIO: Descomentar para producción interactiva
    ! read *, choice_solver
    choice_solver = 1

    select case(choice_solver)
    case(1)
      riemann_solver_id = RS_HLLE
      solver_name = 'hlle'
    case(2)
      if (is_mhd) then
        riemann_solver_id = RS_HLLD
        solver_name = 'hlld'
        print *, ">> INFO: MHD activado. Cambiando HLLC por HLLD automáticamente."
      else
        riemann_solver_id = RS_HLLC
        solver_name = 'hllc'
      end if
    case default
      print *, "Invalid choice. Using Default: HLLE."
      riemann_solver_id = RS_HLLE
      solver_name = 'hlle'
    end select

    ! Unimos los nombres (Ej. 'weno5_hllc' o 'mp5_hlle')
    scheme_name = trim(scheme_name) // '_' // trim(solver_name)

    print *, "=========================================="
    print *, "Selected Architecture: ", trim(scheme_name)
    print *, "Ghost cells depth set to: ", nghost
  end subroutine choose_numerical_architecture


  ! Subrutina: allocate_and_grid
  ! Genera la malla computacional utilizando topología 'Cell-Centered'.
  subroutine allocate_and_grid()
    integer :: i, j, k
    real*8 :: min_dist, r_horizon_coord

    ! Límites lógicos (computacionales) dependientes de la malla
    if (use_log_r) then
      x_min = log(r_min)
      x_max = log(r_max)
    else
      x_min = r_min
      x_max = r_max
    end if
    
    ! Prevención de división por cero para dominios degradados (ej. casos 1D o 2D)
    dx = (x_max - x_min) / max(1, nx)
    dy = (y_max - y_min) / max(1, ny)
    dz = (z_max - z_min) / max(1, nz)

    ! Asignación dinámica de memoria geométrica tridimensional
    allocate(x(-nghost:nx+nghost), x_face(0:nx))
    allocate(y(-nghost:ny+nghost), y_face(0:ny))
    allocate(z(-nghost:nz+nghost), z_face(0:nz))

    ! -------------------------------------------------------------
    ! Generación de Malla en X (Radial)
    ! -------------------------------------------------------------
    do i = 0, nx
      x_face(i) = x_min + dble(i) * dx
    end do
    do i = -nghost, nx+nghost
      x(i) = x_min + (dble(i) - 0.5d0) * dx
    end do

    ! -------------------------------------------------------------
    ! Localizador de sondas para Agujeros Negros
    ! -------------------------------------------------------------
    if (trim(metric_type) /= 'Minkowski' .and. bh_mass > 0.0d0) then
      ! Adaptamos la coordenada de búsqueda a la topología elegida
      if (use_log_r) then
         r_horizon_coord = log(2.0d0 * bh_mass)
      else
         r_horizon_coord = 2.0d0 * bh_mass
      end if
      
      min_dist = 1000.0d0
      do i = 1, nx
          if (abs(x(i) - r_horizon_coord) < min_dist) then
              min_dist = abs(x(i) - r_horizon_coord)
              idx_horizon = i
          end if
      end do
      
      if (x(idx_horizon) < r_horizon_coord) then
        idx_probe = idx_horizon + 1 
      else
        idx_probe = idx_horizon
      end if
      
      ! Imprimimos el valor físico real recuperado (exp(x) o x)
      if (use_log_r) then
         print *, "Event Horizon mapped to cell ", idx_horizon, "(r = ", exp(x(idx_horizon)), ")"
         print *, "Accretion Probe mapped to cell ", idx_probe, "(r = ", exp(x(idx_probe)), ")"
      else
         print *, "Event Horizon mapped to cell ", idx_horizon, "(r = ", x(idx_horizon), ")"
         print *, "Accretion Probe mapped to cell ", idx_probe, "(r = ", x(idx_probe), ")"
      end if
    else
      idx_horizon = 0
      idx_probe = 0
    end if
    
    ! -------------------------------------------------------------
    ! Generación de Malla en Y (Polar / Sagital)
    ! -------------------------------------------------------------
    do j = 0, ny
      y_face(j) = y_min + dble(j) * dy
    end do
    do j = -nghost, ny+nghost
      y(j) = y_min + (dble(j) - 0.5d0) * dy
    end do

    ! -------------------------------------------------------------
    ! Generación de Malla en Z (Azimutal / Ecuatorial)
    ! -------------------------------------------------------------
    do k = 0, nz
      z_face(k) = z_min + dble(k) * dz
    end do
    do k = -nghost, nz+nghost
      z(k) = z_min + (dble(k) - 0.5d0) * dz
    end do

    ! -------------------------------------------------------------
    ! ALOJAMIENTO MASIVO DE RAM (4D: [Ecuaciones, X, Y, Z])
    ! -------------------------------------------------------------
    allocate(u(neq, 1:nx, 1:ny, 1:nz))
    allocate(up(neq, 1:nx, 1:ny, 1:nz))
    allocate(s(neq, 1:nx, 1:ny, 1:nz))
    allocate(rhs(neq, 1:nx, 1:ny, 1:nz))
    allocate(p(neq, -nghost:nx+nghost, -nghost:ny+nghost, -nghost:nz+nghost))
    
    ! Etiquetas de exportación para software de visualización (VisIt, Paraview)
    if (.not. allocated(var_names)) allocate(var_names(neq))
    var_names(eq_de) = 'Density'
    var_names(eq_pr) = 'Pressure'
    var_names(eq_vx) = 'Velocity_X'
    var_names(eq_vy) = 'Velocity_Y'
    var_names(eq_vz) = 'Velocity_Z'

    ! ==========================================================
    ! ALOCACIÓN DE LA CACHÉ DE MÉTRICA
    ! ==========================================================
    ! 1. Centros (1:nx, 1:ny, 1:nz)
    allocate(alpha_c(1:nx, 1:ny, 1:nz))
    allocate(beta_c(3, 1:nx, 1:ny, 1:nz))
    allocate(gamma_c(3, 1:nx, 1:ny, 1:nz))
    allocate(gmunu_c(0:3, 0:3, 1:nx, 1:ny, 1:nz))
    allocate(sqrt_gamma_c(1:nx, 1:ny, 1:nz))
    allocate(chris_c(0:3, 0:3, 0:3, 1:nx, 1:ny, 1:nz))
    allocate(dg_c(0:3, 0:3, 1:3, 1:nx, 1:ny, 1:nz))
    allocate(dlna_c(0:3, 1:nx, 1:ny, 1:nz))

    ! 2. Interfaces X (1:nx+1, 1:ny, 1:nz)
    if (nx > 1) then
      allocate(alpha_f_x(1:nx+1, 1:ny, 1:nz))
      allocate(beta_f_x(3, 1:nx+1, 1:ny, 1:nz))
      allocate(gamma_f_x(3, 1:nx+1, 1:ny, 1:nz))
      allocate(sqrt_gamma_f_x(1:nx+1, 1:ny, 1:nz))
    end if

    ! 3. Interfaces Y (1:nx, 1:ny+1, 1:nz)
    if (ny > 1) then
      allocate(alpha_f_y(1:nx, 1:ny+1, 1:nz))
      allocate(beta_f_y(3, 1:nx, 1:ny+1, 1:nz))
      allocate(gamma_f_y(3, 1:nx, 1:ny+1, 1:nz))
      allocate(sqrt_gamma_f_y(1:nx, 1:ny+1, 1:nz))
    end if

    ! 4. Interfaces Z (1:nx, 1:ny, 1:nz+1)
    if (nz > 1) then
      allocate(alpha_f_z(1:nx, 1:ny, 1:nz+1))
      allocate(beta_f_z(3, 1:nx, 1:ny, 1:nz+1))
      allocate(gamma_f_z(3, 1:nx, 1:ny, 1:nz+1))
      allocate(sqrt_gamma_f_z(1:nx, 1:ny, 1:nz+1))
    end if

  end subroutine allocate_and_grid

  subroutine precalculate_metric_cache()
    use variables
    implicit none
    integer :: i, j, k, ii, jj
    real*8  :: a, b(3), g(3), gmunu_loc(0:3,0:3)
    real*8  :: det_local, dlna_local(0:3)
    real*8  :: chris_local(0:3,0:3,0:3), dg_local(0:3,0:3,1:3)

    print *, "--> Pre-calculando Caché de Métrica 3+1 ..."

    ! ========================================================
    ! 1. MÉTRICA EN LOS CENTROS CELULARES
    ! ========================================================
    !$OMP PARALLEL DO PRIVATE(i, j, k, a, b, g, gmunu_loc, det_local, dlna_local, chris_local, dg_local) COLLAPSE(3)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! Llamada estricta por keywords. Extraemos dlnalpha desde aquí.
          call calculate_metric(x(i), y(j), alpha=a, beta=b, gamma=g, gmunu=gmunu_loc, det=det_local, dlnalpha=dlna_local)
          call calculate_christoffel_symbols(x(i), y(j), chris_local)
          call calculate_metric_derivatives(x(i), y(j), dg_local)
          
          alpha_c(i,j,k) = a
          beta_c(:,i,j,k) = b(:)
          gamma_c(:,i,j,k) = g(:)
          gmunu_c(:,:,i,j,k) = gmunu_loc(:,:)

          ! --- CÁLCULO DE LA INVERSA ESPACIAL (gamma^ij) PARA LA CACHÉ ---
          ! Fórmula ADM: gamma^ij = g^ij + (beta^i * beta^j) / alpha^2
          do ii = 1, 3
            do jj = 1, 3
              gamma_inv_c(ii,jj,i,j,k) = gmunu_loc(ii,jj) + &
                                    (b(ii) * b(jj)) / (a**2)
            end do
          end do
          
          sqrt_gamma_c(i,j,k) = sqrt(det_local)
          chris_c(:,:,:,i,j,k) = chris_local(:,:,:)
          dg_c(:,:,:,i,j,k)    = dg_local(:,:,:)
          dlna_c(:,i,j,k)      = dlna_local(:)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    ! ========================================================
    ! 2. MÉTRICA EN LAS INTERFACES X
    ! ========================================================
    if (nx > 1) then
      !$OMP PARALLEL DO PRIVATE(i, j, k, a, b, g, det_local) COLLAPSE(3)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx+1
            call calculate_metric(x_face(i-1), y(j), alpha=a, beta=b, gamma=g, det=det_local)
            
            alpha_f_x(i,j,k) = a
            beta_f_x(:,i,j,k) = b(:)
            gamma_f_x(:,i,j,k) = g(:)
            sqrt_gamma_f_x(i,j,k) = sqrt(det_local)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end if

    ! ========================================================
    ! 3. MÉTRICA EN LAS INTERFACES Y
    ! ========================================================
    if (ny > 1) then
      !$OMP PARALLEL DO PRIVATE(i, j, k, a, b, g, det_local) COLLAPSE(3)
      do k = 1, nz
        do j = 1, ny+1
          do i = 1, nx
            call calculate_metric(x(i), y_face(j-1), alpha=a, beta=b, gamma=g, det=det_local)
            
            alpha_f_y(i,j,k) = a
            beta_f_y(:,i,j,k) = b(:)
            gamma_f_y(:,i,j,k) = g(:)
            sqrt_gamma_f_y(i,j,k) = sqrt(det_local)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end if

    ! ========================================================
    ! 4. MÉTRICA EN LAS INTERFACES Z
    ! ========================================================
    if (nz > 1) then
      !$OMP PARALLEL DO PRIVATE(i, j, k, a, b, g, det_local) COLLAPSE(3)
      do k = 1, nz+1
        do j = 1, ny
          do i = 1, nx
            call calculate_metric(x(i), y(j), alpha=a, beta=b, gamma=g, det=det_local)
            
            alpha_f_z(i,j,k) = a
            beta_f_z(:,i,j,k) = b(:)
            gamma_f_z(:,i,j,k) = g(:)
            sqrt_gamma_f_z(i,j,k) = sqrt(det_local)
          end do
        end do
      end do
      !$OMP END PARALLEL DO
    end if

  end subroutine precalculate_metric_cache

  ! Subrutina: initialize_problem
  ! Punto de entrada principal invocado desde grhd2.f90.
  ! Orquesta la llamada a todas las rutinas de setup en el orden topológico correcto.
  subroutine initialize_problem()

    print *, "=========================================="
    print *, "Select test case:"
    print *, "1. Sod Shock Tube (1D SRHD)"
    print *, "2. Blast Wave (1D SRHD)"
    print *, "3. Shu-Osher (1D SRHD)"
    print *, "4. Kelvin-Helmholtz Instability (2D SRHD)"
    print *, "5. Axisymmetric Jet (2D SRHD)"
    print *, "6. Convergence Test"
    print *, "7. Michel Accretion (1D GRHD)"
    print *, "8. Dust Accretion (1D GRHD)"
    print *, "9. Off-axis Blast Wave (2D GRHD)"
    print *, "10. Fishbone-Moncrief Torus (2D Equatorial GRHD)"
    print *, "11. Fishbone-Moncrief Torus (2D Sagital GRHD)"
    print *, "12. Fishbone-Moncrief Torus (3D GRHD)"
    
    ! RECORDATORIO: Descomentar para producción interactiva
    ! read *, case_id
    case_id = 10

    ! 1. CONFIGURACIÓN DEL PROBLEMA (Cada rutina define sus propios parámetros y geometría)
    select case(case_id)
    case(1)
      call setup_sod()
      print *, "1. Sod Shock Tube (1D SRHD) selected"
    case(2)
      call setup_strong_shock()
      print *, "2. Blast Wave (1D SRHD) selected"
    case(3)
      call setup_shu_osher()
      print *, "3. Shu-Osher (1D SRHD) selected"
    case(4)
      call setup_kelvin_helmholtz()
      print *, "4. Kelvin-Helmholtz Instability (2D SRHD) selected"
    case(5)
      call setup_axisymmetric_jet()
      print *, "5. Axisymmetric Jet (2D SRHD) selected"
    case(6)
      call setup_convergence_test()
      print *, "6. Convergence Test selected"
      print *, 'Sensor de choques apagado'
    case(7)
      call setup_michel_accretion()
      print *, "7. Michel Accretion (1D GRHD) selected"
    case(8)
      call setup_dust_accretion()
      print *, "8. Dust Accretion (1D GRHD) selected"
    case(9)
      call setup_offaxis_blast()
      print *, "9. Off-axis Blast Wave (2D GRHD) selected"
    case(10)
      call setup_fishbone_moncrief_equatorial()
      print *, "10. Fishbone-Moncrief Torus (2D Equatorial GRHD) selected"
    case(11)
      call setup_fishbone_moncrief_sagital()
      print *, "11. Fishbone-Moncrief Torus (2D Sagital GRHD) selected"
    case(12)
      call setup_fishbone_moncrief()
      print *, "12. Fishbone-Moncrief Torus (3D GRHD) selected"
    case default
      print *, "Invalid selection. Exiting."
      stop
    end select

    ! 2. ASIGNACIÓN DE MÉTRICA
    call set_metric_type()

    ! Inicialización de constantes termodinámicas universales
    g1 = adb_idx / (adb_idx - 1.0d0)
    
    rho_floor = 1.0d-10
    p_floor   = rho_floor * 1.0d-3
    D_floor   = rho_floor
    tau_floor = p_floor / (g1 - 1.0d0)

    ! 3. ARQUITECTURA NUMÉRICA (Define 'scheme_name' y 'nghost')
    call choose_numerical_architecture() 

    ! 4. CREACIÓN DE DIRECTORIOS
    call system('mkdir -p ' // trim(output_folder) //  '/' // trim(scheme_name))

    ! 5. ALOJAMIENTO DE MEMORIA, CONDICIONES INICIALES Y CALCULO DE EIGENVALORES
    call allocate_and_grid()
    call precalculate_metric_cache()
    call set_initial_conditions()
    call init_wavespeed_solver()

    ! 6. INICIALIZACIÓN DE DIAGNÓSTICOS MULTIMENSAJERO
    ! Solo se abren los archivos si la gravedad curvo-espacial está activada
    if (trim(metric_type) /= 'Minkowski' .and. bh_mass > 0.0d0) then
        block
          character(len=150) :: gw_filename, m_dot_filename
          if (do_gw_extraction) then
            write(gw_filename, '(A,A,A,A,A)') trim(output_folder), '/', trim(scheme_name), '/', 'GW_signal.dat'
            open(30, file=gw_filename, status='replace')
            write(30, '(A)') '# Time            h_plus            h_cross'
            close(30)
          end if
          
          if (do_mdot_extraction) then
            write(m_dot_filename, '(A,A,A,A,A)') trim(output_folder), '/', trim(scheme_name), '/', 'm_dot.dat'
            open(88, file=m_dot_filename, status='replace')
            write(88, '(A)') '# Time            M_dot'
            close(88)
          end if
        end block
    end if
  end subroutine initialize_problem

  ! ===================================================================================
  ! RUTINAS DE CONFIGURACIÓN DE PARÁMETROS (Tiempos, Resoluciones y EoS)
  ! ===================================================================================

  subroutine setup_sod()
    case_name = 'Sod'
    geom_type = 'Cartesian'
    metric_type = 'Minkowski'
    use_shock_sensor = .false.
    do_mdot_extraction = .false.
    do_gw_extraction = .false.

    adb_idx = 1.4d0 
    
    nx = 200 ; r_min = 0.0d0 ; r_max = 1.0d0
    ny = 1   ; y_min = 0.0d0 ; y_max = 1.0d0
    nz = 1   ; z_min = 0.0d0 ; z_max = 1.0d0
    
    final_time = 0.4d0
    CFL = 0.4d0
    save_interval = 0.1d0
    
    output_prefix = 'Sod'
    output_folder = 'ShockTubes'
  end subroutine setup_sod

  subroutine setup_strong_shock()
    case_name = 'Strong'
    geom_type = 'Cartesian'
    metric_type = 'Minkowski'
    use_shock_sensor = .false.
    do_mdot_extraction = .false.
    do_gw_extraction = .false.

    adb_idx = 5.0d0/3.0d0 
    
    nx = 400 ; r_min = 0.0d0 ; r_max = 1.0d0
    ny = 1   ; y_min = 0.0d0 ; y_max = 1.0d0
    nz = 1   ; z_min = 0.0d0 ; z_max = 1.0d0
    
    final_time = 0.4d0
    CFL = 0.4d0
    save_interval = 0.1d0
    
    output_prefix = 'Strong'
    output_folder = 'ShockTubes'
  end subroutine setup_strong_shock

  subroutine setup_shu_osher()
    case_name = 'Shu'
    geom_type = 'Cartesian'
    metric_type = 'Minkowski'
    use_shock_sensor = .false.
    do_mdot_extraction = .false.
    do_gw_extraction = .false.

    adb_idx = 5.0d0/3.0d0
    
    nx = 400 ; r_min = 0.0d0 ; r_max = 1.0d0
    ny = 1   ; y_min = 0.0d0 ; y_max = 1.0d0
    nz = 1   ; z_min = 0.0d0 ; z_max = 1.0d0
    
    final_time = 0.4d0
    CFL = 0.4d0
    save_interval = 0.1d0
    
    output_prefix = 'Shu'
    output_folder = 'ShockTubes'
  end subroutine setup_shu_osher

  subroutine setup_kelvin_helmholtz()
    case_name = 'KHI'
    geom_type = 'Cartesian'
    metric_type = 'Minkowski'
    use_shock_sensor = .false.
    do_mdot_extraction = .false.
    do_gw_extraction = .false.

    adb_idx = 5.0d0 / 3.0d0 ! Gamma = 5/3 (Monoatómico estándar para astrofísica)
    
    nx = 400 ; r_min = -0.5d0 ; r_max = 0.5d0
    ny = 1   ; y_min = 0.0d0  ; y_max = 1.0d0
    nz = 400 ; z_min = -0.5d0 ; z_max = 0.5d0

    final_time = 3.0d0 
    CFL = 0.25d0       ! CFL conservador para domar inestabilidades cruzadas
    save_interval = 0.1d0
    
    output_prefix = 'KHI'
    output_folder = 'KelvinHelmholtz'
  end subroutine setup_kelvin_helmholtz

  subroutine setup_axisymmetric_jet()
    case_name = 'Jet'
    geom_type = 'Cylindrical'
    metric_type = 'Minkowski'
    use_shock_sensor = .true.
    do_mdot_extraction = .false.
    do_gw_extraction = .false.

    adb_idx = 5.0d0 / 3.0d0
    
    ! Dominio de Del Zanna (Astrofísica Clásica)
    nx = 320 ; r_min = 0.0d0 ; r_max = 16.0d0
    ny = 1   ; y_min = 0.0d0 ; y_max = 1.0d0
    nz = 800 ; z_min = 0.0d0 ; z_max = 40.0d0

    final_time = 120.0d0 
    CFL = 0.25d0        
    save_interval = 1.0d0
    
    output_prefix = 'JET'
    output_folder = 'AxisymmetricJet_extreme'
  end subroutine setup_axisymmetric_jet

  subroutine setup_convergence_test()
    case_name = 'Conv'
    geom_type = 'Cartesian'
    metric_type = 'Minkowski'
    use_shock_sensor = .false.
    do_mdot_extraction = .false.
    do_gw_extraction = .false.

    adb_idx = 4.0d0 / 3.0d0
    
    nx = 2048 ; r_min = 0.0d0 ; r_max = 1.0d0
    ny = 1    ; y_min = 0.0d0 ; y_max = 1.0d0
    nz = 1    ; z_min = 0.0d0 ; z_max = 1.0d0
    
    final_time = 0.7d0
    CFL = 0.10d0
    save_interval = 0.1d0
    
    output_prefix = 'Conv'
    output_folder = 'ConvergenceTest_N2048'
  end subroutine setup_convergence_test

  subroutine setup_michel_accretion()
    case_name = 'MichelA'
    geom_type = 'Spherical'
    metric_type = 'Eddington-Finkelstein'
    use_shock_sensor = .false.
    do_mdot_extraction = .true.
    do_gw_extraction = .false.

    adb_idx = 4.0d0 / 3.0d0 ! Fluido dominado por radiación
    bh_mass = 1.0d0
    
    nx = 400 ; r_min = 1.01d0 ; r_max = 51.0d0
    ny = 1   ; y_min = 0.0d0 ; y_max = pi
    nz = 1   ; z_min = 0.0d0 ; z_max = 2.0d0 * pi

    ! Activa esto cuando quieras máxima resolución cerca del horizonte
    ! use_log_r = .true. 
    
    final_time = 20.0d0
    CFL = 0.5d0
    save_interval = 1.0d0
    
    output_prefix = 'MA'
    output_folder = 'Michel_test'
  end subroutine setup_michel_accretion

  subroutine setup_dust_accretion()
    case_name = 'Dust'
    geom_type = 'Spherical'
    metric_type = 'Eddington-Finkelstein'
    use_shock_sensor = .true.
    do_mdot_extraction = .true.
    do_gw_extraction = .false.

    adb_idx = 4.0d0 / 3.0d0
    bh_mass = 1.0d0
    
    nx = 500 ; r_min = 1.0d0 ; r_max = 51.0d0
    ny = 1   ; y_min = 0.0d0 ; y_max = pi
    nz = 4   ; z_min = 0.0d0 ; z_max = 2.0d0 * pi 
    
    final_time = 500.0d0
    CFL = 0.5d0
    save_interval = 50.0d0
    
    output_prefix = 'D'
    output_folder = 'Dust2D'
  end subroutine setup_dust_accretion

  subroutine setup_offaxis_blast()
    case_name = 'OffAxis'
    geom_type = 'Spherical'
    metric_type = 'Eddington-Finkelstein'
    use_shock_sensor = .true.
    do_mdot_extraction = .false.
    do_gw_extraction = .true.

    adb_idx = 4.0d0 / 3.0d0
    bh_mass = 1.0d0
    
    nx = 400 ; r_min = 1.0d0 ; r_max = 121.0d0
    ny = 1   ; y_min = 0.0d0 ; y_max = pi
    nz = 200 ; z_min = 0.0d0 ; z_max = 2.0*pi

    final_time = 2000.0d0
    CFL = 0.4d0
    save_interval = 25.0d0
    
    output_prefix = 'OR'
    output_folder = 'OffAxis'
  end subroutine setup_offaxis_blast

  subroutine setup_fishbone_moncrief_equatorial()
    case_name = 'FishMoncEqu'
    geom_type = 'Spherical'
    metric_type = 'Eddington-Finkelstein'
    use_shock_sensor = .true.
    do_mdot_extraction = .true.
    do_gw_extraction = .true.
    use_log_r = .true.

    adb_idx = 4.0d0 / 3.0d0
    bh_mass = 1.0d0
    
    ! Caso original Ecuatorial (1D/2D en phi)
    nx = 400 ; r_min = 2.0d0 ; r_max = 200.0d0
    ny = 1   ; y_min = 0     ; y_max = pi     ! Fijo en el ecuador
    nz = 200 ; z_min = 0.0d0 ; z_max = 2.0*pi

    final_time = 15000.0d0
    CFL = 0.4d0            
    save_interval = 10.0d0   
    
    output_prefix = 'FM_Equ'
    output_folder = 'FishboneMoncrief_Equatorial_test'
  end subroutine setup_fishbone_moncrief_equatorial

  subroutine setup_fishbone_moncrief_sagital()
    case_name = 'FishMoncSag'
    geom_type = 'Spherical'
    metric_type = 'Eddington-Finkelstein'
    use_shock_sensor = .true.
    do_mdot_extraction = .true.
    do_gw_extraction = .true.
    use_log_r = .true.

    adb_idx = 4.0d0 / 3.0d0
    bh_mass = 1.0d0
    
    nx = 400 ; r_min = 2.0d0 ; r_max = 200.0d0
    ny = 200 ; y_min = 0.0d0 ; y_max = pi      ! De Polo a Polo
    nz = 1   ; z_min = 0.0d0 ; z_max = 2.0d0 * pi ! Simetría Azimutal

    final_time = 0.0d0
    CFL = 0.4d0            
    save_interval = 1.0d0   
    
    output_prefix = 'FM_Sag'
    output_folder = 'FishboneMoncrief_Sagital'
  end subroutine setup_fishbone_moncrief_sagital

  subroutine setup_fishbone_moncrief()
    case_name = 'FishMonc3D'
    geom_type = 'Spherical'
    metric_type = 'Eddington-Finkelstein'
    use_shock_sensor = .true.
    do_mdot_extraction = .true.
    do_gw_extraction = .true.
    ! use_log_r = .true.

    adb_idx = 4.0d0 / 3.0d0
    bh_mass = 1.0d0
    
    nx = 100 ; r_min = 1.0d0 ; r_max = 121.0d0
    ny = 100 ; y_min = 0.0d0 ; y_max = pi
    nz = 100 ; z_min = 0.0d0 ; z_max = 2.0d0 * pi

    final_time = 35.0d0
    CFL = 0.4d0            
    save_interval = 1.0d0   
    
    output_prefix = 'FM_3D'
    output_folder = 'FishboneMoncrief_3D'
  end subroutine setup_fishbone_moncrief
  
end module initialization