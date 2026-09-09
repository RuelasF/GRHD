! =================================================================================
! Módulo: variables
! Descripción: Contenedor global de las variables de estado, parámetros físicos, 
! arreglos geométricos y banderas de control del simulador. 
!
! Soporte Multi-Geometría:
! Las variables espaciales x, y, z representan las dimensiones lógicas x1, x2, x3.
! Su significado físico (Cartesianas, Cilíndricas o Esféricas) es dictaminado 
! dinámicamente por el módulo de métricas.
!
! Convención Geométrica Estándar (CFD):
! - Celdas Físicas (Centros): Índices de [1 : nx] y [1 : nz]
! - Interfaces (Caras):       Índices de [0 : nx] y [0 : nz]
! - Celdas Fantasma (Ghosts): Índices [-nghost : 0] y [nx+1 : nx+nghost]
! =================================================================================

module variables
    implicit none
    public

    ! ==========================================================
    ! ÍNDICES DE ECUACIONES (Vectores de Estado)
    ! ==========================================================
    integer, parameter :: neq = 5   
    
    integer :: eq_de = 1 ! D:  Densidad de masa en reposo conservada (rho * W * sqrt(gamma))
    integer :: eq_pr = 2 ! tau: Densidad de energía conservada (excluyendo masa en reposo)
    integer :: eq_vx = 3 ! S_1: Momento en la 1ra dimensión (ej. x, r)
    integer :: eq_vy = 4 ! S_2: Momento en la 2da dimensión (ej. y, theta)
    integer :: eq_vz = 5 ! S_3: Momento en la 3ra dimensión (ej. z, phi)

    character(len=20), allocatable :: var_names(:)

    ! ==========================================================
    ! DOMINIO NUMÉRICO Y CONTROL TEMPORAL
    ! ==========================================================
    integer :: nx, ny, nz                   ! Número estricto de celdas FÍSICAS activas
    real*8  :: x_min, x_max                 ! Límites físicos de la 1ra dimensión lógica
    real*8  :: r_min, r_max                 ! Límites físicos de la 1ra dimensión física (r)
    real*8  :: y_min, y_max                 ! Límites físicos de la 2da dimensión
    real*8  :: z_min, z_max                 ! Límites físicos de la 3ra dimensión
    real*8  :: dx, dy, dz                   ! Espaciado constante de la malla lógica
    
    real*8  :: t, dt, final_time            ! Variables de integración temporal
    real*8  :: integration_time             ! Tiempo actual de la simulación
    real*8  :: CFL                          ! Número de Courant-Friedrichs-Lewy
    integer :: n_steps = 0                  ! Contador total de pasos evaluados

    ! ==========================================================
    ! ARREGLOS GEOMÉTRICOS DE LA MALLA (1D Arrays)
    ! ==========================================================
    real*8, allocatable :: x(:), y(:), z(:) ! Coordenadas en los centros celulares
    real*8, allocatable :: x_face(:)        ! Coordenadas de las interfaces 1D
    real*8, allocatable :: y_face(:)        ! Coordenadas de las interfaces 2D
    real*8, allocatable :: z_face(:)        ! Coordenadas de las interfaces 3D

    ! ==========================================================
    ! ARREGLOS DE ESTADO (HIDRODINÁMICA)
    ! ==========================================================
    real*8, allocatable :: u(:,:,:,:)       ! Estado conservativo paso 'n'
    real*8, allocatable :: up(:,:,:,:)      ! Estado intermedio/final paso 'n+1'
    real*8, allocatable :: s(:,:,:,:)       ! Términos fuente
    real*8, allocatable :: rhs(:,:,:,:)     ! Lado Derecho (RHS)
    real*8, allocatable :: p(:,:,:,:)       ! Variables Primitivas

    ! ==========================================================
    ! ARREGLOS DE MÉTRICA PRE-CALCULADA (CACHÉ 3+1)
    ! ==========================================================
    ! 1. Evaluados en los CENTROS celulares (Para fuentes y conversión U -> P)
    real*8, allocatable :: alpha_c(:,:,:)           ! Lapse function
    real*8, allocatable :: beta_c(:,:,:,:)          ! Shift vector [1:3, x, y, z]
    real*8, allocatable :: gamma_c(:,:,:,:,:)       ! Métrica espacial [1:3, 1:3, x, y, z]
    real*8, allocatable :: gamma_inv_c(:,:,:,:,:)   ! Inversa métrica espacial [1:3, 1:3, x, y, z]
    real*8, allocatable :: gmunu_c(:,:,:,:,:)       ! Inversa métrica espacio-temporal
    real*8, allocatable :: sqrt_gamma_c(:,:,:)      ! Raíz del determinante
    real*8, allocatable :: chris_c(:,:,:,:,:,:)     ! Símbolos de Christoffel (0:3, 0:3, 0:3, x, y, z)
    real*8, allocatable :: dg_c(:,:,:,:,:,:)        ! Derivadas de la métrica (0:3, 0:3, 1:3, x, y, z)
    real*8, allocatable :: dlna_c(:,:,:,:)          ! Derivada del logaritmo del lapso (0:3, x, y, z)

    ! 2. Evaluados en las INTERFACES X (Para flujos F^x)
    real*8, allocatable :: alpha_f_x(:,:,:)
    real*8, allocatable :: beta_f_x(:,:,:,:)
    real*8, allocatable :: gamma_f_x(:,:,:,:,:)
    real*8, allocatable :: gmunu_f_x(:,:,:,:,:)
    real*8, allocatable :: sqrt_gamma_f_x(:,:,:)

    ! 3. Evaluados en las INTERFACES Y (Para flujos F^y)
    real*8, allocatable :: alpha_f_y(:,:,:)
    real*8, allocatable :: beta_f_y(:,:,:,:)
    real*8, allocatable :: gamma_f_y(:,:,:,:,:)
    real*8, allocatable :: gmunu_f_y(:,:,:,:,:)
    real*8, allocatable :: sqrt_gamma_f_y(:,:,:)

    ! 4. Evaluados en las INTERFACES Z (Para flujos F^z)
    real*8, allocatable :: alpha_f_z(:,:,:)
    real*8, allocatable :: beta_f_z(:,:,:,:)
    real*8, allocatable :: gamma_f_z(:,:,:,:,:)
    real*8, allocatable :: gmunu_f_z(:,:,:,:,:)
    real*8, allocatable :: sqrt_gamma_f_z(:,:,:)

    ! ==========================================================
    ! DIRECCIONES DE BARRIDO (Dimensional Splitting)
    ! ==========================================================
    integer, parameter :: DIR_X = 1
    integer, parameter :: DIR_Y = 2
    integer, parameter :: DIR_Z = 3

    ! ==========================================================
    ! PARÁMETROS FÍSICOS Y TERMODINÁMICOS
    ! ==========================================================
    real*8 :: adb_idx                       ! Índice adiabático (Gamma)
    real*8 :: g1                            ! Gamma / (Gamma - 1.0)
    real*8 :: pi = acos(-1.0d0) 
    
    real*8 :: bh_mass = 0.0d0               ! Masa del Agujero Negro
    real*8 :: a_spin  = 0.0d0               ! Parámetro de Espín de Kerr (a = J/M)

    integer :: idx_horizon                  ! Índice del horizonte
    integer :: idx_probe                    ! Índice celda monitoreo

    ! ==========================================================
    ! BANDERAS DE CONTROL GLOBAL
    ! ==========================================================
    character(len=20) :: geom_type               ! 'Cartesian', 'Cylindrical', 'Spherical'
    logical :: use_log_r = .false.               ! .true. comprime 1ra dimensión

    character(len=30) :: metric_type             ! 'Minkowski', 'Eddington-Finkelstein', 'Kerr'
    
    logical :: is_mhd = .false.                  ! Campos B (Futuro)
    logical :: use_shock_sensor = .true.         ! Fallback reconstrucción

    logical :: do_mdot_extraction = .false.      ! Extracción de tasa de acreción
    logical :: do_gw_extraction = .false.        ! Extracción de ondas gravitacionales
    logical :: do_ppi_diagnostics = .false.      ! Modos azimutales y estado global

    ! Control reproducible de perturbaciones para toros de Fishbone--Moncrief.
    integer, parameter :: PERT_NONE = 0
    integer, parameter :: PERT_PRESSURE_NOISE = 1
    integer, parameter :: PERT_DENSITY_NOISE = 2
    integer, parameter :: PERT_DENSITY_MODE = 3
    integer, parameter :: PPI_MAX_MODE = 4

    logical :: apply_perturbation = .false.
    logical :: perturbation_applied = .false.
    integer :: perturbation_type = PERT_PRESSURE_NOISE
    integer :: perturbation_seed = 3435
    integer :: diagnostic_stride = 100
    real*8 :: perturbation_time = 1000.0d0
    real*8 :: perturbation_amplitude = 0.01d0
    real*8 :: perturbation_mode = 1.0d0

    ! ==========================================================
    ! LÍMITES NUMÉRICOS (Atmósferas y Tolerancias)
    ! ==========================================================
    real*8, parameter :: tol_v = 1.0d-5 
    real*8, parameter :: v_max = 1.0d0 - tol_v 

    real*8 :: rho_floor, p_floor
    real*8 :: D_floor, tau_floor

    ! ==========================================================
    ! PARÁMETROS DE LAS CONDICIONES INICIALES
    ! ==========================================================
    ! Tubos de choque.
    real*8 :: sod_interface = 0.5d0
    real*8 :: sod_rho_left = 1.0d0, sod_rho_right = 0.125d0
    real*8 :: sod_pressure_left = 1.0d0, sod_pressure_right = 0.1d0

    real*8 :: strong_interface = 0.5d0
    real*8 :: strong_rho_left = 1.0d0, strong_rho_right = 1.0d0
    real*8 :: strong_pressure_left = 1000.0d0, strong_pressure_right = 0.01d0

    real*8 :: shu_interface = 0.5d0
    real*8 :: shu_rho_left = 5.0d0, shu_pressure_left = 50.0d0
    real*8 :: shu_rho_right = 2.0d0, shu_rho_amplitude = 0.3d0
    real*8 :: shu_wave_number = 50.0d0, shu_pressure_right = 5.0d0

    ! Kelvin--Helmholtz.
    real*8 :: khi_half_width = 0.25d0
    real*8 :: khi_rho_inner = 2.0d0, khi_rho_outer = 1.0d0
    real*8 :: khi_vx_inner = 0.5d0, khi_vx_outer = -0.5d0
    real*8 :: khi_pressure = 2.5d0
    real*8 :: khi_perturbation_amplitude = 0.01d0
    real*8 :: khi_perturbation_wave_number = 10.0d0

    ! Jet axisimétrico.
    real*8 :: jet_ambient_density = 10.0d0
    real*8 :: jet_ambient_pressure = 0.01d0
    real*8 :: jet_ambient_velocity = 0.0d0
    real*8 :: jet_nozzle_radius = 1.0d0, jet_nozzle_length = 1.0d0
    real*8 :: jet_density = 0.1d0, jet_pressure = 0.01d0
    real*8 :: jet_velocity = 0.99d0

    ! Onda suave para convergencia. wave_number está en múltiplos de pi.
    real*8 :: advected_wave_density = 1.0d0
    real*8 :: advected_wave_amplitude = 0.1d0
    real*8 :: advected_wave_speed = 0.5d0
    real*8 :: advected_wave_pressure = 1.0d0
    real*8 :: advected_wave_number = 2.0d0

    ! Acreción y explosión fuera del eje.
    real*8 :: michel_critical_radius = 100.0d0
    real*8 :: michel_critical_density = 0.1d0
    real*8 :: dust_accretion_constant = -0.5d0
    real*8 :: offaxis_radial_center = 15.0d0
    real*8 :: offaxis_phi_center = acos(-1.0d0)
    real*8 :: offaxis_width = 1.5d0
    real*8 :: offaxis_background_density = 0.1d0
    real*8 :: offaxis_background_pressure = 0.1d0
    real*8 :: offaxis_density_amplitude = 1.0d0
    real*8 :: offaxis_pressure_amplitude = 1.0d0

    ! Fishbone--Moncrief: radios expresados en unidades de la masa M.
    real*8 :: fm_inner_radius = 6.25d0
    real*8 :: fm_pressure_max_radius = 9.25d0
    real*8 :: fm_polytropic_constant = 0.0015d0

    ! ==========================================================
    ! ARQUITECTURA NUMÉRICA (Reconstrucción)
    ! ==========================================================
    integer, parameter :: REC_GODUNOV = 1
    integer, parameter :: REC_TVD     = 2
    integer, parameter :: REC_WENO3   = 3
    integer, parameter :: REC_MP5     = 4
    integer, parameter :: REC_WENO5   = 5

    integer, parameter :: RS_HLLE = 1
    integer, parameter :: RS_HLLC = 2
    integer, parameter :: RS_HLLD = 3
    integer :: riemann_solver_id = RS_HLLE

    integer, parameter :: LIM_MINMOD   = 1
    integer, parameter :: LIM_SUPERBEE = 2
    integer, parameter :: LIM_MC       = 3

    integer :: rec_method_id = REC_WENO5
    integer :: tvd_limiter_id = LIM_MC
    character(len=20) :: scheme_name 
    character(len=10) :: solver_name 
    integer :: nghost 
    
    ! ==========================================================
    ! ENTRADA/SALIDA (I/O) Y CHECKPOINTS
    ! ==========================================================
    integer :: case_id 
    character(len=50) :: case_name 
    character(len=50) :: output_folder 
    character(len=50) :: output_prefix 

    ! El estado siempre se calcula en las coordenadas de la metrica. Esta
    ! opcion solo controla la geometria cartesiana escrita en los VTK.
    integer, parameter :: VTK_MAP_PHYSICAL = 1
    integer, parameter :: VTK_MAP_UNTWISTED = 2
    integer :: vtk_mapping_id = VTK_MAP_PHYSICAL

    real*8 :: save_interval 
    real*8 :: next_save_time 

    real*8, allocatable :: michel_injector(:,:,:,:) 

    real*8 :: checkpoint_interval = 100.0d0 
    real*8 :: next_checkpoint 
    logical :: do_restart = .false. 
    character(len=100) :: restart_file = 'FishboneMoncrief_Equatorial_test/checkpoint_weno5_hlle_08800.rst' 

end module variables
