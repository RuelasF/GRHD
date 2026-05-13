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
    real*8, allocatable :: gamma_c(:,:,:,:)         ! Métrica espacial [1:3, 1:3, x, y, z]
    real*8, allocatable :: gmunu_c(:,:,:,:,:)       ! Inversa métrica espacial
    real*8, allocatable :: sqrt_gamma_c(:,:,:)      ! Raíz del determinante
    real*8, allocatable :: chris_c(:,:,:,:,:,:)     ! Símbolos de Christoffel (0:3, 0:3, 0:3, x, y, z)
    real*8, allocatable :: dg_c(:,:,:,:,:,:)        ! Derivadas de la métrica (0:3, 0:3, 1:3, x, y, z)
    real*8, allocatable :: dlna_c(:,:,:,:)        ! Derivada del logaritmo del lapso (0:3, x, y, z)

    ! 2. Evaluados en las INTERFACES X (Para flujos F^x)
    real*8, allocatable :: alpha_f_x(:,:,:)
    real*8, allocatable :: beta_f_x(:,:,:,:)
    real*8, allocatable :: gamma_f_x(:,:,:,:)
    real*8, allocatable :: gmunu_f_x(:,:,:,:,:)
    real*8, allocatable :: sqrt_gamma_f_x(:,:,:)

    ! 3. Evaluados en las INTERFACES Y (Para flujos F^y)
    real*8, allocatable :: alpha_f_y(:,:,:)
    real*8, allocatable :: beta_f_y(:,:,:,:)
    real*8, allocatable :: gamma_f_y(:,:,:,:)
    real*8, allocatable :: gmunu_f_y(:,:,:,:,:)
    real*8, allocatable :: sqrt_gamma_f_y(:,:,:)

    ! 4. Evaluados en las INTERFACES Z (Para flujos F^z)
    real*8, allocatable :: alpha_f_z(:,:,:)
    real*8, allocatable :: beta_f_z(:,:,:,:)
    real*8, allocatable :: gamma_f_z(:,:,:,:)
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

    ! ==========================================================
    ! LÍMITES NUMÉRICOS (Atmósferas y Tolerancias)
    ! ==========================================================
    real*8, parameter :: tol_v = 1.0d-5 
    real*8, parameter :: v_max = 1.0d0 - tol_v 

    real*8 :: rho_floor, p_floor 
    real*8 :: D_floor, tau_floor 

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
    integer :: riemann_solver_id 

    integer, parameter :: LIM_MINMOD   = 1
    integer, parameter :: LIM_SUPERBEE = 2
    integer, parameter :: LIM_MC       = 3

    integer :: rec_method_id 
    integer :: tvd_limiter_id 
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

    real*8 :: save_interval 
    real*8 :: next_save_time 

    real*8, allocatable :: michel_injector(:,:,:,:) 

    real*8 :: checkpoint_interval = 100.0d0 
    real*8 :: next_checkpoint 
    logical :: do_restart = .false. 
    character(len=100) :: restart_file = 'FishboneMoncrief_Equatorial3/checkpoint_weno5_hlle_00100.rst' 

end module variables