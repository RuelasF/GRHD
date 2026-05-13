! =======================================================================================
! Módulo: metrics (Geometría del Espacio-Tiempo y Gravedad)
! ---------------------------------------------------------------------------------------
! Este módulo define la métrica del espacio-tiempo. Proporciona los factores geométricos 
! (alpha, beta, gamma) de la descoposición 3+1 (ADM), los símbolos de Christoffel (fuerzas 
! gravitacionales) y las derivadas de la métrica necesarias para los términos fuente del RHS.
!
! Soporte de Métricas:
! 1. Minkowski (Cartesiano 1D/2D Plano)
! 2. Minkowski Cilíndrico (Axisimétrico sin gravedad)
! 3. Schwarzschild en coordenadas de Eddington-Finkelstein (EF).
!    ¿Por qué EF?: Las coordenadas de Schwarzschild estándar tienen una singularidad 
!    de coordenadas en el horizonte de eventos (r = 2M). Las coordenadas EF son 
!    "penetrantes", permitiendo que el flujo hidrodinámico cruze suavemente el horizonte.
! =======================================================================================

module metrics
  use variables
  implicit none
  private
  
  ! ========================================================================
  ! EXPOSICIÓN PÚBLICA
  ! Las subrutinas matemáticas son punteros dinámicos. Esto evita usar 'if's 
  ! dentro de los bucles de integración temporal, ahorrando tiempo de CPU.
  ! ========================================================================
  public :: calculate_metric, calculate_christoffel_symbols, &
            calculate_metric_derivatives, calculate_stress_energy_tensor, &
            set_metric_type

  ! ========================================================================
  ! INTERFACES ABSTRACTAS (Plantillas para los punteros de función)
  ! ========================================================================
  abstract interface
    subroutine i_metric(x_pos, y_pos, alpha, beta, gamma, gmunu, det, dlnalpha)
      implicit none
      real*8, intent(in) :: x_pos, y_pos
      real*8, intent(out), optional :: alpha, beta(3), gamma(3), det
      real*8, intent(out), optional :: gmunu(0:3, 0:3), dlnalpha(0:3)
    end subroutine i_metric

    subroutine i_christoffel(x_pos, y_pos, chris)
      implicit none
      real*8, intent(in) :: x_pos, y_pos
      real*8, intent(out) :: chris(0:3,0:3,0:3)
    end subroutine i_christoffel

    subroutine i_metric_derivs(x_pos, y_pos, dg)
      implicit none
      real*8, intent(in) :: x_pos, y_pos
      real*8, intent(out) :: dg(0:3, 0:3, 1:3)
    end subroutine i_metric_derivs
  end interface

  ! ========================================================================
  ! PUNTEROS A PROCEDIMIENTOS GLOBALES
  ! Estos son los nombres genéricos que llama el integrador principal.
  ! ========================================================================
  procedure(i_metric), pointer :: calculate_metric => null()
  procedure(i_christoffel), pointer :: calculate_christoffel_symbols => null()
  procedure(i_metric_derivs), pointer :: calculate_metric_derivatives => null()

contains

  ! ========================================================================
  ! INICIALIZADOR DE LA MÉTRICA (LLAMADO UNA VEZ EN EL ARRANQUE)
  ! Configura los punteros según la métrica y topología elegidas.
  ! ========================================================================
  subroutine set_metric_type()
    implicit none
    
    select case(trim(metric_type))
        
      ! ---------------------------------------------------------
      ! ESPACIO PLANO (Sin Gravedad)
      ! ---------------------------------------------------------
      case ('Minkowski')
        
        select case(trim(geom_type))
          case ('Cartesian')
            calculate_metric => metric_minkowski
            calculate_christoffel_symbols => christoffel_minkowski
            calculate_metric_derivatives => metric_derivs_minkowski
            print *, "============================================"
            print *, " METRIC ASSIGNED: Flat Minkowski (Cartesian)"
            print *, "============================================"
            
          case ('Cylindrical')
            calculate_metric => metric_cylindrical
            calculate_christoffel_symbols => christoffel_cylindrical
            calculate_metric_derivatives => metric_derivs_cylindrical
            print *, "============================================"
            print *, " METRIC ASSIGNED: Flat Cylindrical (rho, z) "
            print *, "============================================"
            
          case ('Spherical')
            ! Placeholder por si en el futuro agregas espacio plano en esféricas
            print *, "CRITICAL ERROR: Flat Spherical metric not yet implemented."
            stop
            
          case default
            print *, "CRITICAL ERROR: Unrecognized geom_type for Minkowski."
            stop
        end select

      ! ---------------------------------------------------------
      ! AGUJEROS NEGROS SIN ROTACIÓN
      ! ---------------------------------------------------------
      case ('Eddington-Finkelstein')
        
        if (trim(geom_type) /= 'Spherical') then
          print *, "CRITICAL ERROR: EF Metric strictly requires Spherical geometry."
          stop
        end if
        
        if (use_log_r) then
          calculate_metric => metric_ef_log
          calculate_christoffel_symbols => christoffel_ef_log
          calculate_metric_derivatives => metric_derivs_ef_log
          print *, "============================================"
          print *, " METRIC ASSIGNED: EF (Logarithmic Spherical)"
          print *, "============================================"
        else
          calculate_metric => metric_ef_phys
          calculate_christoffel_symbols => christoffel_ef_phys
          calculate_metric_derivatives => metric_derivs_ef_phys
          print *, "============================================"
          print *, " METRIC ASSIGNED: EF (Physical Spherical)   "
          print *, "============================================"
        end if
        
      case default
        print *, "CRITICAL ERROR: Unrecognized metric_type in metric assignment."
        stop
        
    end select
  end subroutine set_metric_type


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! ESPACIO-TIEMPO DE MINKOWSKI (1D/2D CARTESIANO PLANO)
  ! Geometría euclidiana trivial sin curvatura ni términos fuente geométricos.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine metric_minkowski(x_pos, y_pos, alpha, beta, gamma, gmunu, det, dlnalpha)
    implicit none
    real*8, intent(in) :: x_pos, y_pos
    real*8, intent(out), optional :: alpha, beta(3), gamma(3), det
    real*8, intent(out), optional :: gmunu(0:3, 0:3), dlnalpha(0:3)
    real*8 :: dummy
    dummy = x_pos; dummy = y_pos

    if (present(alpha)) alpha = 1.0d0
    if (present(beta))  beta  = 0.0d0
    if (present(gamma)) gamma = 1.0d0
    if (present(det))   det   = 1.0d0
    if (present(dlnalpha)) dlnalpha = 0.0d0

    if (present(gmunu)) then
      gmunu = 0.0d0
      gmunu(0,0) = -1.0d0
      gmunu(1,1) =  1.0d0
      gmunu(2,2) =  1.0d0
      gmunu(3,3) =  1.0d0
    end if
  end subroutine metric_minkowski

  subroutine christoffel_minkowski(x_pos, y_pos, chris)
    implicit none
    real*8, intent(in) :: x_pos, y_pos
    real*8, intent(out) :: chris(0:3,0:3,0:3)
    real*8 :: dummy
    dummy = x_pos; dummy = y_pos
    chris = 0.0d0 ! Cero absoluto: Sin curvatura, no hay aceleración geodésica.
  end subroutine christoffel_minkowski

  subroutine metric_derivs_minkowski(x_pos, y_pos, dg)
    implicit none
    real*8, intent(in) :: x_pos, y_pos
    real*8, intent(out) :: dg(0:3, 0:3, 1:3)
    real*8 :: dummy
    dummy = x_pos; dummy = y_pos
    dg = 0.0d0
  end subroutine metric_derivs_minkowski


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! ESPACIO-TIEMPO DE MINKOWSKI CILÍNDRICO AXISIMÉTRICO (r, z, phi)
  ! Mapeo Lógico: x1 = r (rho), x2 = z, x3 = phi
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subroutine metric_cylindrical(x_pos, y_pos, alpha, beta, gamma, gmunu, det, dlnalpha)
    implicit none
    real*8, intent(in) :: x_pos, y_pos
    real*8, intent(out), optional :: alpha, beta(3), gamma(3), det
    real*8, intent(out), optional :: gmunu(0:3, 0:3), dlnalpha(0:3)
    real*8 :: dummy
    dummy = y_pos

    if (present(alpha)) alpha = 1.0d0
    if (present(beta))  beta  = 0.0d0
    
    ! Mapeo Estandar: x1 = r, x2 = z, x3 = phi
    if (present(gamma)) then
      gamma(1) = 1.0d0     ! g_rr = 1
      gamma(2) = 1.0d0     ! g_zz = 1
      gamma(3) = x_pos**2  ! g_phiphi = r^2
    end if
    
    if (present(det)) det = x_pos**2
    
    if (present(dlnalpha)) dlnalpha = 0.0d0

    if (present(gmunu)) then
      gmunu = 0.0d0
      gmunu(0,0) = -1.0d0
      gmunu(1,1) =  1.0d0
      gmunu(2,2) =  1.0d0
      gmunu(3,3) =  1.0d0 / max(x_pos**2, 1.0d-20) 
    end if
  end subroutine metric_cylindrical

  subroutine christoffel_cylindrical(x_pos, y_pos, chris)
    implicit none
    real*8, intent(in) :: x_pos, y_pos
    real*8, intent(out) :: chris(0:3,0:3,0:3)
    real*8 :: dummy
    dummy = y_pos
    
    chris = 0.0d0 
    
    ! Símbolos de Christoffel (x1=r, x2=z, x3=phi)
    ! Gamma^r_{phi phi} = -r
    chris(1, 3, 3) = -x_pos
    
    chris(3, 1, 3) = 1.0d0 / max(x_pos, 1.0d-15)
    chris(3, 3, 1) = 1.0d0 / max(x_pos, 1.0d-15)
  end subroutine christoffel_cylindrical

  subroutine metric_derivs_cylindrical(x_pos, y_pos, dg)
    implicit none
    real*8, intent(in) :: x_pos, y_pos
    real*8, intent(out) :: dg(0:3, 0:3, 1:3)
    real*8 :: dummy
    
    dummy = y_pos
    dg = 0.0d0 
    
    ! d/dr (g_phiphi) = d/dr (r^2) = 2r
    ! Inyecta la fuerza centrífuga en S_phi
    dg(3,3,1) = 2.0d0 * x_pos
  end subroutine metric_derivs_cylindrical


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! ESPACIO-TIEMPO DE SCHWARZSCHILD (EDDINGTON-FINKELSTEIN) 2D
  ! Describe un Agujero Negro no rotante con coordenadas penetrantes.
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! ==============================================================================
  ! MÉTRICA FÍSICA ORIGINAL (r, theta, phi)
  ! ==============================================================================
  subroutine metric_ef_phys(r_phys, y_pos, alpha, beta, gamma, gmunu, det, dlnalpha)
    implicit none
    real*8, intent(in) :: r_phys, y_pos
    real*8, intent(out), optional :: alpha, beta(3), gamma(3), det
    real*8, intent(out), optional :: gmunu(0:3, 0:3), dlnalpha(0:3)

    real*8 :: g11, g22, g33, temp_alpha, temp_det
    real*8 :: temp_beta(3), temp_gamma(3)
    real*8 :: temp_gmunu(0:3, 0:3), temp_dlna(0:3)
    real*8 :: sin_th

    sin_th = sin(y_pos)

    temp_alpha     = 1.0d0 / sqrt(1.0d0 + 2.0d0 * bh_mass / r_phys)
    
    temp_beta(:)   = 0.0d0
    temp_beta(1)   = 2.0d0 * bh_mass / r_phys / (1.0d0 + 2.0d0 * bh_mass / r_phys)

    g11 = 1.0d0 + 2.0d0 * bh_mass / r_phys
    g22 = r_phys**2
    g33 = (r_phys**2) * (sin_th**2)

    temp_det = g11 * g22 * g33
    temp_gamma = (/ g11, g22, g33 /)

    temp_gmunu(:,:) = 0.0d0
    temp_gmunu(0,0) = -1.0d0 / (temp_alpha**2)
    temp_gmunu(0,1:3) = temp_beta / temp_alpha**2
    temp_gmunu(1:3,0) = temp_gmunu(0,1:3)

    temp_gmunu(1,1) = 1.0d0/g11 - (temp_beta(1)**2) / (temp_alpha**2)
    temp_gmunu(2,2) = 1.0d0/g22 - (temp_beta(2)**2) / (temp_alpha**2)
    temp_gmunu(3,3) = 1.0d0/g33 - (temp_beta(3)**2) / (temp_alpha**2)

    temp_dlna = 0.0d0
    temp_dlna(1) = bh_mass / (r_phys**2 + 2.0d0 * bh_mass * r_phys) 

    if (present(alpha))    alpha    = temp_alpha
    if (present(beta))     beta     = temp_beta
    if (present(gamma))    gamma    = temp_gamma
    if (present(gmunu))    gmunu    = temp_gmunu
    if (present(det))      det      = temp_det
    if (present(dlnalpha)) dlnalpha = temp_dlna
  end subroutine metric_ef_phys

  subroutine christoffel_ef_phys(r_phys, y_pos, chris)
    implicit none
    real*8, intent(in) :: r_phys, y_pos
    real*8, intent(out) :: chris(0:3,0:3,0:3) 
    real*8 :: sin_th, cos_th, safe_sin
    
    chris = 0.0d0 
    sin_th = sin(y_pos)
    cos_th = cos(y_pos)
    
    ! Blindaje polar preservando el signo para celdas fantasma que cruzan el eje
    safe_sin = sign(max(abs(sin_th), 1.0d-20), sin_th)
    
    chris(0,0,0) = 2.0d0 * bh_mass**2 / r_phys**3
    chris(0,0,1) = bh_mass * (1.0d0 + 2.0d0 * bh_mass / r_phys) / r_phys**2
    chris(0,1,0) = chris(0,0,1) 
    chris(0,1,1) = 2.0d0 * bh_mass * (1.0d0 + bh_mass / r_phys) / r_phys**2
    chris(0,2,2) = -2.0d0 * bh_mass
    chris(0,3,3) = chris(0,2,2) * (sin_th**2)
    
    chris(1,0,0) = bh_mass / r_phys**2 * (1.0d0 - 2.0d0 * bh_mass / r_phys)
    chris(1,0,1) = -chris(0,0,0) 
    chris(1,1,0) = chris(1,0,1)
    chris(1,1,1) = -chris(0,0,1)
    chris(1,2,2) = 2.0d0 * bh_mass - r_phys
    chris(1,3,3) = chris(1,2,2) * (sin_th**2)
    
    chris(2,1,2) = 1.0d0 / r_phys
    chris(2,2,1) = chris(2,1,2)
    chris(2,3,3) = -sin_th * cos_th
    
    chris(3,1,3) = 1.0d0 / r_phys
    chris(3,3,1) = chris(3,1,3)
    chris(3,2,3) = cos_th / safe_sin
    chris(3,3,2) = chris(3,2,3)
  end subroutine christoffel_ef_phys

  subroutine metric_derivs_ef_phys(r_phys, y_pos, dg)
    implicit none
    real*8, intent(in) :: r_phys, y_pos
    real*8, intent(out) :: dg(0:3, 0:3, 1:3) 
    real*8 :: sin_th, cos_th
    
    dg = 0.0d0 
    sin_th = sin(y_pos)
    cos_th = cos(y_pos)
    
    dg(0,0,1) = -2.0d0 * bh_mass / r_phys**2
    dg(0,1,1) = -2.0d0 * bh_mass / r_phys**2
    dg(1,0,1) = dg(0,1,1)
    dg(1,1,1) = -2.0d0 * bh_mass / r_phys**2
    dg(2,2,1) = 2.0d0 * r_phys
    dg(3,3,1) = 2.0d0 * r_phys * (sin_th**2)
    
    dg(3,3,2) = r_phys**2 * 2.0d0 * sin_th * cos_th
    ! SANEAMIENTO 1D: Matar el ruido de punto flotante en el ecuador
    if (ny == 1) dg(:,:,2) = 0.0d0 
  end subroutine metric_derivs_ef_phys

  subroutine metric_ef_log(x_pos, y_pos, alpha, beta, gamma, gmunu, det, dlnalpha)
    implicit none
    real*8, intent(in) :: x_pos, y_pos
    real*8, intent(out), optional :: alpha, beta(3), gamma(3), det
    real*8, intent(out), optional :: gmunu(0:3, 0:3), dlnalpha(0:3)

    real*8 :: g11, g22, g33, temp_alpha, temp_det
    real*8 :: temp_beta(3), temp_gamma(3)
    real*8 :: temp_gmunu(0:3, 0:3), temp_dlna(0:3)
    real*8 :: sin_th, r_phys

    r_phys = exp(x_pos)
    sin_th = sin(y_pos)

    ! Función de lapso (Lapse function - Escalar, no cambia)
    temp_alpha     = 1.0d0 / sqrt(1.0d0 + 2.0d0 * bh_mass / r_phys)
    
    ! Vector de desplazamiento LÓGICO (beta^x = beta^r / r)
    temp_beta(:)   = 0.0d0
    temp_beta(1)   = (2.0d0 * bh_mass / r_phys / (1.0d0 + 2.0d0 * bh_mass / r_phys)) / r_phys

    ! Componentes espaciales LÓGICAS (gamma_xx = gamma_rr * r^2)
    g11 = (1.0d0 + 2.0d0 * bh_mass / r_phys) * (r_phys**2)
    g22 = r_phys**2
    g33 = (r_phys**2) * (sin_th**2)

    temp_det = g11 * g22 * g33
    temp_gamma = (/ g11, g22, g33 /)

    ! Métrica espacio-temporal completa 4x4 (g^mu^nu contravariante)
    ! Al usar beta^x y gamma_xx, la relación algebraica estándar se mantiene intacta.
    temp_gmunu(:,:) = 0.0d0
    temp_gmunu(0,0) = -1.0d0 / (temp_alpha**2)
    temp_gmunu(0,1:3) = temp_beta / temp_alpha**2
    temp_gmunu(1:3,0) = temp_gmunu(0,1:3)

    temp_gmunu(1,1) = 1.0d0/g11 - (temp_beta(1)**2) / (temp_alpha**2)
    temp_gmunu(2,2) = 1.0d0/g22 - (temp_beta(2)**2) / (temp_alpha**2)
    temp_gmunu(3,3) = 1.0d0/g33 - (temp_beta(3)**2) / (temp_alpha**2)

    ! Derivada logarítmica del lapso (d_x ln(alpha) = r * d_r ln(alpha))
    temp_dlna = 0.0d0
    temp_dlna(1) = bh_mass / (r_phys + 2.0d0 * bh_mass) 

    if (present(alpha))    alpha    = temp_alpha
    if (present(beta))     beta     = temp_beta
    if (present(gamma))    gamma    = temp_gamma
    if (present(gmunu))    gmunu    = temp_gmunu
    if (present(det))      det      = temp_det
    if (present(dlnalpha)) dlnalpha = temp_dlna
  end subroutine metric_ef_log

  subroutine christoffel_ef_log(x_pos, y_pos, chris)
    implicit none
    real*8, intent(in) :: x_pos, y_pos
    real*8, intent(out) :: chris(0:3,0:3,0:3) 
    real*8 :: sin_th, cos_th, r_phys, safe_sin
    
    r_phys = exp(x_pos)
    chris = 0.0d0 
    sin_th = sin(y_pos)
    cos_th = cos(y_pos)
    
    ! Blindaje polar preservando el signo para celdas fantasma que cruzan el eje
    safe_sin = sign(max(abs(sin_th), 1.0d-20), sin_th)
    
    ! Símbolos de Christoffel para EF en coordenadas lógicas (x = ln r)
    chris(0,0,0) = 2.0d0 * bh_mass**2 / r_phys**3
    chris(0,0,1) = bh_mass * (1.0d0 + 2.0d0 * bh_mass / r_phys) / r_phys
    chris(0,1,0) = chris(0,0,1) 
    chris(0,1,1) = 2.0d0 * bh_mass * (1.0d0 + bh_mass / r_phys)
    chris(0,2,2) = -2.0d0 * bh_mass
    chris(0,3,3) = chris(0,2,2) * (sin_th**2)
    
    chris(1,0,0) = bh_mass / r_phys**3 * (1.0d0 - 2.0d0 * bh_mass / r_phys)
    chris(1,0,1) = -2.0d0 * bh_mass**2 / r_phys**3
    chris(1,1,0) = chris(1,0,1)
    
    ! Incluye el término no inercial d^2r/dx^2 de la transformación de coordenadas
    chris(1,1,1) = 1.0d0 - bh_mass / r_phys - 2.0d0 * bh_mass**2 / r_phys**2
    chris(1,2,2) = 2.0d0 * bh_mass / r_phys - 1.0d0
    chris(1,3,3) = chris(1,2,2) * (sin_th**2)
    
    ! Las componentes transversales se simplifican enormemente
    chris(2,1,2) = 1.0d0
    chris(2,2,1) = chris(2,1,2)
    chris(2,3,3) = -sin_th * cos_th
    
    chris(3,1,3) = 1.0d0
    chris(3,3,1) = chris(3,1,3)
    chris(3,2,3) = cos_th / safe_sin
    chris(3,3,2) = chris(3,2,3)
  end subroutine christoffel_ef_log

  subroutine metric_derivs_ef_log(x_pos, y_pos, dg)
    implicit none
    real*8, intent(in) :: x_pos, y_pos
    real*8, intent(out) :: dg(0:3, 0:3, 1:3) 
    real*8 :: sin_th, cos_th, r_phys
    
    dg = 0.0d0 
    r_phys = exp(x_pos)
    sin_th = sin(y_pos)
    cos_th = cos(y_pos)
    
    ! Derivadas respecto a x (dirección 1). Recordar d_x = r * d_r
    ! g_00_log = -1 + 2M/r
    dg(0,0,1) = -2.0d0 * bh_mass / r_phys
    
    ! g_01_log = r * g_01_phys = r * (2M/r) = 2M (¡Es constante en x!)
    dg(0,1,1) = 0.0d0
    dg(1,0,1) = dg(0,1,1)
    
    ! g_11_log = r^2 * g_11_phys = r^2 * (1 + 2M/r) = r^2 + 2Mr
    dg(1,1,1) = 2.0d0 * r_phys**2 + 2.0d0 * bh_mass * r_phys
    
    ! g_22_log = r^2
    dg(2,2,1) = 2.0d0 * r_phys**2
    
    ! g_33_log = r^2 * sin^2(th)
    dg(3,3,1) = 2.0d0 * (r_phys**2) * (sin_th**2)
    
    ! Derivadas respecto a theta (dirección 2)
    dg(3,3,2) = 2.0d0 * (r_phys**2) * sin_th * cos_th
    ! SANEAMIENTO 1D: Matar el ruido de punto flotante en el ecuador
    if (ny == 1) dg(:,:,2) = 0.0d0
  end subroutine metric_derivs_ef_log

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! ESPACIO-TIEMPO DE KERR (BOYER-LINDQUIST)
  ! Describe un Agujero Negro rotante. Depende de r (x_pos) y theta (y_pos).
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! subroutine metric_kerr(x_pos, y_pos, alpha, beta, gamma, gmunu, det, dlnalpha)
  !   implicit none
  !   real*8, intent(in) :: x_pos    ! Coordenada radial (r)
  !   real*8, intent(in) :: y_pos    ! Coordenada polar (theta)
  !   real*8, intent(out), optional :: alpha, beta(3), gamma(3), det
  !   real*8, intent(out), optional :: gmunu(0:3, 0:3), dlnalpha(0:3)

  !   real*8 :: a_spin, r, th, sin_th, cos_th
  !   real*8 :: Sigma, Delta, A_func
  !   real*8 :: g00, g03, g11, g22, g33
  !   real*8 :: temp_alpha, temp_det
  !   real*8 :: temp_beta(3), temp_gamma(3)
  !   real*8 :: temp_gmunu(0:3, 0:3), temp_dlna(0:3)

  !   ! Renombramos variables para claridad con la física
  !   r = x_pos
  !   th = y_pos
  !   a_spin = bh_spin ! Asegúrate de tener esta variable global (espín a) declarada en tu módulo
    
  !   sin_th = sin(th)
  !   cos_th = cos(th)

  !   ! 1. Funciones fundamentales de Kerr
  !   Sigma  = r**2 + (a_spin * cos_th)**2
  !   Delta  = r**2 - 2.0d0 * bh_mass * r + a_spin**2
  !   A_func = (r**2 + a_spin**2)**2 - a_spin**2 * Delta * (sin_th**2)

  !   ! 2. Componentes covariantes de la métrica (Boyer-Lindquist)
  !   g00 = -(1.0d0 - 2.0d0 * bh_mass * r / Sigma)
  !   g03 = -2.0d0 * bh_mass * r * a_spin * (sin_th**2) / Sigma
  !   g11 = Sigma / Delta
  !   g22 = Sigma
  !   ! Usamos la forma simplificada y optimizada para g_phiphi usando A_func
  !   g33 = A_func * (sin_th**2) / Sigma

  !   ! 3. Componentes espaciales puras (Formalismo 3+1)
  !   temp_gamma = 0.0d0
  !   temp_gamma(1) = g11
  !   temp_gamma(2) = g22
  !   temp_gamma(3) = g33

  !   temp_det = g11 * g22 * g33

  !   ! 4. Vector de desplazamiento (Shift vector - Arrastre de marco)
  !   temp_beta = 0.0d0
  !   temp_beta(3) = g03 / g33  ! beta^phi

  !   ! 5. Función de lapso explícita (Lapse function)
  !   ! Evitamos la resta numérica usando la solución analítica exacta
  !   if (A_func > 0.0d0 .and. Delta > 0.0d0) then
  !     temp_alpha = sqrt((Sigma * Delta) / A_func)
  !   else
  !     ! Protección numérica en caso de que la malla toque el horizonte coordenado
  !     ! donde Delta = 0, para evitar raíces imaginarias o divisiones por cero.
  !     temp_alpha = 1.0d-10 
  !   end if

  !   ! 6. Métrica espacio-temporal completa 4x4 (g_mu_nu covariante)
  !   temp_gmunu = 0.0d0
  !   temp_gmunu(0,0) = g00
  !   temp_gmunu(0,3) = g03
  !   temp_gmunu(3,0) = g03
  !   temp_gmunu(1,1) = g11
  !   temp_gmunu(2,2) = g22
  !   temp_gmunu(3,3) = g33

  !   ! 7. Derivada logarítmica del lapso
  !   ! NOTA: Ahora que es 2D, temp_dlna(1) es d(ln_a)/dr y temp_dlna(2) es d(ln_a)/dtheta.
  !   ! Puedes dejarlo en 0.0d0 temporalmente y calcular los gradientes numéricamente 
  !   ! en tu solver, o usar Maple para exportar la derivada analítica más adelante.
  !   temp_dlna = 0.0d0 

  !   ! Asignación final a variables de salida opcionales
  !   if (present(alpha))    alpha    = temp_alpha
  !   if (present(beta))     beta     = temp_beta
  !   if (present(gamma))    gamma    = temp_gamma
  !   if (present(gmunu))    gmunu    = temp_gmunu
  !   if (present(det))      det      = temp_det
  !   if (present(dlnalpha)) dlnalpha = temp_dlna
  ! end subroutine metric_kerr


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TENSOR DE ENERGÍA-MOMENTO (T_mu_nu)
  ! Computa el tensor de esfuerzo para un fluido perfecto relativista.
  ! Es completamente agnóstico de la métrica (funciona para todas).
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine calculate_stress_energy_tensor(prim_state, i, j, k, Tmunu)
    real*8, intent(in) :: prim_state(:)
    integer, intent(in) :: i, j, k
    real*8, intent(out) :: Tmunu(0:3,0:3)
    
    real*8 :: alpha, beta(3), g(3), gmunu(0:3,0:3)
    real*8 :: v_sq, h, W, u_four(0:3)
    real*8 :: v1, v2, v3
    integer :: mu, nu

    ! 1. LEER DE LA CACHÉ GLOBAL (Centros)
    alpha      = alpha_c(i,j,k)
    beta(:)    = beta_c(:,i,j,k)
    g(:)       = gamma_c(:,i,j,k)
    gmunu(:,:) = gmunu_c(:,:,i,j,k)

    v1 = prim_state(eq_vx)
    v2 = prim_state(eq_vy)
    v3 = prim_state(eq_vz)

    ! Velocidad cuadrada
    v_sq = v1**2 * g(1) + v2**2 * g(2) + v3**2 * g(3)
    
    ! Limitador relativista de seguridad (previene superar c)
    if (v_sq >= v_max) then
      v1 = v1 * sqrt(v_max / v_sq)
      v2 = v2 * sqrt(v_max / v_sq)
      v3 = v3 * sqrt(v_max / v_sq)
      v_sq = v_max
    end if

    ! Entalpía específica (h) y Factor de Lorentz (W)
    h = 1.0d0 + g1 * prim_state(eq_pr) / prim_state(eq_de)
    W = 1.0d0 / sqrt(1.0d0 - v_sq)

    ! Cuadrivector de velocidad u^mu
    u_four = 0.0d0
    u_four(0) = W / alpha
    u_four(1) = W * (v1 - beta(1) / alpha) 
    u_four(2) = W * (v2 - beta(2) / alpha) 
    u_four(3) = W * (v3 - beta(3) / alpha) 

    ! Ensamblaje del Tensor de Energía-Momento contravariante
    Tmunu = 0.0d0 !T^mu^nu = rho * h * u^mu * u^nu + p * g^mu^nu
    do mu = 0, 3
      do nu = mu, 3
        Tmunu(mu,nu) = (prim_state(eq_de) * h * u_four(mu) * u_four(nu)) + (prim_state(eq_pr) * gmunu(mu,nu))
        Tmunu(nu,mu) = Tmunu(mu,nu) ! Espejear por simetría
      end do
    end do

  end subroutine calculate_stress_energy_tensor

end module metrics