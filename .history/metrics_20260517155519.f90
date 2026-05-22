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
      real*8, intent(out), optional :: alpha, beta(3), gamma(3,3), det
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

      
      ! ---------------------------------------------------------
      ! AGUJEROS NEGROS CON ROTACIÓN
      ! ---------------------------------------------------------
      case ('Kerr-Schild')
        
        if (trim(geom_type) /= 'Spherical') then
          print *, "CRITICAL ERROR: KS Metric strictly requires Spherical geometry."
          stop
        end if
        
        if (use_log_r) then
          calculate_metric => metric_ks_log
          calculate_christoffel_symbols => christoffel_ks_log
          calculate_metric_derivatives => metric_derivs_ks_log
          print *, "============================================"
          print *, " METRIC ASSIGNED: KS (Logarithmic Spherical)"
          print *, "============================================"
        else
          calculate_metric => metric_ks_phys
          calculate_christoffel_symbols => christoffel_ks_phys
          calculate_metric_derivatives => metric_derivs_ks_phys
          print *, "============================================"
          print *, " METRIC ASSIGNED: KS (Physical Spherical)   "
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
    real*8, intent(out), optional :: alpha, beta(3), gamma(3,3), det
    real*8, intent(out), optional :: gmunu(0:3, 0:3), dlnalpha(0:3)
    real*8 :: dummy
    
    dummy = x_pos; dummy = y_pos

    if (present(alpha)) alpha = 1.0d0
    if (present(beta))  beta  = 0.0d0
    
    ! Matriz espacial diagonal (Identidad)
    if (present(gamma)) then
      gamma = 0.0d0
      gamma(1,1) = 1.0d0
      gamma(2,2) = 1.0d0
      gamma(3,3) = 1.0d0
    end if
    
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
    real*8, intent(out), optional :: alpha, beta(3), gamma(3,3), det
    real*8, intent(out), optional :: gmunu(0:3, 0:3), dlnalpha(0:3)
    real*8 :: dummy
    
    dummy = y_pos

    if (present(alpha)) alpha = 1.0d0
    if (present(beta))  beta  = 0.0d0
    
    ! Mapeo Estándar: x1 = r, x2 = z, x3 = phi
    if (present(gamma)) then
      gamma = 0.0d0
      gamma(1,1) = 1.0d0     ! g_rr = 1
      gamma(2,2) = 1.0d0     ! g_zz = 1
      gamma(3,3) = x_pos**2  ! g_phiphi = r^2
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
    real*8, intent(out), optional :: alpha, beta(3), gamma(3,3), det
    real*8, intent(out), optional :: gmunu(0:3, 0:3), dlnalpha(0:3)

    real*8 :: g11, g22, g33, temp_alpha, temp_det
    real*8 :: temp_beta(3), temp_gamma(3,3)
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
    
    ! Asignación diagonal para la métrica espacial
    temp_gamma = 0.0d0
    temp_gamma(1,1) = g11
    temp_gamma(2,2) = g22
    temp_gamma(3,3) = g33

    temp_gmunu(:,:) = 0.0d0
    temp_gmunu(0,0) = -1.0d0 / (temp_alpha**2)
    temp_gmunu(0,1:3) = temp_beta / temp_alpha**2
    temp_gmunu(1:3,0) = temp_gmunu(0,1:3)

    temp_gmunu(1,1) = 1.0d0/g11 - (temp_beta(1)**2) / (temp_alpha**2)
    temp_gmunu(2,2) = 1.0d0/g22 - (temp_beta(2)**2) / (temp_alpha**2)
    temp_gmunu(3,3) = 1.0d0/g33 - (temp_beta(3)**2) / (temp_alpha**2)

    if (abs(gamma(1,2)) > 1.0d-14 .or. abs(gamma(1,3)) > 1.0d-14) then
      print *, "¡ALERTA Métrica! Términos cruzados espaciales NO son cero en EF."
      print *, "gamma(1,2) =", gamma(1,2)
      stop
    end if

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
    ! gamma actualizado a matriz 3x3
    real*8, intent(out), optional :: alpha, beta(3), gamma(3,3), det
    real*8, intent(out), optional :: gmunu(0:3, 0:3), dlnalpha(0:3)

    real*8 :: g11, g22, g33, temp_alpha, temp_det
    real*8 :: temp_beta(3), temp_gamma(3,3)
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
    
    ! Asignación diagonal para la métrica espacial lógica
    temp_gamma = 0.0d0
    temp_gamma(1,1) = g11
    temp_gamma(2,2) = g22
    temp_gamma(3,3) = g33

    ! Métrica espacio-temporal completa 4x4 (g^mu^nu contravariante)
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
  ! ESPACIO-TIEMPO DE KERR-SCHILD (KS) 3D ENTRANTE
  ! Métrica Física Pura (r_phys, theta, phi)
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine metric_ks_phys(r_phys, y_pos, alpha, beta, gamma, gmunu, det, dlnalpha)
    implicit none
    real*8, intent(in) :: r_phys, y_pos
    real*8, intent(out), optional :: alpha, beta(3), gamma(3,3), det
    real*8, intent(out), optional :: gmunu(0:3, 0:3), dlnalpha(0:3)

    real*8 :: a, M, Sigma, H, sin_th, cos_th, sin2, cos2, Delta
    real*8 :: temp_alpha, temp_det
    real*8 :: temp_beta(3), temp_gamma(3,3)
    real*8 :: temp_gmunu(0:3, 0:3), temp_dlna(0:3)

    M = bh_mass
    a = a_spin ! Asegúrate de que a_spin esté definida en tu módulo global

    sin_th = sin(y_pos)
    cos_th = cos(y_pos)
    sin2 = sin_th**2
    cos2 = cos_th**2

    Sigma = r_phys**2 + (a**2) * cos2
    H = M * r_phys / Sigma
    Delta = r_phys**2 - 2.0d0 * M * r_phys + a**2

    ! 1. Factor de Lapso (alpha)
    ! Obtenido de g_tt = -(1 - 2H) = -alpha^2 + beta_i beta^i
    temp_alpha = 1.0d0 / sqrt(1.0d0 + 2.0d0 * H)

    ! 2. Vector de Desplazamiento contravariante (beta^i)
    temp_beta(:) = 0.0d0
    temp_beta(1) = 2.0d0 * H / (1.0d0 + 2.0d0 * H)
    ! beta^theta = 0.0d0
    ! beta^phi   = 0.0d0 

    ! 3. Métrica espacial covariante (gamma_ij)
    ! Es exactamente el bloque 3x3 inferior derecho de tu matriz de Maple
    temp_gamma = 0.0d0
    temp_gamma(1,1) = 1.0d0 + 2.0d0 * H
    temp_gamma(1,3) = -a * (1.0d0 + 2.0d0 * H) * sin2
    temp_gamma(3,1) = temp_gamma(1,3)
    temp_gamma(2,2) = Sigma
    temp_gamma(3,3) = (r_phys**2 + a**2 + 2.0d0 * H * a**2 * sin2) * sin2

    ! Determinante espacial de gamma_ij
    temp_det = (Sigma**2) * sin2 * (1.0d0 + 2.0d0 * H)

    ! 4. Métrica espacio-temporal contravariante (g^munu)
    temp_gmunu = 0.0d0
    temp_gmunu(0,0) = -(1.0d0 + 2.0d0 * H)
    temp_gmunu(0,1) = 2.0d0 * H
    temp_gmunu(1,0) = temp_gmunu(0,1)
    
    temp_gmunu(1,1) = Delta / Sigma
    temp_gmunu(1,3) = a / Sigma
    temp_gmunu(3,1) = temp_gmunu(1,3)
    
    temp_gmunu(2,2) = 1.0d0 / Sigma
    
    ! Blindaje polar estricto para g^phiphi para evitar divisiones por cero en el eje
    if (abs(sin_th) > 1.0d-15) then
        temp_gmunu(3,3) = 1.0d0 / (Sigma * sin2)
    else
        temp_gmunu(3,3) = 1.0d0 / (Sigma * 1.0d-30)
    end if

    ! 5. d(ln alpha) / dx^i
    ! Usando la regla de la cadena sobre alpha = (1 + 2H)^(-1/2)
    temp_dlna = 0.0d0
    block
        real*8 :: dH_dr, dH_dth
        dH_dr = M * (a**2 * cos2 - r_phys**2) / (Sigma**2)
        dH_dth = 2.0d0 * M * r_phys * a**2 * sin_th * cos_th / (Sigma**2)
        
        temp_dlna(1) = - dH_dr / (1.0d0 + 2.0d0 * H)
        temp_dlna(2) = - dH_dth / (1.0d0 + 2.0d0 * H)
    end block

    ! 6. Mapeo final a las variables de salida (Opcionales)
    if (present(alpha))    alpha    = temp_alpha
    if (present(beta))     beta     = temp_beta
    if (present(gamma))    gamma    = temp_gamma
    if (present(gmunu))    gmunu    = temp_gmunu
    if (present(det))      det      = temp_det
    if (present(dlnalpha)) dlnalpha = temp_dlna

  end subroutine metric_ks_phys

  ! ==============================================================================
  ! DERIVADAS DE LA MÉTRICA FÍSICA (KS 3D ENTRANTE)
  ! dg(mu, nu, k) donde k=1 es r, k=2 es theta
  ! ==============================================================================
  subroutine metric_derivs_ks_phys(r_phys, y_pos, dg)
    implicit none
    real*8, intent(in) :: r_phys, y_pos
    real*8, intent(out) :: dg(0:3, 0:3, 1:3) 
    
    real*8 :: a, M, Sigma, H, sin_th, cos_th, sin2, cos2
    real*8 :: dS_dr, dS_dth, dH_dr, dH_dth
    
    M = bh_mass
    a = a_spin
    
    sin_th = sin(y_pos)
    cos_th = cos(y_pos)
    sin2 = sin_th**2
    cos2 = cos_th**2
    
    Sigma = r_phys**2 + a**2 * cos2
    H = M * r_phys / Sigma
    
    ! Derivadas auxiliares de Sigma y H
    dS_dr = 2.0d0 * r_phys
    dS_dth = -2.0d0 * a**2 * sin_th * cos_th
    
    dH_dr = M * (a**2 * cos2 - r_phys**2) / (Sigma**2)
    dH_dth = 2.0d0 * M * r_phys * a**2 * sin_th * cos_th / (Sigma**2)
    
    dg = 0.0d0
    
    ! ---------------------------------------------------------
    ! Derivadas respecto a R (k=1)
    ! ---------------------------------------------------------
    dg(0,0,1) = 2.0d0 * dH_dr
    
    dg(0,1,1) = 2.0d0 * dH_dr
    dg(1,0,1) = dg(0,1,1)
    
    dg(0,3,1) = -2.0d0 * a * sin2 * dH_dr
    dg(3,0,1) = dg(0,3,1)
    
    dg(1,1,1) = 2.0d0 * dH_dr
    
    dg(1,3,1) = -a * sin2 * (2.0d0 * dH_dr)
    dg(3,1,1) = dg(1,3,1)
    
    dg(2,2,1) = dS_dr
    dg(3,3,1) = 2.0d0 * r_phys * sin2 + 2.0d0 * a**2 * (sin2**2) * dH_dr
    
    ! ---------------------------------------------------------
    ! Derivadas respecto a THETA (k=2)
    ! ---------------------------------------------------------
    dg(0,0,2) = 2.0d0 * dH_dth
    
    dg(0,1,2) = 2.0d0 * dH_dth
    dg(1,0,2) = dg(0,1,2)
    
    dg(0,3,2) = -2.0d0 * a * (sin2 * dH_dth + 2.0d0 * H * sin_th * cos_th)
    dg(3,0,2) = dg(0,3,2)
    
    dg(1,1,2) = 2.0d0 * dH_dth
    
    dg(1,3,2) = -a * ( (1.0d0 + 2.0d0 * H) * 2.0d0 * sin_th * cos_th + sin2 * 2.0d0 * dH_dth )
    dg(3,1,2) = dg(1,3,2)
    
    dg(2,2,2) = dS_dth
    
    dg(3,3,2) = 2.0d0 * (r_phys**2 + a**2) * sin_th * cos_th + &
                2.0d0 * a**2 * (sin2**2) * dH_dth + &
                8.0d0 * a**2 * H * sin_th * cos_th * sin2
    
    ! =========================================================
    ! SANEAMIENTO ECUATORIAL (BLINDAJE NUMÉRICO)
    ! Si la malla es puramente ecuatorial (ny == 1), TODAS las 
    ! derivadas en theta deben ser matemáticamente cero. Esto 
    ! mata el residuo de cos(pi/2) ~ 1e-17.
    ! =========================================================
    if (ny == 1) dg(:,:,2) = 0.0d0 

  end subroutine metric_derivs_ks_phys

  ! ==============================================================================
  ! SÍMBOLOS DE CHRISTOFFEL FÍSICOS (Esquema Híbrido Optimizado)
  ! ==============================================================================
  subroutine christoffel_ks_phys(r_phys, y_pos, chris)
    implicit none
    real*8, intent(in) :: r_phys, y_pos
    real*8, intent(out) :: chris(0:3,0:3,0:3)
    
    real*8 :: gmunu_up(0:3, 0:3)
    real*8 :: dg(0:3, 0:3, 1:3)
    integer :: lambda, mu, nu, rho
    real*8 :: d_mu_g, d_nu_g, d_rho_g
    
    ! 1. Obtenemos el contravariante analítico (g^munu)
    call metric_ks_phys(r_phys, y_pos, gmunu=gmunu_up)
    
    ! 2. Obtenemos las derivadas analíticas exactas (d_k g_munu)
    call metric_derivs_ks_phys(r_phys, y_pos, dg)
    
    chris = 0.0d0
    
    ! 3. Contracción tensorial optimizada (Solo 40 cálculos en lugar de 64)
    ! Gamma^lambda_{mu nu} = 0.5 * g^{lambda rho} * ( d_mu g_{rho nu} + d_nu g_{rho mu} - d_rho g_{mu nu} )
    do lambda = 0, 3
      do mu = 0, 3
        ! OJO AQUÍ: El bucle 'nu' empieza en 'mu' para aprovechar la simetría
        do nu = mu, 3
          
          ! Sumatoria sobre el índice mudo 'rho'
          do rho = 0, 3
            d_mu_g = 0.0d0
            d_nu_g = 0.0d0
            d_rho_g = 0.0d0
            
            ! Solo existen derivadas en las direcciones espaciales r (1) y theta (2)
            if (mu == 1 .or. mu == 2) d_mu_g = dg(rho, nu, mu)
            if (nu == 1 .or. nu == 2) d_nu_g = dg(rho, mu, nu)
            if (rho == 1 .or. rho == 2) d_rho_g = dg(mu, nu, rho)
            
            chris(lambda, mu, nu) = chris(lambda, mu, nu) + &
              0.5d0 * gmunu_up(lambda, rho) * (d_mu_g + d_nu_g - d_rho_g)
          end do
          
          ! 4. Clonamos inmediatamente el resultado al índice simétrico
          if (mu /= nu) then
            chris(lambda, nu, mu) = chris(lambda, mu, nu)
          end if
          
        end do
      end do
    end do
  end subroutine christoffel_ks_phys


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! ESPACIO-TIEMPO DE KERR-SCHILD (KS) 3D ENTRANTE - MALLA LOGARÍTMICA
  ! Métrica Lógica (x_pos = ln(r), y_pos = theta, phi)
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine metric_ks_log(x_pos, y_pos, alpha, beta, gamma, gmunu, det, dlnalpha)
    implicit none
    real*8, intent(in) :: x_pos, y_pos
    ! OJO: gamma es matriz 3x3 por el término cruzado x-phi (gamma_13)
    real*8, intent(out), optional :: alpha, beta(3), gamma(3,3), det
    real*8, intent(out), optional :: gmunu(0:3, 0:3), dlnalpha(0:3)

    real*8 :: a, M, Sigma, H, sin_th, cos_th, sin2, cos2, Delta, r_phys
    real*8 :: temp_alpha, temp_det
    real*8 :: temp_beta(3), temp_gamma(3,3)
    real*8 :: temp_gmunu(0:3, 0:3), temp_dlna(0:3)

    M = bh_mass
    a = a_spin 
    r_phys = exp(x_pos)

    sin_th = sin(y_pos)
    cos_th = cos(y_pos)
    sin2 = sin_th**2
    cos2 = cos_th**2

    Sigma = r_phys**2 + (a**2) * cos2
    H = M * r_phys / Sigma
    Delta = r_phys**2 - 2.0d0 * M * r_phys + a**2

    ! 1. Factor de Lapso (alpha) - Invariante ante transformaciones espaciales
    temp_alpha = 1.0d0 / sqrt(1.0d0 + 2.0d0 * H)

    ! 2. Vector de Desplazamiento LÓGICO contravariante (beta^x = beta^r / r)
    temp_beta(:) = 0.0d0
    temp_beta(1) = (2.0d0 * H / (1.0d0 + 2.0d0 * H)) / r_phys

    ! 3. Métrica espacial LÓGICA covariante (gamma_ij_log)
    ! Aplicando las reglas de transformación: 
    ! gamma_xx = r^2 * gamma_rr  |  gamma_xphi = r * gamma_rphi
    temp_gamma = 0.0d0
    temp_gamma(1,1) = (1.0d0 + 2.0d0 * H) * r_phys**2
    temp_gamma(1,3) = -a * (1.0d0 + 2.0d0 * H) * sin2 * r_phys
    temp_gamma(3,1) = temp_gamma(1,3)
    temp_gamma(2,2) = Sigma
    temp_gamma(3,3) = (r_phys**2 + a**2 + 2.0d0 * H * a**2 * sin2) * sin2

    ! Determinante espacial LÓGICO (det_log = r^2 * det_phys)
    temp_det = (Sigma**2 * sin2 * (1.0d0 + 2.0d0 * H)) * r_phys**2

    ! 4. Métrica espacio-temporal completa LÓGICA contravariante (g^munu_log)
    ! Los componentes contravariantes se transforman con derivadas inversas (multiplicando por r)
    temp_gmunu = 0.0d0
    temp_gmunu(0,0) = -(1.0d0 + 2.0d0 * H)
    
    temp_gmunu(0,1) = 2.0d0 * H * r_phys  ! g^xt = r * g^rt
    temp_gmunu(1,0) = temp_gmunu(0,1)
    
    temp_gmunu(1,1) = (Delta / Sigma) * r_phys**2  ! g^xx = r^2 * g^rr
    
    temp_gmunu(1,3) = (a / Sigma) * r_phys  ! g^xphi = r * g^rphi
    temp_gmunu(3,1) = temp_gmunu(1,3)
    
    temp_gmunu(2,2) = 1.0d0 / Sigma
    
    ! Blindaje polar para evitar divisiones por cero en el eje
    if (abs(sin_th) > 1.0d-15) then
        temp_gmunu(3,3) = 1.0d0 / (Sigma * sin2)
    else
        temp_gmunu(3,3) = 1.0d0 / (Sigma * 1.0d-30)
    end if

    ! 5. d(ln alpha) / dx^i en la malla lógica (d_x = r * d_r)
    temp_dlna = 0.0d0
    block
        real*8 :: dH_dr, dH_dth
        dH_dr = M * (a**2 * cos2 - r_phys**2) / (Sigma**2)
        dH_dth = 2.0d0 * M * r_phys * a**2 * sin_th * cos_th / (Sigma**2)
        
        temp_dlna(1) = - (dH_dr * r_phys) / (1.0d0 + 2.0d0 * H)
        temp_dlna(2) = - dH_dth / (1.0d0 + 2.0d0 * H)
    end block

    if (present(alpha))    alpha    = temp_alpha
    if (present(beta))     beta     = temp_beta
    if (present(gamma))    gamma    = temp_gamma
    if (present(gmunu))    gmunu    = temp_gmunu
    if (present(det))      det      = temp_det
    if (present(dlnalpha)) dlnalpha = temp_dlna
  end subroutine metric_ks_log

  ! ==============================================================================
  ! DERIVADAS DE LA MÉTRICA LÓGICA (KS 3D ENTRANTE)
  ! dg_log(mu, nu, k) donde k=1 es x=ln(r), k=2 es theta
  ! ==============================================================================
  subroutine metric_derivs_ks_log(x_pos, y_pos, dg)
    implicit none
    real*8, intent(in) :: x_pos, y_pos
    real*8, intent(out) :: dg(0:3, 0:3, 1:3) 
    
    real*8 :: a, M, Sigma, H, sin_th, cos_th, sin2, cos2, r_phys
    real*8 :: dS_dth, dH_dr, dH_dth
    real*8 :: g_tr_phys, g_rr_phys, g_rphi_phys
    
    M = bh_mass
    a = a_spin
    r_phys = exp(x_pos)
    
    sin_th = sin(y_pos)
    cos_th = cos(y_pos)
    sin2 = sin_th**2
    cos2 = cos_th**2
    
    Sigma = r_phys**2 + a**2 * cos2
    H = M * r_phys / Sigma
    
    dS_dth = -2.0d0 * a**2 * sin_th * cos_th
    dH_dr = M * (a**2 * cos2 - r_phys**2) / (Sigma**2)
    dH_dth = 2.0d0 * M * r_phys * a**2 * sin_th * cos_th / (Sigma**2)
    
    ! Componentes físicas puras requeridas para derivar los productos (Regla de Leibniz)
    g_tr_phys = 2.0d0 * H
    g_rr_phys = 1.0d0 + 2.0d0 * H
    g_rphi_phys = -2.0d0 * a * sin2 * H
    
    dg = 0.0d0
    
    ! ---------------------------------------------------------
    ! Derivadas respecto a X (k=1). 
    ! Aplicando: d_x ( f(r) ) = r * d_r ( f(r) )
    ! ---------------------------------------------------------
    
    ! g_tt_log = g_tt_phys
    dg(0,0,1) = r_phys * (2.0d0 * dH_dr)
    
    ! d_x(r * g_tr_phys) = r * g_tr_phys + r^2 * d_r(g_tr_phys)
    dg(0,1,1) = r_phys * g_tr_phys + r_phys**2 * (2.0d0 * dH_dr)
    dg(1,0,1) = dg(0,1,1)
    
    ! g_tphi_log = g_tphi_phys
    dg(0,3,1) = r_phys * (-2.0d0 * a * sin2 * dH_dr)
    dg(3,0,1) = dg(0,3,1)
    
    ! d_x(r^2 * g_rr_phys) = 2r^2 * g_rr_phys + r^3 * d_r(g_rr_phys)
    dg(1,1,1) = 2.0d0 * r_phys**2 * g_rr_phys + r_phys**3 * (2.0d0 * dH_dr)
    
    ! d_x(r * g_rphi_phys) = r * g_rphi_phys + r^2 * d_r(g_rphi_phys)
    dg(1,3,1) = r_phys * g_rphi_phys + r_phys**2 * (-2.0d0 * a * sin2 * dH_dr)
    dg(3,1,1) = dg(1,3,1)
    
    ! d_x(Sigma) = r * d_r(Sigma) = r * (2r) = 2r^2
    dg(2,2,1) = 2.0d0 * r_phys**2
    
    ! g_phiphi_log = g_phiphi_phys
    dg(3,3,1) = r_phys * (2.0d0 * r_phys * sin2 + 2.0d0 * a**2 * sin2**2 * dH_dr)
    
    ! ---------------------------------------------------------
    ! Derivadas respecto a THETA (k=2)
    ! d_theta ( f(r, theta) ) de los componentes mapeados lógicamente
    ! ---------------------------------------------------------
    dg(0,0,2) = 2.0d0 * dH_dth
    
    dg(0,1,2) = r_phys * (2.0d0 * dH_dth)
    dg(1,0,2) = dg(0,1,2)
    
    dg(0,3,2) = -2.0d0 * a * (sin2 * dH_dth + 2.0d0 * H * sin_th * cos_th)
    dg(3,0,2) = dg(0,3,2)
    
    dg(1,1,2) = r_phys**2 * (2.0d0 * dH_dth)
    
    dg(1,3,2) = r_phys * ( -a * ( (1.0d0 + 2.0d0 * H) * 2.0d0 * sin_th * cos_th + sin2 * 2.0d0 * dH_dth ) )
    dg(3,1,2) = dg(1,3,2)
    
    dg(2,2,2) = dS_dth
    
    dg(3,3,2) = 2.0d0 * (r_phys**2 + a**2) * sin_th * cos_th + &
                2.0d0 * a**2 * (sin2**2) * dH_dth + &
                8.0d0 * a**2 * H * sin_th * cos_th * sin2
    
    ! =========================================================
    ! SANEAMIENTO ECUATORIAL (BLINDAJE NUMÉRICO)
    ! =========================================================
    if (ny == 1) dg(:,:,2) = 0.0d0 

  end subroutine metric_derivs_ks_log


  ! ==============================================================================
  ! SÍMBOLOS DE CHRISTOFFEL LÓGICOS (Esquema Híbrido Optimizado)
  ! ==============================================================================
  subroutine christoffel_ks_log(x_pos, y_pos, chris)
    implicit none
    real*8, intent(in) :: x_pos, y_pos
    real*8, intent(out) :: chris(0:3,0:3,0:3)
    
    real*8 :: gmunu_up(0:3, 0:3)
    real*8 :: dg(0:3, 0:3, 1:3)
    integer :: lambda, mu, nu, rho
    real*8 :: d_mu_g, d_nu_g, d_rho_g
    
    ! Extraemos los tensores lógicos
    call metric_ks_log(x_pos, y_pos, gmunu=gmunu_up)
    call metric_derivs_ks_log(x_pos, y_pos, dg)
    
    chris = 0.0d0
    
    ! Contracción tensorial (Cálculo idéntico, la geometría diferencial se encarga de la lógica)
    do lambda = 0, 3
      do mu = 0, 3
        do nu = mu, 3
          
          do rho = 0, 3
            d_mu_g = 0.0d0
            d_nu_g = 0.0d0
            d_rho_g = 0.0d0
            
            if (mu == 1 .or. mu == 2) d_mu_g = dg(rho, nu, mu)
            if (nu == 1 .or. nu == 2) d_nu_g = dg(rho, mu, nu)
            if (rho == 1 .or. rho == 2) d_rho_g = dg(mu, nu, rho)
            
            chris(lambda, mu, nu) = chris(lambda, mu, nu) + &
              0.5d0 * gmunu_up(lambda, rho) * (d_mu_g + d_nu_g - d_rho_g)
          end do
          
          if (mu /= nu) then
            chris(lambda, nu, mu) = chris(lambda, mu, nu)
          end if
          
        end do
      end do
    end do
  end subroutine christoffel_ks_log


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TENSOR DE ENERGÍA-MOMENTO (T_mu_nu)
  ! Computa el tensor de esfuerzo para un fluido perfecto relativista.
  ! Es completamente agnóstico de la métrica (funciona para todas).
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine calculate_stress_energy_tensor(prim_state, i, j, k, Tmunu)
    implicit none
    real*8, intent(in) :: prim_state(:)
    integer, intent(in) :: i, j, k
    real*8, intent(out) :: Tmunu(0:3,0:3)
    
    real*8 :: alpha, beta(3), g(3,3), gmunu(0:3,0:3)
    real*8 :: v_sq, h, W, u_four(0:3)
    real*8 :: v(3)
    integer :: mu, nu, ii, jj  ! ii, jj agregados para proteger i, j

    ! 1. LEER DE LA CACHÉ GLOBAL (Centros)
    alpha      = alpha_c(i,j,k)
    beta(:)    = beta_c(:,i,j,k)
    g(:,:)     = gamma_c(:,:,i,j,k)
    gmunu(:,:) = gmunu_c(:,:,i,j,k)

    v(1) = prim_state(eq_vx)
    v(2) = prim_state(eq_vy)
    v(3) = prim_state(eq_vz)

    ! 2. Velocidad cuadrada (Usando los índices locales ii, jj)
    v_sq = 0.0d0
    do ii = 1, 3
      do jj = 1, 3
        v_sq = v_sq + g(ii,jj) * v(ii) * v(jj)
      end do
    end do
    
    ! Limitador relativista de seguridad (previene superar c)
    if (v_sq >= v_max) then
      v(1) = v(1) * sqrt(v_max / v_sq)
      v(2) = v(2) * sqrt(v_max / v_sq)
      v(3) = v(3) * sqrt(v_max / v_sq)
      v_sq = v_max
    end if

    ! Entalpía específica (h) y Factor de Lorentz (W)
    h = 1.0d0 + g1 * prim_state(eq_pr) / prim_state(eq_de)
    W = 1.0d0 / sqrt(1.0d0 - v_sq)

    ! Cuadrivector de velocidad u^mu
    u_four = 0.0d0
    u_four(0) = W / alpha
    u_four(1) = W * (v(1) - beta(1) / alpha) 
    u_four(2) = W * (v(2) - beta(2) / alpha) 
    u_four(3) = W * (v(3) - beta(3) / alpha) 

    ! Ensamblaje del Tensor de Energía-Momento contravariante
    Tmunu = 0.0d0 ! T^mu^nu = rho * h * u^mu * u^nu + p * g^mu^nu
    do mu = 0, 3
      do nu = mu, 3
        Tmunu(mu,nu) = (prim_state(eq_de) * h * u_four(mu) * u_four(nu)) + (prim_state(eq_pr) * gmunu(mu,nu))
        if (mu /= nu) then
            Tmunu(nu,mu) = Tmunu(mu,nu) ! Espejear por simetría
        end if
      end do
    end do

  end subroutine calculate_stress_energy_tensor

end module metrics