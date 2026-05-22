! =======================================================================================
! Módulo: fluxes (Solucionadores de Riemann y Velocidades Características)
! ---------------------------------------------------------------------------------------
! Propósito: Resuelve el Problema de Riemann local en cada interfaz de la malla. Dada una 
! discontinuidad entre un estado izquierdo (q_L) y uno derecho (q_R) arrojados por el 
! reconstructor espacial, este módulo calcula el flujo numérico neto que cruza la pared.
!
! Arquitectura: Funciona como un "Hub" o director. Dependiendo de la variable global 
! riemann_solver_id, enruta el cálculo hacia esquemas ultra-robustos (HLLE) o esquemas 
! de altísima nitidez para la onda de contacto (HLLC).
! =======================================================================================
module fluxes
  use equations
  use variables
  use metrics
  implicit none
  private
  
  ! Hacemos público el enrutador de Riemann y el inicializador del puntero
  public :: resolve_riemann_problem, init_wavespeed_solver

  ! ====================================================================
  ! INTERFAZ ABSTRACTA Y DECLARACIÓN DEL PUNTERO
  ! ====================================================================
  abstract interface
    subroutine wavespeed_interface(prim_state, sweep_dir, alpha, beta, g, ,sqg, wave_min, wave_max)
      implicit none
      real*8, intent(in)  :: prim_state(:)
      integer, intent(in) :: sweep_dir
      real*8, intent(in)  :: alpha, beta(3), g(3,3)
      real*8, intent(out) :: wave_min, wave_max
    end subroutine wavespeed_interface
  end interface

  ! Este es el puntero que se comportará como una subrutina normal
  procedure(wavespeed_interface), pointer, public :: calc_wavespeeds => null()

contains

  ! ====================================================================
  ! ENRUTAMIENTO Y ASIGNACIÓN (Llamar al inicio del programa)
  ! ====================================================================
  subroutine init_wavespeed_solver()
    if (trim(metric_type) == 'Minkowski' .and. trim(geom_type) == 'Cartesian') then
      calc_wavespeeds => calc_wavespeeds_srhd
    else
      calc_wavespeeds => calc_wavespeeds_grhd
    end if
  end subroutine init_wavespeed_solver


  ! ====================================================================
  ! SOLUCIONADORES DE RIEMANN
  ! ====================================================================
  subroutine resolve_riemann_problem(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
    real*8, intent(in)  :: q_L(:), q_R(:)
    integer, intent(in) :: sweep_dir
    real*8, intent(in)  :: alpha, beta(3), g(3,3), sqg
    real*8, intent(out) :: flux_out(:)

    select case(riemann_solver_id)
    case(RS_HLLE)
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
    case(RS_HLLC)
      call calc_hllc_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
    case(RS_HLLD)
      call calc_hlld_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out) ! Futuro MHD
    case default
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
    end select
  end subroutine resolve_riemann_problem


  ! Subrutina auxiliar para convertir Q -> U usando estrictamente la métrica de la interfaz
  subroutine prim_to_cons_face(prim_state, cons_state, g, sqg)
    implicit none
    real*8, intent(in)  :: prim_state(:), g(3,3), sqg
    real*8, intent(out) :: cons_state(:)
    
    real*8 :: v_sq, h, W, v(3), v_cov(3)
    integer :: ii, jj

    v(1) = prim_state(eq_vx); v(2) = prim_state(eq_vy); v(3) = prim_state(eq_vz)
    
    ! Cuadrado de la velocidad Euleriana tridimensional (Tensorial Universal)
    v_sq = 0.0d0
    do ii = 1, 3
      do jj = 1, 3
        v_sq = v_sq + g(ii,jj) * v(ii) * v(jj)
      end do
    end do
    
    ! --- ESCUDO DE CAUSALIDAD NUMÉRICA ---
    if (v_sq >= v_max) then
      v(1) = v(1) * sqrt(v_max / v_sq)
      v(2) = v(2) * sqrt(v_max / v_sq)
      v(3) = v(3) * sqrt(v_max / v_sq)
      v_sq = v_max
    end if

    ! Calcular la velocidad covariante (v_i = gamma_ij * v^j)
    v_cov = 0.0d0
    do ii = 1, 3
      do jj = 1, 3
        v_cov(ii) = v_cov(ii) + g(ii,jj) * v(jj)
      end do
    end do

    h = 1.0d0 + g1 * prim_state(eq_pr) / prim_state(eq_de)
    W = 1.0d0 / sqrt(1.0d0 - v_sq)

    cons_state(eq_de) = prim_state(eq_de) * W
    cons_state(eq_pr) = prim_state(eq_de) * h * W**2 - prim_state(eq_pr) - prim_state(eq_de) * W
    
    ! Momentos espaciales covariantes S_i = rho * h * W^2 * v_i
    cons_state(eq_vx) = prim_state(eq_de) * h * W**2 * v_cov(1)
    cons_state(eq_vy) = prim_state(eq_de) * h * W**2 * v_cov(2)
    cons_state(eq_vz) = prim_state(eq_de) * h * W**2 * v_cov(3)

    cons_state = cons_state * sqg
  end subroutine prim_to_cons_face


  subroutine calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
    real*8, intent(in)  :: q_L(:), q_R(:)
    integer, intent(in) :: sweep_dir
    real*8, intent(in)  :: alpha, beta(3), g(3,3), sqg
    real*8, intent(out) :: flux_out(:)
    
    real*8 :: u_L(neq), u_R(neq), f_L(neq), f_R(neq)
    real*8 :: a_plus, a_minus, lambda_L_min, lambda_L_max, lambda_R_min, lambda_R_max

    ! 1. Evaluamos los flujos F(U) usando la métrica de la interfaz
    call calc_fluxes(q_L, f_L, sweep_dir, alpha, beta, g, sqg)
    call calc_fluxes(q_R, f_R, sweep_dir, alpha, beta, g, sqg)

    ! 2. Transformamos variables primitivas a conservadas
    call prim_to_cons_face(q_L, u_L, g, sqg)
    call prim_to_cons_face(q_R, u_R, g, sqg)

    ! 3. Cálculo de velocidades de onda
    call calc_wavespeeds(q_L, sweep_dir, alpha, beta, g, sqg, lambda_L_min, lambda_L_max)
    call calc_wavespeeds(q_R, sweep_dir, alpha, beta, g, sqg, lambda_R_min, lambda_R_max)

    ! 4. Velocidades de onda máxima y mínima en la interfaz
    a_plus  = max(0.0d0, lambda_L_max, lambda_R_max)
    a_minus = min(0.0d0, lambda_L_min, lambda_R_min)

    ! 5. Construcción del flujo numérico de Harten-Lax-van Leer (HLLE)
    if (abs(a_plus - a_minus) < 1.0d-14) then
      flux_out = 0.0d0
    else
      flux_out = (a_plus * f_L - a_minus * f_R + a_plus * a_minus * (u_R - u_L)) / (a_plus - a_minus)
    end if
  end subroutine calc_hlle_fluxes

  ! ====================================================================
  ! SOLUCIONADOR HLLC (Harten-Lax-van Leer-Contact)
  ! Captura la onda intermedia de contacto (restaura resolución de corte)
  ! ====================================================================
  
  subroutine calc_hllc_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
    real*8, intent(in)  :: q_L(:), q_R(:)
    integer, intent(in) :: sweep_dir
    real*8, intent(in)  :: alpha, beta(3), g(3,3), sqg
    real*8, intent(out) :: flux_out(:)
    
    real*8 :: u_L(neq), u_R(neq), f_L(neq), f_R(neq)
    real*8 :: u_hll(neq), f_hll(neq), D_star(neq)
    real*8 :: a_plus, a_minus, lambda_L_min, lambda_L_max, lambda_R_min, lambda_R_max
    real*8 :: E_hll, M_hll, FE_hll, FM_hll
    real*8 :: A_quad, B_quad, C_quad, discriminant, lambda_star, p_star

    ! 1. Evaluamos los flujos F(U) usando la métrica de la interfaz
    call calc_fluxes(q_L, f_L, sweep_dir, alpha, beta, g, sqg)
    call calc_fluxes(q_R, f_R, sweep_dir, alpha, beta, g, sqg)
    
    ! 2. Transformamos variables primitivas a conservadas (rutina de interfaz)
    call prim_to_cons_face(q_L, u_L, g, sqg)
    call prim_to_cons_face(q_R, u_R, g, sqg)

    ! 3. LLAMADA AL PUNTERO DE VELOCIDADES DE ONDA
    call calc_wavespeeds(q_L, sweep_dir, alpha, beta, g, sqg, lambda_L_min, lambda_L_max)
    call calc_wavespeeds(q_R, sweep_dir, alpha, beta, g, sqg, lambda_R_min, lambda_R_max)

    a_plus  = max(0.0d0, lambda_L_max, lambda_R_max)
    a_minus = min(0.0d0, lambda_L_min, lambda_R_min)

    if (abs(a_plus - a_minus) < 1.0d-14) then
      flux_out = 0.0d0
      return
    end if

    u_hll = (a_plus * u_R - a_minus * u_L + f_L - f_R) / (a_plus - a_minus)
    f_hll = (a_plus * f_L - a_minus * f_R + a_plus * a_minus * (u_R - u_L)) / (a_plus - a_minus)

    E_hll = u_hll(eq_pr) + u_hll(eq_de)
    FE_hll = f_hll(eq_pr) + f_hll(eq_de)
    
    if (sweep_dir == DIR_X) then
      M_hll = u_hll(eq_vx) ; FM_hll = f_hll(eq_vx)
    else if (sweep_dir == DIR_Y) then
      M_hll = u_hll(eq_vy) ; FM_hll = f_hll(eq_vy)
    else
      M_hll = u_hll(eq_vz) ; FM_hll = f_hll(eq_vz)
    end if

    A_quad = FE_hll
    B_quad = -(E_hll + FM_hll)
    C_quad = M_hll
    discriminant = B_quad**2 - 4.0d0 * A_quad * C_quad

    if (abs(A_quad) < 1.0d-14) then
      lambda_star = -C_quad / B_quad
    else
      if (discriminant < 0.0d0) then
        flux_out = f_hll 
        return
      end if
      lambda_star = (-B_quad - sign(1.0d0, B_quad) * sqrt(discriminant)) / (2.0d0 * A_quad)
      if (lambda_star < a_minus .or. lambda_star > a_plus) then
        lambda_star = (-B_quad + sign(1.0d0, B_quad) * sqrt(discriminant)) / (2.0d0 * A_quad)
      end if
    end if

    p_star = FM_hll - lambda_star * M_hll

    D_star = 0.0d0
    D_star(eq_pr) = lambda_star
    if (sweep_dir == DIR_X) D_star(eq_vx) = 1.0d0
    if (sweep_dir == DIR_Y) D_star(eq_vy) = 1.0d0
    if (sweep_dir == DIR_Z) D_star(eq_vz) = 1.0d0

    if (lambda_star >= 0.0d0) then
      flux_out = (lambda_star * (a_minus * u_L - f_L) + a_minus * p_star * D_star) / (a_minus - lambda_star)
    else
      flux_out = (lambda_star * (a_plus * u_R - f_R) + a_plus * p_star * D_star) / (a_plus - lambda_star)
    end if
  end subroutine calc_hllc_fluxes


  ! ====================================================================
  ! Esqueleto Reservado: Solucionador HLLD (Miyoshi & Kusano 2005 / Mignone 2009)
  ! Diseñado para Magnetohidrodinámica (MHD). Resuelve 5 ondas incluyendo 
  ! las discontinuidades de Alfvén. 
  ! ====================================================================
  
  subroutine calc_hlld_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
    real*8, intent(in)  :: q_L(:), q_R(:)
    integer, intent(in) :: sweep_dir
    real*8, intent(in)  :: alpha, beta(3), g(3,3), sqg
    real*8, intent(out) :: flux_out(:)
    real*8 :: dummy
    integer :: dummy_int
    
    ! Silenciar warnings del compilador con las nuevas variables
    dummy = q_L(1)
    dummy = q_R(1)
    dummy = alpha
    dummy = beta(1)
    dummy = g(1,1)
    dummy = sqg
    flux_out(1) = dummy
    dummy_int = sweep_dir
    
    ! Prevención de seguridad: Detener el código si se llama por accidente
    ! antes de ser programado.
    write(*,*) "=========================================================="
    write(*,*) " FATAL ERROR: El solucionador HLLD aun no ha sido "
    write(*,*) " implementado en este codigo."
    write(*,*) " Por favor, seleccione HLLE o HLLC para simulaciones."
    write(*,*) "=========================================================="
    stop

    ! Aquí irá en el futuro:
    ! 1. Cálculo de velocidades de onda rápidas (lambda_L, lambda_R)
    ! 2. Cálculo de la velocidad de contacto (lambda_star)
    ! 3. Cálculo de las velocidades de Alfvén (lambda_A_L, lambda_A_R)
    ! 4. Ensamblaje de las 4 regiones intermedias (L*, L**, R**, R*)
    
  end subroutine calc_hlld_fluxes

  ! Opción A: Relatividad Especial Pura (SRHD) - Minkowski
  subroutine calc_wavespeeds_srhd(prim_state, sweep_dir, alpha, beta, g, sqg, wave_min, wave_max)
    implicit none
    real*8, intent(in)  :: prim_state(:)
    integer, intent(in) :: sweep_dir
    real*8, intent(in)  :: alpha, beta(3), g(3,3), sqg
    real*8, intent(out) :: wave_min, wave_max
    
    real*8 :: cs, cs_sq, term1, term2
    real*8 :: eig1, eig2, eig3, enthalpy, v_sq
    real*8 :: v_norm, v(3)
    integer :: ii, jj
    
    v(1) = prim_state(eq_vx); v(2) = prim_state(eq_vy); v(3) = prim_state(eq_vz)

    ! Cuadrado de la velocidad Euleriana tridimensional
    v_sq = 0.0d0
    do ii = 1, 3
      do jj = 1, 3
        v_sq = v_sq + g(ii,jj) * v(ii) * v(jj)
      end do
    end do
    
    if (v_sq >= v_max) then
      v(1) = v(1) * sqrt(v_max / v_sq)
      v(2) = v(2) * sqrt(v_max / v_sq)
      v(3) = v(3) * sqrt(v_max / v_sq)
      v_sq = v_max
    end if

    enthalpy = 1.0d0 + g1 * prim_state(eq_pr) / prim_state(eq_de)
    cs_sq = adb_idx * prim_state(eq_pr) / (prim_state(eq_de) * enthalpy)
    cs = sqrt(cs_sq)

    if (sweep_dir == DIR_X) then
      v_norm = v(1)
    else if (sweep_dir == DIR_Y) then
      v_norm = v(2)
    else 
      v_norm = v(3)
    end if

    ! Ecuaciones simplificadas para fondo plano
    term1 = 1.0d0 / (1.0d0 - v_sq * cs_sq)
    term2 = cs * sqrt((1.0d0 - v_sq) * (1.0d0/term1 - v_norm**2 * (1.0d0 - cs_sq)))
    
    eig1 = v_norm
    eig2 = term1 * (v_norm * (1.0d0 - cs_sq) + term2)
    eig3 = term1 * (v_norm * (1.0d0 - cs_sq) - term2)

    wave_max = max(0.0d0, eig1, eig2, eig3)
    wave_min = min(0.0d0, eig1, eig2, eig3)
  end subroutine calc_wavespeeds_srhd


  ! Opción B: Relatividad General (GRHD) - Curvo
  subroutine calc_wavespeeds_grhd(prim_state, sweep_dir, alpha, beta, g, sqg, wave_min, wave_max)
    implicit none
    real*8, intent(in)  :: prim_state(:)
    integer, intent(in) :: sweep_dir
    real*8, intent(in)  :: alpha, beta(3), g(3,3), sqg
    real*8, intent(out) :: wave_min, wave_max
    
    real*8 :: cs, cs_sq, term1, term2
    real*8 :: eig1, eig2, eig3, enthalpy, v_sq
    real*8 :: v_norm, beta_norm, gamma_up_ii
    real*8 :: v(3)
    integer :: ii, jj

    v(1) = prim_state(eq_vx); v(2) = prim_state(eq_vy); v(3) = prim_state(eq_vz)

    ! Cuadrado de la velocidad con métrica espacial completa
    v_sq = 0.0d0
    do ii = 1, 3
      do jj = 1, 3
        v_sq = v_sq + g(ii,jj) * v(ii) * v(jj)
      end do
    end do
    
    if (v_sq >= v_max) then
      v(1) = v(1) * sqrt(v_max / v_sq)
      v(2) = v(2) * sqrt(v_max / v_sq)
      v(3) = v(3) * sqrt(v_max / v_sq)
      v_sq = v_max
    end if

    enthalpy = 1.0d0 + g1 * prim_state(eq_pr) / prim_state(eq_de)
    cs_sq = adb_idx * prim_state(eq_pr) / (prim_state(eq_de) * enthalpy)
    cs = sqrt(cs_sq)

    ! Extracción analítica del componente contravariante diagonal (gamma^ii) 
    ! usando determinantes de submatrices (Cofactores / det_g). det_g = sqg^2
    if (sweep_dir == DIR_X) then
      v_norm = v(1) ; beta_norm = beta(1) 
      gamma_up_ii = (g(2,2)*g(3,3) - g(2,3)*g(3,2)) / (sqg**2)
    else if (sweep_dir == DIR_Y) then
      v_norm = v(2) ; beta_norm = beta(2)
      gamma_up_ii = (g(1,1)*g(3,3) - g(1,3)*g(3,1)) / (sqg**2)
    else
      v_norm = v(3) ; beta_norm = beta(3)
      gamma_up_ii = (g(1,1)*g(2,2) - g(1,2)*g(2,1)) / (sqg**2)
    end if

    ! Ecuaciones exactas acopladas con el lapso, shift y métrica inversa
    term1 = alpha / (1.0d0 - v_sq * cs_sq)
    term2 = cs * sqrt( (1.0d0 - v_sq) * (gamma_up_ii * (1.0d0 - v_sq * cs_sq) - v_norm**2 * (1.0d0 - cs_sq)) )
    
    eig1 = alpha * v_norm - beta_norm
    eig2 = term1 * (v_norm * (1.0d0 - cs_sq) + term2) - beta_norm
    eig3 = term1 * (v_norm * (1.0d0 - cs_sq) - term2) - beta_norm

    wave_max = max(0.0d0, eig1, eig2, eig3)
    wave_min = min(0.0d0, eig1, eig2, eig3)
  end subroutine calc_wavespeeds_grhd

end module fluxes