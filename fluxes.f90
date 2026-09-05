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
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  
  ! Hacemos público el enrutador de Riemann y el inicializador del puntero
  public :: resolve_riemann_problem, init_wavespeed_solver

  ! ====================================================================
  ! INTERFAZ ABSTRACTA Y DECLARACIÓN DEL PUNTERO
  ! ====================================================================
  abstract interface
    subroutine wavespeed_interface(prim_state, sweep_dir, alpha, beta, g, sqg, wave_min, wave_max)
      implicit none
      real*8, intent(in)  :: prim_state(:)
      integer, intent(in) :: sweep_dir
      real*8, intent(in)  :: alpha, beta(3), g(3,3), sqg
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
  !
  ! El HLLC relativista de Mignone & Bodo (2005) esta formulado en una
  ! base cartesiana ortonormal. En GRHD no se puede aplicar su cuadratica
  ! directamente a S_i densitizados ni a velocidades coordenadas. Para
  ! cada cara se construye aqui un triedro ortonormal cuya primera base es
  ! normal a x^d=cte, se resuelve el problema SRHD local y se transforma
  ! el flujo nuevamente a la base coordenada.
  !
  ! Una cara fija en coordenadas se mueve en ese marco local con
  !   w = beta^d / (alpha * sqrt(gamma^dd)).
  ! Por eso se transforma F_hat-w*U_hat y se compara w con las ondas del
  ! abanico. Ante degeneraciones se vuelve al HLLE GRHD robusto.
  ! ====================================================================

  subroutine calc_hllc_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
    real*8, intent(in)  :: q_L(:), q_R(:)
    integer, intent(in) :: sweep_dir
    real*8, intent(in)  :: alpha, beta(3), g(3,3), sqg
    real*8, intent(out) :: flux_out(:)

    real*8 :: q_L_hat(neq), q_R_hat(neq)
    real*8 :: u_L(neq), u_R(neq), f_L(neq), f_R(neq)
    real*8 :: u_hll(neq), f_hll(neq), pressure_vector(neq)
    real*8 :: u_star_L(neq), u_star_R(neq), f_star_L(neq), f_star_R(neq)
    real*8 :: u_selected(neq), f_selected(neq)
    real*8 :: lambda_L_min, lambda_L_max, lambda_R_min, lambda_R_max
    real*8 :: lambda_left, lambda_right
    real*8 :: E_hll, M_hll, FE_hll, FM_hll
    real*8 :: A_quad, B_quad, C_quad, discriminant, lambda_star, p_star
    real*8 :: A_left, A_right, B_left, B_right, p_star_left, p_star_right
    real*8 :: frame(3,3), coframe(3,3), gamma_up_dd, face_speed
    real*8 :: denom, disc_scale, coefficient_scale, root_alt, root_tol
    real*8 :: pressure_tolerance, state_scale
    logical :: ok, root_is_valid, quadratic_case

    if (sweep_dir < DIR_X .or. sweep_dir > DIR_Z .or. alpha <= 0.0d0 .or. &
        sqg <= 0.0d0 .or. .not. ieee_is_finite(alpha) .or. &
        .not. ieee_is_finite(sqg)) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if

    ! 1. Triedro espacial local: e_(1) es normal a x^d=cte.
    call build_face_orthonormal_frame(g, sweep_dir, frame, coframe, gamma_up_dd, ok)
    if (.not. ok) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if

    call primitives_to_local_frame(q_L, coframe, q_L_hat)
    call primitives_to_local_frame(q_R, coframe, q_R_hat)

    ! 2. Estados, flujos y ondas acusticas en el marco minkowskiano local.
    call calc_srhd_local_state(q_L_hat, u_L, f_L, lambda_L_min, lambda_L_max, ok)
    if (.not. ok) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if
    call calc_srhd_local_state(q_R_hat, u_R, f_R, lambda_R_min, lambda_R_max, ok)
    if (.not. ok) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if

    lambda_left  = min(lambda_L_min, lambda_R_min)
    lambda_right = max(lambda_L_max, lambda_R_max)
    face_speed = beta(sweep_dir) / (alpha * sqrt(gamma_up_dd))

    if (.not. ieee_is_finite(face_speed) .or. &
        lambda_right - lambda_left <= 1.0d-13) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if

    ! Si la cara esta fuera del abanico, el flujo es puramente upwind.
    if (face_speed <= lambda_left) then
      call local_flux_to_coordinate(u_L, f_L, face_speed, alpha, sqg, &
                                    gamma_up_dd, coframe, flux_out, ok)
      if (.not. ok) call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    else if (face_speed >= lambda_right) then
      call local_flux_to_coordinate(u_R, f_R, face_speed, alpha, sqg, &
                                    gamma_up_dd, coframe, flux_out, ok)
      if (.not. ok) call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if

    ! 3. Estado y flujo HLL no densitizados en la base ortonormal.
    denom = lambda_right - lambda_left
    u_hll = (lambda_right * u_R - lambda_left * u_L + f_L - f_R) / denom
    f_hll = (lambda_right * f_L - lambda_left * f_R + &
             lambda_right * lambda_left * (u_R - u_L)) / denom

    if (.not. all(ieee_is_finite(u_hll)) .or. .not. all(ieee_is_finite(f_hll))) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if

    E_hll = u_hll(eq_pr) + u_hll(eq_de)
    FE_hll = f_hll(eq_pr) + f_hll(eq_de)
    M_hll = u_hll(eq_vx)
    FM_hll = f_hll(eq_vx)

    A_quad = FE_hll
    B_quad = -(E_hll + FM_hll)
    C_quad = M_hll
    discriminant = B_quad**2 - 4.0d0 * A_quad * C_quad

    disc_scale = max(tiny(1.0d0), B_quad**2, abs(4.0d0 * A_quad * C_quad))
    if (.not. ieee_is_finite(discriminant) .or. discriminant < -1.0d-12 * disc_scale) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if
    discriminant = max(0.0d0, discriminant)

    ! La raiz con signo menos es la rama fisica de Mignone--Bodo. Se usa
    ! una forma equivalente cuando la resta produciria cancelacion severa.
    coefficient_scale = max(tiny(1.0d0), abs(B_quad), abs(C_quad))
    quadratic_case = abs(A_quad) > 1.0d-13 * coefficient_scale
    if (.not. quadratic_case) then
      if (abs(B_quad) <= 1.0d-14 * max(tiny(1.0d0), abs(C_quad))) then
        call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
        return
      end if
      lambda_star = -C_quad / B_quad
    else
      denom = -B_quad - sqrt(discriminant)
      if (abs(denom) > 1.0d-14 * &
          max(tiny(1.0d0), abs(B_quad), sqrt(discriminant))) then
        lambda_star = denom / (2.0d0 * A_quad)
      else
        denom = -B_quad + sqrt(discriminant)
        if (abs(denom) <= 1.0d-14 * &
            max(tiny(1.0d0), abs(B_quad), sqrt(discriminant))) then
          call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
          return
        end if
        lambda_star = (2.0d0 * C_quad) / denom
      end if
    end if

    root_tol = 1.0d-10 * max(1.0d0, abs(lambda_left), abs(lambda_right))
    root_is_valid = ieee_is_finite(lambda_star) .and. &
                    lambda_star >= lambda_left - root_tol .and. &
                    lambda_star <= lambda_right + root_tol .and. &
                    abs(lambda_star) <= 1.0d0 + root_tol

    ! Por redondeo excepcional se prueba la otra raiz antes de abandonar HLLC.
    if (.not. root_is_valid .and. quadratic_case) then
      root_alt = (-B_quad + sqrt(discriminant)) / (2.0d0 * A_quad)
      if (ieee_is_finite(root_alt) .and. root_alt >= lambda_left - root_tol .and. &
          root_alt <= lambda_right + root_tol .and. abs(root_alt) <= 1.0d0 + root_tol) then
        lambda_star = root_alt
        root_is_valid = .true.
      end if
    end if

    if (.not. root_is_valid) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if
    lambda_star = min(lambda_right, max(lambda_left, lambda_star))

    ! Ecuacion (17) de Mignone--Bodo, evaluada desde ambos lados:
    !   p* = (A lambda* - B)/(1-lambda lambda*)
    ! con A=lambda E-m_n y B=m_n(lambda-v_n)-p.
    A_left = lambda_left * (u_L(eq_pr) + u_L(eq_de)) - u_L(eq_vx)
    B_left = u_L(eq_vx) * (lambda_left - q_L_hat(eq_vx)) - q_L_hat(eq_pr)
    A_right = lambda_right * (u_R(eq_pr) + u_R(eq_de)) - u_R(eq_vx)
    B_right = u_R(eq_vx) * (lambda_right - q_R_hat(eq_vx)) - q_R_hat(eq_pr)

    denom = 1.0d0 - lambda_left * lambda_star
    if (abs(denom) <= 1.0d-13) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if
    p_star_left = (A_left * lambda_star - B_left) / denom

    denom = 1.0d0 - lambda_right * lambda_star
    if (abs(denom) <= 1.0d-13) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if
    p_star_right = (A_right * lambda_star - B_right) / denom

    state_scale = max(tiny(1.0d0), &
                      abs(u_L(eq_pr) + u_L(eq_de)), abs(u_R(eq_pr) + u_R(eq_de)), &
                      abs(u_L(eq_vx)), abs(u_R(eq_vx)))
    pressure_tolerance = 1.0d-8 * &
      max(tiny(1.0d0), abs(q_L_hat(eq_pr)), abs(q_R_hat(eq_pr)), &
          abs(p_star_left), abs(p_star_right)) + &
      1.0d2 * epsilon(1.0d0) * state_scale
    if (.not. ieee_is_finite(p_star_left) .or. .not. ieee_is_finite(p_star_right) .or. &
        abs(p_star_left - p_star_right) > pressure_tolerance) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if
    p_star = 0.5d0 * (p_star_left + p_star_right)
    if (.not. ieee_is_finite(p_star) .or. &
        abs(lambda_left - lambda_star) <= 1.0d-13 .or. &
        abs(lambda_right - lambda_star) <= 1.0d-13) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if

    ! N* expresa el termino de presion en el orden (D,tau,S_n,S_t1,S_t2).
    pressure_vector = 0.0d0
    pressure_vector(eq_pr) = lambda_star
    pressure_vector(eq_vx) = 1.0d0

    u_star_L = (lambda_left * u_L - f_L + p_star * pressure_vector) / &
               (lambda_left - lambda_star)
    u_star_R = (lambda_right * u_R - f_R + p_star * pressure_vector) / &
               (lambda_right - lambda_star)
    f_star_L = f_L + lambda_left * (u_star_L - u_L)
    f_star_R = f_R + lambda_right * (u_star_R - u_R)

    if (.not. conservative_state_is_admissible(u_star_L) .or. &
        .not. conservative_state_is_admissible(u_star_R) .or. &
        .not. all(ieee_is_finite(f_star_L)) .or. &
        .not. all(ieee_is_finite(f_star_R))) then
      call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
      return
    end if

    ! 4. Estado que cruza una cara fija de la malla.
    if (face_speed <= lambda_star) then
      u_selected = u_star_L
      f_selected = f_star_L
    else
      u_selected = u_star_R
      f_selected = f_star_R
    end if

    call local_flux_to_coordinate(u_selected, f_selected, face_speed, alpha, sqg, &
                                  gamma_up_dd, coframe, flux_out, ok)
    if (.not. ok) call calc_hlle_fluxes(q_L, q_R, sweep_dir, alpha, beta, g, sqg, flux_out)
  end subroutine calc_hllc_fluxes

  ! Construye una base espacial ortonormal e_(a)^i y su co-base e^(a)_i.
  ! La primera direccion es normal a la superficie coordenada x^d=cte:
  !   e_(1)^i = gamma^(id) / sqrt(gamma^dd).
  subroutine build_face_orthonormal_frame(g, sweep_dir, frame, coframe, gamma_up_dd, ok)
    implicit none
    real*8, intent(in) :: g(3,3)
    integer, intent(in) :: sweep_dir
    real*8, intent(out) :: frame(3,3), coframe(3,3), gamma_up_dd
    logical, intent(out) :: ok

    real*8 :: ginv(3,3), det_g, seed(3), projection, norm_sq
    integer :: tangent_dirs(2), basis_idx, i, j, a

    ok = .false.
    frame = 0.0d0
    coframe = 0.0d0
    gamma_up_dd = 0.0d0

    if (.not. all(ieee_is_finite(g))) return

    det_g = g(1,1)*(g(2,2)*g(3,3) - g(2,3)*g(3,2)) &
          - g(1,2)*(g(2,1)*g(3,3) - g(2,3)*g(3,1)) &
          + g(1,3)*(g(2,1)*g(3,2) - g(2,2)*g(3,1))
    if (.not. ieee_is_finite(det_g) .or. det_g <= 1.0d-30) return

    ginv(1,1) =  (g(2,2)*g(3,3) - g(2,3)*g(3,2)) / det_g
    ginv(1,2) =  (g(1,3)*g(3,2) - g(1,2)*g(3,3)) / det_g
    ginv(1,3) =  (g(1,2)*g(2,3) - g(1,3)*g(2,2)) / det_g
    ginv(2,1) =  (g(2,3)*g(3,1) - g(2,1)*g(3,3)) / det_g
    ginv(2,2) =  (g(1,1)*g(3,3) - g(1,3)*g(3,1)) / det_g
    ginv(2,3) =  (g(1,3)*g(2,1) - g(1,1)*g(2,3)) / det_g
    ginv(3,1) =  (g(2,1)*g(3,2) - g(2,2)*g(3,1)) / det_g
    ginv(3,2) =  (g(1,2)*g(3,1) - g(1,1)*g(3,2)) / det_g
    ginv(3,3) =  (g(1,1)*g(2,2) - g(1,2)*g(2,1)) / det_g
    if (.not. all(ieee_is_finite(ginv))) return

    gamma_up_dd = ginv(sweep_dir, sweep_dir)
    if (gamma_up_dd <= 1.0d-30 .or. .not. ieee_is_finite(gamma_up_dd)) return
    frame(:,1) = ginv(:,sweep_dir) / sqrt(gamma_up_dd)

    select case(sweep_dir)
    case(DIR_X)
      tangent_dirs = (/ DIR_Y, DIR_Z /)
    case(DIR_Y)
      tangent_dirs = (/ DIR_Z, DIR_X /)
    case(DIR_Z)
      tangent_dirs = (/ DIR_X, DIR_Y /)
    case default
      return
    end select

    ! Gram--Schmidt con el producto interno gamma_ij.
    do basis_idx = 2, 3
      seed = 0.0d0
      seed(tangent_dirs(basis_idx-1)) = 1.0d0
      do a = 1, basis_idx - 1
        projection = metric_dot(seed, frame(:,a), g)
        seed = seed - projection * frame(:,a)
      end do
      norm_sq = metric_dot(seed, seed, g)
      if (.not. ieee_is_finite(norm_sq) .or. norm_sq <= 1.0d-28) return
      frame(:,basis_idx) = seed / sqrt(norm_sq)
    end do

    do a = 1, 3
      do i = 1, 3
        do j = 1, 3
          coframe(a,i) = coframe(a,i) + g(i,j) * frame(j,a)
        end do
      end do
    end do

    if (.not. all(ieee_is_finite(frame)) .or. .not. all(ieee_is_finite(coframe))) return
    ok = .true.
  end subroutine build_face_orthonormal_frame


  real*8 function metric_dot(vector_a, vector_b, g)
    implicit none
    real*8, intent(in) :: vector_a(3), vector_b(3), g(3,3)
    integer :: i, j

    metric_dot = 0.0d0
    do i = 1, 3
      do j = 1, 3
        metric_dot = metric_dot + g(i,j) * vector_a(i) * vector_b(j)
      end do
    end do
  end function metric_dot


  subroutine primitives_to_local_frame(prim_coordinate, coframe, prim_local)
    implicit none
    real*8, intent(in) :: prim_coordinate(:), coframe(3,3)
    real*8, intent(out) :: prim_local(:)
    integer :: a, i

    prim_local = prim_coordinate
    prim_local(eq_vx:eq_vz) = 0.0d0
    do a = 1, 3
      do i = 1, 3
        prim_local(eq_vx + a - 1) = prim_local(eq_vx + a - 1) + &
          coframe(a,i) * prim_coordinate(eq_vx + i - 1)
      end do
    end do
  end subroutine primitives_to_local_frame


  ! Variables conservativas, flujo normal y velocidades extremas SRHD en
  ! una base ortonormal. El orden permanece (D,tau,S_n,S_t1,S_t2).
  subroutine calc_srhd_local_state(prim_local, cons_state, flux_state, wave_min, wave_max, ok)
    implicit none
    real*8, intent(in) :: prim_local(:)
    real*8, intent(out) :: cons_state(:), flux_state(:), wave_min, wave_max
    logical, intent(out) :: ok

    real*8 :: rho, pressure, velocity(3), v_sq, h, W, cs_sq
    real*8 :: D, tau, momentum(3), denom, radicand, acoustic_term

    ok = .false.
    cons_state = 0.0d0
    flux_state = 0.0d0
    wave_min = 0.0d0
    wave_max = 0.0d0

    rho = prim_local(eq_de)
    pressure = prim_local(eq_pr)
    velocity = prim_local(eq_vx:eq_vz)
    if (.not. ieee_is_finite(rho) .or. .not. ieee_is_finite(pressure) .or. &
        .not. all(ieee_is_finite(velocity)) .or. rho <= 0.0d0 .or. pressure < 0.0d0) return

    v_sq = sum(velocity**2)
    if (.not. ieee_is_finite(v_sq) .or. v_sq < 0.0d0) return
    ! Los estados reconstruidos no causales se delegan a HLLE. No se
    ! modifican silenciosamente porque q_hat tambien entra en p*.
    if (v_sq >= v_max) return

    h = 1.0d0 + g1 * pressure / rho
    if (.not. ieee_is_finite(h) .or. h <= 0.0d0) return
    cs_sq = adb_idx * pressure / (rho * h)
    if (.not. ieee_is_finite(cs_sq) .or. cs_sq < 0.0d0 .or. cs_sq >= 1.0d0) return

    W = 1.0d0 / sqrt(1.0d0 - v_sq)
    D = rho * W
    tau = rho * h * W**2 - pressure - D
    momentum = rho * h * W**2 * velocity

    cons_state(eq_de) = D
    cons_state(eq_pr) = tau
    cons_state(eq_vx:eq_vz) = momentum

    flux_state(eq_de) = D * velocity(1)
    flux_state(eq_pr) = momentum(1) - D * velocity(1)
    flux_state(eq_vx) = momentum(1) * velocity(1) + pressure
    flux_state(eq_vy) = momentum(2) * velocity(1)
    flux_state(eq_vz) = momentum(3) * velocity(1)

    denom = 1.0d0 - v_sq * cs_sq
    radicand = (1.0d0 - v_sq) * &
      (1.0d0 - v_sq * cs_sq - velocity(1)**2 * (1.0d0 - cs_sq))
    if (denom <= 1.0d-14 .or. radicand < -1.0d-12) return
    acoustic_term = sqrt(cs_sq) * sqrt(max(0.0d0, radicand))
    wave_min = (velocity(1) * (1.0d0 - cs_sq) - acoustic_term) / denom
    wave_max = (velocity(1) * (1.0d0 - cs_sq) + acoustic_term) / denom

    if (.not. all(ieee_is_finite(cons_state)) .or. &
        .not. all(ieee_is_finite(flux_state)) .or. &
        .not. ieee_is_finite(wave_min) .or. .not. ieee_is_finite(wave_max)) return
    ok = .true.
  end subroutine calc_srhd_local_state


  logical function conservative_state_is_admissible(cons_state)
    implicit none
    real*8, intent(in) :: cons_state(:)
    real*8 :: total_energy, momentum_sq, admissibility_scale

    conservative_state_is_admissible = .false.
    if (.not. all(ieee_is_finite(cons_state))) return
    if (cons_state(eq_de) <= 0.0d0) return

    total_energy = cons_state(eq_pr) + cons_state(eq_de)
    momentum_sq = sum(cons_state(eq_vx:eq_vz)**2)
    ! Criterio relativo: debe seguir siendo valido en la atmosfera, donde
    ! D y E son muchos ordenes de magnitud menores que uno.
    admissibility_scale = max(tiny(1.0d0), total_energy**2, &
                              cons_state(eq_de)**2 + momentum_sq)
    if (total_energy <= 0.0d0) return
    if (total_energy**2 - cons_state(eq_de)**2 - momentum_sq <= &
        1.0d-13 * admissibility_scale) return

    conservative_state_is_admissible = .true.
  end function conservative_state_is_admissible


  ! Convierte sqrt(gamma) alpha sqrt(gamma^dd) (F_hat-w*U_hat)
  ! al flujo de las variables de Valencia en la base coordenada.
  subroutine local_flux_to_coordinate(cons_local, flux_local, face_speed, alpha, sqg, &
                                      gamma_up_dd, coframe, flux_coordinate, ok)
    implicit none
    real*8, intent(in) :: cons_local(:), flux_local(:), face_speed
    real*8, intent(in) :: alpha, sqg, gamma_up_dd, coframe(3,3)
    real*8, intent(out) :: flux_coordinate(:)
    logical, intent(out) :: ok

    real*8 :: moving_flux(neq), factor
    integer :: i, a

    ok = .false.
    flux_coordinate = 0.0d0
    moving_flux = flux_local - face_speed * cons_local
    factor = sqg * alpha * sqrt(gamma_up_dd)
    if (.not. ieee_is_finite(factor) .or. .not. all(ieee_is_finite(moving_flux))) return

    flux_coordinate(eq_de) = factor * moving_flux(eq_de)
    flux_coordinate(eq_pr) = factor * moving_flux(eq_pr)
    do i = 1, 3
      do a = 1, 3
        flux_coordinate(eq_vx + i - 1) = flux_coordinate(eq_vx + i - 1) + &
          factor * coframe(a,i) * moving_flux(eq_vx + a - 1)
      end do
    end do

    if (.not. all(ieee_is_finite(flux_coordinate))) return
    ok = .true.
  end subroutine local_flux_to_coordinate


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
    
    real*8 :: cs, cs_sq, term1, term2, det_g
    real*8 :: eig1, eig2, eig3, enthalpy, v_sq
    real*8 :: v_norm, beta_norm, gamma_up_ii
    real*8 :: v(3)
    integer :: ii, jj
    real*8 :: gamma_up_old, error_ondas

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

    ! --- CÁLCULO EXACTO DEL DETERMINANTE (Evita la pérdida de bits de sqg**2) ---
    det_g = g(1,1)*(g(2,2)*g(3,3) - g(2,3)*g(3,2)) &
          - g(1,2)*(g(2,1)*g(3,3) - g(2,3)*g(3,1)) &
          + g(1,3)*(g(2,1)*g(3,2) - g(2,2)*g(3,1))

    ! Extracción analítica del componente contravariante diagonal (gamma^ii) 
    if (sweep_dir == DIR_X) then
      v_norm = v(1) ; beta_norm = beta(1) 
      gamma_up_ii = (g(2,2)*g(3,3) - g(2,3)*g(3,2)) / det_g
    else if (sweep_dir == DIR_Y) then
      v_norm = v(2) ; beta_norm = beta(2)
      gamma_up_ii = (g(1,1)*g(3,3) - g(1,3)*g(3,1)) / det_g
    else
      v_norm = v(3) ; beta_norm = beta(3)
      gamma_up_ii = (g(1,1)*g(2,2) - g(1,2)*g(2,1)) / det_g
    end if

    ! Ecuaciones exactas acopladas con el lapso, shift y métrica inversa
    term1 = alpha / (1.0d0 - v_sq * cs_sq)
    
    ! Blindaje max(0.0d0) agregado contra ruido de coma flotante extremo
    term2 = cs * sqrt( max(0.0d0, (1.0d0 - v_sq) * (gamma_up_ii * (1.0d0 - v_sq * cs_sq) - v_norm**2 * (1.0d0 - cs_sq))) )
    
    eig1 = alpha * v_norm - beta_norm
    eig2 = term1 * (v_norm * (1.0d0 - cs_sq) + term2) - beta_norm
    eig3 = term1 * (v_norm * (1.0d0 - cs_sq) - term2) - beta_norm

    wave_max = max(0.0d0, eig1, eig2, eig3)
    wave_min = min(0.0d0, eig1, eig2, eig3)


    ! --- SONDA DE DIAGNÓSTICO DE ONDAS ---
    if (n_steps == 0 .and. sweep_dir == DIR_X) then
      gamma_up_old = 1.0d0 / g(1,1)
      error_ondas = abs(gamma_up_ii - gamma_up_old)
      
      if (error_ondas > 1.0d-13) then
          print *, ">>> ¡ERROR EN VELOCIDADES DE ONDA!"
          print *, "gamma^11 Cramer =", gamma_up_ii
          print *, "gamma^11 Viejo  =", gamma_up_old
          print *, "g_tphi          =", g(1,3) ! Revisar arrastre
          stop
      end if
    end if
  end subroutine calc_wavespeeds_grhd

end module fluxes
