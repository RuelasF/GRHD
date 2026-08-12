! =======================================================================================
! Módulo: reconstruction (Reconstrucción Espacial y Limitadores de Flujo)
! ---------------------------------------------------------------------------------------
! Propósito: Interpolar las variables primitivas desde los centros de las celdas (i) 
! hacia las interfaces (i+1/2). Estos estados interpolados (q_L y q_R) son el "input" 
! que el Solucionador de Riemann necesita para resolver el problema de discontinuidad 
! local y calcular los flujos.
!
! Física/Matemática: Implementa esquemas TVD (Total Variation Diminishing) y 
! métodos de alto orden (WENO, MP5) para garantizar precisión espectral en zonas 
! suaves del fluido (turbulencia) y evitar oscilaciones térmicas no físicas 
! cerca de choques fuertes (Fenómeno de Gibbs).
! =======================================================================================
module reconstruction
  use variables
  implicit none
  private
  public :: reconstruct_1d_core

contains

  ! ===================================================================================
  ! LIMITADORES DE PENDIENTE (Slope Limiters)
  ! ===================================================================================
  ! Conjunto de funciones matemáticas para limitar los gradientes espaciales y 
  ! forzar la monotonicidad estricta en las interfaces.

  real*8 function minmod_2(a, b)
    implicit none
    real*8, intent(in) :: a, b
    minmod_2 = 0.5d0 * (sign(1.0d0, a) + sign(1.0d0, b)) * min(abs(a), abs(b))
  end function minmod_2

  real*8 function minmod_3(a, b, c)
    implicit none
    real*8, intent(in) :: a, b, c
    minmod_3 = sign(1.0d0, a) * max(0.0d0, min(abs(a), sign(1.0d0, a)*b, sign(1.0d0, a)*c))
  end function minmod_3

  real*8 function minmod_4(a, b, c, d)
    implicit none
    real*8, intent(in) :: a, b, c, d
    minmod_4 = 0.125d0 * (sign(1.0d0, a) + sign(1.0d0, b)) * abs(sign(1.0d0, a) + sign(1.0d0, c)) &
             * abs(sign(1.0d0, a) + sign(1.0d0, d)) * min(abs(a), abs(b), abs(c), abs(d))
  end function minmod_4


  ! ===================================================================================
  ! ESQUEMAS DE INTERPOLACIÓN DE ALTO ORDEN
  ! ===================================================================================

  ! Subrutina: weno3_interface (WENO3 con Epsilon Dinámico de Castro 2011)
  ! Combina 2 esténciles lineales para lograr 3er orden en zonas suaves.
  ! Utiliza un epsilon dinámico escalado con el tamaño de malla (dx) para 
  ! garantizar la convergencia teórica en puntos críticos.

  ! original
  ! function weno3_interface(q) result(weno_val)
  !   implicit none
  !   real*8, intent(in) :: q(-1:1)
  !   real*8 :: weno_val
  !   real*8 :: p0, p1
  !   real*8 :: beta_0, beta_1, tau_3
  !   real*8 :: alpha_0, alpha_1, weight_0, weight_1
    
  !   ! Parámetros constantes
  !   real*8, parameter :: eps_z = 1.0d-12
  !   real*8, parameter :: d0 = 2.0d0/3.0d0, d1 = 1.0d0/3.0d0

  !   ! 1. Polinomios candidatos de los sub-esténciles
  !   p0 =  0.5d0*q(0) + 0.5d0*q(1)
  !   p1 = -0.5d0*q(-1) + 1.5d0*q(0)

  !   ! 2. Indicadores de suavidad local (Jiang-Shu)
  !   beta_0 = (q(0) - q(1))**2
  !   beta_1 = (q(-1) - q(0))**2

  !   ! 3. Indicador de suavidad global (Magia de WENO-Z)
  !   tau_3 = abs(beta_0 - beta_1)

  !   ! 4. Pesos no lineales penalizados (WENO-Z)
  !   alpha_0 = d0 * (1.0d0 + (tau_3 / (beta_0 + eps_z)))
  !   alpha_1 = d1 * (1.0d0 + (tau_3 / (beta_1 + eps_z)))

  !   ! 5. Normalización estándar
  !   weight_0 = alpha_0 / (alpha_0 + alpha_1)
  !   weight_1 = alpha_1 / (alpha_0 + alpha_1)

  !   ! 6. Ensamblaje final
  !   weno_val = weight_0 * p0 + weight_1 * p1
  ! end function weno3_interface

  ! weno-m
  ! function weno3_interface(q) result(weno_val)
  !   implicit none
  !   real*8, intent(in) :: q(-1:1)
  !   real*8 :: weno_val
  !   real*8 :: p0, p1
  !   real*8 :: beta_0, beta_1
  !   real*8 :: alpha_0, alpha_1, w0_js, w1_js
  !   real*8 :: w0_map, w1_map, sum_map
    
  !   real*8, parameter :: eps = 1.0d-12
  !   real*8, parameter :: d0 = 2.0d0/3.0d0, d1 = 1.0d0/3.0d0

  !   ! 1. Polinomios candidatos
  !   p0 =  0.5d0*q(0) + 0.5d0*q(1)
  !   p1 = -0.5d0*q(-1) + 1.5d0*q(0)

  !   ! 2. Betas (WENO-JS)
  !   beta_0 = (q(0) - q(1))**2
  !   beta_1 = (q(-1) - q(0))**2

  !   ! 3. Pesos Clásicos Jiang-Shu
  !   alpha_0 = d0 / (beta_0 + eps)**2
  !   alpha_1 = d1 / (beta_1 + eps)**2
  !   w0_js = alpha_0 / (alpha_0 + alpha_1)
  !   w1_js = alpha_1 / (alpha_0 + alpha_1)

  !   ! =========================================================
  !   ! 4. LA MAGIA DE WENO-M (Mapped WENO)
  !   ! Esta función amplifica los pesos hacia d_k en zonas suaves
  !   ! g(w, d) = w * (d + d^2 - 3*d*w + w^2) / (d^2 + w*(1 - 2*d))
  !   ! =========================================================
  !   w0_map = w0_js * (d0 + d0**2 - 3.0d0*d0*w0_js + w0_js**2) / &
  !           (d0**2 + w0_js*(1.0d0 - 2.0d0*d0) + eps)

  !   w1_map = w1_js * (d1 + d1**2 - 3.0d0*d1*w1_js + w1_js**2) / &
  !           (d1**2 + w1_js*(1.0d0 - 2.0d0*d1) + eps)

  !   ! 5. Normalización final
  !   sum_map = w0_map + w1_map
  !   w0_map = w0_map / sum_map
  !   w1_map = w1_map / sum_map

  !   ! 6. Ensamblaje
  !   weno_val = w0_map * p0 + w1_map * p1
  ! end function weno3_interface

  ! weno-dinamico
  function weno3_interface(q) result(weno_val)
    implicit none
    real*8, intent(in) :: q(-1:)
    real*8 :: weno_val
    real*8 :: p0, p1
    real*8 :: beta_0, beta_1
    real*8 :: alpha_0, alpha_1, weight_0, weight_1
    
    ! --- Variables para tu Epsilon Híbrido ---
    real*8 :: shock_sensor, eps_hybrid
    real*8, parameter :: eps_min = 1.0d-14
    real*8, parameter :: d0 = 2.0d0/3.0d0, d1 = 1.0d0/3.0d0

    ! 1. Polinomios candidatos
    p0 =  0.5d0*q(0) + 0.5d0*q(1)
    p1 = -0.5d0*q(-1) + 1.5d0*q(0)

    ! 2. Indicadores de suavidad local
    beta_0 = (q(0) - q(1))**2
    beta_1 = (q(-1) - q(0))**2

    ! =========================================================
    ! 3. EL EPSILON ADAPTATIVO
    ! =========================================================
    ! Sensor adimensional: 0 (Suave) ---> 1 (Choque fuerte)
    shock_sensor = abs(beta_0 - beta_1) / (max(beta_0, beta_1) + eps_min)
    
    ! Válvula exponencial: Si es suave da dx, si es choque da eps_min
    eps_hybrid = eps_min + (q(-1)**2 + q(0)**2 + q(1)**2) * exp(-10.0d0 * shock_sensor)

    ! 4. Pesos Clásicos pero con el Epsilon Dinámico Inteligente
    alpha_0 = d0 / (eps_hybrid + beta_0)**2
    alpha_1 = d1 / (eps_hybrid + beta_1)**2

    ! 5. Normalización
    weight_0 = alpha_0 / (alpha_0 + alpha_1)
    weight_1 = alpha_1 / (alpha_0 + alpha_1)

    ! 6. Ensamblaje
    weno_val = weight_0 * p0 + weight_1 * p1
  end function weno3_interface

  ! Subrutina: weno5_interface (WENO5-Z)
  ! Esquema de 5to orden mejorado. Utiliza el indicador de suavidad global (tau5)
  ! para evitar la degradación del orden cerca de extremos suaves, superando
  ! al WENO5 clásico de Jiang-Shu.
  function weno5_interface(q) result(weno_val)
    implicit none
    real*8, intent(in) :: q(-2:)
    real*8 :: weno_val
    real*8 :: b0, b1, b2, tau5, a0, a1, a2, w0, w1, w2
    real*8 :: p0, p1, p2
    real*8, parameter :: eps = 1.0d-14
    real*8, parameter :: d0 = 0.3d0, d1 = 0.6d0, d2 = 0.1d0

    b0 = 13.0d0/12.0d0*(q(0) - 2.0d0*q(1) + q(2))**2 + 0.25d0*(3.0d0*q(0) - 4.0d0*q(1) + q(2))**2
    b1 = 13.0d0/12.0d0*(q(-1) - 2.0d0*q(0) + q(1))**2 + 0.25d0*(q(-1) - q(1))**2
    b2 = 13.0d0/12.0d0*(q(-2) - 2.0d0*q(-1) + q(0))**2 + 0.25d0*(q(-2) - 4.0d0*q(-1) + 3.0d0*q(0))**2

    tau5 = abs(b0 - b2)

    a0 = d0 * (1.0d0 + (tau5 / (eps + b0)))
    a1 = d1 * (1.0d0 + (tau5 / (eps + b1)))
    a2 = d2 * (1.0d0 + (tau5 / (eps + b2)))

    w0 = a0 / (a0 + a1 + a2)
    w1 = a1 / (a0 + a1 + a2)
    w2 = a2 / (a0 + a1 + a2)

    p0 = (1.0d0/3.0d0)*q(0)  + (5.0d0/6.0d0)*q(1) - (1.0d0/6.0d0)*q(2)
    p1 = -(1.0d0/6.0d0)*q(-1) + (5.0d0/6.0d0)*q(0) + (1.0d0/3.0d0)*q(1)
    p2 = (1.0d0/3.0d0)*q(-2) - (7.0d0/6.0d0)*q(-1) + (11.0d0/6.0d0)*q(0)

    weno_val = w0*p0 + w1*p1 + w2*p2
  end function weno5_interface


  ! Subrutina: mp5_interface
  ! Monotonicity Preserving de 5to Orden. Basado en polinomios que no cruzan
  ! límites geométricos locales. Extraordinariamente nítido en choques, aunque
  ! puede generar ligeras inestabilidades en la relatividad acoplada.
  function mp5_interface(q) result(mp5_val)
    implicit none
    real*8, intent(in) :: q(-2:)
    real*8 :: mp5_val, val_L, val_mp, d_jm1, d_j, d_jp1, dm4_jph, dm4_jmh
    real*8 :: val_ul, val_av, val_md, val_lc, val_min, val_max
    real*8, parameter :: alpha_mp5 = 4.0d0 
    real*8, parameter :: eps = 1.0d-10     

    val_L = (2.0d0*q(-2) - 13.0d0*q(-1) + 47.0d0*q(0) + 27.0d0*q(1) - 3.0d0*q(2)) / 60.0d0
    
    val_mp = q(0) + minmod_2(q(1)-q(0), alpha_mp5*(q(0)-q(-1)))

    if ((val_L - q(0))*(val_L - val_mp) < eps) then
      mp5_val = val_L
    else
      d_jm1 = q(-2) - 2.0d0*q(-1) + q(0)
      d_j   = q(-1) - 2.0d0*q(0) + q(1)
      d_jp1 = q(0) - 2.0d0*q(1) + q(2)

      dm4_jph = minmod_4(4.0d0*d_j - d_jp1, 4.0d0*d_jp1 - d_j, d_j, d_jp1)
      dm4_jmh = minmod_4(4.0d0*d_j - d_jm1, 4.0d0*d_jm1 - d_j, d_j, d_jm1)

      val_ul = q(0) + alpha_mp5*(q(0) - q(-1))
      val_av = 0.5d0*(q(0) + q(1))
      val_md = val_av - 0.5d0*dm4_jph
      val_lc = q(0) + 0.5d0*(q(0) - q(-1)) + 4.0d0*dm4_jmh/3.0d0

      val_min = max(min(q(0), q(1), val_md), min(q(0), val_ul, val_lc))
      val_max = min(max(q(0), q(1), val_md), max(q(0), val_ul, val_lc))

      mp5_val = val_L + minmod_2(val_min - val_L, val_max - val_L)
    end if
  end function mp5_interface


  ! Función: tvd_slope
  ! Selector dinámico (numérico) de limitadores de segundo orden.
  function tvd_slope(a, b, limiter_id) result(slope)
    real*8, intent(in) :: a, b
    integer, intent(in) :: limiter_id
    real*8 :: slope

    select case(limiter_id)
    case(LIM_MINMOD)
      slope = minmod_2(a, b)
    case(LIM_SUPERBEE)
      if (a*b <= 0.0d0) then
        slope = 0.0d0
      else
        slope = sign(1.0d0, a) * max(min(2.0d0*abs(a), abs(b)), min(abs(a), 2.0d0*abs(b)))
      end if
    case(LIM_MC)
      slope = minmod_3(2.0d0*b, 0.5d0*(a+b), 2.0d0*a)
    case default
      slope = minmod_2(a, b)
    end select
  end function tvd_slope


  ! ===================================================================================
  ! NÚCLEO PRINCIPAL DE RECONSTRUCCIÓN (Arquitectura Híbrida Vectorizada)
  ! ===================================================================================

  subroutine reconstruct_1d_core(prim_1d, nodes, q_L, q_R, sweep_dir)
    implicit none
    real*8, intent(in)  :: prim_1d(neq, -nghost:nodes+nghost)
    integer, intent(in) :: nodes, sweep_dir
    real*8, intent(out) :: q_L(neq, 0:nodes), q_R(neq, 0:nodes)
    integer :: i, k
    real*8  :: s_min_L, s_max_L, s_min_R, s_max_R
    real*8  :: p_jump, p_min, p_max
    real*8  :: rho_jump, rho_min, rho_max
    logical :: shock_sensor
    real*8, parameter :: atm_factor = 1.0d2 ! 100 veces el piso se considera atmósfera

    ! =====================================================================
    ! CONFIGURACIÓN DEL ESCUDO TÉRMICO (FALLBACK)
    ! Descomentar el conjunto que se desee utilizar
    ! =====================================================================

    ! --- CONJUNTO A: Agresivo  ---
    real*8  :: p_tol_tvd   = 0.80d0  !  80% var. para bajar a TVD
    real*8  :: rho_tol_tvd = 0.50d0  !  50% var. para bajar a TVD
    real*8  :: p_tol_weno  = 0.40d0  !  40% var. para bajar a WENO3
    real*8  :: rho_tol_weno= 0.15d0  !  15% var. para bajar a WENO3

    ! --- CONJUNTO B: Extremo  ---
    ! real*8  :: p_tol_tvd   = 1.50d0    ! 150% var. para bajar a TVD
    ! real*8  :: rho_tol_tvd = 0.90d0    !  90% var. para bajar a TVD
    ! real*8  :: p_tol_weno  = 0.80d0    !  80% var. para bajar a WENO3
    ! real*8  :: rho_tol_weno= 0.40d0    !  40% var. para bajar a WENO3
    ! =====================================================================

    ! --- CONJUNTO C: Acantilado de Vacío ---
    real*8  :: p_tol_godunov   = 10.0d0  ! Salto > 1000% -> Godunov puro
    real*8  :: rho_tol_godunov = 10.0d0

    select case(rec_method_id)
    
    case(REC_GODUNOV)
    
    !$OMP PARALLEL DO PRIVATE(i)
      do i = 0, nodes
        q_L(:, i) = prim_1d(:, i)
        q_R(:, i) = prim_1d(:, i+1)
      end do
    !$OMP END PARALLEL DO

    case(REC_TVD)

    if (use_shock_sensor) then
      !$OMP PARALLEL DO PRIVATE(k,i,s_min_L,s_max_L,s_min_R,s_max_R)
        do k = 1, neq
          do i = 0, nodes
            if (prim_1d(eq_de, i) <= atm_factor * rho_floor .or. prim_1d(eq_de, i+1) <= atm_factor * rho_floor) then
              q_L(k, i) = prim_1d(k, i)
              q_R(k, i) = prim_1d(k, i+1)
            else
              s_min_L = prim_1d(k, i) - prim_1d(k, i-1)
              s_max_L = prim_1d(k, i+1) - prim_1d(k, i)
              q_L(k, i) = prim_1d(k, i) + 0.5d0 * tvd_slope(s_min_L, s_max_L, tvd_limiter_id)

              s_min_R = prim_1d(k, i+1) - prim_1d(k, i)
              s_max_R = prim_1d(k, i+2) - prim_1d(k, i+1)
              q_R(k, i) = prim_1d(k, i+1) - 0.5d0 * tvd_slope(s_min_R, s_max_R, tvd_limiter_id)
            end if
          end do
        end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(k,i,s_min_L,s_max_L,s_min_R,s_max_R)
        do k = 1, neq
          do i = 0, nodes
            s_min_L = prim_1d(k, i) - prim_1d(k, i-1)
            s_max_L = prim_1d(k, i+1) - prim_1d(k, i)
            q_L(k, i) = prim_1d(k, i) + 0.5d0 * tvd_slope(s_min_L, s_max_L, tvd_limiter_id)

            s_min_R = prim_1d(k, i+1) - prim_1d(k, i)
            s_max_R = prim_1d(k, i+2) - prim_1d(k, i+1)
            q_R(k, i) = prim_1d(k, i+1) - 0.5d0 * tvd_slope(s_min_R, s_max_R, tvd_limiter_id)
          end do
        end do
      !$OMP END PARALLEL DO
    end if


    case(REC_WENO3) 

    if (use_shock_sensor) then
      !$OMP PARALLEL DO PRIVATE(k, i, s_min_L, s_max_L, s_min_R, s_max_R, shock_sensor, p_jump, rho_jump)
        do i = 0, nodes
          
          p_jump = abs(prim_1d(eq_pr, i+1) - prim_1d(eq_pr, i-1)) / &
                  max(min(prim_1d(eq_pr, i+1), prim_1d(eq_pr, i-1)), p_floor)

          rho_jump = abs(prim_1d(eq_de, i+1) - prim_1d(eq_de, i-1)) / &
                    max(min(prim_1d(eq_de, i+1), prim_1d(eq_de, i-1)), rho_floor)

          if (p_jump > p_tol_tvd .or. rho_jump > rho_tol_tvd) then
              shock_sensor = .true.
          else
              shock_sensor = .false.
          end if

          do k = 1, neq
            
            if (prim_1d(eq_de, i) <= atm_factor * rho_floor .or. prim_1d(eq_de, i+1) <= atm_factor * rho_floor) then
              q_L(k, i) = prim_1d(k, i)
              q_R(k, i) = prim_1d(k, i+1)

            else if (shock_sensor) then
              s_min_L = prim_1d(k, i) - prim_1d(k, i-1)
              s_max_L = prim_1d(k, i+1) - prim_1d(k, i)
              q_L(k, i) = prim_1d(k, i) + 0.5d0 * tvd_slope(s_min_L, s_max_L, LIM_MC)

              s_min_R = prim_1d(k, i+1) - prim_1d(k, i)
              s_max_R = prim_1d(k, i+2) - prim_1d(k, i+1)
              q_R(k, i) = prim_1d(k, i+1) - 0.5d0 * tvd_slope(s_min_R, s_max_R, LIM_MC)

            else
              q_L(k, i) = weno3_interface(prim_1d(k, i-1 : i+1))
              q_R(k, i) = weno3_interface(prim_1d(k, i+2 : i : -1))
            end if

          end do
        end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(k, i)
        do i = 0, nodes
          do k = 1, neq
            ! Fronteras Geométricas Estrictas (Ahora incluye protección polar)
            if (trim(metric_type) /= 'Minkowski' .and. sweep_dir == DIR_X .and. i <= 0) then
              q_L(k, i) = prim_1d(k, i)
              q_R(k, i) = prim_1d(k, i+1)
            else 
              ! MP5 Puro sin paracaídas (Aplica a WENO3 en este caso)
              q_L(k, i) = weno3_interface(prim_1d(k, i-1 : i+1))
              q_R(k, i) = weno3_interface(prim_1d(k, i+2 : i : -1))
            end if
          end do
        end do
      !$OMP END PARALLEL DO
    end if


    case(REC_MP5) 

    if (use_shock_sensor) then
      !$OMP PARALLEL DO PRIVATE(k, i, s_min_L, s_max_L, s_min_R, s_max_R, p_max, p_min, p_jump, rho_max, rho_min, rho_jump)
        do i = 0, nodes

          p_max = max(prim_1d(eq_pr, i-2), prim_1d(eq_pr, i-1), prim_1d(eq_pr, i), &
                      prim_1d(eq_pr, i+1), prim_1d(eq_pr, i+2))
          p_min = min(prim_1d(eq_pr, i-2), prim_1d(eq_pr, i-1), prim_1d(eq_pr, i), &
                      prim_1d(eq_pr, i+1), prim_1d(eq_pr, i+2))

          rho_max = max(prim_1d(eq_de, i-2), prim_1d(eq_de, i-1), prim_1d(eq_de, i), &
                        prim_1d(eq_de, i+1), prim_1d(eq_de, i+2))
          rho_min = min(prim_1d(eq_de, i-2), prim_1d(eq_de, i-1), prim_1d(eq_de, i), &
                        prim_1d(eq_de, i+1), prim_1d(eq_de, i+2))

          ! Evaluar el salto solo si hay material real
          if (p_max > 1.0d3 * p_floor) then
              p_jump = (p_max - p_min) / max(p_min, p_floor)
          else
              p_jump = 0.0d0
          end if

          if (rho_max > 1.0d3 * rho_floor) then
              rho_jump = (rho_max - rho_min) / max(rho_min, rho_floor)
          else
              rho_jump = 0.0d0
          end if

          do k = 1, neq

            if (trim(metric_type) /= 'Minkowski' .and. sweep_dir == DIR_X .and. i <= 0) then
              q_L(k, i) = prim_1d(k, i)
              q_R(k, i) = prim_1d(k, i+1)

            else if (prim_1d(eq_de, i) <= atm_factor * rho_floor .or. prim_1d(eq_de, i+1) <= atm_factor * rho_floor) then
              q_L(k, i) = prim_1d(k, i)
              q_R(k, i) = prim_1d(k, i+1)

            else if (p_jump >= p_tol_tvd .or. rho_jump >= rho_tol_tvd) then
              s_min_L = prim_1d(k, i) - prim_1d(k, i-1)
              s_max_L = prim_1d(k, i+1) - prim_1d(k, i)
              q_L(k, i) = prim_1d(k, i) + 0.5d0 * tvd_slope(s_min_L, s_max_L, LIM_MC)

              s_min_R = prim_1d(k, i+1) - prim_1d(k, i)
              s_max_R = prim_1d(k, i+2) - prim_1d(k, i+1)
              q_R(k, i) = prim_1d(k, i+1) - 0.5d0 * tvd_slope(s_min_R, s_max_R, LIM_MC)

            else if (p_jump >= p_tol_weno .or. rho_jump >= rho_tol_weno .or. &
                    (trim(metric_type) /= 'Minkowski' .and. sweep_dir == DIR_X .and. i <= 3) .or. &
                    (trim(metric_type) /= 'Minkowski' .and. sweep_dir == DIR_Y .and. (i <= 2 .or. i >= nodes - 2))) then
              q_L(k, i) = weno3_interface(prim_1d(k, i-1 : i+1))
              q_R(k, i) = weno3_interface(prim_1d(k, i+2 : i : -1))

            else 
              q_L(k, i) = mp5_interface(prim_1d(k, i-2 : i+2))
              q_R(k, i) = mp5_interface(prim_1d(k, i+3 : i-1 : -1))
            end if

          end do
        end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(k, i)
        do i = 0, nodes
          do k = 1, neq
            ! Fronteras Geométricas Estrictas
            if (trim(metric_type) /= 'Minkowski' .and. sweep_dir == DIR_X .and. i <= 0) then
              q_L(k, i) = prim_1d(k, i)
              q_R(k, i) = prim_1d(k, i+1)
            else 
              ! MP5 Puro sin paracaídas
              q_L(k, i) = mp5_interface(prim_1d(k, i-2 : i+2))
              q_R(k, i) = mp5_interface(prim_1d(k, i+3 : i-1 : -1))
            end if
          end do
        end do
      !$OMP END PARALLEL DO
    end if

    case(REC_WENO5) 

    if (use_shock_sensor) then
      !$OMP PARALLEL DO PRIVATE(k, i, s_min_L, s_max_L, s_min_R, s_max_R, p_max, p_min, p_jump, rho_max, rho_min, rho_jump)
        do i = 0, nodes

          p_max = max(prim_1d(eq_pr, i-2), prim_1d(eq_pr, i-1), prim_1d(eq_pr, i), &
                      prim_1d(eq_pr, i+1), prim_1d(eq_pr, i+2))
          p_min = min(prim_1d(eq_pr, i-2), prim_1d(eq_pr, i-1), prim_1d(eq_pr, i), &
                      prim_1d(eq_pr, i+1), prim_1d(eq_pr, i+2))

          rho_max = max(prim_1d(eq_de, i-2), prim_1d(eq_de, i-1), prim_1d(eq_de, i), &
                        prim_1d(eq_de, i+1), prim_1d(eq_de, i+2))
          rho_min = min(prim_1d(eq_de, i-2), prim_1d(eq_de, i-1), prim_1d(eq_de, i), &
                        prim_1d(eq_de, i+1), prim_1d(eq_de, i+2))

          ! Evaluar el salto solo si hay material real
          if (p_max > 1.0d3 * p_floor) then
              p_jump = (p_max - p_min) / max(p_min, p_floor)
          else
              p_jump = 0.0d0
          end if

          if (rho_max > 1.0d3 * rho_floor) then
              rho_jump = (rho_max - rho_min) / max(rho_min, rho_floor)
          else
              rho_jump = 0.0d0
          end if

          do k = 1, neq

            ! 1. Frontera Radial Interna Estricta (DIR_X)
            if (trim(metric_type) /= 'Minkowski' .and. sweep_dir == DIR_X .and. i <= 0) then
              q_L(k, i) = prim_1d(k, i)
              q_R(k, i) = prim_1d(k, i+1)

            ! 2. Atmósfera profunda absoluta
            else if (prim_1d(eq_de, i) <= atm_factor * rho_floor .or. prim_1d(eq_de, i+1) <= atm_factor * rho_floor) then
              q_L(k, i) = prim_1d(k, i)
              q_R(k, i) = prim_1d(k, i+1)

            ! 3. ESCUDO DE ACANTILADO (El antídoto para el ruido en r=6.5)
            ! Si la diferencia relativa es astronómica, hay una discontinuidad tipo vacío.
            ! Forzamos Godunov para matar el overshoot (Fenómeno de Gibbs).
            else if (p_jump >= p_tol_godunov .or. rho_jump >= rho_tol_godunov) then
              q_L(k, i) = prim_1d(k, i)
              q_R(k, i) = prim_1d(k, i+1)

            ! 4. Choques normales -> TVD MC
            else if (p_jump >= p_tol_tvd .or. rho_jump >= rho_tol_tvd) then
              s_min_L = prim_1d(k, i) - prim_1d(k, i-1)
              s_max_L = prim_1d(k, i+1) - prim_1d(k, i)
              q_L(k, i) = prim_1d(k, i) + 0.5d0 * tvd_slope(s_min_L, s_max_L, LIM_MC)

              s_min_R = prim_1d(k, i+1) - prim_1d(k, i)
              s_max_R = prim_1d(k, i+2) - prim_1d(k, i+1)
              q_R(k, i) = prim_1d(k, i+1) - 0.5d0 * tvd_slope(s_min_R, s_max_R, LIM_MC)

            ! 5. Transición suave -> WENO3
            else if (p_jump >= p_tol_weno .or. rho_jump >= rho_tol_weno) then
              q_L(k, i) = weno3_interface(prim_1d(k, i-1 : i+1))
              q_R(k, i) = weno3_interface(prim_1d(k, i+2 : i : -1))

            ! 6. Disco y zonas suaves -> WENO5
            else 
              q_L(k, i) = weno5_interface(prim_1d(k, i-2 : i+2))
              q_R(k, i) = weno5_interface(prim_1d(k, i+3 : i-1 : -1))
            end if

          end do

        end do
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(k, i)
        do i = 0, nodes
          do k = 1, neq
            if (trim(metric_type) /= 'Minkowski' .and. sweep_dir == DIR_X .and. i <= 0) then
              q_L(k, i) = prim_1d(k, i)
              q_R(k, i) = prim_1d(k, i+1)
            else 
              ! WENO5 Puro sin paracaídas
              q_L(k, i) = weno5_interface(prim_1d(k, i-2 : i+2))
              q_R(k, i) = weno5_interface(prim_1d(k, i+3 : i-1 : -1))
            end if
          end do
        end do
      !$OMP END PARALLEL DO
    end if

    end select

    ! Escudo de positividad absoluta
    !$OMP PARALLEL DO PRIVATE(i)
    do i = 0, nodes
      q_L(eq_de, i) = max(q_L(eq_de, i), rho_floor)
      q_R(eq_de, i) = max(q_R(eq_de, i), rho_floor)
      q_L(eq_pr, i) = max(q_L(eq_pr, i), p_floor)
      q_R(eq_pr, i) = max(q_R(eq_pr, i), p_floor)
    end do
    !$OMP END PARALLEL DO
    
  end subroutine reconstruct_1d_core

end module reconstruction