! =======================================================================================
! Módulo: conditions (Condiciones Iniciales y de Frontera)
! ---------------------------------------------------------------------------------------
! Este módulo define el "escenario" físico del experimento. Contiene las rutinas para
! inicializar la malla con soluciones analíticas exactas (Toro de Fishbone-Moncrief,
! Acreción de Michel, Tubos de choque) y las reglas topológicas que rigen los bordes
! del dominio computacional (Fronteras periódicas, sumideros relativistas, outflow).
! =======================================================================================
module conditions
  use variables
  use metrics
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: set_initial_conditions, set_boundary_conditions, &
    inject_pressure_noise, setup_michel_accretion_initial, &
    inject_density_mode, inject_density_noise, &
    prepare_fishbone_moncrief, evaluate_fishbone_moncrief_state, &
    evaluate_dust_accretion_state, prepare_case_boundary_data, &
    set_boundary_conditions_accelerator

contains

  ! Subrutina: set_initial_conditions
  ! Propósito: Actúa como el "Big Bang" de la simulación. Asigna a cada celda de la
  ! matriz de primitivas su estado termodinámico y cinemático en T=0.
  subroutine set_initial_conditions()
    integer :: i, j, k
    real*8  :: sinc_factor, sinc_argument
    real*8  :: dist_sq, exponent, gauss_prof, r_loc
    real*8  :: x_pos, z_pos, interface_distance, perturbation_envelope
    real*8  :: dust_state(neq)

    ! Limpieza de RAM: Previene que celdas fantasma contengan "basura" de ejecuciones previas
    p = 0.0d0

    ! Se soporta tanto el nombre clásico como el ID numérico en formato string
    select case(trim(case_name))

    case('Sod', '1') ! Tubo de Choque de Sod (Prueba térmica base)
      do i = -nghost, nx+nghost
        if (x(i) < sod_interface) then
          p(i,:,:,eq_de) = sod_rho_left
          p(i,:,:,eq_pr) = sod_pressure_left
        else
          p(i,:,:,eq_de) = sod_rho_right
          p(i,:,:,eq_pr) = sod_pressure_right
        end if
      end do

    case('Strong', '2') ! Blast Wave (Captura de choques extremos)
      do i = -nghost, nx+nghost
        if (x(i) < strong_interface) then
          p(i,:,:,eq_de) = strong_rho_left
          p(i,:,:,eq_pr) = strong_pressure_left
        else
          p(i,:,:,eq_de) = strong_rho_right
          p(i,:,:,eq_pr) = strong_pressure_right
        end if
      end do

    case('Shu', '3') ! Shu-Osher Relativista (Interacción Choque-Turbulencia)
      do i = -nghost, nx+nghost
        if (x(i) < shu_interface) then
          p(i,:,:,eq_de) = shu_rho_left
          p(i,:,:,eq_pr) = shu_pressure_left
        else
          p(i,:,:,eq_de) = shu_rho_right + &
                               shu_rho_amplitude * sin(shu_wave_number * x(i))
          p(i,:,:,eq_pr) = shu_pressure_right
        end if
      end do

    case('KHI', '4') ! Inestabilidad de Kelvin-Helmholtz (Transición a la turbulencia)
      print *, "-> Inicializando Kelvin-Helmholtz (Setup analitico)..."

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            x_pos = x(i)
            z_pos = z(k)

            ! Perfil Base Discontinuo (Sharp Interface)
            if (abs(z_pos) < khi_half_width) then
              p(i,j,k,eq_de) = khi_rho_inner
              p(i,j,k,eq_vx) = khi_vx_inner
            else
              p(i,j,k,eq_de) = khi_rho_outer
              p(i,j,k,eq_vx) = khi_vx_outer
            end if

            p(i,j,k,eq_pr) = khi_pressure
            p(i,j,k,eq_vy) = 0.0d0

            ! La perturbación debe ser transversal al flujo de cizalla. Se
            ! localiza suavemente alrededor de las interfaces z=+-half_width.
            interface_distance = abs(abs(z_pos) - khi_half_width)
            perturbation_envelope = exp(-0.5d0 * &
              (interface_distance/khi_perturbation_width)**2)
            p(i,j,k,eq_vz) = khi_perturbation_amplitude * &
              sin(khi_perturbation_wave_number*pi*x_pos) * perturbation_envelope
          end do
        end do
      end do

    case('Jet', '5') ! Jet Relativista Centrado para Depuración
      print *, "-> Depurando: Jet centrado en el dominio con Outflow total..."

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            ! Estado base: Ambiente denso y estático
            p(i,j,k,eq_de) = jet_ambient_density
            p(i,j,k,eq_pr) = jet_ambient_pressure
            p(i,j,k,eq_vz) = jet_ambient_velocity

            ! "Tobera" de fluido relativista en el origen
            if (x(i) <= jet_nozzle_radius .and. z(k) <= jet_nozzle_length) then
              p(i,j,k,eq_de) = jet_density
              p(i,j,k,eq_pr) = jet_pressure
              p(i,j,k,eq_vz) = jet_velocity
            end if
          end do
        end do
      end do

    case('Conv', '6') ! Test de Convergencia Espacial L1/L_inf
      print *, '-> Test de Convergencia inicializado'

      sinc_argument = 0.5d0 * advected_wave_number * pi * dx
      if (abs(sinc_argument) > epsilon(1.0d0)) then
        sinc_factor = sin(sinc_argument) / sinc_argument
      else
        sinc_factor = 1.0d0
      end if

      do i = 1, nx
        p(i,:,:,eq_de) = advected_wave_density + &
                            advected_wave_amplitude * sinc_factor * &
                            sin(advected_wave_number * pi * x(i))
      end do

      p(:,:,:,eq_vx) = advected_wave_speed
      p(:,:,:,eq_pr) = advected_wave_pressure

    case('MichelA', '7') ! Acreción de Michel (Caída esférica de gas en un Agujero Negro)
      print *, "-> Inicializando Acrecion de Michel..."
      call setup_michel_accretion_initial(p(1:nx+nghost,1:ny,1:nz,:), &
                                           michel_critical_radius, &
                                           michel_critical_density)
      call cache_michel_injector()
    case('Dust', '8') ! Caída libre estacionaria desde el reposo en el infinito.
      do i = 1, nx+nghost
        call evaluate_dust_accretion_state(x(i), dust_state)
        do k = -nghost, nz+nghost
          do j = -nghost, ny+nghost
            p(i,j,k,:) = dust_state
          end do
        end do
      end do
      do i = -nghost, 0
        p(i,:,:,:) = p(1,:,:,:)
      end do

    case('OffAxis', '9') ! Off-axis Blast Wave (Explosión asimétrica)
      do k = -nghost, nz+nghost
        do j = 1, ny
          do i = -nghost, nx+nghost
            if (use_log_r) then
              r_loc = exp(x(i))
            else
              r_loc = x(i)
            end if
            dist_sq = r_loc**2 + offaxis_radial_center**2 - &
                      2.0d0 * offaxis_radial_center * r_loc * &
                      cos(z(k) - offaxis_phi_center)
            exponent = -dist_sq / (offaxis_width**2)

            if (exponent < -40.0d0) then
              gauss_prof = 0.0d0
            else
              gauss_prof = exp(exponent)
            end if

            p(i,j,k,eq_de) = offaxis_background_density + &
                                 offaxis_density_amplitude * gauss_prof
            p(i,j,k,eq_pr) = offaxis_background_pressure + &
                                 offaxis_pressure_amplitude * gauss_prof
          end do
        end do
      end do

    case('FishMoncEqu', '10') ! Toro de Fishbone-Moncrief en corte ecuatorial (2D r-phi)
      print *, "-> Inicializando Toro de Fishbone-Moncrief (corte ecuatorial 2D)..."
      call setup_fishbone_moncrief_initial(p(1:nx,1:ny,1:nz,:), &
        fm_inner_radius * bh_mass, fm_pressure_max_radius * bh_mass, &
        fm_polytropic_constant)

    case('FishMoncSag', '11') ! Toro de Fishbone-Moncrief en corte sagital (2D r-theta)
      print *, "-> Inicializando Toro de Fishbone-Moncrief (corte sagital 2D)..."
      call setup_fishbone_moncrief_initial(p(1:nx,1:ny,1:nz,:), &
        fm_inner_radius * bh_mass, fm_pressure_max_radius * bh_mass, &
        fm_polytropic_constant)

    case('FishMonc3D', '12')
      print *, "-> Inicializando Toro de Fishbone-Moncrief completo (3D)..."
      call setup_fishbone_moncrief_initial(p(1:nx,1:ny,1:nz,:), &
        fm_inner_radius * bh_mass, fm_pressure_max_radius * bh_mass, &
        fm_polytropic_constant)

    case default
      print *, "CRITICAL ERROR: case_name no valido en initial_conditions."
      stop
    end select
  end subroutine set_initial_conditions


  subroutine prepare_case_boundary_data()
    ! Reconstruye datos auxiliares que no forman parte del estado conservado.
    ! En particular, el inyector de Michel debe existir también al reiniciar.
    select case(trim(case_name))
    case('MichelA', '7')
      call setup_michel_accretion_initial(p(1:nx+nghost,1:ny,1:nz,:), &
        michel_critical_radius, michel_critical_density)
      call cache_michel_injector()
    end select
  end subroutine prepare_case_boundary_data


  subroutine cache_michel_injector()
    integer :: i, j, k

    if (allocated(michel_injector)) then
      if (size(michel_injector,2) /= nghost .or. &
          size(michel_injector,3) /= ny .or. &
          size(michel_injector,4) /= nz) deallocate(michel_injector)
    end if
    if (.not. allocated(michel_injector)) allocate(michel_injector(nghost,ny,nz,neq))

    do k = 1, nz
      do j = 1, ny
        do i = 1, nghost
          michel_injector(i,j,k,:) = p(nx+i,j,k,:)
        end do
      end do
    end do
  end subroutine cache_michel_injector


  subroutine evaluate_dust_accretion_state(radial_coordinate, prim_state)
    implicit none
    real*8, intent(in) :: radial_coordinate
    real*8, intent(out) :: prim_state(neq)
    real*8 :: r_loc, radial_jacobian, free_fall_speed

    if (use_log_r) then
      r_loc = exp(radial_coordinate)
      radial_jacobian = r_loc
    else
      r_loc = radial_coordinate
      radial_jacobian = 1.0d0
    end if
    if (r_loc <= 0.0d0 .or. bh_mass <= 0.0d0 .or. &
        dust_accretion_constant >= 0.0d0) then
      write(*,*) 'CRITICAL ERROR: invalid stationary dust-accretion parameters.'
      error stop 1
    end if

    free_fall_speed = sqrt(2.0d0*bh_mass/r_loc)
    prim_state = 0.0d0
    prim_state(eq_de) = -dust_accretion_constant / (r_loc**2*free_fall_speed)
    ! El código evoluciona un gas ideal: p_floor aproxima el límite frío p=0
    ! sin entregar un estado inadmisible al recuperador de primitivas.
    prim_state(eq_pr) = p_floor
    prim_state(eq_vx) = -1.0d0 / &
      (sqrt(1.0d0 + r_loc/(2.0d0*bh_mass)) * &
       (1.0d0 + free_fall_speed + 2.0d0*bh_mass/r_loc))
    prim_state(eq_vx) = prim_state(eq_vx) / radial_jacobian
  end subroutine evaluate_dust_accretion_state


  ! ========================================================================
  ! ACRECIÓN DE MICHEL (Solución analítica esférica 1D mapeada a 3D)
  ! ========================================================================
  subroutine setup_michel_accretion_initial(prim_state, r_critical, rho_critical, amp)
    implicit none
    real*8, intent(in) :: r_critical, rho_critical
    real*8, intent(in), optional :: amp
    real*8, intent(out) :: prim_state(:,:,:,:)

    real*8 :: u_c, ut_c, ur, ux, r_loc
    real*8 :: rho, pres, K_poly, h_spec, rho_final
    real*8 :: C1, C2, Vc2
    real*8 :: f_res, df_res, eps_tol = 1.0d-6
    real*8 :: f_p, h_p, rho_p, p_p, ur_p, ut_p, ut_tmp
    real*8 :: f_m, h_m, rho_m, p_m, ur_m, ut_m
    real*8 :: W_lorentz, g(3,3), beta(3), alpha
    integer :: local_nx, local_ny, local_nz
    integer :: i, j, k, iter, ii, jj

    ! Variables para la solución cuadrática universal
    real*8 :: g_tt, g_tx, A_quad, B_quad, C_quad, u_t_contra

    local_nx = size(prim_state, 1)
    local_ny = size(prim_state, 2)
    local_nz = size(prim_state, 3)

    ! Constantes de la solución asintótica
    u_c = -sqrt(bh_mass/(2.0d0*r_critical))
    ut_c = sqrt(1.0d0 - 2.0d0*bh_mass/r_critical + u_c*u_c)
    Vc2 = (u_c/ut_c)**2
    pres = (Vc2 * rho_critical) / (adb_idx - Vc2*g1)
    K_poly = pres / rho_critical**adb_idx
    h_spec = 1.0d0 + g1 * pres/rho_critical

    C1 = r_critical**2 * rho_critical * u_c
    C2 = r_critical**2 * rho_critical * h_spec * u_c * ut_c

    do i = 1, local_nx
      ! Extracción segura de la coordenada
      if (use_log_r) then
        r_loc = exp(x(i))
      else
        r_loc = x(i)
      end if

      ur = -sqrt(2.0d0*bh_mass/r_loc)

      ! Buscador de raíces de Newton-Raphson para el fluido
      do iter = 1, 40
        ut_tmp = sqrt(1.0d0 - 2.0d0*bh_mass/r_loc + ur**2)
        rho = C1/(r_loc**2 * ur)
        pres = K_poly * rho**adb_idx
        h_spec = 1.0d0 + g1 * pres/rho
        f_res = h_spec*ut_tmp - C2/C1

        ur_p = ur + eps_tol
        ur_m = ur - eps_tol

        ut_p = sqrt(1.0d0 - 2.0d0*bh_mass/r_loc + ur_p**2)
        rho_p = C1/(r_loc**2 * ur_p)
        p_p = K_poly * rho_p**adb_idx
        h_p = 1.0d0 + g1 * p_p/rho_p
        f_p = h_p*ut_p - C2/C1

        ut_m = sqrt(1.0d0 - 2.0d0*bh_mass/r_loc + ur_m**2)
        rho_m = C1/(r_loc**2 * ur_m)
        p_m = K_poly * rho_m**adb_idx
        h_m = 1.0d0 + g1 * p_m/rho_m
        f_m = h_m*ut_m - C2/C1

        df_res = (f_p - f_m)/(2.0d0*eps_tol)

        ur = ur - f_res/df_res
        if(abs(f_res) < 1.0d-12) exit
      end do

      rho = C1/(r_loc**2 * ur)
      pres = K_poly * rho**adb_idx

      rho_final = rho
      if (present(amp)) rho_final = rho + amp * exp(-0.5d0*((r_loc - 30.0d0)/2.0d0)**2)

      do k = 1, local_nz
        do j = 1, local_ny

          ! 1. OBTENCIÓN DE MÉTRICA SEGURA PARA CELDAS FANTASMA
          ! (Ignoramos la caché local aquí para evitar el Out-Of-Bounds en i > nx)
          call calculate_metric(x(i), y(j), alpha=alpha, beta=beta, gamma=g)

          if (use_log_r) then
            ux = ur / r_loc
          else
            ux = ur
          end if

          ! =========================================================
          ! 2. SOLUCIÓN CUADRÁTICA UNIVERSAL GRHD PARA u^t
          ! Resolvemos exactamente g_mu_nu u^mu u^nu = -1
          ! =========================================================
          g_tt = -alpha**2
          do ii = 1, 3
            do jj = 1, 3
              g_tt = g_tt + g(ii,jj) * beta(ii) * beta(jj)
            end do
          end do

          g_tx = 0.0d0
          do jj = 1, 3
            g_tx = g_tx + g(1,jj) * beta(jj)
          end do

          A_quad = g_tt
          B_quad = 2.0d0 * g_tx * ux
          C_quad = g(1,1) * ux**2 + 1.0d0

          ! Extracción de la raíz físicamente válida (u^t > 0)
          u_t_contra = (-B_quad - sqrt(B_quad**2 - 4.0d0 * A_quad * C_quad)) / (2.0d0 * A_quad)
          W_lorentz = alpha * u_t_contra

          ! 3. ASIGNACIÓN FINAL AL ESTADO PRIMITIVO
          prim_state(i,j,k,eq_de) = rho_final
          prim_state(i,j,k,eq_pr) = pres
          prim_state(i,j,k,eq_vx) = ux / W_lorentz + beta(1) / alpha
          prim_state(i,j,k,eq_vy) = beta(2) / alpha
          prim_state(i,j,k,eq_vz) = beta(3) / alpha

        end do
      end do
    end do
  end subroutine setup_michel_accretion_initial

  ! ========================================================================
  ! TORO DE FISHBONE-MONCRIEF
  ! ========================================================================
  subroutine setup_fishbone_moncrief_initial(prim_state, r_in, r_max_dens, K_poly)
    implicit none
    real*8, intent(inout) :: prim_state(:,:,:,:)
    real*8, intent(in)    :: r_in, r_max_dens, K_poly
    integer :: local_nx, local_ny, local_nz
    integer :: i, j, k
    real*8  :: l_ang, W_in, r_loc
    real*8  :: alpha, beta(3), g(3,3)
    logical :: model_is_valid, cell_is_torus

    local_nx = size(prim_state, 1)
    local_ny = size(prim_state, 2)
    local_nz = size(prim_state, 3)

    ! La solución física se prepara una sola vez. La dimensionalidad del caso
    ! queda determinada únicamente por local_ny y local_nz: ny=1 fija el corte
    ! ecuatorial, nz=1 fija el corte sagital y ambos >1 producen el toro 3D.
    call prepare_fishbone_moncrief(r_in, r_max_dens, l_ang, W_in, model_is_valid)
    if (.not. model_is_valid) then
      write(*,*) 'CRITICAL ERROR: invalid Fishbone-Moncrief model parameters.'
      error stop 1
    end if

    !!$OMP PARALLEL DO PRIVATE(i, j, k, r_loc, alpha, beta, g, cell_is_torus)
    do k = 1, local_nz
      do j = 1, local_ny
        do i = 1, local_nx

          if (use_log_r) then
            r_loc = exp(x(i))
          else
            r_loc = x(i)
          end if

          ! Usamos la caché de la métrica
          alpha = alpha_c(i,j,k)
          beta(:) = beta_c(:,i,j,k)
          g(:,:) = gamma_c(:,:,i,j,k)

          call evaluate_fishbone_moncrief_state(r_loc, r_in, K_poly, l_ang, W_in, &
                                                alpha, beta, g, &
                                                prim_state(i,j,k,:), cell_is_torus)
        end do
      end do
    end do
    !!$OMP END PARALLEL DO

    print *, ">>> Toro FM configurado (Métrica Completa + Kerr Spin a=", a_spin, ")"
    print *, "    Momento Angular (l) =", l_ang, " W_in =", W_in

  end subroutine setup_fishbone_moncrief_initial


  ! Prepara las dos constantes compartidas por cualquier representación del toro:
  ! corte ecuatorial, corte sagital o dominio tridimensional completo.
  subroutine prepare_fishbone_moncrief(r_in, r_max_dens, l_ang, W_in, ok)
    implicit none
    real*8, intent(in) :: r_in, r_max_dens
    real*8, intent(out) :: l_ang, W_in
    logical, intent(out) :: ok
    real*8 :: alpha, beta(3), g(3,3)
    real*8 :: coord_in, sq_Mr, numerator_l, denominator_l
    real*8 :: g_tt, g_tphi, numerator_ut, denominator_ut, ut_sq

    ok = .false.
    l_ang = 0.0d0
    W_in = 0.0d0

    if (bh_mass <= 0.0d0 .or. r_in <= 0.0d0 .or. &
        r_max_dens <= r_in) return

    sq_Mr = sqrt(bh_mass * r_max_dens)
    numerator_l = sq_Mr * (r_max_dens**2 - 2.0d0*a_spin*sq_Mr + a_spin**2)
    denominator_l = r_max_dens**2 - 2.0d0*bh_mass*r_max_dens + a_spin*sq_Mr
    if (.not. ieee_is_finite(denominator_l) .or. &
        abs(denominator_l) <= tiny(1.0d0)) return
    l_ang = numerator_l / denominator_l
    if (.not. ieee_is_finite(l_ang)) return

    if (use_log_r) then
      coord_in = log(r_in)
    else
      coord_in = r_in
    end if
    call calculate_metric(coord_in, pi/2.0d0, alpha=alpha, beta=beta, gamma=g)
    call fishbone_metric_components(alpha, beta, g, g_tt, g_tphi)

    numerator_ut = g_tphi**2 - g_tt*g(3,3)
    denominator_ut = g(3,3) + 2.0d0*l_ang*g_tphi + l_ang**2*g_tt
    if (numerator_ut <= 0.0d0 .or. denominator_ut <= 0.0d0) return

    ut_sq = numerator_ut / denominator_ut
    if (.not. ieee_is_finite(ut_sq) .or. ut_sq <= 0.0d0) return
    W_in = 0.5d0 * log(ut_sq)
    ok = ieee_is_finite(W_in)
  end subroutine prepare_fishbone_moncrief


  ! Evalúa la misma solución analítica en una celda, con independencia de cuántas
  ! direcciones tenga activas la malla. La métrica se recibe desde la caché para no
  ! recalcularla durante la inicialización de dominios 3D.
  subroutine evaluate_fishbone_moncrief_state(r_loc, r_in, K_poly, l_ang, W_in, &
                                              alpha, beta, g, prim_cell, is_torus)
    implicit none
    real*8, intent(in) :: r_loc, r_in, K_poly, l_ang, W_in
    real*8, intent(in) :: alpha, beta(3), g(3,3)
    real*8, intent(out) :: prim_cell(:)
    logical, intent(out) :: is_torus
    real*8 :: g_tt, g_tphi, numerator_ut, denominator_ut, ut_sq
    real*8 :: W_r, enthalpy, density, pressure, Omega, omega_denominator

    prim_cell = 0.0d0
    prim_cell(eq_de) = rho_floor
    prim_cell(eq_pr) = p_floor
    is_torus = .false.

    if (r_loc < r_in .or. K_poly <= 0.0d0 .or. alpha <= 0.0d0) return

    call fishbone_metric_components(alpha, beta, g, g_tt, g_tphi)
    numerator_ut = g_tphi**2 - g_tt*g(3,3)
    denominator_ut = g(3,3) + 2.0d0*l_ang*g_tphi + l_ang**2*g_tt
    if (numerator_ut <= 0.0d0 .or. denominator_ut <= 0.0d0) return

    ut_sq = numerator_ut / denominator_ut
    if (.not. ieee_is_finite(ut_sq) .or. ut_sq <= 0.0d0) return
    W_r = 0.5d0 * log(ut_sq)
    enthalpy = exp(W_in - W_r)
    if (.not. ieee_is_finite(enthalpy) .or. enthalpy <= 1.0d0) return

    density = ((enthalpy - 1.0d0) * (adb_idx - 1.0d0) / &
               (K_poly * adb_idx))**(1.0d0 / (adb_idx - 1.0d0))
    pressure = K_poly * density**adb_idx
    ! El problema discreto usa una atmósfera finita, no vacío. Una celda
    ! analítica por debajo de cualquiera de los pisos debe conservar el estado
    ! atmosférico completo inicializado arriba (incluidas velocidades nulas),
    ! igual que hará posteriormente la recuperación conservativa-primitiva.
    if (.not. ieee_is_finite(density) .or. .not. ieee_is_finite(pressure) .or. &
        density <= rho_floor .or. pressure <= p_floor) return

    omega_denominator = g(3,3) + l_ang*g_tphi
    if (.not. ieee_is_finite(omega_denominator) .or. &
        abs(omega_denominator) <= tiny(1.0d0)) return
    Omega = -(g_tphi + l_ang*g_tt) / omega_denominator
    if (.not. ieee_is_finite(Omega)) return

    prim_cell(eq_de) = density
    prim_cell(eq_pr) = pressure
    prim_cell(eq_vx) = beta(1) / alpha
    prim_cell(eq_vy) = beta(2) / alpha
    prim_cell(eq_vz) = (Omega + beta(3)) / alpha
    if (.not. all(ieee_is_finite(prim_cell))) then
      prim_cell = 0.0d0
      prim_cell(eq_de) = rho_floor
      prim_cell(eq_pr) = p_floor
      return
    end if
    is_torus = .true.
  end subroutine evaluate_fishbone_moncrief_state


  subroutine fishbone_metric_components(alpha, beta, g, g_tt, g_tphi)
    implicit none
    real*8, intent(in) :: alpha, beta(3), g(3,3)
    real*8, intent(out) :: g_tt, g_tphi

    g_tt = -alpha**2 + dot_product(beta, matmul(g, beta))
    g_tphi = dot_product(g(3,:), beta)
  end subroutine fishbone_metric_components

  ! Inicializa de manera portable y reproducible el generador intrínseco.
  subroutine set_reproducible_random_seed(seed_value)
    implicit none
    integer, intent(in) :: seed_value
    integer :: seed_size, idx
    integer, allocatable :: seed_values(:)

    call random_seed(size=seed_size)
    allocate(seed_values(seed_size))
    do idx = 1, seed_size
      seed_values(idx) = modulo(seed_value + 104729 * (idx - 1), huge(1) - 1)
      if (seed_values(idx) <= 0) seed_values(idx) = idx
    end do
    call random_seed(put=seed_values)
    deallocate(seed_values)
  end subroutine set_reproducible_random_seed

  ! Ruido blanco multiplicativo en presión. Es la perturbación PPI predeterminada
  ! porque excita simultáneamente varios modos sin imponer uno dominante.
  subroutine inject_pressure_noise(prim_state, amplitude, seed_value)
    implicit none
    real*8, intent(inout) :: prim_state(-nghost:nx+nghost, &
                                        -nghost:ny+nghost, -nghost:nz+nghost, neq)
    real*8, intent(in) :: amplitude
    integer, intent(in) :: seed_value
    integer :: i, j, k
    real*8  :: rand_val, rho_max_local, factor

    if (amplitude < 0.0d0 .or. amplitude >= 1.0d0) then
      write(*,*) 'Invalid pressure-noise amplitude; require 0 <= A < 1.'
      error stop 1
    end if

    print *, ">>> Inyectando ruido blanco en presion. A =", amplitude, &
             " semilla =", seed_value

    rho_max_local = maxval(prim_state(1:nx,1:ny,1:nz,eq_de))
    call set_reproducible_random_seed(seed_value)

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! Solo perturbar el cuerpo del disco
          if (prim_state(i,j,k,eq_de) > 0.05d0 * rho_max_local) then
            call random_number(rand_val)
            rand_val = 2.0d0 * rand_val - 1.0d0 ! Entre -1 y 1
            factor = 1.0d0 + amplitude * rand_val
            prim_state(i,j,k,eq_pr) = prim_state(i,j,k,eq_pr) * factor
          end if
        end do
      end do
    end do
  end subroutine inject_pressure_noise

  ! Subrutina: inject_density_mode
  ! Inyecta un modo de Fourier específico para estudios controlados
  ! subroutine inject_density_mode(prim_state, mode_m)
  !   implicit none
  !   ! Cambiamos a la dimensionalidad 4D estándar
  !   real*8, intent(inout) :: prim_state(:,:,:,:)
  !   real*8, intent(in)    :: mode_m ! Ahora puedes elegir si es 2.0d0 o 3.0d0 desde el programa principal
  !   integer :: i, j, k
  !   integer :: local_nx, local_ny, local_nz
  !   real*8  :: phi_loc, rho_max_local, envelope, factor

  !   local_nx = size(prim_state, 2)
  !   local_ny = size(prim_state, 3)
  !   local_nz = size(prim_state, 4)

  !   print *, ">>> ATENCION: Inyectando modo azimutal ISENTROPICO m =", mode_m, " en la densidad..."

  !   ! Extraemos la densidad máxima de toda la malla
  !   rho_max_local = maxval(prim_state(:,:,:,eq_de))

  !   do k = 1, local_nz
  !     ! PRECAUCIÓN: Asegúrate de que z(k) es tu coordenada azimutal phi
  !     phi_loc = z(k)

  !     do j = 1, local_ny
  !       do i = 1, local_nx

  !         ! Blindaje del vacío: Solo perturbamos el corazón del toroide
  !         if (prim_state(i,j,k,eq_de) > 0.05d0 * rho_max_local) then

  !           envelope = cos(mode_m * phi_loc)

  !           ! Factor de perturbación del 1%
  !           factor = 1.0d0 + 0.01d0 * envelope

  !           ! 1. Perturbar la densidad
  !           prim_state(i,j,k,eq_de) = prim_state(i,j,k,eq_de) * factor

  !           ! 2. Recalcular la presión isentrópica (CORRECCIÓN CRÍTICA)
  !           ! En lugar de hardcodear K = 0.015, usamos la relación P_nueva = P_vieja * (Rho_nueva / Rho_vieja)^Gamma
  !           ! Esto hace que la subrutina sea universal y a prueba de errores.
  !           prim_state(i,j,k,eq_pr) = prim_state(i,j,k,eq_pr) * (factor**adb_idx)

  !         end if
  !       end do
  !     end do
  !   end do
  ! end subroutine inject_density_mode

  subroutine inject_density_mode(prim_state, mode_m, amplitude)
    implicit none
    real*8, intent(inout) :: prim_state(-nghost:nx+nghost, &
                                        -nghost:ny+nghost, -nghost:nz+nghost, neq)
    real*8, intent(in)    :: mode_m
    real*8, intent(in)    :: amplitude
    integer :: i, j, k
    real*8  :: phi_loc, rho_max_local, envelope, factor

    if (amplitude < 0.0d0 .or. amplitude >= 1.0d0 .or. mode_m < 1.0d0) then
      write(*,*) 'Invalid coherent density mode; require 0 <= A < 1 and m >= 1.'
      error stop 1
    end if

    print *, ">>> Inyectando modo isentropico puro m =", mode_m, " A =", amplitude

    rho_max_local = maxval(prim_state(1:nx,1:ny,1:nz,eq_de))

    do k = 1, nz
      phi_loc = z(k)

      do j = 1, ny
        do i = 1, nx

          ! Blindaje del vacío: Solo perturbamos el corazón del toroide (> 5% de la densidad máxima)
          if (prim_state(i,j,k,eq_de) > 0.05d0 * rho_max_local) then

            envelope = amplitude * cos(mode_m * phi_loc)

            factor = 1.0d0 + envelope

            ! 3. Aplicar a la densidad
            prim_state(i,j,k,eq_de) = prim_state(i,j,k,eq_de) * factor

            ! 4. Recalcular presión isentrópica rigurosamente
            prim_state(i,j,k,eq_pr) = prim_state(i,j,k,eq_pr) * (factor**adb_idx)

          end if
        end do
      end do
    end do
  end subroutine inject_density_mode

  subroutine inject_density_noise(prim_state, amplitude, seed_value)
    implicit none
    real*8, intent(inout) :: prim_state(-nghost:nx+nghost, &
                                        -nghost:ny+nghost, -nghost:nz+nghost, neq)
    real*8, intent(in) :: amplitude
    integer, intent(in) :: seed_value
    integer :: i, j, k
    real*8  :: rho_max_local, factor, rand_val

    if (amplitude < 0.0d0 .or. amplitude >= 1.0d0) then
      write(*,*) 'Invalid density-noise amplitude; require 0 <= A < 1.'
      error stop 1
    end if

    print *, ">>> Inyectando ruido blanco isentropico. A =", amplitude, &
             " semilla =", seed_value

    rho_max_local = maxval(prim_state(1:nx,1:ny,1:nz,eq_de))
    call set_reproducible_random_seed(seed_value)

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! Solo perturbar el cuerpo del disco
          if (prim_state(i,j,k,eq_de) > 0.05d0 * rho_max_local) then

            call random_number(rand_val)
            rand_val = 2.0d0 * rand_val - 1.0d0 ! Entre -1 y 1

            factor = 1.0d0 + amplitude * rand_val

            prim_state(i,j,k,eq_de) = prim_state(i,j,k,eq_de) * factor
            prim_state(i,j,k,eq_pr) = prim_state(i,j,k,eq_pr) * (factor**adb_idx)

          end if
        end do
      end do
    end do
  end subroutine inject_density_noise


  ! Subrutina: set_boundary_conditions
  ! Propósito: Rellena las "celdas fantasma" preservando la topología del dominio en 3D.
  subroutine set_boundary_conditions(q_state)
    implicit none
    real*8, intent(inout) :: q_state(-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost,neq)
    integer :: i, j, k, g, k_shift
    real*8 :: dust_state(neq)

    select case(trim(case_name))
    case('OffAxis', '9', 'FishMoncEqu', '10', 'FishMoncSag', '11', 'FishMonc3D', '12')
      if (.not. allocated(alpha_c) .or. .not. allocated(beta_c)) then
        write(*,*) 'CRITICAL ERROR: radial outflow boundaries require the metric cache.'
        error stop 1
      end if
      if (size(alpha_c,1) < nx .or. size(alpha_c,2) < ny .or. size(alpha_c,3) < nz .or. &
          size(beta_c,2) < nx .or. size(beta_c,3) < ny .or. size(beta_c,4) < nz) then
        write(*,*) 'CRITICAL ERROR: metric cache dimensions do not match the fluid grid.'
        error stop 1
      end if
    end select

    select case(trim(case_name))

    case('Sod', '1', 'Strong', '2', 'Shu', '3') ! Cajas 1D - Gradiente Cero (Outflow libre)
      do g = 1, nghost
        q_state(1-g,:,:,:) = q_state(1,:,:,:)
        q_state(nx+g,:,:,:) = q_state(nx,:,:,:)
      end do

      do g = 1, nghost
        q_state(:,1-g,:,:) = q_state(:,ny + 1 - g,:,:)
        q_state(:,ny+g,:,:) = q_state(:,g,:,:)

        q_state(:,:,1-g,:) = q_state(:,:,nz + 1 - g,:)
        q_state(:,:,nz+g,:) = q_state(:,:,g,:)
      end do

    case('KHI', '4', 'Conv', '6') ! Fronteras periódicas
      do g = 1, nghost
        q_state(1-g,:,:,:) = q_state(nx + 1 - g,:,:,:)
        q_state(nx+g,:,:,:) = q_state(g,:,:,:)

        q_state(:,1-g,:,:) = q_state(:,ny + 1 - g,:,:)
        q_state(:,ny+g,:,:) = q_state(:,g,:,:)

        q_state(:,:,1-g,:) = q_state(:,:,nz + 1 - g,:)
        q_state(:,:,nz+g,:) = q_state(:,:,g,:)
      end do

    case('Jet', '5') ! Slab Jet Cartesiano

      ! Reseteamos fantasmas a ambiente
      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = -nghost, nz+nghost
        do j = -nghost, ny+nghost
          do i = -nghost, nx+nghost
            if (i < 1 .or. i > nx .or. j < 1 .or. j > ny .or. k < 1 .or. k > nz) then
                q_state(i,j,k,eq_de) = jet_ambient_density
                q_state(i,j,k,eq_pr) = jet_ambient_pressure
                q_state(i,j,k,eq_vx:eq_vz) = 0.0d0
                q_state(i,j,k,eq_vz) = jet_ambient_velocity
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! 1. X -> REFLECTANTES
      !$OMP PARALLEL DO PRIVATE(j, k, g)
      do k = -nghost, nz+nghost
        do j = -nghost, ny+nghost
          do g = 1, nghost
            q_state(1-g,j,k,:) = q_state(g,j,k,:)
            q_state(1-g,j,k,eq_vx) = -q_state(g,j,k,eq_vx)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! 2. Z -> TOBERA + OUTFLOW
      !$OMP PARALLEL DO PRIVATE(i, j, g)
      do j = -nghost, ny+nghost
        do i = -nghost, nx+nghost
          do g = 1, nghost
            q_state(i,j,nz+g,:) = q_state(i,j,nz,:)
            q_state(i,j,nz+g,eq_vz) = max(0.0d0, q_state(i,j,nz,eq_vz))

            if (x(i) <= jet_nozzle_radius) then
              q_state(i,j,1-g,eq_de) = jet_density
              q_state(i,j,1-g,eq_pr) = jet_pressure
              q_state(i,j,1-g,eq_vx:eq_vy) = 0.0d0
              q_state(i,j,1-g,eq_vz) = jet_velocity
            else
              q_state(i,j,1-g,:) = q_state(i,j,g,:)
              q_state(i,j,1-g,eq_vz) = -q_state(i,j,g,eq_vz)
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

    case('MichelA', '7')
      ! 1. Radial (X)
      !$OMP PARALLEL DO PRIVATE(j, k, g)
      do k = 1, nz
        do j = 1, ny
          do g = 1, nghost
            q_state(1-g,j,k,:) = q_state(1,j,k,:)
            q_state(nx+g,j,k,:) = michel_injector(g,j,k,:)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! 2. Polar (Y) y Azimutal (Z)
      do g = 1, nghost
        q_state(:,1-g,:,:) = q_state(:,ny + 1 - g,:,:)
        q_state(:,ny+g,:,:) = q_state(:,g,:,:)
        q_state(:,:,1-g,:) = q_state(:,:,nz + 1 - g,:)
        q_state(:,:,nz+g,:) = q_state(:,:,g,:)
      end do

    case('Dust', '8')
      ! El borde interno es de salida hacia el agujero negro. En el borde
      ! externo se mantiene la solución analítica estacionaria como inyector.
      do g = 1, nghost
        call evaluate_dust_accretion_state(x(nx+g), dust_state)
        do k = 1, nz
          do j = 1, ny
            q_state(1-g,j,k,:) = q_state(1,j,k,:)
            q_state(1-g,j,k,eq_de) = max(q_state(1,j,k,eq_de), rho_floor)
            q_state(1-g,j,k,eq_pr) = max(q_state(1,j,k,eq_pr), p_floor)
            q_state(1-g,j,k,eq_vx) = min(0.0d0, q_state(1,j,k,eq_vx))
            q_state(nx+g,j,k,:) = dust_state
          end do
        end do
      end do

      do g = 1, nghost
        q_state(:,1-g,:,:) = q_state(:,ny + 1 - g,:,:)
        q_state(:,ny+g,:,:) = q_state(:,g,:,:)
        q_state(:,:,1-g,:) = q_state(:,:,nz + 1 - g,:)
        q_state(:,:,nz+g,:) = q_state(:,:,g,:)
      end do

    case('OffAxis', '9', 'FishMoncEqu', '10')

      ! 1. Radial
      !$OMP PARALLEL DO PRIVATE(j, k, g)
      do k = 1, nz
        do j = 1, ny
          do g = 1, nghost
            q_state(1-g,j,k,:) = q_state(1,j,k,:)
            q_state(1-g,j,k,eq_de) = max(q_state(1,j,k,eq_de), rho_floor)
            q_state(1-g,j,k,eq_pr) = max(q_state(1,j,k,eq_pr), p_floor)
            q_state(1-g,j,k,eq_vx) = min(beta_c(1,1,j,k)/alpha_c(1,j,k), &
                                                     q_state(1,j,k,eq_vx))
            if (q_state(1,j,k,eq_de) < 1.0d-6) q_state(1-g,j,k,eq_vy:eq_vz) = 0.0d0

            q_state(nx+g,j,k,eq_de) = max(q_state(nx,j,k,eq_de), rho_floor)
            q_state(nx+g,j,k,eq_pr) = max(q_state(nx,j,k,eq_pr), p_floor)
            q_state(nx+g,j,k,eq_vx) = max(beta_c(1,nx,j,k)/alpha_c(nx,j,k), &
                                                      q_state(nx,j,k,eq_vx))
            q_state(nx+g,j,k,eq_vy) = q_state(nx,j,k,eq_vy)
            q_state(nx+g,j,k,eq_vz) = q_state(nx,j,k,eq_vz)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! 2. Transversales
      do g = 1, nghost
        q_state(:,1-g,:,:) = q_state(:,ny + 1 - g,:,:)
        q_state(:,ny+g,:,:) = q_state(:,g,:,:)
        q_state(:,:,1-g,:) = q_state(:,:,nz + 1 - g,:)
        q_state(:,:,nz+g,:) = q_state(:,:,g,:)
      end do

    case('FishMoncSag', '11')

      ! 1. Llenado Radial (X) - Físico
      !$OMP PARALLEL DO PRIVATE(j, k, g)
      do k = 1, nz
        do j = 1, ny
          do g = 1, nghost
            ! --- FRONTERA INTERNA (Agujero Negro) ---
            q_state(1-g,j,k,:) = q_state(1,j,k,:)
            q_state(1-g,j,k,eq_de) = max(q_state(1,j,k,eq_de), rho_floor)
            q_state(1-g,j,k,eq_pr) = max(q_state(1,j,k,eq_pr), p_floor)
            ! No entra fluido desde la excisión: alpha*v^r-beta^r <= 0.
            q_state(1-g,j,k,eq_vx) = min(beta_c(1,1,j,k)/alpha_c(1,j,k), &
                                                     q_state(1,j,k,eq_vx))

            ! --- FRONTERA EXTERNA ---
            q_state(nx+g,j,k,:) = q_state(nx,j,k,:)
            q_state(nx+g,j,k,eq_de) = max(q_state(nx,j,k,eq_de), rho_floor)
            q_state(nx+g,j,k,eq_pr) = max(q_state(nx,j,k,eq_pr), p_floor)
            ! No entra fluido desde el exterior: alpha*v^r-beta^r >= 0.
            q_state(nx+g,j,k,eq_vx) = max(beta_c(1,nx,j,k)/alpha_c(nx,j,k), &
                                                      q_state(nx,j,k,eq_vx))
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! 2. Llenado Polar (Y) - Sobre X Expandido
      !$OMP PARALLEL DO PRIVATE(i, k, g)
      do k = 1, nz
        do i = -nghost, nx+nghost
          do g = 1, nghost
            q_state(i,1-g,k,:) = q_state(i,g,k,:)
            q_state(i,1-g,k,eq_vy) = -q_state(i,g,k,eq_vy)
            q_state(i,ny+g,k,:) = q_state(i,ny+1-g,k,:)
            q_state(i,ny+g,k,eq_vy) = -q_state(i,ny+1-g,k,eq_vy)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! 3. Llenado Azimutal (Z) - Global
      do g = 1, nghost
        q_state(:,:,1-g,:) = q_state(:,:,nz + 1 - g,:)
        q_state(:,:,nz+g,:) = q_state(:,:,g,:)
      end do

    case('FishMonc3D', '12')

      ! Escudo de Paridad
      if (mod(nz, 2) /= 0) then
         print *, "CRITICAL ERROR: nz must be even for 3D Trans-Polar boundary conditions."
         stop
      end if

      ! ========================================================================
      ! 1. Llenado Radial (X) - Sobre dominio Físico Y, Z
      ! ========================================================================
      !$OMP PARALLEL DO PRIVATE(j, k, g)
      do k = 1, nz
        do j = 1, ny
          do g = 1, nghost
            ! Interna (Horizonte de EF)
            q_state(1-g,j,k,:) = q_state(1,j,k,:)
            q_state(1-g,j,k,eq_de) = max(q_state(1,j,k,eq_de), rho_floor)
            q_state(1-g,j,k,eq_pr) = max(q_state(1,j,k,eq_pr), p_floor)
            q_state(1-g,j,k,eq_vx) = min(beta_c(1,1,j,k)/alpha_c(1,j,k), &
                                                     q_state(1,j,k,eq_vx))

            ! Externa (Outflow)
            q_state(nx+g,j,k,:) = q_state(nx,j,k,:)
            q_state(nx+g,j,k,eq_de) = max(q_state(nx,j,k,eq_de), rho_floor)
            q_state(nx+g,j,k,eq_pr) = max(q_state(nx,j,k,eq_pr), p_floor)
            q_state(nx+g,j,k,eq_vx) = max(beta_c(1,nx,j,k)/alpha_c(nx,j,k), &
                                                      q_state(nx,j,k,eq_vx))
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! ========================================================================
      ! 2. Llenado Polar (Y): Pole-Crossing - Sobre dominio extendido X (-nghost:nx+nghost)
      ! ========================================================================
      !$OMP PARALLEL DO PRIVATE(i, k, g, k_shift)
      do k = 1, nz
        k_shift = k + (nz / 2)
        if (k_shift > nz) k_shift = k_shift - nz

        do i = -nghost, nx+nghost
          do g = 1, nghost
            ! Polo Norte
            q_state(i,1-g,k,:) = q_state(i,g,k_shift,:)
            q_state(i,1-g,k,eq_vy) = -q_state(i,g,k_shift,eq_vy)
            q_state(i,1-g,k,eq_vz) =  q_state(i,g,k_shift,eq_vz)

            ! Polo Sur
            q_state(i,ny+g,k,:) = q_state(i,ny+1-g,k_shift,:)
            q_state(i,ny+g,k,eq_vy) = -q_state(i,ny+1-g,k_shift,eq_vy)
            q_state(i,ny+g,k,eq_vz) =  q_state(i,ny+1-g,k_shift,eq_vz)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! ========================================================================
      ! 3. Llenado Azimutal (Z): Periódico - Sobre matriz completa
      ! ========================================================================
      do g = 1, nghost
        q_state(:,:,1-g,:) = q_state(:,:,nz + 1 - g,:)
        q_state(:,:,nz+g,:)  = q_state(:,:,g,:)
      end do

    case default
      write(*,*) 'CRITICAL ERROR: case_name no valido en set_boundary_conditions.'
      error stop 1

    end select
  end subroutine set_boundary_conditions

  ! Condiciones de frontera sin secciones vectoriales de Fortran. Cada fase se
  ! ejecuta en el dispositivo y respeta el mismo orden topológico que la ruta CPU.
  subroutine set_boundary_conditions_accelerator(q_state)
    real*8, intent(inout) :: q_state(-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost,neq)
    integer :: i,j,k,g,m,k_shift

    select case(case_id)
    case(1,2,3)
      !$acc parallel loop collapse(3) gang vector present(q_state)
      do k=1,nz; do j=1,ny; do g=1,nghost
        do m=1,neq
          q_state(1-g,j,k,m)=q_state(1,j,k,m)
          q_state(nx+g,j,k,m)=q_state(nx,j,k,m)
        end do
      end do; end do; end do
      !$acc end parallel loop
    case(4,6)
      !$acc parallel loop collapse(3) gang vector present(q_state)
      do k=1,nz; do j=1,ny; do g=1,nghost
        do m=1,neq
          q_state(1-g,j,k,m)=q_state(nx+1-g,j,k,m)
          q_state(nx+g,j,k,m)=q_state(g,j,k,m)
        end do
      end do; end do; end do
      !$acc end parallel loop
    case(5)
      !$acc parallel loop collapse(3) gang vector present(q_state)
      do k=-nghost,nz+nghost; do j=-nghost,ny+nghost; do i=-nghost,nx+nghost
        if (i < 1 .or. i > nx .or. j < 1 .or. j > ny .or. k < 1 .or. k > nz) then
          q_state(i,j,k,eq_de)=jet_ambient_density
          q_state(i,j,k,eq_pr)=jet_ambient_pressure
          q_state(i,j,k,eq_vx)=0.0d0
          q_state(i,j,k,eq_vy)=0.0d0
          q_state(i,j,k,eq_vz)=jet_ambient_velocity
        end if
      end do; end do; end do
      !$acc end parallel loop
      !$acc parallel loop collapse(3) gang vector present(q_state)
      do k=-nghost,nz+nghost; do j=-nghost,ny+nghost; do g=1,nghost
        do m=1,neq
          q_state(1-g,j,k,m)=q_state(g,j,k,m)
        end do
        q_state(1-g,j,k,eq_vx)=-q_state(g,j,k,eq_vx)
      end do; end do; end do
      !$acc end parallel loop
      !$acc parallel loop collapse(3) gang vector present(q_state,x)
      do j=-nghost,ny+nghost; do i=-nghost,nx+nghost; do g=1,nghost
        do m=1,neq
          q_state(i,j,nz+g,m)=q_state(i,j,nz,m)
        end do
        q_state(i,j,nz+g,eq_vz)=max(0.0d0,q_state(i,j,nz,eq_vz))
        if (x(i) <= jet_nozzle_radius) then
          q_state(i,j,1-g,eq_de)=jet_density
          q_state(i,j,1-g,eq_pr)=jet_pressure
          q_state(i,j,1-g,eq_vx)=0.0d0
          q_state(i,j,1-g,eq_vy)=0.0d0
          q_state(i,j,1-g,eq_vz)=jet_velocity
        else
          do m=1,neq
            q_state(i,j,1-g,m)=q_state(i,j,g,m)
          end do
          q_state(i,j,1-g,eq_vz)=-q_state(i,j,g,eq_vz)
        end if
      end do; end do; end do
      !$acc end parallel loop
      return
    case(7)
      !$acc parallel loop collapse(3) gang vector present(q_state,michel_injector)
      do k=1,nz; do j=1,ny; do g=1,nghost
        do m=1,neq
          q_state(1-g,j,k,m)=q_state(1,j,k,m)
          q_state(nx+g,j,k,m)=michel_injector(g,j,k,m)
        end do
      end do; end do; end do
      !$acc end parallel loop
    case(8)
      !$acc parallel loop collapse(3) gang vector present(q_state)
      do k=1,nz; do j=1,ny; do g=1,nghost
        do m=1,neq
          q_state(1-g,j,k,m)=q_state(1,j,k,m)
        end do
        q_state(1-g,j,k,eq_de)=max(q_state(1,j,k,eq_de),rho_floor)
        q_state(1-g,j,k,eq_pr)=max(q_state(1,j,k,eq_pr),p_floor)
        q_state(1-g,j,k,eq_vx)=min(0.0d0,q_state(1,j,k,eq_vx))
      end do; end do; end do
      !$acc end parallel loop
    case(9,10,11,12)
      !$acc parallel loop collapse(3) gang vector present(q_state,alpha_c,beta_c)
      do k=1,nz; do j=1,ny; do g=1,nghost
        do m=1,neq
          q_state(1-g,j,k,m)=q_state(1,j,k,m)
          q_state(nx+g,j,k,m)=q_state(nx,j,k,m)
        end do
        q_state(1-g,j,k,eq_de)=max(q_state(1,j,k,eq_de),rho_floor)
        q_state(1-g,j,k,eq_pr)=max(q_state(1,j,k,eq_pr),p_floor)
        q_state(1-g,j,k,eq_vx)=min(beta_c(1,1,j,k)/alpha_c(1,j,k),q_state(1,j,k,eq_vx))
        q_state(nx+g,j,k,eq_de)=max(q_state(nx,j,k,eq_de),rho_floor)
        q_state(nx+g,j,k,eq_pr)=max(q_state(nx,j,k,eq_pr),p_floor)
        q_state(nx+g,j,k,eq_vx)=max(beta_c(1,nx,j,k)/alpha_c(nx,j,k),q_state(nx,j,k,eq_vx))
        if ((case_id == 9 .or. case_id == 10) .and. q_state(1,j,k,eq_de) < 1.0d-6) then
          q_state(1-g,j,k,eq_vy)=0.0d0
          q_state(1-g,j,k,eq_vz)=0.0d0
        end if
      end do; end do; end do
      !$acc end parallel loop
    end select

    if (case_id == 11) then
      !$acc parallel loop collapse(3) gang vector present(q_state)
      do k=1,nz; do i=-nghost,nx+nghost; do g=1,nghost
        do m=1,neq
          q_state(i,1-g,k,m)=q_state(i,g,k,m)
          q_state(i,ny+g,k,m)=q_state(i,ny+1-g,k,m)
        end do
        q_state(i,1-g,k,eq_vy)=-q_state(i,g,k,eq_vy)
        q_state(i,ny+g,k,eq_vy)=-q_state(i,ny+1-g,k,eq_vy)
      end do; end do; end do
      !$acc end parallel loop
    else if (case_id == 12) then
      !$acc parallel loop collapse(3) gang vector present(q_state)
      do k=1,nz; do i=-nghost,nx+nghost; do g=1,nghost
        k_shift=k+nz/2
        if (k_shift > nz) k_shift=k_shift-nz
        do m=1,neq
          q_state(i,1-g,k,m)=q_state(i,g,k_shift,m)
          q_state(i,ny+g,k,m)=q_state(i,ny+1-g,k_shift,m)
        end do
        q_state(i,1-g,k,eq_vy)=-q_state(i,g,k_shift,eq_vy)
        q_state(i,ny+g,k,eq_vy)=-q_state(i,ny+1-g,k_shift,eq_vy)
      end do; end do; end do
      !$acc end parallel loop
    else
      !$acc parallel loop collapse(3) gang vector present(q_state)
      do k=1,nz; do i=-nghost,nx+nghost; do g=1,nghost
        do m=1,neq
          q_state(i,1-g,k,m)=q_state(i,ny+1-g,k,m)
          q_state(i,ny+g,k,m)=q_state(i,g,k,m)
        end do
      end do; end do; end do
      !$acc end parallel loop
    end if

    !$acc parallel loop collapse(3) gang vector present(q_state)
    do j=-nghost,ny+nghost; do i=-nghost,nx+nghost; do g=1,nghost
      do m=1,neq
        q_state(i,j,1-g,m)=q_state(i,j,nz+1-g,m)
        q_state(i,j,nz+g,m)=q_state(i,j,g,m)
      end do
    end do; end do; end do
    !$acc end parallel loop
  end subroutine set_boundary_conditions_accelerator

end module conditions
