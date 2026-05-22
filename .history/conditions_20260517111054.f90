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
  implicit none
  private
  public :: set_initial_conditions, set_boundary_conditions, &
    inject_pressure_noise, setup_michel_accretion_initial, &
    inject_density_mode, inject_density_noise

contains

  ! Subrutina: set_initial_conditions
  ! Propósito: Actúa como el "Big Bang" de la simulación. Asigna a cada celda de la 
  ! matriz de primitivas su estado termodinámico y cinemático en T=0.
  subroutine set_initial_conditions()
    integer :: i, j, k
    real*8  :: C_const = -0.5d0
    real*8  :: sinc_factor
    real*8  :: x_center, z_center, dist_sq, width, exponent, gauss_prof
    real*8  :: x_pos, z_pos

    ! --- Variables Físicas Locales (Inicializadas en 0 para evitar warnings) ---
    real*8  :: r_critical   = 0.0d0  ! Radio sónico para acreción de Michel
    real*8  :: rho_critical = 0.0d0  ! Densidad en el punto crítico (o infinito)
    real*8  :: r_in         = 0.0d0  ! Radio interno del disco/toro
    real*8  :: r_max_dens   = 0.0d0  ! Radio de densidad y presión máxima
    real*8  :: K_poly       = 0.0d0  ! Constante politrópica (P = K * rho^Gamma)

    ! Limpieza de RAM: Previene que celdas fantasma contengan "basura" de ejecuciones previas
    p = 0.0d0
    sinc_factor = sin(pi * dx) / (pi * dx)

    ! Se soporta tanto el nombre clásico como el ID numérico en formato string
    select case(trim(case_name)) 

    case('Sod', '1') ! Tubo de Choque de Sod (Prueba térmica base)
      do i = -nghost, nx+nghost
        if (x(i) < 0.5d0) then
          p(eq_de, i, :, :) = 1.0d0
          p(eq_pr, i, :, :) = 1.0d0
        else
          p(eq_de, i, :, :) = 0.125d0
          p(eq_pr, i, :, :) = 0.1d0
        end if
      end do

    case('Strong', '2') ! Blast Wave (Captura de choques extremos)
      do i = -nghost, nx+nghost
        p(eq_de, i, :, :) = 1.0d0
        if (x(i) < 0.5d0) then
          p(eq_pr, i, :, :) = 1000.0d0   
        else
          p(eq_pr, i, :, :) = 0.01d0
        end if
      end do

    case('Shu', '3') ! Shu-Osher Relativista (Interacción Choque-Turbulencia)
      do i = -nghost, nx+nghost
        if (x(i) < 0.5d0) then
          p(eq_de, i, :, :) = 5.0d0
          p(eq_pr, i, :, :) = 50.0d0
        else
          p(eq_de, i, :, :) = 2.0d0 + 0.3d0 * sin(50.0d0 * x(i))
          p(eq_pr, i, :, :) = 5.0d0
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
            if (abs(z_pos) < 0.25d0) then
              p(eq_de, i, j, k) = 2.0d0
              p(eq_vx, i, j, k) = 0.5d0
            else
              p(eq_de, i, j, k) = 1.0d0
              p(eq_vx, i, j, k) = -0.5d0
            end if
            
            p(eq_pr, i, j, k) = 2.5d0
            p(eq_vy, i, j, k) = 0.0d0 
            
            ! Perturbación de velocidad para detonar los vórtices
            p(eq_vx, i, j, k) = p(eq_vx, i, j, k) * (1.0d0 + 0.01d0 * cos(10.0d0 * pi * x_pos) * cos(10.0d0 * pi * z_pos))
          end do
        end do
      end do

    case('Jet', '5') ! Jet Relativista Centrado para Depuración
      print *, "-> Depurando: Jet centrado en el dominio con Outflow total..."
      
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            ! Estado base: Ambiente denso y estático
            p(eq_de, i, j, k) = 10.0d0
            p(eq_pr, i, j, k) = 0.01d0
            p(eq_vz, i, j, k) = 0.0d0
            
            ! "Tobera" de fluido relativista en el origen
            if (x(i) <= 1.0d0 .and. z(k) <= 1.0d0) then
              p(eq_de, i, j, k) = 0.1d0
              p(eq_vz, i, j, k) = 0.99d0
            end if
          end do
        end do
      end do

    case('Conv', '6') ! Test de Convergencia Espacial L1/L_inf
      print *, '-> Test de Convergencia inicializado'

      do i = 1, nx
        p(eq_de, i, :, :) = 1.0d0 + 0.1d0 * sinc_factor * sin(2.0d0 * pi * x(i))
      end do
      
      p(eq_vx, :, :, :) = 0.5d0
      p(eq_pr, :, :, :) = 1.0d0

    case('MichelA', '7') ! Acreción de Michel (Caída esférica de gas en un Agujero Negro)
      print *, "-> Inicializando Acrecion de Michel..."
      r_critical   = 100.0d0
      rho_critical = 1.0d-1
      
      call setup_michel_accretion_initial(p(:, 1:nx+nghost, 1:ny, 1:nz), r_critical, rho_critical)
      
      ! Guardamos el estado asintótico en el borde externo para usarlo como inyector
      if (.not. allocated(michel_injector)) allocate(michel_injector(neq, nghost, ny, nz))
        
      do k = 1, nz
        do j = 1, ny
          do i = 1, nghost
            michel_injector(:, i, j, k) = p(:, nx+i, j, k)
          end do
        end do
      end do

    case('Dust', '8') ! Dust accretion (Acreción de polvo en caída libre, Presión = 0)
      p(eq_de, :, :, :) = -C_const/(x_max**2 * sqrt(2.0d0*bh_mass/x_max))
      p(eq_vx, :, :, :) = -1.0d0/(sqrt(1.0d0+x_max/(2.0d0*bh_mass))*(1.0d0 + sqrt(2.0d0*bh_mass/x_max) + 2.0d0*bh_mass/x_max))

    case('OffAxis', '9') ! Off-axis Blast Wave (Explosión asimétrica)
      x_center = 15.0d0
      z_center = pi
      width = 1.5d0      

      do k = -nghost, nz+nghost
        do j = 1, ny
          do i = -nghost, nx+nghost
            dist_sq = x(i)**2 + x_center**2 - 2.0d0 * x_center * x(i) * cos(z(k) - z_center)
            exponent = -dist_sq / (width**2)
            
            if (exponent < -40.0d0) then
              gauss_prof = 0.0d0
            else
              gauss_prof = exp(exponent)
            end if
            
            p(eq_de, i, j, k) = 0.1d0 + 1.0d0 * gauss_prof
            p(eq_pr, i, j, k) = 0.1d0 + 1.0d0 * gauss_prof
          end do
        end do
      end do

    case('FishMoncEqu', '10') ! Toro de Fishbone-Moncrief en corte ecuatorial (2D r-phi)
      print *, "-> Inicializando Toro de Fishbone-Moncrief (corte ecuatorial 2D)..."
      r_in       = 6.25d0 * bh_mass      
      r_max_dens = 9.25d0 * bh_mass  
      K_poly     = 0.0015d0            
      
      call setup_fishbone_moncrief_initial(p(:, 1:nx, 1:ny, 1:nz), r_in, r_max_dens, K_poly)

    case('FishMoncSag', '11') ! Toro de Fishbone-Moncrief en corte sagital (2D r-theta)
      print *, "-> Inicializando Toro de Fishbone-Moncrief (corte sagital 2D)..."
      r_in       = 6.25d0 * bh_mass      
      r_max_dens = 9.25d0 * bh_mass  
      K_poly     = 0.0015d0
      
      call setup_fishbone_moncrief_initial(p(:, 1:nx, 1:ny, 1:nz), r_in, r_max_dens, K_poly)

    case('FishMonc3D', '12')
      print *, "-> Inicializando Toro de Fishbone-Moncrief completo (3D)..."
      r_in       = 6.0d0 * bh_mass      
      r_max_dens = 12.0d0 * bh_mass  
      K_poly     = 0.015d0
      
      call setup_fishbone_moncrief_initial(p(:, 1:nx, 1:ny, 1:nz), r_in, r_max_dens, K_poly)

    case default
      print *, "CRITICAL ERROR: case_name no valido en initial_conditions."
      stop
    end select
  end subroutine set_initial_conditions


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

    local_nx = size(prim_state, 2)
    local_ny = size(prim_state, 3)
    local_nz = size(prim_state, 4)

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
          prim_state(eq_de, i, j, k) = rho_final
          prim_state(eq_pr, i, j, k) = pres
          prim_state(eq_vx, i, j, k) = ux / W_lorentz + beta(1) / alpha
          prim_state(eq_vy, i, j, k) = beta(2) / alpha
          prim_state(eq_vz, i, j, k) = beta(3) / alpha
          
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
    integer :: i, j, k, ii, jj
    real*8  :: l_ang, W_in, W_r, enthalpy, ut_sq, r_loc
    real*8  :: alpha, beta(3), g(3,3)
    real*8  :: g_tt, g_tphi, g_phiphi, Omega, coord_in

    local_nx = size(prim_state, 2)
    local_ny = size(prim_state, 3)
    local_nz = size(prim_state, 4)
    
    ! Momento angular específico (constante de integración)
    l_ang = sqrt(bh_mass * r_max_dens) / (1.0d0 - 2.0d0 * bh_mass / r_max_dens)

    ! CÁLCULO DE W_in (Usando el puntero métrico general)
    if (use_log_r) then
      coord_in = log(r_in)
    else
      coord_in = r_in
    end if
    
    call calculate_metric(coord_in, pi/2.0d0, alpha=alpha, beta=beta, gamma=g)
    
    ! Construcción Universal de g_tt y g_tphi vía ADM
    g_tt = -alpha**2
    do ii = 1, 3
      do jj = 1, 3
        g_tt = g_tt + g(ii,jj) * beta(ii) * beta(jj)
      end do
    end do
    
    g_tphi = 0.0d0
    do jj = 1, 3
      g_tphi = g_tphi + g(3,jj) * beta(jj)
    end do
    g_phiphi = g(3,3)

    ut_sq = ((g_tphi)**2 - g_tt * g_phiphi) / (g_phiphi + 2.0d0 * l_ang * g_tphi + l_ang**2 * g_tt)
    W_in  = 0.5d0 * log(ut_sq)

    !!$OMP PARALLEL DO PRIVATE(i, j, k, r_loc, alpha, beta, g, ii, jj, g_tt, g_tphi, g_phiphi, ut_sq, Omega, W_r, enthalpy)
    do k = 1, local_nz
      do j = 1, local_ny
        do i = 1, local_nx
          
          if (use_log_r) then
            r_loc = exp(x(i))
          else
            r_loc = x(i)
          end if
          
          alpha   = alpha_c(i,j,k)
          beta(:) = beta_c(:,i,j,k)
          g(:,:)  = gamma_c(:,:,i,j,k)

          ! Blindaje de frontera interna
          if (r_loc < r_in) then
            prim_state(eq_de, i, j, k) = rho_floor
            prim_state(eq_pr, i, j, k) = p_floor
            prim_state(eq_vx:eq_vz, i, j, k) = 0.0d0
            cycle
          end if

          ! Construcción Local Universal ADM
          g_tt = -alpha**2
          do ii = 1, 3
            do jj = 1, 3
              g_tt = g_tt + g(ii,jj) * beta(ii) * beta(jj)
            end do
          end do
          
          g_tphi = 0.0d0
          do jj = 1, 3
            g_tphi = g_tphi + g(3,jj) * beta(jj)
          end do
          g_phiphi = g(3,3)

          ! Solución Analítica General
          ut_sq = ((g_tphi)**2 - g_tt * g_phiphi) / (g_phiphi + 2.0d0 * l_ang * g_tphi + l_ang**2 * g_tt)
          
          if (ut_sq > 0.0d0) then
            W_r = 0.5d0 * log(ut_sq)
            enthalpy = exp(W_in - W_r)

            if (enthalpy > 1.0d0) then
              ! Velocidad Angular de la órbita
              Omega = -(g_tphi + l_ang * g_tt) / (g_phiphi + l_ang * g_tphi)
              
              prim_state(eq_de, i, j, k) = ((enthalpy - 1.0d0) * (adb_idx - 1.0d0) / (K_poly * adb_idx))**(1.0d0 / (adb_idx - 1.0d0))
              prim_state(eq_pr, i, j, k) = K_poly * prim_state(eq_de, i, j, k)**adb_idx

              ! Velocidades Eulerianas (v^i = u^i / (alpha u^t) + beta^i / alpha)
              ! Como u^r = 0 y u^theta = 0, el cálculo es directo y bellísimo
              prim_state(eq_vx, i, j, k) = beta(1) / alpha 
              prim_state(eq_vy, i, j, k) = beta(2) / alpha 
              prim_state(eq_vz, i, j, k) = Omega / alpha + beta(3) / alpha 
            else
              ! Vacío exterior (Atmósfera)
              prim_state(eq_de, i, j, k) = rho_floor
              prim_state(eq_pr, i, j, k) = p_floor
              prim_state(eq_vx:eq_vz, i, j, k) = 0.0d0
            end if
          else
            ! Zonas causalmente desconectadas
            prim_state(eq_de, i, j, k) = rho_floor
            prim_state(eq_pr, i, j, k) = p_floor
            prim_state(eq_vx:eq_vz, i, j, k) = 0.0d0
          end if
        end do
      end do
    end do
    !!$OMP END PARALLEL DO
    
    print *, ">>> Toro FM Universal ADM configurado. Angular Mom. =", l_ang, "W_in =", W_in
  end subroutine setup_fishbone_moncrief_initial

  ! Subrutina: inject_pressure_noise
  ! Rompe la simetría perfecta para detonar inestabilidades multidimensionales.
  subroutine inject_pressure_noise(prim_state)
    implicit none
    ! Cambiamos a la dimensionalidad 4D estándar
    real*8, intent(inout) :: prim_state(:,:,:,:)
    integer :: i, j, k
    integer :: local_nx, local_ny, local_nz
    real*8  :: rand_val, rho_max_local, factor

    local_nx = size(prim_state, 2)
    local_ny = size(prim_state, 3)
    local_nz = size(prim_state, 4)

    print *, ">>> ATENCION: Inyectando 4% de ruido aleatorio en la presion..."

    rho_max_local = maxval(prim_state(eq_de, :, :, :))
    call random_seed()

    do k = 1, local_nz
      do j = 1, local_ny
        do i = 1, local_nx
          ! Solo perturbar el cuerpo del disco
          if (prim_state(eq_de, i, j, k) > 0.05d0 * rho_max_local) then
            call random_number(rand_val) 
            rand_val = 2.0d0 * rand_val - 1.0d0 ! Entre -1 y 1
            factor = 1.0d0 + 0.04d0 * rand_val
            prim_state(eq_pr, i, j, k) = prim_state(eq_pr, i, j, k) * factor
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
  !   rho_max_local = maxval(prim_state(eq_de, :, :, :))

  !   do k = 1, local_nz
  !     ! PRECAUCIÓN: Asegúrate de que z(k) es tu coordenada azimutal phi
  !     phi_loc = z(k) 
      
  !     do j = 1, local_ny
  !       do i = 1, local_nx

  !         ! Blindaje del vacío: Solo perturbamos el corazón del toroide
  !         if (prim_state(eq_de, i, j, k) > 0.05d0 * rho_max_local) then
            
  !           envelope = cos(mode_m * phi_loc)
            
  !           ! Factor de perturbación del 1%
  !           factor = 1.0d0 + 0.01d0 * envelope
            
  !           ! 1. Perturbar la densidad
  !           prim_state(eq_de, i, j, k) = prim_state(eq_de, i, j, k) * factor
            
  !           ! 2. Recalcular la presión isentrópica (CORRECCIÓN CRÍTICA)
  !           ! En lugar de hardcodear K = 0.015, usamos la relación P_nueva = P_vieja * (Rho_nueva / Rho_vieja)^Gamma
  !           ! Esto hace que la subrutina sea universal y a prueba de errores.
  !           prim_state(eq_pr, i, j, k) = prim_state(eq_pr, i, j, k) * (factor**adb_idx)
            
  !         end if                    
  !       end do
  !     end do
  !   end do
  ! end subroutine inject_density_mode

  subroutine inject_density_mode(prim_state, mode_m)
    implicit none
    real*8, intent(inout) :: prim_state(:,:,:,:)
    real*8, intent(in)    :: mode_m
    integer :: i, j, k
    integer :: local_nx, local_ny, local_nz
    real*8  :: phi_loc, rho_max_local, envelope, factor, rand_val

    local_nx = size(prim_state, 2)
    local_ny = size(prim_state, 3)
    local_nz = size(prim_state, 4)

    print *, ">>> ATENCION: Inyectando modo m =", mode_m, " + Ruido Blanco Isentrópico..."
    
    rho_max_local = maxval(prim_state(eq_de, :, :, :))

    ! Inicializar la semilla aleatoria para que el ruido cambie
    call random_seed()

    do k = 1, local_nz
      ! Nota: Asegúrate de que z(k) o y(j) corresponda a tu coordenada azimutal phi
      phi_loc = z(k) 
      
      do j = 1, local_ny
        do i = 1, local_nx

          ! Blindaje del vacío: Solo perturbamos el corazón del toroide (> 5% de la densidad máxima)
          if (prim_state(eq_de, i, j, k) > 0.05d0 * rho_max_local) then
            
            ! 1. Generar ruido aleatorio (rand_val estará entre 0.0 y 1.0)
            call random_number(rand_val)
            
            ! Escalar el ruido para que oscile entre -1.0 y 1.0
            rand_val = 2.0d0 * rand_val - 1.0d0
            
            ! 2. Súper-posición: 1% de semilla coherente (m=2, m=3) + 4% de ruido blanco aleatorio
            envelope = 0.01d0 * cos(mode_m * phi_loc) + 0.04d0 * rand_val
            
            factor = 1.0d0 + envelope
            
            ! 3. Aplicar a la densidad
            prim_state(eq_de, i, j, k) = prim_state(eq_de, i, j, k) * factor
            
            ! 4. Recalcular presión isentrópica rigurosamente
            prim_state(eq_pr, i, j, k) = prim_state(eq_pr, i, j, k) * (factor**adb_idx)
            
          end if                    
        end do
      end do
    end do
  end subroutine inject_density_mode

  subroutine inject_density_noise(prim_state)
    implicit none
    real*8, intent(inout) :: prim_state(:,:,:,:)
    integer :: i, j, k
    integer :: local_nx, local_ny, local_nz
    real*8  :: rho_max_local, factor, rand_val

    local_nx = size(prim_state, 2)
    local_ny = size(prim_state, 3)
    local_nz = size(prim_state, 4)

    print *, ">>> ATENCION: Inyectando 4% de Ruido Blanco Isentrópico en T=0..."
    
    rho_max_local = maxval(prim_state(eq_de, :, :, :))
    call random_seed()

    do k = 1, local_nz
      do j = 1, local_ny
        do i = 1, local_nx
          ! Solo perturbar el cuerpo del disco
          if (prim_state(eq_de, i, j, k) > 0.05d0 * rho_max_local) then
            
            call random_number(rand_val)
            rand_val = 2.0d0 * rand_val - 1.0d0 ! Entre -1 y 1
            
            factor = 1.0d0 + 0.04d0 * rand_val
            
            prim_state(eq_de, i, j, k) = prim_state(eq_de, i, j, k) * factor
            prim_state(eq_pr, i, j, k) = prim_state(eq_pr, i, j, k) * (factor**adb_idx)
            
          end if                    
        end do
      end do
    end do
  end subroutine inject_density_noise


  ! Subrutina: set_boundary_conditions
  ! Propósito: Rellena las "celdas fantasma" preservando la topología del dominio en 3D.
  subroutine set_boundary_conditions(q_state)
    implicit none
    real*8, intent(inout) :: q_state(neq, -nghost:nx+nghost, -nghost:ny+nghost, -nghost:nz+nghost)
    integer :: i, j, k, g, k_shift

    select case(trim(case_name))

    case('Sod', '1', 'Strong', '2', 'Shu', '3') ! Cajas 1D - Gradiente Cero (Outflow libre)
      do g = 1, nghost
        q_state(:, 1-g, :, :) = q_state(:, 1, :, :)
        q_state(:, nx+g, :, :) = q_state(:, nx, :, :)
      end do

      do g = 1, nghost
        q_state(:, :, 1-g, :) = q_state(:, :, ny + 1 - g, :)
        q_state(:, :, ny+g, :) = q_state(:, :, g, :)
        
        q_state(:, :, :, 1-g) = q_state(:, :, :, nz + 1 - g)
        q_state(:, :, :, nz+g) = q_state(:, :, :, g)
      end do

    case('KHI', '4', 'Conv', '6') ! Fronteras periódicas
      do g = 1, nghost
        q_state(:, 1-g, :, :) = q_state(:, nx + 1 - g, :, :)
        q_state(:, nx+g, :, :) = q_state(:, g, :, :)
        
        q_state(:, :, 1-g, :) = q_state(:, :, ny + 1 - g, :)
        q_state(:, :, ny+g, :) = q_state(:, :, g, :)
        
        q_state(:, :, :, 1-g) = q_state(:, :, :, nz + 1 - g)
        q_state(:, :, :, nz+g) = q_state(:, :, :, g)
      end do

    case('Jet', '5') ! Slab Jet Cartesiano
      
      ! Reseteamos fantasmas a ambiente
      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = -nghost, nz+nghost
        do j = -nghost, ny+nghost
          do i = -nghost, nx+nghost
            if (i < 1 .or. i > nx .or. j < 1 .or. j > ny .or. k < 1 .or. k > nz) then
                q_state(eq_de, i, j, k) = 10.0d0
                q_state(eq_pr, i, j, k) = 0.01d0
                q_state(eq_vx:eq_vz, i, j, k) = 0.0d0
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
            q_state(:, 1-g, j, k) = q_state(:, g, j, k)
            q_state(eq_vx, 1-g, j, k) = -q_state(eq_vx, g, j, k) 
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! 2. Z -> TOBERA + OUTFLOW
      !$OMP PARALLEL DO PRIVATE(i, j, g)
      do j = -nghost, ny+nghost
        do i = -nghost, nx+nghost
          do g = 1, nghost
            q_state(:, i, j, nz+g) = q_state(:, i, j, nz)
            q_state(eq_vz, i, j, nz+g) = max(0.0d0, q_state(eq_vz, i, j, nz))

            if (x(i) <= 1.0d0) then
              q_state(eq_de, i, j, 1-g) = 0.1d0
              q_state(eq_pr, i, j, 1-g) = 0.01d0
              q_state(eq_vx:eq_vy, i, j, 1-g) = 0.0d0
              q_state(eq_vz, i, j, 1-g) = 0.99d0
            else
              q_state(:, i, j, 1-g) = q_state(:, i, j, g)
              q_state(eq_vz, i, j, 1-g) = -q_state(eq_vz, i, j, g) 
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
            q_state(:, 1-g, j, k) = q_state(:, 1, j, k)
            q_state(:, nx+g, j, k) = michel_injector(:, g, j, k)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! 2. Polar (Y) y Azimutal (Z)
      do g = 1, nghost
        q_state(:, :, 1-g, :) = q_state(:, :, ny + 1 - g, :)
        q_state(:, :, ny+g, :) = q_state(:, :, g, :)
        q_state(:, :, :, 1-g) = q_state(:, :, :, nz + 1 - g)
        q_state(:, :, :, nz+g) = q_state(:, :, :, g)
      end do

    case('OffAxis', '9', 'FishMoncEqu', '10') 
      
      ! 1. Radial
      !$OMP PARALLEL DO PRIVATE(j, k, g)
      do k = 1, nz
        do j = 1, ny
          do g = 1, nghost
            q_state(:, 1-g, j, k) = q_state(:, 1, j, k) 
            q_state(eq_de, 1-g, j, k) = max(q_state(eq_de, 1, j, k), rho_floor)
            q_state(eq_pr, 1-g, j, k) = max(q_state(eq_pr, 1, j, k), p_floor)
            q_state(eq_vx, 1-g, j, k) = min(0.0d0, q_state(eq_vx, 1-g, j, k))
            if (q_state(eq_de, 1, j, k) < 1.0d-6) q_state(eq_vy:eq_vz, 1-g, j, k) = 0.0d0

            q_state(eq_de, nx+g, j, k) = max(q_state(eq_de, nx, j, k), rho_floor)
            q_state(eq_pr, nx+g, j, k) = max(q_state(eq_pr, nx, j, k), p_floor)
            q_state(eq_vx, nx+g, j, k) = max(0.0d0, q_state(eq_vx, nx, j, k)) 
            q_state(eq_vy, nx+g, j, k) = q_state(eq_vy, nx, j, k)
            q_state(eq_vz, nx+g, j, k) = q_state(eq_vz, nx, j, k)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! 2. Transversales
      do g = 1, nghost
        q_state(:, :, 1-g, :) = q_state(:, :, ny + 1 - g, :)
        q_state(:, :, ny+g, :) = q_state(:, :, g, :)
        q_state(:, :, :, 1-g) = q_state(:, :, :, nz + 1 - g)
        q_state(:, :, :, nz+g) = q_state(:, :, :, g)
      end do

    case('FishMoncSag', '11')
      
      ! 1. Llenado Radial (X) - Físico
      !$OMP PARALLEL DO PRIVATE(j, k, g)
      do k = 1, nz
        do j = 1, ny
          do g = 1, nghost
            ! --- FRONTERA INTERNA (Agujero Negro) ---
            q_state(:, 1-g, j, k) = q_state(:, 1, j, k) 
            q_state(eq_de, 1-g, j, k) = max(q_state(eq_de, 1, j, k), rho_floor)
            q_state(eq_pr, 1-g, j, k) = max(q_state(eq_pr, 1, j, k), p_floor)
            ! CORRECCIÓN VITAL: El gas solo puede CAER hacia adentro (vx <= 0)
            q_state(eq_vx, 1-g, j, k) = min(0.0d0, q_state(eq_vx, 1, j, k))

            ! --- FRONTERA EXTERNA ---
            q_state(:, nx+g, j, k) = q_state(:, nx, j, k)
            q_state(eq_de, nx+g, j, k) = max(q_state(eq_de, nx, j, k), rho_floor)
            q_state(eq_pr, nx+g, j, k) = max(q_state(eq_pr, nx, j, k), p_floor)
            ! El gas solo puede SALIR del dominio (vx >= 0)
            q_state(eq_vx, nx+g, j, k) = max(0.0d0, q_state(eq_vx, nx, j, k)) 
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! 2. Llenado Polar (Y) - Sobre X Expandido
      !$OMP PARALLEL DO PRIVATE(i, k, g)
      do k = 1, nz
        do i = -nghost, nx+nghost
          do g = 1, nghost
            q_state(:, i, 1-g, k) = q_state(:, i, g, k)
            q_state(eq_vy, i, 1-g, k) = -q_state(eq_vy, i, g, k) 
            q_state(:, i, ny+g, k) = q_state(:, i, ny+1-g, k)
            q_state(eq_vy, i, ny+g, k) = -q_state(eq_vy, i, ny+1-g, k) 
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! 3. Llenado Azimutal (Z) - Global
      do g = 1, nghost
        q_state(:, :, :, 1-g) = q_state(:, :, :, nz + 1 - g)
        q_state(:, :, :, nz+g) = q_state(:, :, :, g)
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
            q_state(:, 1-g, j, k) = q_state(:, 1, j, k) 
            q_state(eq_de, 1-g, j, k) = max(q_state(eq_de, 1, j, k), rho_floor)
            q_state(eq_pr, 1-g, j, k) = max(q_state(eq_pr, 1, j, k), p_floor)
            
            ! Externa (Outflow)
            q_state(:, nx+g, j, k) = q_state(:, nx, j, k)
            q_state(eq_de, nx+g, j, k) = max(q_state(eq_de, nx, j, k), rho_floor)
            q_state(eq_pr, nx+g, j, k) = max(q_state(eq_pr, nx, j, k), p_floor)
            q_state(eq_vx, nx+g, j, k) = max(0.0d0, q_state(eq_vx, nx, j, k)) 
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
            q_state(:, i, 1-g, k) = q_state(:, i, g, k_shift)
            q_state(eq_vy, i, 1-g, k) = -q_state(eq_vy, i, g, k_shift) 
            q_state(eq_vz, i, 1-g, k) =  q_state(eq_vz, i, g, k_shift)

            ! Polo Sur
            q_state(:, i, ny+g, k) = q_state(:, i, ny+1-g, k_shift)
            q_state(eq_vy, i, ny+g, k) = -q_state(eq_vy, i, ny+1-g, k_shift) 
            q_state(eq_vz, i, ny+g, k) =  q_state(eq_vz, i, ny+1-g, k_shift)
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! ========================================================================
      ! 3. Llenado Azimutal (Z): Periódico - Sobre matriz completa
      ! ========================================================================
      do g = 1, nghost
        q_state(:, :, :, 1-g) = q_state(:, :, :, nz + 1 - g)
        q_state(:, :, :, nz+g)  = q_state(:, :, :, g)
      end do

    end select
  end subroutine set_boundary_conditions

end module conditions