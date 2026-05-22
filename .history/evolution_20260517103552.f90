! =======================================================================================
! Módulo: evolution (Integración Espacial y Ensamblaje del RHS)
! ---------------------------------------------------------------------------------------
! Propósito: Actúa como el puente entre la geometría de la malla y la física local. 
! Realiza el barrido dimensional (Dimensional Splitting) a lo largo de los ejes X y Z.
! En cada eje: 
!   1. Reconstruye los estados en las interfaces.
!   2. Resuelve el problema de Riemann (HLLE) para obtener los flujos.
!   3. Calcula la divergencia de los flujos y suma los términos fuente gravitacionales.
! =======================================================================================
module evolution
  use variables
  use metrics
  use reconstruction
  use fluxes
  use equations
  implicit none
  private
  public :: calc_rhs, update_adaptive_dt

contains

  ! ===================================================================================
  ! ENSAMBLADOR DEL LADO DERECHO (RHS) Y ESTABILIDAD NUMÉRICA
  ! ===================================================================================

  ! Subrutina: calc_rhs
  ! Ensambla el lado derecho (RHS) de las ecuaciones de Navier-Stokes Relativistas:
  subroutine calc_rhs(p_in, rhs_out)
    use variables
    implicit none
    real*8, intent(in)  :: p_in(neq, -nghost:nx+nghost, -nghost:ny+nghost, -nghost:nz+nghost)
    real*8, intent(out) :: rhs_out(neq, 1:nx, 1:ny, 1:nz) 
    
    integer :: i, j, k
    real*8, allocatable :: prim_1d(:,:), q_L_strip(:,:), q_R_strip(:,:), flux_strip(:,:)
    real*8 :: src_term(neq)
    
    ! Variables temporales para pasar la métrica a las interfaces
    real*8 :: a, b(3), g(3,3), sqg

    rhs_out = 0.0d0

    ! =====================================================================
    ! 1. BARRIDO EN DIRECCIÓN X (Radial / Horizontal)
    ! =====================================================================
    !$OMP PARALLEL PRIVATE(j, k, i, prim_1d, q_L_strip, q_R_strip, flux_strip, src_term, a, b, g, sqg)
    
    ! Alojamiento dinámico de tiras 1D de memoria por cada hilo OMP
    allocate(prim_1d(neq, -nghost:nx+nghost))
    allocate(q_L_strip(neq, 0:nx))
    allocate(q_R_strip(neq, 0:nx))
    allocate(flux_strip(neq, 0:nx))
    
    !$OMP DO COLLAPSE(2)
    do k = 1, nz 
      do j = 1, ny
        ! A. Extraer tira unidimensional a lo largo de r
        prim_1d(:, -nghost:nx+nghost) = p_in(:, -nghost:nx+nghost, j, k)
        
        ! B. Reconstrucción Espacial (Devuelve estados en interfaces 0 a nx)
        call reconstruct_1d_core(prim_1d, nx, q_L_strip, q_R_strip, DIR_X)
        
        ! C. Cálculo de Flujos en las Interfaces (0 a nx)
        do i = 0, nx 
          ! Leer de la caché de la interfaz X (Mapeo: índice 0 físico -> índice 1 del arreglo de interfaz)
          a    = alpha_f_x(i+1, j, k)
          b(:) = beta_f_x(:, i+1, j, k)
          g(:,:) = gamma_f_x(:, :, i+1, j, k)
          sqg  = sqrt_gamma_f_x(i+1, j, k)
          
          call resolve_riemann_problem(q_L_strip(:,i), q_R_strip(:,i), DIR_X, a, b, g, sqg, flux_strip(:,i))
        end do
        
        ! D. Ensamblaje del RHS en los Centros Celulares (1 a nx)
        do i = 1, nx
          call calc_sources(p_in(:,i,j,k), src_term, i, j, k)
          rhs_out(:,i,j,k) = -(flux_strip(:,i) - flux_strip(:,i-1)) / dx + src_term(:)
        end do
      end do
    end do
    !$OMP END DO

    deallocate(prim_1d, q_L_strip, q_R_strip, flux_strip)
    !$OMP END PARALLEL


    ! =====================================================================
    ! 2. BARRIDO EN DIRECCIÓN Y (Polar / Sagital)
    ! =====================================================================
    if (ny > 1) then
      !$OMP PARALLEL PRIVATE(i, j, k, prim_1d, q_L_strip, q_R_strip, flux_strip, a, b, g, sqg)
      
      allocate(prim_1d(neq, -nghost:ny+nghost))
      allocate(q_L_strip(neq, 0:ny))
      allocate(q_R_strip(neq, 0:ny))
      allocate(flux_strip(neq, 0:ny))
      
      !$OMP DO COLLAPSE(2)
      do k = 1, nz 
        do i = 1, nx 
          ! A. Extraer tira unidimensional a lo largo de theta
          prim_1d(:, -nghost:ny+nghost) = p_in(:, i, -nghost:ny+nghost, k)
          
          ! B. Reconstrucción
          call reconstruct_1d_core(prim_1d, ny, q_L_strip, q_R_strip, DIR_Y)

          ! C. Cálculo de Flujos (las interfaces caen en los bordes de celda en theta)
          do j = 0, ny
            a    = alpha_f_y(i, j+1, k)
            b(:) = beta_f_y(:, i, j+1, k)
            g(:,:) = gamma_f_y(:, :, i, j+1, k)
            sqg  = sqrt_gamma_f_y(i, j+1, k)
            
            call resolve_riemann_problem(q_L_strip(:,j), q_R_strip(:,j), DIR_Y, a, b, g, sqg, flux_strip(:,j))
          end do

          ! D. Acumulación de Flujos Y en el RHS
          do j = 1, ny 
            rhs_out(:,i,j,k) = rhs_out(:,i,j,k) - (flux_strip(:,j) - flux_strip(:,j-1)) / dy
          end do
        end do
      end do
      !$OMP END DO
      
      deallocate(prim_1d, q_L_strip, q_R_strip, flux_strip)
      !$OMP END PARALLEL
    end if


    ! =====================================================================
    ! 3. BARRIDO EN DIRECCIÓN Z (Azimutal / Vertical)
    ! =====================================================================
    if (nz > 1) then
      !$OMP PARALLEL PRIVATE(i, j, k, prim_1d, q_L_strip, q_R_strip, flux_strip, a, b, g, sqg)
      
      allocate(prim_1d(neq, -nghost:nz+nghost))
      allocate(q_L_strip(neq, 0:nz))
      allocate(q_R_strip(neq, 0:nz))
      allocate(flux_strip(neq, 0:nz))
      
      !$OMP DO COLLAPSE(2)
      do j = 1, ny
        do i = 1, nx 
          ! A. Extraer tira unidimensional a lo largo de phi
          prim_1d(:, -nghost:nz+nghost) = p_in(:, i, j, -nghost:nz+nghost)
          
          ! B. Reconstrucción
          call reconstruct_1d_core(prim_1d, nz, q_L_strip, q_R_strip, DIR_Z)

          ! C. Cálculo de Flujos
          do k = 0, nz
            a    = alpha_f_z(i, j, k+1)
            b(:) = beta_f_z(:, i, j, k+1)
            g(:,:) = gamma_f_z(:, :, i, j, k+1)
            sqg  = sqrt_gamma_f_z(i, j, k+1)
            
            call resolve_riemann_problem(q_L_strip(:,k), q_R_strip(:,k), DIR_Z, a, b, g, sqg, flux_strip(:,k))
          end do

          ! D. Acumulación de Flujos Z en el RHS
          do k = 1, nz 
            rhs_out(:,i,j,k) = rhs_out(:,i,j,k) - (flux_strip(:,k) - flux_strip(:,k-1)) / dz
          end do
        end do
      end do
      !$OMP END DO
      
      deallocate(prim_1d, q_L_strip, q_R_strip, flux_strip)
      !$OMP END PARALLEL
    end if

  end subroutine calc_rhs

  ! ====================================================================
  ! Cálculo del Paso de Tiempo Dinámico (CFL Adaptativo)
  ! Evalúa la velocidad de la luz coordenada y el transporte del fluido
  ! para garantizar la estabilidad estricta del integrador RK3.
  ! ====================================================================
  subroutine update_adaptive_dt()
    use variables
    implicit none
    
    real*8 :: max_vx, max_vy, max_vz
    real*8 :: local_vx, local_vy, local_vz, alpha, beta(3), ginv(3,3)
    integer :: i, j, k
    
    max_vx = 1.0d-10
    max_vy = 1.0d-10
    max_vz = 1.0d-10

    ! Buscamos las velocidades máximas
    !$OMP PARALLEL DO PRIVATE(i, j, k, alpha, beta, ginv, local_vx, local_vy, local_vz) &
    !$OMP REDUCTION(max:max_vx, max_vy, max_vz) COLLAPSE(3)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          
          ! Leer métrica INVERSA (contravariante) en el centro de la celda
          alpha    = alpha_c(i,j,k)
          beta(:)  = beta_c(:,i,j,k)
          ginv(:,:) = gamma_inv_c(:,:,i,j,k) ! gamma^ij
          
          ! Límite Causal Absoluto Físico (Cono de Luz Coordenado: ds^2 = 0)
          ! La velocidad coordenada de la luz en Kerr usa el componente gamma^ii
          local_vx = alpha * sqrt(ginv(1,1)) + abs(beta(1))
          local_vy = alpha * sqrt(ginv(2,2)) + abs(beta(2))
          local_vz = alpha * sqrt(ginv(3,3)) + abs(beta(3))
          
          max_vx = max(max_vx, local_vx)
          max_vy = max(max_vy, local_vy)
          max_vz = max(max_vz, local_vz)
        end do
      end do
    end do
    !$OMP END PARALLEL DO

    ! Aplicamos la condición Courant-Friedrichs-Lewy (CFL) evaluando las 3 dimensiones
    if (ny > 1 .and. nz > 1) then
      dt = CFL * min(min(dx / max_vx, dy / max_vy), dz / max_vz)
    else if (ny > 1) then
      dt = CFL * min(dx / max_vx, dy / max_vy)
    else if (nz > 1) then
      dt = CFL * min(dx / max_vx, dz / max_vz)
    else
      dt = CFL * (dx / max_vx)
    end if
    
  end subroutine update_adaptive_dt

end module evolution