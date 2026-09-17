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
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: calc_rhs, update_adaptive_dt, initialize_rhs_accelerator, finalize_rhs_accelerator
  public :: calc_rhs_cpu, calc_rhs_accelerator

  real*8, allocatable :: qlx(:,:,:,:), qrx(:,:,:,:), fluxx(:,:,:,:)
  real*8, allocatable :: qly(:,:,:,:), qry(:,:,:,:), fluxy(:,:,:,:)
  real*8, allocatable :: qlz(:,:,:,:), qrz(:,:,:,:), fluxz(:,:,:,:)

contains

  subroutine calc_rhs(p_in, rhs_out)
    real*8, intent(in)  :: p_in(-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost,neq)
    real*8, intent(out) :: rhs_out(1:nx,1:ny,1:nz,neq)

    if (accelerator_enabled) then
      call calc_rhs_accelerator(p_in,rhs_out)
    else
      call calc_rhs_cpu(p_in,rhs_out)
    end if
  end subroutine calc_rhs

  ! ===================================================================================
  ! ENSAMBLADOR DEL LADO DERECHO (RHS) Y ESTABILIDAD NUMÉRICA
  ! ===================================================================================

  ! Subrutina: calc_rhs
  ! Ensambla el lado derecho (RHS) de las ecuaciones de Navier-Stokes Relativistas:
  subroutine calc_rhs_cpu(p_in, rhs_out)
    use variables
    implicit none
    real*8, intent(in)  :: p_in(-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost,neq)
    real*8, intent(out) :: rhs_out(1:nx,1:ny,1:nz,neq)

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
        do i=-nghost,nx+nghost
          prim_1d(:,i)=p_in(i,j,k,:)
        end do

        ! B. Reconstrucción Espacial (Devuelve estados en interfaces 0 a nx)
        call reconstruct_1d_core(prim_1d, nx, q_L_strip, q_R_strip, DIR_X)

        ! C. Cálculo de Flujos en las Interfaces (0 a nx)
        do i = 0, nx
          ! La cara rho=0 del sistema cilíndrico tiene área nula y una métrica
          ! coordenada degenerada. El flujo densitizado se fija exactamente a cero.
          if (trim(geom_type) == 'Cylindrical' .and. &
              abs(x_face(i)) <= 64.0d0*epsilon(1.0d0)) then
            flux_strip(:,i) = 0.0d0
          else
            ! Leer de la caché de la interfaz X (índice de cara i -> caché i+1).
            a    = alpha_f_x(i+1, j, k)
            b(:) = beta_f_x(:, i+1, j, k)
            g(:,:) = gamma_f_x(:, :, i+1, j, k)
            sqg  = sqrt_gamma_f_x(i+1, j, k)

            call resolve_riemann_problem(q_L_strip(:,i), q_R_strip(:,i), &
                                         DIR_X, a, b, g, sqg, flux_strip(:,i))
          end if
        end do

        ! D. Ensamblaje del RHS en los Centros Celulares (1 a nx)
        do i = 1, nx
          call calc_sources(p_in(i,j,k,:), src_term, i, j, k)
          rhs_out(i,j,k,:) = -(flux_strip(:,i) - flux_strip(:,i-1)) / dx + src_term(:)
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
          do j=-nghost,ny+nghost
            prim_1d(:,j)=p_in(i,j,k,:)
          end do

          ! B. Reconstrucción
          call reconstruct_1d_core(prim_1d, ny, q_L_strip, q_R_strip, DIR_Y)

          ! C. Cálculo de Flujos (las interfaces caen en los bordes de celda en theta).
          ! En coordenadas polares las caras theta=0,pi tienen área física nula:
          ! el flujo densitizado sqrt(gamma) F^theta es exactamente cero. No se debe
          ! llamar al solver de Riemann allí porque gamma_ij es degenerada en el eje.
          do j = 0, ny
            if (is_polar_geometry(geom_type) .and. &
                abs(sin(y_face(j))) <= 64.0d0*epsilon(1.0d0)) then
              flux_strip(:,j) = 0.0d0
            else
              a      = alpha_f_y(i, j+1, k)
              b(:)   = beta_f_y(:, i, j+1, k)
              g(:,:) = gamma_f_y(:, :, i, j+1, k)
              sqg    = sqrt_gamma_f_y(i, j+1, k)

              call resolve_riemann_problem(q_L_strip(:,j), q_R_strip(:,j), &
                                           DIR_Y, a, b, g, sqg, flux_strip(:,j))
            end if
          end do

          ! D. Acumulación de Flujos Y en el RHS
          do j = 1, ny
            rhs_out(i,j,k,:) = rhs_out(i,j,k,:) - (flux_strip(:,j) - flux_strip(:,j-1)) / dy
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
          do k=-nghost,nz+nghost
            prim_1d(:,k)=p_in(i,j,k,:)
          end do

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
            rhs_out(i,j,k,:) = rhs_out(i,j,k,:) - (flux_strip(:,k) - flux_strip(:,k-1)) / dz
          end do
        end do
      end do
      !$OMP END DO

      deallocate(prim_1d, q_L_strip, q_R_strip, flux_strip)
      !$OMP END PARALLEL
    end if

  end subroutine calc_rhs_cpu

  subroutine initialize_rhs_accelerator(force_for_test)
    logical, intent(in), optional :: force_for_test
    logical :: force
    force=.false.
    if (present(force_for_test)) force=force_for_test
    if (.not. accelerator_enabled .and. .not. force) return

    allocate(qlx(0:nx,1:ny,1:nz,neq),qrx(0:nx,1:ny,1:nz,neq), &
             fluxx(0:nx,1:ny,1:nz,neq))
    !$acc enter data create(qlx,qrx,fluxx)
    if (ny > 1) then
      allocate(qly(1:nx,0:ny,1:nz,neq),qry(1:nx,0:ny,1:nz,neq), &
               fluxy(1:nx,0:ny,1:nz,neq))
      !$acc enter data create(qly,qry,fluxy)
    end if
    if (nz > 1) then
      allocate(qlz(1:nx,1:ny,0:nz,neq),qrz(1:nx,1:ny,0:nz,neq), &
               fluxz(1:nx,1:ny,0:nz,neq))
      !$acc enter data create(qlz,qrz,fluxz)
    end if
  end subroutine initialize_rhs_accelerator

  subroutine finalize_rhs_accelerator(force_for_test)
    logical, intent(in), optional :: force_for_test
    logical :: force
    force=.false.
    if (present(force_for_test)) force=force_for_test
    if (.not. accelerator_enabled .and. .not. force) return
    if (allocated(qlx)) then
      !$acc exit data delete(qlx,qrx,fluxx)
      deallocate(qlx,qrx,fluxx)
    end if
    if (allocated(qly)) then
      !$acc exit data delete(qly,qry,fluxy)
      deallocate(qly,qry,fluxy)
    end if
    if (allocated(qlz)) then
      !$acc exit data delete(qlz,qrz,fluxz)
      deallocate(qlz,qrz,fluxz)
    end if
  end subroutine finalize_rhs_accelerator

  ! Backend por caras: cada iteración reconstruye y resuelve una interfaz.
  ! No hay arreglos de longitud nx/ny/nz privados por hilo, como en OpenMP.
  subroutine calc_rhs_accelerator(p_in,rhs_out)
    real*8, intent(in)  :: p_in(-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost,neq)
    real*8, intent(out) :: rhs_out(1:nx,1:ny,1:nz,neq)
    integer :: i,j,k,m,offset
    real*8 :: stencil(neq,-2:3),q_left(neq),q_right(neq)
    real*8 :: alpha,beta(3),g(3,3),sqg,src_term(neq),flux_local(neq)

    !$acc parallel loop collapse(3) gang vector present(p_in,qlx,qrx) &
    !$acc& private(stencil,q_left,q_right)
    do k=1,nz
      do j=1,ny
        do i=0,nx
          do offset=-2,2
            do m=1,neq
              stencil(m,offset)=p_in(i+offset,j,k,m)
            end do
          end do
          if (rec_method_id == REC_MP5 .or. rec_method_id == REC_WENO5) then
            do m=1,neq
              stencil(m,3)=p_in(i+3,j,k,m)
            end do
          else
            stencil(:,3)=0.0d0
          end if
          call reconstruct_face_state(stencil,i,nx,DIR_X,rec_method_id,tvd_limiter_id, &
               use_shock_sensor,curved_metric,rho_floor,p_floor,q_left,q_right)
          do m=1,neq
            qlx(i,j,k,m)=q_left(m)
            qrx(i,j,k,m)=q_right(m)
          end do
        end do
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop collapse(3) gang vector present(qlx,qrx,fluxx,alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x,x_face) &
    !$acc& private(q_left,q_right,alpha,beta,g,sqg,flux_local)
    do k=1,nz
      do j=1,ny
        do i=0,nx
          if (cylindrical_geometry .and. abs(x_face(i)) <= 64.0d0*epsilon(1.0d0)) then
            flux_local=0.0d0
          else
            do m=1,neq
              q_left(m)=qlx(i,j,k,m)
              q_right(m)=qrx(i,j,k,m)
            end do
            alpha=alpha_f_x(i+1,j,k)
            beta=beta_f_x(:,i+1,j,k)
            g=gamma_f_x(:,:,i+1,j,k)
            sqg=sqrt_gamma_f_x(i+1,j,k)
            call resolve_riemann_problem(q_left,q_right,DIR_X,alpha,beta,g,sqg,flux_local)
          end if
          do m=1,neq
            fluxx(i,j,k,m)=flux_local(m)
          end do
        end do
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop collapse(3) gang vector &
    !$acc& present(p_in,rhs_out,fluxx,alpha_c,beta_c,gamma_c,gmunu_c) &
    !$acc& present(sqrt_gamma_c,chris_c,dg_c,dlna_c) private(q_left,src_term)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          do m=1,neq
            q_left(m)=p_in(i,j,k,m)
          end do
          call calc_sources(q_left,src_term,i,j,k)
          do m=1,neq
            rhs_out(i,j,k,m)=-(fluxx(i,j,k,m)-fluxx(i-1,j,k,m))/dx+src_term(m)
          end do
        end do
      end do
    end do
    !$acc end parallel loop

    if (ny > 1) then
      !$acc parallel loop collapse(3) gang vector present(p_in,qly,qry) &
      !$acc& private(stencil,q_left,q_right)
      do k=1,nz
        do j=0,ny
          do i=1,nx
            do offset=-2,2
              do m=1,neq
                stencil(m,offset)=p_in(i,j+offset,k,m)
              end do
            end do
            if (rec_method_id == REC_MP5 .or. rec_method_id == REC_WENO5) then
              do m=1,neq
                stencil(m,3)=p_in(i,j+3,k,m)
              end do
            else
              stencil(:,3)=0.0d0
            end if
            call reconstruct_face_state(stencil,j,ny,DIR_Y,rec_method_id,tvd_limiter_id, &
                 use_shock_sensor,curved_metric,rho_floor,p_floor,q_left,q_right)
            do m=1,neq
              qly(i,j,k,m)=q_left(m)
              qry(i,j,k,m)=q_right(m)
            end do
          end do
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop collapse(3) gang vector present(qly,qry,fluxy,alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y,y_face) &
      !$acc& private(q_left,q_right,alpha,beta,g,sqg,flux_local)
      do k=1,nz
        do j=0,ny
          do i=1,nx
            if (polar_geometry .and. abs(sin(y_face(j))) <= 64.0d0*epsilon(1.0d0)) then
              flux_local=0.0d0
            else
              do m=1,neq
                q_left(m)=qly(i,j,k,m)
                q_right(m)=qry(i,j,k,m)
              end do
              alpha=alpha_f_y(i,j+1,k)
              beta=beta_f_y(:,i,j+1,k)
              g=gamma_f_y(:,:,i,j+1,k)
              sqg=sqrt_gamma_f_y(i,j+1,k)
              call resolve_riemann_problem(q_left,q_right,DIR_Y,alpha,beta,g,sqg,flux_local)
            end if
            do m=1,neq
              fluxy(i,j,k,m)=flux_local(m)
            end do
          end do
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop collapse(3) gang vector present(rhs_out,fluxy)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            do m=1,neq
              rhs_out(i,j,k,m)=rhs_out(i,j,k,m)-(fluxy(i,j,k,m)-fluxy(i,j-1,k,m))/dy
            end do
          end do
        end do
      end do
      !$acc end parallel loop
    end if

    if (nz > 1) then
      !$acc parallel loop collapse(3) gang vector present(p_in,qlz,qrz) &
      !$acc& private(stencil,q_left,q_right)
      do k=0,nz
        do j=1,ny
          do i=1,nx
            do offset=-2,2
              do m=1,neq
                stencil(m,offset)=p_in(i,j,k+offset,m)
              end do
            end do
            if (rec_method_id == REC_MP5 .or. rec_method_id == REC_WENO5) then
              do m=1,neq
                stencil(m,3)=p_in(i,j,k+3,m)
              end do
            else
              stencil(:,3)=0.0d0
            end if
            call reconstruct_face_state(stencil,k,nz,DIR_Z,rec_method_id,tvd_limiter_id, &
                 use_shock_sensor,curved_metric,rho_floor,p_floor,q_left,q_right)
            do m=1,neq
              qlz(i,j,k,m)=q_left(m)
              qrz(i,j,k,m)=q_right(m)
            end do
          end do
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop collapse(3) gang vector present(qlz,qrz,fluxz,alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z) &
      !$acc& private(q_left,q_right,alpha,beta,g,sqg,flux_local)
      do k=0,nz
        do j=1,ny
          do i=1,nx
            do m=1,neq
              q_left(m)=qlz(i,j,k,m)
              q_right(m)=qrz(i,j,k,m)
            end do
            alpha=alpha_f_z(i,j,k+1)
            beta=beta_f_z(:,i,j,k+1)
            g=gamma_f_z(:,:,i,j,k+1)
            sqg=sqrt_gamma_f_z(i,j,k+1)
            call resolve_riemann_problem(q_left,q_right,DIR_Z,alpha,beta,g,sqg,flux_local)
            do m=1,neq
              fluxz(i,j,k,m)=flux_local(m)
            end do
          end do
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop collapse(3) gang vector present(rhs_out,fluxz)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            do m=1,neq
              rhs_out(i,j,k,m)=rhs_out(i,j,k,m)-(fluxz(i,j,k,m)-fluxz(i,j,k-1,m))/dz
            end do
          end do
        end do
      end do
      !$acc end parallel loop
    end if
  end subroutine calc_rhs_accelerator

  ! ====================================================================
  ! Cálculo del Paso de Tiempo Dinámico (CFL Adaptativo)
  ! Evalúa la velocidad de la luz coordenada y el transporte del fluido
  ! para garantizar la estabilidad estricta del integrador RK3.
  ! ====================================================================
  subroutine update_adaptive_dt()
    use variables
    implicit none

    real*8 :: max_vx, max_vy, max_vz, inverse_dt
    real*8 :: local_vx, local_vy, local_vz, alpha, beta(3), ginv(3,3)
    integer :: i, j, k, invalid_cfl_cells

    max_vx = 0.0d0
    max_vy = 0.0d0
    max_vz = 0.0d0
    invalid_cfl_cells = 0

    if (.not. ieee_is_finite(dx) .or. dx <= 0.0d0 .or. &
        (ny > 1 .and. (.not. ieee_is_finite(dy) .or. dy <= 0.0d0)) .or. &
        (nz > 1 .and. (.not. ieee_is_finite(dz) .or. dz <= 0.0d0))) then
      write(*,*) 'CRITICAL ERROR: invalid grid spacing in multidimensional CFL condition.'
      error stop 1
    end if

    ! Cotas causales máximas en cada dirección. El cono de luz es una cota
    ! segura tanto para GRHD como para las futuras ondas rápidas de GRMHD.
    !$OMP PARALLEL DO PRIVATE(i, j, k, alpha, beta, ginv, local_vx, local_vy, local_vz) &
    !$OMP REDUCTION(max:max_vx, max_vy, max_vz) REDUCTION(+:invalid_cfl_cells) COLLAPSE(3)
    !$ACC PARALLEL LOOP COLLAPSE(3) REDUCTION(MAX:max_vx,max_vy,max_vz) &
    !$ACC& REDUCTION(+:invalid_cfl_cells) PRESENT(alpha_c,beta_c,gamma_inv_c) &
    !$ACC& PRIVATE(alpha,beta,ginv,local_vx,local_vy,local_vz)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

          ! Leer métrica INVERSA (contravariante) en el centro de la celda
          alpha    = alpha_c(i,j,k)
          beta(:)  = beta_c(:,i,j,k)
          ginv(:,:) = gamma_inv_c(:,:,i,j,k) ! gamma^ij

          if (.not. ieee_is_finite(alpha) .or. alpha <= 0.0d0 .or. &
              .not. all(ieee_is_finite(beta)) .or. &
              .not. ieee_is_finite(ginv(1,1)) .or. ginv(1,1) <= 0.0d0 .or. &
              (ny > 1 .and. (.not. ieee_is_finite(ginv(2,2)) .or. ginv(2,2) <= 0.0d0)) .or. &
              (nz > 1 .and. (.not. ieee_is_finite(ginv(3,3)) .or. ginv(3,3) <= 0.0d0))) then
            invalid_cfl_cells = invalid_cfl_cells + 1
            cycle
          end if

          ! Límite Causal Absoluto Físico (Cono de Luz Coordenado: ds^2 = 0)
          ! La velocidad coordenada de la luz en Kerr usa el componente gamma^ii
          local_vx = alpha * sqrt(ginv(1,1)) + abs(beta(1))
          max_vx = max(max_vx, local_vx)
          if (ny > 1) then
            local_vy = alpha * sqrt(ginv(2,2)) + abs(beta(2))
            max_vy = max(max_vy, local_vy)
          end if
          if (nz > 1) then
            local_vz = alpha * sqrt(ginv(3,3)) + abs(beta(3))
            max_vz = max(max_vz, local_vz)
          end if
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    !$ACC END PARALLEL LOOP

    if (invalid_cfl_cells > 0) then
      write(*,*) 'CRITICAL ERROR: invalid metric data in CFL cells: ', invalid_cfl_cells
      error stop 1
    end if

    ! Para un operador multidimensional, las contribuciones de estabilidad se
    ! suman. Usar min(dx/lambda_i) sólo aplica la restricción unidimensional y
    ! sobreestima dt por hasta el número de direcciones activas.
    inverse_dt = max_vx / dx
    if (ny > 1) inverse_dt = inverse_dt + max_vy / dy
    if (nz > 1) inverse_dt = inverse_dt + max_vz / dz

    if (.not. ieee_is_finite(inverse_dt) .or. inverse_dt <= 0.0d0) then
      write(*,*) 'CRITICAL ERROR: invalid multidimensional CFL propagation rate.'
      error stop 1
    end if
    dt = CFL / inverse_dt

  end subroutine update_adaptive_dt

end module evolution
