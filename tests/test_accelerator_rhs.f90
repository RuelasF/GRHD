program test_accelerator_rhs
  use variables
  use evolution
  use fluxes
  implicit none
  real*8, allocatable :: rhs_cpu(:,:,:,:),rhs_acc(:,:,:,:)
  integer :: i,j,k,m,method,solver,geometry_mode,failures,loc(4)

  failures=0
  nx=9; ny=1; nz=8; nghost=3
  dx=0.1d0; dy=1.0d0; dz=0.2d0
  adb_idx=4.0d0/3.0d0; g1=adb_idx/(adb_idx-1.0d0)
  rho_floor=1.0d-10; p_floor=1.0d-13
  metric_type='Minkowski'; geom_type='Cartesian'
  riemann_solver_id=RS_HLLE; tvd_limiter_id=LIM_MC
  use_shock_sensor=.true.

  allocate(p(-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost,neq))
  allocate(rhs_cpu(nx,ny,nz,neq),rhs_acc(nx,ny,nz,neq))
  allocate(x_face(0:nx),y_face(0:ny),z_face(0:nz))
  allocate(alpha_c(nx,ny,nz),beta_c(3,nx,ny,nz),gamma_c(3,3,nx,ny,nz))
  allocate(gamma_inv_c(3,3,nx,ny,nz),gmunu_c(0:3,0:3,nx,ny,nz))
  allocate(sqrt_gamma_c(nx,ny,nz),chris_c(0:3,0:3,0:3,nx,ny,nz))
  allocate(dg_c(0:3,0:3,1:3,nx,ny,nz),dlna_c(0:3,nx,ny,nz))
  allocate(alpha_f_x(nx+1,ny,nz),beta_f_x(3,nx+1,ny,nz))
  allocate(gamma_f_x(3,3,nx+1,ny,nz),sqrt_gamma_f_x(nx+1,ny,nz))
  allocate(alpha_f_z(nx,ny,nz+1),beta_f_z(3,nx,ny,nz+1))
  allocate(gamma_f_z(3,3,nx,ny,nz+1),sqrt_gamma_f_z(nx,ny,nz+1))

  alpha_c=1.0d0; beta_c=0.0d0; gamma_c=0.0d0; gamma_inv_c=0.0d0
  gmunu_c=0.0d0; sqrt_gamma_c=1.0d0; chris_c=0.0d0; dg_c=0.0d0; dlna_c=0.0d0
  alpha_f_x=1.0d0; beta_f_x=0.0d0; gamma_f_x=0.0d0; sqrt_gamma_f_x=1.0d0
  alpha_f_z=1.0d0; beta_f_z=0.0d0; gamma_f_z=0.0d0; sqrt_gamma_f_z=1.0d0
  do m=1,3
    gamma_c(m,m,:,:,:)=1.0d0; gamma_inv_c(m,m,:,:,:)=1.0d0
    gamma_f_x(m,m,:,:,:)=1.0d0; gamma_f_z(m,m,:,:,:)=1.0d0
  end do
  do i=0,nx; x_face(i)=dble(i)*dx; end do
  y_face=0.0d0
  do k=0,nz; z_face(k)=dble(k)*dz; end do

  do k=-nghost,nz+nghost; do j=-nghost,ny+nghost; do i=-nghost,nx+nghost
    p(i,j,k,eq_de)=1.0d0+0.03d0*sin(0.3d0*dble(i))+0.02d0*cos(0.4d0*dble(k))
    p(i,j,k,eq_pr)=0.5d0+0.02d0*cos(0.2d0*dble(i-k))
    p(i,j,k,eq_vx)=0.08d0*sin(0.1d0*dble(i))
    p(i,j,k,eq_vy)=0.0d0
    p(i,j,k,eq_vz)=0.05d0*cos(0.2d0*dble(k))
  end do; end do; end do

  call initialize_rhs_accelerator(.true.)
  do geometry_mode=0,1
    if (geometry_mode == 0) then
      metric_type='Minkowski'; geom_type='Cartesian'
    else
      metric_type='Eddington-Finkelstein'; geom_type='Spherical'
    end if
    call init_wavespeed_solver()
    do solver=RS_HLLE,RS_HLLC
      riemann_solver_id=solver
      do method=REC_GODUNOV,REC_WENO5
        rec_method_id=method
        call calc_rhs_cpu(p,rhs_cpu)
        call calc_rhs_accelerator(p,rhs_acc)
        if (maxval(abs(rhs_cpu-rhs_acc)) > 2.0d-13) then
          failures=failures+1
          write(*,'(A,3(1X,I0),1X,ES14.6)') 'FAIL accelerator RHS geometry/solver/method', &
               geometry_mode,solver,method,maxval(abs(rhs_cpu-rhs_acc))
          loc=maxloc(abs(rhs_cpu-rhs_acc))
          write(*,'(A,4(1X,I0),2(1X,ES14.6))') '  max at',loc, &
               rhs_cpu(loc(1),loc(2),loc(3),loc(4)),rhs_acc(loc(1),loc(2),loc(3),loc(4))
        end if
      end do
    end do
  end do
  call finalize_rhs_accelerator(.true.)

  if (failures /= 0) then
    write(*,'(A,I0)') 'ACCELERATOR RHS TESTS FAILED: ',failures
    error stop 1
  end if
  write(*,'(A)') 'ACCELERATOR RHS TESTS PASSED'
end program test_accelerator_rhs
