! Persistent OpenACC data environment for NVIDIA GPUs. The hydrodynamic state
! remains on the device between RK stages; host transfers occur only for I/O,
! diagnostics, checkpoints, and the one-time perturbation event.
module accelerator
  use variables
  implicit none
  private
  public :: initialize_accelerator, finalize_accelerator, &
            accelerator_update_host, accelerator_update_device, wall_time_seconds

contains

  real*8 function wall_time_seconds()
    integer :: count,count_rate
    call system_clock(count,count_rate)
    wall_time_seconds=dble(count)/dble(count_rate)
  end function wall_time_seconds

  subroutine initialize_accelerator()
    if (.not. accelerator_enabled) return

    print *, 'OpenACC backend enabled: FP64 state and arithmetic.'

    !$acc enter data copyin(nx,ny,nz,nghost,dx,dy,dz,adb_idx,g1,rho_floor,p_floor, &
    !$acc& rec_method_id,tvd_limiter_id,riemann_solver_id,use_shock_sensor, &
    !$acc& use_srhd_wavespeeds,flat_cartesian_sources,curved_metric, &
    !$acc& cylindrical_geometry,polar_geometry,case_id, &
    !$acc& jet_ambient_density,jet_ambient_pressure,jet_ambient_velocity, &
    !$acc& jet_nozzle_radius,jet_density,jet_pressure,jet_velocity)
    !$acc enter data copyin(x,y,z,x_face,y_face,z_face)
    !$acc enter data copyin(u,up,p,rhs)
    !$acc enter data copyin(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c, &
    !$acc& sqrt_gamma_c,chris_c,dg_c,dlna_c)

    !$acc enter data copyin(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
    if (ny > 1) then
      !$acc enter data copyin(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
    end if
    if (nz > 1) then
      !$acc enter data copyin(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
    end if
    if (allocated(michel_injector)) then
      !$acc enter data copyin(michel_injector)
    end if
  end subroutine initialize_accelerator

  subroutine accelerator_update_host()
    if (.not. accelerator_enabled) return
    !$acc update host(u,up,p)
  end subroutine accelerator_update_host

  subroutine accelerator_update_device()
    if (.not. accelerator_enabled) return
    !$acc update device(u,up,p)
  end subroutine accelerator_update_device

  subroutine finalize_accelerator()
    if (.not. accelerator_enabled) return
    if (allocated(michel_injector)) then
      !$acc exit data delete(michel_injector)
    end if
    if (nz > 1) then
      !$acc exit data delete(alpha_f_z,beta_f_z,gamma_f_z,sqrt_gamma_f_z)
    end if
    if (ny > 1) then
      !$acc exit data delete(alpha_f_y,beta_f_y,gamma_f_y,sqrt_gamma_f_y)
    end if
    !$acc exit data delete(alpha_f_x,beta_f_x,gamma_f_x,sqrt_gamma_f_x)
    !$acc exit data delete(alpha_c,beta_c,gamma_c,gamma_inv_c,gmunu_c, &
    !$acc& sqrt_gamma_c,chris_c,dg_c,dlna_c)
    !$acc exit data delete(u,up,p,rhs)
    !$acc exit data delete(x,y,z,x_face,y_face,z_face)
    !$acc exit data delete(nx,ny,nz,nghost,dx,dy,dz,adb_idx,g1,rho_floor,p_floor, &
    !$acc& rec_method_id,tvd_limiter_id,riemann_solver_id,use_shock_sensor, &
    !$acc& use_srhd_wavespeeds,flat_cartesian_sources,curved_metric, &
    !$acc& cylindrical_geometry,polar_geometry,case_id, &
    !$acc& jet_ambient_density,jet_ambient_pressure,jet_ambient_velocity, &
    !$acc& jet_nozzle_radius,jet_density,jet_pressure,jet_velocity)
  end subroutine finalize_accelerator

end module accelerator
