program test_multidimensional_cfl
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use variables
  use evolution
  implicit none

  integer :: failures

  failures = 0
  CFL = 0.4d0

  call check_cfl_case('1D', 4, 1, 1, 0.5d0, 1.0d0, 1.0d0, &
                      1.0d0, (/ 0.0d0, 0.0d0, 0.0d0 /), &
                      (/ 1.0d0, 1.0d0, 1.0d0 /), &
                      CFL/(1.0d0/0.5d0), failures)

  call check_cfl_case('2D', 4, 3, 1, 0.5d0, 0.25d0, 1.0d0, &
                      1.0d0, (/ 0.0d0, 0.0d0, 0.0d0 /), &
                      (/ 1.0d0, 4.0d0, 1.0d0 /), &
                      CFL/(1.0d0/0.5d0 + 2.0d0/0.25d0), failures)

  call check_cfl_case('3D shifted', 4, 3, 2, 0.5d0, 0.25d0, 0.125d0, &
                      0.8d0, (/ 0.1d0, -0.2d0, 0.05d0 /), &
                      (/ 1.0d0, 4.0d0, 9.0d0 /), &
                      CFL/(0.9d0/0.5d0 + 1.8d0/0.25d0 + 2.45d0/0.125d0), failures)

  if (failures /= 0) then
    write(*,'(A,I0)') 'MULTIDIMENSIONAL CFL TESTS FAILED: ', failures
    error stop 1
  end if
  write(*,'(A)') 'MULTIDIMENSIONAL CFL TESTS PASSED'

contains

  subroutine check_cfl_case(label, nx_in, ny_in, nz_in, dx_in, dy_in, dz_in, &
                            lapse, shift, inverse_diagonal, expected_dt, failures)
    character(len=*), intent(in) :: label
    integer, intent(in) :: nx_in, ny_in, nz_in
    real*8, intent(in) :: dx_in, dy_in, dz_in, lapse
    real*8, intent(in) :: shift(3), inverse_diagonal(3), expected_dt
    integer, intent(inout) :: failures
    integer :: direction

    nx = nx_in
    ny = ny_in
    nz = nz_in
    dx = dx_in
    dy = dy_in
    dz = dz_in

    allocate(alpha_c(nx,ny,nz), beta_c(3,nx,ny,nz), gamma_inv_c(3,3,nx,ny,nz))
    alpha_c = lapse
    beta_c = 0.0d0
    gamma_inv_c = 0.0d0
    do direction = 1, 3
      beta_c(direction,:,:,:) = shift(direction)
      gamma_inv_c(direction,direction,:,:,:) = inverse_diagonal(direction)
    end do

    call update_adaptive_dt()
    if (.not. ieee_is_finite(dt) .or. &
        abs(dt-expected_dt) > 64.0d0*epsilon(1.0d0)*max(1.0d0,abs(expected_dt))) then
      failures = failures + 1
      write(*,'(A,1X,A,2(1X,ES16.8))') 'FAIL', trim(label), dt, expected_dt
    end if

    deallocate(alpha_c, beta_c, gamma_inv_c)
  end subroutine check_cfl_case

end program test_multidimensional_cfl
