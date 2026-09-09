program test_polar_boundaries
  use variables
  use conditions
  implicit none

  real*8, allocatable :: q(:,:,:,:)
  integer :: failures

  failures = 0
  nghost = 2
  nx = 3
  ny = 4
  nz = 6
  rho_floor = 1.0d-10
  p_floor = 1.0d-13
  case_name = 'FishMonc3D'

  allocate(q(neq,-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost))
  call fill_interior(q)
  call set_boundary_conditions(q)
  call check_three_dimensional_poles(q, failures)
  deallocate(q)

  ny = 4
  nz = 1
  case_name = 'FishMoncSag'
  allocate(q(neq,-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost))
  call fill_interior(q)
  call set_boundary_conditions(q)
  call check_axisymmetric_poles(q, failures)
  deallocate(q)

  if (failures /= 0) then
    write(*,'(A,I0)') 'POLAR BOUNDARY TESTS FAILED: ', failures
    error stop 1
  end if
  write(*,'(A)') 'POLAR BOUNDARY TESTS PASSED'

contains

  subroutine fill_interior(state)
    real*8, intent(out) :: state(neq,-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost)
    integer :: i, j, k, component

    state = -huge(1.0d0)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          do component = 1, neq
            state(component,i,j,k) = 1.0d5*component + 1.0d3*i + 10.0d0*j + k
          end do
        end do
      end do
    end do
  end subroutine fill_interior


  subroutine check_three_dimensional_poles(state, failures)
    real*8, intent(in) :: state(neq,-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost)
    integer, intent(inout) :: failures
    integer :: g, k, k_shift, component
    real*8 :: expected

    do k = 1, nz
      k_shift = modulo(k-1+nz/2,nz) + 1
      do g = 1, nghost
        do component = 1, neq
          expected = state(component,1,g,k_shift)
          if (component == eq_vy) expected = -expected
          call assert_equal('3D north pole', state(component,1,1-g,k), expected, failures)

          expected = state(component,1,ny+1-g,k_shift)
          if (component == eq_vy) expected = -expected
          call assert_equal('3D south pole', state(component,1,ny+g,k), expected, failures)
        end do
      end do
    end do

    do g = 1, nghost
      call assert_vector('3D lower periodic phi', state(:,1,2,1-g), &
                         state(:,1,2,nz+1-g), failures)
      call assert_vector('3D upper periodic phi', state(:,1,2,nz+g), &
                         state(:,1,2,g), failures)
    end do
  end subroutine check_three_dimensional_poles


  subroutine check_axisymmetric_poles(state, failures)
    real*8, intent(in) :: state(neq,-nghost:nx+nghost,-nghost:ny+nghost,-nghost:nz+nghost)
    integer, intent(inout) :: failures
    integer :: g, component
    real*8 :: expected

    do g = 1, nghost
      do component = 1, neq
        expected = state(component,1,g,1)
        if (component == eq_vy) expected = -expected
        call assert_equal('sagittal north pole', state(component,1,1-g,1), expected, failures)

        expected = state(component,1,ny+1-g,1)
        if (component == eq_vy) expected = -expected
        call assert_equal('sagittal south pole', state(component,1,ny+g,1), expected, failures)
      end do
    end do
  end subroutine check_axisymmetric_poles


  subroutine assert_equal(label, actual, expected, failures)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: actual, expected
    integer, intent(inout) :: failures

    if (actual /= expected) then
      failures = failures + 1
      write(*,'(A,1X,A,2(1X,ES14.6))') 'FAIL', trim(label), actual, expected
    end if
  end subroutine assert_equal


  subroutine assert_vector(label, actual, expected, failures)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: actual(:), expected(:)
    integer, intent(inout) :: failures

    if (any(actual /= expected)) then
      failures = failures + 1
      write(*,'(A,1X,A)') 'FAIL', trim(label)
    end if
  end subroutine assert_vector

end program test_polar_boundaries
