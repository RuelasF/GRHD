program test_accelerator_reconstruction
  use variables
  use reconstruction
  implicit none

  integer, parameter :: nodes=8
  real*8 :: line(neq,-3:nodes+3),ql_ref(neq,0:nodes),qr_ref(neq,0:nodes)
  real*8 :: stencil(neq,-2:3),ql_face(neq),qr_face(neq)
  integer :: method,sensor,face,component,offset,failures

  failures=0
  nghost=3
  rho_floor=1.0d-10
  p_floor=1.0d-13
  tvd_limiter_id=LIM_MC
  curved_metric=.false.
  metric_type='Minkowski'

  do face=-3,nodes+3
    line(eq_de,face)=1.0d0+0.02d0*dble(face)+0.08d0*sin(0.7d0*dble(face))
    line(eq_pr,face)=0.4d0+0.01d0*dble(face)+0.03d0*cos(0.9d0*dble(face))
    line(eq_vx,face)=0.05d0*sin(0.3d0*dble(face))
    line(eq_vy,face)=0.02d0*cos(0.4d0*dble(face))
    line(eq_vz,face)=-0.01d0*sin(0.2d0*dble(face))
  end do

  do sensor=0,1
    use_shock_sensor=sensor == 1
    do method=REC_GODUNOV,REC_WENO5
      rec_method_id=method
      call reconstruct_1d_core(line,nodes,ql_ref,qr_ref,DIR_X)
      do face=0,nodes
        do offset=-2,3
          do component=1,neq
            stencil(component,offset)=line(component,face+offset)
          end do
        end do
        call reconstruct_face_state(stencil,face,nodes,DIR_X,method,tvd_limiter_id, &
             use_shock_sensor,curved_metric,rho_floor,p_floor,ql_face,qr_face)
        if (maxval(abs(ql_face-ql_ref(:,face))) > 5.0d-15 .or. &
            maxval(abs(qr_face-qr_ref(:,face))) > 5.0d-15) then
          failures=failures+1
          write(*,'(A,3(1X,I0),2(1X,ES12.4))') 'FAIL face reconstruction', &
               method,sensor,face,maxval(abs(ql_face-ql_ref(:,face))), &
               maxval(abs(qr_face-qr_ref(:,face)))
        end if
      end do
    end do
  end do

  if (failures /= 0) then
    write(*,'(A,I0)') 'ACCELERATOR RECONSTRUCTION TESTS FAILED: ',failures
    error stop 1
  end if
  write(*,'(A)') 'ACCELERATOR RECONSTRUCTION TESTS PASSED'
end program test_accelerator_reconstruction
