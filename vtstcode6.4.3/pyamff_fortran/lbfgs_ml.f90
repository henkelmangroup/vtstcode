MODULE lbfgs_ml
  USE nnType
  USE trainType
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: d, y, s, q_, r, ro, al, flat_grad, prev_flat_grad
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: old_dirs, old_stps
  DOUBLE PRECISION :: step, H_diag
  INTEGER :: nbiases
  INTEGER, DIMENSION(:), ALLOCATABLE :: nweights
  ! INTEGER :: hist_size

  CONTAINS

  SUBROUTINE lbfgs_init(hist_size)
    IMPLICIT NONE
    INTEGER :: history_size, i, l
    INTEGER, OPTIONAL :: hist_size

    IF (PRESENT(hist_size)) THEN
      history_size=hist_size
    ELSE
      history_size=20
    END IF

    ALLOCATE(nweights(nelements))

    ! compute number of biases and weights
    nbiases=nelements*(SUM(nhidneurons)+1)
    nweights=0
    DO i=1,nelements
      nweights(i)=nweights(i)+nGs(i)*nhidneurons(1)
      DO l=1,nhidlayers-1
        nweights(i)=nweights(i)+nhidneurons(l)*nhidneurons(l+1)
      END DO
      nweights(i)=nweights(i)+nhidneurons(nhidlayers)
    END DO
    ! print*, 'nbiases=', nbiases
    ! print*, 'nweights=', nweights

    ALLOCATE(d(nbiases+SUM(nweights)))
    ALLOCATE(y(nbiases+SUM(nweights)))
    ALLOCATE(s(nbiases+SUM(nweights)))
    ALLOCATE(q_(nbiases+SUM(nweights)))
    ALLOCATE(r(nbiases+SUM(nweights)))
    ALLOCATE(flat_grad(nbiases+SUM(nweights)))
    ALLOCATE(prev_flat_grad(nbiases+SUM(nweights)))
    ALLOCATE(old_dirs(nbiases+SUM(nweights),history_size))
    ALLOCATE(old_stps(nbiases+SUM(nweights),history_size)) 
    ALLOCATE(ro(history_size)) 
    ALLOCATE(al(history_size))

    ! zero all tensors
    d=0.
    y=0.
    s=0.
    q_=0.
    r=0.
    ro=0.
    al=0.
    old_dirs=0.
    old_stps=0.
    flat_grad=0.
    prev_flat_grad=0.

    ! print*, 'initialized lbfgs'

  END SUBROUTINE

  SUBROUTINE lbfgs_step(epoch,learning_rate,tol_grad,tol_change,hist_size)
    ! line search currently not implemented
    ! reference: https://pytorch.org/docs/stable/generated/torch.optim.LBFGS.html#torch.optim.LBFGS
    IMPLICIT NONE
    ! inputs
    DOUBLE PRECISION, OPTIONAL :: learning_rate, tol_grad, tol_change
    DOUBLE PRECISION :: lr, tolerance_grad, tolerance_change
    INTEGER, OPTIONAL :: hist_size
    INTEGER :: history_size
    INTEGER, INTENT(IN) :: epoch
    ! intermediate variables
    DOUBLE PRECISION :: ys, be_i, gtd
    INTEGER :: num_old, i, j

    ! set parameters if not specified
    IF (PRESENT(learning_rate)) THEN
      lr=learning_rate 
    ELSE
      lr=0.01d0 ! made to be consistent with default settings in PyAMFF config.yaml
    END IF
    IF (PRESENT(tol_grad)) THEN
      tolerance_grad=tol_grad
    ELSE
      tolerance_grad=1e-10
    END IF
    IF (PRESENT(tol_change)) THEN
      tolerance_change=tol_change
    ELSE
      tolerance_change=1e-10
    END IF
    IF (PRESENT(hist_size)) THEN
      history_size=hist_size
    ELSE
      history_size=20
    END IF

    ! print*, 'called fortran lbfgs step'
    ! print*, 'beginning epoch', epoch

    ! gather gradient information
    CALL gather_flat_grad
    ! print*, 'flat grad=', flat_grad
    
    ! to do? add check for first order optimality with tolerance_grad 

    ! compute direction
    IF (epoch == 1) THEN
      d = -flat_grad
    ELSE 
      y = flat_grad-prev_flat_grad
      ! print*, 'y=', y
      s = d*step
      ! print*, 's=', s
      ys = SUM(y*s)
      ! print*, 'ys=', ys
      IF (ys > 1e-10) THEN
        ! print*, 'passed curvature condition'
      
        ! update memory
        IF (epoch-1 <= history_size) THEN ! if we have enough memory
          old_dirs(:,epoch-1) = y
          old_stps(:,epoch-1) = s
          ro(epoch-1) = (1/ys)
        ELSE ! if we need to replace old memory
          DO j = 1, history_size-1
            old_dirs(:,j) = old_dirs(:,j+1)
            old_stps(:,j) = old_stps(:,j+1)
            ro(j) = ro(j+1)
          END DO
          old_dirs(:,history_size) = y
          old_stps(:,history_size) = s
          ro(history_size) = (1/ys)
        END IF

        ! update scale of initial Hessian approx
        H_diag = ys / SUM(y*y)

      END IF

      IF (history_size > epoch-1) THEN
        num_old = epoch-1
      ELSE
        num_old = history_size
      END IF

      ! compute Hessian*gradient product
      q_ = -flat_grad

      DO i=num_old, 1, -1
        al(i) = SUM(old_stps(:,i)*q_) * ro(i)
        q_ = q_ - old_dirs(:,i)*al(i)
      END DO
      ! print*, 'alphas for epoch #', epoch, al
      ! print*, 'q for epoch #', epoch, q_
      ! print*, 'H_diag=', H_diag
      r = q_*H_diag
      DO i=1, num_old
        be_i = SUM(old_dirs(:,i)*r) * ro(i)
        r = r + old_stps(:,i)*(al(i)-be_i)
      END DO
      d = r
    END IF
    ! print*, 'd=', d

    ! compute step length
    IF (epoch == 1) THEN
      step = MIN(1.0, (1.0d0 / SUM(ABS(flat_grad)))) * lr ! scaled first step to be consistent with pytorch
    ELSE
      step = lr
    END IF
    ! print*, 'step=', step

    ! calculate "directional derivative" - remove for now, add back if we implement an internal loop
    ! gtd = SUM(flat_grad*d) ! dot product
    ! IF (gtd > -tolerance_change) THEN
    ! ...

    ! update weights and biases
    ! if line search was implemented, we would call it here
    CALL update_params

    ! torch recalculates loss/gradients here, but we do that in training.f90

    ! save gradients for the next run
    prev_flat_grad = flat_grad

    ! TODO: add more break criteria (lack of loss progress, etc)

  END SUBROUTINE

  SUBROUTINE gather_flat_grad
    IMPLICIT NONE
    ! intermediate variables
    DOUBLE PRECISION, DIMENSION(nbiases) :: bias_flat_grad
    DOUBLE PRECISION, DIMENSION(SUM(nweights)) :: weight_flat_grad
    INTEGER :: i, j, l, n, offset

    ! print *, 'bias_grad=', bias_grad
    ! print *, 'weight_grad=', weight_grad

    ! takes gradient tensors and flattens them into one nonzero row vector
    ! start with bias grads
    offset=0
    DO i=1,nelements
      DO l=1,nhidlayers
        bias_flat_grad(offset+1:offset+nhidneurons(l))=bias_grad(1:nhidneurons(l),l,i)
        offset=offset+nhidneurons(l)
      END DO
      bias_flat_grad(offset+1)=bias_grad(1,nhidlayers+1,i)
      offset=offset+1
    END DO
    ! print*, 'bias_flat_grad=', bias_flat_grad
  
    ! next get weight grads
    offset=0
    DO i=1,nelements ! for each subnetwork...
      DO n=1,nhidneurons(1) ! strip out the first layer
        weight_flat_grad(offset+1:offset+nGs(i))=weight_grad(1:nGs(i),n,1,i)
        offset=offset+nGs(i)
      END DO
      DO l=1,nhidlayers-1 ! strip out the other layers
        DO n=1,nhidneurons(l+1) ! go through each row
          weight_flat_grad(offset+1:offset+nhidneurons(l))=weight_grad(1:nhidneurons(l),n,l+1,i)
          offset=offset+nhidneurons(l)
        END DO
      END DO ! grab the outputs
      weight_flat_grad(offset+1:offset+nhidneurons(nhidlayers))=weight_grad(1:nhidneurons(nhidlayers),1,nhidlayers+1,i)
      offset=offset+nhidneurons(nhidlayers) 
    END DO
    ! print*, 'weight_flat_grad=', weight_flat_grad 

    ! populate the final gradient row vector
    flat_grad(1:SIZE(bias_flat_grad)) = bias_flat_grad(:)
    flat_grad(SIZE(bias_flat_grad)+1:SIZE(bias_flat_grad)+SIZE(weight_flat_grad)) = weight_flat_grad(:)

  END SUBROUTINE 

  SUBROUTINE update_params
    ! probably a cleaner way to write this involves saving indices from gather_flat_grad, but this works
    IMPLICIT NONE
    ! other variables
    INTEGER :: offset, i, l, k

    ! update input biases
    offset=0
    DO i=1, nelements
      DO k=1, natoms_arr(i)        
        in_biases(k,1:nhidneurons(1),i)=in_biases(k,1:nhidneurons(1),i)&
        +step*d((offset+1):(offset+nhidneurons(1)))
        offset=offset+SUM(nhidneurons)+1
      END DO
      
    END DO

    ! update input weights
    offset = nbiases ! reset offset!
    DO i = 1, nelements
      in_weights(1:nGs(i),1:nhidneurons(1),i)=in_weights(1:nGs(i),1:nhidneurons(1),i)&
      +step*RESHAPE(d(offset+1:offset+nGs(i)*nhidneurons(1)), (/nGs(i),nhidneurons(1)/)) ! reshaping is necessary
      offset=offset+nweights(i)
    END DO

    ! update hidden biases
    offset = 0
    DO i=1, nelements
      DO l=1, nhidlayers-1
        offset=offset+nhidneurons(l)
        DO k=1, natoms_arr(i)        
          hid_biases(k,1:nhidneurons(l+1),l,i)=hid_biases(k,1:nhidneurons(l+1),l,i)&
          +step*d(offset+1:offset+nhidneurons(l+1)) ! no need to reshape here
        END DO
      END DO
      offset=offset+nhidneurons(nhidlayers)+1
    END DO

    ! update hidden weights
    offset = nbiases
    DO i=1, nelements
      DO l=1, nhidlayers-1
        IF (l == 1) THEN
          offset = offset+nGs(i)*nhidneurons(1)
        ELSE
          offset = offset+nhidneurons(l-1)*nhidneurons(l)
        END IF
        hid_weights(1:nhidneurons(l),1:nhidneurons(l+1),l,i) = hid_weights(1:nhidneurons(l),1:nhidneurons(l+1),l,i)&
        + step*RESHAPE(d(offset+1:offset+nhidneurons(l)*nhidneurons(l+1)), (/nhidneurons(l),nhidneurons(l+1)/))
      END DO
      offset = offset+nhidneurons(nhidlayers-1)*nhidneurons(nhidlayers)+nhidneurons(nhidlayers)
    END DO

    ! update output biases
    offset = 0
    DO i=1, nelements
      offset = offset+SUM(nhidneurons)+1
      DO k=1, natoms_arr(i)
        out_biases(k,1,i) = out_biases(k,1,i) + step*d(offset)
      END DO
    END DO

    ! update output weights
    offset = nbiases-nhidneurons(nhidlayers)
    DO i=1, nelements
      offset = offset+nweights(i)
      out_weights(1:nhidneurons(nhidlayers),1,i) = out_weights(1:nhidneurons(nhidlayers),1,i)&
      + step*d(offset+1:offset+2)
    END DO

    ! print*, 'updated parameters:'
    ! print*, 'in_biases=', in_biases
    ! print*, 'in_weights=', in_weights
    ! print*, 'out_biases=', out_biases
    ! print*, 'out_weights=', out_weights

  END SUBROUTINE

  SUBROUTINE lbfgs_cleanup
    IMPLICIT NONE

    DEALLOCATE(d)
    DEALLOCATE(y)
    DEALLOCATE(s)
    DEALLOCATE(q_)
    DEALLOCATE(r)
    DEALLOCATE(ro)
    DEALLOCATE(al)
    DEALLOCATE(nweights)
    DEALLOCATE(old_dirs)
    DEALLOCATE(old_stps) 
    DEALLOCATE(flat_grad)
    DEALLOCATE(prev_flat_grad)

    ! print*, 'cleaned up lbfgs'
    
  END SUBROUTINE

END MODULE
