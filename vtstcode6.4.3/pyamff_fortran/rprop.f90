MODULE rprop
  USE nnType
  USE trainType
  IMPLICIT NONE
  ! define global variables
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: prev_bias_grad, bias_etas, bias_steps, signed_bias_grad
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: prev_weight_grad, weight_etas, weight_steps, signed_weight_grad

  CONTAINS

  SUBROUTINE rprop_init
    IMPLICIT NONE

    ! allocate global variables
    ALLOCATE(signed_bias_grad(MAXVAL(nhidneurons),nhidlayers+1,nelements))
    ALLOCATE(prev_bias_grad(MAXVAL(nhidneurons),nhidlayers+1,nelements))
    ALLOCATE(bias_etas(MAXVAL(nhidneurons),nhidlayers+1,nelements))
    ALLOCATE(bias_steps(MAXVAL(nhidneurons),nhidlayers+1,nelements))
    ALLOCATE(signed_weight_grad(max(MAXVAL(nGs),MAXVAL(nhidneurons)),MAXVAL(nhidneurons),nhidlayers+1,nelements))
    ALLOCATE(prev_weight_grad(max(MAXVAL(nGs),MAXVAL(nhidneurons)),MAXVAL(nhidneurons),nhidlayers+1,nelements))
    ALLOCATE(weight_etas(max(MAXVAL(nGs),MAXVAL(nhidneurons)),MAXVAL(nhidneurons),nhidlayers+1,nelements))
    ALLOCATE(weight_steps(max(MAXVAL(nGs),MAXVAL(nhidneurons)),MAXVAL(nhidneurons),nhidlayers+1,nelements))

    ! zero all tensors
    prev_bias_grad = 0.
    bias_etas = 0.
    bias_steps = 0.
    signed_bias_grad = 0.
    prev_weight_grad = 0.
    weight_etas = 0.
    weight_steps = 0.
    signed_weight_grad = 0.

  END SUBROUTINE

  SUBROUTINE rprop_step(epoch,etaplus_in,etaminus_in,step_size_min_in,step_size_max_in,lr_in)
    ! reference: https://pytorch.org/docs/stable/generated/torch.optim.Rprop.html
    IMPLICIT NONE
    ! input variables
    INTEGER, INTENT(IN) :: epoch
    DOUBLE PRECISION, OPTIONAL :: etaplus_in, etaminus_in, step_size_min_in, step_size_max_in, lr_in
    DOUBLE PRECISION :: etaplus, etaminus, step_size_min, step_size_max, lr
    ! intermediate variables
    INTEGER :: i, l

    ! set defaults to be consistent with python implementation
    IF (PRESENT(lr_in)) THEN
      lr=lr_in 
    ELSE
      lr=0.01d0
    END IF
    IF (PRESENT(etaminus_in)) THEN
      etaminus=etaminus_in
    ELSE
      etaminus=0.5d0
    END IF
    IF (PRESENT(etaplus_in)) THEN
      etaplus=etaplus_in
    ELSE
      etaplus=1.2d0
    END IF
    IF (PRESENT(step_size_min_in)) THEN
      step_size_min=step_size_min_in
    ELSE
      step_size_min=1e-6
    END IF
    IF (PRESENT(step_size_max_in)) THEN
      step_size_max=step_size_max_in
    ELSE
      step_size_max=50d0
    END IF

    ! print*, 'called step for epoch #', epoch
    ! print*, 'bias grad=', bias_grad
    ! print*, 'weight grad=', weight_grad
    ! print*, 'in biases=', in_biases
    ! print*, 'in weights=', in_weights
    ! print*, 'hid biases=', hid_biases
    ! print*, 'hid weights=', hid_weights
    ! print*, 'out biases=', out_biases
    ! print*, 'out weights=', out_weights

    ! assign eta values to each parameter
    bias_etas = bias_grad*prev_bias_grad
    weight_etas = weight_grad*prev_weight_grad
    CALL assign_etas(etaplus,etaminus)

    ! set initial step sizes equal to the learning rate
    IF (epoch == 1) THEN
      do i=1, nelements
        bias_steps(1,nhidlayers+1,i) = lr
        weight_steps(1:nhidneurons(nhidlayers),1,nhidlayers+1,i) = lr
        do l=1, nhidlayers
          bias_steps(1:nhidneurons(l),l,i) = lr
          if (l == 1) then
            weight_steps(1:nGs(i),1:nhidneurons(l),l,i) = lr
          else
            weight_steps(1:nhidneurons(l-1),1:nhidneurons(l),l,i) = lr
          end if
        end do
      end do
    END IF

    ! update step sizes
    bias_steps = bias_steps*bias_etas
    weight_steps = weight_steps*weight_etas
    CALL enforce_bounds(step_size_min,step_size_max)
    ! print*, 'bias steps=', bias_steps
    ! print*, 'weight steps=', weight_steps

    ! apply sign operator to gradient tensors
    CALL get_sign
    ! print*, 'sign(bias grad)=', signed_bias_grad
    ! print*, 'sign(weight grad)=', signed_weight_grad

    ! update parameters
    DO i=1, nelements

      in_weights(1:nGs(i),1:nhidneurons(1),i) = in_weights(1:nGs(i),1:nhidneurons(1),i)&
      - weight_steps(1:nGs(i),1:nhidneurons(1),1,i)*signed_weight_grad(1:nGs(i),1:nhidneurons(1),1,i)

      in_biases(1:nhidneurons(1),i) = in_biases(1:nhidneurons(1),i)&
      - bias_steps(1:nhidneurons(1),1,i)*signed_bias_grad(1:nhidneurons(1),1,i)

      DO l=1, nhidlayers-1
        
        hid_weights(1:nhidneurons(l),1:nhidneurons(l+1),l,i) = &
        hid_weights(1:nhidneurons(l),1:nhidneurons(l+1),l,i) - &
        weight_steps(1:nhidneurons(l),1:nhidneurons(l+1),l+1,i) * &
        signed_weight_grad(1:nhidneurons(l),1:nhidneurons(l+1),l+1,i)

        hid_biases(1:nhidneurons(l+1),l,i) = hid_biases(1:nhidneurons(l+1),l,i) - &
        bias_steps(1:nhidneurons(l+1),l+1,i)*signed_bias_grad(1:nhidneurons(l+1),l+1,i)

      END DO
      
      out_weights(1:nhidneurons(nhidlayers),1,i) = out_weights(1:nhidneurons(nhidlayers),1,i) - &
      weight_steps(1:nhidneurons(nhidlayers),1,nhidlayers+1,i) * &
      signed_weight_grad(1:nhidneurons(nhidlayers),1,nhidlayers+1,i)

      out_biases(i) = out_biases(i) - bias_steps(1,nhidlayers+1,i)*signed_bias_grad(1,nhidlayers+1,i)
    
    END DO
    ! print*, 'in biases=', in_biases
    ! print*, 'in weights=', in_weights
    ! print*, 'hid biases=', hid_biases
    ! print*, 'hid weights=', hid_weights
    ! print*, 'out biases=', out_biases
    ! print*, 'out weights=', out_weights
   
    ! save gradient for the next run
    prev_bias_grad = bias_grad
    prev_weight_grad = weight_grad
    
  END SUBROUTINE

  SUBROUTINE assign_etas(etaplus,etaminus)
    IMPLICIT NONE
    ! inputs
    DOUBLE PRECISION, INTENT(IN) :: etaplus, etaminus
    ! other variables
    INTEGER :: i, j, k, l

    do i=1, nelements
      ! input weights
      do j=1,nGs(i) ! for each column
        do k=1,nhidneurons(1) ! for each row
          if (weight_etas(j,k,1,i) > 0.) then
            weight_etas(j,k,1,i) = etaplus
          elseif (weight_etas(j,k,1,i) < 0.) then
            weight_etas(j,k,1,i) = etaminus
            weight_grad(j,k,1,i) = 0.
          else
            weight_etas(j,k,1,i) = 1d0
          end if
        end do
      end do
  
      ! input biases
      do k=1,nhidneurons(1)
        if (bias_etas(k,1,i) > 0.) then
          bias_etas(k,1,i) = etaplus
        elseif (bias_etas(k,1,i) < 0.) then
          bias_etas(k,1,i) = etaminus
          bias_grad(k,1,i) = 0.
        else
          bias_etas(k,1,i) = 1d0
        end if
      end do
    
      do l=1,nhidlayers-1
        ! hidden weights
        do j=1,nhidneurons(l)
          do k=1,nhidneurons(l+1)
            if (weight_etas(j,k,l+1,i) > 0.) then
              weight_etas(j,k,l+1,i) = etaplus
            elseif (weight_etas(j,k,l+1,i) < 0.) then
              weight_etas(j,k,l+1,i) = etaminus
              weight_grad(j,k,l+1,i) = 0.
            else
              weight_etas(j,k,l+1,i) = 1d0
            end if
          end do
        end do

        ! hidden biases
        do k=1,nhidneurons(l+1)
          if (bias_etas(k,l+1,i) > 0.) then
            bias_etas(k,l+1,i) = etaplus
          elseif (bias_etas(k,l+1,i) < 0.) then
            bias_etas(k,l+1,i) = etaminus
            bias_grad(k,l+1,i) = 0.
          else
            bias_etas(k,l+1,i) = 1d0
          end if
        end do
      end do

      ! output weights
      do j=1,nhidneurons(nhidlayers)
        if (weight_etas(j,1,nhidlayers+1,i) > 0) then
          weight_etas(j,1,nhidlayers+1,i) = etaplus
        elseif (weight_etas(j,1,nhidlayers+1,i) < 0) then
          weight_etas(j,1,nhidlayers+1,i) = etaminus
          weight_grad(j,1,nhidlayers+1,i) = 0.
        else
          weight_etas(j,1,nhidlayers+1,i) = 1d0
        end if
      end do

      ! output bias
      if (bias_etas(1,nhidlayers+1,i) > 0.) then
        bias_etas(1,nhidlayers+1,i) = etaplus
      elseif (bias_etas(1,nhidlayers+1,i) < 0.) then
        bias_etas(1,nhidlayers+1,i) = etaminus
        bias_grad(1,nhidlayers+1,i) = 0.
      else
        bias_etas(1,nhidlayers+1,i) = 1d0
      end if
    end do

  END SUBROUTINE

  SUBROUTINE enforce_bounds(step_size_min,step_size_max)
    IMPLICIT NONE
    ! inputs
    DOUBLE PRECISION, INTENT(IN) :: step_size_min, step_size_max
    ! other variables
    INTEGER :: i, j, k, l

    do i=1, nelements
      ! input weights
      do j=1,nGs(i) ! for each column
        do k=1,nhidneurons(1) ! for each row
          if (weight_steps(j,k,1,i) > step_size_max) then
            weight_steps(j,k,1,i) = step_size_max
          elseif (weight_steps(j,k,1,i) < step_size_min) then
            weight_steps(j,k,1,i) = step_size_min
          end if
        end do
      end do
    
      ! input biases
      do k=1,nhidneurons(1)
        if (bias_steps(k,1,i) > step_size_max) then
          bias_steps(k,1,i) = step_size_max
        elseif (bias_steps(k,1,i) < step_size_min) then
          bias_steps(k,1,i) = step_size_min
        end if
      end do
    
      do l=1,nhidlayers-1
        ! hidden weights
        do j=1,nhidneurons(l)
          do k=1,nhidneurons(l+1)
            if (weight_steps(j,k,l+1,i) > step_size_max) then
              weight_steps(j,k,l+1,i) = step_size_max
            elseif (weight_steps(j,k,l+1,i) < step_size_min) then
              weight_steps(j,k,l+1,i) = step_size_min
            end if
          end do
        end do

        ! hidden biases
        do k=1,nhidneurons(l+1)
          if (bias_steps(k,l+1,i) > step_size_max) then
            bias_steps(k,l+1,i) = step_size_max
          elseif (bias_steps(k,l+1,i) < step_size_min) then
            bias_steps(k,l+1,i) = step_size_min
          end if
        end do
      end do

      ! output weights
      do j=1,nhidneurons(nhidlayers)
        if (weight_steps(j,1,nhidlayers+1,i) > step_size_max) then
          weight_steps(j,1,nhidlayers+1,i) = step_size_max
        elseif (weight_steps(j,1,nhidlayers+1,i) < step_size_min) then
          weight_steps(j,1,nhidlayers+1,i) = step_size_min
        end if
      end do

      ! output bias
      if (bias_steps(1,nhidlayers+1,i) > step_size_max) then
        bias_steps(1,nhidlayers+1,i) = step_size_max
      elseif (bias_steps(1,nhidlayers+1,i) < step_size_min) then
        bias_steps(1,nhidlayers+1,i) = step_size_min
      end if
    end do

  END SUBROUTINE

  SUBROUTINE get_sign
    IMPLICIT NONE
    ! other variables
    INTEGER :: i, j, k, l

    ! reset sign tensors for safety
    signed_bias_grad = 0.
    signed_weight_grad = 0.

    do i=1, nelements
      ! input weights
      do j=1,nGs(i) ! for each column
        do k=1,nhidneurons(1) ! for each row
          if (weight_grad(j,k,1,i) > 0.) then
            signed_weight_grad(j,k,1,i) = 1d0
          elseif (weight_grad(j,k,1,i) < 0.) then
            signed_weight_grad(j,k,1,i) = -1d0
          else
            signed_weight_grad(j,k,1,i) = 0.
          end if
        end do
      end do
  
      ! input biases
      do k=1,nhidneurons(1)
        if (bias_grad(k,1,i) > 0.) then
          signed_bias_grad(k,1,i) = 1d0
        elseif (bias_grad(k,1,i) < 0.) then
          signed_bias_grad(k,1,i) = -1d0
        else
          signed_bias_grad(k,1,i) = 0.
        end if
      end do
    
      do l=1,nhidlayers-1
        ! hidden weights
        do j=1,nhidneurons(l)
          do k=1,nhidneurons(l+1)
            if (weight_grad(j,k,l+1,i) > 0.) then
              signed_weight_grad(j,k,l+1,i) = 1d0
            elseif (weight_grad(j,k,l+1,i) < 0.) then
              signed_weight_grad(j,k,l+1,i) = -1d0
            else
              signed_weight_grad(j,k,l+1,i) = 0.
            end if
          end do
        end do

        ! hidden biases
        do k=1,nhidneurons(l+1)
          if (bias_grad(k,l+1,i) > 0.) then
            signed_bias_grad(k,l+1,i) = 1d0
          elseif (bias_grad(k,l+1,i) < 0.) then
            signed_bias_grad(k,l+1,i) = -1d0
          else
            signed_bias_grad(k,l+1,i) = 0.
          end if
        end do
      end do

      ! output weights
      do j=1,nhidneurons(nhidlayers)
        if (weight_grad(j,1,nhidlayers+1,i) > 0) then
          signed_weight_grad(j,1,nhidlayers+1,i) = 1d0
        elseif (weight_grad(j,1,nhidlayers+1,i) < 0) then
          signed_weight_grad(j,1,nhidlayers+1,i) = -1d0
        else
          signed_weight_grad(j,1,nhidlayers+1,i) = 0.
        end if
      end do

      ! output bias
      if (bias_grad(1,nhidlayers+1,i) > 0.) then
        signed_bias_grad(1,nhidlayers+1,i) = 1d0
      elseif (bias_grad(1,nhidlayers+1,i) < 0.) then
        signed_bias_grad(1,nhidlayers+1,i) = -1d0
      else
        signed_bias_grad(1,nhidlayers+1,i) = 0.
      end if
    end do

  END SUBROUTINE

  SUBROUTINE rprop_cleanup
    IMPLICIT NONE

    DEALLOCATE(signed_weight_grad)
    DEALLOCATE(signed_bias_grad)
    DEALLOCATE(prev_weight_grad)
    DEALLOCATE(prev_bias_grad)
    DEALLOCATE(weight_steps)
    DEALLOCATE(bias_steps)
    DEALLOCATE(weight_etas)
    DEALLOCATE(bias_etas)

  END SUBROUTINE

END MODULE
