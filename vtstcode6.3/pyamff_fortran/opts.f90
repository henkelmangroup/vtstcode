MODULE opts
  USE adam
  USE lbfgs_ml
  USE rprop
  IMPLICIT NONE
 
  CONTAINS

  SUBROUTINE opt_init(opt_type)
    IMPLICIT NONE
    !Inputs
    CHARACTER(*) :: opt_type

    IF ((opt_type == 'adam') .or. (opt_type == 'ADAM')) THEN 
      CALL adam_init
    ELSE IF ((opt_type == 'lbfgs') .or. (opt_type == 'LBFGS')) THEN
      CALL lbfgs_init
    ELSE IF ((opt_type == 'rprop') .or. (opt_type == 'RPROP')) THEN
      CALL rprop_init
    ELSE  
      print *, 'Fortran error: Not available optimizer type:', opt_type
      print*, 'vasp fails here'
      STOP
    END IF

  END SUBROUTINE
  
  SUBROUTINE opt_step(opt_type,epoch)
    !All default values follow pytorch
    IMPLICIT NONE
    !Inputs
    CHARACTER(*) :: opt_type
    INTEGER :: epoch
    !Variables
    !DOUBLE PRECISION :: beta1, beta2, learning_rate, weight_decay, eps, time

    IF ((opt_type == 'adam') .or. (opt_type == 'ADAM')) THEN
      !Call adam step
      CALL adam_step(epoch)
    ELSE IF ((opt_type == 'lbfgs') .or. (opt_type == 'LBFGS')) THEN
      CALL lbfgs_step(epoch)
    ELSE IF ((opt_type == 'rprop') .or. (opt_type == 'RPROP')) THEN
      CALL rprop_step(epoch)   
    ELSE
      print *, 'Fortran error: Not available optimizer type:', opt_type
      STOP
    END IF

  END SUBROUTINE

  SUBROUTINE opt_cleanup(opt_type) 
    IMPLICIT NONE
    !Inputs
    CHARACTER(*) :: opt_type

    IF (opt_type == 'adam' .or. (opt_type == 'ADAM')) THEN
      CALL adam_cleanup
    ELSE IF ((opt_type == 'lbfgs') .or. (opt_type == 'LBFGS')) THEN
      CALL lbfgs_cleanup
    ELSE IF ((opt_type == 'rprop') .or. (opt_type == 'RPROP')) THEN
      CALL rprop_cleanup
    ELSE  
      print *, 'Fortran error: Not available optimizer type:', opt_type
      STOP
    END IF

  END SUBROUTINE
END MODULE
