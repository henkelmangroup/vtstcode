MODULE neuralnetwork
    USE fpType
    USE nnType
    USE fbp, only: find_loc_char
    USE trainType!, only: energy_training, all_neurons, img_idx, input_fps, natomsE,&
    !force_training, input_dfps, natomsF, inputF, targetF,& 
    !symbols_img, nneighbors_img, neighborlists_img

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: loadfNNParas, init, forward, nncleanup, nncleanup_ase, nncleanup_optim,&
              nncleanup_atom

    CONTAINS

    SUBROUTINE loadfNNParas(nelement, MAX_FPS, uniq_elements)

        IMPLICIT NONE
        !input
        INTEGER :: nelement, MAX_FPS
        CHARACTER*2, DIMENSION(nelement) :: uniq_elements

        !Variables
        INTEGER :: i, j, k, myid
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: flatten_inweights, flatten_hidweights
        CHARACTER(LEN=30) :: model_type, line
        CHARACTER*2 :: atom_type
        
        ! Set common variables
        nelements = nelement
        !total_natoms = natoms
        
        ! Read weights and biases
        OPEN (11, FILE='mlff.pyamff', status='old')
        DO WHILE (line .NE. "#MachineLearning")
            READ (11,*) line
        ENDDO
        READ (11,*) model_type
        READ (11,*) !Skip #Activation function type
        READ (11,*) actfuncId
        READ (11,*) !Skip command line
        READ (11,*) nhidlayers
        ALLOCATE(nhidneurons(nhidlayers))
        READ (11,*) nhidneurons
        ALLOCATE(flatten_inweights(MAXVAL(nGs)*nhidneurons(1)))
        ALLOCATE(flatten_hidweights(MAXVAL(nhidneurons)*MAXVAL(nhidneurons)))
        !ALLOCATE(weights(max(MAXVAL(nGs),MAXVAL(nhidneurons)),MAXVAL(nhidneurons),nhidlayers+1,nelements))
        !ALLOCATE(biases(1,MAXVAL(nhidneurons),nhidlayers+1,nelements))
        ALLOCATE (in_weights(MAXVAL(nGs),nhidneurons(1),nelements))
        ALLOCATE (in_biases(nhidneurons(1),nelements))
        ALLOCATE (hid_weights(MAXVAL(nhidneurons),MAXVAL(nhidneurons),nhidlayers-1,nelements))
        ALLOCATE (hid_biases(MAXVAL(nhidneurons),nhidlayers-1,nelements))
        ALLOCATE (out_weights(MAXVAL(nhidneurons),1,nelements))
        ALLOCATE (out_biases(nelements)) 
               
        READ (11,*) !Skip #Model structure
        DO i = 1, nelements
            ! Read atom type skipping the first symbol #.
            READ (11,*) atom_type
            ! Find index of corresponding atom_type
            myid=find_loc_char(uniq_elements, atom_type, SIZE(uniq_elements))
            ! Input weights, biases
            READ (11,*) flatten_inweights(1:nGs(myid)*nhidneurons(1))
            in_weights(1:nGs(myid),1:nhidneurons(1),myid)=&
            reshape(flatten_inweights(1:nGs(myid)*nhidneurons(1)),(/nGs(myid),nhidneurons(1)/)) 
            READ (11,*) in_biases(1:nhidneurons(1),myid)
            ! Hidden weights, biases
            DO j = 1, nhidlayers-1
                READ (11,*) flatten_hidweights(1:nhidneurons(j)*nhidneurons(j+1))
                hid_weights(1:nhidneurons(j),1:nhidneurons(j+1),j,myid) = &
                reshape(flatten_hidweights(1:nhidneurons(j)*nhidneurons(j+1)),(/nhidneurons(j), nhidneurons(j+1)/))
                READ (11,*) hid_biases(1:nhidneurons(j+1),j,myid)
            END DO
            ! Out weights, biases
            READ (11,*) out_weights(1:nhidneurons(nhidlayers),1,myid)
            READ (11,*) out_biases(myid)
        END DO
        READ (11,*) !skip command
        READ (11,*) slope, intercept
        CLOSE (11)

    END SUBROUTINE

    SUBROUTINE init
        IMPLICIT NONE
        INTEGER :: i,j,k
        ! allocate gradients, forces according to atomic info
        ALLOCATE (in_gradients(MAXVAL(natoms_arr),MAXVAL(nGs), nhidneurons(1),nelements))
        ALLOCATE (hid_gradients(MAXVAL(natoms_arr),MAXVAL(nhidneurons),MAXVAL(nhidneurons),nhidlayers-1,nelements))
        ALLOCATE (out_gradients(MAXVAL(natoms_arr),nhidneurons(nhidlayers),1,nelements))
        ALLOCATE (forces(3,total_natoms))

    END SUBROUTINE

    SUBROUTINE forward(nneighbors, max_nneighbors, neighborlists, symbols, in_atomidx, &
    ordered_fps, fp_primes, max_nGs, max_natoms_arr, max_hidneurons)

      IMPLICIT NONE
      !f2py INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(aux) :: nGs
      !f2py INTEGER, INTENT(aux) :: total_natoms
      !f2py INTEGER, INTENT(aux) :: nelements
      ! Input parameters
      INTEGER :: max_nneighbors, max_nGs, max_natoms_arr, max_hidneurons
      INTEGER, DIMENSION(max_natoms_arr,nelements) :: in_atomidx
      INTEGER, DIMENSION(total_natoms) :: nneighbors, symbols
      INTEGER, DIMENSION(total_natoms, max_nneighbors) :: neighborlists
      !f2py INTEGER, DIMENSION(total_natoms, MAX(nneighbors)) , INTENT(inout) :: neighborlists
      !DOUBLE PRECISION, DIMENSION(total_natoms,max_nGs) :: in_fps
      DOUBLE PRECISION, DIMENSION(max_natoms_arr,max_nGs,nelements) :: ordered_fps
      DOUBLE PRECISION, DIMENSION(total_natoms,max_nneighbors+1,3,max_nGs) :: fp_primes
      ! Variables
      INTEGER :: i, j, k, m, l, myid, nid
      INTEGER :: neigh, i2, k2
      DOUBLE PRECISION, DIMENSION(nelements) :: Es_pElement
      DOUBLE PRECISION, DIMENSION(total_natoms,max_nGs) :: dEdGs
      DOUBLE PRECISION, DIMENSION(max_natoms_arr,max_hidneurons,nhidlayers+1,nelements) :: fps, fps2

      Etotal = 0.d0
      forces = 0.d0 
      
      DO i = 1, nelements
          ! input to 1st hidden layer
          fps(1:natoms_arr(i),1:nhidneurons(1),1,i)=&
          MATMUL(ordered_fps(1:natoms_arr(i),1:nGs(i),i),in_weights(1:nGs(i),1:nhidneurons(1),i))
          DO k=1, natoms_arr(i)
            fps(k,1:nhidneurons(1),1,i)=fps(k,1:nhidneurons(1),1,i)+in_biases(1:nhidneurons(1),i)
          END DO
          fps(1:natoms_arr(i),1:nhidneurons(1),1,i)=actfunc(&
          fps(1:natoms_arr(i),1:nhidneurons(1),1,i),nhidneurons(1),natoms_arr(i),.FALSE.)

          IF (actfuncId == 'silu') THEN
              ! Get sigmoid values
              fps2(1:natoms_arr(i),1:nhidneurons(1),1,i)=&
              MATMUL(ordered_fps(1:natoms_arr(i),1:nGs(i),i),in_weights(1:nGs(i),1:nhidneurons(1),i))
              DO k=1, natoms_arr(i)
                  fps2(k,1:nhidneurons(1),1,i)=fps2(k,1:nhidneurons(1),1,i)+in_biases(1:nhidneurons(1),i)
              END DO
              fps2(1:natoms_arr(i),1:nhidneurons(1),1,i)=actfunc(&
              fps2(1:natoms_arr(i),1:nhidneurons(1),1,i),nhidneurons(1),natoms_arr(i),.TRUE.)

              ! gradient matrix of the first layer 
              in_gradients(1:natoms_arr(i),1:nGs(i),1:nhidneurons(1),i) = &
              grad(fps(1:natoms_arr(i),1:nhidneurons(1),1,i),in_weights(1:nGs(i),1:nhidneurons(1),i),&
              natoms_arr(i),nGs(i),nhidneurons(1),fps2(1:natoms_arr(i),1:nhidneurons(1),1,i))
              
          ELSE   
              ! gradient matrix of the first layer
              in_gradients(1:natoms_arr(i),1:nGs(i),1:nhidneurons(1),i) = &
              grad(fps(1:natoms_arr(i),1:nhidneurons(1),1,i),in_weights(1:nGs(i),1:nhidneurons(1),i), &
              natoms_arr(i),nGs(i),nhidneurons(1))
          END IF

          !If force training is true 
          IF (force_training .EQV. .TRUE.) THEN
              in_gradients_img(1:natoms_arr(i),1:nGs(i),1:nhidneurons(1),i,epoch_img_idx) = &
              in_gradients(1:natoms_arr(i),1:nGs(i),1:nhidneurons(1),i)
          END IF

          ! from the 1st hidden layer to the last hidden layer
          DO j = 1, nhidlayers-1
              fps(1:natoms_arr(i),1:nhidneurons(j+1),j+1,i)=&
              MATMUL(fps(1:natoms_arr(i),1:nhidneurons(j),j,i),hid_weights(1:nhidneurons(j),1:nhidneurons(j+1),j,i))
              DO k=1, natoms_arr(i)
                  fps(k,1:nhidneurons(j+1),j+1,i)=fps(k,1:nhidneurons(j+1),j+1,i)+hid_biases(1:nhidneurons(j+1),j,i)
              END DO
              fps(1:natoms_arr(i),1:nhidneurons(j+1),j+1,i)=actfunc(&
              fps(1:natoms_arr(i),1:nhidneurons(j+1),j+1,i),nhidneurons(j+1),natoms_arr(i),.FALSE.)

              IF (actfuncId == 'silu') THEN 
                  ! Get sigmoid values
                  fps2(1:natoms_arr(i),1:nhidneurons(j+1),j+1,i)=&
                  MATMUL(fps(1:natoms_arr(i),1:nhidneurons(j),j,i),hid_weights(1:nhidneurons(j),1:nhidneurons(j+1),j,i))
                  DO k=1, natoms_arr(i)
                      fps2(k,1:nhidneurons(j+1),j+1,i)=&
                      fps2(k,1:nhidneurons(j+1),j+1,i)+hid_biases(1:nhidneurons(j+1),j,i)
                  END DO
                  fps2(1:natoms_arr(i),1:nhidneurons(j+1),j+1,i)=actfunc(&
                  fps2(1:natoms_arr(i),1:nhidneurons(j+1),j+1,i),nhidneurons(j+1),natoms_arr(i),.TRUE.)

                  hid_gradients(1:natoms_arr(i),1:nhidneurons(j),1:nhidneurons(j+1),j,i)= &
                  grad(fps(1:natoms_arr(i),1:nhidneurons(j+1),j+1,i),hid_weights(1:nhidneurons(j),1:nhidneurons(j+1),j,i),&
                  natoms_arr(i),nhidneurons(j),nhidneurons(j+1),fps2(1:natoms_arr(i),1:nhidneurons(j+1),j+1,i))
              ELSE          
                  ! gradient matrix from a hidden layer to the next hidden layer
                  hid_gradients(1:natoms_arr(i),1:nhidneurons(j),1:nhidneurons(j+1),j,i) = &
                  grad(fps(1:natoms_arr(i),1:nhidneurons(j+1),j+1,i),hid_weights(1:nhidneurons(j),1:nhidneurons(j+1),j,i), &
                  natoms_arr(i),nhidneurons(j),nhidneurons(j+1))

              END IF
              IF (force_training .EQV. .TRUE.) THEN
                  hid_gradients_img(1:natoms_arr(i),1:nhidneurons(j),1:nhidneurons(j+1),j,i,epoch_img_idx) = &
                  hid_gradients(1:natoms_arr(i),1:nhidneurons(j),1:nhidneurons(j+1),j,i)
              END IF 
          END DO
          ! The last hidden layer to output layer
          fps(1:natoms_arr(i),:1,nhidlayers+1,i)=&
          MATMUL(fps(1:natoms_arr(i),1:nhidneurons(nhidlayers),nhidlayers,i),out_weights(1:nhidneurons(nhidlayers),:1,i))
          DO k=1, natoms_arr(i)
              fps(k,1,nhidlayers+1,i)=fps(k,1,nhidlayers+1,i)+out_biases(i)
          END DO
          
          ! gradient matrix from the last hidden layer to output layer
          DO l = 1, natoms_arr(i)
              out_gradients(l,1:nhidneurons(nhidlayers),1,i) = out_weights(1:nhidneurons(nhidlayers),1,i)
              
              IF (force_training .EQV. .TRUE.) THEN 
                  out_gradients_img(l,1:nhidneurons(nhidlayers),1,i,epoch_img_idx) = &
                  out_gradients(l,1:nhidneurons(nhidlayers),1,i)
              END IF
          END DO
          ! Summing output layer's values over natoms
          Es_pElement(i) = SUM(fps(1:natoms_arr(i),1,nhidlayers+1,i))
          Etotal = Etotal + Es_pElement(i)

          DO k = 1, natoms_arr(i)
              IF (nhidlayers .EQ. 1) THEN
                  dEdGs(in_atomidx(k,i),1:nGs(i)) = &
                  MATMUL(in_gradients(k,1:nGs(i),1:nhidneurons(1),i),out_gradients(k,1:nhidneurons(1),1,i))
              ELSE
                  dEdGs(in_atomidx(k,i),1:nGs(i)) = &
                  chainrule(in_gradients(k,1:nGs(i),1:nhidneurons(1),i),hid_gradients(k,:,:,:,i),out_gradients(k,:,:1,i),i)
              END IF
          END DO
      END DO

      ! Store neuron values of multiple images in memory if training is requested
      IF ((energy_training .EQV. .TRUE.) .OR. (force_training .EQV. .TRUE.)) THEN
          all_neurons(1:max_natoms_arr,:,:,:,epoch_img_idx)=fps
      END IF
      ! Forces
      DO i = 1, total_natoms
          myid = symbols(i)
          forces(1:3,i) = MATMUL(fp_primes(i,1,1:3,1:nGs(myid)),dEdGs(i,1:nGs(myid)))
          DO j = 1, nneighbors(i)
              m = neighborlists(i,j)
              nid = symbols(m)
              forces(1:3,i) = forces(1:3,i) + MATMUL(fp_primes(i,j+1,1:3,1:nGs(nid)),dEdGs(m,1:nGs(nid)))
          END DO
      END DO

      ! Consider slope (Etotal, forces) and intercept (Etotal)
      !Note: an 'IF-THEN' statement needed for different Scaler
      Etotal = Etotal*slope + intercept
      !Etotal = Etotal*slope + intercept * total_natoms !only for STDScaler
      forces = forces*(-slope)

      ! Mai edit :  Cohesive Energy
      IF (use_cohesive_energy .EQV. .TRUE.) THEN
          DO i = 1, nelements
              Etotal = Etotal + natoms_arr(i)*coheEs(i)
          END DO
          !Etotal = Etotal / total_natoms
          !forces = forces * (total_natoms)
      END IF

      !PRINT *, 'Etotal'
      !PRINT *, Etotal
      !DO i = 1, total_natoms
      !  PRINT *, 'Forces at ', i
      !  PRINT *, Forces(1:3,i)
      !END DO

      CONTAINS

      FUNCTION actfunc(tensor, nneuron, natom_,sigmoid_flag) RESULT(values)
          IMPLICIT NONE
          LOGICAL :: sigmoid_flag
          INTEGER :: nneuron, i, j, natom_
          DOUBLE PRECISION, DIMENSION(natom_, nneuron) :: tensor, values

          IF ((actfuncId == 'sigmoid') .OR. (sigmoid_flag .EQV. .TRUE.)) THEN
              DO i = 1, nneuron
                  DO j = 1, natom_
                      values(j,i) = 1.d0/(1.d0+exp(-tensor(j,i)))
                  END DO
              END DO
          ELSE IF (actfuncId == 'tanh') THEN
              DO i = 1, nneuron
                  DO j = 1, natom_
                      values(j,i) = (exp(tensor(j,i))-exp(-tensor(j,i)))/(exp(tensor(j,i))+exp(-tensor(j,i)))
                  END DO
              END DO
          ELSE IF (actfuncId == 'softplus') THEN
              DO i = 1, nneuron
                  DO j = 1, natom_
                      values(j,i) = log(1.d0+exp(tensor(j,i)))
                  END DO
              END DO
          ELSE IF (actfuncId == 'relu') THEN
              DO i = 1, nneuron
                  DO j = 1, natom_
                      values(j,i) = max(0.d0,tensor(j,i))
                  END DO
              END DO
          ELSE IF (actfuncId == 'silu') THEN
              DO i = 1, nneuron
                  DO j = 1, natom_
                      values(j,i) = tensor(j,i)/(1.d0+exp(-tensor(j,i)))
                  END DO
              END DO
          ELSE
              !If unknown activation function is set, exit the program 
              print *, 'Unknown activation function type:',actfuncId, 'Check your mlff.pyamff file'
              stop
          END IF

    END FUNCTION

    FUNCTION grad(actval, weight, i, j, k, actval2) RESULT (gradient)
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(i, k) :: actval
        DOUBLE PRECISION, DIMENSION(i, k),OPTIONAL :: actval2
        DOUBLE PRECISION, DIMENSION(k, k) :: diagonal
        DOUBLE PRECISION, DIMENSION(j, k) :: weight
        DOUBLE PRECISION, DIMENSION(i,j,k) :: gradient
        INTEGER :: a, b, i, j, k 

        IF (actfuncId == 'sigmoid') THEN
            DO a = 1, i
                diagonal = 0.0d0
                DO b = 1, k
                    diagonal(b,b) = actval(a,b)*(1.0d0 - actval(a,b))
                END DO
                gradient(a,1:j,1:k) = MATMUL(weight,diagonal)
            END DO
        ELSE IF (actfuncId == 'tanh') THEN
            DO a = 1, i
                diagonal = 0.0d0
                DO b = 1, k
                    diagonal(b,b) = 1.0d0-actval(a,b)**2 
                END DO
                gradient(a,1:j,1:k) = MATMUL(weight,diagonal)
            END DO
        ELSE IF (actfuncId == 'softplus') THEN
            DO a = 1, i
                diagonal = 0.0d0
                DO b = 1, k
                    diagonal(b,b) = 1.0d0-(1.0d0/exp(actval(a,b)))
                END DO
                gradient(a,1:j,1:k) = MATMUL(weight,diagonal)
            END DO
        ELSE IF (actfuncId == 'relu') THEN
            DO a = 1, i
                diagonal = 0.0d0
                DO b = 1, k
                    IF (actval(a,b) > 0) THEN
                        diagonal(b,b) = 1.0d0
                    ELSE
                        diagonal(b,b) = 0.0d0
                    END IF
                END DO
                gradient(a,1:j,1:k) = MATMUL(weight,diagonal)
            END DO
        ELSE IF (actfuncId == 'silu') THEN
            DO a = 1, i
                diagonal = 0.0d0
                DO b = 1, k
                    diagonal(b,b) = actval(a,b)+actval2(a,b)*(1.0d0-actval(a,b))
                END DO
                gradient(a,1:j,1:k) = MATMUL(weight,diagonal)
            END DO
        END IF

    END FUNCTION

    FUNCTION chainrule(in_gradient, hid_gradient_in, out_gradient, idx) RESULT (dEdG)
        IMPLICIT NONE
        INTEGER j, idx
        DOUBLE PRECISION, DIMENSION(nGs(idx),nhidneurons(1)) :: in_gradient
        DOUBLE PRECISION, DIMENSION(nhidneurons(nhidlayers)) :: out_gradient
        DOUBLE PRECISION, DIMENSION(MAXVAL(nhidneurons),MAXVAL(nhidneurons)) :: hid_gradient
        DOUBLE PRECISION, DIMENSION(MAXVAL(nhidneurons),MAXVAL(nhidneurons),nhidlayers-1) :: hid_gradient_in
        DOUBLE PRECISION, DIMENSION(nGs(idx)) :: dEdG

        hid_gradient(1:nhidneurons(1),1:nhidneurons(2))=&
        hid_gradient_in(1:nhidneurons(1),1:nhidneurons(2),1)
        IF (nhidlayers-2 .EQ. 0) THEN
            CONTINUE
        ELSE
            DO j = 1, nhidlayers-2
                hid_gradient = MATMUL(hid_gradient(1:nhidneurons(1),1:nhidneurons(j+1)),&
                hid_gradient_in(1:nhidneurons(j+1),1:nhidneurons(j+2),j+1)) 
            END DO
        END IF
        dEdG = MATMUL(in_gradient,MATMUL(hid_gradient(1:nhidneurons(1),1:nhidneurons(nhidlayers)),out_gradient))

    END FUNCTION

    END SUBROUTINE

    SUBROUTINE TruncatedNormal
        IMPLICIT NONE
        
    END SUBROUTINE

    SUBROUTINE nncleanup() BIND(C,name='nncleanup')
        USE, INTRINSIC :: iso_c_binding
        ! This subroutine is called only once at the very end of calculation
        IMPLICIT NONE
        !DEALLOCATE (biases)
        !DEALLOCATE (weights)
        DEALLOCATE (nhidneurons)
        DEALLOCATE (natoms_arr)
        DEALLOCATE (nGs)
        DEALLOCATE (in_weights)
        DEALLOCATE (in_biases)
        DEALLOCATE (in_gradients)
        DEALLOCATE (out_weights)
        DEALLOCATE (out_biases)
        DEALLOCATE (out_gradients)
        DEALLOCATE (hid_weights)
        DEALLOCATE (hid_biases)
        DEALLOCATE (hid_gradients)
        DEALLOCATE (forces)
        DEALLOCATE (atom_idx)
        DEALLOCATE (fpminvs)
        DEALLOCATE (fpmaxvs)
        DEALLOCATE (diffs)
        DEALLOCATE (interceptScale)
        DEALLOCATE (magnitude)
        DEALLOCATE (coheEs)
    END SUBROUTINE
    
    SUBROUTINE nncleanup_ase
        IMPLICIT NONE
        DEALLOCATE (natoms_arr)
        DEALLOCATE (in_weights)
        DEALLOCATE (in_biases)
        DEALLOCATE (in_gradients)
        DEALLOCATE (out_weights)
        DEALLOCATE (out_biases)
        DEALLOCATE (out_gradients)
        DEALLOCATE (hid_weights)
        DEALLOCATE (hid_biases)
        DEALLOCATE (hid_gradients)
        DEALLOCATE (forces)
        DEALLOCATE (atom_idx)
        DEALLOCATE (interceptScale)
        DEALLOCATE (magnitude)
    END SUBROUTINE
    
    SUBROUTINE nncleanup_atom
        IMPLICIT NONE
        DEALLOCATE (natoms_arr)
        !print *, 'natoms_arr clean'
        DEALLOCATE (in_gradients)
        !print *, 'ingrad clean'
        DEALLOCATE (out_gradients)
        !print *, 'outgrad'
        DEALLOCATE (hid_gradients)
        !print *, 'hidgrad'
        DEALLOCATE (interceptScale)
        !print *, 'intercept scale'
        DEALLOCATE (magnitude) 
        !print *, 'mag'
        DEALLOCATE (forces)
        !print *, 'force'
    END SUBROUTINE
    
    SUBROUTINE nncleanup_optim
        IMPLICIT NONE
        !DEALLOCATE (biases)
        !DEALLOCATE (weights)
        DEALLOCATE (nhidneurons)
        DEALLOCATE (nGs)
        DEALLOCATE (fpminvs)
        DEALLOCATE (fpmaxvs)
        DEALLOCATE (diffs)
        DEALLOCATE (coheEs)
    END SUBROUTINE

END MODULE

