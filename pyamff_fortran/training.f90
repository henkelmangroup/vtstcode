MODULE training
    USE pyamffType
    USE nnType
    USE trainType
    USE neuralnetwork
    USE fnnmodule
    USE lossgrad
    USE opts
    USE normalize

    IMPLICIT NONE
    
    CONTAINS

    SUBROUTINE train_init(nAtoms,nelement,max_fps,uniqElems,filename,seedval)
        USE fpCalc
        IMPLICIT NONE
        !Inputs
        INTEGER :: nAtoms, nelement, max_fps
        CHARACTER*3, DIMENSION(nelement) :: uniqElems
        !Optional inputs
        CHARACTER(*), OPTIONAL :: filename
        INTEGER, OPTIONAL :: seedval
        !Variables
        INTEGER :: i, j, istat
        CHARACTER*3, DIMENSION(92) :: elementArray
        INTEGER, DIMENSION(nAtoms) :: symbols
        DOUBLE PRECISION, DIMENSION(nelement) :: coeh

        ! Try opening the input filename 
        OPEN(55, FILE=filename, STATUS='old', IOSTAT=istat)
        IF (istat == 6 .OR. istat == 29) THEN 
            !If file open is failed, load deafult fingerprints  
            ! 11/16/22: Loading default mlff is broken currently. It should be fixed.  
            PRINT *, 'Warning: PyAMFF cannot find the input file, ',filename, &
            ', specified in INCAR. Hence default fpParas will be used!'
            CALL load_default_mlff(nelement, max_fps, uniqElems, seedval)
            !PRINT *, 'Error: Input file missing! Current version requires an input file,', filename 
            !STOP
        ! Input file (filename) is found 
        ELSE
          !Read *.pyamff 
          CALL read_mlff(nelement, max_fps, filename, seedval)
        END IF
 
        !backward initiation
        CALL init_backward
        
    END SUBROUTINE

    SUBROUTINE trainExec(nAtoms,pos_car,cell,symbols,nelement,uniq_elements,&
                max_natarr,maxfps,opt_type,max_epoch,force_coeff,&
                conv_method,energy_tol,force_tol,grad_tol,&
                newImg,update_idx)
        USE nlist
        USE fpCalc
        USE normalize
        IMPLICIT NONE
        !Inputs
        LOGICAL :: newImg
        CHARACTER(*) :: opt_type, conv_method
        INTEGER :: nAtoms, max_epoch, nelement, maxfps, max_natarr
        INTEGER, OPTIONAL :: update_idx
        CHARACTER*3, DIMENSION(nelement) :: uniq_elements
        INTEGER, DIMENSION(nAtoms) :: symbols
        DOUBLE PRECISION, DIMENSION(nAtoms,3) :: pos_car
        DOUBLE PRECISION, DIMENSION(3,3) :: cell
        DOUBLE PRECISION :: force_coeff, energy_tol, force_tol, grad_tol
        !Variables
        INTEGER :: i, j, img, maxneighs
        INTEGER, PARAMETER ::forceEngine = 1, max_neighs=100
        DOUBLE PRECISION, DIMENSION(nAtoms, maxfps) :: fps
        DOUBLE PRECISION, DIMENSION(max_natarr, maxfps, nelement) :: ordered_fps
        DOUBLE PRECISION, DIMENSION(nAtoms, max_neighs, 3, maxfps) :: temp_dfps
        DOUBLE PRECISION, DIMENSION(nAtoms, nAtoms, 3, maxfps) :: dfps
        INTEGER, DIMENSION(nAtoms, MAX_NEIGHS) :: neighs
        INTEGER, DIMENSION(nAtoms, nAtoms) :: sub_neighs
        INTEGER, DIMENSION(nAtoms) :: num_neigh, sub_num_neigh
        
        REAL :: start, finish

        temp_dfps = 0.0
        dfps = 0.0
        fps = 0.0

        !print *, 'trainExec is called in pyamff side'
        
        !If we want to skip fp calc for the previous image, folllowing code
        !works.
        !fingerprint calculation of each image
        !IF (newImg .EQV. .TRUE.) THEN 
        !  IF (img_idx >= update_idx) THEN 
        !    CALL calcfps(nAtoms, pos_car, cell, symbols, MAXFPs, nelement, forceEngine, &
        !         fps, dfps, neighs, num_neigh)
        !  ELSE 
            !Do nothing (fingerprints of old images are already calculated)
        !    print *, 'img ', img_idx, 'fp calculation is skipped'
        !    GOTO 33  
        !  END IF
        !ELSE !All of initial images' fingerprints should be calculated 
        !  CALL calcfps(nAtoms, pos_car, cell, symbols, MAXFPs, nelement, forceEngine, &
        !       fps, dfps, neighs, num_neigh)
        !END IF
        CALL calcfps(nAtoms, pos_car, cell, symbols, maxfps, nelement, forceEngine, &
        fps, temp_dfps, neighs, num_neigh)
        
        !print *, 'fps calced'
        CALL ghost_dfps_correct(nelement, nAtoms, maxfps, MAX_NEIGHS, num_neigh, neighs, &
        sub_num_neigh, sub_neighs, temp_dfps, dfps)

        !Cleanup the ghost parts
        CALL atomsCleanup

        !Store nneighbors, neighborlists in TrainImg
        !TrainImg(img_idx)%nneighbors=num_neigh
        !TrainImg(img_idx)%neighborlists=neighs
        TrainImg(img_idx)%nneighbors=sub_num_neigh
        TrainImg(img_idx)%neighborlists=sub_neighs
        
        !Update fprange based on fingerprints
        CALL normalizeParas2(nelement, nAtoms, MAXFPS, TrainImg(img_idx)%natoms_arr, &
        MAXVAL(TrainImg(img_idx)%natoms_arr), symbols, TrainImg(img_idx)%atom_idx, fps, ordered_fps)
        
        !Store ordered fps and dfps in TrainImg
        TrainImg(img_idx)%input_fps=ordered_fps 
        !TrainImg(img_idx)%calc_fps=ordered_fps !before normalized
        TrainImg(img_idx)%input_dfps(:,1:MAXVAL(sub_num_neigh)+1,:,:)=&
        dfps(:,1:MAXVAL(sub_num_neigh)+1,:,:)
        !TrainImg(img_idx)%calc_dfps(:,1:MAXVAL(num_neigh)+1,:,:)=dfps(:,1:MAXVAL(num_neigh)+1,:,:) !before normalized

  !33    CONTINUE
        
        !Copy not normalized fps/dfps to input_fps/dfps before normalization
        !When fp calc is skipped, it copies not normalized fps to input before
        !normalization 
        !TrainImg(img_idx)%input_fps=TrainImg(img_idx)%calc_fps
        !TrainImg(img_idx)%input_dfps=TrainImg(img_idx)%calc_dfps
 
        ! Call trainer when computations of all images are done
        IF (img_idx == nimages) THEN
            DO img=1, nimages
              epoch_img_idx=img
              !Update atomic info of each image for NN 
              CALL update_atomInfo2(TrainImg(img)%natoms,nelement,maxfps,TrainImg(img)%natoms_arr)
              !print *, 'atomic info updated'
              !Normalize fingerprints over all images
              CALL normalizeFPs2(nelement,TrainImg(img)%natoms,uniq_elements,maxfps,&
              MAXVAL(TrainImg(img)%nneighbors),MAXVAL(TrainImg(img)%natoms_arr),&
              TrainImg(img)%nneighbors,TrainImg(img)%neighborlists(:,1:MAXVAL(TrainImg(img)%nneighbors)), &
              TrainImg(img)%symbols,TrainImg(img)%input_fps,&
              TrainImg(img)%input_dfps(:,1:MAXVAL(TrainImg(img)%nneighbors)+1,:,:))
              !print *, 'normalized'
              !Calculate energy and forces of each image
              CALL forward(TrainImg(img)%nneighbors,MAXVAL(TrainImg(img)%nneighbors),&
              TrainImg(img)%neighborlists(:,1:MAXVAL(TrainImg(img)%nneighbors)),&
              TrainImg(img)%symbols,TrainImg(img)%atom_idx,&
              TrainImg(img)%input_fps,&
              TrainImg(img)%input_dfps(:,1:MAXVAL(TrainImg(img)%nneighbors)+1,:,:),&
              MAXVAL(nGs),MAXVAL(trainImg(img)%natoms_arr),MAXVAL(nhidneurons))
              !print *, 'forward done'
              !Set values for arrays required for training
              CALL prepTrain(img)
              !print *, 'train prepared'
              CALL nncleanup_atom
              !print *, 'atomic info is cleanup'
            END DO   
            !Timings
            CALL cpu_time(start)
            CALL Trainer(opt_type,max_epoch,MAXVAL(nGs),MAXVAL(nhidneurons),uniq_elements,&
                force_coeff,conv_method,energy_tol,force_tol,grad_tol)
            CALL cpu_time(finish)
            print *, 'Training Time: ', finish-start, "seconds"
        END IF
    END SUBROUTINE
   
    SUBROUTINE Trainer(opt_type,maxepochs,max_nGs,max_hidneurons,uniq_elements,&
    fconst,conv_method,etol,ftol,gtol,&
    learningRate)
        IMPLICIT NONE
        !Inputs 
        CHARACTER(*) :: opt_type, conv_method 
        INTEGER, INTENT(IN) :: maxepochs, max_nGs, max_hidneurons
        DOUBLE PRECISION :: fconst, etol, ftol, gtol
        DOUBLE PRECISION, OPTIONAL :: learningRate
        CHARACTER*3, DIMENSION(nelements) :: uniq_elements 
        !Must be input variables eventually
        DOUBLE PRECISION :: beta1, beta2, eps, weight_decay
        !variables
        LOGICAL :: exist_flag, model_converge
        INTEGER :: epoch, time, i, totnatoms
        DOUBLE PRECISION :: lr
        DOUBLE PRECISION :: energyloss, forceloss, loss
        DOUBLE PRECISION :: energyRMSE, forceRMSE, gradnorm
 
        IF (PRESENT(learningRate)) THEN
            lr=learningRate
        ELSE  
            lr=0.01
        END IF   

        ! Initiate model_converge flag
        model_converge=.FALSE.
 
        ! Check uncertainty of model
        CALL calc_model_conv(conv_method,etol,ftol,gtol,fconst,0,&
        energyloss,forceloss,energyRMSE,forceRMSE,gradnorm,model_converge)

        IF (model_converge) GOTO 30
 
        DO epoch=1, maxepochs
            curr_epoch=epoch
            !print *, '*********************************'
            !print *, 'epoch=', epoch
            !print *, '*********************************'

            ! If conv_method is RMSE, backward propagation           
            IF (conv_method == 'RMSE') CALL backward(fconst)

            IF (epoch == 1) THEN
                CALL opt_init(opt_type)
                !INQUIRE(file='pyamff.log',exist=exist_flag) 
                !IF (exist_flag) THEN 
                !  OPEN(unit=10, file='pyamff.log', status='old',position='append')
                !ELSE 
                !  OPEN(unit=10, file='pyamff.log', status='new')
                !END IF
                !WRITE(10,'(a)') 'epoch    lossValue   EnergyRMSE    ForceRMSE'
            END IF
             
            ! Take an optimization step       
            ! TODO: how to set parameters of each optimizers
            CALL opt_step(opt_type,epoch)

            ! Compute energy and forces with updated parameters
            DO i=1, nimages
              CALL update_atomInfo2(TrainImg(i)%natoms,nelements,MAX_nGs,TrainImg(i)%natoms_arr)
              epoch_img_idx=i
              totnatoms=TrainImg(i)%natoms
              !!! Some of forward variables should be changed to not common variables
              CALL forward(TrainImg(i)%nneighbors,MAXVAL(TrainImg(i)%nneighbors),&
              TrainImg(i)%neighborlists,TrainImg(i)%symbols,TrainImg(i)%atom_idx,&
              TrainImg(i)%input_fps,TrainImg(i)%input_dfps(:,1:MAXVAL(TrainImg(i)%nneighbors)+1,:,:),&
              max_nGs, MAXVAL(TrainImg(i)%natoms_arr), max_hidneurons)
                   
              !Update calculated energy of image i
              IF (energy_training) THEN
                inputE(i)=Etotal
              END IF
              !Update calculated force of image i
              IF (force_training) THEN
                IF (i==1) THEN
                  inputF(1:3,1:totnatoms)=forces(1:3,1:totnatoms)*TrainImg(i)%Free(1:3,1:TrainImg(i)%natoms)
                ELSE
                  inputF(1:3,nAtimg_ptr(i)+1:nAtimg_ptr(i)+totnatoms)=&
                  forces(1:3,1:totnatoms)*TrainImg(i)%Free(1:3,1:TrainImg(i)%natoms)
                END IF  
              END IF
              !cleanup atomic info
              CALL nncleanup_atom
              !print *, 'atomic info cleanup for this epoch'
            END DO
             
            ! Check conv of model with updated parameters
            CALL calc_model_conv(conv_method,etol,ftol,gtol,fconst,epoch,&
            energyloss,forceloss,energyRMSE,forceRMSE,gradnorm,model_converge)

            IF (model_converge) GOTO 30

        END DO

        PRINT *, 'Maximum number of epochs reached but minimization NOT converged.'
        final_fRMSE=forceRMSE
        final_eRMSE=energyRMSE
        IF (conv_method == 'GRADNORM') final_gradnorm=gradnorm
        ! Print loss in the same unit (not squared)
        final_loss=energyRMSE+forceRMSE*fconst
        final_epoch=epoch-1

   30   CONTINUE
 
    END SUBROUTINE 
    
    SUBROUTINE calc_model_conv(conv_method,etol,ftol,gtol,fconst,epoch,&
                eloss,floss,eRMSE,fRMSE,gradnorm,model_converge)
      IMPLICIT NONE
      ! Inputs
      INTEGER :: epoch
      CHARACTER(*) :: conv_method
      DOUBLE PRECISION :: etol, ftol, gtol, fconst
      ! Outputs
      LOGICAL :: model_converge
      DOUBLE PRECISION :: eloss, floss, loss, eRMSE, fRMSE 
      DOUBLE PRECISION :: gradnorm

      ! Initiate model_converge flag
      model_converge=.FALSE.

      ! Calculate loss before epoch starts
      CALL LossFunction(fconst,eloss,floss,loss)
      ! Calculate RMSE values 
      eRMSE=sqrt(eloss/nimages)
      fRMSE=sqrt(floss/nimages)
  
      IF (conv_method == 'RMSE') THEN
          IF (eRMSE < etol .and. fRMSE < ftol) THEN
              final_fRMSE=fRMSE
              final_eRMSE=eRMSE
              final_loss=loss
              final_epoch=epoch
              PRINT *, 'Minimization converged'
              model_converge=.TRUE.
          END IF
      ELSE IF (conv_method == 'GRADNORM') THEN
          ! Backpropagation
          CALL backward(fconst)
          ! Calculate gradient magnitude (norm)
          CALL calc_gradnorm(gradnorm)
          ! Check if the convergence criterion is met before taking an optimization step
          IF (gradnorm < gtol) THEN
              final_fRMSE=fRMSE
              final_eRMSE=eRMSE
              final_gradnorm=gradnorm
              ! Print loss in the same unit (not squared)
              final_loss=eRMSE+fRMSE*fconst
              final_epoch=epoch
              PRINT *, 'Minimization converged'
              model_converge=.TRUE.
          END IF
      END IF

    END SUBROUTINE calc_model_conv 
   
    SUBROUTINE calc_gradnorm(gradnorm)
    ! Calculate gradient magnitude to quantify uncertainty of our model.
        IMPLICIT NONE
        DOUBLE PRECISION :: gradnorm

        gradnorm=0.d0
        gradnorm=SUM(bias_grad*bias_grad)+SUM(weight_grad*weight_grad)
        gradnorm=SQRT(gradnorm)

    END SUBROUTINE calc_gradnorm
 
    SUBROUTINE write_mlff(nelement, uniq_elements)
      !------------------------------------------------------------------------!
      !This is currently temporary format. Only prints out Model parameters    !
      !Eventually this should go to the fingerprints to write the fingerprint  !
      !parameters as well                                                      !
      !------------------------------------------------------------------------!
      IMPLICIT NONE
      !Inputs
      INTEGER :: nelement
      CHARACTER*3, DIMENSION(nelement) :: uniq_elements
      !Variables
      INTEGER :: i, l

      !TODO: file name should beTrainImg(img)%natoms!TODO: make more formatted 
      OPEN(unit=11, file='trained.pyamff', status='unknown')   
      WRITE(11,'(a)') '#Model Parameters'
      DO i=1, nelement
        WRITE(11,'(A2)') uniq_elements(i)
        WRITE(11,*) in_weights(1:nGs(i),1:nhidneurons(1),i), '#inputLayer weight'
        WRITE(11,*) in_biases(1:nhidneurons(1),i), '#inputLayer bias'
        DO l=1, nhidlayers-1
          WRITE(11,*) hid_weights(1:nhidneurons(l),1:nhidneurons(l+1),l,i),&
                                 '#hiddenLayer_',l,'weight' 
          WRITE(11,*) hid_biases(1:nhidneurons(l+1),l,i),&
                                 '#hiddenLayer_',l,'bias'
        END DO
        WRITE(11,*) out_weights(1:nhidneurons(nhidlayers),1,i), '#outputLayer weight'  
        WRITE(11,*) out_biases(i), '#outputLayer bias'
      END DO
      WRITE(11,'(a)') '#Energy Scaling Parameters'
      WRITE(11,'(a)') scaler_type
      WRITE(11,*) slope, intercept       
      CLOSE(11) 
    END SUBROUTINE

    SUBROUTINE prepTrain(img)
    !--------------------------------------------------------------------!
    !This subroutine sets common variables required for backward.        !
    !--------------------------------------------------------------------!
      IMPLICIT NONE
      !Inputs
      INTEGER :: img
      !Variables
      INTEGER :: max_nneighbors, max_nGs, totnat
      
      !print *, 'prep train is called'
      max_nneighbors=MAXVAL(TrainImg(img)%nneighbors)
      max_nGs=MAXVAL(nGs)
      totnat=TrainImg(img)%natoms
      !print *, 'required variable is set ' 
      IF ((energy_training .EQV. .TRUE.) .OR. (force_training .EQV. .TRUE.)) THEN
        input_fps(:,:,:,img)=TrainImg(img)%input_fps
        natomsE(img)=TrainImg(img)%natoms
        IF (energy_training .EQV. .TRUE.) THEN
            targetE(img) = TrainImg(img)%targetE
            inputE(img) = Etotal
        END IF
        IF (force_training .EQV. .TRUE.) THEN
            input_dfps(1:totnat,1:max_nneighbors+1,1:3,1:max_nGs,img)=&
            TrainImg(img)%input_dfps(1:totnat,1:max_nneighbors+1,1:3,1:max_nGs) 
          IF (img==1) THEN
            natomsF(1:totnat)=TrainImg(img)%natoms
            targetF(1:3,1:totnat)=TrainImg(img)%targetF(1:3,1:totnat)
            !inputF(1:3,1:totnat)=forces(1:3,1:totnat)
            !Make sure Free atom forces to be zeros
            inputF(1:3,1:totnat)=forces(1:3,1:totnat)*TrainImg(img)%Free(1:3,1:totnat)
            !Initialize the first image's pointer 
            nAtimg_ptr(img)=1
            !Update the next image's pointer if the next image is present 
            IF (img+1 <= nimages) THEN
              nAtimg_ptr(img+1)=TrainImg(img)%natoms
            END IF
          ELSE
            natomsF(nAtimg_ptr(img)+1:nAtimg_ptr(img)+totnat)=TrainImg(img)%natoms
            targetF(1:3,nAtimg_ptr(img)+1:nAtimg_ptr(img)+totnat)=&
            TrainImg(img)%targetF(1:3,1:TrainImg(img)%natoms)
            inputF(1:3,nAtimg_ptr(img)+1:nAtimg_ptr(img)+totnat)=&
            forces(1:3,1:totnat)*TrainImg(img)%Free(1:3,1:TrainImg(img)%natoms)

            !Update the next image's pointer if the next image is present
            IF (img+1 <= nimages) THEN
              nAtimg_ptr(img+1)=nAtimg_ptr(img)+TrainImg(img)%natoms
            END IF
          END IF
        END IF 
      END IF

    END SUBROUTINE

    SUBROUTINE traincleanup(opt_type)
        USE fpCalc
        IMPLICIT NONE
        CHARACTER(*) :: opt_type

        !CALL nncleanup_atom
        CALL opt_cleanup(opt_type)
        CALL backcleanup

    END SUBROUTINE

END MODULE

