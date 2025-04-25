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
    INTEGER :: global_dft_call =0 
    CONTAINS

    SUBROUTINE train_init(nAtoms,nelement,uniqElems,filename,seedval) !train_init(nAtoms,nelement,max_fps,uniqElems,filename,seedval)
        USE fpCalc
        IMPLICIT NONE
        !Inputs
        INTEGER :: nAtoms, nelement !, max_fps
        CHARACTER*2, DIMENSION(nelement) :: uniqElems
        !Optional inputs
        CHARACTER(*), OPTIONAL :: filename
        INTEGER, OPTIONAL :: seedval
        !Variables
        INTEGER :: i, j, istat
        CHARACTER*3, DIMENSION(92) :: elementArray
        INTEGER, DIMENSION(nAtoms) :: symbols
        DOUBLE PRECISION, DIMENSION(nelement) :: coeh
        !print *,'line 29 training.f90'
        ! Try opening the input filename 
        OPEN(55, FILE=filename, STATUS='old', IOSTAT=istat)
        IF (istat == 6 .OR. istat == 29) THEN 
            !If file open is failed, load deafult fingerprints  
            !PRINT *, 'Warning: PyAMFF cannot find the input file, ',filename, &
            !', specified in INCAR. Hence default fpParas will be used!'
            !CALL load_default_mlff(nelement, max_fps, uniqElems, seedval)
            print *, 'Error: Input file missing! Current version requires an input file,', filename 
            STOP
        ! Input file (filename) is found 
        ELSE
          !Read *.pyamff
          !print *, 'line 42 training.f90' 
          CALL read_mlff(nelement, filename) !read_mlff(nelement, max_fps, filename, seedval)
        END IF
        !print *, 'line 45 training.f90'
        !backward initiation
        CALL init_backward
        !print *, 'line 48 END of train_init training.f90'
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
        CHARACTER*2, DIMENSION(nelement) :: uniq_elements
        INTEGER, DIMENSION(nAtoms) :: symbols
        DOUBLE PRECISION, DIMENSION(nAtoms,3) :: pos_car
        DOUBLE PRECISION, DIMENSION(3,3) :: cell
        DOUBLE PRECISION :: force_coeff, energy_tol, force_tol, grad_tol
        !Variables
        INTEGER :: i, j, img !, maxneighs
        INTEGER, PARAMETER ::forceEngine = 1, max_neighs=100
        ! DOUBLE PRECISION, DIMENSION(nAtoms, maxfps) :: fps
        DOUBLE PRECISION, DIMENSION(max_natarr, maxfps, nelement) :: ordered_fps
        !DOUBLE PRECISION, DIMENSION(nAtoms, max_neighs, 3, maxfps) :: temp_dfps
        ! DOUBLE PRECISION, DIMENSION(nAtoms, nAtoms, 3, maxfps) :: dfps
        INTEGER, DIMENSION(nAtoms, MAX_NEIGHS) :: neighs
        INTEGER, DIMENSION(nAtoms, nAtoms) :: sub_neighs
        INTEGER, DIMENSION(nAtoms) :: num_neigh, sub_num_neigh
        
        REAL :: start, finish
        !print *, 'line 81 nAtimg_ptr: in training.f90 ', nAtimg_ptr

        !temp_dfps = 0.0
        !dfps = 0.0
        !fps = 0.0

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
        !
        !SUBROUTINE calcfps(nAtoms, pos_car, cell, symbols, nelement,forceEngine, &
        !                 nneigh_incell, neighs_incell, num_neigh, num_cells, neighsDefined)
        !
        
        !print *, 'line 108 training.f90'
        CALL calcfps(nAtoms, pos_car, cell, symbols, nelement, forceEngine, sub_num_neigh,neighs, num_neigh)
        !print *, 'line 110 training.f90'
        !CALL calcfps(nAtoms, pos_car, cell, symbols, maxfps, nelement, forceEngine, &fps, temp_dfps, neighs, num_neigh)
        !print *, 'supersymbols: ',supersymbols
        
        !print *, 'fps calced'
        ! CALL ghost_dfps_correct(nelement, nAtoms, maxfps, MAX_NEIGHS, num_neigh, neighs, &
        ! sub_num_neigh, sub_neighs, temp_dfps, dfps)

        !print *, 'line 119 training.f90'
        !print *, 'supersymbols: ',supersymbols

        !Store nneighbors, neighborlists in TrainImg
        !TrainImg(img_idx)%nneighbors=num_neigh
        !TrainImg(img_idx)%neighborlists=neighs
        TrainImg(img_idx)%nneighbors=sub_num_neigh
        !!!TrainImg(img_idx)%neighborlists=sub_neighs  !!!OLD
        TrainImg(img_idx)%neighborlists=neighs  !!!edit
        
        !Update fprange based on fingerprints
        !CALL normalizeParas(nelement, nAtoms, MAXFPS, TrainImg(img_idx)%natoms_arr, &
        !MAXVAL(TrainImg(img_idx)%natoms_arr), symbols, TrainImg(img_idx)%atom_idx, fps, ordered_fps)
        !print *, 'line 130 training.f90'
        !print *, 'nAtoms: ',nAtoms
        !print *, 'natoms_arr +132 training.f90: ',TrainImg(img_idx)%natoms_arr  
        !natoms_arr(i) = TrainImg(img_idx)%natoms_arr
        !print *, 'symbols in training.f90 before normalizeParas: ',symbols
        !!!CALL normalizeParas(nelement) !!! original 
        CALL update_atomInfo(nAtoms,nelement,MAX_FPS,symbols,TrainImg(img_idx)%natoms_arr)
        !print *, 'line 141 training.f90'
        !print*, 'atom_idx line 141  training.f90: ',atom_idx
        CALL normalizeFPs(nelement,nAtoms,uniq_elements,maxfps,&
        MAXVAL(TrainImg(img_idx)%nneighbors),&
        TrainImg(img_idx)%nneighbors,TrainImg(img_idx)%neighborlists(:,1:MAXVAL(TrainImg(img_idx)%nneighbors)))
        !print *, 'atom_idx line 144  training.f90: ',atom_idx
        TrainImg(img_idx)%atom_idx = atom_idx
        DO i=1,nelement
            !DO j=1,natoms_arr(i) !orignal
            DO j=1,TrainImg(img_idx)%natoms_arr(i)
                !ordered_fps(j,1:nGs(i),i) = fps(atom_idx(j,i),1:nGs(i)) !original
                ordered_fps(j,1:nGs(i),i) = fps(TrainImg(img_idx)%atom_idx(j,i),1:nGs(i))
            END DO
        END DO
        TrainImg(img_idx)%tnAtoms = tnAtoms
        IF (ALLOCATED(TrainImg(img_idx)%supersymbols)) THEN
            DEALLOCATE(TrainImg(img_idx)%supersymbols)
        END IF
        ALLOCATE(TrainImg(img_idx)%supersymbols(tnAtoms))
        
        TrainImg(img_idx)%supersymbols=supersymbols
        !print *, 'line 152 training.f90'
        !Store ordered fps and dfps in TrainImg
        !Cleanup the ghost parts
        CALL atomsCleanup
        TrainImg(img_idx)%input_fps=ordered_fps 
        !TrainImg(img_idx)%calc_fps=ordered_fps !before normalized
        TrainImg(img_idx)%input_dfps(:,1:MAXVAL(sub_num_neigh)+1,:,:)=&
        dfps(:,1:MAXVAL(sub_num_neigh)+1,:,:)
        !IF (ALLOCATED(TrainImg(img_idx)%atom_idx)) THEN
        !    DEALLOCATE (TrainImg(img_idx)%atom_idx)
        !END IF
        !ALLOCATE (TrainImg(img_idx)%atom_idx(nelement))
        !print *, 'line 150 training.f90'
        !TrainImg(img_idx)%calc_dfps(:,1:MAXVAL(num_neigh)+1,:,:)=dfps(:,1:MAXVAL(num_neigh)+1,:,:) !before normalized

  !33    CONTINUE
        
        !Copy not normalized fps/dfps to input_fps/dfps before normalization
        !When fp calc is skipped, it copies not normalized fps to input before
        !normalization 
        !TrainImg(img_idx)%input_fps=TrainImg(img_idx)%calc_fps
        !TrainImg(img_idx)%input_dfps=TrainImg(img_idx)%calc_dfps
 
        ! Call trainer when computations of all images are done
        !print *,'line 155 training.f90'
        !print *, 'nimages training.f90: ',nimages 
        !print *, 'line 179 nAtimg_ptr: in training.f90 ', nAtimg_ptr
       IF (img_idx == nimages) THEN
            DO img=1, nimages
              epoch_img_idx=img
              !Update atomic info of each image for NN 
              !!!CALL update_atomInfo(TrainImg(img)%natoms,nelement,maxfps,TrainImg(img)%natoms_arr) ! original
              !print *, 'img training.f90: ',img
              CALL update_atomInfo(TrainImg(img)%natoms,nelement,maxfps,symbols,TrainImg(img)%natoms_arr)
              !print*, 'atom_idx line 190 training.f90: ',atom_idx
              atom_idx = TrainImg(img)%atom_idx
              !print*, 'atom_idx line 197 training.f90: ',atom_idx
              !print *, 'atomic info updated'
              !Normalize fingerprints over all images
              !print *, 'line 178 training.f90'
              !print *, 'TrainImg(img)%input_dfps line 175 training.f90: ',TrainImg(img)%input_dfps(:,1:MAXVAL(TrainImg(img)%nneighbors)+1,:,:)


              !!!! ORIGINAL FUNCTION CALL WITH ORDERED NORMALIZATION  - currently works (Feb 6/25)
              !CALL normalizeFPs_ordered(nelement,TrainImg(img)%natoms,uniq_elements,maxfps,&
              !MAXVAL(TrainImg(img)%nneighbors),MAXVAL(TrainImg(img)%natoms_arr),&
              !TrainImg(img)%nneighbors,TrainImg(img)%neighborlists(:,1:MAXVAL(TrainImg(img)%nneighbors)), &
              !TrainImg(img)%symbols,TrainImg(img)%input_fps,&
              !TrainImg(img)%input_dfps(:,1:MAXVAL(TrainImg(img)%nneighbors)+1,:,:))
                
              !print *, 'normalized training.f90 done'
              !Calculate energy and forces of each image
              !!SUBROUTINE forward(nneighbors, max_nneighbors, neighborlists, &
              !!!ordered_fps, max_nGs, max_natoms_arr, max_hidneurons)
              IF (ALLOCATED(supersymbols)) THEN
                  DEALLOCATE(supersymbols)
              END IF
              ALLOCATE(supersymbols(trainImg(img)%tnAtoms))
              supersymbols = trainImg(img)%supersymbols
              !print*, 'line 213 training.f90'
              CALL allocate_outputs(trainImg(img)%nAtoms,maxfps,MAXVAL(trainImg(img)%nneighbors))
              !!!
              !fps = TrainImg(img)%input_fps
              dfps = TrainImg(img)%input_dfps(:,1:MAXVAL(TrainImg(img)%nneighbors)+1,:,:)
              !CALL normalizeFPs(nelement,TrainImg(img)%natoms,uniq_elements,maxfps,&
              !!MAXVAL(TrainImg(img)%nneighbors),&
              !TrainImg(img)%nneighbors,TrainImg(img)%neighborlists(:,1:MAXVAL(TrainImg(img)%nneighbors)))
              !print*, 'line 217 training.f90 atom_idx: ',atom_idx
              !print *,'line 214'
              !print *, 'line 215 nAtimg_ptr: in training.f90 ', nAtimg_ptr
              !print*, 'line 223 training.f90 atom_idx: ',atom_idx
              CALL forward(trainImg(img)%nneighbors,MAXVAL(trainImg(img)%nneighbors),&
              trainImg(img)%neighborlists(:,1:MAXVAL(TrainImg(img)%nneighbors)),&
              TrainImg(img)%input_fps, MAXVAL(nGs), MAXVAL(TrainImg(img)%natoms_arr), MAXVAL(nhidneurons))
              !print *, 'line 218'

              !CALL forward_old(TrainImg(img)%nneighbors,MAXVAL(TrainImg(img)%nneighbors),&
              !TrainImg(img)%neighborlists(:,1:MAXVAL(TrainImg(img)%nneighbors)),&
              !TrainImg(img)%symbols,TrainImg(img)%atom_idx,&
              !TrainImg(img)%input_fps,&
              !TrainImg(img)%input_dfps(:,1:MAXVAL(TrainImg(img)%nneighbors)+1,:,:),&
              !MAXVAL(nGs),MAXVAL(trainImg(img)%natoms_arr),MAXVAL(nhidneurons))

              !print *, 'forward done'
              !Set values for arrays required for training
              CALL prepTrain(img)
              !print *, 'train prepared'
              !CALL nncleanup_atom
              !!!!CALL nncleanup !orignal, deallocates nGs and  nHidneurons
              !!!!CALL nncleanup_ase
              !print *, 'line 234 nAtimg_ptr(2): in training.f90 ', nAtimg_ptr(2)
              IF (img .LT. nimages) THEN   !!!if not on last img. clean up all data
                  !print *, 'in if img < nimages: if statement training.f90 +208'
                  !print *, 'nGs line 207 training.f90: ',nGs
                  CALL nncleanup_ase
              ELSE   !!!except if on the last image, keep the weights
                  !print *, 'in else block of img < nimages training.f90'
                  CALL nncleanup_atom
              END IF
              !print *, 'atomic info is cleanup'
              !print *, 'nnclean +205 training.f90'
              !print *, 'line 245 nAtimg_ptr(2): in training.f90 ', nAtimg_ptr(2)
              CALL deallocate_outputs()
              DEALLOCATE(supersymbols)
              !print *, 'line 248 nAtimg_ptr(2): in training.f90 ', nAtimg_ptr(2)
            END DO   
            !Timings
            !print *, 'line 208 training.f90'
            CALL cpu_time(start)
            !print  *,'line 210 training.f90'
            !print *, 'opt_type: ',opt_type
            !print *, 'max_epoch: ',max_epoch
            !print *, 'force coefficient: ',force_coeff
            !print *, 'conv_method: ',conv_method
            !print *, 'energy_tol: ',energy_tol
            !print *, 'force_tol: ',force_tol
            !print *, 'grad_tol: ',grad_tol
            !print *, 'uniq_elements: ',uniq_elements
            !print *, 'MAXVAL(nhidneurons): ',MAXVAL(nhidneurons)
            !print *, 'nGs: ',nGs
            !print *, 'Max Ngs: ',MAXVAL(nGs)
            CALL Trainer(opt_type,max_epoch,MAXVAL(nGs),MAXVAL(nhidneurons),uniq_elements,&
                force_coeff,conv_method,energy_tol,force_tol,grad_tol)
            CALL cpu_time(finish)
            !!! CALLING NNcleanup 
            CALL nncleanup_weights
            !CALL nncleanup_optim
            !print *, 'Training Time: ', finish-start, "seconds"
        END IF
    !print *, 'END of trainExec +183 training.f90'
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
        CHARACTER*2, DIMENSION(nelements) :: uniq_elements 
        !Must be input variables eventually
        DOUBLE PRECISION :: beta1, beta2, eps, weight_decay
        !variables
        LOGICAL :: exist_flag, model_converge
        INTEGER :: epoch, time, i, totnatoms
        DOUBLE PRECISION :: lr
        DOUBLE PRECISION :: energyloss, forceloss, loss
        DOUBLE PRECISION :: energyRMSE, forceRMSE, gradnorm
        !print *, 'line 237 training.f90'
        IF (PRESENT(learningRate)) THEN
            lr=learningRate
        ELSE  
            lr=0.01
        END IF   
        !print *, 'training.f90 line 242'
        ! Initiate model_converge flag
        model_converge=.FALSE.
 
        ! Check uncertainty of model
        CALL calc_model_conv(conv_method,etol,ftol,gtol,fconst,0,&
        energyloss,forceloss,energyRMSE,forceRMSE,gradnorm,model_converge)
        !print *, 'line 249 training.f90'
        IF (model_converge) GOTO 30
        !print *, 'line 251 training.f90'
        DO epoch=1, maxepochs
            curr_epoch=epoch
            !print *, '*********************************'
            !print *, 'epoch=', epoch
            !print *, '*********************************'

            ! If conv_method is RMSE, backward propagation           
            IF (conv_method == 'RMSE') CALL backward(fconst)
            !print *, 'conv_method: ',conv_method
            IF (epoch == 1) THEN
                !print *, 'line 283 training.f90'
                CALL opt_init(opt_type)
                !print *, 'line 285 training.f90'
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
            !print *, 'line 297 training.f90'
            CALL opt_step(opt_type,epoch)
            !print *, 'line 299 training.f90'

            ! Compute energy and forces with updated parameters
            !print *, 'line 302 training.f90'
            DO i=1, nimages
              !CALL update_atomInfo(TrainImg(i)%natoms,nelements,MAX_nGs,TrainImg(i)%natoms_arr)
              CALL update_atomInfo(TrainImg(i)%natoms,nelements,MAX_nGs,TrainImg(i)%symbols,TrainImg(i)%natoms_arr)
              epoch_img_idx=i
              totnatoms=TrainImg(i)%natoms
              !print *, 'line 308 training.f90'
              !!! Some of forward variables should be changed to not common variables
              !CALL forward(TrainImg(i)%nneighbors,MAXVAL(TrainImg(i)%nneighbors),&
              !TrainImg(i)%neighborlists,TrainImg(i)%symbols,TrainImg(i)%atom_idx,&
              !TrainImg(i)%input_fps,TrainImg(i)%input_dfps(:,1:MAXVAL(TrainImg(i)%nneighbors)+1,:,:),&
              !max_nGs, MAXVAL(TrainImg(i)%natoms_arr), max_hidneurons)
              !print *, 'line 312 training.f90'
              !!!CALL forward(TrainImg(i)%nneighbors,MAXVAL(TrainImg(i)%nneighbors),&
              !!!TrainImg(i)%neighborlists,&
              !!!TrainImg(i)%input_fps,&
              !!!max_nGs, MAXVAL(TrainImg(i)%natoms_arr), max_hidneurons)
              CALL forward_old(TrainImg(i)%nneighbors,MAXVAL(TrainImg(i)%nneighbors),&
              TrainImg(i)%neighborlists(:,1:MAXVAL(TrainImg(i)%nneighbors)),&
              TrainImg(i)%symbols,TrainImg(i)%atom_idx,&
              TrainImg(i)%input_fps,&
              TrainImg(i)%input_dfps(:,1:MAXVAL(TrainImg(i)%nneighbors)+1,:,:),&
              MAXVAL(nGs),MAXVAL(trainImg(i)%natoms_arr),MAXVAL(nhidneurons))
              !print *, 'line 325 training.f90'

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
              !CALL nncleanup
              !print *, 'atomic info cleanup for this epoch'
            END DO
            !print *, 'line 308 training.f90'
            ! Check conv of model with updated parameters
            CALL calc_model_conv(conv_method,etol,ftol,gtol,fconst,epoch,&
            energyloss,forceloss,energyRMSE,forceRMSE,gradnorm,model_converge)

            IF (model_converge) GOTO 30

        END DO
        !print *,'line 316 training.f90'
        !PRINT *, 'Maximum number of epochs reached but minimization NOT converged.'
        final_fRMSE=forceRMSE
        final_eRMSE=energyRMSE
        IF (conv_method == 'GRADNORM') final_gradnorm=gradnorm
        ! Print loss in the same unit (not squared)
        final_loss=energyRMSE+forceRMSE*fconst
        final_epoch=epoch-1

   30   CONTINUE
    !print *, 'END OF TRAINER SUBROUTINE, training.f90 +326'
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
      !print *, 'line 343 training.f90'
      ! Calculate loss before epoch starts
      CALL LossFunction(fconst,eloss,floss,loss)
      !print *, 'line 346 training.f90 after LossFunction'
      ! Calculate RMSE values 
      eRMSE=sqrt(eloss/nimages)
      fRMSE=sqrt(floss/nimages)
      
      IF (conv_method == 'RMSE') THEN
          IF (eRMSE < etol .and. fRMSE < ftol) THEN
              final_fRMSE=fRMSE
              final_eRMSE=eRMSE
              final_loss=loss
              final_epoch=epoch
              !PRINT *, 'Minimization converged'
              model_converge=.TRUE.
          END IF
      ELSE IF (conv_method == 'GRADNORM') THEN
          ! Backpropagation
          !print *, 'line 375 training.f90'
          CALL backward(fconst)
          !print *, 'line 377 training.f90'
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
    !print *, 'END OF SUBROUTINE calc_model_conv line 377 training.f90'
    END SUBROUTINE calc_model_conv 
   
    SUBROUTINE calc_gradnorm(gradnorm)
    ! Calculate gradient magnitude to quantify uncertainty of our model.
        IMPLICIT NONE
        DOUBLE PRECISION :: gradnorm

        gradnorm=0.d0
        gradnorm=SUM(bias_grad*bias_grad)+SUM(weight_grad*weight_grad)
        gradnorm=SQRT(gradnorm)

    END SUBROUTINE calc_gradnorm
 
    SUBROUTINE write_mlff(nelement, uniq_elements,num_ml_model)
    !SUBROUTINE write_mlff(nelement, uniq_elements)
      !------------------------------------------------------------------------!
      !This is currently temporary format. Only prints out Model parameters    !
      !Eventually this should go to the fingerprints to write the fingerprint  !
      !parameters as well                                                      !
      !------------------------------------------------------------------------!
      IMPLICIT NONE
      !Inputs
      INTEGER :: nelement
      INTEGER :: num_params
      CHARACTER*2, DIMENSION(nelement) :: uniq_elements
      CHARACTER*50 :: param_format
      !Variables
      INTEGER :: i, l
      INTEGER,OPTIONAL :: num_ml_model
      CHARACTER*20 :: file_name
      IF (PRESENT(num_ml_model)) THEN
         WRITE (file_name,'(A,I0)') 'trained.pyamff',num_ml_model
      ELSE
         file_name='trained.pyamff'
      END IF
      num_params = 80
      !!! Switch to 14F.6 format for floats
      !TODO: file name should beTrainImg(img)%natoms!TODO: make more formatted 
      
      !OPEN(unit=11, file='trained.pyamff', status='unknown')   
      OPEN(unit=11, file=file_name, status='unknown')   
      WRITE(11,'(a)') '#Model Parameters'
      DO i=1, nelement
        WRITE(11,'(A2)') uniq_elements(i)
        !WRITE(11,*) in_weights(1:nGs(i),1:nhidneurons(1),i), '#inputLayer weight'
        !WRITE(11,'(1000(F14.6," "))',advance='NO') in_weights(1:nGs(i),1:nhidneurons(1),i)
        !WRITE(11,*) '#inputLayer weight'
        WRITE(param_format,'(A,I0,A)') '(',num_params,'(1X,E17.11)," #inputLayer weight")'
        WRITE(11,param_format) in_weights(1:nGs(i),1:nhidneurons(1),i)
        !!!WRITE(11,*) '#inputLayer weight'
        !WRITE(11,*) in_biases(1:nhidneurons(1),i), '#inputLayer bias'
        !!!WRITE(11,*) in_biases(1:natoms_arr(i),1:nhidneurons(1),i), '#inputLayer bias'
        WRITE(param_format,'(A,I0,A)') '(',natoms_arr(i)*nhidneurons(1),'(1X,E17.11)," #inputLayer bias")'
        WRITE(11,param_format) in_biases(1:natoms_arr(i),1:nhidneurons(1),i)
        DO l=1, nhidlayers-1
          !WRITE(11,*) hid_weights(1:nhidneurons(l),1:nhidneurons(l+1),l,i),&
          !                      '#hiddenLayer_',l,'weight'
          WRITE(param_format,"(A,I0,A,I0,A)") '(',nhidneurons(l)*nhidneurons(l+1),"(1X,E17.11),' #hiddenLayer_",l, "weight')"
          !print *, 'param_format line 469: ',param_format
          WRITE(11,param_format) hid_weights(1:nhidneurons(l),1:nhidneurons(l+1),l,i)
          !                      '#hiddenLayer_',l,'weight'
          WRITE(param_format,'(A,I0,A,I0,A)') '(',natoms_arr(i)*nhidneurons(l+1),"(1X,E17.11),' #hiddenLayer_",l, "bias')"
          WRITE(11,param_format) hid_biases(1:natoms_arr(i),1:nhidneurons(l+1),l,i)
          !                       '#hiddenLayer_',l,'bias'
        END DO
        WRITE(param_format,'(A,I0,A)') '(',nhidneurons(nhidlayers),'(1X,E17.11)," #outputLayer weight")'
        WRITE(11,param_format) out_weights(1:nhidneurons(nhidlayers),1,i)
        WRITE(param_format,'(A,I0,A)') '(',1,'(1X,E17.11)," #outputLayer bias")'
        WRITE(11,param_format) out_biases(natoms_arr(1),1,i)
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
        !print *, 'line 539 nAtimg_ptr(2): ',nAtimg_ptr(2)
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
            !print *, 'line 555 nAtimg_ptr(2): ',nAtimg_ptr(2)
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
            !print *, 'line 567 nAtimg_ptr(2): ',nAtimg_ptr(2)
          END IF
        END IF
        !print *, 'line 570 nAtimg_ptr(2): ',nAtimg_ptr(2) 
      END IF

    END SUBROUTINE

    SUBROUTINE traincleanup(opt_type)
        USE fpCalc
        IMPLICIT NONE
        CHARACTER(*) :: opt_type
        CHARACTER(2) :: elements
        elements = 'Cu'
        CALL write_mlff(1,elements,global_dft_call)
        global_dft_call =global_dft_call+1
        !CALL nncleanup_atom
        CALL opt_cleanup(opt_type)
        CALL backcleanup

    END SUBROUTINE

END MODULE

