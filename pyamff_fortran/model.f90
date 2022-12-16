MODULE PyAMFF
    USE pyamffType
    USE training
    IMPLICIT NONE
    CHARACTER*2, DIMENSION(:), ALLOCATABLE :: uniqElems
    CONTAINS
    
    SUBROUTINE PyAMFF_init(NAts, nelement, nimg, &
    in_uniqElems, in_natarr, in_maxNimg, filename, seedval, uniqueNrs, atomicNrs)
    !---------------------------------------------------------------------------------!
    !Initiate PyAMFF by allocating ImgInfo data type with maximum setup values.       !
    !Maximum number of images(max_nimg) is set to allow addition of images on the fly !
    !The number of input images(nimg) is the number of initial training images        !
    !This subroutine sets variables (symbols, natoms_arr, uniq_elements,etc)          !
    !of initial training images.                                                      ! 
    !---------------------------------------------------------------------------------!
      IMPLICIT NONE
      !Inputs
      INTEGER :: nimg, nelement
      INTEGER, DIMENSION(nimg) :: NAts
      !Optional inputs
      !TODO: we should change optional to not optional
      INTEGER, OPTIONAL :: in_maxNimg, seedval
      CHARACTER*2, DIMENSION(nelement), OPTIONAL :: in_uniqElems
      INTEGER, DIMENSION(nelement),OPTIONAL :: uniqueNrs
      INTEGER, DIMENSION(nelement, nimg), OPTIONAL :: in_natarr
      INTEGER, DIMENSION(MAXVAL(NAts),nimg),OPTIONAL :: atomicNrs
      CHARACTER(*), OPTIONAL :: filename
      !Variables
      INTEGER, PARAMETER :: max_fps=500
      INTEGER, PARAMETER :: max_neighs=100
      INTEGER :: maxNimg
      INTEGER :: i, j, img, max_natarr, ptr, ptr2
      INTEGER, DIMENSION(nelement) :: idx_arr
      CHARACTER*3, DIMENSION(92) :: elementArray
       
      DATA elementArray / "H","He","Li","Be","B","C","N","O", &
                 "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc", &
                 "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se", &
                 "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag", &
                 "Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd", &
                 "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta", &
                 "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", &
                 "Fr","Ra","Ac","Th","Pa","U" /

      IF (PRESENT(in_maxNimg)) THEN 
        maxNimg=in_maxNimg
      ELSE
        maxNimg=nimg !Default of maximum number of images = number of initial input images
      END IF

      !Allocate arrays 
      ALLOCATE(TrainImg(maxNimg)) !use max_nimg for addition of images later
      ALLOCATE(uniqElems(nelement)) 
      DO img=1, nimg
        ALLOCATE(TrainImg(img)%natoms_arr(nelement))
        ALLOCATE(TrainImg(img)%symbols(NAts(img)))
        ALLOCATE(TrainImg(img)%poscar(NAts(img),3))
        ALLOCATE(TrainImg(img)%cell(3,3))
        ALLOCATE(TrainImg(img)%nneighbors(NAts(img)))
        ALLOCATE(TrainImg(img)%neighborlists(NAts(img),max_neighs))
        ALLOCATE(TrainImg(img)%targetF(3,NAts(img)))
        ALLOCATE(TrainImg(img)%predF(3,NAts(img)))
        ALLOCATE(TrainImg(img)%Free(3,NAts(img)))
      END DO

      !Set unique elemnt array
      DO j=1, nelement
        IF (PRESENT(uniqueNrs)) THEN
          uniqElems(j)=elementArray(uniqueNrs(j))
        ELSE IF (PRESENT(in_uniqElems)) THEN
          uniqElems(j)=in_uniqElems(j)
        ELSE
          PRINT *, 'Information on unique elements is missing. Job is killed by PyAMFF_init'
          STOP
        END IF
      END DO
      
      !Set atomic information of initial images
      DO img=1, nimg
        TrainImg(img)%natoms=NAts(img)
        IF (PRESENT(atomicNrs) .and. PRESENT(uniqueNrs)) THEN 
          DO i=1, nAts(img)
            DO j=1, nelement
              IF (atomicNrs(i,img) == uniqueNrs(j)) THEN 
                TrainImg(img)%symbols(i)=j
                TrainImg(img)%natoms_arr(j)=TrainImg(img)%natoms_arr(j)+1
              END IF
            END DO
          END DO
        ELSE IF (PRESENT(in_natarr)) THEN  !VASP-friendly
          ptr2=0
          DO j=1, nelement
            TrainImg(img)%natoms_arr(j)=in_natarr(j,img)
            TrainImg(img)%symbols(ptr2+1:ptr2+in_natarr(j,img))=j
            ptr2=ptr2+in_natarr(j,img)
          END DO 
          !print *, 'symbols=', TrainImg(img)%symbols 
        ELSE 
          PRINT *, 'Information on atomic numbers is missing. Job is killed by PyAMFF_init'
          STOP
        END IF
        
        !Allocate and set atom_idx 
        ALLOCATE(TrainImg(img)%atom_idx(maxval(TrainImg(img)%natoms_arr),nelement))
        idx_arr=1
        DO i=1, nAts(img)
          ptr=TrainImg(img)%symbols(i)
          TrainImg(img)%atom_idx(idx_arr(ptr),ptr)=i
          idx_arr(ptr)=idx_arr(ptr)+1
        END DO
      END DO
     
      !Store values in main memory
      !total number of elements
      nelements=nelement
      !total number of initial images
      nimages=nimg
      !total natoms of all initial images
      nAtimg=SUM(nAts(1:nimg)) 
      
      !TODO:
      !This should be called maybe in step but only once. 
      !It depends on if we pass the INCAR info in init or step
      !Initiate training (reading mlff, allocate arrays related to backward prop)
      IF (PRESENT(filename)) THEN
        CALL train_init(MAXVAL(nAts),nelement, max_fps,uniqElems,filename,seedval)
      ELSE  
        CALL train_init(MAXVAL(nAts),nelement, max_fps,uniqElems)
      END IF
      !print *, 'train_init is done in pyamff'

      !Allocate input_fps and input_dfps of TrainImg 
      DO img=1, nimg
        max_natarr=MAXVAL(TrainImg(img)%natoms_arr)
        ALLOCATE(TrainImg(img)%input_fps(max_natarr,MAXVAL(nGs),nelement))
        ALLOCATE(TrainImg(img)%input_dfps(NAts(img),NAts(img),3,MAXVAL(nGs)))
      END DO
      !Allocate elementwise magnitude and intercept scale
      ALLOCATE(elem_mag(MAXVAL(nGs),nelement))
      ALLOCATE(elem_intercept(MAXVAL(nGs),nelement))
      elem_mag=0.d0
      elem_intercept=0.d0

      !print *, 'pyamff_init is done in pyamff '

    END SUBROUTINE

    SUBROUTINE PyAMFF_step(posions,trueEs,trueFs,boxes,nelement,in_nimg,totnimg,&
    max_epoch,max_natoms,opt_type,force_coeff,newImg,eflag,fflag,Frees,&
    !energy_tol,force_tol,&
    conv_method,energy_tol,force_tol,grad_tol,&
    new_nat,in_natarr)
    !TODO: Adding multiple images to the training dataset
    !---------------------------------------------------------------------!
    !This subroutine manage one time of whole training process:           !
    !1. Calculate neighborlists and fingerprints of all images            !
    !2. Calculate energy and forces of all images                         !
    !3. Calculate loss gradients and update weights and biases            !
    !4. Go back to step 2 and repete until minimization                   !
    !---------------------------------------------------------------------! 
      IMPLICIT NONE
      !Inputs
      INTEGER :: nelement, in_nimg, totnimg, max_epoch, max_natoms
      CHARACTER(*) :: opt_type, conv_method
      DOUBLE PRECISION :: force_coeff, energy_tol, force_tol, grad_tol
      LOGICAL :: newImg, eflag, fflag
      DOUBLE PRECISION, DIMENSION(in_nimg) :: trueEs
      DOUBLE PRECISION, DIMENSION(3,max_natoms,in_nimg) :: trueFs, Frees
      DOUBLE PRECISION, DIMENSION(3,max_natoms,in_nimg) :: posions
      DOUBLE PRECISION, DIMENSION(9,in_nimg) :: boxes
      !Optional inputs for adding one or multiple new images
      INTEGER, DIMENSION(in_nimg), OPTIONAL :: new_nat
      INTEGER, DIMENSION(nelement,in_nimg), INTENT(IN), OPTIONAL :: in_natarr
      !INTEGER, OPTIONAL :: new_nat
      !INTEGER, DIMENSION(nelement), INTENT(IN), OPTIONAL :: in_natarr
      !INTEGER, DIMENSION(new_nat), INTENT(IN), OPTIONAL :: atomicNrs
      !Variables
      INTEGER :: img, update_idx, i, img_ptr
    

      !flags for training !need better way
      energy_training=eflag
      force_training=fflag

      !Get the initial index of new images for fingerprints
      !This will be used to avoid repeat of computations
      IF (newImg .EQV. .TRUE.) THEN
        !print *, 'new image is added' 
        update_idx=nimages+1
        !print *, 'update_idx=', update_idx
      ELSE !initial
        update_idx=1 
      END IF

      !Update TrainImg for new images
      !print *, 'total nimages=', totnimg
      !TODO: we should be able to add multiple new images to the training set
      img_ptr=1
      DO img=update_idx, totnimg
        IF (newImg .EQV. .TRUE.) THEN
          !Allocate TrainImage key arrays for a new image (similar to pyamff_init but for only new images)
          IF (PRESENT(in_natarr)) THEN !VTST-ML inputs
            CALL UpdateTrainImage(img, new_nat(img_ptr), nelement, in_natarr(:,img_ptr))
            !print *, 'new image training image is updated'
            !print *, 'nimages =', nimages
          ELSE 
            PRINT *, 'new image atomic information is missing. Job is killed by pyamff_step'
            STOP
          END IF
          
          !Reset arrays for backward propagation  
          IF (img == totnimg) CALL init_backward
          
          !Set cells of image img
          TrainImg(img)%cell(1,1)=boxes(1,img_ptr)
          TrainImg(img)%cell(1,2)=boxes(2,img_ptr)
          TrainImg(img)%cell(1,3)=boxes(3,img_ptr)
          TrainImg(img)%cell(2,1)=boxes(4,img_ptr)
          TrainImg(img)%cell(2,2)=boxes(5,img_ptr)
          TrainImg(img)%cell(2,3)=boxes(6,img_ptr)
          TrainImg(img)%cell(3,1)=boxes(7,img_ptr)
          TrainImg(img)%cell(3,2)=boxes(8,img_ptr)
          TrainImg(img)%cell(3,3)=boxes(9,img_ptr)  
          !Set positions  
          DO i=1, TrainImg(img)%natoms
            TrainImg(img)%poscar(i,1:3)=posions(1:3,i,img_ptr)
          END DO
          !Set target E and F
          TrainImg(img)%targetE=trueEs(img_ptr)
          TrainImg(img)%targetF(1:3,1:TrainImg(img)%natoms)=trueFs(1:3,1:TrainImg(1)%natoms,img_ptr)
          TrainImg(img)%Free(1:3,1:TrainImg(img)%natoms)=Frees(1:3,1:TrainImg(1)%natoms,img_ptr)
          img_ptr=img_ptr+1
          !print *, 'new image cell, poscar, true E and F are updated'
        
        !Initial images
        ELSE !Initial image(s)
          !Set cells of image img
          TrainImg(img)%cell(1,1)=boxes(1,img)
          TrainImg(img)%cell(1,2)=boxes(2,img)
          TrainImg(img)%cell(1,3)=boxes(3,img)
          TrainImg(img)%cell(2,1)=boxes(4,img)
          TrainImg(img)%cell(2,2)=boxes(5,img)
          TrainImg(img)%cell(2,3)=boxes(6,img)
          TrainImg(img)%cell(3,1)=boxes(7,img)
          TrainImg(img)%cell(3,2)=boxes(8,img)
          TrainImg(img)%cell(3,3)=boxes(9,img)
          !Set positions  
          DO i=1, TrainImg(img)%natoms
            TrainImg(img)%poscar(i,1:3)=posions(1:3,i,img)
          END DO
          !Set target E and F
          TrainImg(img)%targetE=trueEs(img)
          TrainImg(img)%targetF(1:3,1:TrainImg(img)%natoms)=trueFs(1:3,1:TrainImg(img)%natoms,img)
          TrainImg(img)%Free(1:3,1:TrainImg(img)%natoms)=Frees(1:3,1:TrainImg(img)%natoms,img)
          !print *, 'initial image ', img, 'box, poscar, true E and F are updated'
          !print *, 'initial nimages are ', nimages
        END IF
       
      END DO
      
      !Execute fingerprints, forward, and trainer.  
      !Fingerprints and neighborlists should be skipped if they are already calculated
      !update_idx will be used as an initial pointer of neighborlists, fingerprints calculation
      !forward should be redone with the updated parameters
      !backward only occurs when img=nimg
      DO img=1, totnimg
        img_idx=img
        CALL trainExec(TrainImg(img)%natoms,TrainImg(img)%poscar,&
        TrainImg(img)%cell,TrainImg(img)%symbols,nelement,uniqElems,&
        MAXVAL(TrainImg(img)%natoms_arr),MAXVAL(nGs),opt_type,&
        !max_epoch,force_coeff,energy_tol,force_tol,newImg,update_idx)
        max_epoch,force_coeff,conv_method,energy_tol,force_tol,grad_tol,&
        newImg,update_idx)
      END DO

    END SUBROUTINE

    SUBROUTINE UpdateTrainImage(img, NAt, nelement, in_natarr)
    !----------------------------------------------------------------!
    !Similar to pyamff_init but only for a new image                 !
    !----------------------------------------------------------------!
      IMPLICIT NONE
      !Inputs
      INTEGER, INTENT(IN) :: img, NAt, nelement
      INTEGER, DIMENSION(nelement), INTENT(IN) :: in_natarr
      !Variables
      INTEGER, PARAMETER :: max_fps=500
      INTEGER, PARAMETER :: max_neighs=100 
      INTEGER :: curr_len, i, j, max_natarr, ptr, ptr2
      INTEGER, DIMENSION(nelement) :: idx_arr
      CHARACTER*3, DIMENSION(92) :: elementArray

      DATA elementArray / "H","He","Li","Be","B","C","N","O", &
                 "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc", &
                 "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se", &
                 "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag", &
                 "Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd", &
                 "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta", &
                 "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", &
                 "Fr","Ra","Ac","Th","Pa","U" /  
       
      !Allocate arrays for a new image
      curr_len=SIZE(TrainImg(:))
      !print *, 'current length of TrainImg=', curr_len
      IF (img > curr_len) THEN !need allocation
        !print *, 'TrainImg length is not enough. the length should be doubled'
        CALL Reallocate_TrainImg(curr_len) !reallocate TrainImg with the two times of the current length
        !print *, 'TrainImg is backed up with double length:', SIZE(TrainImg(:))
      END IF  
      ALLOCATE(TrainImg(img)%natoms_arr(nelement))
      !print *, 'natoms_arr allocated'
      ALLOCATE(TrainImg(img)%symbols(NAt))
      !print *, 'symbols allocated'
      ALLOCATE(TrainImg(img)%poscar(NAt,3))
      !print *, 'poscar allocated'
      ALLOCATE(TrainImg(img)%cell(3,3))
      !print *, 'cell allocated'
      ALLOCATE(TrainImg(img)%nneighbors(NAt))
      !print *, 'nneighbors allocated'
      ALLOCATE(TrainImg(img)%neighborlists(NAt,max_neighs))
      !print *, 'neighborlists allocated'
      ALLOCATE(TrainImg(img)%targetF(3,NAt))
      !print *, 'targetF allocated'
      ALLOCATE(TrainImg(img)%predF(3,NAt))
      !print *, 'predF allocated'
      ALLOCATE(TrainImg(img)%Free(3,NAt))
      !Set atomic information of a new image
      TrainImg(img)%natoms=NAt
      ptr2=0
      DO j=1, nelement
        TrainImg(img)%natoms_arr(j)=in_natarr(j)
        TrainImg(img)%symbols(ptr2+1:ptr2+in_natarr(j))=j
        ptr2=ptr2+in_natarr(j)
      END DO
      !print *, 'TrainImg natoms_arr, symbols are updated'

      !Allocate and set atom_idx
      ALLOCATE(TrainImg(img)%atom_idx(maxval(TrainImg(img)%natoms_arr),nelement))
      !print *, 'atom_idx allocated'
      idx_arr=1
      DO i=1, NAt
        ptr=TrainImg(img)%symbols(i)
        TrainImg(img)%atom_idx(idx_arr(ptr),ptr)=i
        idx_arr(ptr)=idx_arr(ptr)+1
      END DO
      !print *, 'atom_idx is updated' 
      !Update values in main memory
      !total number of elements
      nelements=nelement
      !total number of current images
      nimages=nimages+1
      !total natoms of all current images
      nAtimg=nAtimg+NAt !This should work
      !print *, 'nimages=', nimages 
      !Allocate input_fps and input_dfps of new image
      max_natarr=MAXVAL(TrainImg(img)%natoms_arr)
      ALLOCATE(TrainImg(img)%input_fps(max_natarr,MAXVAL(nGs),nelement))
      ALLOCATE(TrainImg(img)%input_dfps(NAt,NAt,3,MAXVAL(nGs)))
      
    END SUBROUTINE

    SUBROUTINE Reallocate_TrainImg(curr_nimg)
    !---------------------------------------------------------------------!
    ! Double the length of TrainImg.                                      !
    !---------------------------------------------------------------------!
      IMPLICIT NONE
      !Inputs
      INTEGER, INTENT(IN) :: curr_nimg
      !Variables
      INTEGER, PARAMETER :: max_fps=500
      INTEGER, PARAMETER :: max_neighs=100
      INTEGER :: img, nelement, nat, max_natarr
      TYPE(ImgInfo), DIMENSION(:), ALLOCATABLE :: bk_TrainImg
      
      !Allocate for backups and the doubled ones
      ALLOCATE(bk_TrainImg(curr_nimg))
      nelement=SIZE(TrainImg(img)%natoms_arr)
      DO img=1, curr_nimg
        ALLOCATE(bk_TrainImg(img)%natoms_arr(nelement))
        nat=TrainImg(img)%natoms
        ALLOCATE(bk_TrainImg(img)%symbols(nat))
        ALLOCATE(bk_TrainImg(img)%poscar(nat,3))
        ALLOCATE(bk_TrainImg(img)%cell(3,3))
        ALLOCATE(bk_TrainImg(img)%nneighbors(nat))
        ALLOCATE(bk_TrainImg(img)%neighborlists(nat,max_neighs))
        ALLOCATE(bk_TrainImg(img)%targetF(3,nat))
        ALLOCATE(bk_TrainImg(img)%predF(3,nat))
        max_natarr=MAXVAL(TrainImg(img)%natoms_arr)
        ALLOCATE(bk_TrainImg(img)%input_fps(max_natarr,MAXVAL(nGs),nelement)) !nGs must be known
        ALLOCATE(bk_TrainImg(img)%input_dfps(nat,nat,3,MAXVAL(nGs)))
      END DO

      !Copy all of the data of the current TrainImg
      bk_TrainImg=TrainImg

      !Deallocate TrainImg 
      DEALLOCATE(TrainImg)

      !Double the length of TrainImg
      ALLOCATE(TrainImg(2*curr_nimg))
      !And reallocate the doubled TrainImg 
      DO img=1, curr_nimg
        ALLOCATE(TrainImg(img)%natoms_arr(nelement))
        nat=TrainImg(img)%natoms
        ALLOCATE(TrainImg(img)%symbols(nat))
        ALLOCATE(TrainImg(img)%poscar(nat,3))
        ALLOCATE(TrainImg(img)%cell(3,3))
        ALLOCATE(TrainImg(img)%nneighbors(nat))
        ALLOCATE(TrainImg(img)%neighborlists(nat,max_neighs))
        ALLOCATE(TrainImg(img)%targetF(3,nat))
        ALLOCATE(TrainImg(img)%predF(3,nat))
        max_natarr=MAXVAL(bk_TrainImg(img)%natoms_arr)
        ALLOCATE(TrainImg(img)%input_fps(max_natarr,MAXVAL(nGs),nelement))
        ALLOCATE(TrainImg(img)%input_dfps(nat,nat,3,MAXVAL(nGs)))
      END DO

      !Copy backup data to the doubled one 
      DO img=1, curr_nimg
        TrainImg(img)=bk_TrainImg(img)
      END DO

      !Deallocate backup 
      DEALLOCATE(bk_TrainImg)

    END SUBROUTINE

    SUBROUTINE PyAMFF_calc(nAtoms,posion,box,in_natarr,max_natarr,nelement,uniq_elements,&
    predE,predF)
    !----------------------------------------------------------------------!
    ! Calculate energy and forces on the input model                       !
    !----------------------------------------------------------------------!
      use fpCalc
      IMPLICIT NONE
      !Inputs
      INTEGER :: nAtoms, nelement,max_natarr
      INTEGER, DIMENSION(nelement) :: in_natarr
      DOUBLE PRECISION, DIMENSION(3,nAtoms) :: posion
      DOUBLE PRECISION, DIMENSION(9) :: box
      CHARACTER*2, DIMENSION(nelement) :: uniq_elements
      !Outputs
      DOUBLE PRECISION :: predE, predF(3,nAtoms)
      !Variables
      INTEGER, PARAMETER :: max_neighs=100
      INTEGER :: i, j, ptr, ptr2
      INTEGER, PARAMETER :: forceEngine = 1
      INTEGER, DIMENSION(nelement) :: idx_arr
      INTEGER, DIMENSION(nAtoms, max_neighs) :: neighs
      INTEGER, DIMENSION(nAtoms, nAtoms) :: sub_neighs
      INTEGER, DIMENSION(nAtoms) :: num_neigh, symbols, sub_num_neigh
      INTEGER, DIMENSION(max_natarr,nelement) :: in_atomidx
      DOUBLE PRECISION, DIMENSION(nAtoms,3) :: pos_car
      DOUBLE PRECISION, DIMENSION(3,3) :: cell
      !Note nGs should be known to use this subroutine
      DOUBLE PRECISION, DIMENSION(nAtoms, MAXVAL(nGs)) :: fps
      DOUBLE PRECISION, DIMENSION(max_natarr, MAXVAL(nGs), nelement) :: ordered_fps
      DOUBLE PRECISION, DIMENSION(nAtoms, max_neighs, 3, MAXVAL(nGs)) :: temp_dfps
      DOUBLE PRECISION, DIMENSION(nAtoms, nAtoms, 3, MAXVAL(nGs)) :: dfps

      !print *, 'vasp calculator is called'
     
      ptr2=0
      DO i=1, nelement
        symbols(ptr2+1:ptr2+in_natarr(i))=i
        ptr2=ptr2+in_natarr(i)
      END DO
      idx_arr=1
      DO i=1, nAtoms
        ptr=symbols(i)
        in_atomidx(idx_arr(ptr),ptr)=i
        idx_arr(ptr)=idx_arr(ptr)+1
      END DO
      !Set pos_car and cell
      DO i=1, nAtoms
        pos_car(i,1:3)=posion(1:3,i)
      END DO
      cell(1,1:3)=box(1:3)
      cell(2,1:3)=box(4:6)
      cell(3,1:3)=box(7:9)
      
      !Set training flag to be false  
      energy_training=.FALSE.
      force_training=.FALSE.
      
      fps=0.
      dfps=0.
      temp_dfps=0.

      CALL calcfps(nAtoms, pos_car, cell, symbols, MAXVAL(nGs), nelement, forceEngine, &
                   fps, temp_dfps, neighs, num_neigh)
      !print *, 'fps=', fps
      !print *, 'fp is calculated'
      CALL ghost_dfps_correct(nelement, nAtoms, MAXVAL(nGs), MAX_NEIGHS, num_neigh, neighs, &
           sub_num_neigh, sub_neighs, temp_dfps, dfps)

      !Cleanup ghost parts
      CALL atomsCleanup
      
      CALL update_atomInfo2(nAtoms, nelement, MAXVAL(nGs), in_natarr)
      !print *, 'atomic info is updated'

      !Reorder fingerprints
      DO i = 1, nelements
        DO j = 1, natoms_arr(i)
          ordered_fps(j,1:nGs(i),i) = fps(in_atomidx(j,i),1:nGs(i))
        END DO
      END DO
      CALL normalizeFPs2(nelement, nAtoms, uniq_elements, MAXVAL(nGs), MAXVAL(sub_num_neigh), &
      max_natarr,sub_num_neigh,sub_neighs,symbols,ordered_fps, dfps(:,1:MAXVAL(sub_num_neigh)+1,:,:))

      CALL forward(sub_num_neigh, MAXVAL(sub_num_neigh),sub_neighs(:,1:MAXVAL(sub_num_neigh)),&
      symbols, in_atomidx, ordered_fps, dfps(:,1:MAXVAL(sub_num_neigh)+1,:,:), &
      MAXVAL(nGs), max_natarr, MAXVAL(nhidneurons))
      !print *, 'forward propagation is done'
      !print *, 'Etotal=', Etotal
      !print *, 'forces=', forces
      predE=Etotal
      predF=forces
      CALL nncleanup_atom
     
    END SUBROUTINE

    SUBROUTINE PyAMFF_clean(opt_type,opt_flag)
    !Deallocate all of allocated arrays
      IMPLICIT NONE
      CHARACTER(*) :: opt_type
      LOGICAL :: opt_flag

      IF (opt_flag) THEN
        CALL traincleanup(opt_type)
      ELSE
        CALL backcleanup
      END IF  
    END SUBROUTINE

END MODULE
