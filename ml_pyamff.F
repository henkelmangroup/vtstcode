#include "symbol.inc"
!**********************************************************************
!
! This module implements the PyAMFF ML method.  For more information refer
! to the web page:
!   http://theory.cm.utexas.edu/vtsttools/
! 
!
!**********************************************************************

  MODULE ML_PyAMFF
    USE prec
    USE lattice
    USE PyAMFF !pyamff source code 
    USE qm !to update posions
    !USE fire
    USE main_mpi
    
    IMPLICIT NONE
    SAVE
    PRIVATE
    PUBLIC :: ML_PyAMFF_train, ML_PyAMFF_step, ML_PyAMFF_init, ML_PyAMFF_cleanup, ml_pes
    
    INTEGER :: nions, nelement, iu0, iu6, mdout, pesout, totnimg
    INTEGER :: images_local
    INTEGER, DIMENSION(:), ALLOCATABLE :: NAts
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: nat_arr !2nd dim is for consistency with PyAMFF structure
    LOGICAL :: initImg, newImg, eflag, fflag, global
    INTEGER :: maxEpochs
    CHARACTER*8 :: optType, conv_method
    CHARACTER*2, DIMENSION(:), ALLOCATABLE :: uniqEls
    REAL(q) :: force_coeff, force_tol, energy_tol, gradtol
    
    !for qm step calc
    REAL(q) :: maxmove, dt, stepsize

    !NEB
    LOGICAL :: cell_flag
    INTEGER :: max_iter
    REAL(q), ALLOCATABLE :: posion_all(:,:,:),latt_a_all(:,:,:),latt_b_all(:,:,:)
    REAL(q), ALLOCATABLE :: boxes_all(:,:), toten_all(:), force_all(:,:,:)
    REAL(q), ALLOCATABLE :: Free_all(:,:,:) 

!**********************************************************************
! Machine learning (ML) method
!**********************************************************************
    CONTAINS

    SUBROUTINE ML_PyAMFF_init(T_INFO,IO) !posion info, MLvariable info
      USE base
      USE poscar
      !USE main_mpi !NEB
      IMPLICIT NONE
      !Inputs
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO
      
      !Variables
      INTEGER IDUM,IERR,N
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      INTEGER :: i, in_maxnimg, seedval
      CHARACTER*20 :: filename
 
      nions=T_INFO%nions
      nelement=T_INFO%ntyp
      iu0=IO%IU0
      iu6=IO%IU6
      !mdout=IO%MDOUT ! for debug, eventually it will be gone
      !pesout=IO%PESOUT ! for debug, eventually it will be gone
      initImg=.TRUE.

      IF (images == 0) THEN 
        images_local=1
      ELSE
        images_local=images
      END IF  
      totnimg=images_local

      !Allocate common variables
      ALLOCATE(NAts(images_local))
      ALLOCATE(uniqEls(nelement)) 
      ALLOCATE(nat_arr(nelement,images_local))
      dt=0.1_q
      maxmove=0.01_q
     
      !Allocate posions, latts, boxes, energy, force, Frees
      ALLOCATE(posion_all(3,nions,images_local))
      ALLOCATE(latt_a_all(3,3,images_local))
      ALLOCATE(latt_b_all(3,3,images_local))
      ALLOCATE(boxes_all(9,images_local))
      ALLOCATE(toten_all(images_local))
      ALLOCATE(force_all(3,nions,images_local))
      ALLOCATE(Free_all(3,nions,images_local)) 
 
      ! Check cell_flag in NEB
      cell_flag=.FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LNEBCELL','=','#',';','L', &
       &            IDUM,RDUM,CDUM,cell_flag,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LNEBCELL'' from file INCAR.'
        STOP
      ENDIF
        
      !Check global flag 
      global=.TRUE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LGLOBAL','=','#',';','L', &
     &            IDUM,RDUM,CDUM,global,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
        IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''LGLOBAL'' from file INCAR.'
        STOP
      ENDIF

      !read variables used for the ML optimzier
      !Note: RDATAB subroutine is in /lib/drdatab.F
     
      ! Check uncertainty method 
      conv_method='GRADNORM'
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_CONV','=','#',';','S', &
      &            IDUM,RDUM,CDUM,LDUM,conv_method,N,LEN(filename),IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
          IF (iu0>=0) WRITE(iu0,*)'Error reading item ''PYAMFF_CONV'' from file INCAR.'
          STOP
      ENDIF
      
      IF (conv_method == 'RMSE') THEN  
          energy_tol=0.001_q
          CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_ETOL','=','#',';','F', &
                  IDUM,energy_tol,CDUM,LDUM,CHARAC,N,1,IERR)
          IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
              IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_ETOL'' from file INCAR.'
              STOP
          ENDIF
             
          force_tol=0.010_q
          CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_FTOL','=','#',';','F', &
                  IDUM,force_tol,CDUM,LDUM,CHARAC,N,1,IERR)
          IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
             IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_FTOL'' from file INCAR.'
             STOP
          ENDIF
      ELSE IF (conv_method == 'GRADNORM') THEN 
          ! This is for gradnorm as a model convergence criteria.
          gradtol=0.001_q
          CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_TOL','=','#',';','F', &
                  IDUM,gradtol,CDUM,LDUM,CHARAC,N,1,IERR)
          IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
              IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_TOL'' from file INCAR.'
              STOP
          ENDIF
      ELSE 
          IF (iu0>=0) WRITE(iu0,*) 'Value Error of ''PYAMFF_CONV'' from file INCAR.'
          IF (iu0>=0) WRITE(iu0,*) 'Available values: ''GRADNORM, RMSE'
      END IF

      force_coeff=1.00_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_FCOEFF','=','#',';','F', &
                  IDUM,force_coeff,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_FCOEFF'' from file INCAR.'
         STOP
      ENDIF

      maxEpochs=2000
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_MAXEPOCH','=','#',';','I', &
                  maxEpochs,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_MAXEPOCH'' from file INCAR.'
         STOP
      ENDIF

      filename='mlff.pyamff'
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_MODEL','=','#',';','S', &
      &            IDUM,RDUM,CDUM,LDUM,filename,N,LEN(filename),IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
          IF (iu0>=0) WRITE(iu0,*)'Error reading item ''PYAMFF_MODEL'' from file INCAR.'
          STOP
      ENDIF

      optType='RPROP'
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_OPT','=','#',';','S', &
                  IDUM,RDUM,CDUM,LDUM,optType,N,LEN(optType),IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
        IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_OPT'' from file INCAR.'
        STOP
      ENDIF

      seedval=123456789
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_SEED','=','#',';','I', &
                  seedval,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_SEED'' from file INCAR.'
         STOP
      ENDIF

      !Assign number of atoms to the array type
      NAts(1:images_local)=nions
      !Prepare inputs of PyAMFF_init
      DO i=1, nelement
        uniqEls(i)=T_INFO%type(i)
        nat_arr(i,1:images_local)=T_INFO%nityp(i)
      END DO
      
      !Initiate pyamff in pyamff source code (/path-to-pyamff/pyamff/fortran/model.f90) 
      !All of the required variables will be allocated in pyamff side
      !TEMP: set maxnimg (memory issue on frontera)
      in_maxnimg=150
      CALL PyAMFF_init(NAts(1:images_local), nelement, images_local, & 
      uniqEls(1:nelement), nat_arr(1:nelement,1:images_local),in_maxnimg,filename,seedval)
      !print *, 'pyamff_init is done in vasp side'

    END SUBROUTINE    
   
    SUBROUTINE ML_PyAMFF_train(posion,toten,force,latt_a,latt_b,Free,&
               in_eflag,in_fflag)

      IMPLICIT NONE
      ! Inputs
      REAL(q) :: posion(3,nions), toten, force(3,nions), Free(3,nions)
      REAL(q) :: latt_a(3,3), latt_b(3,3)
      ! Optional Inputs
      LOGICAL, OPTIONAL :: in_eflag, in_fflag
      ! Variables
      INTEGER :: node, nimg, ni, nj
      REAL(q) :: ftemp
      DOUBLE PRECISION, DIMENSION(1) :: trueEs
      DOUBLE PRECISION, DIMENSION(9,1) :: boxes
      DOUBLE PRECISION, DIMENSION(3,nions,1) :: posions, trueFs, Frees
      ! Calc_pes
      REAL(q) :: pes_posion(3,nions)

      !Set flags whether it is initial images or not and update total number of images 
      IF (Initimg .EQV. .TRUE.) THEN
        Initimg=.FALSE. !For the next time to call step function
        newImg=.FALSE.

        IF (PRESENT(in_eflag)) THEN
          eflag=in_eflag
        ELSE
          eflag=.TRUE.
        END IF

        IF (PRESENT(in_fflag)) THEN
          fflag=in_fflag
        ELSE
          fflag=.TRUE.
        END IF
      ELSE !not initial images
        newImg=.TRUE.
        totnimg=totnimg+images_local
      END IF

      !Prepare inputs of PyAMFF_step
      !Convert the position into Cartesian coordinates
      CALL dirkar(nions, posion, latt_a)
      IF (images>0) THEN
#if defined(MPI) || defined (MPI_CHAIN)
        node=comm_chain%node_me
        !Store all NEB images' energy, force, positions
        posion_all(:,:,1:images)=0._q
        posion_all(:,:,node)=posion
        CALLMPI_C( M_SUM_d( comm_chain, posion_all(1,1,1), nions*3*images))

        IF (cell_flag) THEN
          latt_a_all(:,:,1:images)=0._q
          latt_a_all(:,:,node)=latt_a
          latt_b_all(:,:,1:images)=0._q
          latt_b_all(:,:,node)=latt_b
          CALLMPI_C( M_SUM_d( comm_chain, latt_a_all(1,1,1), 3*3*images))
          CALLMPI_C( M_SUM_d( comm_chain, latt_b_all(1,1,1), 3*3*images))
        END IF

        boxes_all(:,1:images)=0._q
        boxes_all(1:3,node)=latt_a(1:3,1)
        boxes_all(4:6,node)=latt_a(1:3,2)
        boxes_all(7:9,node)=latt_a(1:3,3)
        CALLMPI_C( M_SUM_d( comm_chain, boxes_all(1,1), 9*images))

        toten_all=0._q
        toten_all(node)=toten
        CALLMPI_C( M_SUM_d (comm_chain,toten_all(1),images))

        force_all=0._q
        force_all(1:3,1:nions,node)=force(1:3,1:nions)
        CALLMPI_C( M_SUM_d (comm_chain,force_all(1,1,1),nions*3*images))      
        Free_all=0._q
        Free_all(1:3,1:nions,node)=Free(1:3,1:nions)
        CALLMPI_C( M_SUM_d (comm_chain,Free_all(1,1,1),nions*3*images))
#endif
      ELSE
        posion_all(:,:,:)=0._q
        posion_all(:,:,1)=posion

        boxes_all(:,:)=0._q
        boxes_all(1:3,1)=latt_a(1:3,1)
        boxes_all(4:6,1)=latt_a(1:3,2)
        boxes_all(7:9,1)=latt_a(1:3,3)

        toten_all=0._q
        toten_all(1)=toten
        
        force_all(:,:,:)=0._q
        force_all(:,:,1)=force

        Free_all(:,:,:)=0._q
        Free_all(:,:,1)=Free 
      END IF
     
      ! TODO: Pass conv_method  
      !Perform PyAMFF training
      IF (newImg) THEN
        CALL PyAMFF_step(posion_all,toten_all,force_all,boxes_all,nelement,&
          images_local,totnimg,maxEpochs,MAXVAL(nat_arr),optType,force_coeff,newImg,&
           eflag,fflag,Free_all,conv_method,energy_tol,force_tol,gradtol,&
           NAts,nat_arr)
      ELSE
        CALL PyAMFF_step(posion_all,toten_all,force_all,boxes_all,nelement,&
          images_local,totnimg,maxEpochs,MAXVAL(nat_arr),optType,force_coeff,newImg,&
           eflag,fflag,Free_all,conv_method,energy_tol,force_tol,gradtol)
      END IF

      IF (iu6>=0) WRITE(iu6,'(A11,A)') 'ML-PyAMFF:',' Training information '
      IF (iu6>=0) WRITE(iu6,'(A)') '          epoch      lossValue       EnergyRMSE       ForceRMSE'
      IF (iu6>=0) WRITE(iu6,'(A)') '        ---------------------------------------------------------'  
      IF (iu6>=0) WRITE(iu6,'(A7,I6,F16.8,F16.8,F16.8)') 'Final: ',final_epoch, final_loss, final_eRMSE, final_fRMSE

      ! Pd111_H: evolution of PES
      !pes_posion=posion
      !CALL calc_pes(pes_posion,boxes_all(:,:1),1.0_q,21) 
        
      !Convert posion back to direct
      CALL kardir(nions,posion,latt_b)

    END SUBROUTINE 
 
    SUBROUTINE ML_PyAMFF_step(optflag,posion,toten,force,hstress,latt_a,latt_b,Free)
     
      IMPLICIT NONE
      !Inputs
      REAL(q) :: posion(3,nions), toten, force(3,nions), Free(3,nions)
      REAL(q) :: hstress(3,3), latt_a(3,3), latt_b(3,3)
      LOGICAL :: optflag

      ! Take qm-step on ML-PES
      CALL qm_step(optflag,posion,toten,force,hstress,latt_a,latt_b)
      !CALL fire_step(optflag,posion,toten,force,hstress,latt_a,latt_b)
  
    END SUBROUTINE
   
    SUBROUTINE ML_PyAMFF_cleanup
      IMPLICIT NONE
    
      IF (final_epoch /= 0) THEN
        CALL PyAMFF_clean(optType,.TRUE.)
      ELSE
        CALL PyAMFF_clean(optType,.FALSE.)
      END IF

    END SUBROUTINE

    SUBROUTINE ml_pes(posion,latt_a,latt_b,Free,predE,predF)
      IMPLICIT NONE
      !Inputs
      REAL(q) :: posion(3,nions), latt_a(3,3), latt_b(3,3)
      REAL(q) :: Free(3,nions) 
      !Outputs
      REAL(q) :: predE, predF(3,nions)
      !Variables
      DOUBLE PRECISION :: boxes(9,1)

      !current dims of boxes are 2 but it is tentative 
      boxes(1:3,1)=latt_a(1:3,1)
      boxes(4:6,1)=latt_a(1:3,2)
      boxes(7:9,1)=latt_a(1:3,3)

      !Convert the position into Cartesian coordinates
      CALL dirkar(nions, posion, latt_a)
      !Run calculator to get energy and force of the posion
      CALL PyAMFF_calc(nions,posion,boxes(1:9,1),nat_arr(:,1),MAXVAL(nat_arr),&
      nelement,uniqEls,predE,predF)
      ! For frozen atoms
      predF=predF*Free
      ! convert position back to direct coordinates
      CALL kardir(nions,posion,latt_b)

    END SUBROUTINE 

!    SUBROUTINE write_mlff_vasp
!      !------------------------------------------------------------------------!
!      !This is currently temporary format. Only prints out Model parameters    !
!      !Eventually this should go to the fingerprints to write the fingerprint  !
!      !parameters as well                                                      !
!      !------------------------------------------------------------------------!
!      IMPLICIT NONE
!      !Variables
!      INTEGER :: i, l
!
!      !TODO: file name should beTrainImg(img)%natoms!TODO: make more formatted
!      OPEN(mdout, file="trained.mlff", status='unknown') 
!      WRITE(mdout,'(a)') '#Model Parameters'
!      DO i=1, nelement
!        WRITE(mdout,'(A2)') uniqEls(i)
!        WRITE(mdout,*) in_weights(1:nGs(i),1:nhidneurons(1),i), '#inputLayer weight'
!        WRITE(mdout,*) in_biases(1:nhidneurons(1),i), '#inputLayer bias'
!        DO l=1, nhidlayers-1
!          WRITE(mdout,*) hid_weights(1:nhidneurons(l),1:nhidneurons(l+1),l,i),&
!                                 '#hiddenLayer_',l,'weight'
!          WRITE(mdout,*) hid_biases(1:nhidneurons(l+1),l,i),&
!                                 '#hiddenLayer_',l,'bias'
!        END DO
!        WRITE(mdout,*) out_weights(1:nhidneurons(nhidlayers),1,i), '#outputLayer weight'
!        WRITE(mdout,*) out_biases(i), '#outputLayer bias'
!      END DO
!      WRITE(mdout,'(a)') '#Energy Scaling Parameters'
!      WRITE(mdout,'(a)') scaler_type
!      WRITE(mdout,*) slope, intercept
!      !CLOSE(mdout)
!    END SUBROUTINE
!
!    SUBROUTINE calc_pes(pes_posion,boxes,maxz,move_idx)
!    !SUBROUTINE calc_pes(pes_posion,latt_a,maxy,move_idx)
!      IMPLICIT NONE
!      !Inputs
!      REAL(q) :: maxz, pes_posion(3,nions)
!      DOUBLE PRECISION, DIMENSION(9,1) :: boxes
!      INTEGER :: move_idx 
!      !Variables
!      INTEGER :: i
!      REAL(q) :: z, step
!      REAL(q) :: predE, predF(3,nions)
!
!      !print *, 'calc pes called'
!      step=maxz/50.
!      OPEN(pesout, file="pes.dat", status='unknown') 
!      DO i=1, 50
!        !Move atom "move_idx" from 0 to maxy
!        z=18.0+(i-1)*step
!        !z=5.904+(i-1)*step
!        pes_posion(3,move_idx)=z
!        !pes_posion(1,move_idx)=z
!        !Run on trained model
!        CALL PyAMFF_calc(nions,pes_posion,boxes(1:9,1),nat_arr(:,1),MAXVAL(nat_arr),&
!        nelement,uniqEls,predE,predF)
!        IF (i==1) THEN 
!          WRITE(pesout,'(a)') "#step   z   predE " !Pd111_h
!          !WRITE(pesout,'(a)') "#step   x   predE "
!        END IF
!        WRITE(pesout,'(I6,F16.8,F16.8)') i, z, predE
!      END DO
!    END SUBROUTINE

  END MODULE
