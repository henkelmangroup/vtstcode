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
    
    INTEGER :: nions, nelement, iu0, iu6, totnimg
    INTEGER :: images_local
    INTEGER, DIMENSION(:), ALLOCATABLE :: NAts
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: nat_arr !2nd dim is for consistency with PyAMFF structure
    LOGICAL :: initImg, newImg, eflag, fflag, global
    INTEGER :: maxEpochs
    CHARACTER*8 :: optType
    CHARACTER*2, DIMENSION(:), ALLOCATABLE :: uniqEls
    REAL(q) :: force_coeff, force_tol, energy_tol, gradtol
    !Stopping criteria test variables
    !TEST #1
    !REAL(q) :: init_maxforce, old_maxforce
    !REAL(q), ALLOCATABLE :: fmaxatom(:), ftot(:), dmaxatom(:) 
    !TEST #2
    !REAL(q) :: trustR, ml_minmove
    !REAL(q), ALLOCATABLE :: init_R(:,:,:), change_in_R(:,:,:)
    !REAL(q) :: EDIFFG_local

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

      INTEGER :: i, in_maxnimg
      CHARACTER*20 :: filename
 
      nions=T_INFO%nions
      nelement=T_INFO%ntyp
      iu0=IO%IU0
      iu6=IO%IU6
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

      ! This is for gradnorm as a model convergence criteria.
      gradtol=0.001_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_TOL','=','#',';','F', &
                  IDUM,gradtol,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_TOL'' from file INCAR.'
         STOP
      ENDIF 
 
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
      uniqEls(1:nelement), nat_arr(1:nelement,1:images_local),in_maxnimg,filename)
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
      
      !Perform PyAMFF training
      IF (newImg) THEN
        CALL PyAMFF_step(posion_all,toten_all,force_all,boxes_all,nelement,&
          images_local,totnimg,maxEpochs,MAXVAL(nat_arr),optType,force_coeff,newImg,&
          eflag,fflag,Free_all,energy_tol,force_tol,&
          !eflag,fflag,Free_all,gradtol,&   !This is for gradnorm as model convergence criteria
          NAts,nat_arr)
      ELSE
        CALL PyAMFF_step(posion_all,toten_all,force_all,boxes_all,nelement,&
          images_local,totnimg,maxEpochs,MAXVAL(nat_arr),optType,force_coeff,newImg,&
          eflag,fflag,Free_all,energy_tol,force_tol)
          !eflag,fflag,Free_all,gradtol)    !This is for gradnorm as model convergence criteria
      END IF

      IF (iu6>=0) WRITE(iu6,'(A11,A)') 'ML-PyAMFF:',' Training information '
      IF (iu6>=0) WRITE(iu6,'(A)') 'epoch    lossValue   EnergyRMSE    ForceRMSE'
      IF (iu6>=0) WRITE(iu6,'(I6,F16.8,F16.8,F16.8)') final_epoch, final_loss, final_eRMSE, final_fRMSE

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

    SUBROUTINE ml_pes(posion,latt_a,latt_b,Free,predF)
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

  END MODULE
