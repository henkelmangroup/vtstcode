#include "symbol.inc"
!**********************************************************************
!
! Module which controls the running of four methods, the nudged
! elastic band (neb.F), the dimer (dimer.F), lanczos (lanczos.F), and
! the dynamical matrix method (dynmat.F).  The purpose of this module
! is to determine which of these methods should be run, and to do so at
! each ionic step.
!
! A set of force-based optimers are also included for use with normal
! optimizations, or transition state calculations.  The optimizers are
! steepest-descent, quick-min, conjugate-gradients, and LBFGS.
!
! NOTE: the vasp folks have now added their own dynamical matrix
! code.  The only advantage of this dynamical matrix implementation
! is that you can combine forces from multiple vasp runs.  This
! makes it easier to separate a large calculation into several jobs,
! and systematically check for convergence of the normal modes, or
! prefactors of reactions.
!
! For more information see: http://henkelmanlab.org/vtsttools/
!
! Contributers:
!   Andri Arnaldsson
!   Sam Chill
!   Graeme Henkelman
!   Hannes Jonsson
!   Jiyoung Lee
!   Daniel Sheppard
!   Akksay Singh
!   Blas Uberuaga
!   Eboni Williams
!   Lijun Xu
!   Liang Zhang
!
! Email: henkelman@utexas.edu
!
!**********************************************************************

  MODULE chain
    USE prec
    USE main_mpi
    USE poscar
    USE lattice
    USE neb
    USE dynmat
    USE dimer
    USE lanczos
    USE bbm
    USE instanton
    USE opt
    USE ML_PyAMFF, only: ML_PyAMFF_train, ML_PyAMFF_step, ML_PyAMFF_cleanup, ml_pes
    USE PyAMFF, only: final_epoch

    IMPLICIT NONE
    SAVE
!    PRIVATE
    PUBLIC :: chain_force, chain_init, chain_stress, parallel_tempering
    PUBLIC :: Sum_Chain, And_Chain, LHYPER_NUDGE, ML_chain_criteria
    INTEGER mpmd_client_rank
#ifdef EAM
    PRIVATE :: EAMForce
#endif
    INTEGER :: ICHAIN, IOPT, ISIF_local
    LOGICAL :: optflag, fconverge, ftot_flag, LINTERACT, LMPMD
    LOGICAL :: sconverge, cell_flag, twodim_flag
    REAL(q) :: EDIFFG_local, ftot_val, jacobian, PSTRESS_local
    REAL(q),ALLOCATABLE :: Free(:,:)

#ifdef EAM
    ! Variables used for calling the EAM potential (testing)
    LOGICAL :: eaminit
    REAL(q),ALLOCATABLE :: Ream(:,:), Feam(:,:)
    REAL(q),ALLOCATABLE :: Rvasp(:,:), Fvasp(:,:)
    REAL(q) :: Uvasp, Ueam
#endif

    REAL(q),ALLOCATABLE,SAVE :: posion_all(:,:,:)
    REAL(q) :: spring=10
    INTEGER :: nions_max

!!! Tempering variables

! tag VCAIMAGES allows to perform two MD's with e.g. different POTCAR files
! allowing to do force averaging
    LOGICAL :: LVCAIMAGES
! weight of first images, weight of second image is 1-VCAIMAGES
! average forces over two images
! it is recommended to combine this with the VCA type calculations
! where the POSCAR/POTCAR/INCAR etc. are strictly identical
! and the VCA tag is set for one type to 0 and 1 respectively
! in the two POTCAR files
! forces and energies are averaged over the two nodes
! IMAGES in mpi_main.F is forced to two in this case
    REAL(q) :: VCAIMAGES


! the tag LTEMPER allows to perform parallel tempering runs
! IMAGES images of VASP are kicked up and the 
! temperatures between the images (TEBEG) are swapped using
! a Monte-Carlo algorithm (with corresponding temperature scaling)
    LOGICAL, SAVE :: LTEMPER = .FALSE.
    INTEGER, SAVE :: NTEMPER       !< replica exchange every NTEMPER steps
    !> random variable between [0,NTEMPER*2]
    !> will be decremented and when 0 is reached
    !> replica exchange is attempted
    INTEGER, SAVE :: NTEMPER_NEXT=0 

    INTEGER, ALLOCATABLE, SAVE :: ATTEMPTS(:)
    INTEGER, ALLOCATABLE, SAVE :: SUCCESS(:)
!!!
! ML-PYAMFF
    INTEGER :: ml_chain_iter, max_iter
    REAL(q) :: force_tol, SWFTOL !Switch force tolerance from ML to LBFGS
    REAL(q) :: trustR, minmove
    LOGICAL :: ml_turnoff, ml_chain
!**********************************************************************
!  General force routine for any method using the multiple image mode.
!  The variable ICHAIN determines which method to use (see chain_init)
!**********************************************************************

  CONTAINS
    SUBROUTINE chain_force(nions, posion, toten, force, stress, a, b, IU6)
      INTEGER :: nions, ni, nj, IU6, I
      REAL(q) :: ftot, frms, fmaxatom, ftemp, fmaxdim, toten
      REAL(q),DIMENSION(3,nions) :: posion, force
      REAL(q),DIMENSION(3,nions) :: posion_vasp, force_vasp, force_dimlan
      REAL(q),DIMENSION(3,3) :: a, b
      REAL(q),DIMENSION(3,3) :: stress, hstress, stress_vasp
      REAL(q),DIMENSION(3,3) :: sdA, sdB, sdA2
      REAL(q) :: stot, smaxdim
      REAL(q) :: omega
      LOGICAL :: stopcar_exists, newcar_exists
      INTEGER NIOND, NIONPD, NTYPPD, NTYPD, ierr
      TYPE (latt):: LATT_CUR
      TYPE (TYPE_info) :: T_I
      TYPE (dynamics) :: DYN
      TYPE(nh_chains) :: NHC
      TYPE (in_struct) :: IO

      INTEGER node
      REAL(q),ALLOCATABLE,SAVE :: force_all(:,:,:)

      ! ML-optimizer
      LOGICAL :: ml_cleanup
      REAL(q) :: dft_fmaxatom, ml_fmaxatom, toten_ml
      REAL(q), DIMENSION(3,nions) :: force_mlreturn, force_ml, init_R, posion_ml

#if defined(MPI) || defined (MPI_CHAIN)
      node = comm_chain%node_me

!======================================================================
! Parallel tempering return
!======================================================================
      IF (LTEMPER) THEN
         RETURN
!======================================================================
! VCA average forces over two images
!======================================================================
      ELSE IF (LVCAIMAGES) THEN
         ALLOCATE(force_all(3,nions,2))
         force_all(:,:,1:2) = 0
         force_all(:,:,node) = force
#ifndef old_vca
         ! first core in image merges the data
         IF (comm%node_me==1) THEN
            CALLMPI_C( M_sum_d( comm_chain, force_all(1,1,1), nions*3*2))
            force(:,:) = force_all(:,:,1)*VCAIMAGES + force_all(:,:,2)*(1-VCAIMAGES)
         ELSE
            force=0
         ENDIF
         ! now sum over all cores in one image
         CALLMPI_C( M_sum_d( comm, force(1,1), nions*3))
#else
         CALLMPI_C( M_sum_d( comm_chain, force_all(1,1,1), nions*3*2))
         force=0
         force(:,:) = force_all(:,:,1)*VCAIMAGES + force_all(:,:,2)*(1-VCAIMAGES)
#endif
         DEALLOCATE(force_all)
         RETURN
      ENDIF
#endif
!!!

      ! Interactive Mode only have rank 0 process write data
!      IF (LINTERACT .AND. IU6>0) THEN
      IF (LINTERACT) THEN
        IF (IU6>0) THEN
          ! Write the force-energy file.
          OPEN(UNIT = 1, FILE = "FU")
          WRITE(1, *) toten
          DO ni = 1,nions
            WRITE(1, *) force(1, ni), force(2, ni), force(3, ni)
          ENDDO
          CLOSE(1)
          ! Wait for the NEWCAR or STOPCAR.
          WRITE(*,*) 'LINTERACT: Waiting for NEWCAR.'
        ENDIF
        DO
          INQUIRE(FILE = "STOPCAR", EXIST = stopcar_exists)
          INQUIRE(FILE = "NEWCAR", EXIST = newcar_exists)
          IF (stopcar_exists .OR. newcar_exists) exit
          CALL Sleep(1)
        ENDDO
        ! Read the NEWCAR in.
        IF (newcar_exists) THEN
          CALL RD_POSCAR_HEAD(LATT_CUR, T_I, NIOND, NIONPD, NTYPD, NTYPPD, IO%IU0, IO%IU6)
          CALL RD_POSCAR(LATT_CUR, T_I, DYN, NHC, NIOND, NIONPD, NTYPD, NTYPPD, IO%IU0, IO%IU6)
          posion = DYN%POSION
        ENDIF
#if defined(MPI) || defined(MPI_CHAIN)
        CALLMPI(MPI_Barrier(comm_chain%mpi_comm, ierr))
#endif
        IF (IU6>0) THEN       
          CALL unlink("NEWCAR")
        ENDIF
      ENDIF

#if defined(MPI) || defined(MPI_CHAIN)
          CALLMPI(MPI_Barrier(comm_chain%mpi_comm, ierr))
#endif

      ! /Interactive Mode

      ! MPMD mode
#ifdef VASP_MPMD
      IF (LMPMD .AND. IU6>0) THEN
        CALL mpmd_send(nions, toten, force)
        CALL mpmd_recv(nions, posion, b)
      ENDIF
#endif
      CALLMPI(MPI_Barrier(comm_chain%mpi_comm, ierr))
      ! End MPMD mode

      ! save the original stress tensor

      stress_vasp = stress

      ! if we are not using vasp optimizers we need to add 
      ! the pressure to the stress tensor
      IF (IOPT .NE. 0) THEN
        ! find volume of cell
        omega =  a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2))
        omega = omega - a(1,2)*(a(2,1)*a(3,3) - a(3,1)*a(2,3))
        omega = omega + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
        DO I=1,3
          stress(I,I) = stress(I,I) - omega*PSTRESS_local
        ENDDO
      ENDIF

      ! for dimer/lanczos, save the force and posion for the vasp stop criteria
      IF (ICHAIN==2 .OR. ICHAIN==3) THEN
        force_vasp = force
        posion_vasp = posion
      ENDIF

#ifdef EAM
      IF (.NOT. eaminit) THEN
        eaminit = .TRUE.
        Ream = posion
        Rvasp = posion
        Fvasp = force
        Uvasp = toten
      ENDIF
      posion = Ream
      IF (IU6>0) WRITE(IU6,*) 'Calling EAMForce'
      IF (IU6>0) WRITE(IU6,*) 'posion'
      IF (IU6>0) WRITE(IU6,'(3F14.6)') posion
      CALL EAMForce(IU6, nions, posion, Feam, Ueam, a, b)
      force = Feam
      toten = Ueam
      IF (IU6>0) WRITE(IU6,*) 'force'
      IF (IU6>0) WRITE(IU6,'(3F14.6)') force
      IF (IU6>0) WRITE(IU6,*) 'U ',Ueam
#endif

      ! optflag indicates who has control
      !  true: optimizer is active
      !  false: chain method is active (for dimer/lanczos)
     
      ! Initialize ml_chain and ml_chain_iter
      ml_chain_iter = 0
      ml_chain = .TRUE.
 
      ! Initialize ML-optimizer flags 
      IF ( IOPT == 8 ) THEN
        ml_turnoff = .FALSE.
        ! Save initial posion for the user-defiend trust radius 
        init_R = posion
      ELSE
        ml_turnoff = .TRUE.  
      END IF

      DO WHILE (ml_chain)
        ! Set ml_chain_iter and ml_chain
        ! ml_chain_iter >= 1 for ml-optimizer (iopt = 8) 
        ! ml_chain_iter = 0 otherwise 
        IF ((IOPT == 8) .AND. (.NOT. ml_turnoff)) THEN 
          ml_chain_iter = ml_chain_iter + 1
        ELSE 
          ml_chain = .FALSE.
        END IF

        ! ML-optimizer: Save true force and posion for training before taking chain-steps
        IF ((IOPT == 8) .AND. (.NOT. ml_turnoff)) THEN
          IF (ml_chain_iter == 1) THEN 
            force_vasp = force
            posion_vasp = posion
          ELSE
            ! Evaluate force on ML-PES when ml_chain_iter > 1
            CALL ml_pes(posion, a, b, Free, toten, force)
            ! For dimer and lanczos, store true force on ML-PES
            force_ml = force
            posion_ml = posion
          END IF
        END IF

        IF (ICHAIN == 0) THEN
          CALL neb_step(optflag, posion, toten, force, stress, a, b)
        ELSEIF (ICHAIN == 1) THEN
          CALL dynmat_step(optflag, posion, toten, force, a, b)
        ELSEIF (ICHAIN == 2) THEN
          IF (IOPT == 8) THEN
            IF (ml_chain_iter > 1) THEN
              ! Keep turning on optflag until ml-chain becomes false
              optflag = .TRUE.
            END IF  
          END IF
          CALL dimer_step(optflag, posion, toten, force, a, b)
          force_dimlan = force
        ELSEIF (ICHAIN == 3) THEN
          IF (IOPT == 8) THEN
            IF (ml_chain_iter > 1) THEN
              ! Keep turning on optflag until ml-chain becomes false
              optflag = .TRUE.
            END IF
          END IF
          CALL lanczos_step(optflag, posion, toten, force, a, b)
          force_dimlan = force
        ELSEIF (ICHAIN == 6) THEN
          CALL bbm_step(optflag, posion, toten, force, a, b)
!INS_BEGIN
        ELSEIF (ICHAIN == 4) THEN
          CALL instanton_step(optflag, posion, toten, force, a, b)
!INS_END
        ENDIF

        ! Write stress matrix for the first ml-chain loop only
        IF (ml_chain_iter <= 1) THEN 
          IF (IU6>0) THEN
            WRITE(IU6,*) "stress matrix after NEB project (eV)"
            WRITE(IU6,'(3F13.5)') stress
          ENDIF
        END IF

        ! 2D calculations (only in x,y)
        IF (twodim_flag) THEN
          stress(1,3) = 0._q
          stress(2,3) = 0._q
          stress(3,3) = 0._q
        ENDIF

        ! zero out any added forces on frozen atoms
        ! atomic constraints in direct coordinates, thanks to Christopher Stihl

        CALL KARDIR(NIONS, force, B)
        force = force*Free
        CALL DIRKAR(NIONS, force, A)

        ! for dimer/lanczos, check force criteria using the vasp force
        IF (ICHAIN == 2 .OR. ICHAIN == 3) THEN
          force_dimlan = force
          IF ((IOPT == 8) .AND. (ml_chain_iter > 1)) THEN
            force = force_ml
          ELSE
            force = force_vasp
          END IF
          !ML-optimizer only trains until optflag becomes true
          IF (IOPT == 8 .AND. (.NOT. optflag)) ml_chain = .FALSE. 
        ENDIF

        ftot = 0._q
        fmaxatom = 0._q
        fmaxdim = 0._q
        DO ni = 1,nions
          ftemp = 0._q
          DO nj = 1,3
            ftot = ftot + force(nj,ni)**2
            ftemp = ftemp + force(nj,ni)**2
            IF (fmaxdim .LT. ABS(force(nj,ni))) THEN
              fmaxdim = ABS(force(nj,ni))
            ENDIF
          ENDDO
          IF (ftemp .GT. fmaxatom) fmaxatom = ftemp
        ENDDO
        frms = SQRT(ftot/REAL(nions))
        fmaxatom = SQRT(fmaxatom)

        ! Check if ML-optimizer should be turned on or not 
        ! This check should happen when ml_chain_iter == 1 (when force is DFT chain force)
        IF (IOPT == 8) THEN
          IF ((ml_chain_iter == 1) .AND. (fmaxatom .LT. SWFTOL)) THEN
            ml_turnoff = .TRUE.
            ml_chain = .FALSE.
          ELSE
            ml_turnoff = .FALSE.
          END IF
#if defined(MPI) || defined(MPI_CHAIN)
        CALLMPI_C(and_chain(ml_turnoff))
        CALLMPI_C(and_chain(ml_chain))
#endif
        END IF

        ! Train ML-PES when ml_chain_iter=1 and ml_turnoff is false
        IF ( (.NOT. ml_turnoff) .AND. (ml_chain_iter == 1)) THEN 
          CALL ML_PyAMFF_train(posion_vasp, toten, force_vasp, a, b, Free) 
          ml_cleanup = .TRUE.
        END IF
        !For sanity check 
        !IF ( (IOPT == 8) .AND. (.NOT. ml_turnoff) .AND. (ml_chain_iter==1)) THEN
        !  Evaluate new osion's forces on ML-PES
        !  CALL ml_pes(posion,a,b,Free,toten_ml,force_ml)
        !END IF

        !ml_fmaxatom = 0._q
        !DO ni = 1,nions
        !  ftemp = 0._q
        !  DO nj = 1,3
        !    ftemp = ftemp + force_ml(nj,ni)**2
        !  ENDDO
        !  IF (ftemp .GT. ml_fmaxatom) ml_fmaxatom = ftemp
        !ENDDO
        !ml_fmaxatom = SQRT(ml_fmaxatom)

        ! Write fmaxatom and frms when ml_chain_iter <= 1 
        IF (ml_chain_iter <= 1) THEN 
          !dft_fmaxatom=fmaxatom 
          IF (IU6>=0) WRITE(IU6, 4693) fmaxatom, frms
 4693     FORMAT(1x,' FORCES: max atom, RMS ', 2f12.6)
          IF (IU6>=0) WRITE(IU6, 4694) SQRT(ftot), fmaxdim
 4694     FORMAT(1x,' FORCE total and by dimension', 2f12.6)
          !IF (IOPT == 8 .AND. IU6>=0) WRITE(IU6,4696) ml_fmaxatom

        ELSE
          IF (IU6>=0) WRITE(IU6, 4696) fmaxatom
 4696     FORMAT(1x,' ML-PyAMFF:',' ML-Forces: max atom ', 2f12.6)
          !IF (ICHAIN == 0) THEN 
          !  IF (IU6 >= 0) WRITE(IU6,*) 'ML-predicted true force fmaxatom=', ml_fmaxatom
          !END IF
        END IF

        ! added to monitor the stress vector 
        stot = 0._q
        smaxdim = 0._q
        DO ni = 1,3
          DO nj = 1,3
            stot = stot + stress(nj,ni)**2
            IF(smaxdim .LT. ABS(stress(nj,ni))) THEN
              smaxdim = ABS(stress(nj,ni))
            ENDIF
          ENDDO
        ENDDO
 4695   FORMAT(1x,' Stress total and by dimension', 2f12.6)

        ! Write stot, smaxdim when ml_chain_iter <= 1
        IF (ml_chain_iter <= 1) THEN
          IF (IU6>=0) WRITE(IU6,4695) SQRT(stot), smaxdim
        END IF

        ! for dimer/lanczos, use projected force in optimizer
        IF (ICHAIN==2 .OR. ICHAIN==3) THEN
          force = force_dimlan
        ENDIF

        ! Skip the following IF statement when ml_chain_iter > 1 
        IF (ml_chain_iter > 1) GOTO 147

        ! stops based on the Magnitude of the Force
        IF(ftot_flag) THEN
          fconverge = (SQRT(ftot) .LT. ftot_val)
          CALLMPI_C(and_chain(fconverge))
          IF (fconverge) THEN
            IF (IU6>=0) WRITE(IU6,*) 'CONVERGED based on Magnitude of Force'
            STOP
          ENDIF
        ENDIF

  147   CONTINUE

        IF (IOPT .NE. 0) THEN
          hstress = stress
          CALL sdotA(hstress,a)
          ! freeze out rotation
          hstress(2,1) = 0._q
          hstress(3,1) = 0._q
          hstress(3,2) = 0._q
          hstress = hstress/jacobian

          ! When ml_chain_iter > 1, fconverge is determined based on ML-PES
          fconverge = (fmaxatom .LT. ABS(EDIFFG_local))
          CALLMPI_C(and_chain(fconverge))
          IF (.NOT. fconverge) THEN
            ! our own optimizers (optflag: do or do not optimize)  
            CALL opt_step(optflag, posion, toten, force, hstress, a, b, &
                          ml_turnoff, Free, ml_chain_iter)
            ! When ml_chain_iter=1, save force for the vasp stop criteria 
            IF (ml_chain_iter == 1) THEN
              ! For dimer/lanczos, force should be true-force
              IF (ICHAIN==2 .OR. ICHAIN==3) THEN 
                force_mlreturn = force_vasp 
              ELSE
                force_mlreturn = force
              END IF
            ! Else, check if ml-pes should be trusted or not
            ELSE
              CALL ML_chain_criteria(nions, init_R, posion, A, iu6)
            END IF 
               
          ELSE
            IF (cell_flag) THEN
              sconverge = (smaxdim .LT. ABS(EDIFFG_local)*dble(nions))
              CALLMPI_C(and_chain(sconverge))
            ELSE
              sconverge = .TRUE.
            ENDIF
            IF (.NOT. sconverge) THEN
              CALL opt_step(optflag, posion, toten, force, hstress, a, b, &
                            ml_turnoff, Free, ml_chain_iter)
              ! The force returned to main.F should be based on DFT
              IF (ml_chain_iter == 1) THEN 
                ! For dimer/lanczos, force should be true-force
                IF (ICHAIN == 2 .OR. ICHAIN == 3) THEN
                  force_mlreturn = force_vasp
                ELSE
                  force_mlreturn = force
                END IF
              ! Else, check if ml-pes should be trustable or not
              ELSE
                CALL ML_chain_criteria(nions, init_R, posion, A, iu6)
              END IF
  
            ELSE
              IF (IU6>=0) WRITE(IU6,*) 'OPT: skip step - force has converged'
              ! Turn off ml-optimizer to avoid problem of deallocating not allocated arrays
              ml_turnoff = .TRUE.
              ml_chain = .FALSE.
              ! GH: remove any finite difference steps taken by dimer or lan
              IF (ICHAIN == 2 .OR. ICHAIN == 3) THEN
                IF (IOPT == 8 .AND. (ml_chain_iter > 1)) THEN 
                  posion = posion_ml
                ELSE
                  posion = posion_vasp
                END IF
              END IF            
              IF (ICHAIN == 2) THEN
                IF (IOPT == 8 .AND. (ml_chain_iter > 1)) THEN 
                  CONTINUE
                ELSE
                  CALL dimer_fin(posion,a,b)
                END IF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      END DO

#ifdef EAM
      Ream = posion
      ! copy saved variables back for vasp
      posion = Rvasp
      force = Fvasp
      toten = Uvasp
#endif

    IF (IOPT == 8 .AND. (ml_chain_iter > 1)) THEN
      ! Write ML-chain usage info 
      IF (iu6>=0) WRITE(iu6,'(A11,A29,I6)') 'ML-PyAMFF:', ' Number of ML-potential used: ',ml_chain_iter-1

      ! Return to the correct force
      force = force_mlreturn
      ! Reset fconverge to FALSE. This step is essential to check fconverge based on DFT not ML-PES
      fconverge = .FALSE.
    ELSE
      ! for dimer/lanczos, return the true force to vasp
      IF (ICHAIN == 2 .OR. ICHAIN == 3) THEN
        force = force_vasp
      END IF
    END IF

    IF (IOPT==8 .AND. (ml_cleanup)) CALL ML_PyAMFF_cleanup

    ! return unprojected stress vasp
    stress = stress_vasp

    END SUBROUTINE chain_force

    !>
    !> averaging the stress tensor when using
    !> VCAIMAGES.
    !>
    !> Function is only called for DYN%ISIF > 0
    !>@param stress stress tensor which will be averaged over
    !> the two molecular dynamics streams
    !>
    SUBROUTINE chain_stress( stress )

       implicit none
       REAL(q),intent( inout ) :: stress(3,3)
       REAL(q)                 :: stress_all(3,3,2)
       INTEGER                 :: node

       IF ( images == 0 ) RETURN
       IF ( spring == -1000 ) RETURN

#if defined(MPI) || defined (MPI_CHAIN)
       node = comm_chain%node_me

!======================================================================
! Parallel tempering return
!======================================================================
       IF ( LTEMPER ) THEN
          RETURN
!======================================================================
! VCA average stresses over two images
!======================================================================
       ELSE IF ( LVCAIMAGES ) THEN
          stress_all(:,:,1:2)  = 0
          stress_all(:,:,node) = stress
#ifndef old_vca
          ! first core in image merges the data
          IF ( comm%node_me == 1 ) THEN
             !! 3x3 stress components and two molecular dynamics streams
             CALLMPI_C( M_sum_d( comm_chain, stress_all(1,1,1), 3*3*2))
             stress(:,:) = stress_all(:,:,1)*VCAIMAGES + stress_all(:,:,2)*(1-VCAIMAGES)
          ELSE
             stress=0
          ENDIF
          ! now sum over all cores in one image
          CALLMPI_C( M_sum_d( comm, stress( 1, 1 ), 3*3 ) )
#else
          !! 3x3 stress components and two molecular dynamics streams
          CALLMPI_C( M_sum_d( comm_chain, stress_all( 1, 1, 1 ), 3*3*2 ) )
          stress = 0
          stress( : , : ) = stress_all( :, :, 1 )*VCAIMAGES + stress_all( :, :, 2 )*( 1 - VCAIMAGES )
          RETURN
#endif
       END IF
#endif
    END SUBROUTINE chain_stress

!**********************************************************************
!
! there are several points where a global sum between chains
! is required in order to get correct results
! terms like the total energy / or the kinetic energy
! make only sense in a global and not per image sense
!
!**********************************************************************

    SUBROUTINE sum_chain( value )
      REAL(q) :: value
      REAL(q) :: value_all(images)
      INTEGER node

      IF (images==0    ) RETURN
      IF (spring==-1000) RETURN
      IF (LTEMPER) RETURN
#if defined(MPI) || defined (MPI_CHAIN)
      ! VCAIMAGES returns average value for energies etc.
      IF (LVCAIMAGES) THEN
         node=comm_chain%node_me
         value_all=0
         value_all(node)=value
#ifndef old_vca
         ! first core in image merges data
         IF (comm%node_me==1) THEN
            CALLMPI_C( M_sum_d( comm_chain, value_all, 2))
            value=value_all(1)*VCAIMAGES+value_all(2)*(1-VCAIMAGES)
         ELSE
            value=0
         ENDIF
         ! now sum over all cores in one image
         CALLMPI_C( M_sum_d( comm, value, 1))
#else
         CALLMPI_C( M_sum_d( comm_chain, value_all, 2))
         value=value_all(1)*VCAIMAGES+value_all(2)*(1-VCAIMAGES)
#endif
      ELSE
         CALLMPI_C( M_sum_d( comm_chain, value, 1 ))
      ENDIF
#endif

    END SUBROUTINE sum_chain

!**********************************************************************
!
! also the logical break conditionions must be
!
!**********************************************************************

    SUBROUTINE PARALLEL_TEMPERING(NSTEP, NIONS, POSION, VEL, TOTEN, FORCE, TEBEG, TEEND, &
        A, B, IU6)
      USE constant
      IMPLICIT NONE
      INTEGER :: NSTEP           ! step
      INTEGER :: NIONS           ! number of ions
      INTEGER :: IU6             ! std-out  (unit for OUTCAR file)
      REAL(q) :: POSION(3,nions) ! position of ions
      REAL(q) :: VEL(3,nions)    ! velocity of ions
      REAL(q) :: FORCE(3,nions)
      REAL(q) :: TOTEN           ! total energy
      REAL(q) :: TEBEG           ! temperature of ensemble
      REAL(q) :: TEEND           ! temperature of ensemble (must equal TEBEG)
      REAL(q) :: A(3,3),B(3,3)   ! lattice and reciprocal lattice
! local
      INTEGER,SAVE :: SWAP_START=1
      REAL(q),ALLOCATABLE :: TEBEG_old(:), TEBEG_new(:), TOTEN_all(:)
      REAL(q) :: E1, E2
      INTEGER,ALLOCATABLE :: ID(:)
      INTEGER :: SWAP
      INTEGER :: NODE, I, J
      REAL(q) :: TEBEG_new_local, E
      REAL(q) :: value

#if defined(MPI) || defined (MPI_CHAIN)
      IF (.NOT. LTEMPER .OR. NTEMPER==0) RETURN

      IF (NTEMPER_NEXT>0) THEN
         NTEMPER_NEXT=NTEMPER_NEXT-1
      ELSE
         CALL RANDOM_NUMBER(value)
         ! broadcast to all nodes
         value=value*NTEMPER*2-1
         CALLMPI_C( M_bcast_d( comm_world, value, 1))
         NTEMPER_NEXT=value
         IF (IU6>=0) WRITE(IU6,'(A,I8,A,I8,A,I8)') ' parallel tempering entered NSTEP= ',NSTEP,' NTEMPER= ', NTEMPER,' NTEMPER_NEXT = ',NTEMPER_NEXT
      ENDIF

      IF (NTEMPER_NEXT == 0) THEN
         IF (IU6>=0) WRITE(IU6,'(A,3I8)') ' parallel tempering performing NSTEP= ', NSTEP

         ALLOCATE(TEBEG_old(images), TEBEG_new(images), TOTEN_all(images), ID(images))

         ! my id, i.e. which subdir I am running in
         node=comm_chain%node_me

         TEBEG_old=0
         TEBEG_old(node)=TEBEG
         TOTEN_all=0
         TOTEN_all(node)=TOTEN
         ID=0
         ID(node)=node

         CALLMPI_C( M_sum_d( comm_chain, TEBEG_old,  images))
         CALLMPI_C( M_sum_d( comm_chain, TOTEN_all,  images))
         CALLMPI_C( M_sum_i( comm_chain, ID,  images))

! sort the temperatures ascendingly and store the corresponding node id 
         CALL SORT_ASC_REAL(images, TEBEG_old, ID)
!======================================================================
! now select swaps randomly
! images/2 swaps 
! the array SWAP stores the list of images to be swapped 
! 0 no swap
! 1        two lowest temperatures are swapped
! 2        next two are swapped and so on
! images-1 swap last two
!======================================================================
! old code, one swap attempt each step
!         CALL RANDOM_NUMBER(value)  ! values between [0,1[
!         CALLMPI_C( M_bcast_d( comm_world, value, 1))
!         SWAP=(images-1)*value+1    ! create a random number [1,images-1]
!         CALLMPI_C( M_bcast_i( comm_chain, SWAP,  1))
! new version
       TEBEG_new=TEBEG_old

       SWAP_START=3-SWAP_START ! start swap either at first or second image

       DO SWAP=SWAP_START,images-1,2

         IF (IU6>=0) WRITE(IU6,'(A,16I10)')   ' attempting swapping', SWAP

         ATTEMPTS(SWAP)=ATTEMPTS(SWAP)+1

         E=(1.0_q/TEBEG_old(SWAP)-1.0_q/TEBEG_old(SWAP+1))/ BOLKEV * (TOTEN_all(ID(SWAP))-TOTEN_all(ID(SWAP+1)))
         E1=(1.0_q/TEBEG_old(SWAP)-1.0_q/TEBEG_old(SWAP+1))/ BOLKEV
         E2=(TOTEN_all(ID(SWAP))-TOTEN_all(ID(SWAP+1)))
         E=EXP(E)
         CALL RANDOM_NUMBER(value) ! values between [0,1[
         CALLMPI_C( M_bcast_d( comm_world, value, 1))

         IF (IU6>=0) WRITE(IU6,'(A,16F10.4)') '  1/T1-1/T2         ', E1
         IF (IU6>=0) WRITE(IU6,'(A,16F10.4)') '  E1  -E2           ', E2
         IF (IU6>=0) WRITE(IU6,'(A,16F10.7)')   '            random  ', value

         IF (value>E) THEN
            ! Metropolis forbids swap
            IF (IU6>=0) WRITE(IU6,'(A,16I10)')   '          noswapping', SWAP
         ELSE
            SUCCESS(SWAP)=SUCCESS(SWAP)+1
            TEBEG_new(SWAP)  =TEBEG_old(SWAP+1)
            TEBEG_new(SWAP+1)=TEBEG_old(SWAP)
            IF (IU6>=0) WRITE(IU6,'(A,16I10)')   '            swapping', SWAP
         ENDIF

      ENDDO

         IF (IU6>=0) WRITE(IU6,'(A,99F14.7)') ' parallel tempering old TOTEN ', (TOTEN_all(ID(I)),I=1,images)
         IF (IU6>=0) WRITE(IU6,'(A,99F14.7)') ' parallel tempering old TEBEG ', TEBEG_old
         IF (IU6>=0) WRITE(IU6,'(A,99F14.7)') ' parallel tempering new TEBEG ', TEBEG_new
         IF (IU6>=0) WRITE(IU6,'(A,99F14.7)') ' Acceptance ratio for swaps          ', REAL(SUCCESS,q)/MAX(1,ATTEMPTS)

         DO I=1,images
            IF ( ID(I)==node) THEN
               IF (TEBEG /= TEBEG_old(I)) THEN
                  WRITE(*,*) 'internal error in PARALELL_TEMPERING:', I,ID(I), TEBEG, TEBEG_old(I)
               ENDIF
               TEBEG_new_local=TEBEG_new(I)
            ENDIF
         END DO

         VEL(:,:)=VEL(:,:)*SQRT(TEBEG_new_local/TEBEG)

         TEBEG=TEBEG_new_local
         TEEND=TEBEG_new_local

         IF (IU6>=0) WRITE(IU6,*)
         IF (IU6>=0) WRITE(IU6,"('   TEBEG  = ',F6.1,';   TEEND  =',F6.1)") TEBEG,TEEND
         IF (IU6>=0) WRITE(IU6,*)

         DEALLOCATE(TEBEG_old, TEBEG_new, TOTEN_all, ID)
      ENDIF

#endif
    END SUBROUTINE PARALLEL_TEMPERING


!**********************************************************************
! sorts RA in descending order, and rearanges an index array RB
! seems to be a quicksort, but I am not sure
! subroutine writen by Florian Kirchhof
!**********************************************************************

  SUBROUTINE SORT_ASC_REAL(N, RA, RB)
    REAL(q) :: RA(N)
    INTEGER :: RB(N)
    REAL(q) :: RRA
    INTEGER :: RRB
    INTEGER :: N, L, IR, J, I

    IF (N<=1) RETURN

    L=N/2+1
    IR=N
10  CONTINUE
    IF(L.GT.1)THEN
       L=L-1
       RRA=RA(L)
       RRB=RB(L)
    ELSE
       RRA=RA(IR)
       RRB=RB(IR)
       RA(IR)=RA(1)
       RB(IR)=RB(1)
       IR=IR-1
       IF(IR.EQ.1)THEN
          RA(1)=RRA
          RB(1)=RRB
          RETURN
       ENDIF
    ENDIF
    I=L
    J=L+L
20  IF(J.LE.IR)THEN
       IF(J.LT.IR)THEN
          IF(RA(J).GT.RA(J+1))J=J+1
       ENDIF
       IF(RRA.GT.RA(J))THEN
          RA(I)=RA(J)
          RB(I)=RB(J)
          I=J
          J=J+J
       ELSE
          J=IR+1
       ENDIF
       GO TO 20
    ENDIF
    RA(I)=RRA
    RB(I)=RRB
    GO TO 10
  END SUBROUTINE SORT_ASC_REAL

!**********************************************************************
!
! also the logical break conditionions must be
!
!**********************************************************************

    SUBROUTINE and_chain( value )
      LOGICAL :: value
      REAL(q) :: sum

      IF (images==0) RETURN
      IF (spring==-1000) RETURN
#if defined(MPI) || defined(MPI_CHAIN)
! if one node is .FALSE., .FALSE. is returned on all nodes
      IF (value) THEN
         sum=0
      ELSE
         sum=1
      ENDIF
      IF (LVCAIMAGES) THEN
#ifndef old_vca
         ! first core in image merges data
         IF (comm%node_me==1) THEN
            CALLMPI_C( M_sum_d( comm_chain, sum, 1 ))
         ELSE
            sum=0
         ENDIF
         ! now sum over all cores in one image
         CALLMPI_C( M_sum_d( comm, sum, 1))
#else
         CALLMPI_C( M_sum_d( comm_chain, sum, 1 ))
#endif
      ELSE
         CALLMPI_C( M_sum_d( comm_chain, sum, 1 ))
      ENDIF

      IF (sum>=1) THEN
         value=.FALSE.
      ELSE
         value=.TRUE.
      ENDIF
#endif
    END SUBROUTINE and_chain

    SUBROUTINE or_chain( value )
      LOGICAL :: value
      REAL(q) :: sum

      IF (images==0) RETURN
      IF (spring==-1000) RETURN
#if defined(MPI) || defined(MPI_CHAIN)
! if one node is .FALSE., .FALSE. is returned on all nodes
      IF (value) THEN
         sum=1
      ELSE
         sum=0
      ENDIF
      IF (LVCAIMAGES) THEN
#ifndef old_vca
         ! first core in image merges data
         IF (comm%node_me==1) THEN
            CALLMPI_C( M_sum_d( comm_chain, sum, 1 ))
         ELSE
            sum=0
         ENDIF
         ! now sum over all cores in one image
         CALLMPI_C( M_sum_d( comm, sum, 1))
#else
         CALLMPI_C( M_sum_d( comm_chain, sum, 1 ))
#endif
      ELSE
         CALLMPI_C( M_sum_d( comm_chain, sum, 1 ))
      ENDIF

      IF (sum>=1) THEN
         value=.TRUE.
      ELSE
         value=.FALSE.
      ENDIF
#endif
    END SUBROUTINE or_chain

!**************** SUBROUTINE s2ts ************************************
! transform stress tensor to force on vectors
! true_stress = (stress dot B)^T
! ) reciprocal lattice  (BASIS must be equal to A direct lattice)
!***********************************************************************

      SUBROUTINE s2ts(V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER N
      DIMENSION V(3,3),BASIS(3,3)

      DO N=1,3
        V1 = V(N,1)*BASIS(1,1) + V(N,2)*BASIS(2,1) + V(N,3)*BASIS(3,1)
        V2 = V(N,1)*BASIS(1,2) + V(N,2)*BASIS(2,2) + V(N,3)*BASIS(3,2)
        V3 = V(N,1)*BASIS(1,3) + V(N,2)*BASIS(2,3) + V(N,3)*BASIS(3,3)
        V(N,1) = V1
        V(N,2) = V2
        V(N,3) = V3
      ENDDO

      RETURN
      END SUBROUTINE s2ts

!**************** SUBROUTINE ts2s ************************************
! transforms force on vectors to stress tensor
!  stress=true_stress^T dot A
!  direct lattice      (BASIS must be equal to A direct lattice)
! to cartesian coordinates
!***********************************************************************

      SUBROUTINE ts2s(V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER N
      DIMENSION V(3,3),BASIS(3,3)

      DO N=1,3
        V1 = V(N,1)*BASIS(1,1) + V(N,2)*BASIS(1,2) + V(N,3)*BASIS(1,3)
        V2 = V(N,1)*BASIS(2,1) + V(N,2)*BASIS(2,2) + V(N,3)*BASIS(2,3)
        V3 = V(N,1)*BASIS(3,1) + V(N,2)*BASIS(3,2) + V(N,3)*BASIS(3,3)
        V(N,1) = V1
        V(N,2) = V2
        V(N,3) = V3
      ENDDO

      RETURN
      END SUBROUTINE ts2s

      SUBROUTINE sdotA(V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER I,J,K
      DIMENSION V(3,3),BASIS(3,3),AC(3,3)
      AC = 0
      DO J = 1,3
         DO I = 1,3
            !A(I,J) = AC(I,J)
            DO K = 1,3
               AC(I,J) = AC(I,J) + V(I,K)*BASIS(K,J)
            ENDDO
         ENDDO
      ENDDO
      V = AC

      RETURN
      END SUBROUTINE sdotA

      SUBROUTINE sdotB(V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER I,J,K
      DIMENSION V(3,3),BASIS(3,3),AC(3,3)
      AC = 0
      DO J = 1,3
         DO I = 1,3
            !A(I,J) = AC(I,J)
            DO K = 1,3
               AC(I,J) = AC(I,J) + V(I,K)*BASIS(J,K)
            ENDDO
         ENDDO
      ENDDO
      V = AC

      RETURN
      END SUBROUTINE sdotB

!======================================================================
! Returns a unit vector along v1
!======================================================================
      FUNCTION return_unit(v1,nmax)
      INTEGER :: nmax
      REAL(q) :: v1(3,nmax)
      REAL(q),DIMENSION(3,nmax) :: return_unit
      IF (SUM(v1*v1) .NE. 0._q) THEN
        return_unit = v1*(1._q/SQRT(SUM(v1*v1)))
      ELSE
        return_unit = v1*0._q
      ENDIF
      END FUNCTION return_unit

#ifdef EAM
!**********************************************************************
! Empirical potential, for testing purposes only
!**********************************************************************

      SUBROUTINE EAMForce(IU6, NIONS, POSION, Feam, Ueam, A, B)
      INTEGER :: NIONS, IU6
      REAL(q) :: POSION(3,NIONS)
      REAL(q) :: Reamcar(3,NIONS), Feam(3,NIONS)
      REAL(q) :: A(3,3), B(3,3)
      INTEGER :: NI, NJ

      REAL*8 :: Rtmp(3*NIONS), Ftmp(3*NIONS)
      REAL*8 :: Ueam
      REAL*8 :: Ax,Ay,Az

      EXTERNAL EAM_FORCE

      Ax = A(1,1)
      Ay = A(2,2)
      Az = A(3,3)

      IF (IU6>0) WRITE(IU6,*) 'IN EAMForce'
      ! Convert positions to Cartesian
      Reamcar = POSION
      CALL DIRKAR(NIONS, Reamcar, A)
      IF (IU6>0) WRITE(IU6,*) 'Reamcar'
      IF (IU6>0) WRITE(IU6,'(3F14.8)') Reamcar(:,:)
      ! 
      DO NI = 1,NIONS
         DO NJ = 1,3
            Rtmp(3*(NI-1)+NJ) = Reamcar(NJ,NI)
         ENDDO
      ENDDO
      call EAM_FORCE(Nions, Rtmp, Ftmp, Ueam, Ax, Ay, Az)
      DO NI = 1,NIONS
         DO NJ = 1,3
            Feam(NJ,NI) = Ftmp(3*(NI-1)+NJ)
         ENDDO
      ENDDO
      ! Apply constraints in direct coordinates
      CALL KARDIR(NIONS, Feam, B)
      Feam = Feam*Free
      CALL DIRKAR(NIONS, Feam, A)

      IF (IU6>0) WRITE(IU6,*) 'EAMForce'
      IF (IU6>0) WRITE(IU6,'(3F14.8)') FEAM(:,:)
      IF (IU6>0) WRITE(IU6,*) 'Ueam: ',Ueam
      END SUBROUTINE EAMForce
#endif


!**********************************************************************
! Initialize the chain (repeated image mode) and determine which of the
! three possible methods to use based on the ICHAIN variable:
!   ICHAIN==0: nudged elastic band (default)
!   ICHAIN==1: dynamical matrix
!   ICHAIN==2: dimer method
!   ICHAIN==3: lanczos method
!   ICHAIN==6: bond-boost method
!**********************************************************************

    SUBROUTINE chain_init(T_INFO, IO)
      USE base
      USE reader_tags
      USE tutor, ONLY: vtutor

      IMPLICIT NONE

      TYPE (in_struct) :: IO
      TYPE (type_info) :: T_INFO

! needed only temporarily
      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      TYPE (latt)      :: LATT_CUR
      TYPE (type_info) :: T_I
      TYPE (dynamics)  :: DYN
      TYPE(nh_chains)  :: NHC
      INTEGER     IDUM,IERR,N,idir,node
      CHARACTER (1)   CHARAC
      COMPLEX(q)  CDUM  ; LOGICAL  LDUM
      REAL(q) :: RDUM

!!!
      INTEGER :: NI, NJ, IU0 ,IU6
      REAL(q) :: omega
      REAL, PARAMETER :: EVTOJ=1.60217733E-19_q

#ifdef EAM
      EXTERNAL EAM_POTINIT
#endif
      IU0 = IO%IU0
      IU6 = IO%IU6

      ! write the version number
      IF (IU6>=0) WRITE(IU6,'(/,A,/)') ' VTST: version 6.6.0, (04/01/26)'

      IF (IU6>=0) WRITE(IU6,*) 'CHAIN: initializing optimizer'

      ! initialize optimizer
      CALL opt_init(T_INFO, IO)
      optflag = .FALSE.

      IF (IOPT == 8) THEN
        ml_turnoff = .FALSE.
      ELSE
        ml_turnoff = .TRUE.
      END IF
  
      ! initialize chain based method
      ICHAIN = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'ICHAIN','=','#',';','I', &
     &            ICHAIN,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''ICHAIN'' from file INCAR.'
        STOP
      ENDIF

      LINTERACT = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LINTERACT','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LINTERACT,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LINTERACT'' from file INCAR.'
        STOP
      ENDIF

      LMPMD = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LMPMD','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LMPMD,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LMPMD'' from file INCAR.'
        STOP
      ENDIF

      IOPT = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'IOPT','=','#',';','I', &
     &            IOPT,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''IOPT'' from file INCAR.'
        STOP
      ENDIF

      !IF ML-PyAMFF is used for optmizer, read input variables for chain.F
      IF (IOPT == 8) THEN
        force_tol=0.010_q
        CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_FTOL','=','#',';','F', &
                  IDUM,force_tol,CDUM,LDUM,CHARAC,N,1,IERR)
        IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
          IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_FTOL'' from file INCAR.'
          STOP
        ENDIF

        ! By default, ML-optimizer is switched over to L-BFGS when fmaxatom < 0.050 
        SWFTOL=0.050_q
        CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_SWFTOL','=','#',';','F', &
                    IDUM,SWFTOL,CDUM,LDUM,CHARAC,N,1,IERR)
        IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
          IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_SWFTOL'' from file INCAR.'
          STOP
        END IF

        trustR=0.5_q
        CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_MAXMOVE','=','#',';','F', &
                  IDUM,trustR,CDUM,LDUM,CHARAC,N,1,IERR)
        IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
          IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_MAXMOVE'' from file INCAR.'
         STOP
        ENDIF

        minmove=0.001_q
        CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_MINMOVE','=','#',';','F', &
                  IDUM,minmove,CDUM,LDUM,CHARAC,N,1,IERR)
        IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
          IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_MINMOVE'' from file INCAR.'
          STOP
        ENDIF

        max_iter=30
        CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PYAMFF_MAXITER','=','#',';','I', &
                  max_iter,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
        IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
          IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''PYAMFF_MAXITER'' from file INCAR.'
          STOP
        ENDIF

      END IF

      EDIFFG_local = 0.1_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'EDIFFG','=','#',';','F', &
     &            IDUM,EDIFFG_local,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''EDIFFG'' from file INCAR.'
        STOP
      ENDIF

!     determines if cell should change in NEB
      cell_flag = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LNEBCELL','=','#',';','L', &
     &            IDUM,RDUM,CDUM,cell_flag,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LNEBCELL'' from file INCAR.'
        STOP
      ENDIF

!     sets stress tensor to two dimensions 
      twodim_flag = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LTWODIM','=','#',';','L', &
     &            IDUM,RDUM,CDUM,twodim_flag,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LTWODIM'' from file INCAR.'
        STOP
      ENDIF

!     USED to converge based on Magnitude of the force
      ftot_flag = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FMAGFLAG','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ftot_flag,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''FMAGFLAG'' from file INCAR.'
        STOP
      ENDIF

      ftot_val = 0.01_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FMAGVAL','=','#',';','F', &
     &            IDUM,ftot_val,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''FMAGVAL'' from file INCAR.'
        STOP
      ENDIF

      ! make sure that convergence is force based when using IOPT
      IF((IOPT .NE. 0) .AND. (EDIFFG_local .GT. 0.0_q)) THEN
        IF(IU6>=0) WRITE(IU6,*) 'Must set  EDIFFG < 0 when using IOPT > 0'
        STOP
      ENDIF

      ! check that solid state NEB uses our optimizers 
      IF((cell_flag .EQV. .TRUE.) .AND. ((IOPT .NE. 3) .AND. (IOPT .NE. 7))) THEN
        IF(IU6>=0) WRITE(IU6,*) 'Must set  IOPT = 3 or 7 when using LNEBCELL=.TRUE.'
        STOP
      ENDIF

      PSTRESS_local = 0.0_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'PSTRESS','=','#',';','F', &
     &            IDUM,PSTRESS_local,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''PSTRESS'' from file INCAR.'
        STOP
      ENDIF
      PSTRESS_local = PSTRESS_local/(EVTOJ*1E22_q)

      ISIF_local = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'ISIF','=','#',';','I', &
     &            ISIF_local,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''ISIF'' from file INCAR.'
        STOP
      ENDIF

      ! check that if ISIF is set, either vasp optimizers are used or IOPT=3,7
      IF ((IOPT>0) .AND. (ISIF_local>2)) THEN
        IF (ISIF_local/=3) THEN
          IF(IU6>=0) WRITE(IU6,*) 'Only ISIF = 3 is supported with IOPT = 3 or 7'
          STOP
        ELSE IF ((IOPT.NE.3) .AND. (IOPT.NE.7)) THEN
          IF (IU6>=0) WRITE(IU6,*) 'Must set  IOPT = 3 or 7 when using ISIF = 3'
          IF (IU6>=0) WRITE(IU6,*) 'Alternatively, use the vasp optimizers'
          STOP
        ENDIF
      ENDIF

      IF (IU6>=0) WRITE(IU6,*) 'CHAIN: Read ICHAIN ',ICHAIN

      IF (IMAGES==0) THEN
          CALL RD_POSCAR_HEAD(LATT_CUR, T_I, &
               NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)
          CALL RD_POSCAR(LATT_CUR, T_I, DYN, NHC, &
               NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)

          IF (T_I%NIONS /= T_INFO%NIONS) THEN
             IF (IU0>=0) WRITE(IU0,*)'ERROR: image mode number of ions wrong'
          STOP
          ENDIF
          omega = LATT_CUR%OMEGA
          jacobian = (omega/DBLE(T_INFO%nions))**(1._q/3._q)*sqrt(DBLE(T_INFO%nions))
      ENDIF

      IF (ICHAIN==0) THEN
        IF (IMAGES>0) THEN
          jacobian = 0.0_q
          IF (IU6>=0) WRITE(IU6,*) 'CHAIN: Running the NEB'
        ENDIF
        CALL neb_init(T_INFO, IO, jacobian)
      ELSEIF (ICHAIN==1) THEN
        IF (IU6>=0) WRITE(IU6,*) 'CHAIN: Running the Dynamical Matrix'
        CALL dynmat_init(T_INFO, IO)
      ELSEIF (ICHAIN==2) THEN
        IF (IU6>=0) WRITE(IU6,*) 'CHAIN: Running the Dimer method'
        CALL dimer_init(T_INFO, IO)
        IF (IOPT==0) THEN
          IF (IU6>=0) WRITE(IU6,*) 'CHAIN: Must set IOPT>0 to use the Dimer method'
          STOP
        ENDIF
      ELSEIF (ICHAIN==3) THEN
        IF (IU6>=0) WRITE(IU6,*) 'CHAIN: Running the Lanczos method'
        CALL lanczos_init(T_INFO,IO)
        IF (IOPT==0) THEN
          IF (IU6>=0) WRITE(IU6,*) 'CHAIN: Must set IOPT>0 to use the Lanczos method'
          STOP
        ENDIF
      ELSEIF (ICHAIN==4) THEN
        IF (IO%IU6>=0) WRITE(IU6,*) 'CHAIN: Running the Instanton method'
        CALL instanton_init(T_INFO,IO)
      ELSEIF (ICHAIN==6) THEN
      IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Running the bond-boost method'
      CALL bbm_init(T_INFO,IO)
      ENDIF

      ! Make vector to zero out force on frozen atoms 
      ALLOCATE(Free(3, T_INFO%nions))
      Free(:,:) = 1._q
      IF (T_INFO%LSDYN) THEN
        DO NI = 1,T_INFO%nions
          DO NJ = 1,3
            IF (.NOT.T_INFO%LSFOR(NJ,NI)) Free(NJ,NI) = 0._q
          ENDDO
        ENDDO
      ENDIF

#ifdef VASP_MPMD
      IF (LMPMD .EQV. .TRUE.) THEN
        IF (IO%IU6>=0) WRITE(IU6,*) 'CHAIN: MPMD Mode'
        mpmd_client_rank = comm_chain%mpmd_client_rank
      ENDIF
#endif

#ifdef EAM
      CALL EAM_POTINIT()
      eaminit = .FALSE.
      ALLOCATE(Feam(3,T_INFO%nions), Ream(3,T_INFO%nions))
      ALLOCATE(Fvasp(3,T_INFO%nions), Rvasp(3,T_INFO%nions))
#endif

!!!
! quick return, if we are not running in image mode
      IF (images==0) RETURN
#if defined(MPI) || defined(MPI_CHAIN)
! read the VCAIMAGES
      VCAIMAGES=-1

      CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'VCAIMAGES', VCAIMAGES, IERR, WRITEXMLINCAR)

      IF (VCAIMAGES==-1) THEN
         LVCAIMAGES=.FALSE.
      ELSE
         LVCAIMAGES=.TRUE.
      ENDIF


! LTEMPER -- use subspace diagonalization or not (default is TRUE):
      LTEMPER=.FALSE.
      CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'LTEMPER', LTEMPER, IERR, WRITEXMLINCAR)

      IF (LTEMPER) THEN
         ! read NTEMPER
         NTEMPER=200
         CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'NTEMPER', NTEMPER, IERR, WRITEXMLINCAR)

         ALLOCATE(ATTEMPTS(images-1))
         ALLOCATE(SUCCESS(images-1))
         ATTEMPTS=0
         SUCCESS=0
      ENDIF

      IF (LVCAIMAGES .OR. LTEMPER) THEN
         nions_max=T_INFO%NIONS
#ifndef old_vca
         IF (comm%node_me==1) THEN
            CALL M_max_i(comm_chain, nions_max, 1 )
         ELSE
            nions_max=0
         ENDIF
         ! now sum over all cores in one image
         CALLMPI_C( M_sum_i( comm, nions_max, 1))
#else
         CALL M_max_i(comm_chain, nions_max, 1 )
#endif
         IF (T_INFO%NIONS /= nions_max) THEN
            CALL vtutor%error("ERROR: image mode number of ions wrong")
         ENDIF
      ELSE

! allocate work array for all positions posion_all
         ALLOCATE(posion_all(3,t_info%nions,0:images+1))
! default spring constant
         spring=-5.
! read the spring constant
         CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'SPRING', SPRING, IERR, WRITEXMLINCAR)

! read the start and end points to posion_all
! read 00/POSCAR file, a little bit of fiddling is required
         idir=0
         CALL MAKE_DIR_APP(idir)

         CALL RD_POSCAR_HEAD(LATT_CUR, T_I, &
              NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)
         CALL RD_POSCAR(LATT_CUR, T_I, DYN, NHC, &
              NIOND,NIONPD, NTYPD,NTYPPD, &
              IO%IU0, IO%IU6)

         IF (T_I%NIONS /= T_INFO%NIONS) THEN
            CALL vtutor%error("ERROR: image mode number of ions wrong")
         ENDIF
         posion_all(:,:,0)= DYN%POSION
! read images+1/POSCAR file
         idir=images+1
         CALL MAKE_DIR_APP(idir)

         CALL RD_POSCAR_HEAD(LATT_CUR, T_I, &
              NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)
         CALL RD_POSCAR(LATT_CUR, T_I, DYN, NHC, &
              NIOND,NIONPD, NTYPD,NTYPPD, &
              IO%IU0, IO%IU6)

         IF (T_I%NIONS /= T_INFO%NIONS) THEN
            CALL vtutor%error("ERROR: image mode number of ions wrong")
         ENDIF
         posion_all(:,:,images+1)= DYN%POSION
      ENDIF

      node=comm_chain%node_me
      CALL MAKE_DIR_APP(node)
#endif
    END SUBROUTINE chain_init

!**********************************************************************
!
! returns true if hyper nudged elastic band method is used
!
!**********************************************************************

    FUNCTION LHYPER_NUDGE()
      LOGICAL LHYPER_NUDGE
      IF (images==0 .OR. spring /= 0 ) THEN
         LHYPER_NUDGE=.FALSE.
      ELSE
         LHYPER_NUDGE=.TRUE.
      ENDIF

    END FUNCTION LHYPER_NUDGE


!**********************************************************************
!
!> This subroutine automatically synchronizes some tag values over all
!> images.
!
!**********************************************************************

    SUBROUTINE SYNC_VCAIMAGES(DYN, LBLOCK_HELP, COMM)
       USE tutor, ONLY: vtutor
       USE poscar_struct_def, ONLY: dynamics
       IMPLICIT NONE
       TYPE (dynamics), INTENT(INOUT) :: DYN
       INTEGER, INTENT(INOUT)         :: LBLOCK_HELP
       TYPE (communic), INTENT(IN)    :: COMM

       ! Compute minimum of LBLOCK_HELP and NBLOCK over all images.
       CALLMPI_C( M_min_i(COMM, LBLOCK_HELP, 1) )
       CALLMPI_C( M_min_i(COMM, DYN%NBLOCK, 1) )

    END SUBROUTINE SYNC_VCAIMAGES

    SUBROUTINE ML_chain_criteria(nions, init_R, posion, A, iu6)
      ! ML-optimizer criteira of relaxing geomtery on ML-potential
      ! This subroutine determines ml_chain 

      ! Inputs
      INTEGER :: nions, iu6
      REAL(q),DIMENSION(3,3) :: A
      REAL(q), DIMENSION(3,nions) :: init_R, posion
      ! Variables
      INTEGER :: ni, nj
      REAL(q) :: dmaxatom, dtemp, fmaxatom, ftemp
      REAL(q), DIMENSION(3,nions) :: change_in_R
      LOGICAL :: not_trustable

      ! Initiate logical
      not_trustable = .FALSE.

      ! Calculate maximum move per atom
      ! Here, positions are all in direct 
      dmaxatom = 0._q
      change_in_R=posion-init_R
      ! change direct to cartesian
      CALL DIRKAR(NIONS, change_in_R, A)
      DO ni = 1, nions
        dtemp = 0._q
        DO nj = 1,3
          dtemp = dtemp + change_in_R(nj,ni)**2
        END DO
        IF (dtemp .GT. dmaxatom) dmaxatom = dtemp
      END DO
      dmaxatom = SQRT(dmaxatom)

      ! If any of dmaxatom is greater than user-defined trustR, ml-chain is false
      IF (dmaxatom .GT. trustR) not_trustable = .TRUE.
      CALLMPI_C(or_chain(not_trustable)) 

      !IF (iu6>=0) WRITE(iu6,*) 'dmaxatom=', dmaxatom  
      ! If dmaxatom at the first step is too small, we will turn off ml-opt
      ! The frms value becomes SWFTOL value 
      !IF ( (ml_chain_iter == 2) .AND. (dmaxatom .LT. minmove) ) THEN
      !  IF (iu6>=0) WRITE(iu6,*) 'The first step on the ML-PES is less than minmove =', minmove 
      !  ml_minmove = .TRUE. 
      !END IF
      !CALLMPI_C(and_chain(ml_minmove))

      ! Determine ml_chain value
      IF ( (not_trustable) .OR. &
           (ml_chain_iter-1 .GE. max_iter) ) THEN
        ml_chain = .FALSE.
      END IF

      CALLMPI_C(or_chain(ml_chain)) 

    END SUBROUTINE ML_chain_criteria
    
END MODULE chain
