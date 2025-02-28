MODULE fnnmodule
    USE fpType
    USE neuralnetwork
    USE nnType
    USE atomsProp
    USE normalize
    !USE fpCalc
    IMPLICIT NONE
    PUBLIC

    CONTAINS

    SUBROUTINE prepfNN(natoms, nelement, atomicNumbers, uniqueNrs) BIND(C,name='prepfNN')
        USE, INTRINSIC :: iso_c_binding
        IMPLICIT NONE
        ! Inputs
        INTEGER(c_long) :: natoms, i, j
        INTEGER(c_int) :: nelement !max_nneighs = MAXVAL(num_neigh)
        INTEGER(c_int), DIMENSION(nAtoms) :: atomicNumbers
        INTEGER(c_int), DIMENSION(nelement) :: uniqueNrs
        INTEGER :: nAt
        INTEGER, DIMENSION(natoms) :: symbols
        CHARACTER*2, DIMENSION(nelement) :: uniq_elements

        CHARACTER*2, DIMENSION(92) :: elementArray

        DATA elementArray / "H","He","Li","Be","B","C","N","O", &
                 "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc", &
                 "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se", &
                 "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag", &
                 "Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd", &
                 "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta", &
                 "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", &
                 "Fr","Ra","Ac","Th","Pa","U" /

        DO i = 1, nelement
            ! Read from pos.con
            uniq_elements(i) = elementArray(uniqueNrs(i))
        END DO

        DO i = 1, nAtoms
            DO j = 1, nelement
                IF (atomicNumbers(i) == uniqueNrs(j)) THEN
                    symbols(i) = j
                END IF
            END DO
        END DO

        !Define common variable for fNN
        total_natoms=nAtoms
        ! Load fNN parameters by reading mlff.pyamff
        CALL loadfNNParas(nelement, uniq_elements)
        CALL init
    END SUBROUTINE

!    SUBROUTINE prepfNN_ase(nelement, MAX_FPS, uniq_elements)
!        ! Read neural network parameters from mlff.pyamff for ASE Fortran calculator.
!        ! This subroutine is also desgined to be called only once.
!        IMPLICIT NONE
!        ! Inputs
!        INTEGER :: nelement, MAX_FPS 
!        CHARACTER*2, DIMENSION(nelement) :: uniq_elements
!  
!        ! Load fNN parameters by reading mlff.pyamff
!        CALL loadfNNParas(nelement, MAX_FPS, uniq_elements)
!        !CALL init
!
!    END SUBROUTINE
    
    SUBROUTINE update_atomInfo(nAtoms, nelement, MAX_FPS, symbols,in_natomsarr) 
        !This subroutine updates atomic information and ASE fortran calculator uses this subroutine to update atomic information
        !of each image. 
        
        ! Inputs
        INTEGER :: natoms, MAX_FPS, nelement
        INTEGER, DIMENSION(natoms) :: symbols
        INTEGER, DIMENSION(nelement), OPTIONAL :: in_natomsarr
        !Variables
        INTEGER :: i, j, k
        !print *, 'total_natoms : line 80 nnmodule.f90 ',total_natoms
        ! Update common variable, total_natoms, for fNN for each image
        total_natoms=nAtoms
        !print *, 'line 83 nnmodule.f90'
        IF (.NOT. ALLOCATED(natoms_arr)) THEN
             ! Allocate natoms_arr and update the array from symbols array for fNN
             ALLOCATE(natoms_arr(nelement))
        END IF
        IF (PRESENT(in_natomsarr)) THEN
            natoms_arr=in_natomsarr
        ELSE
            DO k = 1, nelement
                natoms_arr(k) = COUNT(symbols .EQ. k)
            END DO
            ! Allocate atom_idx (atomic indices)
        END IF
        IF (.NOT. ALLOCATED(atom_idx)) THEN
            ALLOCATE(atom_idx(MAXVAL(natoms_arr),nelement))
        END IF
        IF (.NOT. ALLOCATED(magnitude)) THEN !if magnitude is not allocated, interceptscale won't be either ; should be able to use only 1 if statement
            ! Allocate magnitudes and interceptScale array (for Normalization) here
            ALLOCATE(magnitude(MAX_FPS, MAXVAL(natoms_arr), nelement))
            ALLOCATE(interceptScale(MAX_FPS, MAXVAL(natoms_arr), nelement))
        END IF 
        magnitude = 0.0
        interceptScale = 0.0
        ! Store magnitude interceptScale in memory
        CALL normalizeParas(nelement)
        
        ! Prepare weights, biases with updated atomic information 
        ! This process is create weights, bias arrays for matmul computations
        CALL init
    !print *, 'END of SUBROUTINE UpdateAtomInfo'
    END SUBROUTINE
      
    ! SUBROUTINE update_atomInfo2(nAtoms, nelement, MAX_FPS, in_natomsarr)
    !     !This subroutine is used for training module
    !     IMPLICIT NONE
    !     !Inputs
    !     INTEGER :: nAtoms, nelement, MAX_FPS
    !     INTEGER, DIMENSION(nelement) :: in_natomsarr
    !     !Variables
    !     INTEGER :: i, k

    !     ! Update common variable, total_natoms, for fNN for each image
    !     total_natoms=nAtoms
    !     ! Allocate natoms_arr and update the array from input
    !     IF (ALLOCATED(natoms_arr) .EQV. .FALSE.) THEN 
    !       ALLOCATE(natoms_arr(nelement))
    !     END IF
    !     natoms_arr=in_natomsarr

    !     ! Allocate magnitudes and interceptScale array (for Normalization) here
    !     IF (ALLOCATED(magnitude) .EQV. .FALSE.) THEN
    !       ALLOCATE(magnitude(MAX_FPS, MAXVAL(natoms_arr), nelement))
    !       ALLOCATE(interceptScale(MAX_FPS, MAXVAL(natoms_arr), nelement))
    !     END IF
    !     !Load magnitude and intercept scale from elementwise magnitude and intercept scale
    !     magnitude = 0.0
    !     interceptScale = 0.0
    !     DO i=1, nelement
    !         DO k=1, natoms_arr(i)
    !           magnitude(1:nGs(i),k,i) = elem_mag(1:nGs(i),i)
    !           interceptScale(1:nGs(i),k,i) = elem_intercept(1:nGs(i),i)
    !         END DO
    !     END DO    
    !     ! Prepare weights, biases with updated atomic information 
    !     ! This process is create weights, bias arrays for matmul computations
    !     CALL init
 

    ! END SUBROUTINE
    
    !SUBROUTINE forwardCalc(num_neigh, max_num_neigh, neighs, symbols, &
    !                       fps, dfps, max_nGs, max_natoms_arr, max_hidneurons)
    SUBROUTINE forwardCalc(num_neigh, MAX_NEIGHS, neighs, &
                           max_nGs, max_natoms_arr, max_hidneurons)
                           !fps, dfps, max_nGs, max_natoms_arr, max_hidneurons
        IMPLICIT NONE
        !f2py INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(aux) :: nGs
        !f2py INTEGER, INTENT(aux) :: total_natoms
        INTEGER :: MAX_NEIGHS, max_nGs, max_natoms_arr, max_hidneurons
        !INTEGER, DIMENSION(total_natoms) :: num_neigh, symbols
        INTEGER, DIMENSION(total_natoms) :: num_neigh
        !INTEGER, DIMENSION(total_natoms,total_natoms) :: neighs
        INTEGER, DIMENSION(total_natoms, MAX_NEIGHS) :: neighs
        !DOUBLE PRECISION, DIMENSION(total_natoms, max_nGs) :: fps
        !DOUBLE PRECISION, DIMENSION(total_natoms, MAX_NEIGHS, 3, max_nGs) :: dfps
        !Variables
        INTEGER :: i, j, k
        DOUBLE PRECISION, DIMENSION(max_natoms_arr,max_nGs,nelements) :: ordered_fps
       
        ! Reorder calculated fps
        DO i = 1, nelements
            DO j = 1, natoms_arr(i)
                ordered_fps(j,1:nGs(i),i) = fps(atom_idx(j,i),1:nGs(i))
            END DO
            !print *, 'ordered_fps'
            !print *, ordered_fps(1:natoms_arr(i), 1:nGs(i),i)
        END DO
        !CALL forward(num_neigh, max_num_neigh, neighs, symbols, ordered_fps, dfps, max_nGs, max_natoms_arr, max_hidneurons)
        CALL forward(num_neigh, MAX_NEIGHS, neighs, ordered_fps, max_nGs, max_natoms_arr, max_hidneurons)
    
    END SUBROUTINE

END MODULE
