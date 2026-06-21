MODULE fpCalc
    USE nlist
    USE fbp
    USE fpType
    USE atomsProp
    USE normalize, only: normalizeParas, loadnormalizeParas
    USE nnType !, only: natoms_arr, nGs, fpminvs, fpmaxvs, diffs, magnitude, interceptScale, coheEs, use_cohesive_energy, atom_idx
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: read_fpParas, read_mlffParas, calcg2s_n, calcfps, cleanup, read_mlff, calcg1s, &
              atomsCleanup, cleanup_ase, clean_mlff, &
              update_mlff_model, update_mlff, write_mlff
    TYPE (fingerprints),DIMENSION(:),ALLOCATABLE,SAVE :: fpParas
    DOUBLE PRECISION, SAVE :: max_rcut   !fetch rcut from fpParas
    CHARACTER*2, DIMENSION(:), ALLOCATABLE :: uniq_elements
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: rmins
    INTEGER :: nG1, nG2    
    INTEGER, PUBLIC :: max_fps, maxneighs
    CONTAINS

!!!!!!-----------------------------------------------------------------------------------!
!!!!!! loopfps: calculate fingerprints for each images based on type of fingerprints
!!!!!!           i.e., G1(H, H)
!!!!!!          symbols: [0, 1, 1, 1]
!!!!!!-----------------------------------------------------------------------------------!
    SUBROUTINE read_fpParas(filename, nelement, coeh)

        ! inputs
        CHARACTER*100 :: filename
        INTEGER :: nelement
!f2py   INTENT(IN) :: filename
!f2py   INTENT(IN) :: nelement
        CHARACTER*2 :: G_type
        CHARACTER*2 :: center, neigh1, neigh2
        INTEGER*4 :: i, j, k, accN

        ! variables
        INTEGER :: cidx, tempidx, nidx1, nidx2 !nG1, nG2
        INTEGER :: currIndex
        INTEGER :: nFPs
        DOUBLE PRECISION :: eta, junk, Rs, zeta, lambda, thetas, rcut
        CHARACTER*10 :: fp_type
        !CHARACTER*2, DIMENSION(nelement) :: uniq_elements

        ! outputs
        DOUBLE PRECISION, DIMENSION(nelement) :: coeh
!f2py   INTENT(OUT) :: coeh
        !COMMON fpParas, max_rcut

        ! Read fp Parameters from fpParas.dat
        !print *, 'reading', filename
        IF(ALLOCATED(rmins)) THEN
          DEALLOCATE(rmins)
        END IF
        IF(ALLOCATED(uniq_elements)) THEN
          DEALLOCATE(uniq_elements)
        END IF
        ALLOCATE(uniq_elements(nelement))
        max_rcut = 0.d0
        !Initialize rmins
        ALLOCATE(rmins(nelement, nelement))
        DO i = 1, nelement
            DO j = 1, nelement
                rmins(i, j) = 100.0
                !print*, nelement, i, j
            END DO
        END DO

        OPEN (12, FILE = filename, STATUS = 'old')
        READ (12,*)
        READ (12,*) fp_type
        READ (12,*)
        READ (12,*) uniq_elements
        READ (12,*) coeh
        READ (12,*)
        READ (12,*) nG1, nG2
        !print *, nG1, nG2
        IF(ALLOCATED(fpParas)) THEN
          DEALLOCATE(fpParas)
        END IF
        ALLOCATE(fpParas(nelement))
        DO i = 1, nelement
            ALLOCATE(fpParas(i)%g1s(nelement))
            ALLOCATE(fpParas(i)%g2s(nelement, nelement))
            fpParas(i)%tnFPs = 0
            fpParas(i)%g1_startpoint = 0
            fpParas(i)%g1_endpoint = 0
            fpParas(i)%g2_startpoint = 0
            fpParas(i)%g2_endpoint = 0
            DO j = 1, nelement
                fpParas(i)%g1s(j)%startpoint = 0
                fpParas(i)%g1s(j)%endpoint = 0
                fpParas(i)%g1s(j)%nFPs = 0
                fpParas(i)%g1s(j)%currIndex = 0
                DO k = 1, nelement
                    fpParas(i)%g2s(j, k)%startpoint = 0
                    fpParas(i)%g2s(j, k)%endpoint = 0
                    fpParas(i)%g2s(j, k)%nFPs = 0
                    fpParas(i)%g2s(j, k)%currIndex = 0
                END DO
            END DO
        END DO
        !print *, 'allocated'
        IF (nG1 .GT. 0) THEN
            READ (12,*)
            DO i = 1, nG1
                READ (12,*) G_type, center, neigh1, junk, junk, junk
                cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                nidx1 = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                !print *, 'nG1 index', center, neigh1, cidx, nidx1
                fpParas(cidx)%g1s(nidx1)%nFPs = fpParas(cidx)%g1s(nidx1)%nFPs + 1
            END DO
            !print *,'nG1', fpParas(1)%g1s(1)%nFPs
        END IF

        ! Count nG2 if there are G2s
        IF (nG2 .GT. 0) THEN
            READ (12,*)
            DO i = 1, nG2
                READ (12,*) G_type, center, neigh1, neigh2, junk, junk, junk, junk, junk
                cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                nidx1 = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                nidx2 = find_loc_char(uniq_elements, neigh2, SIZE(uniq_elements))
                fpParas(cidx)%g2s(nidx1, nidx2)%nFPs = fpParas(cidx)%g2s(nidx1, nidx2)%nFPs + 1
                IF (nidx1 .NE. nidx2) THEN
                    fpParas(cidx)%g2s(nidx2, nidx1)%nFPs = fpParas(cidx)%g2s(nidx2, nidx1)%nFPs + 1
                END IF
            END DO
            !print *, fpParas(1)%g2s(1,1)%nFPs
        END IF

        ! ALLOCATE fpParas
        DO i = 1, nelement
            DO j = 1, nelement
                nFPs = fpParas(i)%g1s(j)%nFPs
                ALLOCATE(fpParas(i)%g1s(j)%etas(nFPs))
                ALLOCATE(fpParas(i)%g1s(j)%rss(nFPs))
                ALLOCATE(fpParas(i)%g1s(j)%r_cuts(nFPs))
                DO k = 1, nelement
                    nFPs = fpParas(i)%g2s(j,k)%nFPs
                    IF (nFPs .GT. 0) THEN
                        ALLOCATE(fpParas(i)%g2s(j,k)%etas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%gammas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%lambdas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%zetas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%r_cuts(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%theta_ss(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%g2(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%dg2(3,3,nFPs)) ! atoms w.r.t * xyz * nFPs
                    END IF
                END DO
            END DO
        END DO

        ! Rewind the fpParas.dat file
        REWIND 12
        ! read fp Paras dat again and put data into the right place 
        READ (12,*) !skip first line
        READ (12,*) !skip second line
        READ (12,*) !skip third line
        READ (12,*) !skip fourth line
        READ (12,*) !skip fifth line
        READ (12,*) !skip sixth line
        READ (12,*) !skip seventh line
        IF (nG1 .GT. 0) THEN
            READ (12,*)
            DO i = 1, nG1
                READ (12,*) G_type, center, neigh1, eta, Rs, rcut
                ! ideally have dynamic rcut for different interactions
                IF (rcut .GT. max_rcut) THEN
                    IF (max_rcut > 0.d0) THEN
                        PRINT *, 'Error: multiple rcuts are not supported yet. Use a single rcut value'
                        STOP
                    ELSE
                    max_rcut = rcut
                    END IF
                ENDIF
                !print *, G_type, center, neigh1, eta, Rs, rcut
                cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                nidx1 = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))

                fpParas(cidx)%g1s(nidx1)%fp_type = G_type
                fpParas(cidx)%g1s(nidx1)%species1 = center
                fpParas(cidx)%g1s(nidx1)%species1_code = cidx
                fpParas(cidx)%g1s(nidx1)%species2 = neigh1
                fpParas(cidx)%g1s(nidx1)%currIndex = fpParas(cidx)%g1s(nidx1)%currIndex + 1

                currIndex = fpParas(cidx)%g1s(nidx1)%currIndex
                fpParas(cidx)%g1s(nidx1)%etas(currIndex) = eta
                fpParas(cidx)%g1s(nidx1)%rss(currIndex) = Rs
                fpParas(cidx)%g1s(nidx1)%r_cuts(currIndex) = rcut
                fpParas(cidx)%tnFPs = fpParas(cidx)%tnFPs + 1

                IF (nidx1 .EQ. 1) THEN
                    fpParas(cidx)%g1s(nidx1)%startpoint = 1
                    fpParas(cidx)%g1s(nidx1)%endpoint = fpParas(cidx)%g1s(nidx1)%endpoint + 1
                ELSE
                    fpParas(cidx)%g1s(nidx1)%startpoint = fpParas(cidx)%g1s(nidx1-1)%endpoint + 1
                    fpParas(cidx)%g1s(nidx1)%endpoint = fpParas(cidx)%g1s(nidx1)%startpoint + fpParas(cidx)%g1s(nidx1)%nFPs - 1
                END IF
            END DO
        END IF
        !print *, 'nG1 reading done'
        IF (nG2 .GT. 0) THEN
            READ (12,*) !skip #type
            DO k = 1, nelement
                fpParas(k)%g1_endpoint = fpParas(k)%tnFPs
                fpParas(k)%g2_startpoint = fpParas(k)%g1_endpoint + 1
            END DO

            DO i = 1, nG2
                READ (12,*) G_type, center, neigh1, neigh2, eta, zeta, lambda, thetas, rcut

                ! ideally have dynamic rcut for different interactions
                IF (rcut .GT. max_rcut) THEN
                    IF (max_rcut > 0.d0) THEN
                        PRINT *, 'Error: multiple rcuts are not supported yet. Use a single rcut value'
                        STOP
                    ELSE
                        max_rcut = rcut
                    END IF
                END IF

  !!!!!!!        cidx = FINDLOC(uniq_elements, VALUE=center, DIM=1)
                cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                nidx1 = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                nidx2 = find_loc_char(uniq_elements, neigh2, SIZE(uniq_elements))
                !print *, 'nG2 index', center, neigh1, neigh2, cidx, nidx1, nidx2
                IF (nidx1 .GT. nidx2) THEN
                    tempidx = nidx1
                    nidx1 = nidx2
                    nidx2 = tempidx
                END IF
                fpParas(cidx)%g2s(nidx1, nidx2)%fp_type = G_type
                fpParas(cidx)%g2s(nidx1, nidx2)%species1 = center
                fpParas(cidx)%g2s(nidx1, nidx2)%species1_code = cidx
                fpParas(cidx)%g2s(nidx1, nidx2)%species2 = neigh1
                fpParas(cidx)%g2s(nidx1, nidx2)%currIndex = fpParas(cidx)%g2s(nidx1, nidx2)%currIndex + 1
                currIndex = fpParas(cidx)%g2s(nidx1, nidx2)%currIndex

                !print *, 'currIndex', currIndex
                fpParas(cidx)%g2s(nidx1, nidx2)%etas(currIndex)     = eta
                fpParas(cidx)%g2s(nidx1, nidx2)%zetas(currIndex)    = zeta
                fpParas(cidx)%g2s(nidx1, nidx2)%lambdas(currIndex)  = lambda
                fpParas(cidx)%g2s(nidx1, nidx2)%theta_ss(currIndex) = thetas
                fpParas(cidx)%g2s(nidx1, nidx2)%r_cuts(currIndex)   = rcut
                fpParas(cidx)%tnFPs = fpParas(cidx)%tnFPs + 1

                !fpParas(cidx)%g2s(nidx1, nidx2)%nFPs = fpParas(cidx)%g2s(nidx1, nidx2)%nFPs + 1
            END DO
            !print *, 'half done'

            DO i = 1, nelement
                DO j = 1, nelement
                    DO k = j, nelement
                       IF (fpParas(i)%g2s(j, k)%nFPs .EQ. 0) THEN
                            accN = 0
                        ELSE
                            accN = 1
                        END IF
                        IF ((j .EQ. 1) .AND. (k .EQ. 1)) THEN
                            fpParas(i)%g2s(j, k)%startpoint = fpParas(i)%g1_endpoint + accN
                        ELSE
                            IF (j .EQ. k) THEN
                                fpParas(i)%g2s(j, k)%startpoint = fpParas(i)%g2s(j-1, nelement)%endpoint + accN
                            ELSE
                                fpParas(i)%g2s(j, k)%startpoint = fpParas(i)%g2s(j, k-1)%endpoint + accN
                            END IF
                        END IF
                        IF (accN .EQ. 0) THEN
                            fpParas(i)%g2s(j, k)%endpoint = fpParas(i)%g2s(j, k)%startpoint
                        ELSE
                            fpParas(i)%g2s(j, k)%endpoint = fpParas(i)%g2s(j, k)%startpoint + fpParas(i)%g2s(j, k)%nFPs - 1
                        END IF
                    END DO
                END DO
            END DO
            !print *, fpParas(2)%g2s(2,2)%startpoint, fpParas(2)%g2s(2,2)%endpoint
            DO i = 1, nelement
                DO j = 1, nelement
                    DO k = 1, j-1
                        !print *, i, j, k
                        fpParas(i)%g2s(j, k)%etas       =  fpParas(i)%g2s(k, j)%etas
                        fpParas(i)%g2s(j, k)%zetas      =  fpParas(i)%g2s(k, j)%zetas
                        fpParas(i)%g2s(j, k)%lambdas    =  fpParas(i)%g2s(k, j)%lambdas
                        fpParas(i)%g2s(j, k)%theta_ss   =  fpParas(i)%g2s(k, j)%theta_ss
                        fpParas(i)%g2s(j, k)%r_cuts     =  fpParas(i)%g2s(k, j)%r_cuts
                        fpParas(i)%g2s(j, k)%startpoint =  fpParas(i)%g2s(k, j)%startpoint
                        fpParas(i)%g2s(j, k)%endpoint   =  fpParas(i)%g2s(k, j)%endpoint
                    END DO
                END DO
            END DO
            !print*, 'read fpparas'
            !i=1
            !DO j = 1, nelement
            !    DO k = 1, nelement
            !        print *, i, j, k
            !        print*,fpParas(i)%g2s(j, k)%etas       
            !        print*,fpParas(i)%g2s(j, k)%zetas      
            !        print*,fpParas(i)%g2s(j, k)%lambdas    
            !        print*,fpParas(i)%g2s(j, k)%theta_ss   
            !        print*,fpParas(i)%g2s(j, k)%r_cuts     
            !        print*,fpParas(i)%g2s(j, k)%startpoint
            !        print*,fpParas(i)%g2s(j, k)%endpoint   
            !    END DO
            !END DO
            !write(*,*) 'Parameter reading done'
        END IF
        !print*, 'maxcutoff',max_cutoff(1), coeh
        CLOSE(12)
        !
        !Read mlff.pyamff to obtain fpminvs, and fpmaxvs
        !CALL loadnormalizeParas(nAtoms, nelement, MAX_FPS, symbols, uniq_element)
        !ALLOCATE(magnitude(MAX_FPS, MAXVAL(natoms_arr), nelement))
        !ALLOCATE(interceptScale(MAX_FPS, MAXVAL(natoms_arr), nelement))
        !Store magnitude interceptScale in memory
        !CALL normalizeParas(nelement)

    END SUBROUTINE

    !SUBROUTINE read_mlffParas(nAtoms, nelement, MAX_FPS, atomicNumbers, uniqueNrs) BIND(C,name='read_mlffParas')
    SUBROUTINE read_mlffParas(nAtoms, nelement, atomicNumbers, uniqueNrs) BIND(C,name='read_mlffParas')
        !This subroutine is for EON-PyAMFF interface.
        USE, INTRINSIC :: iso_c_binding
        ! Inputs
        INTEGER(c_long) :: nAtoms
        INTEGER(c_int) :: nelement
        INTEGER(c_int), DIMENSION(nAtoms) :: atomicNumbers
        INTEGER(c_int), DIMENSION(nelement) :: uniqueNrs
        INTEGER, DIMENSION(nAtoms) :: symbols
        !CHARACTER*2, DIMENSION(nelement) :: uniq_element

        ! Variables
        INTEGER :: nFPs, cidx, nidx, nidx1, nidx2, tempidx
        INTEGER :: i, j, k, accN, currIndex, m, n
        INTEGER, DIMENSION(nelement) :: fprange_idx 
        DOUBLE PRECISION :: djunk, fpmin, fpmax, eta, Rs, rcut, thetas, zeta, lambda
        CHARACTER*3 :: center, neigh1, neigh2, cjunk
        CHARACTER(LEN=30) :: G_type, line
        CHARACTER*3, DIMENSION(92) :: elementArray
        INTEGER, DIMENSION(:), ALLOCATABLE :: numGs, tnumGs

        INTEGER:: myid, ios
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: flatten_inweights,flatten_hidweights
        CHARACTER(LEN=30) :: model_type
        CHARACTER*2 :: atom_type
        LOGICAL :: haveV1 = .FALSE.

        DATA elementArray / "H","He","Li","Be","B","C","N","O", &
                 "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc", &
                 "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se", &
                 "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag", &
                 "Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd", &
                 "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta", &
                 "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", &
                 "Fr","Ra","Ac","Th","Pa","U" /

        CHARACTER(len=2), DIMENSION(nelement) :: coheElement
        DOUBLE PRECISION, DIMENSION(nelement) :: temp_coheEs
        
        ALLOCATE(coheEs(nelement))
        use_cohesive_energy = .FALSE.
        max_rcut=0.d0
        !filename = 'fpParas.dat'
        !print*, 'allocating rmins'
        nG1 = 0
        nG2 = 0
        ALLOCATE(rmins(nelement, nelement))
        rmins=0.0

        DO i = 1, nAtoms
            DO j = 1, nelement
                IF (atomicNumbers(i) == uniqueNrs(j)) THEN
                    symbols(i) = j
                END IF
            END DO
        END DO

        IF(ALLOCATED(numGs)) THEN 
          DEALLOCATE(numGs)
        END IF
        IF(ALLOCATED(tnumGs)) THEN 
          DEALLOCATE(tnumGs)
        END IF
        ALLOCATE(numGs(nelement))
        ALLOCATE(tnumGs(nelement))
        numGs = 0
        tnumGs = 0
        ! Outputs
        ALLOCATE (uniq_elements(nelement))
        ALLOCATE(natoms_arr(nelement))
        ALLOCATE(nGs(nelement))
        
        nGs = 0
        fprange_idx = 0

        !Define common variables for fNN
        DO k = 1, nelement
           natoms_arr(k) = COUNT(symbols .EQ. k)
        END DO

        ALLOCATE(fpParas(nelement))
        DO i = 1, nelement
            ALLOCATE(fpParas(i)%g1s(nelement))
            ALLOCATE(fpParas(i)%g2s(nelement, nelement))
            fpParas(i)%tnFPs = 0
            fpParas(i)%g1_startpoint = 0
            fpParas(i)%g1_endpoint = 0
            fpParas(i)%g2_startpoint = 0
            fpParas(i)%g2_endpoint = 0
            DO j = 1, nelement
                fpParas(i)%g1s(j)%startpoint = 0
                fpParas(i)%g1s(j)%endpoint = 0
                fpParas(i)%g1s(j)%nFPs = 0
                fpParas(i)%g1s(j)%currIndex = 0
                DO k = 1, nelement
                    fpParas(i)%g2s(j, k)%startpoint = 0
                    fpParas(i)%g2s(j, k)%endpoint = 0
                    fpParas(i)%g2s(j, k)%nFPs = 0
                    fpParas(i)%g2s(j, k)%currIndex = 0
                END DO
            END DO
        END DO
        WRITE(*,*) 'Reading fingerprint parameters from THE mlff.pyamff'
        OPEN (11, FILE='mlff.pyamff', status='old')
        !print*, 'reading mlparas'
        READ (11,*) !skip #Fingerprint type
        READ (11,*) !fp_type
        READ (11,*) !Rmins
        READ (11,*) uniq_elements
        !print*, 'reading Rmins', uniq_elements
        DO i = 1, nelement
            READ (11,*) (rmins(i,j), j = 1, nelement)
        END DO

        !!!!!! If there is a #Cohesive tag
        !print*, 'Counting G1s'
        DO WHILE (line .NE. "#MachineLearning")
            READ (11,*) line! can be #  type or #Cohesive 
            IF (line .EQ. "#Cohesive") THEN
                use_cohesive_energy = .TRUE.
                READ (11,*) coheElement
                READ (11,*) temp_coheEs
                DO m = 1, nelement
                    n = find_loc_char(coheElement, uniq_elements(m), nelement)
                    coheEs(m) = temp_coheEs(n)
                END DO
                READ (11,*) line! skip # type 
            END IF
            IF (line .EQ. "#MachineLearning") GOTO 40
            READ (11,*,IOSTAT=ios) G_type, numGs
            IF (ios .LT. 0) THEN
                WRITE(*,*) "End-of-file reached while reading 'mlff.pyamff'"
                STOP
            ELSE IF (ios .EQ. 5010) THEN
                ! This is the error code corresponding to "Bad real number in item 3 of list input"
                ! It means that we have 1 number representing the total number instead of a number per element
                IF(.NOT. haveV1) WRITE(*,*) "WARNING: V1 of mlff.pyamff detected. A new mlff.pyamff will be generated."
                haveV1 = .TRUE.
            END IF
            IF (G_type .EQ. 'G1') THEN
                ! nG1 = numGs
                IF (.NOT. haveV1) THEN
                    nG1 = SUM(numGs)
                    READ (11,*) ! skip # center neighbor ...
                ELSE
                    nG1 = numGs(1)
                END IF
                DO i = 1, nG1
                    READ (11,*) center, neigh1, djunk, djunk, djunk, djunk, djunk
                    cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                    nidx = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                    nGs(cidx) = nGs(cidx) + 1
                    fpParas(cidx)%g1s(nidx)%nFPs = fpParas(cidx)%g1s(nidx)%nFPs + 1
                END DO
                IF (haveV1) numGs = nGs
            ELSE IF (G_type .EQ. 'G2') THEN
                ! nG2 = numGs
                IF (.NOT. haveV1) THEN
                    nG2 = SUM(numGs)
                    READ (11,*) ! skip # center neighbor ...
                ELSE
                    nG2 = numGs(1)
                END IF
                DO i = 1, nG2
                    READ (11,*) center, neigh1, neigh2, djunk, djunk, djunk, djunk, djunk, djunk, djunk
                    cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                    nidx1 = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                    nidx2 = find_loc_char(uniq_elements, neigh2, SIZE(uniq_elements))
                    nGs(cidx) = nGs(cidx) + 1
                    fpParas(cidx)%g2s(nidx1, nidx2)%nFPs = fpParas(cidx)%g2s(nidx1, nidx2)%nFPs + 1
                    IF (nidx1 .NE. nidx2) THEN
                        fpParas(cidx)%g2s(nidx2, nidx1)%nFPs = fpParas(cidx)%g2s(nidx2, nidx1)%nFPs + 1
                    END IF
                END DO
                IF (haveV1) numGs = nGs
            END IF
            tnumGs = tnumGs + numGs
        END DO

  40    REWIND 11
        !print*, 'allocating fpParas'
        MAX_FPs = MAXVAL(tnumGs)
        ALLOCATE(fpminvs(MAX_FPS, nelement))
        ALLOCATE(fpmaxvs(MAX_FPS, nelement))
        ALLOCATE(diffs(MAX_FPS, nelement))
        fpminvs = 0
        fpmaxvs = 0
        diffs = 0
        ! ALLOCATE fpParas
        DO i = 1, nelement
            DO j = 1, nelement
                nFPs = fpParas(i)%g1s(j)%nFPs
                ALLOCATE(fpParas(i)%g1s(j)%etas(nFPs))
                ALLOCATE(fpParas(i)%g1s(j)%rss(nFPs))
                ALLOCATE(fpParas(i)%g1s(j)%r_cuts(nFPs))
                DO k = 1, nelement
                    nFPs = fpParas(i)%g2s(j,k)%nFPs
                    IF (nFPs .GT. 0) THEN
                        ALLOCATE(fpParas(i)%g2s(j,k)%etas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%gammas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%lambdas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%zetas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%r_cuts(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%theta_ss(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%g2(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%dg2(3,3,nFPs)) ! atoms w.r.t * 
                        fpParas(i)%g2s(j,k)%etas = 0.0
                        fpParas(i)%g2s(j,k)%gammas = 0.0
                        fpParas(i)%g2s(j,k)%lambdas = 0.0
                        fpParas(i)%g2s(j,k)%zetas = 0.0
                        fpParas(i)%g2s(j,k)%r_cuts = 0.0
                        fpParas(i)%g2s(j,k)%theta_ss = 0.0
                        fpParas(i)%g2s(j,k)%g2 = 0.0
                        fpParas(i)%g2s(j,k)%dg2 = 0.0
                    END IF
                END DO
            END DO
        END DO

        READ (11,*) !skip #Fingerprint type
        READ (11,*) !fp_type 
        READ (11,*) !Rmins
        READ (11,*) !element type
        DO i = 1, nelement
            READ (11,*)   !Rmins
        END DO

         ! Mai edit: Read the cohesive Energy
         !!!!!!!!! If there is a cohesive tag in mlff.pyamff
         ! read the cohesive energies into the array

        READ (11,*) line!  can be #  type or #Cohesive
        IF (line .EQ. "#Cohesive") THEN
            READ (11,*)
            READ (11,*)
            READ (11,*) ! # type
        END IF
        IF (.NOT. haveV1) THEN
            READ (11,*) G_type, numGs
        ELSE
            READ(11,*) G_type, numGs(1)
        END IF
        IF (nG1 .GT. 0) THEN
            READ (11,*) !skip center ...
            DO i = 1, nG1
                READ (11,*) center, neigh1, eta, Rs, rcut, fpmin, fpmax
                ! ideally have dynamic rcut for different interactions
                IF (rcut .GT. max_rcut) THEN
                    IF (max_rcut > 0.d0) THEN
                        print *, 'Error: multiple rcuts are not supported yet. Use a single rcut value'
                        STOP
                    ELSE    
                        max_rcut = rcut
                    END IF
                END IF
                cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                nidx1 = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                fpParas(cidx)%g1s(nidx1)%fp_type = G_type
                fpParas(cidx)%g1s(nidx1)%species1 = center
                fpParas(cidx)%g1s(nidx1)%species1_code = cidx
                fpParas(cidx)%g1s(nidx1)%species2 = neigh1
                fpParas(cidx)%g1s(nidx1)%currIndex = fpParas(cidx)%g1s(nidx1)%currIndex + 1
                currIndex = fpParas(cidx)%g1s(nidx1)%currIndex
                fpParas(cidx)%g1s(nidx1)%etas(currIndex) = eta
                fpParas(cidx)%g1s(nidx1)%rss(currIndex) = Rs
                fpParas(cidx)%g1s(nidx1)%r_cuts(currIndex) = rcut
                fpParas(cidx)%tnFPs = fpParas(cidx)%tnFPs + 1
                
                fpminvs(1+fprange_idx(cidx),cidx) = fpmin
                fpmaxvs(1+fprange_idx(cidx),cidx) = fpmax
                fprange_idx(cidx)=fprange_idx(cidx)+1 
                
                IF (nidx1 .EQ. 1) THEN
                    fpParas(cidx)%g1s(nidx1)%startpoint = 1
                    fpParas(cidx)%g1s(nidx1)%endpoint = fpParas(cidx)%g1s(nidx1)%endpoint + 1
                ELSE
                    fpParas(cidx)%g1s(nidx1)%startpoint = fpParas(cidx)%g1s(nidx1-1)%endpoint + 1
                    fpParas(cidx)%g1s(nidx1)%endpoint = fpParas(cidx)%g1s(nidx1)%startpoint + fpParas(cidx)%g1s(nidx1)%nFPs - 1
                END IF
            END DO
        END IF
        !print*, 'reading G2s'
        IF (nG2 .GT. 0) THEN
            IF (nG1 .GT. 0) THEN
                READ (11,*) !skip type ...
                IF (.NOT. haveV1) THEN
                    READ (11,*) G_type, numGs
                ELSE
                    READ (11,*) G_type, numGs(1)
                END IF
            ELSE
                CONTINUE
            END IF    
            !print *, G_type, numGs
            DO k = 1, nelement
                fpParas(k)%g1_endpoint = fpParas(k)%tnFPs
                fpParas(k)%g2_startpoint = fpParas(k)%g1_endpoint + 1
            END DO
            READ (11,*) !skip # center...
            DO i = 1, nG2
                READ (11,*) center, neigh1, neigh2, eta, zeta, lambda, thetas, rcut, fpmin, fpmax
                ! ideally have dynamic rcut for different interactions
                IF (rcut .GT. max_rcut) THEN
                    IF (max_rcut > 0.d0) THEN
                        print *, 'Error: multiple rcuts are not supported yet. Use a single rcut value'
                        STOP
                    ELSE
                        max_rcut = rcut
                    END IF
                END IF
                cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                nidx1 = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                nidx2 = find_loc_char(uniq_elements, neigh2, SIZE(uniq_elements))
                IF (nidx1 .GT. nidx2) THEN
                    tempidx = nidx1
                    nidx1 = nidx2
                    nidx2 = tempidx
                END IF
                fpParas(cidx)%g2s(nidx1, nidx2)%fp_type = G_type
                fpParas(cidx)%g2s(nidx1, nidx2)%species1 = center
                fpParas(cidx)%g2s(nidx1, nidx2)%species1_code = cidx
                fpParas(cidx)%g2s(nidx1, nidx2)%species2 = neigh1
                fpParas(cidx)%g2s(nidx1, nidx2)%currIndex = fpParas(cidx)%g2s(nidx1, nidx2)%currIndex + 1
                currIndex = fpParas(cidx)%g2s(nidx1, nidx2)%currIndex
 
                fpParas(cidx)%g2s(nidx1, nidx2)%etas(currIndex)     = eta
                fpParas(cidx)%g2s(nidx1, nidx2)%zetas(currIndex)    = zeta
                fpParas(cidx)%g2s(nidx1, nidx2)%lambdas(currIndex)  = lambda
                fpParas(cidx)%g2s(nidx1, nidx2)%theta_ss(currIndex) = thetas
                fpParas(cidx)%g2s(nidx1, nidx2)%r_cuts(currIndex)   = rcut
                fpParas(cidx)%tnFPs = fpParas(cidx)%tnFPs + 1
                
                fpminvs(1+fprange_idx(cidx),cidx) = fpmin
                fpmaxvs(1+fprange_idx(cidx),cidx) = fpmax
                fprange_idx(cidx)=fprange_idx(cidx)+1
            END DO

            DO i = 1, nelement
                DO j = 1, nelement
                    DO k = j, nelement
                        IF (fpParas(i)%g2s(j, k)%nFPs .EQ. 0) THEN
                            accN = 0
                        ELSE
                            accN = 1
                        END IF

                        IF ((j .EQ. 1) .AND. (k .EQ. 1)) THEN
                            fpParas(i)%g2s(j, k)%startpoint = fpParas(i)%g1_endpoint + accN
                        ELSE
                            IF (j .EQ. k) THEN
                                fpParas(i)%g2s(j, k)%startpoint = fpParas(i)%g2s(j-1, nelement)%endpoint+accN
                            ELSE
                                fpParas(i)%g2s(j, k)%startpoint = fpParas(i)%g2s(j, k-1)%endpoint+accN
                            END IF
                        END IF
                        IF (accN .EQ. 0) THEN
                            fpParas(i)%g2s(j, k)%endpoint = fpParas(i)%g2s(j, k)%startpoint
                        ELSE
                            fpParas(i)%g2s(j, k)%endpoint = fpParas(i)%g2s(j, k)%startpoint + fpParas(i)%g2s(j, k)%nFPs - 1
                        END IF
                    END DO
                END DO
            END DO
            !print *, fpParas(2)%g2s(2,2)%startpoint, fpParas(2)%g2s(2,2)%endpoint
            !print*, 'G2 reading done'
            DO i = 1, nelement
                DO j = 1, nelement
                    DO k = 1, j - 1
                        !print *, i, j, k
                        IF (fpParas(i)%g2s(j,k)%nFPs .GT. 0) THEN
                            fpParas(i)%g2s(j, k)%etas       =  fpParas(i)%g2s(k, j)%etas
                            fpParas(i)%g2s(j, k)%zetas      =  fpParas(i)%g2s(k, j)%zetas
                            fpParas(i)%g2s(j, k)%lambdas    =  fpParas(i)%g2s(k, j)%lambdas
                            fpParas(i)%g2s(j, k)%theta_ss   =  fpParas(i)%g2s(k, j)%theta_ss
                            fpParas(i)%g2s(j, k)%r_cuts     =  fpParas(i)%g2s(k, j)%r_cuts
                            fpParas(i)%g2s(j, k)%startpoint =  fpParas(i)%g2s(k, j)%startpoint
                            fpParas(i)%g2s(j, k)%endpoint   =  fpParas(i)%g2s(k, j)%endpoint
                        END IF
                    END DO
                END DO
            END DO
        END IF 
                
        !Allocate atom_idx here 
        ALLOCATE(atom_idx(MAXVAL(natoms_arr),nelement))
        !print*, 'Allocating mags'
        ALLOCATE(magnitude(MAX_FPS, MAXVAL(natoms_arr), nelement))
        ALLOCATE(interceptScale(MAX_FPS, MAXVAL(natoms_arr), nelement))
        magnitude = 0
        interceptScale = 0
        atom_idx = 0
        !Store magnitude interceptScale in memory
        CALL normalizeParas(nelement)

        CLOSE (11)
    END SUBROUTINE

    SUBROUTINE read_mlff(nelement, mlff_file)
    !This subroutine is for ASE Fortran calculator and this is designed to be called only once. 
    !Unlike read_mlffParas, this subroutine does not contain atomic information since 
    !each image can have different number of atoms. 
        ! Inputs
        !INTEGER :: nAtoms
        INTEGER :: nelement

        ! Variables
        INTEGER :: nFPs, cidx, nidx, nidx1, nidx2, tempidx
        INTEGER :: i, j, k, accN, currIndex, m, n
        INTEGER, DIMENSION(:), ALLOCATABLE :: numGs, tnumGs
        INTEGER, DIMENSION(nelement) :: fprange_idx 
        DOUBLE PRECISION :: djunk, fpmin, fpmax, eta, Rs, rcut, thetas, zeta, lambda
        CHARACTER*3 :: center, neigh1, neigh2, cjunk
        CHARACTER*20 :: mlff_file
        CHARACTER(LEN=30) :: G_type, line
        ! CHARACTER*3, DIMENSION(92) :: elementArray
        INTEGER :: myid, mlff_file_unit=44
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: flatten_inweights, flatten_hidweights
        CHARACTER(LEN=30) :: model_type,comment
        CHARACTER*2 :: atom_type

        CHARACTER(len=2), DIMENSION(nelement) :: coheElement
        DOUBLE PRECISION, DIMENSION(nelement) :: temp_coheEs
        INTEGER,PARAMETER :: seedval = 6 !!!Dummy starter (Eboni suggested) ... 
        INTEGER :: ios
        LOGICAL :: haveV2 = .FALSE., haveV1 = .FALSE., haveV2_noparams = .FALSE.

        IF(ALLOCATED(coheEs)) THEN
          DEALLOCATE(coheEs)
        END IF
        ALLOCATE(coheEs(nelement))
        !print *, 'line 700 in read_mlff fingerprints.f90'
        use_cohesive_energy = .FALSE.
        max_rcut=0.d0
        !filename = 'fpParas.dat'
        nG1 = 0
        nG2 = 0
        IF(ALLOCATED(rmins)) THEN 
          DEALLOCATE(rmins)
        END IF
        ALLOCATE(rmins(nelement, nelement))
        rmins = 0.0
        !print *, 'line 711 fingerprints.f90'
        IF(ALLOCATED(uniq_elements)) THEN 
          DEALLOCATE(uniq_elements)
        END IF
        IF(ALLOCATED(nGs)) THEN 
          DEALLOCATE(nGs)
        END IF
        IF(ALLOCATED(numGs)) THEN 
          DEALLOCATE(numGs)
        END IF
        IF(ALLOCATED(tnumGs)) THEN 
          DEALLOCATE(tnumGs)
        END IF
        IF(ALLOCATED(fpminvs)) THEN 
          DEALLOCATE(fpminvs)
        END IF
        IF(ALLOCATED(fpmaxvs)) THEN 
          DEALLOCATE(fpmaxvs)
        END IF
        IF(ALLOCATED(diffs)) THEN 
          DEALLOCATE(diffs)
        END IF
        ALLOCATE(uniq_elements(nelement))
        ALLOCATE(nGs(nelement))
        ALLOCATE(numGs(nelement))
        ALLOCATE(tnumGs(nelement))
        nGs = 0
        tnumGs = 0
        fprange_idx=0
        !print *, 'line 740 fingerprints.f90'

        IF(ALLOCATED(fpParas)) THEN 
          DEALLOCATE(fpParas)
        END IF
        ALLOCATE(fpParas(nelement))
        DO i = 1, nelement
            ALLOCATE(fpParas(i)%g1s(nelement))
            ALLOCATE(fpParas(i)%g2s(nelement, nelement))
            fpParas(i)%tnFPs = 0
            fpParas(i)%g1_startpoint = 0
            fpParas(i)%g1_endpoint = 0
            fpParas(i)%g2_startpoint = 0
            fpParas(i)%g2_endpoint = 0
            DO j = 1, nelement
                fpParas(i)%g1s(j)%startpoint = 0
                fpParas(i)%g1s(j)%endpoint = 0
                fpParas(i)%g1s(j)%nFPs = 0
                fpParas(i)%g1s(j)%currIndex = 0
                DO k = 1, nelement
                    fpParas(i)%g2s(j, k)%startpoint = 0
                    fpParas(i)%g2s(j, k)%endpoint = 0
                    fpParas(i)%g2s(j, k)%nFPs = 0
                    fpParas(i)%g2s(j, k)%currIndex = 0
                END DO
            END DO
        END DO
        !print *, 'line 767 fingerprints.f90'
        !OPEN (11, FILE='mlff.pyamff', status='old')
        OPEN (newunit=mlff_file_unit, FILE=mlff_file, status='old',action='read',iostat=ios)
        IF (ios .NE. 0) THEN
            WRITE(*,*) "Error opening file: ", TRIM(mlff_file)
            STOP
        END IF
        READ (mlff_file_unit,*) !skip #Fingerprint type
        READ (mlff_file_unit,*) !fp_type 
        READ (mlff_file_unit,*) !Rmins
        READ (mlff_file_unit,*) uniq_elements
        !print*, 'reading Rmins'
        DO i = 1, nelement
            READ (mlff_file_unit,*) (rmins(i,j), j = 1, nelement)
        END DO
        !print *, 'line 778 fingerprints.f90'

        ! If there is a #Cohesive tag
        DO WHILE (line .NE. "#MachineLearning")
            READ (mlff_file_unit,*) line! can be #  type or #Cohesive
            IF (line .EQ. "#Cohesive") THEN
                use_cohesive_energy = .TRUE.
                READ (mlff_file_unit,*) coheElement
                READ (mlff_file_unit,*) temp_coheEs
                DO m = 1, nelement
                    n = find_loc_char(coheElement, uniq_elements(m), nelement)
                    coheEs(m) = temp_coheEs(n)
                END DO
                READ (mlff_file_unit,*) line ! skip # type 
            END IF
            IF (line .EQ. "#MachineLearning") GOTO 40
            READ (mlff_file_unit,*,IOSTAT=ios) G_type, numGs
            ! If ios == 0, then things are okay
            IF (ios .LT. 0) THEN
                WRITE(*,*) "End-of-file reached while reading " // TRIM(mlff_file)
                STOP
            ELSE IF (ios .EQ. 5010) THEN 
                ! This is the error code corresponding to "Bad real number in item 3 of list input"
                ! It means that we have 1 number representing the total number instead of a number per element
                IF(.NOT. haveV1) WRITE(*,*) "WARNING: V1 of mlff.pyamff detected. A new mlff.pyamff will be generated."
                haveV1 = .TRUE.
            END IF
            IF (G_type .EQ. 'G1') THEN
                IF (.NOT. haveV1) THEN
                    nG1 = SUM(numGs)
                    READ (mlff_file_unit,*) ! skip # center neighbor ...
                ELSE
                    nG1 = numGs(1)
                END IF
                DO i = 1, nG1
                    READ (mlff_file_unit,*) center, neigh1, djunk, djunk, djunk, djunk, djunk
                    cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                    nidx = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                    nGs(cidx) = nGs(cidx) + 1
                    fpParas(cidx)%g1s(nidx)%nFPs = fpParas(cidx)%g1s(nidx)%nFPs + 1
                END DO
                IF (haveV1) numGs = nGs
            !print *, 'line 806 fingerprints.f90'

            ELSE IF (G_type .EQ. 'G2') THEN
                IF (.NOT. haveV1) THEN
                    nG2 = SUM(numGs)
                    READ (mlff_file_unit,*) ! skip # center neighbor ...
                ELSE
                    nG2 = numGs(1)
                END IF
                DO i = 1, nG2
                    READ (mlff_file_unit,*) center, neigh1, neigh2, djunk, djunk, djunk, djunk, djunk, djunk, djunk
                    cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                    nidx1 = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                    nidx2 = find_loc_char(uniq_elements, neigh2, SIZE(uniq_elements))
                    nGs(cidx) = nGs(cidx) + 1
                    fpParas(cidx)%g2s(nidx1, nidx2)%nFPs = fpParas(cidx)%g2s(nidx1, nidx2)%nFPs + 1

                    IF (nidx1 .NE. nidx2) THEN
                        fpParas(cidx)%g2s(nidx2, nidx1)%nFPs = fpParas(cidx)%g2s(nidx2, nidx1)%nFPs + 1
                    END IF
                END DO
                IF (haveV1) numGs = nGs
            END IF
            tnumGs = tnumGs + numGs
        END DO
        !print *, 'line 825 fingerprints.f90'

  40    REWIND mlff_file_unit
        MAX_FPs = MAXVAL(tnumGs)
        ALLOCATE(fpminvs(MAX_FPS, nelement))
        ALLOCATE(fpmaxvs(MAX_FPS, nelement))
        ALLOCATE(diffs(MAX_FPS, nelement))
        fpminvs = 0
        fpmaxvs = 0
        diffs = 0
        !print *, 'line 835 fingerprints.f90'

        ! ALLOCATE fpParas
        DO i = 1, nelement
            DO j = 1, nelement
                nFPs = fpParas(i)%g1s(j)%nFPs
                ALLOCATE(fpParas(i)%g1s(j)%etas(nFPs))
                ALLOCATE(fpParas(i)%g1s(j)%rss(nFPs))
                ALLOCATE(fpParas(i)%g1s(j)%r_cuts(nFPs))
                DO k = 1, nelement
                    nFPs = fpParas(i)%g2s(j,k)%nFPs
                    IF (nFPs .GT. 0) THEN
                        ALLOCATE(fpParas(i)%g2s(j,k)%etas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%gammas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%lambdas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%zetas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%r_cuts(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%theta_ss(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%g2(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%dg2(3,3,nFPs)) ! atoms w.r.t * xyz * nFPs 
                        fpParas(i)%g2s(j,k)%etas = 0.0
                        fpParas(i)%g2s(j,k)%gammas = 0.0
                        fpParas(i)%g2s(j,k)%lambdas = 0.0
                        fpParas(i)%g2s(j,k)%zetas = 0.0
                        fpParas(i)%g2s(j,k)%r_cuts = 0.0
                        fpParas(i)%g2s(j,k)%theta_ss = 0.0
                        fpParas(i)%g2s(j,k)%g2 = 0.0
                        fpParas(i)%g2s(j,k)%dg2 = 0.0
                    ENDIF
                    !print*, 'fpParas', fpParas(i)%g2s(j,k)%etas
                END DO
            END DO
        END DO
        !print *, 'line 868 fingerprints.f90'

        READ (mlff_file_unit,*) !skip #Fingerprint type
        READ (mlff_file_unit,*) !fp_type
        READ (mlff_file_unit,*) !Rmins
        READ (mlff_file_unit,*) !element type
        DO i = 1, nelement
            READ (mlff_file_unit,*)  !Rmins
        END DO
        !print *, 'line 877 fingerprints.f90'

         ! Mai edit: Read the cohesive Energy
         !!!!!!!!! If there is a #Cohesive tag
         ! read the cohesive energies to the array

         !print*, 'reading g1'
        READ (mlff_file_unit,*) line ! can be #  type or #Cohesive
        !print *, 'line 885 fingerprints.f90'

        IF (line .EQ. "#Cohesive") THEN
            READ (mlff_file_unit,*)
            READ (mlff_file_unit,*)
            READ (mlff_file_unit,*) !skip # type
        END IF
        IF (.NOT. haveV1) THEN
            READ (mlff_file_unit,*) G_type, numGs
        ELSE
            READ(mlff_file_unit,*) G_type, numGs(1)
        END IF
        IF (nG1 .GT. 0) THEN
            READ (mlff_file_unit,*) !skip center ...
            DO i = 1, nG1
                READ (mlff_file_unit,*) center, neigh1, eta, Rs, rcut, fpmin, fpmax
                ! ideally have dynamic rcut for different interactions
                IF (rcut .GT. max_rcut) THEN
                    IF (max_rcut > 0.d0) THEN
                        print *, 'Error: multiple rcuts are not supported yet. Use a single rcut value'
                        STOP
                    ELSE
                        max_rcut = rcut
                    END IF
                END IF
                cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                nidx1 = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                fpParas(cidx)%g1s(nidx1)%fp_type = G_type
                fpParas(cidx)%g1s(nidx1)%species1 = center
                fpParas(cidx)%g1s(nidx1)%species1_code = cidx
                fpParas(cidx)%g1s(nidx1)%species2 = neigh1
                fpParas(cidx)%g1s(nidx1)%currIndex = fpParas(cidx)%g1s(nidx1)%currIndex + 1
                currIndex = fpParas(cidx)%g1s(nidx1)%currIndex
                fpParas(cidx)%g1s(nidx1)%etas(currIndex) = eta
                fpParas(cidx)%g1s(nidx1)%rss(currIndex) = Rs
                fpParas(cidx)%g1s(nidx1)%r_cuts(currIndex) = rcut
                fpParas(cidx)%tnFPs = fpParas(cidx)%tnFPs + 1
      
                fpminvs(1+fprange_idx(cidx),cidx) = fpmin
                fpmaxvs(1+fprange_idx(cidx),cidx) = fpmax
                fprange_idx(cidx)=fprange_idx(cidx)+1

                IF (nidx1 .EQ. 1) THEN
                    fpParas(cidx)%g1s(nidx1)%startpoint = 1
                    fpParas(cidx)%g1s(nidx1)%endpoint = fpParas(cidx)%g1s(nidx1)%endpoint + 1
                ELSE
                    fpParas(cidx)%g1s(nidx1)%startpoint = fpParas(cidx)%g1s(nidx1-1)%endpoint + 1
                    fpParas(cidx)%g1s(nidx1)%endpoint = fpParas(cidx)%g1s(nidx1)%startpoint + fpParas(cidx)%g1s(nidx1)%nFPs - 1
                END IF
            END DO
        END IF
        !print *, 'line 932 fingerprints.f90'

        IF (nG2 .GT. 0) THEN
            IF (nG1 .GT. 0) THEN
                READ (mlff_file_unit,*) !skip type ...
                IF (.NOT. haveV1) THEN
                    READ (mlff_file_unit,*) G_type, numGs
                ELSE
                    READ (mlff_file_unit,*) G_type, numGs(1)
                END IF
            ELSE
                CONTINUE
            END IF    
            !print *, 'line 941 fingerprints.f90'

            !print *, G_type, numGs
            DO k = 1, nelement
                fpParas(k)%g1_endpoint = fpParas(k)%tnFPs
                fpParas(k)%g2_startpoint = fpParas(k)%g1_endpoint + 1
            END DO
            READ (mlff_file_unit,*) !skip # center...
            DO i = 1, nG2
                READ (mlff_file_unit,*) center, neigh1, neigh2, eta, zeta, lambda, thetas, rcut, fpmin, fpmax
                ! ideally have dynamic rcut for different interactions
                IF (rcut .GT. max_rcut) THEN
                    IF (max_rcut > 0.d0) THEN
                        print *, 'Error: multiple rcuts are not supported yet. Use a single rcut value'
                        STOP
                    ELSE
                        max_rcut = rcut
                    END IF
                END IF
                cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                nidx1 = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                nidx2 = find_loc_char(uniq_elements, neigh2, SIZE(uniq_elements))
                IF (nidx1 .GT. nidx2) THEN
                    tempidx = nidx1
                    nidx1 = nidx2
                    nidx2 = tempidx
                END IF
                fpParas(cidx)%g2s(nidx1, nidx2)%fp_type = G_type
                fpParas(cidx)%g2s(nidx1, nidx2)%species1 = center
                fpParas(cidx)%g2s(nidx1, nidx2)%species1_code = cidx
                fpParas(cidx)%g2s(nidx1, nidx2)%species2 = neigh1
                fpParas(cidx)%g2s(nidx1, nidx2)%currIndex = fpParas(cidx)%g2s(nidx1, nidx2)%currIndex + 1
                currIndex = fpParas(cidx)%g2s(nidx1, nidx2)%currIndex

                fpParas(cidx)%g2s(nidx1, nidx2)%etas(currIndex)     = eta
                fpParas(cidx)%g2s(nidx1, nidx2)%zetas(currIndex)    = zeta
                fpParas(cidx)%g2s(nidx1, nidx2)%lambdas(currIndex)  = lambda
                fpParas(cidx)%g2s(nidx1, nidx2)%theta_ss(currIndex) = thetas
                fpParas(cidx)%g2s(nidx1, nidx2)%r_cuts(currIndex)   = rcut
                fpParas(cidx)%tnFPs = fpParas(cidx)%tnFPs + 1

                fpminvs(1+fprange_idx(cidx),cidx)=fpmin
                fpmaxvs(1+fprange_idx(cidx),cidx)=fpmax
                fprange_idx(cidx)=fprange_idx(cidx)+1
            END DO
            !print *, 'line 986 fingerprints.f90'

            DO i = 1, nelement
                DO j = 1, nelement
                    DO k = j, nelement
                        IF (fpParas(i)%g2s(j, k)%nFPs .EQ. 0) THEN
                            accN = 0
                        ELSE
                            accN = 1
                        END IF

                        IF ((j .EQ. 1) .AND. (k .EQ. 1)) THEN
                            fpParas(i)%g2s(j, k)%startpoint = fpParas(i)%g1_endpoint + accN
                        ELSE
                            IF (j .EQ. k) THEN
                                fpParas(i)%g2s(j, k)%startpoint = fpParas(i)%g2s(j-1, nelement)%endpoint+accN
                            ELSE
                                fpParas(i)%g2s(j, k)%startpoint = fpParas(i)%g2s(j, k-1)%endpoint+accN
                            END IF
                        END IF
                        IF (accN .EQ. 0) THEN
                            fpParas(i)%g2s(j, k)%endpoint = fpParas(i)%g2s(j, k)%startpoint
                        ELSE
                            fpParas(i)%g2s(j, k)%endpoint = fpParas(i)%g2s(j, k)%startpoint + fpParas(i)%g2s(j, k)%nFPs - 1
                        END IF
                    END DO
                END DO
            END DO
            !print *, 'line 1014 fingerprints.f90'

            !print *, fpParas(2)%g2s(2,2)%startpoint, fpParas(2)%g2s(2,2)%endpoint
            DO i = 1, nelement
                DO j = 1, nelement
                    DO k = 1, j - 1
                        !print *, i, j, k
                        IF (fpParas(i)%g2s(j,k)%nFPs .GT. 0) THEN
                            fpParas(i)%g2s(j, k)%etas       =  fpParas(i)%g2s(k, j)%etas
                            fpParas(i)%g2s(j, k)%zetas      =  fpParas(i)%g2s(k, j)%zetas
                            fpParas(i)%g2s(j, k)%lambdas    =  fpParas(i)%g2s(k, j)%lambdas
                            fpParas(i)%g2s(j, k)%theta_ss   =  fpParas(i)%g2s(k, j)%theta_ss
                            fpParas(i)%g2s(j, k)%r_cuts     =  fpParas(i)%g2s(k, j)%r_cuts
                            fpParas(i)%g2s(j, k)%startpoint =  fpParas(i)%g2s(k, j)%startpoint
                            fpParas(i)%g2s(j, k)%endpoint   =  fpParas(i)%g2s(k, j)%endpoint
                        END IF
                    END DO
                END DO
            END DO
        END IF
        !print *, 'line 1034 fingerprints.f90'

        ! Read NN parameters 
        ! Set common variables
        nelements = nelement
        !total_natoms = natoms

        ! Read weights and biases
        DO WHILE (line .NE. "#MachineLearning")
            READ (mlff_file_unit,*) line
        ENDDO
        READ (mlff_file_unit,*) model_type
        READ (mlff_file_unit,*) !Skip #Activation function type
        READ (mlff_file_unit,*) actfuncId
        !2025-02 changed so that weights and biases are at the end
        READ (mlff_file_unit,*) comment
        IF ((comment(:6) .EQ. "#Model") .AND. (.NOT. haveV1)) THEN
            WRITE(*,*) "WARNING: V2 of mlff.pyamff detected. A new mlff.pyamff will be generated."
            haveV2 = .TRUE.
        END IF
        IF ((.NOT. haveV2) .AND. (.NOT. haveV1)) THEN
            READ (mlff_file_unit,*) scaler_type
            READ (mlff_file_unit,*) slope, intercept
            READ (mlff_file_unit,*) !Skip command line
        END IF
        READ (mlff_file_unit,*) nhidlayers
        IF(ALLOCATED(nhidneurons)) THEN
          DEALLOCATE(nhidneurons)
        END IF
        ALLOCATE(nhidneurons(nhidlayers))
        READ (mlff_file_unit,*) nhidneurons
        IF(ALLOCATED(flatten_inweights)) THEN
          DEALLOCATE(flatten_inweights)
        END IF
        IF(ALLOCATED(flatten_hidweights)) THEN
          DEALLOCATE(flatten_hidweights)
        END IF
        IF(ALLOCATED(weights)) THEN
          DEALLOCATE(weights)
        END IF
        IF(ALLOCATED(biases)) THEN
          DEALLOCATE(biases)
        END IF
        !print *, 'line 1067 fingerprints.f90'

        ALLOCATE(flatten_inweights(MAXVAL(nGs)*nhidneurons(1)))
        ALLOCATE(flatten_hidweights(MAXVAL(nhidneurons)*MAXVAL(nhidneurons)))
        ALLOCATE(weights(max(MAXVAL(nGs),MAXVAL(nhidneurons)),MAXVAL(nhidneurons),nhidlayers+1,nelements))
        ALLOCATE(biases(1,MAXVAL(nhidneurons),nhidlayers+1,nelements))
        READ (mlff_file_unit,*) !Skip #Model structure
        !print *, 'line 1074 fingerprints.f90'

        DO i = 1, nelements
            ! Read atom type skipping the first symbol #.
            READ (mlff_file_unit,*) atom_type
            ! Find index of corresponding atom_type
            myid=find_loc_char(uniq_elements, atom_type, SIZE(uniq_elements))
            ! Input weights, biases
            IF ((.NOT. haveV2) .AND. (.NOT. haveV1)) READ (mlff_file_unit,*) !inputLayer weights tag
            !print*, 'myid ', myid, nGs(myid)*nhidneurons(1)
            READ (mlff_file_unit,*,IOSTAT=ios) flatten_inweights(1:nGs(myid)*nhidneurons(1))
            IF (ios .GT. 0) THEN ! 59 (intel) or 225 (nvidia) is if have "#Energy Scaling Parameters"
                haveV2_noparams = .TRUE.
                CALL initialize_NN_params(uniq_elements, seedval)
            ELSE IF (ios .LT. 0) THEN ! EOF
                CALL initialize_NN_params(uniq_elements, seedval)
            ELSE IF (ios .EQ. 0) THEN
                weights(1:nGs(myid),1:nhidneurons(1),1,myid) = &
                    reshape(flatten_inweights(1:nGs(myid)*nhidneurons(1)),(/nGs(myid),nhidneurons(1)/))
                IF ((.NOT. haveV2) .AND. (.NOT. haveV1)) READ (mlff_file_unit,*) !inputLayer bias tag
                READ (mlff_file_unit,*) biases(1,1:nhidneurons(1),1,myid)
                ! Hidden weights, biases
                DO j = 1, nhidlayers-1
                    IF ((.NOT. haveV2) .AND. (.NOT. haveV1)) READ (mlff_file_unit,*) !hiddenLayer weights tag
                    READ (mlff_file_unit,*) flatten_hidweights(1:nhidneurons(j)*nhidneurons(j+1))
                    weights(1:nhidneurons(j),1:nhidneurons(j+1),j+1,myid) = &
                        reshape(flatten_hidweights(1:nhidneurons(j)*nhidneurons(j+1)),(/nhidneurons(j), nhidneurons(j+1)/))
                    IF ((.NOT. haveV2) .AND. (.NOT. haveV1)) READ (mlff_file_unit,*) !hiddenLayer bias tag
                    READ (mlff_file_unit,*) biases(1,1:nhidneurons(j+1),j+1,myid)
                END DO
                ! Out weights, biases
                IF ((.NOT. haveV2) .AND. (.NOT. haveV1)) READ (mlff_file_unit,*) !outputLayer weights tag
                READ (mlff_file_unit,*) weights(1:nhidneurons(nhidlayers),1,nhidlayers+1,myid)
                IF ((.NOT. haveV2) .AND. (.NOT. haveV1)) READ (mlff_file_unit,*) !outputLayer bias tag
                READ (mlff_file_unit,*) biases(1,1,nhidlayers+1,myid)
                !!! IF (i ==  nelements) READ (11,*) !
            ! ELSE
            !     WRITE(*,'(A,I4,A)') "Error (",ios,") while reading in weights and biases from '" // TRIM(ADJUSTL(mlff_file)) // "'"
            !     CLOSE(mlff_file_unit)
            !     STOP
            END IF
        END DO
        IF (haveV2 .OR. haveV1) THEN
            IF (.NOT. haveV2_noparams) READ (mlff_file_unit,*) ! skip command
            READ (mlff_file_unit,*) scaler_type
            READ (mlff_file_unit,*) slope, intercept
            CLOSE(mlff_file_unit)
            IF (haveV2) CALL update_mlff(mlff_file)
            IF (haveV1) THEN 
                CALL RENAME(TRIM(mlff_file),TRIM(ADJUSTL(mlff_file)) // ".v1")
                CALL write_mlff(mlff_file)
            END IF
        ELSE
            CLOSE (mlff_file_unit)
        END IF

    END SUBROUTINE

    SUBROUTINE initialize_NN_params(uniqElems, seedval)
        USE mtmod
        IMPLICIT NONE
        !Inputs
        CHARACTER*2, DIMENSION(nelements), INTENT(IN) :: uniqElems
        INTEGER :: seedval
        !Variables
        INTEGER :: i, j, k
        print *, 'NN parameters are randomly generated by fortran'
        DO i=1, nelements
            !print *, 'line 1125 fps.f90 nGs(1): ',nGs(i)
            DO k=1, nGs(i)
                CALL gather_grnd(weights(k,1:nhidneurons(1),1,i),nhidneurons(1),seedval)
            END DO
            !print *, 'line 1128 fingerprints.f90 in_weights'!, in_weights(1:nGs(i),1:nhidneurons(1),i)
            !CALL gather_grnd(in_biases(1:nhidneurons(1),i),nhidneurons(1),seedval)
            CALL gather_grnd(biases(1,1:nhidneurons(1),1,i),nhidneurons(1),seedval)
            !print *, 'line 1131 fingerprints.f90 in_biases'!, in_biases(1:nhidneurons(1),i)
            DO j=1, nhidlayers-1
                DO k=1, nhidneurons(j)
                    CALL gather_grnd(weights(k,1:nhidneurons(j+1),j+1,i),nhidneurons(j+1),seedval)
                END DO
                !print *, 'line 1136 fingerprints.f90 hid_weights'!, hid_weights(1:nhidneurons(j),1:nhidneurons(j+1),j,i)
                CALL gather_grnd(biases(1,1:nhidneurons(j+1),j+1,i),nhidneurons(j+1),seedval)
                !print *, 'line 1138 fingerprints.f90 hid_biases'! hid_biases(1:nhidneurons(j+1),j,i)
            END DO
            CALL gather_grnd(weights(1:nhidneurons(nhidlayers),1,nhidlayers+1,i),nhidneurons(nhidlayers),seedval)
            !print *, 'out_weights=', out_weights(1:nhidneurons(nhidlayers),1,i)
            !out_biases(i)=grnd(seedval)
            CALL gather_grnd(biases(1,1,nhidlayers+1,i),1,seedval)
            !print *, 'out_biases=', out_biases(i)
        END DO

    !print *, 'END OF initalize_NN_params 1147 fingerprints.f90'
    END SUBROUTINE


    SUBROUTINE calcfps(nAtoms, pos_car, cell, symbols, nelement,forceEngine, &
                       nneigh_incell, neighs_incell, num_neigh, num_cells, neighsDefined)
                       !fps, dfps, nneigh_incell, neighs_incell, num_neigh, num_cells, neighsDefined)

        ! Parameters
        INTEGER, PARAMETER :: MAX_NEIGHS  = 100

        ! input values
        LOGICAL, OPTIONAL:: neighsDefined
        INTEGER :: nAtoms, nelement
        INTEGER :: forceEngine
        INTEGER, DIMENSION(3), OPTIONAL:: num_cells ! there's a better way to do this but not right now
        INTEGER*4, DIMENSION(nAtoms) :: symbols
        DOUBLE PRECISION, DIMENSION(nAtoms,3) :: pos_car
        DOUBLE PRECISION, DIMENSION(3,3) :: cell
        !INTEGER :: MAX_FPS  !find max number of fingerprints, element-based
!f2py   INTENT(IN) :: pos_car
!f2py   INTENT(IN) :: cell
!f2py   INTENT(IN) :: symbols
!f2py   INTENT(IN) :: forceEngine
!f2py   INTENT(IN) :: nelement
!!!!f2py    INTENT(IN) :: nAtoms, nelement
!!!f2py   INTENT(IN) :: MAX_FPS  !find max number of fingerprints, element-based

        ! variables for nlist.calc
        INTEGER :: j, k, l,i
        INTEGER :: npairs, npairs_incell, tncells, max_npairs
        INTEGER, DIMENSION(2) :: num_pairs
        INTEGER, DIMENSION(3) :: ncells
        DOUBLE PRECISION, DIMENSION(3, 3) :: supercell
        DOUBLE PRECISION :: start_T, end_T
        LOGICAL :: HAS_NEIGHS  !!!Need this to make sure .AND. works for both condition a and b  ; NOT a or b
        !COMMON fpParas, max_rcut

        ! variables for calcTriplet
        INTEGER :: ntriplets, nFPs
        DOUBLE PRECISION, DIMENSION(3, 3) :: vectsigns
        INTEGER, DIMENSION(nAtoms*MAX_NEIGHS) :: pair_offset   ! offset of neigh atom

        ! variables for DataToMatrix

        ! output values
        !DOUBLE PRECISION, DIMENSION(nAtoms, MAX_FPS) :: fps
        !DOUBLE PRECISION, DIMENSION(nAtoms, MAX_NEIGHS, 3, MAX_FPS) :: dfps
        INTEGER, DIMENSION(nAtoms, MAX_NEIGHS) :: neighs
        INTEGER, DIMENSION(nAtoms, MAX_NEIGHS) :: neighs_incell
        INTEGER, DIMENSION(nAtoms) :: num_neigh, nneigh_incell
!!!f2py   INTENT(OUT) :: fps
!!!f2py   INTENT(OUT) :: dfps
!!!!f2py   INTENT(OUT) :: neighs
!f2py   INTENT(OUT) :: neighs_incell
!f2py   INTENT(OUT) :: num_neigh
!f2py   INTENT(OUT) :: nneigh_incell
        ! CHARACTER(100) :: filename,xform, yform,zform
        ! INTEGER :: fileunit = 33, ii, jj, kk, mg
        !vectsigns = reshape((/1, 1, -1, 1, -1, -1, 1, -1, 1/), [3,3])
        !print *, 'in calcfps line 1214 fingerprints.f90'
        !print *, ' nAtoms: ',nAtoms
        CALL calcCellNum(nAtoms, pos_car, cell, max_rcut, ncells)
        tncells = ncells(1) * ncells(2) * ncells(3)
        !print *, '1217 after calcCellNum +1217 fingerprints.f90'
        !print *, 'tncells line 1150 calcfps: ',tncells
        !print *, 'tncells: ',tncells
        !print *, 'is present',PRESENT(neighsDefined)
        !print *, 'neighsDefined: ',neighsDefined
        !IF ( PRESENT(neighsDefined) .AND. neighsDefined ) THEN
        IF ( PRESENT(neighsDefined)) THEN
            IF (neighsDefined) THEN
            HAS_NEIGHS = .TRUE.
            ! ncells(1) = num_cells(1)
            ! ncells(2) = num_cells(2)
            ! ncells(3) = num_cells(3)
            ! tncells = ncells(1) * ncells(2) * ncells(3)
            !print *,'neighsDefined'
            supercell(1,:) = ncells(1) * cell(1,:)
            supercell(2,:) = ncells(2) * cell(2,:)
            supercell(3,:) = ncells(3) * cell(3,:)
            END IF
        ELSE
        HAS_NEIGHS = .FALSE.
        !print *, 'line 1230 nAtoms: ',nAtoms
        tnAtoms = nAtoms * tncells
        !print *, 'line 1232 fingerprints.f90'
        ALLOCATE(pool_pos_car(tnAtoms,3))
        ALLOCATE(pool_ids(tnAtoms))
        ALLOCATE(supersymbols(tnAtoms))
        ALLOCATE(tneighs(tnAtoms, MAX_NEIGHS))
        END IF
        !print *, 'line: 1233 fingerprints.f90'
        IF ( PRESENT(neighsDefined) .AND. neighsDefined) THEN
            supercell(1,:) = ncells(1) * cell(1,:)
            supercell(2,:) = ncells(2) * cell(2,:)
            supercell(3,:) = ncells(3) * cell(3,:)
        ELSE
        CALL genSupercell(nAtoms, pos_car,symbols, cell, ncells, tncells,  max_rcut, supercell)
                          !supercell, supersymbols, pool_pos_car, pool_ids)
        END IF
        ! print *, 'line 1175 after call genSupercell'
        max_npairs = tnAtoms*MAX_NEIGHS
        ALLOCATE(pairs(max_npairs))
        !ALLOCATE(fcutoff(max_npairs))
        ALLOCATE(pair_info(tnAtoms,MAX_NEIGHS))
        ALLOCATE(num_eachpair(tnAtoms))
        ALLOCATE(pair_indices(2, max_npairs))
        ALLOCATE(pair_global_indices(2, max_npairs))
        ALLOCATE(unitvects_pair(2, max_npairs, 3))
        ALLOCATE(gvects(max_npairs))
        ! ALLOCATE(tneighs(tnAtoms, MAX_NEIGHS))
        ALLOCATE(tneighs_incell(tnAtoms, MAX_NEIGHS))
        ALLOCATE(tnum_neigh(tnAtoms))
        !print *, 'line 1256 fingerprints.f90'
        tneighs_incell = 0
        tnum_neigh = 0
        pairs = 0
        pair_info = 0
        num_eachpair =0
        pair_indices = 0
        unitvects_pair = 0
        gvects = 0
        ! print *, 'line 1196 fingerprints.f90'
        !fps = 0
        !dfps = 0
        !print *,'line 1266 calcfps fingerprints.f90'
        CALL calcNlist(tnAtoms, MAX_NEIGHS, pool_pos_car, supercell, symbols,  max_rcut, nelement, &
                       forceEngine, rmins, tncells, nneigh_incell, &
                       npairs_incell, npairs, tnum_neigh, tneighs, tneighs_incell, (HAS_NEIGHS))
        !print *,'line 1270 calcfps fingerprints.f90'
        ! print *, 'line 1202 after calcNlist in fingerprints.f90'
        num_neigh = tnum_neigh(1:nAtoms)
        maxneighs = MAXVAL(num_neigh)
        !neighs = tneighs(1:nAtoms, :)
        neighs_incell = tneighs_incell(1:nAtoms, :)
        neighs = tneighs_incell(1:nAtoms, :)

        CALL allocate_outputs(nAtoms, max_fps, maxneighs)
        ! print *, 'line 1210 after allocate_outputs in fps.f90'
        !CALL calcg1s(nAtoms, npairs_incell, symbols, maxneighs, MAX_NEIGHS, max_fps, neighs, & 
        CALL calcg1s(nAtoms, npairs_incell, symbols, maxneighs, MAX_NEIGHS, neighs, & 
                     num_neigh, max_npairs)
                     !fps, dfps)
        !print *,'line 1284 fingerprints.f90'
        ! print *, 'line 1215 after calcg1s in fingerprints.f90'
        IF (nG2 .GT. 0) THEN
           ALLOCATE(angles(3, nAtoms*maxneighs*maxneighs))
           ALLOCATE(angles_indices(3, nAtoms*maxneighs*maxneighs))
           CALL calcTriplet(nAtoms, tnAtoms, MAX_NEIGHS, max_rcut, npairs, ntriplets)
           !CALL cpu_time(start_T)
           !CALL calcg2s_n(natoms, npairs, MAX_NPAIRS, ntriplets, maxneighs, MAX_NEIGHS, max_fps, &
           CALL calcg2s_n(natoms, npairs, MAX_NPAIRS, ntriplets, maxneighs, MAX_NEIGHS, &
                        symbols, max_rcut, neighs, num_neigh, &
                          unitvects_pair(1, :, :))
                          !fps, dfps)
           !CALL cpu_time(end_T)
           !print*, 'g2 time used', end_T-start_T
           DEALLOCATE(angles)
           DEALLOCATE(angles_indices)
        END IF
        !print *,'line 1301 fingerprints.f90'
        !print *, 'END of Calc FPS fingerprints.f90'
        ! WRITE(filename,'(A)') "fps.dat"
        ! OPEN(22,file=trim(filename),access="append")
        ! WRITE(22,'(A)') "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        ! WRITE(22,'(A)') "%%%%%%                   FPS                  %%%%%%"
        ! WRITE(22,'(A)') "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        ! WRITE(22,*), 'FINGERPRINTS'
        ! DO j = 1, natoms
        !     WRITE(22,'(T12,i4,": ",100F15.8)'), j,fps(j,:)
        ! ENDDO
        ! CLOSE(22)

        ! WRITE(filename,'(A)') "dfps.dat"
        ! OPEN(fileunit,file=trim(filename),access="append") ! append in case have multiple images
        ! WRITE(fileunit,'(A)') "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        ! WRITE(fileunit,'(A)') "%%%%%%                  DFPS                  %%%%%%"
        ! WRITE(fileunit,'(A)') "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        ! WRITE(fileunit,'(A)') "<atom_no>: (<neigh_no/self>) -> <dir> [<dFPS>]"
        ! WRITE(fileunit,'(A)') "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        ! WRITE(fileunit,'(A)') "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        ! mg = nG1+nG2
        ! WRITE(xform,'(A,I0,A)') '(T10," -> x : [",',mg,'f10.5,"]")'
        ! WRITE(yform,'(A,I0,A)') '(T10," -> y : [",',mg,'f10.5,"]")'
        ! WRITE(zform,'(A,I0,A)') '(T10," -> z : [",',mg,'f10.5,"]")'
        ! FLUSH(fileunit)
        ! DO ii = 1,natoms
        !     kk = num_neigh(ii)+1
        !     ! print self
        !     WRITE(fileunit,'(I3)') ii
        !     WRITE(fileunit,'(T4,A)') "-> self"
        !     WRITE(fileunit,xform) dfps(ii,1,1,:)
        !     WRITE(fileunit,yform) dfps(ii,1,2,:)
        !     WRITE(fileunit,zform) dfps(ii,1,3,:)
        !     DO jj=2,kk
        !         WRITE(fileunit,'(T4,A,I0)') "-> neigh no.",jj-1
        !         WRITE(fileunit,xform) dfps(ii,jj,1,:)
        !         WRITE(fileunit,yform) dfps(ii,jj,2,:)
        !         WRITE(fileunit,zform) dfps(ii,jj,3,:)
        !     END DO
        !     FLUSH(fileunit)
        ! END DO
        ! CLOSE(fileunit)
    END SUBROUTINE

    !SUBROUTINE calcg1s(nAtoms, npairs, symbols, maxneighs, MAX_NEIGHS, max_fps, neighs, & 
    SUBROUTINE calcg1s(nAtoms, npairs, symbols, maxneighs, MAX_NEIGHS, neighs, & 
                       num_neigh, max_npairs)
                       !fps, dfps)

        ! input values
         INTEGER :: nAtoms, npairs, maxneighs, MAX_NEIGHS, max_npairs
        INTEGER, DIMENSION(nAtoms) :: symbols
         INTEGER, DIMENSION(nAtoms, MAX_NEIGHS) :: neighs
         !DOUBLE PRECISION, DIMENSION(max_npairs) :: pairs  
         !INTEGER, DIMENSION(2, max_npairs) :: pair_indices 
         !DOUBLE PRECISION, DIMENSION(2, max_npairs, 3) :: unitvects_pair
        INTEGER :: nfps
        INTEGER, DIMENSION(nAtoms) :: num_neigh
         !INTEGER, DIMENSION(2,max_npairs) :: pair_global_indices   ! offset of neigh atom
         !TYPE (fingerprints), DIMENSION(nelement) :: fps

        ! variables
         INTEGER :: i, k, j, lid
        INTEGER :: center_1, index_1, sp_1, sp_2, ep_1, ep_2
        INTEGER :: center_2, index_2, o_index_1, o_index_2
        INTEGER, DIMENSION(1) :: neighloc_1, neighloc_2
         !TYPE (fingerprintsData), DIMENSION(2, npairs) :: fpdata_temp
        TYPE (fingerprintsData), DIMENSION(2) :: fpdata_temp
        ! outputs
         !DOUBLE PRECISION, DIMENSION(nAtoms, MAX_FPS) :: fps
         !DOUBLE PRECISION, DIMENSION(nAtoms, MAX_NEIGHS, 3, MAX_FPS) :: dfps

         !Loop over all pairs: Calculate fps and common parts in dfps
         !print *, shape(pair_indices)
         !fps = 0.0
         !dfps = 0.0
        DO i = 1, npairs
            !centered on indices_1(i)
            !print *, 'starting i', i, npairs
            index_1 = pair_global_indices(1,i)
            index_2 = pair_global_indices(2,i)
            IF ((index_1>nAtoms).AND.(index_2>nAtoms)) THEN
                print*, 'Warning ghost atoms found in head of pairs'
                CYCLE
            END IF
            o_index_1 = pair_indices(1,i)
            o_index_2 = pair_indices(2,i)
            center_1 = symbols(o_index_1)
            center_2 = symbols(o_index_2)
            IF (center_1 == 0) THEN
                CYCLE
            END IF
             !print*, 'fetching sp ep'
             !print *, 'pair', i, index_1, index_2
             !print *,i, pairs(i)
            IF ((index_1<=nAtoms).AND.(index_2<=nAtoms)) THEN
                sp_1 = fpParas(center_1)%g1s(center_2)%startpoint
                ep_1 = fpParas(center_1)%g1s(center_2)%endpoint
                nFPs = fpParas(center_1)%g1s(center_2)%nFPs
                IF (nFPs .EQ. 0) THEN
                    CYCLE
                END IF
                 !print *,'allocating', center_1, center_2, nFPs
                 !possible optimization: allocate the size based on the pair only not in total number of fps
                ALLOCATE(fpdata_temp(1)%gs(fpParas(center_1)%tnFPs))
                ALLOCATE(fpdata_temp(1)%pre_dgdxs(fpParas(center_1)%tnFPs))
                fpdata_temp(1)%gs = 0
                fpdata_temp(1)%pre_dgdxs = 0
                CALL fg1(pairs(i), nFPs, &
                      fpParas(center_1)%g1s(center_2)%etas, &
                      fpParas(center_1)%g1s(center_2)%rss, &
                      fpParas(center_1)%g1s(center_2)%r_cuts, &
                      fpdata_temp(1)%gs(sp_1:ep_1))
                CALL fdg1(pairs(i), nFPs, &
                      fpParas(center_1)%g1s(center_2)%etas, &
                      fpParas(center_1)%g1s(center_2)%rss, &
                      fpParas(center_1)%g1s(center_2)%r_cuts, &
                      fpdata_temp(1)%pre_dgdxs(sp_1:ep_1))
                 !print*, 'g1 fdg1'
                 !centered on indices_2(i)
                
                IF (center_1 .EQ. center_2) THEN
                    ALLOCATE(fpdata_temp(2)%gs(fpParas(center_2)%tnFPs))
                    ALLOCATE(fpdata_temp(2)%pre_dgdxs(fpParas(center_2)%tnFPs))
                    sp_2 = sp_1
                    ep_2 = ep_1
                    fpdata_temp(2)%gs = fpdata_temp(1)%gs
                    fpdata_temp(2)%pre_dgdxs = fpdata_temp(1)%pre_dgdxs
                ELSE
                    sp_2 = fpParas(center_2)%g1s(center_1)%startpoint
                    ep_2 = fpParas(center_2)%g1s(center_1)%endpoint
                    nFPs = fpParas(center_2)%g1s(center_1)%nFPs
                    ALLOCATE(fpdata_temp(2)%gs(fpParas(center_2)%tnFPs))
                    ALLOCATE(fpdata_temp(2)%pre_dgdxs(fpParas(center_2)%tnFPs))
                    fpdata_temp(2)%gs = 0
                    fpdata_temp(2)%pre_dgdxs = 0
                    CALL fg1(pairs(i), nFPs, &
                          fpParas(center_2)%g1s(center_1)%etas, &
                          fpParas(center_2)%g1s(center_1)%rss, &
                          fpParas(center_2)%g1s(center_1)%r_cuts, &
                          fpdata_temp(2)%gs(sp_2:ep_2))
                    CALL fdg1(pairs(i), nFPs, &
                          fpParas(center_2)%g1s(center_1)%etas, &
                          fpParas(center_2)%g1s(center_1)%rss, &
                          fpParas(center_2)%g1s(center_1)%r_cuts, &
                          fpdata_temp(2)%pre_dgdxs(sp_2:ep_2))
                END IF
                
                fps(index_1, sp_1:ep_1) = fps(index_1, sp_1:ep_1) + &
                                            fpdata_temp(1)%gs(sp_1:ep_1)
                !Only for atoms in original cell not for ghost atoms
                fps(index_2, sp_2:ep_2) = fps(index_2, sp_2:ep_2) + &
                                            fpdata_temp(2)%gs(sp_2:ep_2)
                neighloc_2 = find_loc_int(neighs(index_2,:), index_1, SIZE(neighs(index_2,:))) +1
                neighloc_1 = find_loc_int(neighs(index_1,:), index_2, SIZE(neighs(index_1,:))) + 1
                !derivative
                DO k = 1, 3
                    !center on itself
                    dfps(index_1, 1, k, sp_1:ep_1) = &
                    dfps(index_1, 1, k, sp_1:ep_1) - &
                    fpdata_temp(1)%pre_dgdxs(sp_1:ep_1) * &
                    unitvects_pair(1, i, k)
      
                    dfps(index_2, 1, k,sp_2:ep_2) = &
                    dfps(index_2, 1, k,sp_2:ep_2) - &
                    fpdata_temp(2)%pre_dgdxs(sp_2:ep_2) * &
                    unitvects_pair(2, i, k)
                    !Center on index_2 but w.r.t. index_1
                     dfps(index_1, neighloc_1(1), k,sp_2:ep_2) = dfps(index_1, neighloc_1(1), k,sp_2:ep_2) + &
                    fpdata_temp(2)%pre_dgdxs(sp_2:ep_2) * &
                    unitvects_pair(2, i, k)
                    !Center on index_1 but w.r.t. index_2
                     dfps(index_2, neighloc_2(1), k,sp_1:ep_1) = dfps(index_2, neighloc_2(1), k,sp_1:ep_1) + &
                    fpdata_temp(1)%pre_dgdxs(sp_1:ep_1) * &
                    unitvects_pair(1, i, k)
                     !print*, dfps(index_2, index_1, k,sp_1:ep_1)
                END DO
                DEALLOCATE(fpdata_temp(1)%gs)
                DEALLOCATE(fpdata_temp(1)%pre_dgdxs)
                DEALLOCATE(fpdata_temp(2)%gs)
                DEALLOCATE(fpdata_temp(2)%pre_dgdxs)
            END IF

            IF ((index_1<=nAtoms).AND.(index_2>nAtoms)) THEN
                sp_1 = fpParas(center_1)%g1s(center_2)%startpoint
                ep_1 = fpParas(center_1)%g1s(center_2)%endpoint
                nFPs = fpParas(center_1)%g1s(center_2)%nFPs
                IF (nFPs .EQ. 0) THEN
                    CYCLE
                END IF
                !print *,'allocating', center_1, center_2, nFPs
                !possible optimization: allocate the size based on the pair only not in total number of fps
                ALLOCATE(fpdata_temp(1)%gs(fpParas(center_1)%tnFPs))
                ALLOCATE(fpdata_temp(1)%pre_dgdxs(fpParas(center_1)%tnFPs))
                fpdata_temp(1)%gs = 0
                fpdata_temp(1)%pre_dgdxs = 0
                CALL fg1(pairs(i), nFPs, &
                         fpParas(center_1)%g1s(center_2)%etas, &
                         fpParas(center_1)%g1s(center_2)%rss, &
                         fpParas(center_1)%g1s(center_2)%r_cuts, &
                         fpdata_temp(1)%gs(sp_1:ep_1))
                CALL fdg1(pairs(i), nFPs, &
                          fpParas(center_1)%g1s(center_2)%etas, &
                          fpParas(center_1)%g1s(center_2)%rss, &
                          fpParas(center_1)%g1s(center_2)%r_cuts, &
                          fpdata_temp(1)%pre_dgdxs(sp_1:ep_1))

                sp_2 = fpParas(center_2)%g1s(center_1)%startpoint
                ep_2 = fpParas(center_2)%g1s(center_1)%endpoint
                nFPs = fpParas(center_2)%g1s(center_1)%nFPs
                ALLOCATE(fpdata_temp(2)%gs(fpParas(center_2)%tnFPs))
                ALLOCATE(fpdata_temp(2)%pre_dgdxs(fpParas(center_2)%tnFPs))
                fpdata_temp(2)%gs = 0
                fpdata_temp(2)%pre_dgdxs = 0
                CALL fg1(pairs(i), nFPs, &
                         fpParas(center_2)%g1s(center_1)%etas, &
                         fpParas(center_2)%g1s(center_1)%rss, &
                         fpParas(center_2)%g1s(center_1)%r_cuts, &
                         fpdata_temp(2)%gs(sp_2:ep_2))
                CALL fdg1(pairs(i), nFPs, &
                          fpParas(center_2)%g1s(center_1)%etas, &
                          fpParas(center_2)%g1s(center_1)%rss, &
                          fpParas(center_2)%g1s(center_1)%r_cuts, &
                          fpdata_temp(2)%pre_dgdxs(sp_2:ep_2))
            
                !Only for atoms in original cell not for ghost atoms
                !print*, 'g1 fdg1'
                !centered on indices_2(i)
                fps(index_1, sp_1:ep_1) = fps(index_1, sp_1:ep_1) + &
                                          fpdata_temp(1)%gs(sp_1:ep_1)
                !Only for atoms in original cell not for ghost atoms
                 lid = index_2 - nAtoms * (int((index_2 - 0.1)/nAtoms))
                 IF (o_index_2 == index_1) THEN
                   neighloc_1 = 1
                 ELSE
                   neighloc_1 = find_loc_int(neighs(index_1,:), o_index_2, SIZE(neighs(index_1,:))) + 1
                 END IF
                neighloc_2 = find_loc_int(neighs(o_index_2,:), index_1, SIZE(neighs(o_index_2,:))) +1
                !derivative
                DO k = 1, 3
                    !center on itself
                    dfps(index_1, 1, k, sp_1:ep_1) = &
                    dfps(index_1, 1, k, sp_1:ep_1) - &
                    fpdata_temp(1)%pre_dgdxs(sp_1:ep_1) * &
                    unitvects_pair(1, i, k)

                    !Center on index_2 but w.r.t. index_1
                     dfps(index_1, neighloc_1(1), k,sp_2:ep_2) = dfps(index_1, neighloc_1(1), k,sp_2:ep_2) + &
                    fpdata_temp(2)%pre_dgdxs(sp_2:ep_2) * &
                    unitvects_pair(2, i, k)
                END DO
                DEALLOCATE(fpdata_temp(1)%gs)
                DEALLOCATE(fpdata_temp(1)%pre_dgdxs)
                DEALLOCATE(fpdata_temp(2)%gs)
                DEALLOCATE(fpdata_temp(2)%pre_dgdxs)
            END IF
            !index_1 is always smaller than index_2
            IF ((index_1>nAtoms).AND.(index_2<=nAtoms)) THEN
                print*, 'WARNING INDEX 2 < INDEX 1'
                sp_2 = fpParas(center_2)%g1s(center_1)%startpoint
                ep_2 = fpParas(center_2)%g1s(center_1)%endpoint
                nFPs = fpParas(center_2)%g1s(center_1)%nFPs
                ALLOCATE(fpdata_temp(2)%gs(fpParas(center_2)%tnFPs))
                ALLOCATE(fpdata_temp(2)%pre_dgdxs(fpParas(center_2)%tnFPs))
                fpdata_temp(2)%gs = 0
                fpdata_temp(2)%pre_dgdxs = 0
                CALL fg1(pairs(i), nFPs, &
                         fpParas(center_2)%g1s(center_1)%etas, &
                         fpParas(center_2)%g1s(center_1)%rss, &
                         fpParas(center_2)%g1s(center_1)%r_cuts, &
                         fpdata_temp(2)%gs(sp_2:ep_2))
                CALL fdg1(pairs(i), nFPs, &
                          fpParas(center_2)%g1s(center_1)%etas, &
                          fpParas(center_2)%g1s(center_1)%rss, &
                          fpParas(center_2)%g1s(center_1)%r_cuts, &
                          fpdata_temp(2)%pre_dgdxs(sp_2:ep_2))
                !Only for atoms in original cell not for ghost atoms
                fps(index_2, sp_2:ep_2) = fps(index_2, sp_2:ep_2) + &
                                          fpdata_temp(2)%gs(sp_2:ep_2)
                neighloc_2 = find_loc_int(neighs(index_2,:), index_1, SIZE(neighs(index_2,:))) +1
                neighloc_1 = find_loc_int(neighs(o_index_1,:), index_2, SIZE(neighs(o_index_1,:))) + 1
                 lid = index_1 - nAtoms * (int((index_1 - 0.1)/nAtoms))
                DO k = 1, 3
                     !center on itself
                     dfps(index_2, index_2, k,sp_2:ep_2) = &
                                    dfps(index_2, index_2, k,sp_2:ep_2) - &
                    fpdata_temp(2)%pre_dgdxs(sp_2:ep_2) * &
                    unitvects_pair(2, i, k)
                    !Center on index_2 but w.r.t. index_1
                    IF (gvects(i)>0) THEN
                        dfps(o_index_1, lid, k,sp_2:ep_2) = dfps(o_index_1, lid, k,sp_2:ep_2) + &
                        fpdata_temp(2)%pre_dgdxs(sp_2:ep_2) * &
                        unitvects_pair(2, i, k)
                    END IF
                END DO
                DEALLOCATE(fpdata_temp(2)%gs)
                DEALLOCATE(fpdata_temp(2)%pre_dgdxs)
            END IF
        END DO
         !print *, 'fps'
         !DO i = 1, natoms
         !    print *, i, fpdata(i)%gs
         !END DO
         !print *, '-dfps'
         !DO i = 1, natoms
         !    DO j = 1, 3
         !        !DO k = 1,num_neigh(i)+1
         !        DO k = 1,natoms
         !            print *, i, k, j
         !           !print *, '         ', fpdata(i)%dgdxs(k)%dgdx(j,:)
         !            print *, '         ', dfps(i,k,j,1:2)
         !        END DO
         !    END DO
         !END DO
    END SUBROUTINE

    !2023-04-09
    !SUBROUTINE calcg2s_n(natoms, npairs, MAX_NPAIRS, ntriplets, maxneighs, MAX_NEIGHS, max_fps, symbols, r_cut, neighs, num_neigh, &
    SUBROUTINE calcg2s_n(natoms, npairs, MAX_NPAIRS, ntriplets, maxneighs, MAX_NEIGHS, symbols, r_cut, neighs, num_neigh, &
                       unitvects)
                       !unitvects, fps, dfps)

        ! Define inputs
        INTEGER :: nAtoms, npairs, ntriplets, maxneighs
        INTEGER :: MAX_NPAIRS, MAX_NEIGHS
        DOUBLE PRECISION :: r_cut
        INTEGER, DIMENSION(nAtoms) :: symbols
        INTEGER, DIMENSION(nAtoms, MAX_NEIGHS) :: neighs
        DOUBLE PRECISION, DIMENSION(MAX_NPAIRS, 3) :: unitvects
        DOUBLE PRECISION, DIMENSION(3,3) :: vectsigns
        INTEGER, DIMENSION(nAtoms) :: num_neigh

        ! Define variables
        DOUBLE PRECISION, PARAMETER :: Pi = 4.*ATAN(1.0_8)
        DOUBLE PRECISION, PARAMETER :: radians = Pi/180.0_8 !same as np.radians
        LOGICAL :: calc_dth
        INTEGER :: i, j, k, lid
        INTEGER :: sp, ep, nFPs
        INTEGER, DIMENSION(1) :: neigh1_loc, neigh2_loc
        INTEGER, DIMENSION(2) :: currneighs
        INTEGER, DIMENSION(3) :: center, pairIndex, csymbol
        INTEGER, DIMENSION(3,2) :: neighbor
        INTEGER, DIMENSION(3,3) :: indices, vertices, gcenter
        DOUBLE PRECISION :: fexpt, fct, inv_rcut 
        DOUBLE PRECISION, DIMENSION(3) :: sins, fcuts_p, Rs, inv_sin, cos_theta
        DOUBLE PRECISION, DIMENSION(3,3) :: uvects !(uvect_ij, uvect_ik, uvect_jk)
        DOUBLE PRECISION, DIMENSION(3,3) :: dfc, vect, dexps, dfc_i, dfc_j
        DOUBLE PRECISION, DIMENSION(3,3,3) :: dth_ijk
        DOUBLE PRECISION, DIMENSION(npairs) :: pair_sins, fcuts, fexp
        !DOUBLE PRECISION, DIMENSION(3, max_fps) :: theta_ss, reduced_etas, zetas, lambdas

        ! Define outputs
        !DOUBLE PRECISION, DIMENSION(nAtoms, MAX_FPS) :: fps
        !DOUBLE PRECISION, DIMENSION(nAtoms, MAX_NEIGHS, 3, MAX_FPS) :: dfps

        ! Loop over all angles and calculate common parts
        inv_rcut = 1/r_cut
        CALL commonparts(npairs, pairs, inv_rcut,&
                          pair_sins, fcuts, fexp)
        DO i = 1, ntriplets
            ! fetch the pair index for each pair in triplet
            pairIndex(1) = angles_indices(1, i)
            pairIndex(2) = angles_indices(2, i)
            pairIndex(3) = angles_indices(3, i)
            ! fetch the center index for each angle in triplet
            center(1) = pair_indices(1, pairIndex(1))
            center(2) = pair_indices(2, pairIndex(1))
            center(3) = pair_indices(2, pairIndex(2))

            ! record the indices of vertex of the triplet
            csymbol = symbols(center)
            neighbor(1,:) = [csymbol(2), csymbol(3)] !neighbor(:,1)
            neighbor(2,:) = [csymbol(1), csymbol(3)]
            neighbor(3,:) = [csymbol(1), csymbol(2)]

            gcenter(1, :) = [pair_global_indices(1, pairIndex(1)), &   !i
                            pair_global_indices(2, pairIndex(1)), &    !j
                            pair_global_indices(2, pairIndex(2))]      !k
            IF ((gcenter(1, 1) > nAtoms) .AND. &
                (gcenter(1, 2) > nAtoms) .AND. &
                (gcenter(1, 3) > nAtoms)) THEN
               CYCLE
            END IF

            fexpt = PRODUCT(fexp(pairIndex))
            fct   = PRODUCT(fcuts(pairIndex))

            !Derivative of cutoff functions
            Rs = pairs(pairIndex(:))
            !sins = Pi/(2.*r_cut) * pair_sins(pairIndex(:))
            sins = Pi* 0.5 * inv_rcut * pair_sins(pairIndex(:))
            uvects = unitvects(pairIndex(:), :)
            DO j=1,3
               dfc_i(j, :) = sins(j) * uvects(j,:)  !(x,y,z)
               dfc_j(j, :) = -dfc_i(j, :)
               vect(j,:) = Rs(j) * uvects(j,:)   !(x, y, z)
            END DO
            fcuts_p(1) = fcuts(pairIndex(1)) * fcuts(pairIndex(2))
            fcuts_p(2) = fcuts(pairIndex(2)) * fcuts(pairIndex(3))
            fcuts_p(3) = fcuts(pairIndex(1)) * fcuts(pairIndex(3))
            dfc(1,:) = fcuts_p(2) * dfc_i(1, :) + fcuts_p(3) * dfc_i(2,:) !(x,y,z)
            dfc(2,:) = fcuts_p(2) * dfc_j(1, :) + fcuts_p(1) * dfc_i(3,:) !(x,y,z)
            dfc(3,:) = fcuts_p(3) * dfc_j(2, :) + fcuts_p(1) * dfc_j(3,:) !(x,y,z)
            dexps(1,:) = vect(1,:)  + vect(2,:)
            dexps(2,:) = -vect(1,:)  + vect(3,:)
            dexps(3,:) = -vect(2,:)  - vect(3,:)
            cos_theta = COS(angles(:,i))                    !theta_ijk jik kij
            calc_dth = (abs(cos_theta(1)) < 0.999999)
            IF (calc_dth) THEN
                inv_sin = -1/sqrt(1-cos_theta * cos_theta) !r_ijik r_jijk r_kjki
                dth_ijk(1,2,:) = inv_sin(1) * (uvects(2,:) - uvects(1,:)*cos_theta(1))/Rs(1)    !rij_ik       !wrt j
                dth_ijk(1,3,:) = inv_sin(1) * (uvects(1,:) - uvects(2,:)*cos_theta(1))/Rs(2)       !wrt k
                dth_ijk(1,1,:) = - dth_ijk(1,2,:) - dth_ijk(1,3,:)   !wrt i
                !center on j
                dth_ijk(2,1,:) = inv_sin(2) * (uvects(3,:) + uvects(1,:)*cos_theta(2))/Rs(1)  !rij_jk
                dth_ijk(2,3,:) = -inv_sin(2) * (uvects(1,:) + uvects(3,:)*cos_theta(2))/Rs(3) !rjk_ij
                dth_ijk(2,2,:) = -dth_ijk(2,1,:) - dth_ijk(2,3,:)
                !center on k
                dth_ijk(3,1,:) =  -inv_sin(3) * (uvects(3,:) - uvects(2,:)*cos_theta(3))/Rs(2) !rjk_ik
                dth_ijk(3,2,:) =  -inv_sin(3) * (uvects(2,:) - uvects(3,:)*cos_theta(3))/Rs(3) !rik_jk
                dth_ijk(3,3,:) =  -dth_ijk(3,2,:) - dth_ijk(3,1,:)
            ELSE
                dth_ijk = 0.0
            END IF

            DO j = 1, 3
                currneighs = neighbor(j,:)
                nFPs = fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%nFPs
                IF (nFPs == 0) THEN
                    CYCLE
                END IF
                CALL fg2_dg(nFPs, angles(j,i), inv_rcut, fexpt, fct, &
                             dexps, dth_ijk(j,:,:), dfc, calc_dth,&
                         fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%theta_ss, &
                         fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%etas, &
                         fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%zetas, &
                         fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%lambdas, &
                             fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%g2, &
                             fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%dg2) 
                             !g2, dg2(:,j,:,:))

                sp = fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%startpoint
                ep = fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%endpoint
                IF (gcenter(1,j)<=nAtoms) THEN
                   !print *, 'j: ',j
                   !print *, 'currneighs(1): ',currneighs(1)
                   !print *, 'currneighs(2) +1656 fps.f90: ',currneighs(2)
                   !print *, 'ccsymbol(j)): ',csymbol(j)
                   !print *, 'fpParas term ',fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%dg2(:,j,:)
                   !print *, 'size fpParas term ',size(fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%dg2(:,j,:))
                   !print *, 'size dfps(gcenter(1,j), sp:ep, 1, :): ',size(dfps(gcenter(1,j), 1, :,sp:ep))(gcenter(1,j), 1,:,sp:ep)
                   fps(gcenter(1, j), sp:ep) = fps(gcenter(1,j), sp:ep) + &
                                               fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%g2(:)  
                                               !g2(1:nFPs)
                   dfps(gcenter(1,j),  1,:, sp:ep) = &            !wrt itself
                                     dfps(gcenter(1,j), 1, :,sp:ep) + &
                                     fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%dg2(j,:,:)
                                     !dg2(1:nFPs, j, :) 
                   DO k =1, 3
                      IF (j==k) THEN
                         CYCLE
                    END IF
                      IF (gcenter(1,k) <= nAtoms) THEN
                         neigh1_loc = find_loc_int(neighs(center(k),:), gcenter(1,j), SIZE(neighs(center(k),:))) + 1
                         dfps(center(k), neigh1_loc(1), :,sp:ep) = &
                                           dfps(center(k), neigh1_loc(1), :,sp:ep) + &
                                           fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%dg2(k,:,:)
                                           !dg2(1:nFPs, k, :) 
                    END IF
            END DO
                ELSE
                   DO k =1, 3
                      IF (j==k) THEN
                         CYCLE
                END IF
                      IF (gcenter(1,k) <= nAtoms) THEN
                         IF (center(j) == center(k)) THEN
                           neigh1_loc = 1
                         ELSE
                           neigh1_loc = find_loc_int(neighs(center(k),:), center(j), SIZE(neighs(center(k),:))) + 1
                    END IF
                         dfps(center(k), neigh1_loc(1), :,sp:ep) = &
                                           dfps(center(k), neigh1_loc(1), :,sp:ep) + &
                                           fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%dg2(k,:,:)
                                           !dg2(1:nFPs, k, :) 
                    END IF
                END DO
                    END IF
                END DO
            END DO
    END SUBROUTINE
    
    SUBROUTINE cleanup() BIND(C,name='cleanup')
         USE, INTRINSIC :: iso_c_binding
         IF(ALLOCATED(fpParas)) THEN
         DEALLOCATE(fpParas)
         END IF
         IF(ALLOCATED(uniq_elements)) THEN
         DEALLOCATE(uniq_elements)
         END IF
         IF(ALLOCATED(rmins)) THEN
         DEALLOCATE(rmins)
         END IF
    END SUBROUTINE

    SUBROUTINE atomsCleanup()
         DEALLOCATE(pairs)
         !DEALLOCATE(pair_start)
         !DEALLOCATE(pair_end)
         DEALLOCATE(gvects)
         DEALLOCATE(tneighs_incell)
         DEALLOCATE(num_eachpair)
         DEALLOCATE(pair_indices)
         DEALLOCATE(pair_info)
         DEALLOCATE(pair_global_indices)
         DEALLOCATE(unitvects_pair)
         DEALLOCATE(tneighs)
         DEALLOCATE(tnum_neigh)
         DEALLOCATE(pool_pos_car)
         DEALLOCATE(pool_ids)
         DEALLOCATE(supersymbols)
    END SUBROUTINE
 
    SUBROUTINE cleanup_ase()
        IF(ALLOCATED(cohees)) THEN
          DEALLOCATE(cohees)
        END IF
        IF(ALLOCATED(nGs)) THEN
          DEALLOCATE (nGs)
        END IF
        IF(ALLOCATED(fpminvs)) THEN
          DEALLOCATE (fpminvs)
        END IF
        IF(ALLOCATED(fpmaxvs)) THEN
          DEALLOCATE (fpmaxvs)
        END IF
        IF(ALLOCATED(diffs)) THEN
          DEALLOCATE (diffs)
        END IF
        IF(ALLOCATED(rmins)) THEN
          DEALLOCATE (rmins)
        END IF
        IF(ALLOCATED(uniq_elements)) THEN
          DEALLOCATE (uniq_elements)
        END IF
        IF(ALLOCATED(fpParas)) THEN
          DEALLOCATE (fpParas)
        END IF
        IF(ALLOCATED(nhidneurons)) THEN
          DEALLOCATE (nhidneurons)
        END IF
        IF(ALLOCATED(biases)) THEN
          DEALLOCATE (biases)
        END IF
        IF(ALLOCATED(weights)) THEN
          DEALLOCATE (weights)
        END IF
    END SUBROUTINE

    !for active learning where re-loading mlff is required
    SUBROUTINE clean_mlff() 
        IF(ALLOCATED(cohees)) THEN
          DEALLOCATE(cohees)
        END IF
        DEALLOCATE(rmins)
        DEALLOCATE(uniq_elements)
        DEALLOCATE(nGs)
        DEALLOCATE(fpminvs)
        DEALLOCATE(fpmaxvs)
        DEALLOCATE(diffs)
        DEALLOCATE(fpParas)
        DEALLOCATE(nhidneurons)
        !DEALLOCATE(flatten_inweights)
        !DEALLOCATE(flatten_hidweights)
        DEALLOCATE(weights)
        DEALLOCATE(biases)
    END SUBROUTINE

    SUBROUTINE update_mlff_model(mlff_file,dest_file,prefix)
        !NOTE: This copies the old mlff_file and just updates the 
        !      weights and biases section at the end
        IMPLICIT NONE
        ! Variables
        INTEGER :: nelement, nFPs
        INTEGER :: i, j, k
        INTEGER :: io, d_unit=22, s_unit=33
        CHARACTER(len=50) :: mlff_file, out_file
        CHARACTER(len=50),OPTIONAL :: dest_file ! if don't want the same name
        CHARACTER(len=50),OPTIONAL :: prefix ! if want the same name, what should the prefix be: prefix_
        CHARACTER(len=256) :: line
        CHARACTER(len=25) :: format, tag
        LOGICAL :: exists, haveDest = .FALSE., havePre = .FALSE.
        
        inquire(file=TRIM(ADJUSTL(mlff_file)), exist=exists)
        IF (.NOT. exists) THEN
            WRITE(*,*) "ERROR: '" // TRIM(ADJUSTL(mlff_file)) // "' does not exist. It cannot be updated."
            ! STOP
            !TODO: If it doesn't exist, call the write_mlff function instead of erroring
            CALL write_mlff(mlff_file)
            RETURN
        END IF
        !TODO: Check both be opened
        IF (PRESENT(dest_file)) THEN
            IF (LEN(TRIM(dest_file)) .NE. 0) THEN
                out_file = dest_file
                haveDest = .TRUE.
            END IF
        ELSE IF(PRESENT(prefix)) THEN
            IF (LEN(TRIM(prefix)) .NE. 0) THEN
                tag = TRIM(ADJUSTL(prefix))
                havePre = .TRUE.
            END IF
        END IF
        IF (.NOT. haveDest) THEN
            out_file = mlff_file
            IF (.NOT. havePre) THEN
                tag = "bkp"
            ENDIF
            CALL RENAME(mlff_file,TRIM(ADJUSTL(tag)) // "_" // TRIM(ADJUSTL(mlff_file)))
            mlff_file = TRIM(ADJUSTL(tag)) // "_" // TRIM(ADJUSTL(mlff_file))
        END IF
        OPEN (newunit=s_unit, file=mlff_file, status='old',action='read',iostat=io)
        OPEN (newunit=d_unit, file=out_file, status='replace',action='write',iostat=io)
        !copy until we hit Model Parameters
        DO
                READ(s_unit,'(A)',iostat=io) line
                IF (io .NE. 0) THEN ! error >, EOF <, either case we have a problem
                    WRITE(*,*) "Error while reading file '" // mlff_file // "'"
                    CLOSE(s_unit)
                    CLOSE(d_unit)
                    STOP
                END IF
                WRITE(d_unit,'(A)') TRIM(line)
                IF (line .EQ. "#Model Parameters") THEN
                    CLOSE(s_unit)
                    EXIT 
                END IF
        END DO

        ! Write NN parameters
        DO i = 1, nelements
            ! Write atom type
            WRITE (d_unit,'(A)') uniq_elements(i)
            WRITE(format,'(A)') '(5(4X,ES16.8))' ! per 5
            ! weights, biases
            WRITE(d_unit,'(2X,A)') " #inputLayer weight"
            WRITE(d_unit, format) weights(1:nGs(i),1:nhidneurons(1),1,i) 
            WRITE(d_unit,'(2X,A)') " #inputLayer bias"
            WRITE(d_unit,format) biases(1,1:nhidneurons(1),1,i)
            DO j = 1, nhidlayers-1
                WRITE(d_unit,'(2X,A,I0,A)') " #hiddenLayer_",j," weight"
                WRITE(d_unit, format) weights(1:nhidneurons(j),1:nhidneurons(j+1),j+1,i) 
                WRITE(d_unit,'(2X,A,I0,A)') " #hiddenLayer_",j," bias"
                WRITE(d_unit, format) biases(1,1:nhidneurons(j+1),j+1,i)
            END DO
            WRITE(d_unit,'(2X,A)') " #outputLayer weight"
            WRITE(d_unit, format) weights(1:nhidneurons(nhidlayers),1,nhidlayers+1,i)
            WRITE(d_unit,'(2X,A)') " #outputLayer bias"
            WRITE(d_unit, format) biases(1,1,nhidlayers+1,i)
        END DO
        CLOSE (d_unit)
    END SUBROUTINE

    SUBROUTINE update_mlff(file)
        USE nnType
        IMPLICIT NONE
        !NOTE: This updates V2 mlff.pyamff to the latest (it assumes you have already updated G1 and G2 lines)

        CHARACTER(len=20), OPTIONAL :: file
        ! Variables
        INTEGER :: nelement, nFPs
        INTEGER :: i, j, k
        INTEGER :: io, d_unit=22, s_unit=33
        CHARACTER(len=20) :: mlff_file, out_file
        CHARACTER(len=256) :: line
        CHARACTER(len=30) :: format

        IF(PRESENT(file)) THEN
            IF(LEN(TRIM(file)) .NE. 0) THEN
                mlff_file = file
            ELSE
                mlff_file = 'mlff.pyamff'
            END IF
        ELSE
            mlff_file = 'mlff.pyamff'
        END IF
        out_file = mlff_file

        CALL RENAME(TRIM(mlff_file),TRIM(ADJUSTL(mlff_file)) // ".v2")
        mlff_file = TRIM(ADJUSTL(mlff_file)) // ".v2"
        OPEN (newunit=s_unit, file=mlff_file, status='old',action='read',iostat=io)
        OPEN (newunit=d_unit, file=out_file, status='replace',action='write',iostat=io)
        !copy until we hit Model Structure
        DO
                READ(s_unit,'(A)',iostat=io) line
                IF (io .GT. 0) THEN ! error
                    WRITE(*,*) "Error while reading file '" // TRIM(ADJUSTL(mlff_file)) // "'"
                    CLOSE(s_unit)
                    CLOSE(d_unit)
                    ERROR STOP
                ELSE IF (io .LT. 0) THEN ! EOF
                    EXIT
                END IF
                IF (line .EQ. "#Model Structure") THEN
                    CLOSE(s_unit)
                    EXIT 
                END IF
                WRITE(d_unit,'(A)') TRIM(line)
        END DO
        WRITE(d_unit,'(A)') "#Energy Scaling Parameters"
        WRITE(d_unit,'(A)') scaler_type
        ! verify that this matches the python side
        ! WRITE(d_unit,'(1X,F16.8,4X,F16.8,4X,A)') slope,intercept,"#slope  intercept"
        WRITE(line,'(F16.8,4X,F16.8)') slope,intercept
        WRITE(d_unit,'(A,4X,A)') ADJUSTL(TRIM(line)),"#slope  intercept"
        WRITE(d_unit,'(A)') "#Model Structure"
        WRITE(d_unit,'(I0)') nhidlayers
        WRITE(format,'(A,I0,A)') "(",nhidlayers,"(I0,1X))" 
        WRITE(d_unit,format) nhidneurons
        WRITE(d_unit,'(A)') "#Model Parameters"
        ! Write NN parameters
        DO i = 1, nelements
            ! Write atom type
            WRITE (d_unit,'(A)') uniq_elements(i)
            WRITE(format,'(A)') '(5(4X,ES16.8))' ! per 5
            ! weights, biases
            WRITE(d_unit,'(2X,A)') " #inputLayer weight"
            WRITE(d_unit, format) weights(1:nGs(i),1:nhidneurons(1),1,i) 
            WRITE(d_unit,'(2X,A)') " #inputLayer bias"
            WRITE(d_unit,format) biases(1,1:nhidneurons(1),1,i)
            DO j = 1, nhidlayers-1
                WRITE(d_unit,'(2X,A,I0,A)') " #hiddenLayer_",j," weight"
                WRITE(d_unit, format) weights(1:nhidneurons(j),1:nhidneurons(j+1),j+1,i) 
                WRITE(d_unit,'(2X,A,I0,A)') " #hiddenLayer_",j," bias"
                WRITE(d_unit, format) biases(1,1:nhidneurons(j+1),j+1,i)
            END DO
            WRITE(d_unit,'(2X,A)') " #outputLayer weight"
            WRITE(d_unit, format) weights(1:nhidneurons(nhidlayers),1,nhidlayers+1,i)
            WRITE(d_unit,'(2X,A)') " #outputLayer bias"
            WRITE(d_unit, format) biases(1,1,nhidlayers+1,i)
        END DO
        CLOSE (d_unit)
    END SUBROUTINE

    SUBROUTINE write_mlff(file)
        !NOTE: If this is in vasp, we can just copy the file and rewrite the model information
        USE nnType
        IMPLICIT NONE

        ! Variables
        INTEGER :: nelement, nFPs
        INTEGER :: i, j, k , l
        INTEGER :: io, file_unit=22
        CHARACTER(len=20) :: format,mlff_file
        CHARACTER(len=20),OPTIONAL :: file
        CHARACTER(len=50) :: line
        INTEGER,DIMENSION(nelements) :: indices
        INTEGER,DIMENSION(2,nelements) :: numGs


        IF(PRESENT(file)) THEN
            IF(LEN(TRIM(file)) .NE. 0) THEN
                mlff_file = file
            ELSE
                mlff_file = 'mlff.pyamff'
            END IF
        ELSE
            mlff_file = 'mlff.pyamff'
        END IF

        !Probably some better way to format this but this is working so I'll update later
        OPEN (file_unit, FILE=mlff_file, status='replace')
        WRITE (file_unit,'(A)') "#Fingerprint type"
        WRITE (file_unit,'(A)') "BP"
        WRITE (file_unit,'(A)') "#Rmins"
        WRITE(format,'(A,I0,A)') "(",nelements,"(2X,A2))" 
        WRITE (file_unit,format) uniq_elements
        WRITE(format,'(A,I0,A)') "(",nelements,"(F8.2))" 
        DO i = 1, nelements
            WRITE (file_unit,format) rmins(i,:) 
        END DO

        ! If there is a #Cohesive tag
        IF (use_cohesive_energy) THEN
            WRITE(file_unit,'(A)') "#Cohesive"
            WRITE(format,'(A,I0,A)') "(",nelements,"(2X,A2))"  
            WRITE(file_unit,'(4X,A)') uniq_elements ! This might be wrong . . . should be coheEl so order might be wrong
            WRITE(format,'(A,I0,A)') "(",nelements,"(F8.2))" 
            WRITE(file_unit,format) coheEs
        END IF
        ! This is where things get tricky
        numGs = 0
        indices = 0
        DO i=1,nelements
            numGs(1,i) = fpParas(i)%g1_endpoint - fpParas(i)%g1_startpoint
            numGs(2,i) = fpParas(i)%tnFPs - numGs(1,i) ! what happens if no G1s?
        END DO

        IF (SUM(numGs(1,:)) .GT. 0) THEN
            WRITE(file_unit,'(A)') "#    type   number"
            WRITE(file_unit, '(6X,A)',advance="no") "G1"
            WRITE(format,'(A,I0,A)') "(",nelements,"(6X,I0))" 
            WRITE(file_unit,format) numGs(1,:)
            WRITE(file_unit,'(A,1X,A,4X,A,7X,A,5X,A,8X,A,14X,A)') "#  center","neighbor","eta","Rs","rcut","fpmin","fpmax"

            DO i=1,nelements
                IF(numGs(1,i) .GT. 0) THEN
                    l = 1
                    DO j=1,nelements
                        nFPS = fpParas(i)%g1s(j)%nFPS
                        IF( nFPS .NE. 0) THEN
                            DO k=1,nFPs !TODO: Save fpParas info so don't need to keep indexing, NOTE: Use ES16.8
                                WRITE (file_unit,'(2(6X,A2),3(1X,F8.3),2(1X,E18.10))') uniq_elements(i), uniq_elements(j), &
                                                                        fpParas(i)%g1s(j)%etas(k), fpParas(i)%g1s(j)%rss(k), &
                                                                        fpParas(i)%g1s(j)%r_cuts(k), fpminvs(l,i), fpmaxvs(l,i)
                                l = l + 1
                            END DO
                        END IF
                    END DO
                    indices(i) = l
                END IF
            END DO
        END IF

        IF (SUM(numGs(2,:)) .GT. 0) THEN
            WRITE(file_unit,'(A)') "#    type   number"
            WRITE(file_unit, '(6X,A)',advance="no") "G2"
            WRITE(format,'(A,I0,A)') "(",nelements,"(6X,I0))" 
            WRITE(file_unit,format) numGs(2,:)
            WRITE(file_unit,'(A,2(1X,A),3(4X,A),2(4X,A),8X,A,14X,A)') &
                         "#  center","neighbor1","neighbor2","eta","zeta","lambda","thetas","rcut","fpmin","fpmax"
            DO i=1,nelements
                IF(numGs(1,i) .GT. 0) THEN
                    DO j=1,nelements
                        DO k=1,nelements ! do I start with 1 or j?
                            nFPS = fpParas(i)%g2s(j,k)%nFPS
                            IF( nFPS .NE. 0) THEN
                                DO l=1,nFPs !NOTE: Use ES16.8 !'(2(6X,A2),3(1X,F8.3),2(1X,E20.12))'
                                    WRITE (file_unit,'(3(6X,A2),3X,5(1X,F8.3),2(1X,E18.10))') &
                                                        uniq_elements(i), uniq_elements(j),uniq_elements(k), &
                                                        fpParas(i)%g2s(j,k)%etas(l), fpParas(i)%g2s(j,k)%zetas(l), &
                                                        fpParas(i)%g2s(j,k)%lambdas(l), fpParas(i)%g2s(j,k)%theta_ss(l), &
                                                        fpParas(i)%g2s(j,k)%r_cuts(l), fpminvs(indices(i),i), fpmaxvs(indices(i),i)
                                                        
                                    indices(i) = indices(i) + 1
                                END DO
                            END IF
                        END DO
                    END DO
                END IF
            END DO
        END IF

        WRITE (file_unit,'(A)') "#MachineLearning model type"
        WRITE (file_unit,'(A)') "atomic-NN"
        WRITE (file_unit,'(A)') "#Activation function type"
        WRITE (file_unit,'(A)') actfuncId
        WRITE(file_unit,'(A)') "#Energy Scaling Parameters"
        WRITE(file_unit,'(A)') scaler_type
        WRITE(line,'(F16.8,4X,F16.8)') slope,intercept
        WRITE(file_unit,'(A,4X,A)') ADJUSTL(TRIM(line)),"#slope  intercept"
        ! WRITE(file_unit,'(F16.8,1X,F16.8,4X,A)') slope,intercept,"#slope intercept"
        WRITE(file_unit,'(A)') "#Model Structure"
        WRITE(file_unit,'(I0)') nhidlayers
        WRITE(format,'(A,I0,A)') "(",nhidlayers,"(I0,1X))" 
        WRITE(file_unit,format) nhidneurons
        WRITE(file_unit,'(A)') "#Model Parameters"
        ! Write NN parameters
        DO i = 1, nelements
            ! Write atom type
            WRITE (file_unit,'(A)') uniq_elements(i)
            WRITE(format,'(A)') '(5(4X,ES16.8))' ! per 5
            ! weights, biases
            WRITE(file_unit,'(2X,A)') " #inputLayer weight"
            WRITE(file_unit, format) weights(1:nGs(i),1:nhidneurons(1),1,i) 
            WRITE(file_unit,'(2X,A)') " #inputLayer bias"
            WRITE(file_unit,format) biases(1,1:nhidneurons(1),1,i)
            DO j = 1, nhidlayers-1
                WRITE(file_unit,'(2X,A,I0,A)') " #hiddenLayer_",j," weight"
                WRITE(file_unit, format) weights(1:nhidneurons(j),1:nhidneurons(j+1),j+1,i) 
                WRITE(file_unit,'(2X,A,I0,A)') " #hiddenLayer_",j," bias"
                WRITE(file_unit, format) biases(1,1:nhidneurons(j+1),j+1,i)
            END DO
            WRITE(file_unit,'(2X,A)') " #outputLayer weight"
            WRITE(file_unit, format) weights(1:nhidneurons(nhidlayers),1,nhidlayers+1,i)
            WRITE(file_unit,'(2X,A)') " #outputLayer bias"
            WRITE(file_unit, format) biases(1,1,nhidlayers+1,i)
        END DO
        CLOSE (file_unit)
    END SUBROUTINE
END MODULE
