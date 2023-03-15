MODULE fpCalc
    USE nlist
    USE fbp
    USE fpType
    USE atomsProp
    USE normalize, only: normalizeParas, loadnormalizeParas
    USE nnType !, only: natoms_arr, nGs, fpminvs, fpmaxvs, diffs, magnitude, interceptScale, coheEs, use_cohesive_energy, atom_idx
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: read_fpParas, read_mlffParas, calcg2s, calcfps, cleanup, read_mlff, calcg1s
    PUBLIC :: load_default_mlff, load_default_mlff_gr, atomsCleanup
    TYPE (fingerprints),DIMENSION(:),ALLOCATABLE,SAVE :: fpParas
    DOUBLE PRECISION, SAVE :: max_rcut   !fetch rcut from fpParas
    CHARACTER*2, DIMENSION(:), ALLOCATABLE :: uniq_elements
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: rmins
    INTEGER :: nG1, nG2    
    CONTAINS

!!!!!!-----------------------------------------------------------------------------------!
!!!!!! loopfps: calculate fingerprints for each images based on type of fingerprints
!!!!!!           i.e., G1(H, H)
!!!!!!          symbols: [0, 1, 1, 1]
!!!!!!-----------------------------------------------------------------------------------!
    SUBROUTINE read_fpParas(filename, nelement,coeh)

        ! inputs
        CHARACTER*20 :: filename
        INTEGER :: nelement
!f2py   INTENT(IN) :: filename
!f2py   INTENT(IN) :: nelement
        CHARACTER*2 :: G_type
        CHARACTER*2 :: center, neigh1, neigh2
        INTEGER*4 :: i, j, k, accN

        ! variables
        INTEGER :: cidx, tempidx, nidx1, nidx2 
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
                    ELSE
                        ALLOCATE(fpParas(i)%g2s(j,k)%etas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%gammas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%lambdas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%zetas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%r_cuts(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%theta_ss(nFPs))
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
        READ (12,*) !skip third line
        READ (12,*) !skip third line
        READ (12,*) !skip third line
        READ (12,*) !skip third line
        IF (nG1 .GT. 0) THEN
            READ (12,*)
            DO i = 1, nG1
                READ (12,*) G_type, center, neigh1, eta, Rs, rcut
                ! ideally have dynamic rcut for different interactions
                IF (rcut .GT. max_rcut) THEN
                    max_rcut = rcut
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
                    max_rcut = rcut
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
            !write(*,*) 'Parameter reading done'
        END IF
        !print*, 'maxcutoff',max_cutoff(1), coeh
        !CLOSE(12)
        !
        !Read mlff.pyamff to obtain fpminvs, and fpmaxvs
        !CALL loadnormalizeParas(nAtoms, nelement, MAX_FPS, symbols, uniq_element)
        !ALLOCATE(magnitude(MAX_FPS, MAXVAL(natoms_arr), nelement))
        !ALLOCATE(interceptScale(MAX_FPS, MAXVAL(natoms_arr), nelement))
        !Store magnitude interceptScale in memory
        !CALL normalizeParas(nelement)

    END SUBROUTINE

    SUBROUTINE read_mlffParas(nAtoms, nelement, MAX_FPS, atomicNumbers, uniqueNrs) BIND(C,name='read_mlffParas')
        !This subroutine is for EON-PyAMFF interface.
        USE, INTRINSIC :: iso_c_binding
        ! Inputs
        INTEGER(c_long) :: nAtoms
        INTEGER(c_int) :: nelement, MAX_FPS
        INTEGER(c_int), DIMENSION(nAtoms) :: atomicNumbers
        INTEGER(c_int), DIMENSION(nelement) :: uniqueNrs
        INTEGER, DIMENSION(nAtoms) :: symbols
        CHARACTER*2, DIMENSION(nelement) :: uniq_element

        ! Variables
        INTEGER :: numGs, nFPs, cidx, nidx, nidx1, nidx2, tempidx
        INTEGER :: i, j, k, accN, currIndex, m, n
        INTEGER, DIMENSION(nelement) :: fprange_idx 
        DOUBLE PRECISION :: djunk, fpmin, fpmax, eta, Rs, rcut, thetas, zeta, lambda
        CHARACTER*3 :: center, neigh1, neigh2, cjunk
        CHARACTER(LEN=30) :: G_type, line
        CHARACTER*3, DIMENSION(92) :: elementArray

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
        ALLOCATE(rmins(nelement, nelement))
        rmins=0.0

        DO i = 1, nAtoms
            DO j = 1, nelement
                IF (atomicNumbers(i) == uniqueNrs(j)) THEN
                    symbols(i) = j
                END IF
            END DO
        END DO

        ! Outputs
        !DOUBLE PRECISION, DIMENSION(MAX_FPS, nelement) :: fpminvs, fpmaxvs, diffs

        ALLOCATE (uniq_elements(nelement))
        ALLOCATE(natoms_arr(nelement))
        ALLOCATE(nGs(nelement))
        ALLOCATE(fpminvs(MAX_FPS, nelement))
        ALLOCATE(fpmaxvs(MAX_FPS, nelement))
        ALLOCATE(diffs(MAX_FPS, nelement))
        
        nGs = 0
        fpminvs = 0
        fpmaxvs = 0
        diffs = 0
        fprange_idx=0

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
        OPEN (11, FILE='mlff.pyamff', status='old')
        !print*, 'reading mlparas'
        READ (11,*) !skip #Fingerprint type
        READ (11,*) !fp_type
        READ (11,*) !Rmins
        READ (11,*) uniq_element
        uniq_elements=uniq_element
        !print*, 'reading Rmins'
        DO i = 1, nelement
            READ (11,*) (rmins(i,j), j = 1, nelement)
        END DO

        ! Mai edit
        !!!!!! If there is a #Cohesive tag
        ! Read the cohesives energies into the array
        !print*, 'Counting G1s'
        DO WHILE (line .NE. "#MachineLearning")
            READ (11,*) line! can be #  type or #Cohesive 
            IF (line .EQ. "#Cohesive") THEN
                use_cohesive_energy = .TRUE.
                READ (11,*) coheElement
                READ (11,*) temp_coheEs
                DO m = 1, nelement
                    n = find_loc_char(coheElement, uniq_element(m), nelement)
                    coheEs(m) = temp_coheEs(n)
                END DO
                READ (11,*) line! skip # type 
            END IF
            IF (line .EQ. "#MachineLearning") GOTO 40
            READ (11,*) G_type, numGs
            IF (G_type .EQ. 'G1') THEN
                nG1 = numGs
                READ (11,*) !skip # center neighbor ...
                DO i = 1, nG1
                    READ (11,*) center, neigh1, djunk, djunk, djunk, djunk, djunk
                    cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                    nidx = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                    nGs(cidx) = nGs(cidx) + 1
                    fpParas(cidx)%g1s(nidx)%nFPs = fpParas(cidx)%g1s(nidx)%nFPs + 1
                END DO
            ELSE IF (G_type .EQ. 'G2') THEN
                nG2 = numGs
                READ (11,*) !skip # center neighbor ...
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
            END IF
        END DO

  40    REWIND 11
        !print*, 'allocating fpParas'
        ! ALLOCATE fpParas
        DO i = 1, nelement
            DO j = 1, nelement
                nFPs = fpParas(i)%g1s(j)%nFPs
                ALLOCATE(fpParas(i)%g1s(j)%etas(nFPs))
                ALLOCATE(fpParas(i)%g1s(j)%rss(nFPs))
                ALLOCATE(fpParas(i)%g1s(j)%r_cuts(nFPs))
                DO k = 1, nelement
                    nFPs = fpParas(i)%g2s(j,k)%nFPs
                    !IF (nFPs .GT. 0) THEN
                    ALLOCATE(fpParas(i)%g2s(j,k)%etas(nFPs))
                    ALLOCATE(fpParas(i)%g2s(j,k)%gammas(nFPs))
                    ALLOCATE(fpParas(i)%g2s(j,k)%lambdas(nFPs))
                    ALLOCATE(fpParas(i)%g2s(j,k)%zetas(nFPs))
                    ALLOCATE(fpParas(i)%g2s(j,k)%r_cuts(nFPs))
                    ALLOCATE(fpParas(i)%g2s(j,k)%theta_ss(nFPs))
                    !END IF
                END DO
            END DO
        END DO

        READ (11,*) !skip #Fingerprint type
        READ (11,*) !fp_type 
        READ (11,*) !Rmins
        READ (11,*) !element type
        DO i = 1, nelement
            READ (11,*)
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
        READ (11,*) G_type, numGs
        !print*, 'G_type', G_type, numGs
        IF (nG1 .GT. 0) THEN
            READ (11,*) !skip center ...
            DO i = 1, nG1
                READ (11,*) center, neigh1, eta, Rs, rcut, fpmin, fpmax
                ! ideally have dynamic rcut for different interactions
                IF (rcut .GT. max_rcut) THEN
                    IF (max_rcut > 0.d0) THEN
                        PRINT *, 'Error: multiple rcuts are not supported yet. Use a single rcut value.'
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
                READ (11,*) G_type, numGs
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
                        PRINT *, 'Error: multiple rcuts are not supported yet. Use a single rcut value.'
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
        !Store magnitude interceptScale in memory
        CALL normalizeParas(nelement)
        !print*, 'reading mlParas Done'
    END SUBROUTINE

    SUBROUTINE read_mlff(nelement, MAX_FPS, filename, seedval)
    !This subroutine is for ASE Fortran calculator and this is designed to be called only once. 
    !Unlike read_mlffParas, this subroutine does not contain atomic information since 
    !each image can have different number of atoms. 
        ! Inputs
        !INTEGER :: nAtoms
        INTEGER :: nelement, MAX_FPS
        !Optional inputs
        CHARACTER(*), OPTIONAL :: filename
        INTEGER, OPTIONAL :: seedval
        ! Variables
        INTEGER :: numGs, nFPs, cidx, nidx, nidx1, nidx2, tempidx
        INTEGER :: i, j, k, accN, currIndex, m, n
        INTEGER, DIMENSION(nelement) :: fprange_idx 
        DOUBLE PRECISION :: djunk, fpmin, fpmax, eta, Rs, rcut, thetas, zeta, lambda
        CHARACTER*3 :: center, neigh1, neigh2, cjunk
        CHARACTER(LEN=30) :: G_type, line
        CHARACTER*3, DIMENSION(92) :: elementArray

        INTEGER :: myid, ios
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: flatten_inweights, flatten_hidweights
        CHARACTER(LEN=30) :: model_type
        CHARACTER*2 :: atom_type

        CHARACTER(len=2), DIMENSION(nelement) :: coheElement
        DOUBLE PRECISION, DIMENSION(nelement) :: temp_coheEs

        ALLOCATE(coheEs(nelement))
        use_cohesive_energy = .FALSE.
        max_rcut=0.d0
        !filename = 'fpParas.dat'
        nG1 = 0
        nG2 = 0
        ALLOCATE(rmins(nelement, nelement))
        rmins = 0.0

        ALLOCATE (uniq_elements(nelement))
        ALLOCATE(nGs(nelement))
        ALLOCATE(fpminvs(MAX_FPS, nelement))
        ALLOCATE(fpmaxvs(MAX_FPS, nelement))
        ALLOCATE(diffs(MAX_FPS, nelement))
        
        nGs = 0
        fpminvs = 0
        fpmaxvs = 0
        diffs = 0
        fprange_idx=0

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

        IF (PRESENT(filename)) THEN 
          OPEN(11, FILE=filename, status='old')
        ELSE  
          OPEN (11, FILE='mlff.pyamff', status='old')
        END IF
        
        READ (11,*) !skip #Fingerprint type
        READ (11,*) !fp_type 
        READ (11,*) !Rmins
        READ (11,*) uniq_elements
        !print*, 'reading Rmins'
        DO i = 1, nelement
            READ (11,*) (rmins(i,j), j = 1, nelement)
        END DO

        ! If there is a #Cohesive tag

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
                READ (11,*) line ! skip # type 
            END IF
            IF (line .EQ. "#MachineLearning") GOTO 40
            READ (11,*) G_type, numGs
            IF (G_type .EQ. 'G1') THEN
                nG1 = numGs
                READ (11,*) ! skip # center neighbor ...
                DO i = 1, nG1
                    READ (11,*) center, neigh1, djunk, djunk, djunk, djunk, djunk
                    cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                    nidx = find_loc_char(uniq_elements, neigh1, SIZE(uniq_elements))
                    nGs(cidx) = nGs(cidx) + 1
                    fpParas(cidx)%g1s(nidx)%nFPs = fpParas(cidx)%g1s(nidx)%nFPs + 1
                END DO

            ELSE IF (G_type .EQ. 'G2') THEN
                nG2 = numGs
                READ (11,*) ! skip # center neighbor ...
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
            END IF
        END DO

  40    REWIND 11

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
                        fpParas(i)%g2s(j,k)%etas = 0.0
                        fpParas(i)%g2s(j,k)%gammas = 0.0
                        fpParas(i)%g2s(j,k)%lambdas = 0.0
                        fpParas(i)%g2s(j,k)%zetas = 0.0
                        fpParas(i)%g2s(j,k)%r_cuts = 0.0
                        fpParas(i)%g2s(j,k)%theta_ss = 0.0
                    ENDIF
                    !print*, 'fpParas', fpParas(i)%g2s(j,k)%etas
                END DO
            END DO
        END DO
   
        READ (11,*) !skip #Fingerprint type
        READ (11,*) !fp_type 
        READ (11,*) !Rmins
        READ (11,*) !element type
        DO i = 1, nelement
            READ (11,*)  !Rmins
        END DO

         ! Mai edit: Read the cohesive Energy
         !!!!!!!!! If there is a #Cohesive tag
         ! read the cohesive energies to the array

         !print*, 'reading g1'
        READ (11,*) line ! can be #  type or #Cohesive
        IF (line .EQ. "#Cohesive") THEN
            READ (11,*)
            READ (11,*)
            READ (11,*) !skip # type
        END IF
        READ (11,*) G_type, numGs
        IF (nG1 .GT. 0) THEN
            READ (11,*) !skip center ...
            DO i = 1, nG1
                READ (11,*) center, neigh1, eta, Rs, rcut, fpmin, fpmax
                ! ideally have dynamic rcut for different interactions
                IF (rcut .GT. max_rcut) THEN
                    IF (max_rcut > 0.d0) THEN
                        PRINT *, 'Error: multiple rcuts are not supported yet. Use a single rcut value.'
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

        IF (nG2 .GT. 0) THEN
            IF (nG1 .GT. 0) THEN
                READ (11,*) !skip type ...
                READ (11,*) G_type, numGs
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
                        PRINT *, 'Error: multiple rcuts are not supported yet. Use a single rcut value.'
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

        ! Read NN parameters 
        ! Set common variables
        nelements = nelement
        !total_natoms = natoms

        ! Read weights and biases
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
        ALLOCATE (in_weights(MAXVAL(nGs),nhidneurons(1),nelements))
        ALLOCATE (in_biases(nhidneurons(1),nelements))
        ALLOCATE (hid_weights(MAXVAL(nhidneurons),MAXVAL(nhidneurons),nhidlayers-1,nelements))
        ALLOCATE (hid_biases(MAXVAL(nhidneurons),nhidlayers-1,nelements))
        ALLOCATE (out_weights(MAXVAL(nhidneurons),1,nelements))
        ALLOCATE (out_biases(nelements))
       
        READ (11,*) !Skip #Model Parameters
        DO i = 1, nelements
            ! Read atom type skipping the first symbol #.
            READ (11,*) atom_type
            ! Find index of corresponding atom_type
            myid=find_loc_char(uniq_elements, atom_type, SIZE(uniq_elements))
            ! Input weights, biases
            !READ (11,*) flatten_inweights(1:nGs(myid)*nhidneurons(1))
            READ (11,*,IOSTAT=ios) flatten_inweights(1:nGs(myid)*nhidneurons(1))
            IF (ios == 59) THEN
                ! TODO: 11/16/22: This is broken currently. Each core generates different params. It should be fixed
                CALL initialize_NN_params(uniq_elements, seedval) 
                !PRINT *, 'Error: NN parameters are missing in your file: ', filename
                !STOP
            ELSE
                in_weights(1:nGs(myid),1:nhidneurons(1),myid)=&
                reshape(flatten_inweights(1:nGs(myid)*nhidneurons(1)),(/nGs(myid),nhidneurons(1)/))
                READ (11,*) in_biases(1:nhidneurons(1),myid)
                ! Hidden weights, biases
                DO j = 1, nhidlayers-1
                    READ (11,*) flatten_hidweights(1:nhidneurons(j)*nhidneurons(j+1))
                    hid_weights(1:nhidneurons(j),1:nhidneurons(j+1),j,myid) = &
                    reshape(flatten_hidweights(1:nhidneurons(j)*nhidneurons(j+1)),(/nhidneurons(j),nhidneurons(j+1)/))
                    READ (11,*) hid_biases(1:nhidneurons(j+1),j,myid)
                END DO
                ! Out weights, biases
                READ (11,*) out_weights(1:nhidneurons(nhidlayers),1,myid)
                READ (11,*) out_biases(myid)
                IF (i == nelements) READ (11,*) !skip command
            END IF
        END DO
        READ (11,'(a)') scaler_type
        READ (11,*) slope, intercept
        CLOSE (11)
    END SUBROUTINE

    SUBROUTINE calcfps(nAtoms, pos_car, cell, symbols, max_fps, nelement,forceEngine, &
                       fps, dfps, neighs_incell, num_neigh)

        ! Parameters
        INTEGER, PARAMETER :: MAX_NEIGHS  = 100

        ! input values
        INTEGER :: nAtoms, nelement
        INTEGER :: forceEngine
        INTEGER*4, DIMENSION(nAtoms) :: symbols
        DOUBLE PRECISION, DIMENSION(nAtoms,3) :: pos_car
        DOUBLE PRECISION, DIMENSION(3,3) :: cell
        INTEGER :: MAX_FPS  !find max number of fingerprints, element-based
        !!!CHARACTER*2, INTENT(IN), DIMENSION(nelement) :: uniq_elements
!f2py   INTENT(IN) :: pos_car
!f2py   INTENT(IN) :: cell
!f2py   INTENT(IN) :: symbols
!f2py   INTENT(IN) :: forceEngine
!f2py   INTENT(IN) :: nelement
!!!!f2py    INTENT(IN) :: nAtoms, nelement
!f2py   INTENT(IN) :: MAX_FPS  !find max number of fingerprints, element-based

        ! variables for nlist.calc
        INTEGER :: j, k, l
        INTEGER :: npairs, maxneighs, npairs_incell, tncells, max_npairs
        INTEGER, DIMENSION(2) :: num_pairs
        INTEGER, DIMENSION(3) :: ncells
        DOUBLE PRECISION, DIMENSION(3, 3) :: supercell

        ! variables for calcTriplet
        INTEGER :: ntriplets
        DOUBLE PRECISION, DIMENSION(3, nAtoms*MAX_NEIGHS*MAX_NEIGHS) :: angles
        INTEGER, DIMENSION(3, nAtoms*MAX_NEIGHS*MAX_NEIGHS) :: angle_indices
        DOUBLE PRECISION, DIMENSION(3, 3) :: vectsigns
        INTEGER, DIMENSION(nAtoms*MAX_NEIGHS) :: pair_offset   ! offset of neigh atom

        ! output values
        DOUBLE PRECISION, DIMENSION(nAtoms, MAX_FPS) :: fps
        DOUBLE PRECISION, DIMENSION(nAtoms, MAX_NEIGHS, 3, MAX_FPS) :: dfps
        INTEGER, DIMENSION(nAtoms, MAX_NEIGHS) :: neighs
        INTEGER, DIMENSION(nAtoms, MAX_NEIGHS) :: neighs_incell
        INTEGER, DIMENSION(nAtoms) :: num_neigh
!f2py   INTENT(OUT) :: fps
!f2py   INTENT(OUT) :: dfps
!f2py   INTENT(OUT) :: neighs
!f2py   INTENT(OUT) :: num_neigh
        
        !print *, 'calc fps is called'
        vectsigns = reshape((/1, 1, -1, 1, -1, -1, 1, -1, 1/), [3,3])

        CALL calcCellNum(nAtoms, pos_car, cell, max_rcut, ncells)

        tncells = ncells(1) * ncells(2) * ncells(3)
        tnAtoms = nAtoms * tncells
        ALLOCATE(pool_pos_car(tnAtoms,3))
        ALLOCATE(pool_ids(tnAtoms))
        ALLOCATE(supersymbols(tnAtoms))
        CALL genSupercell(nAtoms, pos_car,symbols, cell, ncells, tncells,  max_rcut, supercell)
        max_npairs = tnAtoms*MAX_NEIGHS
        ALLOCATE(pairs(max_npairs))
        ALLOCATE(pair_info(tnAtoms,MAX_NEIGHS))
        ALLOCATE(num_eachpair(tnAtoms))
        ALLOCATE(pair_indices(2, max_npairs))
        ALLOCATE(pair_global_indices(2, max_npairs))
        ALLOCATE(unitvects_pair(2, max_npairs, 3))
        ALLOCATE(gvects(max_npairs))
        ALLOCATE(tneighs(tnAtoms, MAX_NEIGHS))
        ALLOCATE(tneighs_incell(tnAtoms, MAX_NEIGHS))
        ALLOCATE(tnum_neigh(tnAtoms))

        CALL calcNlist(tnAtoms, MAX_NEIGHS, pool_pos_car, supercell, symbols,  max_rcut, nelement, &
                       forceEngine, rmins, tncells, &
                       npairs_incell, npairs, tnum_neigh, tneighs, tneighs_incell)
        num_neigh = tnum_neigh(1:nAtoms)
        maxneighs = MAXVAL(num_neigh)
        neighs = tneighs(1:nAtoms, :)
        neighs_incell = tneighs_incell(1:nAtoms, :)
      
        CALL calcg1s(nAtoms, npairs_incell, symbols, maxneighs, MAX_NEIGHS, max_fps, neighs, &
                     num_neigh, max_npairs, &
                     fps, dfps)

        !print *, 'dfps g1s'
        !DO j = 1, natoms
        !  DO k = 1, 3
        !    DO l = 1,num_neigh(j)+1
        !       IF (l==1) THEN
        !         print *, j-1, j-1, k-1, dfps(j, l, k, :)
        !       ELSE
        !         IF (neighs(j, l-1) <= natoms) THEN
        !           print *, j-1, neighs(j,l-1)-1, k-1, dfps(j, l, k, :)
        !         END IF
        !       ENDIF
        !    ENDDO
        !  ENDDO
        !ENDDO 

        IF (nG2 .GT. 0) THEN
            CALL calcTriplet(nAtoms, tnAtoms, MAX_NEIGHS, max_rcut, npairs, &
                            ntriplets, angles, angle_indices)
            
            CALL calcg2s(nAtoms, npairs, MAX_NPAIRS, ntriplets, maxneighs, MAX_NEIGHS, max_fps, &
                        symbols, max_rcut, neighs, num_neigh, &
                        unitvects_pair(1, :, :), &
                        vectsigns, angles(:, 1:ntriplets), angle_indices(:, 1:ntriplets), &
                        fps, dfps)
 
        END IF

        !print *, 'dfps main'
        !DO j = 1, natoms
        !  DO k = 1, 3
        !    DO l = 1,num_neigh(j)+1
        !       IF (l==1) THEN
        !         print *, j-1, j-1, k-1, dfps(j, l, k, :)
        !       ELSE
        !         IF (neighs(j, l-1) <= natoms) THEN
        !           print *, j-1, neighs(j,l-1)-1, k-1, dfps(j, l, k, :)
        !         END IF
        !       ENDIF
        !    ENDDO
        !  ENDDO
        !ENDDO
        !write(*,*) 'Fingerprint calculation is done'

    END SUBROUTINE

    SUBROUTINE calcg1s(nAtoms, npairs, symbols, maxneighs, MAX_NEIGHS, max_fps, neighs, &
                       num_neigh, max_npairs, &
                       fps, dfps)

        ! input values
        INTEGER :: nAtoms, npairs, maxneighs, MAX_NEIGHS, max_fps, max_npairs 
        INTEGER, DIMENSION(nAtoms) :: symbols
        INTEGER, DIMENSION(nAtoms, MAXNEIGHS) :: neighs
        INTEGER :: nfps
        INTEGER, DIMENSION(nAtoms) :: num_neigh

        ! variables
        INTEGER :: i, k, j
        INTEGER :: center_1, index_1, sp_1, sp_2, ep_1, ep_2
        INTEGER :: center_2, index_2, o_index_1, o_index_2
        INTEGER, DIMENSION(1) :: neighloc_1, neighloc_2
        TYPE (fingerprintsData), DIMENSION(2) :: fpdata_temp
        ! outputs
        DOUBLE PRECISION, DIMENSION(nAtoms, MAX_FPS) :: fps
        DOUBLE PRECISION, DIMENSION(nAtoms, MAX_NEIGHS, 3, MAX_FPS) :: dfps

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
      
            IF ((index_1<=nAtoms).AND.(index_2<=nAtoms)) THEN
                sp_1 = fpParas(center_1)%g1s(center_2)%startpoint
                ep_1 = fpParas(center_1)%g1s(center_2)%endpoint
                nFPs = fpParas(center_1)%g1s(center_2)%nFPs
                IF (nFPs .EQ. 0) THEN
                    CYCLE
                END IF
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
                    dfps(index_1, neighloc_1(1), k,sp_2:ep_2) = &
                    dfps(index_1, neighloc_1(1), k,sp_2:ep_2) + &
                    fpdata_temp(2)%pre_dgdxs(sp_2:ep_2) * &
                    unitvects_pair(2, i, k)
                    !Center on index_1 but w.r.t. index_2
                    dfps(index_2, neighloc_2(1), k,sp_1:ep_1) = &
                    dfps(index_2, neighloc_2(1), k,sp_1:ep_1) + &
                    fpdata_temp(1)%pre_dgdxs(sp_1:ep_1) * &
                    unitvects_pair(1, i, k)
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
                neighloc_2 = find_loc_int(neighs(o_index_2,:), index_1, SIZE(neighs(o_index_2,:))) +1
                neighloc_1 = find_loc_int(neighs(index_1,:), index_2, SIZE(neighs(index_1,:))) + 1

                !derivative
                DO k = 1, 3
                    !center on itself
                    dfps(index_1, 1, k, sp_1:ep_1) = &
                    dfps(index_1, 1, k, sp_1:ep_1) - &
                    fpdata_temp(1)%pre_dgdxs(sp_1:ep_1) * &
                    unitvects_pair(1, i, k)

                    !Center on index_2 but w.r.t. index_1
                    dfps(index_1, neighloc_1(1), k,sp_2:ep_2) = &
                    dfps(index_1, neighloc_1(1), k,sp_2:ep_2) + &
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

                DO k = 1, 3
                     !center on itself
                    dfps(index_2, 1, k,sp_2:ep_2) = &
                    dfps(index_2, 1, k,sp_2:ep_2) - &
                    fpdata_temp(2)%pre_dgdxs(sp_2:ep_2) * &
                    unitvects_pair(2, i, k)
                    !Center on index_2 but w.r.t. index_1
                    IF (gvects(i)>0) THEN
                        dfps(o_index_1, neighloc_1(1), k,sp_2:ep_2) = &
                        dfps(o_index_1, neighloc_1(1), k,sp_2:ep_2) + &
                        fpdata_temp(2)%pre_dgdxs(sp_2:ep_2) * &
                        unitvects_pair(2, i, k)
                    END IF
                END DO
                DEALLOCATE(fpdata_temp(2)%gs)
                DEALLOCATE(fpdata_temp(2)%pre_dgdxs)
            END IF
        END DO

    END SUBROUTINE

    SUBROUTINE calcg2s(natoms, npairs, MAX_NPAIRS, ntriplets, maxneighs, MAX_NEIGHS, max_fps, symbols,&
                       r_cut, neighs, num_neigh, unitvects, vectsigns, angles, angles_indices, fps, dfps)
        ! Define inputs
        INTEGER :: nAtoms, npairs, ntriplets, maxneighs, max_fps
        INTEGER :: MAX_NPAIRS, MAX_NEIGHS
        DOUBLE PRECISION :: r_cut
        INTEGER, DIMENSION(nAtoms) :: symbols
        INTEGER, DIMENSION(nAtoms, MAX_NEIGHS) :: neighs
        !TYPE (fingerprints), DIMENSION(nelement) :: fps
        !DOUBLE PRECISION, DIMENSION(MAX_NPAIRS) :: pairs
        DOUBLE PRECISION, DIMENSION(MAX_NPAIRS, 3) :: unitvects
        DOUBLE PRECISION, DIMENSION(3,3) :: vectsigns
        !INTEGER, DIMENSION(2, MAX_NPAIRS) :: pair_indices
        !INTEGER, DIMENSION(2, MAX_NPAIRS) :: pair_global_indices
        DOUBLE PRECISION, DIMENSION(3, ntriplets) :: angles
        INTEGER, DIMENSION(3, ntriplets) :: angles_indices
        INTEGER, DIMENSION(nAtoms) :: num_neigh

        ! Define variable
        INTEGER :: i, j, k
        INTEGER :: sp, ep, nFPs
        INTEGER, DIMENSION(1) :: neigh1_loc, neigh2_loc
        INTEGER, DIMENSION(3) :: center, pairIndex, csymbol
        INTEGER, DIMENSION(3,2) :: neighbor
        DOUBLE PRECISION, DIMENSION(3,3) :: uvects
        INTEGER, DIMENSION(2) :: currneighs
        DOUBLE PRECISION :: fexpt, fct
        DOUBLE PRECISION, DIMENSION(npairs) :: pair_sins, fcuts, fexp
        INTEGER, DIMENSION(3,3) :: indices, vertices, gcenter
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: f_thetas_1, f_thetas, f_rs
        DOUBLE PRECISION, DIMENSION(:, :, :), allocatable :: arm
        DOUBLE PRECISION, DIMENSION(:), allocatable :: gs

        ! Define outputs
        DOUBLE PRECISION, DIMENSION(nAtoms, MAX_FPS) :: fps
        DOUBLE PRECISION, DIMENSION(nAtoms, MAX_NEIGHS, 3, MAX_FPS) :: dfps

        !dfps = 0.0
        !fps = 0.0
        ! Loop over all angles and calculate common parts
        CALL commonparts(npairs, pairs, r_cut,&
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
            vertices(1, :) = center
            vertices(2, :) = [center(2), center(3), center(1)]
            vertices(3, :) = [center(3), center(1), center(2)]

            csymbol = symbols(center)
            neighbor(1,:) = [csymbol(2), csymbol(3)]
            neighbor(2,:) = [csymbol(1), csymbol(3)]
            neighbor(3,:) = [csymbol(1), csymbol(2)]

            fexpt = PRODUCT(fexp(pairIndex))
            fct   = PRODUCT(fcuts(pairIndex))

            gcenter(1, :) = [pair_global_indices(1, pairIndex(1)), &
                            pair_global_indices(2, pairIndex(1)), &
                            pair_global_indices(2, pairIndex(2))]
            gcenter(2, :) = [pair_global_indices(2, pairIndex(1)), &
                            pair_global_indices(2, pairIndex(2)), &
                            pair_global_indices(1, pairIndex(1))]
            gcenter(3, :) = [pair_global_indices(2, pairIndex(2)), &
                            pair_global_indices(1, pairIndex(1)), &
                            pair_global_indices(2, pairIndex(1))]

            IF ((gcenter(1, 1) > nAtoms) .AND. &
                (gcenter(1, 2) > nAtoms) .AND. &
                (gcenter(1, 3) > nAtoms)) THEN
               CYCLE
            END IF

            indices(1, :) = pairIndex
            indices(2, :) = pairIndex([3, 1, 2])  !j, k, i
            indices(3, :) = pairIndex([2, 3, 1])  !k, i, j

            DO j = 1, 3
                uvects = unitvects(indices(j, :), :)
                currneighs = neighbor(j,:)
                nFPs = fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%nFPs
                IF (nFPs == 0) THEN
                    CYCLE
                END IF

                sp = fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%startpoint
                ep = fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%endpoint

                ALLOCATE(gs(fpParas(csymbol(j))%tnFPs))
                ALLOCATE(arm(3,3, fpParas(csymbol(j))%tnFPs))
                ALLOCATE(f_thetas_1(nFPs))
                ALLOCATE(f_thetas(nFPs))
                ALLOCATE(f_rs(nFPs))
                arm = 0.0
                f_thetas = 0.0
                f_rs = 0.0
                CALL fg2(nFPs, &
                         fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%theta_ss, &
                         fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%etas, &
                         fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%zetas, &
                         fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%lambdas, &
                         angles(j, i), fexpt, fct, &
                         gs(sp:ep), &
                         f_thetas_1, f_thetas, f_rs)

                CALL fdg2(nFPs, &
                         fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%theta_ss, &
                         fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%etas, &
                         fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%zetas, &
                         fpParas(csymbol(j))%g2s(currneighs(1),currneighs(2))%lambdas, &
                         r_cut, &
                         angles(j, i), &
                         f_thetas_1, &
                         f_thetas, &
                         f_rs, &
                         pairs(indices(j,:)), &
                         fcuts(indices(j, :)), &
                         pair_sins(indices(j, :)), &
                         unitvects(indices(j, :), :), &
                         vectsigns(j,:), &
                         arm(1, :, sp:ep), &
                         arm(2, :, sp:ep), &
                         arm(3, :, sp:ep))
                DEALLOCATE(f_thetas_1)
                DEALLOCATE(f_thetas)
                DEALLOCATE(f_rs)

                IF (gcenter(j,1)<=nAtoms) THEN
                    fps(gcenter(j, 1), sp:ep) = fps(gcenter(j,1), sp:ep) + gs(sp:ep)
                    dfps(gcenter(j,1), 1, :, sp:ep) = &
                    dfps(gcenter(j,1), 1, :, sp:ep) + &
                    arm(1, :, sp:ep)
                    IF (gcenter(j,2) <= nAtoms) THEN
                        neigh1_loc = find_loc_int(neighs(vertices(j, 2),:), gcenter(j,1), SIZE(neighs(vertices(j, 2),:))) + 1
                        dfps(vertices(j, 2), neigh1_loc(1), :, sp:ep) = &
                        dfps(vertices(j, 2),neigh1_loc(1), :, sp:ep) + &
                        arm(2, :, sp:ep)
                    END IF
                    IF (gcenter(j,3) <= nAtoms) THEN
                        neigh2_loc = find_loc_int(neighs(vertices(j, 3),:), gcenter(j,1), SIZE(neighs(vertices(j, 3),:))) +1
                        dfps(vertices(j, 3), neigh2_loc(1), :, sp:ep) = &
                        dfps(vertices(j, 3), neigh2_loc(1), :, sp:ep) + &
                        arm(3, :, sp:ep)
                    END IF
                ELSE
                    IF (gcenter(j,2) <= nAtoms) THEN
                        neigh1_loc = find_loc_int(neighs(vertices(j, 2),:), gcenter(j,1), SIZE(neighs(vertices(j, 2),:))) + 1
                        dfps(vertices(j, 2), neigh1_loc(1), :, sp:ep) = &
                        dfps(vertices(j, 2),neigh1_loc(1), :, sp:ep) + &
                        arm(2, :, sp:ep)
                    END IF
                    IF (gcenter(j,3) <= nAtoms) THEN
                        neigh2_loc = find_loc_int(neighs(vertices(j, 3),:), gcenter(j,1), SIZE(neighs(vertices(j, 3),:))) +1
                        dfps(vertices(j, 3), neigh2_loc(1), :, sp:ep) = &
                        dfps(vertices(j, 3), neigh2_loc(1), :, sp:ep) + &
                        arm(3, :, sp:ep)
                    END IF
                END IF

                DEALLOCATE(gs)
                DEALLOCATE(arm)

            END DO
        END DO

        !print *, 'fps'
        !DO i = 1, natoms
        !    print *, i, fpdata(i)%gs
        !END DO

        !print *, 'dfps g2s'
        !DO i = 1, natoms
        !    DO j = 1, 3
        !        DO k = 1,num_neigh(i)+1
        !            print *, i, k, j
        !            print *, '         ', fpdata(i)%dgdxs(k)%dgdx(j,:)
        !            print *, '         ', dfps(i,k,j,:)
        !        END DO
        !    END DO
        !END DO
        !print *, 'done'
    END SUBROUTINE
        
    SUBROUTINE load_default_mlff(nelement, MAX_FPS, uniqElems, seedval)
    !Ref: AMP  
        !Inputs
        INTEGER, INTENT(IN) :: nelement, MAX_FPS, seedval
        CHARACTER*2, DIMENSION(nelement), INTENT(IN) :: uniqElems
        !Variables
        INTEGER :: i, j, k, nFPs, accN
        CHARACTER(LEN=30) :: model_type
        DOUBLE PRECISION, DIMENSION(4) :: g1_etas
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: flatten_inweights, flatten_hidweights

        print *, 'load default mlff subroutine is called'

        ALLOCATE(coheEs(nelement))
        use_cohesive_energy = .FALSE.

        nG1 = 0
        nG2 = 0
        ALLOCATE(rmins(nelement, nelement))
        rmins = 0.0

        ALLOCATE (uniq_elements(nelement))
        ALLOCATE(nGs(nelement))
        ALLOCATE(fpminvs(MAX_FPS, nelement))
        ALLOCATE(fpmaxvs(MAX_FPS, nelement))
        ALLOCATE(diffs(MAX_FPS, nelement))
        nGs = 0
        fpminvs = 0
        fpmaxvs = 0
        diffs = 0 

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
        
        !Set default values following mlff.pyamff order
        DO i=1, nelement
            !Load unique elements
            uniq_elements(i)=uniqElems(i)
            !Set minimum distance between elements to 0.01
            rmins(i,1:nelement)=0.01
        END DO  
        
        !Set number of fingerprints, g1 and g2s
        DO i=1, nelement
            DO j=1, nelement
                fpParas(i)%g1s(j)%nFPs=4
                nGs(i)=nGs(i)+4
                DO k=1, nelement
                    fpParas(i)%g2s(j,k)%nFPs=4
                    nGs(i)=nGs(i)+4
                END DO
            END DO
        END DO
      
        !Allocate fpParas and set deafult values 
        !(REF: make_default_symmetry_functions in /amp/amp/descriptor/gaussian.py)
        !Allocate first
        DO i=1, nelement
            DO j=1, nelement
                nFPs = fpParas(i)%g1s(j)%nFPs
                ALLOCATE(fpParas(i)%g1s(j)%etas(nFPs))
                ALLOCATE(fpParas(i)%g1s(j)%rss(nFPs))
                ALLOCATE(fpParas(i)%g1s(j)%r_cuts(nFPs))
                DO k=1, nelement
                    nFPs = fpParas(i)%g2s(j,k)%nFPs
                    IF (nFPs .GT. 0) THEN
                        ALLOCATE(fpParas(i)%g2s(j,k)%etas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%gammas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%lambdas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%zetas(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%r_cuts(nFPs))
                        ALLOCATE(fpParas(i)%g2s(j,k)%theta_ss(nFPs))
                        fpParas(i)%g2s(j,k)%etas = 0.0
                        fpParas(i)%g2s(j,k)%gammas = 0.0
                        fpParas(i)%g2s(j,k)%lambdas = 0.0
                        fpParas(i)%g2s(j,k)%zetas = 0.0
                        fpParas(i)%g2s(j,k)%r_cuts = 0.0
                        fpParas(i)%g2s(j,k)%theta_ss = 0.0
                    ENDIF
                END DO
            END DO
        END DO    
        
        !Now, set fpParas
        !Get etas from logspace function
        g1_etas=(/0.050,0.230,1.080,5.00/) !same as log_space(LOG10(0.05d0),LOG10(5.d0),4) in python
        !G1 fingerprint
        DO i=1, nelement
            DO j=1, nelement
                nFPs = fpParas(i)%g1s(j)%nFPs
                fpParas(i)%g1s(j)%fp_type='G1'
                fpParas(i)%g1s(j)%species1=uniq_elements(i)
                fpParas(i)%g1s(j)%species2=uniq_elements(j)
                fpParas(i)%g1s(j)%etas(1:nFPs)=g1_etas
                fpParas(i)%g1s(j)%rss(1:nFPs)=0.0
                fpParas(i)%g1s(j)%r_cuts(1:nFPs)=3.0 !TODO: default?
                fpParas(i)%tnFPs=fpParas(i)%tnFPs+4
                
                IF (j == 1) THEN 
                    fpParas(i)%g1s(j)%startpoint=1
                    fpParas(i)%g1s(j)%endpoint=fpParas(i)%g1s(j)%endpoint+4
                ELSE
                    fpParas(i)%g1s(j)%startpoint=fpParas(i)%g1s(j-1)%endpoint+1
                    fpParas(i)%g1s(j)%endpoint=fpParas(i)%g1s(j)%startpoint+fpParas(i)%g1s(j)%nFPs-1
                END IF
            END DO
        END DO

        !G2 fingerprint
        DO i=1, nelement
            fpParas(i)%g1_endpoint=fpParas(i)%tnFPs
            fpParas(i)%g2_startpoint=fpParas(i)%g1_endpoint+1
        END DO
        DO i=1, nelement
            DO j=1, nelement
                DO k=1, nelement
                    fpParas(i)%g2s(j,k)%fp_type='G2'
                    fpParas(i)%g2s(j,k)%species1=uniq_elements(i)
                    fpParas(i)%g2s(j,k)%species1_code=i
                    fpParas(i)%g2s(j,k)%species2=uniq_elements(j)
                    fpParas(i)%g2s(j,k)%etas = 0.005
                    fpParas(i)%g2s(j,k)%r_cuts = 3.0 !TODO: default?
                    fpParas(i)%g2s(j,k)%lambdas(1:2) = 1.0 !TODO: default?
                    fpParas(i)%g2s(j,k)%lambdas(3:4) = -1.0    
                    fpParas(i)%g2s(j,k)%theta_ss = 0.0 !TODO: default?
                    fpParas(i)%g2s(j,k)%zetas(1) = 1.0
                    fpParas(i)%g2s(j,k)%zetas(2) = 4.0
                    fpParas(i)%g2s(j,k)%zetas(3) = 1.0 
                    fpParas(i)%g2s(j,k)%zetas(4) = 4.0
                END DO
            END DO
            fpParas(i)%tnFPs=fpParas(i)%tnFPs+4
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

        !Set max_rcut
        max_rcut=3.0
 
        !Make sure nelements are known
        nelements=nelement

        !Load default NN 
        model_type='BP'
        actfuncId='tanh' 
        nhidlayers=2
        ALLOCATE(nhidneurons(nhidlayers))
        !nhidneurons=(/5,5/)
        nhidneurons=(/10,10/)
        ALLOCATE (in_weights(MAXVAL(nGs),nhidneurons(1),nelements))
        ALLOCATE (in_biases(nhidneurons(1),nelements))
        ALLOCATE (hid_weights(MAXVAL(nhidneurons),MAXVAL(nhidneurons),nhidlayers-1,nelements))
        ALLOCATE (hid_biases(MAXVAL(nhidneurons),nhidlayers-1,nelements))
        ALLOCATE (out_weights(MAXVAL(nhidneurons),1,nelements))
        ALLOCATE (out_biases(nelements))
 
        !Initialize random NN parameters
        CALL initialize_NN_params(uniqElems, seedval)
        
        !Set scaler type, slope and intercept
        scaler_type='NoScaler'
        slope=1.0
        intercept=0.0

    END SUBROUTINE
    
    SUBROUTINE load_default_mlff_gr(nelement, MAX_FPS, uniqElems, seedval, pos_car, nAtoms)
    ! Calculate g(r) and g(theta) of input configuration(s) and
    ! generate fingerprint parameters and neural nets
        !Inputs
        INTEGER, INTENT(IN) :: nelement, MAX_FPS, seedval, nAtoms
        CHARACTER*2, DIMENSION(nelement), INTENT(IN) :: uniqElems
        DOUBLE PRECISION, DIMENSION(nAtoms,3), INTENT(IN) :: pos_car
        !Variables
        INTEGER :: i, j, k, nFPs
        CHARACTER(LEN=30) :: model_type
        DOUBLE PRECISION, DIMENSION(4) :: g1_etas
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: flatten_inweights, flatten_hidweights

        print *, 'load default mlff gr subroutine is called'

        ALLOCATE(coheEs(nelement))
        use_cohesive_energy = .FALSE.

        nG1 = 0
        nG2 = 0
        ALLOCATE(rmins(nelement, nelement))
        rmins = 0.0

        ALLOCATE (uniq_elements(nelement))
        ALLOCATE(nGs(nelement))
        ALLOCATE(fpminvs(MAX_FPS, nelement))
        ALLOCATE(fpmaxvs(MAX_FPS, nelement))
        ALLOCATE(diffs(MAX_FPS, nelement))
        nGs = 0
        fpminvs = 0
        fpmaxvs = 0
        diffs = 0

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

        ! TODO ----------------------------------------------
        ! Calculate g(r) of input configuration (pos_car)

        ! Calculate g(theta) of input configuration

        ! Set number of fingerprints, g1 (fpParas%g1(s)%nFPs), g2s (fpParas%g2(s)%nFPs, and nGs of each element 

        ! ---------------------------------------------------

        !Make sure nelements are known
        nelements=nelement

        !Load default NN 
        model_type='BP'
        actfuncId='tanh'
        nhidlayers=2
        ALLOCATE(nhidneurons(nhidlayers))
        !nhidneurons=(/5,5/)
        nhidneurons=(/10,10/)
        ALLOCATE (in_weights(MAXVAL(nGs),nhidneurons(1),nelements))
        ALLOCATE (in_biases(nhidneurons(1),nelements))
        ALLOCATE (hid_weights(MAXVAL(nhidneurons),MAXVAL(nhidneurons),nhidlayers-1,nelements))
        ALLOCATE (hid_biases(MAXVAL(nhidneurons),nhidlayers-1,nelements))
        ALLOCATE (out_weights(MAXVAL(nhidneurons),1,nelements))
        ALLOCATE (out_biases(nelements))

        !Initialize random NN parameters
        CALL initialize_NN_params(uniqElems, seedval)

        !Set scaler type, slope and intercept
        scaler_type='NoScaler'
        slope=1.0
        intercept=0.0

    END SUBROUTINE load_default_mlff_gr
    
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
            DO k=1, nGs(i)
                CALL gather_grnd(in_weights(k,1:nhidneurons(1),i),nhidneurons(1),seedval)
            END DO
            !print *, 'in_weights=', in_weights(1:nGs(i),1:nhidneurons(1),i)
            CALL gather_grnd(in_biases(1:nhidneurons(1),i),nhidneurons(1),seedval)
            !print *, 'in_biases=', in_biases(1:nhidneurons(1),i)
            DO j=1, nhidlayers-1
                DO k=1, nhidneurons(j)
                    CALL gather_grnd(hid_weights(k,1:nhidneurons(j+1),j,i),nhidneurons(j+1),seedval)
                END DO
                !print *, 'hid_weights=', hid_weights(1:nhidneurons(j),1:nhidneurons(j+1),j,i)
                CALL gather_grnd(hid_biases(1:nhidneurons(j+1),j,i),nhidneurons(j+1),seedval)
                !print *, 'hid_biases=', hid_biases(1:nhidneurons(j+1),j,i)
            END DO
            CALL gather_grnd(out_weights(1:nhidneurons(nhidlayers),1,i),nhidneurons(nhidlayers),seedval)
            !print *, 'out_weights=', out_weights(1:nhidneurons(nhidlayers),1,i)
            out_biases(i)=grnd(seedval)
            !print *, 'out_biases=', out_biases(i)
        END DO
        
    END SUBROUTINE

!    SUBROUTINE write_mlff(nelement,uniq_elements,trained)
!    !------------------------------------------------------------------------!
!    !Write (not)trained model in a file name, trained.pyamff or mlff.pyamff  !
!    !trained.pyamff contains trained model parameters with fp parameters     !
!    !mlff.pyamff contains not completely trained model parameters but saved  !
!    !for restart.                                                            !
!    !------------------------------------------------------------------------!
!        IMPLICIT NONE
!        !Inputs
!        INTEGER :: nelement
!        LOGICAL :: trained
!        CHARACTER*3, DIMENSION(nelement) :: uniq_elements
!        !Variables
!        INTEGER :: i, l
!
!        !TODO: should be more formatted!!
!        IF (trained) THEN 
!            OPEN(unit=15, file='trained.pyamff', status='unknown')
!        ELSE
!            OPEN(unit=15, file='mlff.pyamff', status='unknown')
!        END IF
!
!        !FP parameters
!        !TODO: finish this!!!
!        WRITE(15,'(a)') '#Fingerprint type'
!        WRITE(15,'(a)') 'BP'
!        WRITE(15,'(a)') '#Rmins'
!        WRITE(15,*) uniq_elements
!        DO i=1, nelement
!          WRITE(15,*) (rmins(i,j),j=1,nelement)
!        END DO
!        IF (use_cohesive_energy) THEN 
!          WRITE(15,'(a)') '#Cohesive'
!          WRITE(15,*) coheElement
!          WRITE(15,*) coheEs
!        END IF
!        IF (nG1 .GT. 0) THEN 
!          WRITE(15,'(a)') '#    type   number'
!          WRITE(15,*) 'G1', nG1
!          WRITE(15,'(a)') '#  center neighbor      eta       Rs     rcut       fpmin        fpmax'
!          DO i=1, nelement
!            DO j=1, nelement
!              DO k=1, fpParas(i)%g1s(j)%nFPs  
!                WRITE(15,*) uniq_elements(i), uniq_element(j), fpParas(i)%g1s(j)%etas(k), &
!                fpParas(i)%g1s(j)%rss(k), fpParas(i)%g1s(j)%r_cuts(k) !min, max
!              END DO
!            END DO
!          END DO
!        END IF
!        IF (nG2 .GT. 0) THEN 
!          WRITE(15,'(a)') '#    type   number'
!          WRITE(15,*) 'G2', nG2
!          WRITE(15,'(a)') '#  center neighbor1 neighbor2      eta     zeta    lambda   thetas     rcut        fpmin     fpmax'
!        END IF 
!        WRITE(15,'(a)') '#Model Parameters'
!        DO i=1, nelement
!          WRITE(15,'(A2)') uniq_elements(i)
!          WRITE(15,*) in_weights(1:nGs(i),1:nhidneurons(1),i), '#inputLayer weight  '
!          WRITE(15,*) in_biases(1:nhidneurons(1),i), '#inputLayer bias'
!          DO l=1, nhidlayers-1
!            WRITE(15,*) hid_weights(1:nhidneurons(l),1:nhidneurons(l+1),l,i),&
!                                   '#hiddenLayer_',l,'weight'
!            WRITE(15,*) hid_biases(1:nhidneurons(l+1),l,i),&
!                                   '#hiddenLayer_',l,'bias'
!          END DO
!          WRITE(15,*) out_weights(1:nhidneurons(nhidlayers),1,i), '#outputLayer weight  '
!          WRITE(15,*) out_biases(i), '#outputLayer bias'
!        END DO
!        WRITE(15,'(a)') '#Energy Scaling Parameters'
!        WRITE(15,'(a)') scaler_type
!        WRITE(15,*) slope, intercept
!        CLOSE(15)
!
! 
!
!    END SUBROUTINE

    SUBROUTINE cleanup() BIND(C,name='cleanup')
         USE, INTRINSIC :: iso_c_binding
         DEALLOCATE(fpParas)
         DEALLOCATE(uniq_elements)
         DEALLOCATE(rmins)
    END SUBROUTINE

    SUBROUTINE atomsCleanup()
         DEALLOCATE(pairs)
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
 
END MODULE
