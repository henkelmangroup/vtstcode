MODULE nlist
    USE atomsProp
    IMPLICIT NONE
    PUBLIC

    CONTAINS

!-----------------------------------------------------------------------------------!
! calc:  calculate the neighborlist based upon positions, pos, and cell, cell
!-----------------------------------------------------------------------------------!

    SUBROUTINE calcCellNum(nAtoms, pos_car, cell, rcut, ncells)
        INTEGER :: nAtoms
        DOUBLE PRECISION :: rcut
        DOUBLE PRECISION,DIMENSION(nAtoms,3) :: pos_car
        DOUBLE PRECISION,DIMENSION(3,3) :: cell

!f2py   INTENT(IN) :: nAtoms
!f2py   INTENT(IN) :: pos_car
!f2py   INTENT(IN) :: cell
!f2py   INTENT(IN) :: rcut
        !Variables
        INTEGER :: i, j, k
        INTEGER, DIMENSION(3) :: sort_ids
        DOUBLE PRECISION, DIMENSION(3) :: latt_const, orth_latt, cross_ab,cross_bc
        DOUBLE PRECISION :: tmp, cross_ab_norm, cross_bc_norm, sin_3

        !Output
        INTEGER, DIMENSION(3) :: ncells
!f2py   INTENT(OUT) :: ncells

        latt_const(1) = norm(cell(1, :))
        latt_const(2) = norm(cell(2, :))
        latt_const(3) = norm(cell(3, :))
        !Rank lattice constant
        !print*, 'lattice', latt_const
        sort_ids = (/1, 2, 3/)
        DO i =1, 3
          DO j = 1, 3-i
             IF (latt_const(j) > latt_const(j+1)) THEN
               k = sort_ids(j)
               sort_ids(j) = sort_ids(j+1)
               sort_ids(j+1) = k
               tmp = latt_const(j)
               latt_const(j) = latt_const(j+1)
               latt_const(j+1) = tmp
             END IF
          END DO
        END DO
        orth_latt(sort_ids(1)) = latt_const(1)
        !print*, 'orth_ids1', orth_latt(sort_ids(1)), sort_ids(1)
        !tmp = abs(cross_product(cell(sort_ids(1), :),cell(sort_ids(2), : )))
        cross_ab = cross_product(cell(sort_ids(1), :), cell(sort_ids(2),:))
        cross_ab_norm = norm(cross_ab)
        IF (cross_ab_norm ==0) THEN
           print*, 'Cell is not 3-D'
           CALL EXIT(1)
        END IF
        orth_latt(sort_ids(2)) = cross_ab_norm/latt_const(1)
        cross_bc = cross_product(cell(sort_ids(3), :),cross_ab)
        cross_bc_norm = norm(cross_bc)
        !print*, 'cross bc', cross_bc
        !print*, 'bc norm', cross_bc_norm, cross_ab_norm
        sin_3 = cross_bc_norm/(cross_ab_norm*latt_const(3))
        !print*, 'angle', sin_3
        IF (sin_3 <= 0.05) THEN   !sin(4_degree)=0.07
          orth_latt(sort_ids(3)) = latt_const(3)
        ELSE
          orth_latt(sort_ids(3)) = cross_bc_norm/cross_ab_norm
        END IF
        !print*, 'orth 3', orth_latt(sort_ids(3))
        !print*,'orth_latt',orth_latt
        !print*, 'rcut', rcut
        ncells = INT(2*rcut /orth_latt+0.5) + 1
        !print*, 'ncells', ncells
        !CALL EXIT(1)
    END SUBROUTINE

    SUBROUTINE genSupercell(nAtoms, pos_car, symbols, cell, ncells, tncells, rcut, supercell)
                            !supercell, supersymbols, pool_pos_car, pool_ids)

        INTEGER :: nAtoms, tncells
        DOUBLE PRECISION :: rcut
        DOUBLE PRECISION,DIMENSION(nAtoms,3) :: pos_car
        INTEGER, DIMENSION(nAtoms) :: symbols
        DOUBLE PRECISION,DIMENSION(3,3) :: cell

        !Variables
        INTEGER :: i, j, k, i_offset, id
        INTEGER, DIMENSION(3) :: ncells
        INTEGER,DIMENSION(tncells,3) :: offset_list  !
        !INTEGER,DIMENSION(tncells*nAtoms) :: pool_global_ids
        !INTEGER,DIMENSION(tncells*nAtoms,3) :: pool_offset
        DOUBLE PRECISION,DIMENSION(tncells,3) :: offset_vector

        !Output
        DOUBLE PRECISION,DIMENSION(3,3) :: supercell
        !INTEGER, DIMENSION(tncells*nAtoms) :: supersymbols
        !DOUBLE PRECISION,DIMENSION(tncells*nAtoms,3) :: pool_pos_car
        !INTEGER,DIMENSION(tncells*nAtoms) :: pool_ids

        i_offset = 1
        DO i = 1, ncells(1)
           DO j = 1, ncells(2)
              DO k = 1, ncells(3)
                 offset_list(i_offset, :) = (/i-1, j-1, k-1/)
                 i_offset = i_offset + 1
              END DO
           END DO
        END DO
        offset_vector = MATMUL(offset_list, cell)
        supercell(1,:) = ncells(1) * cell(1,:)
        supercell(2,:) = ncells(2) * cell(2,:)
        supercell(3,:) = ncells(3) * cell(3,:)
        supersymbols = 0
        DO i = 1, tncells
           DO j = 1, nAtoms
              id = (i-1)*nAtoms+j
              pool_pos_car(id,:) = offset_vector(i,:) + pos_car(j,:)
              !write(11, *)pool_pos_car(id, 1), pool_pos_car(id, 2), pool_pos_car(id, 3)
              pool_ids(id) = j
              supersymbols(id) = symbols(j)
              !pool_global_ids(id) = id
              !pool_offset(id,:) = offset_list(i,:)
           END DO
        END DO
        !print*, 'gen', supersymbols(1:10)

    END SUBROUTINE

    SUBROUTINE calcNlist(nAtoms, MAX_NEIGHS, pos_car, cell, symbols, rcut, nelement, &
                         forceEngine, rmins, tncells, &
                         npairs_incell, num_pairs, num_neigh, neighs, neighs_incell)

        DOUBLE PRECISION,PARAMETER :: PI = 4.*ATAN(1.0_8)
        DOUBLE PRECISION,PARAMETER :: RAD2DEG = 180.0_8/PI

        ! input values
        INTEGER :: nAtoms, MAX_NEIGHS, nelement, forceEngine, tncells
        DOUBLE PRECISION :: rcut
        DOUBLE PRECISION,INTENT(IN),DIMENSION(nAtoms,3) :: pos_car
        DOUBLE PRECISION,INTENT(IN),DIMENSION(3,3) :: cell
        INTEGER, DIMENSION(nAtoms) :: symbols

        ! variables
        INTEGER, DIMENSION(nAtoms) :: ghost_indices 
        DOUBLE PRECISION,DIMENSION(nAtoms,3) :: pos_dir
        DOUBLE PRECISION,DIMENSION(3,3) :: car2dir, dir2car
        DOUBLE PRECISION,DIMENSION(3) :: v_dir, v_car, v_unit
        DOUBLE PRECISION,DIMENSION(3) :: v_car_o, v_dir_o
        DOUBLE PRECISION :: dist, dist2, rcut2, dist2_o
        INTEGER :: i, j, check, natoms_percell, n_ghost

        !output valeus
        INTEGER :: num_pairs, npairs_incell
        INTEGER, DIMENSION(nAtoms) :: num_neigh
        INTEGER,DIMENSION(nAtoms,MAX_NEIGHS) :: neighs
        INTEGER,DIMENSION(nAtoms,MAX_NEIGHS) :: neighs_incell
        DOUBLE PRECISION, DIMENSION(nelement, nelement) :: rmins

        rcut2 = rcut**2  ! should check the rcut is defined
        natoms_percell = nAtoms/tncells 
        dir2car = TRANSPOSE(cell)
        car2dir = inverse(dir2car)
        ! to direct coordinates
        DO i = 1, nAtoms
            !print*, i, pos_car(i,:)
            pos_dir(i,:) = MATMUL(car2dir, pos_car(i,:))
        END DO
        ! pair terms
        num_pairs = 0
        npairs_incell = 0
        num_neigh = 0
        neighs = 0
        pair_indices = 0
        pair_global_indices=0
        unitvects_pair=0
        pair_info = 0
        gvects = 0
        num_eachpair = 0

        DO i = 1, nAtoms
            DO j = i+1, nAtoms
                !check = check +1
                !print *, check
                v_dir = pos_dir(j,:) - pos_dir(i,:)
                v_dir_o = v_dir
                v_dir = MOD(v_dir + 100.5, 1.0) - 0.5
                v_car = MATMUL(dir2car, v_dir)
                dist2 = SUM(v_car*v_car)
 
                IF (dist2<rcut2) THEN
                    num_neigh(i) = num_neigh(i) + 1
                    num_neigh(j) = num_neigh(j) + 1
                    dist = SQRT(dist2)
                    !print*, i, j, dist
                    v_unit = v_car/dist

                    !Record minimum distance of two atoms
                    IF (dist<rmins(symbols(pool_ids(i)), symbols(pool_ids(j)))) THEN
                        IF (forceEngine == 0) THEN
                            rmins(symbols(pool_ids(i)), symbols(pool_ids(j))) = dist
                            rmins(symbols(pool_ids(j)), symbols(pool_ids(i))) = dist
                        ELSE IF (forceEngine == 1) THEN
                            dist = rmins(symbols(pool_ids(i)), symbols(pool_ids(j)))
                        END IF
                    END IF

                    !make neighs to record the index in pairs not origin structure
                    IF (i <= natoms_percell) THEN
                       npairs_incell = npairs_incell+1
                    END IF
                    num_pairs = num_pairs + 1
                    num_eachpair(i) = num_eachpair(i) + 1
                    num_eachpair(j) = num_eachpair(j) + 1
                    neighs(i, num_neigh(i)) = j
                    neighs(j, num_neigh(j)) = i
                    neighs_incell(i, num_neigh(i)) = pool_ids(j)
                    neighs_incell(j, num_neigh(j)) = pool_ids(i)
                    pairs(num_pairs) = dist
                    
                    pair_indices(1, num_pairs) = pool_ids(i)
                    pair_indices(2, num_pairs) = pool_ids(j)

                    unitvects_pair(1, num_pairs,:) = v_unit
                    unitvects_pair(2, num_pairs,:) = -v_unit
                    pair_info(i, num_eachpair(i)) = num_pairs

                    pair_global_indices(1, num_pairs) = i
                    pair_global_indices(2, num_pairs) = j
                    IF (ANY((v_dir*v_dir_o<0.0))) THEN
                       gvects(num_pairs) = -1
                    ELSE
                       gvects(num_pairs) = 1
                    END IF
                END IF
            END DO
        END DO

    END SUBROUTINE

!-----------------------------------------------------------------------------------!
! norm:  return the norm of A(3) vector
!-----------------------------------------------------------------------------------!

    FUNCTION norm(A)

        DOUBLE PRECISION,DIMENSION(3) :: A
        DOUBLE PRECISION :: norm
        norm = SQRT(SUM(A * A))
        RETURN

    END FUNCTION

!-----------------------------------------------------------------------------------!
! cross_product: return the cross product of a(3) and b(3)
!-----------------------------------------------------------------------------------!

    FUNCTION cross_product(a,b)

        DOUBLE PRECISION,DIMENSION(3) :: a,b
        DOUBLE PRECISION,DIMENSION(3) :: cross_product

        cross_product(1) = a(2)*b(3) - a(3)*b(2)
        cross_product(2) = a(3)*b(1) - a(1)*b(3)
        cross_product(3) = a(1)*b(2) - a(2)*b(1)

        RETURN
    END FUNCTION cross_product

!-----------------------------------------------------------------------------------!
! inverse:  return the inverse of A(3,3)
!-----------------------------------------------------------------------------------!

    FUNCTION inverse(A)

        DOUBLE PRECISION,INTENT(IN),DIMENSION(3,3) :: A
        DOUBLE PRECISION,DIMENSION(3,3) :: inverse
        DOUBLE PRECISION :: det

        det = determinant(A)
        IF (det == 0) STOP 'Divide by zero in matrix inverse'
        inverse = adjoint(A)/det

        RETURN
    END FUNCTION inverse

!-----------------------------------------------------------------------------------!
! adjoint: return the adjoint of A(3,3)
!-----------------------------------------------------------------------------------!

    FUNCTION adjoint(A)

        DOUBLE PRECISION,INTENT(IN),DIMENSION(3,3) :: A
        DOUBLE PRECISION,DIMENSION(3,3) :: adjoint

        adjoint(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
        adjoint(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
        adjoint(1,3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)

        adjoint(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
        adjoint(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
        adjoint(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)

        adjoint(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
        adjoint(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
        adjoint(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

        RETURN
    END FUNCTION adjoint

!-----------------------------------------------------------------------------------!
! determinant: of a 3x3 matrix 
!-----------------------------------------------------------------------------------!

    FUNCTION determinant(A)

        DOUBLE PRECISION,INTENT(IN),DIMENSION(3,3) :: A
        DOUBLE PRECISION :: determinant

        determinant = A(1,1)*A(2,2)*A(3,3) &
                    - A(1,1)*A(2,3)*A(3,2) &
                    - A(1,2)*A(2,1)*A(3,3) &
                    + A(1,2)*A(2,3)*A(3,1) &
                    + A(1,3)*A(2,1)*A(3,2) &
                    - A(1,3)*A(2,2)*A(3,1)

        RETURN
    END FUNCTION

END MODULE
