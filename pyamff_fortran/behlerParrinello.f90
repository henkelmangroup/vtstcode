MODULE fbp
    USE nlist
    USE atomsProp
    IMPLICIT NONE
    PUBLIC

    CONTAINS

    SUBROUTINE fg1(Rijs, nFPs, etas,r_ss,r_cuts, outputs)

        IMPLICIT NONE
        !Define input values
        INTEGER :: nFPs
        DOUBLE PRECISION, DIMENSION(nFPs) :: etas, r_ss, r_cuts
        DOUBLE PRECISION :: Rijs

        !Define values needed for calculation
        DOUBLE PRECISION, DIMENSION(nFPs) :: outputs
        DOUBLE PRECISION, PARAMETER :: Pi = 4.*ATAN(1.0_8)

        outputs = EXP(-((Rijs-r_ss)**2/(r_cuts**2))*etas)*(0.5*(COS(Pi*Rijs/r_cuts)+1))

        RETURN
    END SUBROUTINE

    SUBROUTINE fdg1(Rij, nFPs, etas,r_ss,r_cuts, outputs)

        IMPLICIT NONE
        !Define input values
        INTEGER :: nFPs
        DOUBLE PRECISION :: Rij
        DOUBLE PRECISION, DIMENSION(nFPs) :: etas, r_ss, r_cuts
        !Define parameters
        DOUBLE PRECISION, PARAMETER :: Pi = 4.*ATAN(1.0_8)
        !Define outputs
        DOUBLE PRECISION, DIMENSION(nFPs) :: outputs

        outputs = (((-etas*(Rij-r_ss)*(COS(Pi*Rij/r_cuts) + 1))/(r_cuts**2)) &
          - (0.5*(Pi/r_cuts)*SIN(Pi*Rij/r_cuts)))&
          * EXP(-etas*((Rij-r_ss)**2)/(r_cuts**2))

    END SUBROUTINE


    SUBROUTINE commonparts(npairs, pairs, r_cut, &
                           pair_sins, fcuts, fexp)
        IMPLICIT NONE
        !Define input vales
        INTEGER :: npairs
        DOUBLE PRECISION, DIMENSION(npairs) :: pairs
        DOUBLE PRECISION :: r_cut

        !Define output vales
        DOUBLE PRECISION, DIMENSION(npairs) :: pair_sins
        DOUBLE PRECISION, DIMENSION(npairs) :: fcuts, fexp

        !Define parameters
        DOUBLE PRECISION, PARAMETER :: Pi = 4.*ATAN(1.0_8)

        pair_sins = SIN(Pi*pairs/r_cut)
        fcuts = 0.5*(COS(Pi*pairs/r_cut)+1)
        !only for g2
        fexp = EXP(-pairs**2 / r_cut **2)

    END SUBROUTINE

    !SUBROUTINE calcTriplet(nAtoms, tnAtoms, MAX_NEIGHS, rcut, npairs, pairs, pair_indices, pair_global_indices,&
    !                       unitvects, ntriplets, angles, angle_indices)
    SUBROUTINE calcTriplet(nAtoms, tnAtoms, MAX_NEIGHS, rcut, npairs, &
                           ntriplets, angles, angle_indices)
        IMPLICIT NONE

        ! Define input values
        INTEGER :: nAtoms, MAX_NEIGHS   !nAtoms: number of atoms in original cell
        INTEGER :: tnAtoms !number of atoms in the whole supercell
        INTEGER :: npairs !number of pairs in the whole supercell
        DOUBLE PRECISION :: rcut, rcut2

        ! Define variables
        INTEGER :: i , i_neigh_1, j_neigh, k_neigh, ij_pair_loc, ik_pair_loc, jk_pair_loc
        INTEGER :: p_size, i_neigh_2, s_neigh, l_neigh, s_pair_loc, l_pair_loc
        INTEGER, DIMENSION(1) :: jk_index !used to fetch the 3rd side
        DOUBLE PRECISION, DIMENSION(3) :: v1_unit, v2_unit, v3_unit

        ! Define output values
        INTEGER :: ntriplets
        DOUBLE PRECISION, DIMENSION(3, nAtoms*MAX_NEIGHS*MAX_NEIGHS) :: angles 
        INTEGER, DIMENSION(3, nAtoms*MAX_NEIGHS*MAX_NEIGHS) :: angle_indices
        ! Place holders
        DOUBLE PRECISION :: a1, b1, c1, dist_3

        rcut2 = rcut**2
        ntriplets = 0
        angles = 0
        angle_indices = 0
        !print*, 'Triplet nAtoms', nAtoms
        DO i = 1, nAtoms !for each atom 
            !print*, 'center', i, 'neighs',pair_info(i,1:num_eachpair(i))
            DO j_neigh = 1, num_eachpair(i) !for each pair for that atom
                ij_pair_loc = pair_info(i,j_neigh)  
                v1_unit(:)=unitvects_pair(1,ij_pair_loc,:)
                DO k_neigh = j_neigh+1, num_eachpair(i) !for next neighbor atom
                    !only calculate triplet that has all distance in the cutoff range
                    !select the pairs that are led by neighbors of atom i
                    !pair_info(pair_global_indices(ij_pair),:) give index of pairs
                    ik_pair_loc = pair_info(i,k_neigh)
                    s_neigh = pair_global_indices(2, ij_pair_loc)
                    l_neigh = pair_global_indices(2, ik_pair_loc)

                    IF (l_neigh < s_neigh) THEN ! make sure j is always the smaller side
                        CYCLE
                    END IF

                    IF (s_neigh > tnAtoms) THEN
                       CYCLE
                    END IF
                    ! when pair_indices(j) < nAtoms but all pairs belong to this atom have been included in atoms with smaller index,
                    ! then j_ps =0; for pair_indices(j)= nAtoms case hase been exclued in j loop, because ps = 0
                    p_size = SIZE(pair_global_indices(2, pair_info(s_neigh,1:num_eachpair(s_neigh))))
                    IF (p_size .EQ. 0) THEN
                        CYCLE
                    END IF
                    jk_index = find_loc_int(pair_global_indices(2, pair_info(s_neigh,1:num_eachpair(s_neigh))), &
                                            l_neigh, &
                                            p_size)
                    IF (jk_index(1) /= 0) THEN !i.e find_loc_int found an index, i.e atom in pairs array
                        jk_pair_loc = pair_info(s_neigh,jk_index(1))
                        v2_unit(:) = unitvects_pair(1, ik_pair_loc, :)  !vect_ik
                        IF ((gvects(ik_pair_loc) < 0).OR.(gvects(ij_pair_loc)<0)) THEN !!!IF ghost vector is a vector in triangle 
                           a1 = SUM(v1_unit * v2_unit) !!Dot Product, i.e Theta
                           dist_3 = pairs(ij_pair_loc) * pairs(ij_pair_loc) +  &
                                    pairs(ik_pair_loc) * pairs(ik_pair_loc) -  &
                                    2 * pairs(ij_pair_loc) * pairs(ik_pair_loc) * a1

                           v3_unit(:) = (pairs(ik_pair_loc) * v2_unit(:)  - pairs(ij_pair_loc)*v1_unit(:))/sqrt(dist_3) 
                        ELSE !!Otherwise there is not ghost vector in the Triangle; BUT WE NEED TO COMPUTE DIST3 STILL,WE CAN HAVE DIST3>Rc but NO GHOST VECTOR
                          !cosine is even, can do this for sure
                           a1 = SUM(v1_unit*v2_unit)
                           v3_unit(:) = unitvects_pair(1, jk_pair_loc, :) !vect_jk
                           dist_3 = pairs(ij_pair_loc) * pairs(ij_pair_loc) +  &
                                    pairs(ik_pair_loc) * pairs(ik_pair_loc) -  &
                                    2 * pairs(ij_pair_loc) * pairs(ik_pair_loc) * a1
                          !dist_3 = 0.0
                        END IF
                        IF (dist_3 < rcut2) THEN
                           ntriplets = ntriplets + 1
                           IF (ABS(a1) > 0.999999999999) THEN
                               angles(1, ntriplets) = 0.d0
                           ELSE 
                               angles(1, ntriplets) = ACOS(a1)
                           END IF
                           b1 = SUM(-v1_unit*v3_unit)
                           IF (ABS(b1) > 0.999999999999) THEN
                               angles(2, ntriplets) = 0.d0
                           ELSE 
                               angles(2, ntriplets) = ACOS(b1)
                           END IF
                           c1 = SUM(v2_unit*v3_unit)
                           IF (ABS(c1) > 0.999999999999) THEN
                               angles(3, ntriplets) = 0.d0
                           ELSE 
                               angles(3, ntriplets) = ACOS(c1)
                           END IF
                           
                           angle_indices(1, ntriplets) = ij_pair_loc  !j
                           angle_indices(2, ntriplets) = ik_pair_loc  !k
                           angle_indices(3, ntriplets) = jk_pair_loc  !jk_index

                        END IF
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE

    SUBROUTINE fg2(nFPs, theta_ss, etas, zetas, lambdas, theta, fexpt, fct, &
                   outputs, f_thetas_1, f_thetas, f_rs)
        IMPLICIT NONE
        !Define input values
        INTEGER :: nFPs
        DOUBLE PRECISION, DIMENSION(nFPs) :: etas, theta_ss, zetas, lambdas
        DOUBLE PRECISION :: theta
        !DOUBLE PRECISION, DIMENSION(3) :: Rs   !(Rij, Rik, Rjk)
        DOUBLE PRECISION :: fexpt !fexp_ij, fexp_ik, fexp_jk
        DOUBLE PRECISION :: fct   !fc_rij, fc_rik, fc_rjk

        !Define values needed for calculation
        DOUBLE PRECISION, PARAMETER :: Pi = 4.*ATAN(1.0_8)
        DOUBLE PRECISION, PARAMETER :: radians = Pi/180.0_8 !same as np.radians

        !Define output arrays
        DOUBLE PRECISION, DIMENSION(nFPs) :: outputs
        DOUBLE PRECISION, DIMENSION(nFPs) :: f_thetas, f_rs, f_thetas_1

        !print *, 'calcuting g2'
        !fexpt = fexp(1) * fexp(2) * fexp(3) 
        !DO i = 1, 3
        !f_thetas(i) = (1 + lambdas(i,:) * cos(thetas(i) - theta_ss(i,:))) ** zetas(i, :)
        f_thetas_1 = 1 + lambdas * COS(theta - theta_ss)
        !print *, "In fg2 COS(theta - theta_ss):  +278",COS(theta - theta_ss)
        f_thetas = f_thetas_1 ** (zetas-1)
        !print *, f_thetas
        !f_rs(i) = fexpt * etas(i,:)
        f_rs = fexpt ** etas
        !outputs(i) = f_thetas(i) * f_rs(i) * fexpt * 2 ** (1-zetas)
        outputs = f_thetas * f_thetas_1 * f_rs * fct * 2 ** (1-zetas)
        !ENDDO

    END SUBROUTINE

    SUBROUTINE fdg2(nFPs, theta_ss, etas, zetas, lambdas, r_cut, thetas, f_thetas_1, f_thetas, f_rs, Rs, fcs, sins, & 
                    unitvects, vectsigns, &
                    dfp_i, dfp_j, dfp_k)

        IMPLICIT NONE
        ! Define input values
        INTEGER :: nFPs
        DOUBLE PRECISION, DIMENSION(nFPs) :: etas, theta_ss, zetas, lambdas
        DOUBLE PRECISION :: thetas, r_cut
        DOUBLE PRECISION, DIMENSION(nFPs) :: f_thetas_1, f_thetas, f_rs
        DOUBLE PRECISION, DIMENSION(3) :: Rs  !Rij, Rik, Rjk
        DOUBLE PRECISION, DIMENSION(3) :: sins  !sin_ij, sin_ik, sin_jk
        DOUBLE PRECISION, DIMENSION(3) :: fcs !fc_rij, fc_rik, fc_rjk
        !DOUBLE PRECISION :: unitvect_ij, unitvect_ik, unitvect_jk
        DOUBLE PRECISION, DIMENSION(3, 3) :: unitvects
        DOUBLE PRECISION, DIMENSION(3) :: vectsigns

        ! Define variables
        INTEGER :: i
        DOUBLE PRECISION :: pi_rc, inv_Rc2
        DOUBLE PRECISION :: infinity
        DOUBLE PRECISION, DIMENSION(nFPs) :: cv, cos_thetas
        DOUBLE PRECISION, PARAMETER :: Pi = 4.*ATAN(1.0_8)
        DOUBLE PRECISION, PARAMETER :: radians = Pi/180.0_8 !same as np.radians
        DOUBLE PRECISION, DIMENSION(nFPs) :: coeff1, coeff2, coeff3
        DOUBLE PRECISION, DIMENSION(nFPs) :: tcoeff2, tcoeff3
        DOUBLE PRECISION, DIMENSION(nFPs) :: ttcoeff3
        DOUBLE PRECISION, DIMENSION(nFPs) :: i_term1, i_term2, i_term3
        DOUBLE PRECISION, DIMENSION(nFPs) :: j_term1, j_term2, j_term3
        DOUBLE PRECISION, DIMENSION(nFPs) :: k_term1, k_term2, k_term3
        DOUBLE PRECISION, DIMENSION(nFPs) :: dfpi_p1, dfpj_p1, dfpk_p1
        DOUBLE PRECISION, DIMENSION(nFPs) :: dfpi_p2, dfpj_p2, dfpk_p2
        DOUBLE PRECISION, DIMENSION(nFPs) :: dfpj_p3, dfpk_p3
        !COMMON pi_rc, inv_Rc2
        
        ! Define output arrays
        !DOUBLE PRECISION, DIMENSION(nFPs) :: arm1, arm2, arm1_to_neigh, arm2_to_neigh, arm3_to_neigh
        DOUBLE PRECISION, DIMENSION(3, nFPs) :: dfp_i, dfp_j, dfp_k

        !print *, 'debug'
        !print *, thetas
        !print *, f_thetas, f_rs
        pi_rc = Pi/(2.*r_cut)
        inv_Rc2 = 2./(r_cut**2)

        cv = f_thetas*f_rs*(2**(1.-zetas))
        cos_thetas = COS(thetas)
        !IF (thetas == Pi) THEN
        !   cos_thetas = COS(thetas - 0.0001)
        !ELSE IF (thetas == 0.0) THEN
        !   cos_thetas = COS(thetas + 0.0001)
        !ELSE
        !   cos_thetas = COS(thetas)
        !END IF

        coeff1 = pi_rc
        coeff2 = inv_Rc2*etas

        !Following AMP, if angle is 180 degrees, set coeff3 to zero
        !Else if denominator is eqaul to 0 and thetas = theta_ss, set coeff3 to deriv wrt thetas (lambdas*zetas/(lambdas+1))
        !Otherwise, throw errors and stop the program
        
        IF ((thetas /= Pi).AND. (thetas /= 0)) THEN
           coeff3 = lambdas * zetas * SIN(thetas - theta_ss) / &
                   SQRT(1.-cos_thetas**2) 
           ! w.r.t. Atom i
           !   ij arm
           i_term1 = pi_rc * sins(1) * fcs(2)
         
           tcoeff2 = coeff2 * fcs(1) * fcs(2)
           i_term2 = tcoeff2 * Rs(1)
         
           tcoeff3 = coeff3 * fcs(1) * fcs(2)
           i_term3 = tcoeff3 * (cos_thetas/Rs(1) - 1/Rs(2))
         
           dfpi_p1 = f_thetas_1 * (i_term1 + i_term2) + i_term3
         
           !   ik arm
           i_term1 = pi_rc * sins(2) * fcs(1)
         
           i_term2 = tcoeff2 * Rs(2)
         
           i_term3 = tcoeff3 * (cos_thetas/Rs(2) - 1/Rs(1))
         
           dfpi_p2 = f_thetas_1 * (i_term1 + i_term2) + i_term3
         
           ! w.r.t. Atom j and k
           !    ij arm for j; ik arm for k
           j_term1 = pi_rc * sins(1)
           k_term1 = pi_rc * sins(2)
         
           j_term2 = coeff2 * fcs(1) * Rs(1)
           k_term2 = coeff2 * fcs(2) * Rs(2)
         
           !print *, 'coeff3', coeff3
           ttcoeff3 = coeff3 * cos_thetas
           j_term3 = ttcoeff3 * fcs(1) / Rs(1)
           k_term3 = ttcoeff3 * fcs(2) / Rs(2)
         
           !print *, 'k terms'
           !print *, k_term1, k_term2, k_term3
           dfpj_p1 = -(j_term3 + f_thetas_1 * (j_term1 + j_term2)) *  fcs(3)
           dfpk_p2 = -(k_term3 + f_thetas_1 * (k_term1 + k_term2)) * fcs(3)
         
           !print *, 'dfpj_p1', dfpj_p1
           !   ik arm for j; ij arm for k
           tcoeff3 = coeff3 * fcs(3)
           dfpj_p2 = tcoeff3 * fcs(1) / Rs(1)
           dfpk_p1 = tcoeff3 * fcs(2) / Rs(2)
         
           !   jk arm
           coeff1 = pi_rc * sins(3)
           j_term1 = coeff1 * fcs(1)
           k_term1 = coeff1 * fcs(2)
         
           tcoeff2 = coeff2 * Rs(3) * fcs(3)
           j_term2 = tcoeff2 * fcs(1)
           k_term2 = tcoeff2 * fcs(2)
         
           dfpj_p3 = f_thetas_1 * (j_term1 + j_term2)
           dfpk_p3 = - f_thetas_1 * (k_term1 + k_term2)
         
           DO i = 1, 3
               dfp_i(i,:) = (dfpi_p1 * unitvects(1,i) * vectsigns(1) + dfpi_p2 * unitvects(2,i) * vectsigns(2)) * cv * fcs(3)
               dfp_j(i,:) = (dfpj_p1 * unitvects(1,i) * vectsigns(1) + dfpj_p2 * unitvects(2,i) * vectsigns(2) + &
                             dfpj_p3 * unitvects(3,i) * vectsigns(3)) * cv * fcs(2)
               dfp_k(i,:) = (dfpk_p1 * unitvects(1,i) * vectsigns(1) + dfpk_p2 * unitvects(2,i) * vectsigns(2) + &
                             dfpk_p3 * unitvects(3,i) *vectsigns(3)) * cv * fcs(1) 
           END DO

        ELSE 
           !IF (thetas == Pi) THEN
           !   print *, coeff3
           !END IF
           
           ! w.r.t. Atom i
           !   ij arm
           cv = f_thetas_1 * cv
           i_term1 = pi_rc * sins(1) * fcs(2)
         
           tcoeff2 = coeff2 * fcs(1) * fcs(2)
           i_term2 = tcoeff2 * Rs(1)
         
           !tcoeff3 = coeff3 * fcs(1) * fcs(2)
           !i_term3 = tcoeff3 * (cos_thetas/Rs(1) - 1/Rs(2))
         
           dfpi_p1 = i_term1 + i_term2
         
           !   ik arm
           i_term1 = pi_rc * sins(2) * fcs(1)
         
           i_term2 = tcoeff2 * Rs(2)
         
           !i_term3 = tcoeff3 * (cos_thetas/Rs(2) - 1/Rs(1))
         
           dfpi_p2 = i_term1 + i_term2
         
           ! w.r.t. Atom j and k
           !    ij arm for j; ik arm for k
           j_term1 = pi_rc * sins(1)
           k_term1 = pi_rc * sins(2)
         
           j_term2 = coeff2 * fcs(1) * Rs(1)
           k_term2 = coeff2 * fcs(2) * Rs(2)
         
           !print *, 'coeff3', coeff3
           !ttcoeff3 = coeff3 * cos_thetas
           !j_term3 = ttcoeff3 * fcs(1) / Rs(1)
           !k_term3 = ttcoeff3 * fcs(2) / Rs(2)
         
           !print *, 'k terms'
           !print *, k_term1, k_term2, k_term3
           dfpj_p1 = -(j_term1 + j_term2) *  fcs(3)
           dfpk_p2 = -(k_term1 + k_term2) * fcs(3)
         
           !print *, 'dfpj_p1', dfpj_p1
           !   ik arm for j; ij arm for k
           !tcoeff3 = coeff3 * fcs(3)
           !dfpj_p2 = tcoeff3 * fcs(1) / Rs(1)
           !dfpk_p1 = tcoeff3 * fcs(2) / Rs(2)
         
           !   jk arm
           coeff1 = pi_rc * sins(3)
           j_term1 = coeff1 * fcs(1)
           k_term1 = coeff1 * fcs(2)
         
           tcoeff2 = coeff2 * Rs(3) * fcs(3)
           j_term2 = tcoeff2 * fcs(1)
           k_term2 = tcoeff2 * fcs(2)
         
           dfpj_p3 = j_term1 + j_term2
           dfpk_p3 = - k_term1 - k_term2
         
           DO i = 1, 3
               dfp_i(i,:) = (dfpi_p1 * unitvects(1,i) * vectsigns(1) + dfpi_p2 * unitvects(2,i) * vectsigns(2)) * cv * fcs(3)
               dfp_j(i,:) = (dfpj_p1 * unitvects(1,i) * vectsigns(1) + &
                             dfpj_p3 * unitvects(3,i) * vectsigns(3)) * cv * fcs(2)
               dfp_k(i,:) = (dfpk_p2 * unitvects(2,i) * vectsigns(2) + &
                             dfpk_p3 * unitvects(3,i) *vectsigns(3)) * cv * fcs(1) 
           END DO
        END IF
        !print *, dfp_i(1, :)
        !print *, 'end debug'
    END SUBROUTINE

!------------------------------------------------------------------------------------!
! find loc function
!------------------------------------------------------------------------------------!

    FUNCTION find_loc_int(arr, val, size_arr)

        INTEGER,INTENT(IN) :: size_arr   
        INTEGER,INTENT(IN),DIMENSION(size_arr) :: arr
        INTEGER,INTENT(IN) :: val
        INTEGER :: i
        INTEGER :: find_loc_int

        find_loc_int = 0
        DO i = 1, size_arr
            IF (arr(i) .EQ. val) THEN
                find_loc_int = i
                EXIT
            END IF
        END DO

        RETURN
    END FUNCTION

    FUNCTION find_loc_char(arr, val, size_arr)

        INTEGER,INTENT(IN) :: size_arr
        CHARACTER*2, INTENT(IN),DIMENSION(size_arr) :: arr
        CHARACTER*2,INTENT(IN) :: val
        INTEGER :: i
        INTEGER :: find_loc_char

        find_loc_char = 0
        DO i = 1, size_arr
            IF (arr(i) .EQ. val) THEN
                find_loc_char = i
                EXIT
            END IF
        END DO

        RETURN
    END FUNCTION

    FUNCTION find_loc_dbl(arr, val, size_arr)

        INTEGER,INTENT(IN) :: size_arr
        DOUBLE PRECISION,INTENT(IN),DIMENSION(size_arr) :: arr
        DOUBLE PRECISION,INTENT(IN) :: val
        INTEGER :: i
        INTEGER :: find_loc_dbl
        find_loc_dbl = 0
        DO i = 1, size_arr
            IF (arr(i) .EQ. val) THEN
                find_loc_dbl = i
                EXIT
            END IF
        END DO

        RETURN
    END FUNCTION


END MODULE
