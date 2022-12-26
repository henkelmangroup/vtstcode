MODULE normalize
    USE fbp, only: find_loc_char, find_loc_int
    USE nnType, only: natoms_arr, nGs, atom_idx, fpminvs, fpmaxvs, diffs, magnitude, interceptScale,&
                      elem_mag, elem_intercept
    USE trainType, only: img_idx, nimages
    USE atomsProp    

    IMPLICIT NONE
    PUBLIC

    CONTAINS

    SUBROUTINE normalizeFPs(nelement, nAtoms, uniq_elements, MAX_FPS, max_nneighs, &
                            symbols, nneighbors,neighborlists, fps, dfps)
        IMPLICIT NONE
        ! Inputs
        INTEGER :: MAX_FPS, max_nneighs, nelement, nAtoms
        INTEGER, DIMENSION(nAtoms) :: nneighbors 
        INTEGER, DIMENSION(nAtoms) :: symbols
        INTEGER, DIMENSION(nelement) :: idx_arr, idx_arr2
        INTEGER, DIMENSION(nAtoms, max_nneighs) :: neighborlists
        CHARACTER*3, DIMENSION(nelement) :: uniq_elements
        DOUBLE PRECISION, DIMENSION(nAtoms, MAX_FPS) :: fps
        DOUBLE PRECISION, DIMENSION(nAtoms, max_nneighs+1, 3, MAX_FPS) :: dfps

        ! Variables
        INTEGER :: ptr, ptr2, i, j, j2, m, max_natoms
        
        idx_arr = 1
        DO i = 1, nAtoms
            ptr = symbols(i)
            !print *, 'fps before normalized of atom ', i
            !print *, fps(i,1:nGs(ptr))
            fps(i,1:nGs(ptr)) = &
            fps(i,1:nGs(ptr))*magnitude(1:nGs(ptr),idx_arr(ptr),ptr) + interceptScale(1:nGs(ptr),idx_arr(ptr),ptr)
            idx_arr2 = 1
            DO j = 1, nneighbors(i)+1
                IF (j == 1) THEN
                    !print *, 'dfps from itself before normalized'
                    DO j2 = 1, 3
                        !print *, 'j2=', j2
                        !print *, dfps(i,j,j2,1:nGs(ptr))
                        dfps(i,j,j2,1:nGs(ptr)) = dfps(i,j,j2,1:nGs(ptr))*magnitude(1:nGs(ptr),idx_arr(ptr),ptr)
                    END DO
                ELSE  
                    m = neighborlists(i,j-1)
                    ptr2 = symbols(m)
                    !print *, 'dfps from neighbor', m, 'before normalized'
                    DO j2 = 1, 3
                        !print *, 'j2=', j2
                        !print *, dfps(i,j,j2,1:nGs(ptr2))
                        dfps(i,j,j2,1:nGs(ptr2)) = dfps(i,j,j2,1:nGs(ptr2))*magnitude(1:nGs(ptr2),idx_arr2(ptr2),ptr2)
                    END DO
                    idx_arr2(ptr2) = idx_arr2(ptr2) + 1
                END IF
            END DO
            idx_arr(ptr) = idx_arr(ptr) + 1
        END DO

    END SUBROUTINE

    SUBROUTINE loadnormalizeParas(nAtoms, nelement, MAX_FPS, symbols, uniq_elements)
        IMPLICIT NONE

        ! Input
        INTEGER :: nAtoms, nelement, MAX_FPS
        INTEGER, DIMENSION(nAtoms) :: symbols
        CHARACTER*3, DIMENSION(nelement) :: uniq_elements

        ! Variables
        INTEGER :: nG1, nG2, numGs, cidx, i, k
        INTEGER, DIMENSION(nelement) :: g1_endpts, g2_endpts
        DOUBLE PRECISION :: djunk, fpmin, fpmax
        CHARACTER*3 :: center, cjunk
        CHARACTER(LEN=30) :: line, G_type

        ! Output
        !DOUBLE PRECISION, DIMENSION(MAX_FPS, nelement) :: fpminvs, fpmaxvs, diffs

        ALLOCATE (natoms_arr(nelement))
        ALLOCATE (nGs(nelement))
        ALLOCATE (fpminvs(MAX_FPS, nelement))
        ALLOCATE (fpmaxvs(MAX_FPS, nelement))
        ALLOCATE (diffs(MAX_FPS, nelement))
        !print *, 'allocated done'
        g1_endpts = 0
        g2_endpts = 0
        nGs = 0
        fpminvs = 0
        fpmaxvs = 0
        diffs = 0
        DO k = 1, nelement
            natoms_arr(k) = COUNT(symbols .EQ. k)
        END DO

        OPEN (11, FILE='mlff.pyamff', status='old')
        READ (11,*) !skip #Fingerprint type
        READ (11,*) !fp_type 
        DO WHILE (line .NE. "#MachineLearning")
            READ (11,*) line !skip #type number
            IF (line .EQ. "#MachineLearning") GOTO 40
            READ (11,*) G_type, numGs

            IF (G_type .EQ. 'G1') THEN
                nG1 = numGs
                READ (11,*) !skip # center neighbor ...
                DO i = 1, nG1
                    READ (11,*) center, cjunk, djunk, djunk, djunk, fpmin, fpmax
                    cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                    nGs(cidx) = nGs(cidx) + 1
                    IF (cidx .EQ. 1) THEN
                        fpminvs(i,cidx) = fpmin
                        fpmaxvs(i,cidx) = fpmax
                        g1_endpts(cidx) = g1_endpts(cidx)+1
                    ELSE
                        fpminvs(i-g1_endpts(cidx-1), cidx) = fpmin
                        fpmaxvs(i-g1_endpts(cidx-1), cidx) = fpmax
                        g1_endpts(cidx) = g1_endpts(cidx) + 1
                    END IF
                END DO

            ELSE IF (G_type .EQ. 'G2') THEN
                nG2 = numGs
                READ (11,*) !skip # center neighbor ...
                DO i = 1, nG2
                    READ (11,*) center, cjunk, cjunk, djunk, djunk, djunk, djunk, djunk, fpmin, fpmax
                    cidx = find_loc_char(uniq_elements, center, SIZE(uniq_elements))
                    nGs(cidx) = nGs(cidx) + 1
                    IF (cidx .EQ. 1) THEN
                        fpminvs(g1_endpts(cidx)+i,cidx) = fpmin
                        fpmaxvs(g1_endpts(cidx)+i,cidx) = fpmax
                        g2_endpts(cidx) = g2_endpts(cidx) + 1
                    ELSE
                        fpminvs(g1_endpts(cidx)+i-g2_endpts(cidx-1),cidx) = fpmin
                        fpmaxvs(g1_endpts(cidx)+i-g2_endpts(cidx-1),cidx) = fpmax
                        g2_endpts(cidx) = g2_endpts(cidx) + 1
                    END IF
                END DO
            END IF
        END DO
 40     BACKSPACE 11

    END SUBROUTINE

    SUBROUTINE normalizeParas(nelement)
        IMPLICIT NONE
        ! Inputs
        INTEGER :: nelement
        !DOUBLE PRECISION, DIMENSION(MAX_FPS, nelement) :: fpminvs, fpmaxvs, diffs

        ! Variables
        INTEGER :: i, j, k
        
        DO i = 1, nelement
            DO j = 1, nGs(i)
                diffs(j,i) = fpmaxvs(j,i) - fpminvs(j,i)
                IF (diffs(j,i) .LT. 1.0E-8) THEN
                    fpminvs(j,i) = -1.0d0
                    diffs(j,i) = 2.d0
                END IF    
                IF (fpminvs(j,i) == -1.0d0) THEN
                    diffs(j,i) = 2.d0
                END IF   
            END DO
            DO k = 1, natoms_arr(i)
                magnitude(1:nGs(i),k,i) = 2.d0/diffs(1:nGs(i),i)
                interceptScale(1:nGs(i),k,i) = -2.d0*fpminvs(1:nGs(i),i)/diffs(1:nGs(i),i) - 1.0d0
            END DO
        END DO

    END SUBROUTINE

    SUBROUTINE normalizeParas2(nelement,nAtoms,MAX_FPS,in_natomsarr,max_natoms_arr,symbols,in_atomidx,fps,ordered_fps)
        !We might need this subroutine at some point but not used currently
        IMPLICIT NONE
        ! Inputs
        INTEGER :: nelement, nAtoms, MAX_FPS, max_natoms_arr
        DOUBLE PRECISION, DIMENSION(nAtoms, MAX_FPS) :: fps
        INTEGER, DIMENSION(nAtoms) :: symbols
        INTEGER, DIMENSION(nelement) :: in_natomsarr
        INTEGER, DIMENSION(max_natoms_arr,nelement) :: in_atomidx
        ! Variables
        INTEGER :: i, j, k, ptr
        INTEGER, DIMENSION(nelement) :: idx_arr
        !DOUBLE PRECISION, DIMENSION(MAX_FPS,nelement) :: minv, maxv
        !outputs
        DOUBLE PRECISION, DIMENSION(max_natoms_arr,MAX_FPS,nelement) :: ordered_fps

        !print *, 'normalize paras is called'
        idx_arr = 1
        DO i = 1, nAtoms
          !print *, 'i=',i
          ptr = symbols(i)
          idx_arr(ptr) = idx_arr(ptr) + 1
        END DO
      
        DO i = 1, nelement
            !print *, 'in_atomidx=', in_atomidx(:,i)
            DO j = 1, in_natomsarr(i)
                ordered_fps(j,1:nGs(i),i) = fps(in_atomidx(j,i),1:nGs(i))
                !print *, 'ordered_fps of atom', j
                !print *, ordered_fps(j,1:nGs(i),i)
            END DO
        END DO

        !Initialize only once
        IF (img_idx == 1) THEN
          fpminvs=0.
          fpmaxvs=0.
          !print *, 'min max initialized'
        END IF
        DO i = 1, nelement
            DO j = 1, nGs(i)    
                IF (fpminvs(j,i) == 0.) THEN
                    fpminvs(j,i)=MINVAL(ordered_fps(1:in_natomsarr(i),j,i))
                ELSE IF (fpminvs(j,i) > MINVAL(ordered_fps(1:in_natomsarr(i),j,i))) THEN
                    fpminvs(j,i)=MINVAL(ordered_fps(1:in_natomsarr(i),j,i))
                ELSE
                    !Do nothing
                    CONTINUE
                END IF
                IF (fpmaxvs(j,i) == 0.) THEN
                    fpmaxvs(j,i)=MAXVAL(ordered_fps(1:in_natomsarr(i),j,i))
                ELSE IF (fpmaxvs(j,i) < MAXVAL(ordered_fps(1:in_natomsarr(i),j,i))) THEN
                    fpmaxvs(j,i)=MAXVAL(ordered_fps(1:in_natomsarr(i),j,i))
                ELSE
                    !Do nothing
                    CONTINUE
                END IF
            END DO
            !print *, 'fpminvs=', fpminvs(1:nGs(i),i)
            !print *, 'fpmaxvs=', fpmaxvs(1:nGs(i),i)
        END DO
        !print *, 'min max are calculated'

        DO i = 1, nelement
            DO j = 1, nGs(i)
                diffs(j,i)=fpmaxvs(j,i) - fpminvs(j,i)
                IF (img_idx == nimages) THEN
                    IF (diffs(j,i) .LT. 1.0E-8) THEN
                        fpminvs(j,i) = -1.0d0
                        diffs(j,i) = 2.d0
                    END IF
                END IF
            END DO
            !print *, 'diffs=', diffs(1:nGs(i),i)
            elem_mag(1:nGs(i),i) = 2.d0/diffs(1:nGs(i),i)
            elem_intercept(1:nGs(i),i) = -2.d0*fpminvs(1:nGs(i),i)/diffs(1:nGs(i),i) - 1.0d0
        END DO
    END SUBROUTINE

    SUBROUTINE normalizeFPs2(nelement, nAtoms, uniq_elements, maxfps, max_nneighs, &
               max_natoms_arr,nneighbors,neighborlists, symbols, ordered_fps, dfps)
        !This module is for training module. Instead of reading fprange from .pyamff file, 
        !this module uses fprange calculated from actual fingerprints
        IMPLICIT NONE
        ! Inputs
        INTEGER :: maxfps, max_nneighs, nelement, nAtoms, max_natoms_arr
        INTEGER, DIMENSION(nAtoms) :: nneighbors 
        INTEGER, DIMENSION(nAtoms) :: symbols
        INTEGER, DIMENSION(nAtoms, max_nneighs) :: neighborlists
        CHARACTER*2, DIMENSION(nelement) :: uniq_elements
        DOUBLE PRECISION, DIMENSION(max_natoms_arr,maxfps,nelement) :: ordered_fps
        DOUBLE PRECISION, DIMENSION(nAtoms, max_nneighs+1, 3, maxfps) :: dfps

        ! Variables
        INTEGER :: ptr, ptr2, i, j, j2, m
        INTEGER, DIMENSION(nelement) :: idx_arr, idx_arr2
        
        idx_arr = 1
        DO i = 1, nAtoms
            ptr = symbols(i)
            !print *, 'before normalized fps=', ordered_fps(idx_arr(ptr),1:nGs(ptr),ptr)
            ordered_fps(idx_arr(ptr),1:nGs(ptr),ptr) = &
            ordered_fps(idx_arr(ptr),1:nGs(ptr),ptr)*magnitude(1:nGs(ptr),idx_arr(ptr),ptr) &
            + interceptScale(1:nGs(ptr),idx_arr(ptr),ptr)
            !print *, 'ordered_fps=', ordered_fps(idx_arr(ptr),1:nGs(ptr),ptr)
            idx_arr2 = 1
            DO j = 1, nneighbors(i)+1
                IF (j == 1) THEN
                    DO j2 = 1, 3
                        dfps(i,j,j2,:) = dfps(i,j,j2,1:nGs(ptr))*magnitude(1:nGs(ptr),idx_arr(ptr),ptr)
                    END DO
                ELSE  
                    m = neighborlists(i,j-1)
                    ptr2 = symbols(m)
                    DO j2 = 1, 3
                        dfps(i,j,j2,1:nGs(ptr2)) = dfps(i,j,j2,1:nGs(ptr2))*magnitude(1:nGs(ptr2),idx_arr2(ptr2),ptr2)
                    END DO
                    idx_arr2(ptr2) = idx_arr2(ptr2) + 1
                END IF
            END DO
            idx_arr(ptr) = idx_arr(ptr) + 1
        END DO
    END SUBROUTINE

    SUBROUTINE ghost_dfps_correct(nelement, nAtoms, MAX_FPS, max_nneighs, nneighbors, &
    neighborlists, sub_nneighbors, sub_nlist, in_dfps, out_dfps)
        IMPLICIT NONE
        INTEGER :: MAX_FPS, max_nneighs, nelement, nAtoms, len_subnlist
        INTEGER, DIMENSION(nAtoms) :: nneighbors, sub_nneighbors
        INTEGER, DIMENSION(nAtoms, max_nneighs) :: neighborlists
        INTEGER, DIMENSION(nAtoms, nAtoms) :: sub_nlist
        DOUBLE PRECISION, DIMENSION(nAtoms, max_nneighs, 3, MAX_FPS) :: in_dfps
        DOUBLE PRECISION, DIMENSION(nAtoms, nAtoms, 3, MAX_FPS) :: out_dfps
        ! Variables
        INTEGER :: ptr, ptr2, i, j, j2, m, max_natoms, m_loc, idx

        sub_nlist=0
        sub_nneighbors=0
        DO i = 1, nAtoms
          !print *, 'neighborlists of ',i,'=', neighborlists(i,1:nneighbors(i))
          idx=1
          DO j = 1, nneighbors(i)
            IF (idx <= nAtoms) THEN
              m=neighborlists(i,j)
              IF ((ANY(sub_nlist(i,:)==m)) .OR. (m==i)) THEN 
                CONTINUE
              ELSE
                sub_nlist(i,idx)=m
                sub_nneighbors(i)=sub_nneighbors(i)+1
                idx=idx+1
              END IF
            END IF
          END DO
        END DO

        out_dfps=0

        DO i = 1, nAtoms
            ptr = supersymbols(i)
            DO j2 = 1, 3
                DO j = 1, nneighbors(i)+1
                    m = neighborlists(i,j)
                    IF (j == 1) THEN
                        out_dfps(i,j,j2,1:nGs(ptr))=out_dfps(i,j,j2,1:nGs(ptr))+in_dfps(i,j,j2,1:nGs(ptr))
                    ELSE
                        m = neighborlists(i,j-1)
                        IF (m <= nAtoms) THEN
                            IF (m==i) THEN
                                m_loc=1
                            ELSE
                                m_loc = find_loc_int(sub_nlist(i,:),m,nAtoms)+1
                            END IF
                            ptr2 = supersymbols(m)
                            out_dfps(i,m_loc,j2,1:nGs(ptr2))=out_dfps(i,m_loc,j2,1:nGs(ptr2))+in_dfps(i,j,j2,1:nGs(ptr2))
                        END IF
                    END IF
                END DO
            END DO
        END DO

        !DO i=1,nAtoms
        !  ptr = supersymbols(i)
        !  DO j2 = 1, 3
        !    DO j = 1, sub_nneighbors(i)+1
        !      IF (j == 1) THEN
        !        print *, 'fp primes at', i, i,j2
        !        print *, out_dfps(i,j,j2,1:nGs(ptr))
        !      ELSE
        !        m=sub_nlist(i,j-1)
        !        IF (m <= nAtoms) THEN
        !          m_loc=find_loc_int(sub_nlist(i,:),m,nAtoms)+1
        !          ptr2=supersymbols(m)
        !          print *, 'fp primes at', i, m,j2
        !          print *, out_dfps(i,m_loc,j2,1:nGs(ptr2))
        !        END IF
        !      END IF  
        !    END DO
        !  END DO
        !END DO
    END SUBROUTINE

END MODULE
