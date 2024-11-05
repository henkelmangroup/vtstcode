MODULE atomsProp

    IMPLICIT NONE
    PUBLIC
    
    INTEGER :: tnAtoms
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: pool_pos_car
    INTEGER, DIMENSION(:), ALLOCATABLE :: pool_ids
    INTEGER, DIMENSION(:), ALLOCATABLE :: supersymbols

    DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: pairs 
    !INTEGER, DIMENSION(:), ALLOCATABLE :: pair_start, pair_end
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: pair_info
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: pair_indices 
    !INTEGER, DIMENSION(:,:), ALLOCATABLE :: pair_ghosts
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: pair_global_indices
    DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: unitvects_pair
    INTEGER, DIMENSION(:), ALLOCATABLE :: gvects
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: tneighs
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: tneighs_incell
    INTEGER, DIMENSION(:), ALLOCATABLE :: tnum_neigh
    INTEGER, DIMENSION(:), ALLOCATABLE :: num_eachpair

    CONTAINS 

    SUBROUTINE dummynn()
    END SUBROUTINE

END MODULE
