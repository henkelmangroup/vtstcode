MODULE pyamffType

    IMPLICIT NONE
    TYPE :: ImgInfo
        INTEGER :: natoms
        INTEGER, DIMENSION(:), ALLOCATABLE :: nneighbors, symbols, natoms_arr
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: neighborlists, atom_idx
        DOUBLE PRECISION :: predE, targetE
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: poscar, cell
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: predF, targetF, Free
        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: input_fps
        DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: input_dfps
        DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE :: layer_backgrad
        DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: in_gradients
        DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE :: hid_gradients
        DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: out_gradients
    END TYPE

    TYPE(ImgInfo), DIMENSION(:), ALLOCATABLE :: TrainImg
    !for output
    INTEGER :: final_epoch
    DOUBLE PRECISION :: final_fRMSE, final_eRMSE, final_loss, final_gradnorm
    
    CONTAINS 

    SUBROUTINE dummyPyamff()
    END SUBROUTINE

END MODULE
