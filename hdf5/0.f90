PROGRAM main

        USE MPI
        USE HDF5
        IMPLICIT NONE

!    INCLUDE 'mpif.h'

    !----------------------------------------------------------------------------!
    !  localities
    !----------------------------------------------------------------------------!
    INTEGER                 :: info
    INTEGER                 :: nproc, rank
    INTEGER                 :: i
    CHARACTER(len=256)      :: filename
    INTEGER(hid_t)          :: fapl_id                         ! file access identifier
    INTEGER(hid_t)          :: file_id                         ! file identifier
    INTEGER                 :: hdferror, error                 ! HDF hdferror flag
    INTEGER                 :: filename_length

    CALL MPI_Init(info)
    CALL MPI_Comm_Size(MPI_COMM_WORLD, nproc, info)
    CALL MPI_Comm_Rank(MPI_COMM_WORLD, rank, info)

    filename = 'out.h5'
    filename_length = LEN_TRIM(filename)

    CALL h5open_f(error)
    
    IF (rank.EQ.0) THEN

        CALL h5fcreate_f(filename(1:filename_length), H5F_ACC_TRUNC_F, file_id, hdferror)
        IF (hdferror.LT.0) WRITE(*,*) 'file creation failed'
        CALL h5fclose_f(file_id, hdferror)
        IF (hdferror.LT.0) WRITE(*,*) 'file closing failed'

    ENDIF



    DO i = 1, 200

        CALL h5pcreate_f(H5P_FILE_ACCESS_F, fapl_id, hdferror)
        IF (hdferror.LT.0) WRITE(*,*) 'fapl_id creation failed'

        CALL h5pset_fapl_mpio_f(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdferror)
        IF (hdferror.LT.0) WRITE(*,*) 'fapl_id modification failed'

        CALL h5fopen_f(filename(1:filename_length), H5F_ACC_RDWR_F, file_id, hdferror, access_prp = fapl_id)
        IF (hdferror.LT.0) WRITE(*,*) 'file opening failed'
        
        CALL h5pclose_f(fapl_id, hdferror)
        IF (hdferror.LT.0) WRITE(*,*) 'fapl_id closing failed'

        CALL h5fclose_f(file_id, hdferror)
        IF (hdferror.LT.0) WRITE(*,*) 'file closing failed'

        IF (rank.EQ.0) WRITE(*,*) i

    ENDDO

    CALL h5close_f(error)
    CALL MPI_Finalize(info)

END PROGRAM main
