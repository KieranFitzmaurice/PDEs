MODULE output_arrays
CONTAINS

  !****************************************************************************!

  SUBROUTINE output1D_binary(filename,A,n1)
    ! This subroutine writes a 1D array as a binary file
    IMPLICIT NONE

    CHARACTER(LEN = *),INTENT(IN) :: filename
    REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: A
    INTEGER, INTENT(IN) :: n1

    OPEN(UNIT=1,FILE=filename,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')

    WRITE(1) n1
    WRITE(1) A(1:n1)

    CLOSE(1)
  END SUBROUTINE

  !****************************************************************************!

  SUBROUTINE output2D_binary(filename,A,n1,n2)
    ! This subroutine writes a 2D array as a binary file
    IMPLICIT NONE

    CHARACTER(LEN = *),INTENT(IN) :: filename
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: A
    INTEGER, INTENT(IN) :: n1, n2

    OPEN(UNIT=1,FILE=filename,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')

    WRITE(1) n1
    WRITE(1) n2
    WRITE(1) A(1:n1,1:n2)

    CLOSE(1)
  END SUBROUTINE

  !****************************************************************************!

  SUBROUTINE output3D_binary(filename,A,n1,n2,n3)
    ! This subroutine writes a 3D array as a binary file
    IMPLICIT NONE

    CHARACTER(LEN = *),INTENT(IN) :: filename
    REAL(KIND = 8), DIMENSION(:,:,:), INTENT(IN) :: A
    INTEGER, INTENT(IN) :: n1, n2, n3

    OPEN(UNIT=1,FILE=filename,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')

    WRITE(1) n1
    WRITE(1) n2
    WRITE(1) n3
    WRITE(1) A(1:n1,1:n2,1:n3)

    CLOSE(1)
  END SUBROUTINE

  !****************************************************************************!

  SUBROUTINE output4D_binary(filename,A,n1,n2,n3,n4)
    ! This subroutine writes a 4D array as a binary file
    IMPLICIT NONE

    CHARACTER(LEN = *),INTENT(IN) :: filename
    REAL(KIND = 8), DIMENSION(:,:,:,:), INTENT(IN) :: A
    INTEGER, INTENT(IN) :: n1, n2, n3, n4

    OPEN(UNIT=1,FILE=filename,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')

    WRITE(1) n1
    WRITE(1) n2
    WRITE(1) n3
    WRITE(1) n4
    WRITE(1) A(1:n1,1:n2,1:n3,1:n4)

    CLOSE(1)
  END SUBROUTINE

  !****************************************************************************!

END MODULE
