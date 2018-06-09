MODULE array_io
CONTAINS

  !****************************************************************************!

  SUBROUTINE size1D_binary(filename,n1)
    ! This subroutine reads in the length of a 1D array from a binary file.
    ! This binary file should follow the same format as those that are written
    ! by the array output subroutines in this module.

    CHARACTER(LEN = *), INTENT(IN) :: filename
    INTEGER, INTENT(OUT):: n1

    OPEN(UNIT=1,FILE=filename,FORM ='UNFORMATTED',STATUS='OLD',ACTION='READ')
    READ(1) n1
    CLOSE(1)

  END SUBROUTINE size1D_binary

  !****************************************************************************!

  SUBROUTINE size2D_binary(filename,n1,n2)
    ! This subroutine reads in the dimensions of a 2D array from a binary file.
    ! This binary file should follow the same format as those that are written
    ! by the array output subroutines in this module.

    CHARACTER(LEN = *), INTENT(IN) :: filename
    INTEGER, INTENT(OUT):: n1, n2

    OPEN(UNIT=1,FILE=filename,FORM ='UNFORMATTED',STATUS='OLD',ACTION='READ')
    READ(1) n1
    READ(1) n2
    CLOSE(1)

  END SUBROUTINE size2D_binary

  !****************************************************************************!

  SUBROUTINE size3D_binary(filename,n1,n2,n3)
    ! This subroutine reads in the dimensions of a 3D array from a binary file.
    ! This binary file should follow the same format as those that are written
    ! by the array output subroutines in this module.

    CHARACTER(LEN = *), INTENT(IN) :: filename
    INTEGER, INTENT(OUT):: n1, n2, n3

    OPEN(UNIT=1,FILE=filename,FORM ='UNFORMATTED',STATUS='OLD',ACTION='READ')
    READ(1) n1
    READ(1) n2
    READ(1) n3
    CLOSE(1)

  END SUBROUTINE size3D_binary

  !****************************************************************************!

  SUBROUTINE size4D_binary(filename,n1,n2,n3,n4)
    ! This subroutine reads in the dimensions of a 4D array from a binary file.
    ! This binary file should follow the same format as those that are written
    ! by the array output subroutines in this module.

    CHARACTER(LEN = *), INTENT(IN) :: filename
    INTEGER, INTENT(OUT):: n1, n2, n3, n4

    OPEN(UNIT=1,FILE=filename,FORM ='UNFORMATTED',STATUS='OLD',ACTION='READ')
    READ(1) n1
    READ(1) n2
    READ(1) n3
    READ(1) n4
    CLOSE(1)

  END SUBROUTINE size4D_binary

  !****************************************************************************!

  FUNCTION input1D_binary(filename,n1) RESULT(A)
    ! This subroutine reads in a 1D array from a binary file
    CHARACTER(LEN = *), INTENT(IN) :: filename
    INTEGER, INTENT(IN):: n1
    REAL(KIND = 8), DIMENSION(n1) :: A

    OPEN(UNIT=1,FILE=filename,FORM ='UNFORMATTED',STATUS='OLD',ACTION='READ')
    READ(1) ! Skip entry specifying value of n1
    READ(1) A(1:n1)
    CLOSE(1)

  END FUNCTION input1D_binary

  !****************************************************************************!

  FUNCTION input2D_binary(filename,n1,n2) RESULT(A)
    ! This subroutine reads in a 2D array in column order from a binary file
    CHARACTER(LEN = *), INTENT(IN) :: filename
    INTEGER, INTENT(IN):: n1, n2
    REAL(KIND = 8), DIMENSION(n1,n2) :: A

    OPEN(UNIT=1,FILE=filename,FORM ='UNFORMATTED',STATUS='OLD',ACTION='READ')
    READ(1) ! Skip entry specifying value of n1
    READ(1) ! Skip entry specifying value of n2
    READ(1) A(1:n1,1:n2)

    CLOSE(1)

  END FUNCTION input2D_binary

  !****************************************************************************!

  FUNCTION input3D_binary(filename,n1,n2,n3) RESULT(A)
    ! This subroutine reads in a 3D array in column order from a binary file
    CHARACTER(LEN = *), INTENT(IN) :: filename
    INTEGER, INTENT(IN):: n1, n2, n3
    REAL(KIND = 8), DIMENSION(n1,n2,n3) :: A

    OPEN(UNIT=1,FILE=filename,FORM ='UNFORMATTED',STATUS='OLD',ACTION='READ')
    READ(1) ! Skip entry specifying value of n1
    READ(1) ! Skip entry specifying value of n2
    READ(1) ! Skip entry specifying value of n3
    READ(1) A(1:n1,1:n2,1:n3)
    CLOSE(1)

  END FUNCTION input3D_binary

  !****************************************************************************!

  FUNCTION input4D_binary(filename,n1,n2,n3,n4) RESULT(A)
    ! This subroutine reads in a 4D array in column order from a binary file
    CHARACTER(LEN = *), INTENT(IN) :: filename
    INTEGER, INTENT(IN):: n1, n2, n3, n4
    REAL(KIND = 8), DIMENSION(n1,n2,n3,n4) :: A

    OPEN(UNIT=1,FILE=filename,FORM ='UNFORMATTED',STATUS='OLD',ACTION='READ')
    READ(1) ! Skip entry specifying value of n1
    READ(1) ! Skip entry specifying value of n2
    READ(1) ! Skip entry specifying value of n3
    READ(1) ! Skip entry specifying value of n4
    READ(1) A(1:n1,1:n2,1:n3,1:n4)
    CLOSE(1)

  END FUNCTION input4D_binary

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
  END SUBROUTINE output1D_binary

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
  END SUBROUTINE output2D_binary

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
  END SUBROUTINE output3D_binary

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
  END SUBROUTINE output4D_binary

  !****************************************************************************!

END MODULE
