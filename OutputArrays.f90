MODULE OutputArrays
CONTAINS

!******************************************************************************!

  SUBROUTINE output2D_txt(filepath,filename,A,nrows,ncols)
  ! This subroutine writes a 2D array to a .txt file
  IMPLICIT NONE

  CHARACTER(LEN = *),INTENT(IN) :: filepath, filename
  REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: A
  INTEGER, INTENT(IN) :: nrows, ncols
  INTEGER :: i,j
  CHARACTER(LEN = 100) :: full_file
  CHARACTER(LEN = 20) :: frmt

  frmt = '(ES14.7)'
  full_file = TRIM(ADJUSTL(filepath))//TRIM(ADJUSTL(filename))//'.txt'

  OPEN(UNIT=1,FILE=full_file,FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')

  DO j = 1,ncols
    DO i = 1,nrows
      WRITE(1,*) SNGL(A(i,j))
    ENDDO
  ENDDO

  CLOSE(1)
  END SUBROUTINE

!******************************************************************************!

  SUBROUTINE output2D_binary(filepath,filename,A,nrows,ncols)
  ! This subroutine writes a 2D array as a binary file
  IMPLICIT NONE

  CHARACTER(LEN = *),INTENT(IN) :: filepath, filename
  REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: A
  INTEGER, INTENT(IN) :: nrows, ncols
  CHARACTER(LEN = 100) :: full_file

  full_file = TRIM(ADJUSTL(filepath))//TRIM(ADJUSTL(filename))//'.dat'

  OPEN(UNIT=1,FILE=full_file,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')

  WRITE(1) A(1:nrows,1:ncols)

  CLOSE(1)
  END SUBROUTINE

!******************************************************************************!

END MODULE
