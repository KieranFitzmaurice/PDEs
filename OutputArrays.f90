MODULE OutputArrays
CONTAINS

  SUBROUTINE output2D_txt(filepath,filename,A,nx,ny)
  ! This subroutine writes a 2D array to a .txt file
  IMPLICIT NONE

  CHARACTER(LEN = *),INTENT(IN) :: filepath, filename
  REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: A
  INTEGER, INTENT(IN) :: nx, ny
  CHARACTER(LEN = 100) :: full_file

  full_file = TRIM(filepath)//TRIM(filename)//'.txt'

  OPEN(UNIT=1,FILE=full_file,FORM='FORMATTED',STATUS='REPLACE',ACTION='WRITE')

  WRITE(1,*) SNGL(A(1:nx,1:ny))

  CLOSE(1)
  END SUBROUTINE

END MODULE
