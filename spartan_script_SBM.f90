PROGRAM Spartan_Script_SBM
  ! This program demonstrates the level set method of smoothing, and applies the
  ! smoothed-boundary method to solving diffusion in a domain defined by the
  ! Michigan State script logo

  USE array_io
  IMPLICIT NONE

  INTEGER :: i, j, M, N, n_avg
  REAL(KIND = 8) :: h, dtau, dt, gamma, delta, sum_avg
  REAL(KIND = 8), ALLOCATABLE :: image(:,:), phi(:,:)
  REAL(KIND = 8), ALLOCATABLE :: smoothing(:,:,:)
  LOGICAL, ALLOCATABLE :: mymask(:,:)
  CHARACTER(LEN = 125) :: input_file, output_file

  input_file = '/Users/kieranfitzmaurice/Documents/REU_2018/Spartan_Script.dat'
  output_file = '/Users/kieranfitzmaurice/Documents/REU_2018/Spartan_SBM.dat'

  ! Dimensions of input array
  CALL size2D_binary(input_file,M,N)

  ALLOCATE(image(M,N))
  ALLOCATE(phi(M,N))
  ALLOCATE(mymask(M,N))

  image = input2D_binary(input_file,M,N)

  ! Mesh spacing
  h = 1

  ! Level set method parameters (unphysical)
  dtau = 0.1*h**2
  gamma = 3*(h**2)**0.5
  delta = 3*h ! controls width of interface. ~ 3*h is good

  ! Scale so phi varies from -delta/2 to delta/2
  phi = image*delta - delta/2

  ! Level set method (Gererates signed distance function near boundary)
  DO WHILE (MAXVAL(phi) <= 6*h)
     phi = phi + dtau*smoothed_sign(phi,gamma)*(1 - gudonov_upwind(phi,h))
  ENDDO

  ! Diffusion smoothing (Eliminates noise in second derivatives)
  DO i = 1,20
    phi = phi + dtau*laplacian_noflux(phi,h)
  ENDDO

  ! Renormalize
  mymask = ABS(phi).LT.h
  n_avg = COUNT(mymask)
  sum_avg = SUM(ABS(gradient_noflux(phi,h)), MASK = mymask)

  phi = phi/(sum_avg/n_avg)

  CALL output2D_binary(output_file,phi,M,N)

CONTAINS

  !______________Smoothed sign function for level set method_____________________!
  FUNCTION smoothed_sign(phi,gamma) RESULT(S)
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: phi
    REAL(KIND = 8), INTENT(IN) :: gamma
    REAL(KIND = 8), DIMENSION(1:SIZE(phi,1),1:SIZE(phi,2)) :: S

    S = phi/(phi**2 + gamma**2)**0.5

  END FUNCTION smoothed_sign

  !_____________________Returns higest values of two arrays______________________!
  FUNCTION max_array(A,B) RESULT(arr)
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: A, B
    REAL(KIND = 8), DIMENSION(1:SIZE(A,1),1:SIZE(A,2)) :: arr
    INTEGER :: M, N, i, j

    M = SIZE(A,1)
    N = SIZE(A,2)

    DO j = 1,N
       DO i = 1,M

          arr(i,j) = MAX(A(i,j),B(i,j))

       ENDDO
    ENDDO


  END FUNCTION max_array

  !_______________________1st order forward difference___________________________!

  FUNCTION forward_diff(f,h,dim) RESULT(partial_f)
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: f
    REAL(KIND = 8), INTENT(IN) :: h
    INTEGER, INTENT(IN) :: dim
    REAL(KIND = 8), DIMENSION(1:SIZE(f,1),1:SIZE(f,2)) :: partial_f
    INTEGER :: M, N, x, y

    M = SIZE(f,1)
    N = SIZE(f,2)
    x = 2
    y = 1

    IF (dim == x) THEN
       partial_f(:,2:(N - 1)) = (f(:,1:(N - 2)) - f(:,2:(N - 1)))/h
       partial_f(:,1) = 0.0
       partial_f(:,N) = 0.0
    ELSE
       partial_f(2:(M - 1),:) = (f(3:M,:) - f(2:(M - 1),:))/h
       partial_f(1,:) = 0.0 ! No flux boundary conditions
       partial_f(M,:) = 0.0
    ENDIF

  END FUNCTION forward_diff

  !_______________________1st order backward difference___________________________!

  FUNCTION backward_diff(f,h,dim) RESULT(partial_f)
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: f
    REAL(KIND = 8), INTENT(IN) :: h
    INTEGER, INTENT(IN) :: dim
    REAL(KIND = 8), DIMENSION(1:SIZE(f,1),1:SIZE(f,2)) :: partial_f
    INTEGER :: M, N, x, y

    M = SIZE(f,1)
    N = SIZE(f,2)
    x = 2
    y = 1

    IF (dim == x) THEN
       partial_f(:,2:(N - 1)) = (f(:,2:(N - 1)) - f(:,3:N))/h
       partial_f(:,1) = 0.0 ! No flux boundary conditions
       partial_f(:,N) = 0.0
    ELSE
       partial_f(2:(M - 1),:) = (f(2:(M - 1),:) - f(1:(M - 2),:))/h
       partial_f(1,:) = 0.0
       partial_f(M,:) = 0.0
    ENDIF

  END FUNCTION backward_diff

  !_________________1st Order Central Difference (no flux)_____________________!

  FUNCTION gradient_noflux(f,h) RESULT(del_f)
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: f
    REAL(KIND = 8), INTENT(IN) :: h
    REAL(KIND = 8), DIMENSION(1:SIZE(f,1),1:SIZE(f,2)) :: del_f
    INTEGER :: M, N

    M = SIZE(f,1)
    N = SIZE(f,2)

    del_f(2:(M - 1),2:(N - 1)) = (f(2:(M - 1),3:N) - f(2:(M - 1),1:(N - 2)) &
    + f(1:(M - 2),2:(N - 1)) - f(3:M,2:(N - 1)))/(2*h)

    del_f(1,:) = 0.0
    del_f(M,:) = 0.0
    del_f(:,1) = 0.0
    del_f(:,N) = 0.0

  END FUNCTION gradient_noflux

  !_________________2nd Order Central Difference (no flux)_____________________!

  FUNCTION laplacian_noflux(f,h) RESULT(del2_f)
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: f
    REAL(KIND = 8), INTENT(IN) :: h
    REAL(KIND = 8), DIMENSION(1:SIZE(f,1),1:SIZE(f,2)) :: del2_f
    INTEGER :: M, N

    M = SIZE(f,1)
    N = SIZE(f,2)

    del2_f(2:(M - 1),2:(N - 1)) = (f(3:M,2:(N - 1)) + f(1:(M - 2),2:(N - 1)) &
         + f(2:(M - 1),1:(N - 2)) + f(2:(M - 1),3:N) - 4*f(2:(M - 1),2:(N - 1)))/h**2

    ! No flux boundaries
    del2_f(2:(M - 1),1) = (f(1:(M - 2),1) + f(3:M,1) + 2*f(2:(M - 1),2) - 4*f(2:(M - 1),1))/h**2
    del2_f(2:(M - 1),N) = (f(1:(M - 2),N) + f(3:M,N) + 2*f(2:(M - 1),N - 1) - 4*f(2:(M - 1),N))/h**2
    del2_f(1,2:(N - 1)) = (f(1,3:N) + f(1,1:(N - 2)) + 2*f(2,2:(N - 1)) - 4*f(1,2:(N - 1)))/h**2
    del2_f(M,2:(N - 1)) = (f(M,3:N) + f(M,1:(N - 2)) + 2*f(M - 1,2:(N - 1)) - 4*f(M,2:(N - 1)))/h**2
    del2_f(1,1) = (2*f(2,1) + 2*f(1,2) - 4*f(1,1))/h**2
    del2_f(1,N) = (2*f(1,N - 1) + 2*f(2,N) - 4*f(1,N))/h**2
    del2_f(M,N) = (2*f(M,N - 1) + 2*f(M - 1,N) - 4*f(M,N))/h**2
    del2_f(M,1) = (2*f(M,2) + 2*f(M - 1,1) - 4*f(M,1))/h**2

  END FUNCTION laplacian_noflux

  !__________________________Gudonov upwind scheme_____________________________!
  FUNCTION gudonov_upwind(phi,h) RESULT(abs_grad_phi)
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: phi
    REAL(KIND = 8), INTENT(IN) :: h
    REAL(KIND = 8), DIMENSION(1:SIZE(phi,1),1:SIZE(phi,2)) :: zeros
    REAL(KIND = 8), DIMENSION(1:SIZE(phi,1),1:SIZE(phi,2)) :: abs_grad_phi
    INTEGER :: x, y

    x = 2
    y = 1

    ! Array of zeros
    zeros = 0.0

    abs_grad_phi = ( &
         max_array(max_array(backward_diff(phi,h,x),zeros)**2, &
         max_array(-1*forward_diff(phi,h,x),zeros)**2) &
         + &
         max_array(max_array(backward_diff(phi,h,y),zeros)**2, &
         max_array(-1*forward_diff(phi,h,y),zeros)**2) )**0.5

  END FUNCTION gudonov_upwind

END PROGRAM Spartan_Script_SBM
