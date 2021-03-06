PROGRAM Cahn_Hilliard
  ! This program uses an FTCS scheme to model Spinodal Decomposition in 2D

  USE array_io ! Array I/O module
  IMPLICIT NONE

  INTEGER :: i, j, Nx, Ny, t_steps, numframes, interval
  REAL(KIND = 8) :: h, dt, epsilon, M, W, t_0, t_f, t, c_avg, noise
  REAL(KIND = 8), ALLOCATABLE :: c(:,:), concentration(:,:,:), time(:)
  CHARACTER(LEN = 100) :: filename

  filename = '/Users/kieranfitzmaurice/Documents/MATLAB/Spinodal.dat'

  h = 0.1        ! Lattice spacing
  dt = 0.0003    ! Timestep

  epsilon = 0.1  ! Gradient energy coefficient
  M = 1.0        ! Mobility
  W = 1.0        ! Double well potential height parameter

  t_0 = 0.0      ! Elapsed time and number of time steps
  t_f = 5.0

  numframes = 300! Number of frames to save for visualization
  interval = FLOOR((t_f - t_0)/(dt*numframes)) ! Save every nth frame
  t_steps = interval*numframes ! Number of time steps

  Nx = 250       ! Number of lattice points in x-direction
  Ny = 250       ! Number of lattice points in y-direction

  c_avg = 0.5    ! Average concentration in system
  noise = 0.1    ! Strength of random noise

  ! c is current profile, concentration stores history
  ALLOCATE(c(Ny,Nx))
  ALLOCATE(concentration(Ny,Nx,numframes))
  ALLOCATE(time(numframes))

  ! Initial conditions (homogeneous solution with random fluctuations)
  CALL RANDOM_NUMBER(c)
  c = c_avg + (2*c - 1)*noise

  ! Use FTCS scheme to advance system in time according to Cahn-Hilliard Eq.
  t = t_0
  DO i = 1,numframes
    DO j = 1,interval
      c = c + dt*M*central_diff_2D(grad_FE(c,h,W,epsilon),h)
      t = t + dt
    ENDDO
    concentration(:,:,i) = c
    time(i) = t
  ENDDO

  ! Write results to file
  CALL output3D_binary(filename,concentration,Ny,Nx,numframes)


!****************************END OF MAIN PROGRAM*******************************!

CONTAINS

  !______________________CENTRAL DIFFERENCE FUNCTION___________________________!

  FUNCTION central_diff_2D(f,h) RESULT(laplacian)
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: f
    REAL(KIND = 8), INTENT(IN) :: h
    REAL(KIND = 8), DIMENSION(1:SIZE(f,1),1:SIZE(f,2)) :: laplacian

    laplacian = (CSHIFT(f,1,1)+CSHIFT(f,-1,1)+CSHIFT(f,1,2)+CSHIFT(f,-1,2)-4*f)/h**2

  END FUNCTION central_diff_2D

  !______________GRADIENT OF FREE ENERGY DENSITY FUNCTION______________________!

  FUNCTION grad_FE_density(c,W) RESULT(df_dc)
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: c
    REAL(KIND = 8), INTENT(IN) :: W
    REAL(KIND = 8), DIMENSION(1:SIZE(c,1),1:SIZE(c,2)) :: df_dc

    df_dc = 0.5*W*c*(1-c)*(1-2*c)

  END FUNCTION grad_FE_density

  !______________GRADIENT OF SYSTEM FREE ENERGY FUNCTION______________________!

  FUNCTION grad_FE(c,h,W,epsilon) RESULT(dF_dc)

    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: c
    REAL(KIND = 8), INTENT(IN) :: h, W, epsilon
    REAL(KIND = 8), DIMENSION(1:SIZE(c,1),1:SIZE(c,2)) :: dF_dc

    dF_dc = grad_FE_density(c,W) - epsilon**2*central_diff_2D(c,h)

  END FUNCTION grad_FE

END PROGRAM Cahn_Hilliard
