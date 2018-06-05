PROGRAM Crystal_Growth
  ! This program uses an FTCS scheme to model dendritic crystal growth in 2D
  !
  ! Created by Kieran Fitzmaurice on June 4th 2018
  ! Based on model by Ryo Kobayashi published in
  ! Physica D 63 (1993) 410-423

  USE output_arrays
  IMPLICIT NONE

  REAL(KIND = 8), PARAMETER :: PI = 3.14159265359

  INTEGER :: i, j, mode_num, M, N, i_nuc, j_nuc, t_steps, numframes, interval, x, y
  REAL(KIND = 8) :: Lx, Ly, K, delta
  REAL(KIND = 8) :: epsilon_avg, theta, theta_0, tau, alpha, gamma, amp, dt, h
  REAL(KIND = 8) :: temp_0,temp_eq, t, t_0, t_f
  REAL(KIND = 8), ALLOCATABLE :: epsilon_value(:,:), epsilon_prime(:,:), noise(:,:)
  REAL(KIND = 8), ALLOCATABLE :: p(:,:), dpdt(:,:), temp(:,:), dTdt(:,:)
  REAL(KIND = 8), ALLOCATABLE :: dedx(:,:), dedy(:,:), grad_e2(:,:)
  REAL(KIND = 8), ALLOCATABLE :: time(:), temperature(:,:,:), phase(:,:,:)
  CHARACTER(LEN = 125) :: time_file, temp_file, phase_file

  time_file = '/Users/kieranfitzmaurice/Documents/REU_2018/PDE_data/time_file.dat'
  temp_file = '/Users/kieranfitzmaurice/Documents/REU_2018/PDE_data/temp_file.dat'
  phase_file = '/Users/kieranfitzmaurice/Documents/REU_2018/PDE_data/phase_file.dat'

  ! Dimension of array to specify when taking partial derivative w.r.t. x or y
  x = 2
  y = 1

  ! System size
  Lx = 9.0
  Ly = 9.0

  ! Mesh Numbers
  M = 400
  N = 400

  ! Mesh spacing and timestep
  h = Lx/N
  dt = 0.0001

  t_0 = 0.0
  t_f = 1.0

  numframes = 300! Number of frames to save for visualization +1 initial config
  interval = FLOOR((t_f - t_0)/(dt*numframes)) ! Save every nth frame
  t_steps = interval*numframes ! Number of time steps

  ! Fixed Parameters
  epsilon_avg = 0.01  ! Computational parameter for thickness of interface
  tau = 0.0003        ! tau = 1/D where D is diffusion coefficient of p
  alpha = 0.9         ! Parameter influencing m(T)
  gamma = 10.0        ! Parameter influencing m(T)
  amp = 0.01          ! Amplitude of random noise at interface
  temp_0 = 0.2        ! Initial temperature
  temp_eq = 1.0       ! Equilibrium temperature
  theta_0 = PI/4       ! Parameter influencing crystal orientation

  ! Varied Parameters
  K = 1.6             ! Dimensionless latent heat
  mode_num = 6        ! Mode number of anisotropy
  delta = 0.04        ! Strength of anisotropy

  ! location of crystal seed
  i_nuc = FLOOR(REAL(M)/2)
  j_nuc = FLOOR(REAL(N)/2)

  ! Allocate arrays
  ALLOCATE(p(M,N))
  ALLOCATE(dpdt(M,N))
  ALLOCATE(temp(M,N))
  ALLOCATE(dTdt(M,N))
  ALLOCATE(epsilon_value(M,N))
  ALLOCATE(epsilon_prime(M,N))
  ALLOCATE(noise(M,N))
  ALLOCATE(dedx(M,N))
  ALLOCATE(dedy(M,N))
  ALLOCATE(grad_e2(M,N))
  ALLOCATE(time(0:numframes)) ! index 0 for initial config
  ALLOCATE(temperature(1:M,1:N,0:numframes))
  ALLOCATE(phase(1:M,1:N,0:numframes))

  ! Initial Conditions
  DO j = 1,N
     DO i = 1,M
        p(i,j) = 0.0 ! Start as liquid
        temp(i,j) = temp_0 ! Begin at cooling temperature
        theta = ATAN(REAL(j - j_nuc)/(i - i_nuc))
        epsilon_value(i,j) = epsilon_avg*(1 + delta*COS(mode_num*(theta - theta_0)))
        epsilon_prime(i,j) = epsilon_avg*delta*(-1*mode_num*SIN(mode_num*(theta - theta_0)))

     ENDDO
  ENDDO

  ! must do this to avoid creating NaN
  epsilon_value(i_nuc,j_nuc) = 0.0
  epsilon_prime(i_nuc,j_nuc) = 0.0

  dedx = partial_2D(epsilon_value*epsilon_prime,h,x)
  dedy = partial_2D(epsilon_value*epsilon_prime,h,y)
  grad_e2 = gradient_2D(epsilon_value**2,h)

  ! Seed crystal
  p(i_nuc,j_nuc) = 1.0

  ! Record initial state of system
  t = t_0
  time(0) = t
  temperature(:,:,0) = temp
  phase(:,:,0) = p

  ! Use FTCS scheme to advance system in time

  DO i = 1,numframes
     DO j = 1,interval

        ! Calculate dp/dt

        dpdt = (-1*dedx*partial_2D(p,h,y) + dedy*partial_2D(p,h,x) &
        + grad_e2*gradient_2D(p,h) + epsilon_value**2*laplacian_2D(p,h) &
        + p*(1 - p)*(p - 0.5 + m_value(temp,temp_eq,gamma,alpha)))/tau

        ! Add random noise at interface
        CALL INIT_RANDOM_SEED()
        CALL RANDOM_NUMBER(noise)
        noise = amp*p*(1-p)*(noise - 0.5)
        dpdt = dpdt + noise

        ! Calculate dT/dt
        dTdt = laplacian_2D(temp,h) + K*dpdt

        ! Advance in time
        p = p + dt*dpdt
        temp = temp + dt*dTdt
        t = t + dt

     ENDDO
     ! Record state of system
     time(i) = t
     temperature(:,:,i) = temp
     phase(:,:,i) = p
  ENDDO

  ! Write results to file
  CALL output1D_binary(time_file,time,numframes + 1)
  CALL output3D_binary(temp_file,temperature,M,N,numframes + 1)
  CALL output3D_binary(phase_file,phase,M,N,numframes + 1)

CONTAINS

  !______________________1st Order Partial Derivative__________________________!

  FUNCTION partial_2D(f,h,dim) RESULT(partial)
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: f
    REAL(KIND = 8), INTENT(IN) :: h
    INTEGER, INTENT(IN) :: dim
    REAL(KIND = 8), DIMENSION(1:SIZE(f,1),1:SIZE(f,2)) :: partial

    partial = (CSHIFT(f,1,dim) - CSHIFT(f,-1,dim))/(2*h)

  END FUNCTION partial_2D

  !______________________1st Order Central Difference__________________________!

  FUNCTION gradient_2D(f,h) RESULT(del_f)
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: f
    REAL(KIND = 8), INTENT(IN) :: h
    REAL(KIND = 8), DIMENSION(1:SIZE(f,1),1:SIZE(f,2)) :: del_f

    del_f = (CSHIFT(f,1,1) - CSHIFT(f,-1,1) + CSHIFT(f,1,2) - CSHIFT(f,-1,2))/(2*h)

  END FUNCTION gradient_2D

  !______________________2nd Order Central Difference__________________________!

  FUNCTION laplacian_2D(f,h) RESULT(del_squared_f)
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: f
    REAL(KIND = 8), INTENT(IN) :: h
    REAL(KIND = 8), DIMENSION(1:SIZE(f,1),1:SIZE(f,2)) :: del_squared_f

    del_squared_f = (CSHIFT(f,1,1)+CSHIFT(f,-1,1)+CSHIFT(f,1,2)+CSHIFT(f,-1,2) - 4*f)/h**2

  END FUNCTION laplacian_2D

  !_______________Temperature Dependence of Free Energy________________________!

  FUNCTION m_value(temp,temp_eq,gamma,alpha) RESULT(m)
    REAL(KIND = 8), PARAMETER :: PI = 3.14159265359
    REAL(KIND = 8), DIMENSION(:,:), INTENT(IN) :: temp
    REAL(KIND = 8), INTENT(IN) :: temp_eq, gamma, alpha
    REAL(KIND = 8), DIMENSION(1:SIZE(temp,1),1:SIZE(temp,2)) :: m

    m = alpha/PI*ATAN(gamma*(temp_eq - temp))

  END FUNCTION m_value

  !__________Generates a new random seed each time it is invoked_______________!

  ! Taken from GNU Fortran compiler manual
  ! https://gcc.gnu.org/onlinedocs/gcc-4.7.4/gfortran/

  SUBROUTINE INIT_RANDOM_SEED()
    IMPLICIT NONE
    INTEGER, ALLOCATABLE :: seed(:)
    INTEGER :: i, n, un, istat, dt(8), pid, t(2), s
    INTEGER(KIND = 8) :: count, tms

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    ! First try if the OS provides a random number generator
    OPEN(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    IF (istat == 0) THEN
       READ(un) seed
       CLOSE(un)
    ELSE
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       CALL SYSTEM_CLOCK(count)
       IF (count /= 0) THEN
          t = TRANSFER(count, t)
       ELSE
          CALL DATE_AND_TIME(values=dt)
          tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24 * 60 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
          t = TRANSFER(tms, t)
       END IF
       s = IEOR(t(1), t(2))
       pid = getpid() + 1099279 ! Add a prime
       s = IEOR(s, pid)
       IF (n >= 3) THEN
          seed(1) = t(1) + 36269
          seed(2) = t(2) + 72551
          seed(3) = pid
          IF (n > 3) THEN
             seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
          END IF
       ELSE
          seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
       END IF
    END IF
    CALL RANDOM_SEED(put=seed)
  END SUBROUTINE INIT_RANDOM_SEED

END PROGRAM Crystal_Growth
