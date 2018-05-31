PROGRAM Diffusion
! This program uses an FTCS scheme to solve the 1D Diffusion Equation

USE OutputArrays
IMPLICIT NONE

INTERFACE
   SUBROUTINE u_analytic(x,x_0,t,t_0,u)
     REAL(KIND = 8), INTENT(IN) :: x_0, t, t_0
     REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: x
     REAL(KIND = 8), DIMENSION(:), INTENT(OUT) :: u
   END SUBROUTINE u_analytic
END INTERFACE

INTEGER :: i, n, t_steps
REAL(KIND = 8) :: k, l, dx, dt, max_dt,t_0,t_f,x_0
REAL(KIND = 8), ALLOCATABLE :: x(:), u(:), concentration(:,:)
CHARACTER(LEN = 25) :: filepath,filename

filepath = '/Users/kieranfitzmaurice/Documents/GitHub/PDEs/'
filename = 'diffusion_example'

k = 1 ! Diffusion coefficient
l = 20 ! Length of rod
x_0 = l/2 ! Peak of concentration at beginning (used for analytical solution)
n = 100 ! Number of points in rod

dx = l/(n-1)
max_dt = (0.5*dx**2)/k ! Von Neumman Stability
dt = 0.4*max_dt

t_0 = 0.1
t_f = 5

t_steps = FLOOR((t_f - t_0)/dt)

! Array concentration stores how u changes in time
ALLOCATE(x(1:n))
ALLOCATE(u(1:n))
ALLOCATE(concentration(1:n,1:t_steps))

! Store x positions
DO i = 1,n
  x(i) = (i-1)*dx
ENDDO

! Get initial concentration profile from analytical solution
CALL u_analytic(x,x_0,t_0,t_0,u)

! Use FTCS method to numerically solve diffusion
DO i = 1,t_steps
  u = u + dt*(CSHIFT(u,SHIFT = -1) - 2*u + CSHIFT(u,SHIFT = 1))/dx**2
  concentration(:,i) = u
ENDDO

! Write concentration profiles to output file
CALL output2D_txt(filepath,filename,concentration,t_steps,n)

END PROGRAM Diffusion

!******************************************************************************!

SUBROUTINE u_analytic(x,x_0,t,t_0,u)
! Analytical solution to diffusion equation for this system
! Used to generate initial concentration profile

IMPLICIT NONE

REAL(KIND = 8), INTENT(IN) :: x_0, t, t_0
REAL(KIND = 8), DIMENSION(:), INTENT(IN) :: x
REAL(KIND = 8), DIMENSION(:), INTENT(OUT) :: u

u = (t_0/t)**0.5 * EXP(-(x-x_0)**2/(4*t))

END SUBROUTINE u_analytic
