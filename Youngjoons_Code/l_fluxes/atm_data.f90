MODULE atm_data

  IMPLICIT NONE

  ! dp = Positive double which represents the step in p-axis
  ! dx = Positive double which represents the step in x-axis
  ! dt = Positive double which represents the step in time
  ! t = double, time at which the computation occurs
  ! theta : double, parameter for the central-upwind fluxes
  ! DEBUG = integer, if 0 no debug, if 1 debug is on
  ! Np = integer equals to the number of cell for the mesh along the pressure p
  ! Nx = integer equals to the number of cell for the mesh along the x
  ! p = tabular of Np+1 double, p(i) is the i times p-coordinate 
  !		for 1<=i<=Np+1 of the mesh
  ! pm = tabular of Np double, pm(i) is the i times center 
  !		p-coordinate for 1<=i<=Np of the cell i
  ! x = tabular of Nx+1 double, x(i) is the i times x-coordinate 
  !		for 1<=i<=Nx+1 of the mesh
  ! xm = tabular of Nx double, xm(i) is the i times center 
  !		x-coordinate for 1<=i<=Nx of the cell i
  ! A = double of dimension 2*Nx*Np,
  !     A(l,i,j)= value of T' (if l=1) the temperature, or q (if l=2) specific
  !             humidity, on the cell (i,j) of the mesh, i=1...Nx, j=1...Np
  ! DA= double of dimension 2*Nx*Np, DA~dA/dt
  ! Tbar = double of dimension Nx*Np, value of \bar{T} on the mesh
  ! w = double of dimension Nx*(Np+1), value of omega on the alternate mesh
  ! u = double of dimension (Nx+1)*Np, value of the velocity on the alternate mesh
  ! du = double of dimension (Nx+1)*Np, du~du/dt
  ! S = double of DIMENSION 2*Nx*Np which represents the values of 
  !	   transpose(ST',Sq) at time t
  !	   S(l,i,j) is the value of the l-time coordonnate
  !	   on the cell (i,j) for i=1,..,N for the source term at t
  ! Sdt = double of DIMENSION 2*Nx*Np, Sdt~dS/dt(t+dt)
  ! Sdt2 = double of DIMENSION 2*Nx*Np, Sdt2~dS/dt(t+dt/2)		

  DOUBLE PRECISION :: dx, dt, t, deltaT, p0, Ep_dump
  DOUBLE PRECISION :: th1,th2
  INTEGER :: DEBUG
  CHARACTER*10 :: snt
  INTEGER :: nt

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: A, DA
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: S,Sdt,Sdt2
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: S2, S2dt, S2dt2
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Tbar, phi, normal, phix
  INTEGER :: Nx, Np
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: p, x, ac, ac2, ac3, rqs
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: pm, xm, pm2, pm3, xm2, xm3, zm
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: FluxHori, FluxVert
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: w, Ab, DAb,cu,phitemp,phixtemp
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dp, hhori, lda, ldax, ld,ldtemp
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: m, coeff, m2, m3, coeff3
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: dux,qs_temp
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: rS,rW
  
  INTEGER :: bc, scheme
  INTEGER :: PHYSIC

  DOUBLE PRECISION :: R, Rv, Cp, g

  PARAMETER (R=287.d0)
  PARAMETER (Cp=1004.d0)
  PARAMETER (Rv=461.50d0)

CONTAINS

INCLUDE 'atm_functions.f90'

END MODULE atm_data
