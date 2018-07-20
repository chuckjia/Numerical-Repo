!	Goal : Find the FluxVert (North and South flux defined by (3.10) in PDF file.
!	Parameters 
!			FluxVert is (2,Nx,Np) matrix. 
!			FluxVert(1,i,j) = Flux for T(i+1,j+1)=Tij, FluxVert(2,i,j) = Flux for q(i+1,j+1)=qij
!			FluxNorth : North flux i.e. G(i,j+0.5), double vector size 2 for T and q
!   		FluxSouth : South flux i.e. G(i,j-0.5), double vector size 2 for T and q
!			ac(i,j) : area for (i,j) mesh.
!			Nx, Np : the number of elements of x and p, respectively.
SUBROUTINE dfdp()

  IMPLICIT NONE

  INTEGER :: i,j
  DOUBLE PRECISION :: FluxNorth(3), FluxSouth(3)

  Do i=1,Nx
     Do j=1,Np
		! Top and bottom parts of mesh are different.
		
		! T(i+1,2) is related to the bottom of the mesh
		! FluxSouth is zero because w is zero at the bottom
        IF (j .EQ. 1) THEN
           FluxNorth=FluxS(i,j+1)
           FluxSouth=0.d0
        
        ! T(i+1,j+1) is related to the top of the mesh
        ! If we consider imaginary meshes, we can compute this in a same way with below.
        ! If T(i+1,j+1) is not in the boundary, we can find the fluxes using (3.3) in PDF file
        ELSE
           FluxSouth=FluxNorth
           FluxNorth=FluxS(i,j+1)
        END IF
        
        ! FluxVert  = (G(i+1,j+1.5) - G(i+1,j+0.5))/area
        FluxVert(:,i,j)=(FluxNorth-FluxSouth)/ac(i,j)
        !IF (j .EQ. 2) THEN
			!PRINT *, i,j,FluxNorth,FluxSouth,FluxVert(:,i,j)
			!STOP
        !END IF
        !WRITE (*,*) i,j,FluxNorth(1), FluxSouth(1), FluxVert(1,i,j)
     END DO
  END DO

END SUBROUTINE dfdp

  !==========

!	Goal : To calculate G(i,j-0.5)
!	Input : i,j - index of T or q
!	Output : FluxS(i,j) = G(i,j-0.5)	
!  	Parameters 
!		hhori : vector whose size is Nx*(Np+1). the length of x-direction of the meshes.
!		vcheck : vector whose size is Nx*(Np+1). 
!					Defined by (3.4) in PDF file. In the code, we calculate this in central_upwind.f90							
!		A : matrix size (2,Nx+2,Np+2). This parameter means T^n, q^n. We use this for Runge-Kutta method.

FUNCTION FluxS(i,j)
	IMPLICIT NONE

	INTEGER, INTENT(in) :: i,j
	DOUBLE PRECISION, DIMENSION(3) :: FluxS
	DOUBLE PRECISION :: ucheck, wcheck, vcheck

	ucheck =  (rS(i,j,1)*A(3,i+1,j+1)+rS(i,j,2)*A(3,i+1,j))/(rS(i,j,1)+rS(i,j,2))
	wcheck = (rS(i,j,1)*w(i+1,j+1)+rS(i,j,2)*w(i+1,j))/(rS(i,j,1)+rS(i,j,2))
	
!	ucheck = (A(3,i+1,j+1)+A(3,i+1,j))/2.d0
!	wcheck = (w(i+1,j+1)+w(i+1,j))/2.d0
	
	
	vcheck = normal(1,(j-1)*Nx+i)*ucheck+normal(2,(j-1)*Nx+i)*wcheck
	
	! use upwind or downward
	! Index A is followed by index of T
	IF (j .EQ. 1) THEN
		FluxS=0.d0
!	ELSE IF (j .EQ. Np+1) THEN
!		FluxS=0.d0
	ELSE IF (vcheck .GE. 0) THEN
		FluxS=hhori((j-1)*Nx+i)*vcheck*A(:,i+1,j)
	ELSE
		FluxS=hhori((j-1)*Nx+i)*vcheck*A(:,i+1,j+1)
	END IF
	
	!WRITE(22,*) "i,j=",i,j,r,ucheck,wcheck,normal(1,(j-1)*Nx+i),normal(2,(j-1)*Nx+i),hhori((j-1)*Nx+i),vcheck
	!WRITE(22,*) "i,j=",i,j,xm3(i+1,j),xm(i+1,j),xm(i+1,j+1),xm(i+1,j),r,w(i+1,j+1),w(i+1,j),wcheck
	
	
	!PRINT *, hhori((j-1)*Nx+i), vcheck(i,j), A(:,i,j-1)

END FUNCTION  

  !==========

!	Goal : Find the FluxHori (East and West flux defined by (3.9) in PDF file.
!	Parameters 
!			FluxHori is (2,Nx,Np) matrix. 
!			FluxHori(1,i,j) = Flux for T(i+1,j+1), FluxVert(2,i,j) = Flux for q(i+1,j+1)
!			FluxNorth : North flux i.e. F(i+0.5,j+1), double vector size 2 for T and q
!   		FluxSouth : South flux i.e. F(i+1.5,j+1), double vector size 2 for T and q
!			ac(i,j) : area for (i,j) mesh.
!			Nx, Np : the number of elements of x and p, respectively.
!	More or less, follow the similar method with the previous North-South case.
SUBROUTINE dgdx()
  IMPLICIT NONE

  INTEGER :: i,j
  DOUBLE PRECISION :: FluxEast(3), FluxWest(3)

  !print *, 'upwind reconstruction'
  ! Loop preceed from left to right in terms of the cells. i.e. i-direction
  Do j=1,Np
     Do i=1,Nx
		! call FluxW function which means F(i,j) = F(i-0.5,j)
		! On the left part of boundary, we do
        IF ( i .EQ. 1) THEN
           FluxWest=FluxW(i,j)
           FluxEast=FluxW(i+1,j)
        ELSE
           FluxWest=FluxEast
           FluxEast=FluxW(i+1,j)
        END IF
        FluxHori(:,i,j)=(FluxEast-FluxWest)/ac(i,j)
     END DO
  END DO

END SUBROUTINE dgdx

  !==========


!	Goal : To calculate F(i-0.5,j)
!	Input : i,j - index of T or q
!	Output : FluxS(i,j) = F(i-0.5,j)	
!  	Parameters 
!		A : matrix size (2,Nx+2,Np+2). This parameter means T^n, q^n. We use this for Runge-Kutta method.
FUNCTION FluxW(i,j)

	IMPLICIT NONE
	
	INTEGER, INTENT(in) :: i,j
	DOUBLE PRECISION, DIMENSION(3) :: FluxW
	DOUBLE PRECISION, DIMENSION(3) :: Aw
	DOUBLE PRECISION :: ucheck
	
	ucheck =  (rW(i,j,1)*A(3,i+1,j+1)+rW(i,j,2)*A(3,i,j+1))/(rW(i,j,1)+rW(i,j,2))
!	ucheck =  (0.5d0*A(3,i+1,j+1)+0.5d0*A(3,i,j+1))
	! Here, we use upwind & downward scheme.
	IF ( ucheck .GE. 0) THEN
		Aw=A(:,i,j+1)
	ELSE
		Aw=A(:,i+1,j+1)
	END IF
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Former eq. : FluxW = dp(i)*u(i,j)*Aw
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FluxW=dp(i)*ucheck*Aw
  
END FUNCTION
 
  !==========
