MODULE fluxes

  USE central_upwind
  USE atm_data

CONTAINS

  !==========

  ! Res1(method) compute DQ where :
  ! DZ(1) =  dF(Q)/dx - Phi(Q)
  ! 
  !
  ! Shallow Water equations with Coriolis force :
  ! dQ/dt = Res(Q) + S(t)
  !
  ! With :
  ! Q = transpose(h,U,V) (data)
  !	F = transpose(Fh,Fu,Fv)
  !	Phi = transpose(Phih,Phiu,Phiv) (coriolis force which equals to f)
  !	S = transpose(Sh,SU,SV) (source terme)
  !
  ! INPUT :
  !	method = 0 center fluxes, 1 central-upwind fluxes
  !
  
  SUBROUTINE Res1()

    IMPLICIT NONE

    INTEGER :: i, j

    ! Fluxes
    !CALL Fvcheck()
    IF (scheme .EQ. 1) THEN
		CALL dfdp()
		CALL dgdx()
	ELSE IF (scheme .EQ. 2) THEN
		CALL dfdp2()
		CALL dgdx2()
	ELSE !IF (scheme .EQ. 3) THEN
		CALL dfdp3()
		CALL dgdx3()
!	ELSE 
!		CALL dfdp4()
!		CALL dgdx4()
	END IF
	CALL FPHI()
	CALL FPHIX()
	

	! DA(1,:,:) and DA(2,:,:) means -1*second term of (2.2)_1 and (2.2)_2, respectively
	! Here, index of DA followed by the index of FluxVert and FluxHori
    DA= - FluxVert - FluxHori
    
    !write(*,*) SIZE(DA)
    !write(7,*) 'vert',FluxVert
    !write(8,*) 'hori',FluxHori
!    DO i=1,Nx
!		write (9,*) 'i is ',i,FluxVert(1,i,:)
!	end do
!	DO i=1,Nx
!		write (10,*) 'i is ',i,FluxHori(2,i,:)
!	end do
	
	! If PHYSIC == 0
	!		R.H.S of (2.2)_1 and (2.2)_2 are zero.	
	! If PHYSIC == 1
	!		we apply the equations on the R.H.S. of (2.2)_1 and (2.2)_2.
	!PRINT *, "PHYSIC=", physic
	IF (PHYSIC .EQ. 1) THEN	
		DO i=1,Nx
		   DO j=1,Np
		      ! p0 or pB?
			  DA(1,i,j) = DA(1,i,j) + & 
			  w(i+1,j+1)/pm(i+1,j+1)*(R*A(1,i+1,j+1)-delta(i+1,j+1)*L(i+1,j+1)*F(i+1,j+1))/Cp
				
			  DA(2,i,j) = DA(2,i,j) + delta(i+1,j+1)*F(i+1,j+1)*w(i+1,j+1)/pm(i+1,j+1)
			
			  DA(3,i,j) = DA(3,i,j) - phix(i,j)
			  
			  !write(*,*) '1',w(i+1,j+1)/pm(i+1,j+1)*(R*A(1,i+1,j+1)-delta(i+1,j+1)*L(i+1,j+1)*F(i+1,j+1))/Cp
			  
		   END DO
		END DO
!	ELSEIF (PHYSIC .EQ. 2) THEN
!		  DA(3,:,:) = DA(3,:,:) - phix(:,:)
	END IF
	! Here, delta, F, L is defined in atm_functions.f90
	END SUBROUTINE Res1

  !==========

  ! Here, we apply the boundary values for Tij, qij.
  ! We simply use the periodic conditions.
  ! Parameter :
  ! 		A - A(1,i,j) is T(i+1,j+1)
  !			A - A(2,i,j) is q(i+1,j+1)
  
  ! The order is important. We should apply top and bottom and then apply
  ! left and right boundary.
  SUBROUTINE Borders1()
  
  IMPLICIT NONE
  
  INTEGER :: i,j
    
  DOUBLE PRECISION, DIMENSION(3) :: T0
  DOUBLE PRECISION	::OP1,OP2,x_bar,y_bar,M1,M2
  ! First, apply the top and boundary nodes.
  ! This follows the above nodes (bottom boundary case)  
  ! and the right below nodes (top boundary caes)
  ! For example, A(2,1) = A(2,2) => known value, the above nodes.
  !				 A(2,Np+2) = A(2,Np+1) => known value , the below nodes
!  DO i=1,Nx
!	A(:,i+1,1   )=A(:,i+1,2   )
!	A(:,i+1,Np+2)=A(:,i+1,Np+1)
!  END DO
  DO i=1,Nx
	A(:,i+1,1   )=A(:,i+1,2   )
	!A(:,i+1,Np+2)=A(:,i+1,Np+1)
!	IF (  abs( p(i,Np+1) - p(i+1,Np+1) ) < 1E-8) THEN
	A(:,i+1,Np+2)=A(:,i+1,Np+1)
!	ELSE
!	! x_bar, y_bar are the coordinate of the "Neumann" imaginary point
!	! perpendicular to the mountain.
!	! And then we project.
!		T0 = A(:,i+1,Np+1)
		
!		M1 = (pm2(i+1,Np+2)-pm2(i,Np+2))/(xm2(i+1,Np+2)-xm2(i,Np+2))
!		M2 = (xm2(i,Np+2)-xm2(i+1,Np+2))/(pm2(i+1,Np+2)-pm2(i,Np+2))
!		x_bar = (M1*xm2(i,Np+2)-M2*xm(i+1,Np+1)+pm(i+1,Np+1)-pm2(i,Np+2))/(M1-M2)
!		y_bar = (M1*(x_bar-xm2(i,Np+2)))+pm2(i,Np+2)
		
!		OP2 = sqrt((xm(i+1,Np+1)-x_bar)**2 + (pm(i+1,Np+1)-y_bar)**2)
!		!OP1 = (dp(i)+dp(i+1))/4
!		OP1 = sqrt((xm(i+1,Np+2)-xm(i+1,Np+1))**2 + (pm(i+1,Np+2)-pm(i+1,Np+1))**2)
!		A(:,i+1,Np+2) = T0*OP1/OP2
		
!		!write(*,*) 'num1',OP1/OP2
!	END IF

  END DO

  
  ! Now, apply the west and east boundary nodes
	IF (bc .EQ. 0) THEN
	! Periodic case
	  DO j=1,Np+2
		! The left boundary nodes follows the second last nodes from the right boundary
		! Here, we consider the imaginary meshes and periodicity.
		A(:,1   ,j)=A(:,Nx+1,j)
		! Similarly, the right boundary nodes follow the second nodes from the left boundary
		! Here, again, we consider the imaginary meshes and periodicity
		A(:,Nx+2,j)=A(:,2   ,j)
	  END DO
	ELSE IF (bc .EQ. 1) THEN
	! Neumann case
	  DO j=1,Np+2
		A(:,1   ,j)=A(:,2   ,j)
		A(:,Nx+2,j)=A(:,Nx+1,j)
	  END DO
	ELSE
		PRINT *, "Wrong bc"
	END IF
  
  END SUBROUTINE Borders1
  
  !==========

  ! Here, we apply the boundary values for Tij, qij.
  ! We simply use the periodic conditions.
  ! Parameter :
  ! 		A(3) - A(3,i+1,j+1) is u_{i,j}
  !			phi - phi(i+1,j+1) is phi_{i,j}
  ! The order is important. We should apply top and bottom and then apply
  ! left and right boundary.
  SUBROUTINE Borders2()
  
  IMPLICIT NONE
  
  INTEGER :: i,j
  
  ! First, apply the top and boundary nodes.
  ! This follows the above nodes (bottom boundary case)  
  ! and the right below nodes (top boundary caes)
  DO i=1,Nx
	A(3,i+1,1   )=A(3,i+1,2   )
	A(3,i+1,Np+2)=A(3,i+1,Np+1)
	phi(i+1,1   )=phi(i+1,2   )
	phi(i+1,Np+2)=phi(i+1,Np+1)
  END DO
  
  ! Now, apply the left and right boundary nodes
  IF (bc .EQ. 0) THEN
	  DO j=1,Np+2
		! The left boundary nodes follows the second last nodes from the right boundary
		! Here, we consider the imaginary meshes and periodicity.
		A(3,1   ,j)=A(3,Nx+1,j)
		phi(1   ,j)=phi(Nx+1,j)
		! Similarly, the right boundary nodes follow the second nodes from the left boundary
		! Here, again, we consider the imaginary meshes and periodicity
		A(3,Nx+2,j)=A(3,2   ,j)
		phi(Nx+2,j)=phi(2   ,j)
	END DO
  ELSE IF (bc .EQ. 1) THEN
  ! Neumann
	DO j=1,Np+2
		A(3,1   ,j)=A(3,2   ,j)
		phi(1   ,j)=phi(2   ,j)
		A(3,Nx+2,j)=A(3,Nx+1,j)
		phi(Nx+2,j)=phi(Nx+1,j)
	END DO
  ELSE
	PRINT *, "Wrong bc"
  END IF
  
  END SUBROUTINE Borders2
  
  !==========
  ! Here, we apply the boundary values for wij.
  ! We simply use the periodic conditions.
  SUBROUTINE Borders3()
  
  IMPLICIT NONE
  INTEGER :: i,j
  DOUBLE PRECISION	::T0,OP1,OP2,x_bar,y_bar,M1,M2 

!  DO i=1,Nx
!	w(i+1,Np+2)=w(i+1,Np+1)
!  END DO
  
  
  ! New one
  DO i=1,Nx
	w(i+1, 1) = 0.d0
	w(i+1,Np+2)=w(i+1,Np+1)
  END DO
 
  IF (bc .EQ. 0) THEN
	DO j=2,Np+2
		w(1,j) = w(Nx+1,j)
		w(Nx+2,j) = w(2,j) 
	END DO
  ELSE IF (bc .EQ. 1) THEN
	DO j=2,Np+2
		w(1   ,j) = w(2   ,j)
		w(Nx+2,j) = w(Nx+1,j) 
	END DO
  ELSE
	PRINT *, "Wrong bc"
  END IF
  
  END SUBROUTINE Borders3
  
    !==========
  
  SUBROUTINE FPHI()
  
  IMPLICIT NONE

  INTEGER :: i,j

!!! WAY 1
!   DO i=1,Nx
!	  phi(i+1,2)=lda(i)
!	  DO j=2,Np
!		  phi(i+1,j+1) = phi(i+1,j)-R*A(1,i+1,j+1)/p(i+1,j+1)*(dp(i)+dp(i+1))/2
!	  END DO
!   END DO

!!! WAY 2
   DO i=1,Nx
	  !phi(i+1,Np+2)= zm(i+1,Np+2)*g
	  phi(i+1,Np+1)= zm(i+1,Np+2)*g - R*A(1,i+1,Np+1)/pm(i+1,Np+1)*(pm(i+1,Np+1)-pm(i+1,Np+2))
	  DO j=Np-1,1,-1
		  phi(i+1,j+1) = phi(i+1,j+2) - R*(A(1,i+1,j+1)+A(1,i+1,j))/2/pm(i+1,j+1)*(pm(i+1,j+1)-pm(i+1,j+2))
	  END DO
  END DO		
  
  DO i=1,Nx
	phi(i+1,1   )=phi(i+1,2   )
	phi(i+1,Np+2)=phi(i+1,Np+1)
  END DO

! Now, apply the west and east boundary nodes
	IF (bc .EQ. 0) THEN
	! Periodic case
	  DO j=1,Np+2
		phi(1   ,j)=phi(Nx+1,j)
		phi(Nx+2,j)=phi(2   ,j)
	  END DO
	ELSE IF (bc .EQ. 1) THEN
	! Neumann case
	  DO j=1,Np+2
		phi(1   ,j)=phi(2   ,j)
		phi(Nx+2,j)=phi(Nx+1,j)
	  END DO
	ELSE
		PRINT *, "Wrong bc"
	END IF
  
  END SUBROUTINE FPHI
  
  !==========
  
  SUBROUTINE FPHIX()
  
  IMPLICIT NONE
  ! dim(phix) = Nx*Np
  ! dim(phi) = (Nx+2)*(Np+2)
  DOUBLE PRECISION, DIMENSION(Nx+1,Np-1) :: pu, cT
  DOUBLE PRECISION, DIMENSION(Nx,Np-1) :: dTx
  DOUBLE PRECISION :: Gradx1, Gradx2, Gradx3
  DOUBLE PRECISION :: ma,mb,mc,md,AC1,AC2
  INTEGER :: i,j
  DOUBLE PRECISION, DIMENSION(Nx,Np) :: phicopy
  
  
	 DO i=1,Nx+1
		DO j=1,Np-1
			cT(i,j) = coeff(i,j,1)*A(1,i,j+1)+coeff(i,j,2)*A(1,i+1,j+1)&
			+coeff(i,j,3)*A(1,i,j+2)+coeff(i,j,4)*A(1,i+1,j+2)
		END DO
	END DO

	phix(:,:) = 0.d0
	DO i=1,Nx
		phix(i,1)=0.d0
		DO j=1,Np-1
			dTx(i,j)=m(i,j,1)*(cT(i+1,j)-cT(i,j)) + m(i,j,2)&
			*(A(1,i+1,j+2)-A(1,i+1,j+1))
			
			phix(i,j+1)=phix(i,j)-R*dTx(i,j)/pm3(i+1,j+1)*(pm(i+1,j+2)-pm(i+1,j+1))
		END DO
	END DO
  
 
 !!! WAY 3
  
  
  END SUBROUTINE FPHIX
  
    !==========
  
  SUBROUTINE PROJECTION(s)
  
	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(in ) :: s
	DOUBLE PRECISION :: s1,s2,AreconBot,AreconTop
	DOUBLE PRECISION, DIMENSION(Nx-1) :: Gradx
	DOUBLE PRECISION, DIMENSION(Nx) :: lambda
	DOUBLE PRECISION, DIMENSION(Nx) :: ld2
	DOUBLE PRECISION, DIMENSION(5,Nx) :: LUt
	INTEGER :: i,j
	
	!!!! WAY 1
		DO i=1,Nx
		s1=0.d0
		s2=0.d0
		IF (i .EQ. Nx) THEN
			DO j=1,Np
				s2=s2+A(3,i+2,j+1)*(dp(i+1)+dp(i+1))/2.d0
				s1=s1+A(3,i  ,j+1)*(dp(i-1)+dp(i  ))/2.d0
			END DO
		ELSE IF (i .EQ. 1) THEN
			DO j=1,Np
				s2=s2+A(3,i+2,j+1)*(dp(i+1)+dp(i+2))/2.d0
				s1=s1+A(3,i  ,j+1)*(dp(i  )+dp(i  ))/2.d0
			END DO
		ELSE
			DO j=1,Np
				s2=s2+A(3,i+2,j+1)*(dp(i+1)+dp(i+2))/2.d0
				s1=s1+A(3,i  ,j+1)*(dp(i-1)+dp(i  ))/2.d0
			END DO
		END IF
	
		ld(i)=(s2-s1)/(2*dx)/(pm(i+1,Np+2)-pm(i+1,1))
	END DO
	
	ld(Nx)=0.d0

	LUt=LU(Nx,s)
	
	!! Solving U*ld2=ld
	ld2(1)=ld(1)/LUt(1,1)
	DO i=2,Nx-1
		ld2(i)=(ld(i)-ld2(i-1)*LUt(2,i))/LUt(1,i)
	END DO
	ld2(Nx)=ld(Nx)-ld2(Nx-1)*LUt(2,Nx)
	
	DO i=1,Nx-2
		ld2(Nx)=ld2(Nx)-ld2(i)*LUt(3,i)
	END DO
	ld2(Nx)=ld2(Nx)/LUt(1,Nx)

	
	!! Solving L*ld=ld2
	lda(Nx)=ld2(Nx)
	lda(Nx-1)=(ld2(Nx-1)-lda(Nx)*LUt(4,Nx-1))
	DO i=Nx-2,1,-1
		lda(i)=(ld2(i)-lda(i+1)*LUt(4,i)-LUt(5,i)*lda(Nx))
	END DO
	
	
	!!!! WAY 2
	!CALL Corner_u()
	!CALL Gradu()
!	DO i=1,Nx
!		ld(i)= ( A(3,i+1,Np+2)*(p(i+1,Np+1)-p(i,Np+1))/dx &
!			+ (A(3,i+2,1)-A(3,i,1))/(2.d0*dx)*(pm(i+1,2)-pm(i+1,1)) &
!			+ dux(i,Np-1)*(pm(i+1,Np+2)-pm(i+1,Np+1)) )/(pm(i+1,Np+2)-pm(i+1,1))
!		DO j=1,Np-1
!			ld(i) = ld(i) + (dux(i,j)*(pm(i+1,j+2)-pm(i+1,j+1)))/(pm(i+1,Np+2)-pm(i+1,1))
!		END DO
!	END DO
	
	!!!! WAY 3
!	DO i=1,Nx-1
!			AreconBot=(A(3,i+2,1)-A(3,i+1,1))/2.d0
					
!			AreconTop=coeff3(i+1,2,1)*A(3,i+1,2)+coeff3(i+1,2,2)*A(3,i+2,2)+&
!				coeff3(i+1,2,3)*A(3,i+2,3)+coeff3(i+1,2,4)*A(3,i+1,3)
					
!			Gradx(i)=m3(i+1,1,1)*(A(3,i+2,2)-A(3,i+1,2))+m3(i+1,1,2)*(AreconTop-AreconBot)
!		DO j=2,Np-1
!			AreconBot=AreconTop
		
!			AreconTop=coeff3(i+1,j+1,1)*A(3,i+1,j+1)+coeff3(i+1,j+1,2)*A(3,i+2,j+1)+&
!				coeff3(i+1,j+1,3)*A(3,i+2,j+2)+coeff3(i+1,j+1,4)*A(3,i+1,j+2)
					
!			Gradx(i)=Gradx(i)+m3(i+1,j,1)*(A(3,i+2,j+1)-A(3,i+1,j+1))+m3(i+1,j,2)*(AreconTop-AreconBot)
!		END DO
!		AreconBot=AreconTop
		
!		AreconTop=(A(3,i+2,Np+2)-A(3,i+1,Np+2))/2.d0
					
!		Gradx(i)=Gradx(i)+m3(i+1,Np,1)*(A(3,i+2,Np+1)-A(3,i+1,Np+1))+m3(i+1,Np,2)*(AreconTop-AreconBot)
!	END DO
	
!	lambda(1)=0.d0
	
!	DO i=1,Nx-1
!		lambda(i+1)=(((A(3,i+2,Np+2)-A(3,i+1,Np+2))/2.d0)*(pm(i+2,Np+2)-pm(i+1,Np+2))/dx &
!				+ dp(i+1)*Gradx(i))*dx+lambda(i)
!	END DO
	
!	ldax(:) = lambda(:)/(pm(2:Nx+1,Np+2)-pm(2:Nx+1,1))
	
  
  END SUBROUTINE PROJECTION
  
  !==========
  
  !LU(1)=l1 diagonal
  !LU(2)=l2
  !LU(3)=l3
  !LU(4)=u1
  !LU(5)=u2
  
  FUNCTION LU(Nxt,s)
  
	IMPLICIT NONE
	
	INTEGER, INTENT(in) :: Nxt
	DOUBLE PRECISION, INTENT(in) :: s
	DOUBLE PRECISION, DIMENSION(Nxt) :: diag,low,up
	DOUBLE PRECISION, DIMENSION(5,Nxt) :: LU
	INTEGER :: i,j
	
	diag=-2.d0*s/(dx**2)
	low=(-(pm2(2:Nx+1,Np+2)-pm2(1:Nx,Np+2))/(2*dx**2)/(pm(2:Nx+1,Np+2)-pm(2:Nx+1,1))+1/dx**2)*s
	up=((pm2(2:Nx+1,Np+2)-pm2(1:Nx,Np+2))/(2*dx**2)/(pm(2:Nx+1,Np+2)-pm(2:Nx+1,1))+1/dx**2)*s
	
!	write(*,*) 'a',a
!	write(*,*) 'b',b
!	write(*,*) 'c',c
	
	
	IF (bc .EQ. 0) THEN
		LU(1,1)=diag(1)!-2.d0
		LU(4,1)=up(1)/LU(1,1)!1/LU(1,1)
		LU(5,1)=low(1)/LU(1,1)!1/LU(1,1)
		LU(3,1)=up(Nx)!1.d0
		!LU(2,2)=low(2)!1.d0
		
		
		 
		DO i=2,Nxt-2
			LU(2,i)=low(i)
			LU(1,i)=diag(i)-LU(2,i)*LU(4,i-1)
			LU(3,i)=-LU(3,i-1)*LU(4,i-1)
			LU(4,i)=up(i)/LU(1,i)
			LU(5,i)=-(LU(2,i)*LU(5,i-1))/LU(1,i)
		END DO
		
		! Changes lambda_N-1 = 0
		LU(2,Nxt-1)=0.d0
		LU(1,Nxt-1)=1.d0
		LU(4,Nxt-1)=0.d0
		
		LU(2,Nxt)=low(Nx)-LU(3,Nxt-2)*LU(4,Nxt-2)
		
!		LU(2,Nxt-1) = b(Nxt-1)
!		LU(1,Nxt-1) = a(Nxt-1)-LU(2,Nxt-1)*LU(4,Nxt-2)
!		LU(4,Nxt-1) = (c(Nxt-1) - LU(2,Nxt-1)*LU(5,Nxt-2) )/LU(1,Nxt-1)
!		LU(2,Nxt) = b(Nxt) - LU(3,Nxt-2)*LU(4,Nxt-2)
		
		
		LU(1,Nxt)=diag(Nx)
		DO i=1,Nxt-2
			LU(1,Nxt)=LU(1,Nxt)-LU(3,i)*LU(5,i)
		END DO
		LU(1,Nxt)=LU(1,Nxt)-LU(2,Nxt)*LU(4,Nxt-1)
		
		!write(11,*) LU(1,:)
	ELSE
		LU(3,:)=0.d0
		LU(5,:)=0.d0
		LU(1,1)=diag(1)+low(1)
		LU(4,1)=up(1)/LU(1,1)
		LU(2,1) = 0.d0
		DO i=2,Nxt-1
			LU(2,i)=low(i)
			LU(1,i)=diag(i)-LU(2,i)*LU(4,i-1)
			LU(4,i)=up(i)/LU(1,i)
		END DO
		
		LU(2,Nxt)=low(Nxt)
		LU(1,Nxt) = diag(Nxt)+up(Nxt) - LU(2,Nxt)*LU(4,Nxt-1)
		LU(4,Nxt) = 0.d0

		!!!!!
		LU(1,Nxt)=1.0
		LU(2,Nxt)=0.0
		LU(4,Nxt)=0.0
		!!!!


	END IF
	
  END FUNCTION LU
  
  !==========
  
  ! Compute MLU(1)=L and MLU(2)=U
  FUNCTION MLU(Nxt,s)
  
	IMPLICIT NONE
	INTEGER, INTENT(in) :: Nxt
	DOUBLE PRECISION, INTENT(in ) :: s
	INTEGER :: i
	DOUBLE PRECISION, DIMENSION(2,Nxt,Nxt) :: MLU
	DOUBLE PRECISION, DIMENSION (5,Nxt) :: LUt
	
	LUt=LU(Nxt,s)
	MLU=0.d0
	
	DO i=1,Nxt
		MLU(1,i,i)=LUt(1,i)
		MLU(2,i,i)=1.d0
	END DO
	DO i=2,Nxt
		MLU(1,i,i-1)=LUt(2,i)
		MLU(2,i-1,i)=LUt(4,i-1)
	END DO
	DO i=1,Nxt-2
		MLU(1,Nxt,i)=LUt(3,i)
		MLU(2,i,Nxt)=LUt(5,i)
	END DO
  
  END FUNCTION MLU  
  
   !==========
   
   SUBROUTINE PROJECTION2(s)
	
		IMPLICIT NONE
		
		DOUBLE PRECISION, INTENT(in ) :: s
		!DOUBLE PRECISION, DIMENSION(Nx) :: lambdax
		DOUBLE PRECISION :: deltaT,p0,T0,C
		INTEGER :: i,j
		

		
		CALL PROJECTION(s)
		
		!PRINT *,"lambda",lda

		DO i=1,Nx
			IF (i .EQ. 1) THEN
				IF (bc .EQ. 0) THEN
					ldax(i)=(lda(i+1)-lda(Nx))/(2.d0*dx)
				ELSE IF (bc .EQ. 1) THEN
					ldax(i)=(lda(i+1)-lda(i))/(2.d0*dx)
				END IF
			ELSE IF (i .EQ. Nx) THEN
				IF (bc .EQ. 0) THEN
					ldax(i)=(lda(1)-lda(i-1))/(2.d0*dx)
				ELSE IF (bc .EQ. 1) THEN
					ldax(i)=(lda(i)-lda(i-1))/(2.d0*dx)
				END IF
			ELSE
				ldax(i)=(lda(i+1)-lda(i-1))/(2.d0*dx)
			END IF
		END DO
		
		
!		DO i=1,Nx
!			phi(i+1,Np+1)=-R*A(1,i+1,Np+1)/pm(i+1,Np+1)*(dp(i)+dp(i+1))/4.d0+zm(i+1,Np+2)*g
!			DO j=Np-1,1,-1
!				!PRINT *, "i,j=",i,j
!				phi(i+1,j+1)=phi(i+1,j+2)-R*A(1,i+1,j+1)/pm(i+1,j+1)*(dp(i)+dp(i+1))/2.d0
!			END DO
!		END DO
!		DO i=1,Nx
!			phi(i+1,Np+1)=-R*A(1,i+1,Np+1)/pm(i+1,Np+1)*(dp(i)+dp(i+1))/2.d0+zm(i+1,Np+2)*g
!			DO j=Np-1,1,-1
!				phi(i+1,j+1)=phi(i+1,j+2)-R*A(1,i+1,j+1)/pm(i+1,j+1)*(dp(i)+dp(i+1))/2.d0
!			END DO
!		END DO		
		!PRINT *,R
!		DO i=1,Nx
!			phi(i+1,1)=-R*A(1,i+1,Np+1)/pm(i+1,Np+1)*(dp(i)+dp(i+1))/2.d0+zm(i+1,1)*g
!			DO j=2,Np-1
!				phi(i+1,j+1)=phi(i+1,j+2)-R*A(1,i+1,j+1)/pm(i+1,j+1)*(dp(i)+dp(i+1))/2.d0
!			END DO
!		END DO				

		!phi(2:Nx+1,2:Np+1)=0.d0!zm(2:Nx+1,2:Np+1)*g
	
!	deltaT=50.d0
!	p0=1000.d0
!	T0=300.d0	
!	C = R*(T0-deltaT)*log(p0) + R*deltaT
!    phi(2:Nx+1,2:Np+1)= (-R*(T0-deltaT)*log(pm(2:Nx+1,2:Np+1)) - R*deltaT/p0*pm(2:Nx+1,2:Np+1) +C)

		

		
		DO j=1,Np
			!WRITE (9,*) j,A(3,2:Nx+1,j+1)
			A(3,2:Nx+1,j+1)=A(3,2:Nx+1,j+1)-ldax(1:Nx)*s
			!WRITE (9,*) j,A(3,2:Nx+1,j+1)
			!phi(2:Nx+1,j+1)=lda(1:Nx)
			!phi(2:Nx+1,j+1)=phi(2:Nx+1,j+1)+lda(1:Nx)
		END DO
		
		CALL Borders2()
   
   END SUBROUTINE PROJECTION2
   
  !==========
	
  ! Goal : Get the w(i,j)
  SUBROUTINE FOMEGA()

    INTEGER :: i,j
    DOUBLE PRECISION :: dux1, dux2, dux3
	
	CALL Corner_u()
	CALL Gradu()
	
	! we do not apply this for left and right boundary cases.
	! Thus, the loop goes from 1 to Nx.
!    DO i=1,Nx
!		! w == 0 when it's on the bottom
!		w(i+1,1)=0.d0
		
!		!! j=1
!		dux2=(A(3,i+2,2)-A(3,i,2))/(2.d0*dx)
!		w(i+1,2)= - dux2*(pm(i+1,2)-pm(i+1,1))
		
!		!! j=2
!		!dux1=dux(i,1) 
!		w(i+1,3) = w(i+1,2) - ( 3.d0*dux(i,1) - dux2 )*(pm(i+1,3)-pm(i+1,2))/2.d0
!		w(i+1,4) = w(i+1,3) - ( 23.d0*dux(i,2) - 16.d0*dux(i,1) +5.d0*dux2)*(pm(i+1,4)-pm(i+1,3))/12.d0
!		DO j=4,Np
!			w(i+1,j+1) = w(i+1,j) - ( 23.d0*dux(i,j-1) - 16.d0*dux(i,j-2)  +5.d0*dux(i,j-3))*(pm(i+1,j+1)-pm(i+1,j))/12.d0
!			!dux2=dux1
!			!dux1=dux(i,j-1)
!		END DO
!	END DO

    DO i=1,Nx
		! w == 0 when it's on the bottom
		w(i+1,1)=0.d0
		w(i+1,2)= - (A(3,i+2,2)-A(3,i,2))/(2.d0*dx)*(pm(i+1,2)-pm(i+1,1))
		DO j=2,Np
			w(i+1,j+1) = w(i+1,j) - dux(i,j-1)*(pm(i+1,j+1)-pm(i+1,j))
		END DO
	END DO

	
	
	
  END SUBROUTINE FOMEGA

  !==========
  
  
  SUBROUTINE Corner_u()
	
    INTEGER :: i,j
	DO i=1,Nx+1
		DO j=1,Np-1
			cu(i,j) = coeff(i,j,1)*A(3,i,j+1)+coeff(i,j,2)*A(3,i+1,j+1)&
			+coeff(i,j,3)*A(3,i,j+2)+coeff(i,j,4)*A(3,i+1,j+2)
		END DO
	END DO
	
	! In vector calc form
	! cu(:,:) = coeff(:,:,1)*A(3,1:Nx+1,2:Np) + coeff(:,:,2)*A(3,2:Nx+2,2:Np) + 
	!                        coeff(:,:,3)*A(3,1:Nx+1,3:Np+1) + coeff(:,:,4)*A(3,2:Nx+2,3:Np+1)
	
  END SUBROUTINE Corner_u
  
  
  !==========
  
  
  
  
  SUBROUTINE Gradu()

    INTEGER :: i,j
    
	DO i=1,Nx
		DO j=1,Np-1
			dux(i,j) = m(i,j,1)*(cu(i+1,j)-cu(i,j)) + m(i,j,2)*(A(3,i+1,j+2)-A(3,i+1,j+1))
			!du(i,j,2) = m(i,j,3)*(cu(i+1,j)-cu(i,j)) + m(i,j,4)*(A(3,i+1,j+2)-A(3,i+1,j+1))
		END DO
	END DO
	
	! du(:,:,1) = m(:,:,1)*(cu(2:Nx+1,1:Np-1)-cu(1:Nx,1:Np-1)) + m(:,:,2)*(A(3,2:Nx+1,3:Np+1)-A(3,2:Nx+1,2:Np))
	! du(:,:,1) = m(:,:,3)*(cu(2:Nx+1,1:Np-1)-cu(1:Nx,1:Np-1)) + m(:,:,4)*(A(3,2:Nx+1,3:Np+1)-A(3,2:Nx+1,2:Np))
  END SUBROUTINE Gradu

END MODULE fluxes
