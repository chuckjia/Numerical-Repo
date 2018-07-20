! Central Upwind fluxes of order 2 for T,q,u

  !==========



!	Goal : Find the FluxVert (North and South flux defined by (3.10) in PDF file.
!	Parameters 
!			FluxVert is (2,Nx,Np) matrix. 
!			FluxVert(1,i,j) = Flux for T(i+1,j+1)=Tij, FluxVert(2,i,j) = Flux for q(i+1,j+1)=qij
!			FluxNorth : North flux i.e. G(i,j+0.5), double vector size 2 for T and q
!   		FluxSouth : South flux i.e. G(i,j-0.5), double vector size 2 for T and q
!			ac(i,j) : area for (i,j) mesh.
!			Nx, Np : the number of elements of x and p, respectively.
SUBROUTINE dfdp_CM()

  IMPLICIT NONE

  INTEGER :: i,j
  DOUBLE PRECISION :: FluxNorth(3), FluxSouth(3)

  Do i=1,Nx
     Do j=1,Np
		! Top and bottom parts of mesh are different.
		
		! T(i+1,2) is related to the bottom of the mesh
		! FluxSouth is zero because w is zero at the bottom
        IF (j .EQ. 1) THEN
           FluxNorth=FluxS_CM(i,j+1)
           FluxSouth=0.d0
        
        ! T(i+1,j+1) is related to the top of the mesh
        ! If we consider imaginary meshes, we can compute this in a same way with below.
        ! If T(i+1,j+1) is not in the boundary, we can find the fluxes using (3.3) in PDF file
        ELSE
           FluxSouth=FluxNorth
           FluxNorth=FluxS_CM(i,j+1)
        END IF
        
        ! FluxVert  = (G(i+1,j+1.5) - G(i+1,j+0.5))/area
        FluxVert(:,i,j)=(FluxNorth-FluxSouth)/ac(i,j)
        !IF (j .EQ. 2) THEN
			!PRINT *, i,j,FluxNorth,FluxSouth,FluxVert(:,i,j)
			!STOP
        !END IF
        !WRITE (7,*) i,j,FluxNorth, FluxSouth, FluxVert(:,i,j)
     END DO
  END DO

END SUBROUTINE dfdp_CM

  !==========

!	Goal : To calculate G(i,j-0.5)
!	Input : i,j - index of T or q
!	Output : FluxS_CM(i,j) = G(i,j-0.5)	
!  	Parameters 
!		hhori : vector whose size is Nx*(Np+1). the length of x-direction of the meshes.
!		vcheck : vector whose size is Nx*(Np+1). 
!					Defined by (3.4) in PDF file. In the code, we calculate this in central_upwind.f90							
!		A : matrix size (2,Nx+2,Np+2). This parameter means T^n, q^n. We use this for Runge-Kutta method.

FUNCTION FluxS_CM(i,j)
	IMPLICIT NONE

	INTEGER, INTENT(in) :: i,j
	INTEGER :: l
	DOUBLE PRECISION, DIMENSION(3) :: FluxS_CM, Aplus, Aminus,H_x,H_p
	DOUBLE PRECISION :: wplus, wminus, aaplus, aaminus, dpt, bplus,bminus
	DOUBLE PRECISION :: ldmaxplus, ldmaxminus, ldminplus, ldminminus
	
	IF (j .EQ. 1) THEN
		FluxS_CM=0.d0
	ELSE IF (j .EQ. Np+1) THEN
		FluxS_CM=hhori(Np*Nx+i)*normal(2,Np*Nx+i)*w(i+1,Np+2)*A(:,i+1,Np+2)&
				+ hhori(Np*Nx+i)*normal(1,Np*Nx+i)*A(3,i+1,Np+2)*A(:,i+1,Np+2)
	ELSE		
		dpt=(dp(i+1)+dp(i))/2.d0
		
		IF (j .EQ. 2) THEN
			DO l=1,3
				Aplus(l) = A(l,i+1,j  )+minmod_CM((A(l,i+1,j  )-A(l,i+1,j-1))/dpt, &
						                      (A(l,i+1,j+1)-A(l,i+1,j  ))/(dpt/2.d0))*dpt/2.d0
				Aminus(l)= A(l,i+1,j+1)-minmod_CM((A(l,i+1,j+1)-A(l,i+1,j  ))/dpt, &
						                      (A(l,i+1,j+2)-A(l,i+1,j+1))/dpt)*dpt/2.d0
			END DO
				wplus  =w(i+1,j  )+minmod_CM((w(i+1,j  )-w(i+1,j-1))/dpt, &
							             (w(i+1,j+1)-w(i+1,j  ))/(dpt/2.d0))*dpt
				wminus=w(i+1,j+1)-minmod_CM((w(i+1,j+1)-w(i+1,j  ))/dpt, &
							            (w(i+1,j+2)-w(i+1,j+1))/(dpt))*dpt/2.d0
		ELSE IF (j .EQ. Np) THEN
			DO l=1,3
				Aplus(l) = A(l,i+1,j  )+minmod_CM((A(l,i+1,j  )-A(l,i+1,j-1))/dpt, &
						                      (A(l,i+1,j+1)-A(l,i+1,j  ))/dpt)*dpt/2.d0
				Aminus(l)= A(l,i+1,j+1)-minmod_CM((A(l,i+1,j+1)-A(l,i+1,j  ))/(dpt/2.d0), &
						                      (A(l,i+1,j+2)-A(l,i+1,j+1))/dpt)*dpt
			END DO
		
				wplus  =w(i+1,j  )+minmod_CM((w(i+1,j  )-w(i+1,j-1))/dpt, &
						                 (w(i+1,j+1)-w(i+1,j  ))/(dpt))*dpt/2.d0
				wminus=w(i+1,j+1)-minmod_CM((w(i+1,j+1)-w(i+1,j  ))/(dpt/2.d0), &
						                (w(i+1,j+2)-w(i+1,j+1))/(dpt))*dpt
		ELSE
			DO l=1,3
				Aplus(l) = A(l,i+1,j  )+minmod_CM((A(l,i+1,j  )-A(l,i+1,j-1))/dpt, &
						                      (A(l,i+1,j+1)-A(l,i+1,j  ))/dpt)*dpt/2.d0
				Aminus(l)= A(l,i+1,j+1)-minmod_CM((A(l,i+1,j+1)-A(l,i+1,j  ))/dpt, &
						                      (A(l,i+1,j+2)-A(l,i+1,j+1))/dpt)*dpt/2.d0
			END DO
		
				wplus  =w(i+1,j  )+minmod_CM((w(i+1,j  )-w(i+1,j-1))/dpt, &
						                 (w(i+1,j+1)-w(i+1,j  ))/(dpt))*dpt/2.d0
				wminus=w(i+1,j+1)-minmod_CM((w(i+1,j+1)-w(i+1,j  ))/dpt, &
						                (w(i+1,j+2)-w(i+1,j+1))/(dpt))*dpt/2.d0
		END IF
		
		aaplus=2.d0*max(ABS(Aplus(3)),ABS(Aminus(3)))
		
		bplus=max(ABS(wplus),ABS(wminus))
		

		H_x(:) = hhori((j-1)*Nx+i)*normal(1,(j-1)*Nx+i)&
			*(Aplus(3)*Aplus(:)+Aminus(3)*Aminus(:) - aaplus*(Aminus(:)-Aplus(:)))/2.d0

		H_p(:) = hhori((j-1)*Nx+i)*normal(2,(j-1)*Nx+i)&
			*(wplus*Aplus(:)+wminus*Aminus(:)- bplus*(Aminus(:)-Aplus(:)))/2.d0

		
		FluxS_CM(:)=H_x(:) + H_p(:)
		
		!FluxS_CM(:)= hhori((j-1)*Nx+i)*normal(1,(j-1)*Nx+i)*( (aaplus*Aplus(3)*Aplus(:)-aaminus*Aminus(3)*Aminus(:)
		!		+ aaplus*aaminus*(Aplus(:)-Aminus(:)))/(aaplus-aaminus) ) 
		!		+ hhori((j-1)*Nx+i)*normal(2,(j-1)*Nx+i)*( (bplus*wplus*Aplus(:)-bminus*wminus*Aminus(:)
		!		+ bplus*bminus*(Aplus(:)-Aminus(:)))/(bplus-bminus) ) 
							
!		IF ((i .EQ. 3) .and. (j .EQ. 2)) THEN
!			write(10,*) 'omega',w(i+1,:)
!			write(9,*) 'omega',w(i,:)
!			write(*,*) 'omega(4,2), omega(4,3)',w(i+1,j),w(i+1,j+1)
!			write(*,*) 'aaplus, aaminus, bplus, bminus',aaplus,aaminus,bplus,bminus
!			write(*,*) 'up,um,wp,wm',Aplus(3),Aminus(3),wplus,wminus
!		END IF
	END IF
	


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
SUBROUTINE dgdx_CM()
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
           FluxWest=FluxW_CM(i,j)
           FluxEast=FluxW_CM(i+1,j)
        ELSE
           FluxWest=FluxEast
           FluxEast=FluxW_CM(i+1,j)
        END IF
        FluxHori(:,i,j)=(FluxEast-FluxWest)/ac(i,j)
        !WRITE (8,*) i,j,FluxWest, FluxEast, FluxHori(:,i,j)
     END DO
  END DO

END SUBROUTINE dgdx_CM

  !==========


!	Goal : To calculate F(i-0.5,j)
!	Input : i,j - index of T or q
!	Output : FluxS_CM(i,j) = F(i-0.5,j)	
!  	Parameters 
!		A : matrix size (2,Nx+2,Np+2). This parameter means T^n, q^n. We use this for Runge-Kutta method.
FUNCTION FluxW_CM(i,j)

	IMPLICIT NONE
	
	INTEGER, INTENT(in) :: i,j
	INTEGER :: l
	DOUBLE PRECISION, DIMENSION(3) :: FluxW_CM, Aplus, Aminus
	DOUBLE PRECISION, DIMENSION(3) :: AreconTop, AreconBot
	DOUBLE PRECISION, DIMENSION(3) :: Gradx1, Gradx2, Gradx3, Gradp1, Gradp2, Gradp3
	DOUBLE PRECISION :: aaminus, aaplus
	

	IF (i .EQ. 1) THEN
		FluxW_CM=dp(i)*A(3,i,j+1)*A(:,i,j+1)
	ELSE IF (i .EQ. Nx+1) THEN
		FluxW_CM=dp(Nx+1)*A(3,Nx+1,j+1)*A(:,Nx+1,j+1)
	ELSE
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Computation of Grad1
		AreconBot=coeff3(i-1,j  ,1)*A(:,i-1,j)+coeff3(i-1,j  ,2)*A(:,i,j  )+&
		coeff3(i-1,j  ,3)*A(:,i,j+1)+coeff3(i-1,j  ,4)*A(:,i-1,j+1)
		
		AreconTop=coeff3(i-1,j+1,1)*A(:,i-1,j+1)+coeff3(i-1,j+1,2)*A(:,i,j+1)+&
		coeff3(i-1,j+1,3)*A(:,i,j+2)+coeff3(i-1,j+1,4)*A(:,i-1,j+2)
		
		Gradx1=m3(i-1,j,1)*(A(:,i,j+1)-A(:,i-1,j+1))+m3(i-1,j,2)*(AreconTop-AreconBot)
		Gradp1=m3(i-1,j,3)*(A(:,i,j+1)-A(:,i-1,j+1))+m3(i-1,j,4)*(AreconTop-AreconBot)

		! Computation of Grad3
		AreconBot=coeff3(i,j  ,1)*A(:,i,j)+coeff3(i,j  ,2)*A(:,i+1,j  )+&
		coeff3(i,j  ,3)*A(:,i+1,j+1)+coeff3(i,j  ,4)*A(:,i,j+1)
		AreconTop=coeff3(i,j+1,1)*A(:,i,j+1)+coeff3(i,j+1,2)*A(:,i+1,j+1)+&
		coeff3(i,j+1,3)*A(:,i+1,j+2)+coeff3(i,j+1,4)*A(:,i,j+2)
		
		Gradx3=m3(i,j,1)*(A(:,i+1,j+1)-A(:,i,j+1))+m3(i,j,2)*(AreconTop-AreconBot)
		Gradp3=m3(i,j,3)*(A(:,i+1,j+1)-A(:,i,j+1))+m3(i,j,4)*(AreconTop-AreconBot)
		
		! Computation of Grad3
		!Gradx2=m2(i-1,j,1)*(A(:,i+1,j+1)-A(:,i-1,j+1))+m2(i-1,j,2)*(A(:,i,j+2)-A(:,i,j))
		!Gradp2=m2(i-1,j,3)*(A(:,i+1,j+1)-A(:,i-1,j+1))+m2(i-1,j,4)*(A(:,i,j+2)-A(:,i,j))
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		DO l=1,3
			Aplus(l) = A(l,i  ,j+1)+minmod_CM(Gradx1(l),Gradx3(l))*dx/2.d0  &
						+minmod_CM(Gradp1(l),Gradp3(l))*(pm2(i,j+1)-pm(i  ,j+1))
		END DO
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Computation of Grad1
		Gradx1=Gradx3
		Gradp1=Gradp3

		! Computation of Grad3
		AreconBot=coeff3(i+1,j  ,1)*A(:,i+1,j)+coeff3(i+1,j  ,2)*A(:,i+2,j  )+&
		coeff3(i+1,j  ,3)*A(:,i+2,j+1)+coeff3(i+1,j  ,4)*A(:,i+1,j+1)
		AreconTop=coeff3(i+1,j+1,1)*A(:,i+1,j+1)+coeff3(i+1,j+1,2)*A(:,i+2,j+1)+&
		coeff3(i+1,j+1,3)*A(:,i+2,j+2)+coeff3(i+1,j+1,4)*A(:,i+1,j+2)
		
		Gradx3=m3(i+1,j,1)*(A(:,i+2,j+1)-A(:,i+1,j+1))+m3(i+1,j,2)*(AreconTop-AreconBot)
		Gradp3=m3(i+1,j,3)*(A(:,i+2,j+1)-A(:,i+1,j+1))+m3(i+1,j,4)*(AreconTop-AreconBot)
		
		! Computation of Grad3
		!Gradx2=m2(i,j,1)*(A(:,i+2,j+1)-A(:,i,j+1))+m2(i,j,2)*(A(:,i+1,j+2)-A(:,i+1,j))
		!Gradp2=m2(i,j,3)*(A(:,i+2,j+1)-A(:,i,j+1))+m2(i,j,4)*(A(:,i+1,j+2)-A(:,i+1,j))
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		DO l=1,3
			Aminus(l)= A(l,i+1,j+1)-minmod_CM(Gradx1(l),Gradx3(l))*dx/2.d0 &
						+minmod_CM(Gradp1(l),Gradp3(l))*(pm2(i,j+1)-pm(i+1,j+1))
		END DO
		

		
		aaplus=2.d0*max(ABS(Aplus(3)),ABS(Aminus(3)))
		
		

		FluxW_CM(:)=dp(i)*(Aplus(3)*Aplus(:)+Aminus(3)*Aminus(:) &
			- aaplus*(Aminus(:)-Aplus(:)))/2.d0
		
	END IF
	
!	WRITE(7,*)"W",i,j,FluxW
  
END FUNCTION

!==========

! minmod_CM(a,b,c)
!
FUNCTION minmod_CM(a,b)

	IMPLICIT NONE
	
	DOUBLE PRECISION :: minmod_CM
	DOUBLE PRECISION :: a,b
	DOUBLE PRECISION :: minarg, maxarg
	
	minarg = min(a,b)
	maxarg = max(a,b)
	IF (minarg > 0) THEN
	  minmod_CM = minarg
	ELSEIF (maxarg < 0) THEN
	  minmod_CM = maxarg
	ELSE
	  minmod_CM = 0.d0
	END IF

END FUNCTION

!==========

