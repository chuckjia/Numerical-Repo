MODULE time_method

  USE fluxes
  USE atm_data ! module with data for the method

CONTAINS

  SUBROUTINE euler1()

    IMPLICIT NONE

    INTEGER :: i, j, l
    DOUBLE PRECISION, DIMENSION(3,Nx,Np) :: A0
    
    
    CALL Res1()
    A0=A(:,2:Nx+1,2:Np+1)
    A(:,2:Nx+1,2:Np+1)= A0 + dt*(DA + S(:,2:Nx+1,2:Np+1))
    CALL Borders1()
    CALL PROJECTION2(dt)
    CALL FOMEGA()
	CALL Borders3()

    
  
  END SUBROUTINE euler1

  !==========

  ! rk4source(methods)
  ! Compute the Runge Kutta 4 method on all the mesh for the shallow water system.
  !
  ! dA/dt = Res(A) + S(t)
  !
  ! With :
  ! 	A = transpose(Tq,u) (data)
  ! 	S = transpose(ST,Sq,Su) (source terme)
  !
  ! /!\ this function uses the data in the module upwindfluxes_data /!\
  !
  ! INPUT :
  !	method = 0 center fluxes, 1 central-upwind fluxes
  !
  SUBROUTINE rk4source()

    IMPLICIT NONE

    INTEGER :: i, j, l
    DOUBLE PRECISION, DIMENSION(3,Nx,Np) :: A0, k1, k2, k3

    !PRINT *, 'TOTO

    !print*,'ATTENTION : f=0'

    !PRINT *,FluxHori(1,1,:)
	
    ! Computation of k1
    ! T, q
    CALL Res1()
    A0=A(:,2:Nx+1,2:Np+1)
    k1=DA + S(:,2:Nx+1,2:Np+1)
    A(:,2:Nx+1,2:Np+1)= A0 + dt*k1/2.d0
    CALL Borders1()
	! omega
    CALL FOMEGA()
	CALL Borders3()
    

	! Computation of k2
    CALL Res1()
    k2=DA + Sdt2(:,2:Nx+1,2:Np+1)
    A(:,2:Nx+1,2:Np+1)= A0 + dt*k2/2.d0
    CALL Borders1()
	! omega
	CALL FOMEGA()
	CALL Borders3()

    ! Computation of k3
    CALL Res1()
    k3= DA + Sdt2(:,2:Nx+1,2:Np+1)
    A(:,2:Nx+1,2:Np+1) = A0 + dt*k3
    CALL Borders1()
	! omega
    CALL FOMEGA()
	CALL Borders3()
	
    ! Computation of k4
	! Results
    CALL Res1()
    A(:,2:Nx+1,2:Np+1) = A0 + (dt/6.d0)*( k1 + 2.d0*k2 + 2.d0*k3 + DA + Sdt(:,2:Nx+1,2:Np+1))
    CALL Borders1()
	! omega
    CALL FOMEGA()
	CALL Borders3()

    !PRINT *,DZ(1,:)
    !PRINT *,Z(1,:)

  END SUBROUTINE rk4source

  !==========

  SUBROUTINE rk4_2()

    IMPLICIT NONE

    INTEGER :: i,j
    DOUBLE PRECISION, DIMENSION(3,Nx,Np) :: A0, k1, k2, k3

    ! Computation of k1
    CALL Res1()
    A0=A(:,2:Nx+1,2:Np+1)
    k1=DA + S(:,2:Nx+1,2:Np+1)
    A(:,2:Nx+1,2:Np+1)= A0 + dt*k1/2.d0
    CALL Borders1()
    CALL PROJECTION2(dt/2.d0)
    CALL FOMEGA()
	CALL Borders3()

	! Computation of k2
	CALL Res1()
    k2=DA + Sdt2(:,2:Nx+1,2:Np+1)
    A(:,2:Nx+1,2:Np+1)= A0 + dt*k2/2.d0
    CALL Borders1()
    CALL PROJECTION2(dt/2.d0)
	CALL FOMEGA()
	CALL Borders3()

	
    ! Computation of k3
    CALL Res1()
    k3= DA + Sdt2(:,2:Nx+1,2:Np+1)
    A(:,2:Nx+1,2:Np+1) = A0 + dt*k3
    CALL Borders1()
    CALL PROJECTION2(dt)
    CALL FOMEGA()
	CALL Borders3()
	
	
    ! Computation of k4
	! Results
    CALL Res1()
    A(:,2:Nx+1,2:Np+1) = A0 + (dt/6.d0)*( k1 + 2.d0*k2 + 2.d0*k3 + DA + Sdt(:,2:Nx+1,2:Np+1))
    CALL Borders1()
    CALL PROJECTION2(dt)
    CALL FOMEGA()
	CALL Borders3()


  END SUBROUTINE rk4_2

  !==========
  
  SUBROUTINE rk4_3()

    IMPLICIT NONE

    INTEGER :: i,j
    DOUBLE PRECISION, DIMENSION(3,Nx,Np) :: A0, k1, k2, k3

    ! Computation of k1
    CALL Res1()
    A0=A(:,2:Nx+1,2:Np+1)
    k1=DA + S(:,2:Nx+1,2:Np+1)
    A(:,2:Nx+1,2:Np+1)= A0 + dt*k1/2.d0
    CALL Borders1()
 


	! Computation of k2
	CALL Res1()
    k2=DA + Sdt2(:,2:Nx+1,2:Np+1)
    A(:,2:Nx+1,2:Np+1)= A0 + dt*k2/2.d0
    CALL Borders1()


	
    ! Computation of k3
    CALL Res1()
    k3= DA + Sdt2(:,2:Nx+1,2:Np+1)
    A(:,2:Nx+1,2:Np+1) = A0 + dt*k3
    CALL Borders1()

	
    ! Computation of k4
	! Results
    CALL Res1()
    A(:,2:Nx+1,2:Np+1) = A0 + (dt/6.d0)*( k1 + 2.d0*k2 + 2.d0*k3 + DA + Sdt(:,2:Nx+1,2:Np+1))
    CALL Borders1()



  END SUBROUTINE rk4_3

  !==========
  
END MODULE time_method
