MODULE time_method

  USE fluxes
  USE atm_data ! module with data for the method

CONTAINS

  ! rk2source(methods)
  ! Compute the Runge Kutta 2 method on all the mesh for the shallow water system.
  !
  ! dQ/dt = Res(Q) + S(t)
  !
  ! With :
  ! 	Q = transpose(h,U,V) (data)
  ! 	S = transpose(Sh,SU,SV) (source terme)
  !
  ! /!\ this function uses the data in the module upwindfluxes_data /!\
  !
  ! INPUT :
  !	method = 0 center fluxes, 1 central-upwind fluxes
  !
  SUBROUTINE rk2source()

    IMPLICIT NONE

    INTEGER :: i,j
    DOUBLE PRECISION, DIMENSION(2,Nx,Np) :: A0, k1
    DOUBLE PRECISION, DIMENSION(Nx,Np) :: u_ini, ku1


    IF (DEBUG==1) THEN
       write(snt,"(i3.3)") nt
       open(unit=12,file='DEBUG/rk2source_'//trim(snt)//'.log',status="replace")
       WRITE (12,*) "###############################"
       WRITE (12,*) "### RK2Source ..."
       WRITE (12,*) "###############################"
       WRITE (12,*) "Nx, Np =",Nx, Np
       WRITE (12,*) "###############################"
       WRITE (12,*) "### K1 ..."
       WRITE (12,*) "###############################"
       write(snt,"(i3.3)") nt
       snt=trim(snt)//'-K1'
    END IF

    IF (DEBUG==1) THEN
       write(snt,"(i3.3)") nt
       open(unit=13,file='DEBUG/rk2sourceU_'//trim(snt)//'.log',status="replace")
       WRITE (13,*) "###############################"
       WRITE (13,*) "### RK2Source ..."
       WRITE (13,*) "###############################"
       WRITE (13,*) "Nx, Np =",Nx, Np
       WRITE (13,*) "###############################"
       WRITE (13,*) "### Ku1 ..."
       WRITE (13,*) "###############################"
       write(snt,"(i3.3)") nt
       snt=trim(snt)//'-K1'
       CLOSE(13)
    END IF


    !print*,"ATTENTION : f=0"

    !PRINT *,FluxHori(1,1,:)

    CALL Res()

    ! Computation of k1
    DO i=1,Nx
       DO j=1,Np
          A0(1,i,j)=A(1,i,j)
          A0(2,i,j)=A(2,i,j)
          k1(1,i,j)=DA(1,i,j) + S(1,i,j)
          k1(2,i,j)=DA(2,i,j) + S(2,i,j)
          A(1,i,j) = A0(1,i,j) + dt*k1(1,i,j)
          A(2,i,j) = A0(2,i,j) + dt*k1(2,i,j)
          u_ini(i,j)=u(i,j)
          ku1(i,j)=DU(i,j) + S2(i,j)
          u(i,j) = u_ini(i,j) + dt*ku1(i,j)
       END DO
    END DO
    DO j=1,Np
       u(Nx+1,j)=u(1,j)
    END DO

    !CALL FOMEGA()

    ! Computation of k2
    IF (DEBUG==1) THEN
       WRITE (12,*) "-----------------"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "DA(1,i,j)=",DA(1,i,j)
             WRITE (12,*) "DA(2,i,j)=",DA(2,i,j)
          END DO
       END DO
       WRITE (12,*) "-----------------"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "A(1,i,j)=",A(1,i,j)
             WRITE (12,*) "A(2,i,j)=",A(2,i,j)
          END DO
       END DO
       write(snt,"(i3.3)") nt
       snt=trim(snt)//'-K2'
       WRITE (12,*) "###############################"
       WRITE (12,*) "### End K1"
       WRITE (12,*) "###############################"
       WRITE (12,*) "###############################"
       WRITE (12,*) "### K2 ..."
       WRITE (12,*) "###############################"
    END IF


    CALL Res

    IF (DEBUG==1) THEN
       WRITE (12,*) "-----------------"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "DA(1,i,j)=",DA(1,i,j)
             WRITE (12,*) "DA(2,i,j)=",DA(2,i,j)
          END DO
       END DO
       WRITE (12,*) "-----------------"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "A(1,i,j)=",A(1,i,j)
             WRITE (12,*) "A(2,i,j)=",A(2,i,j)
          END DO
       END DO
       WRITE (12,*) "###############################"
       WRITE (12,*) "### End K2"
       WRITE (12,*) "###############################"
       WRITE (12,*) "###############################"
       WRITE (12,*) "### Results"
       WRITE (12,*) "###############################"
    END IF

    ! Results
    DO i=1,Nx
       DO j=1,Np
          A(1,i,j) = A0(1,i,j) + (dt/2.d0)*( k1(1,i,j) + DA(1,i,j) + Sdt(1,i,j))
          A(2,i,j) = A0(2,i,j) + (dt/2.d0)*( k1(2,i,j) + DA(2,i,j) + Sdt(2,i,j))
          u(i,j) = u_ini(i,j) + (dt/2.d0)*( ku1(i,j) + Du(i,j) + S2dt(i,j))
       END DO
    END DO
    DO j=1,Np
       u(Nx+1,j)=u(1,j)
    END DO

    !CALL FOMEGA()

    IF (DEBUG==1) THEN
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "A(1,i,j)=",A(1,i,j)
             WRITE (12,*) "A(2,i,j)=",A(2,i,j)
          END DO
       END DO
       WRITE (12,*) "###############################"
       WRITE (12,*) "### End RK2Source"
       WRITE (12,*) "###############################"
       CLOSE(12)
    END IF

    !PRINT *,DQ(1,2,:)
    !PRINT *,Q(1,2,:)

  END SUBROUTINE rk2source

  !==========

  ! rk4source(methods)
  ! Compute the Runge Kutta 4 method on all the mesh for the shallow water system.
  !
  ! dQ/dt = Res(Q) + S(t)
  !
  ! With :
  ! 	Q = transpose(h,U,V) (data)
  ! 	S = transpose(Sh,SU,SV) (source terme)
  !
  ! /!\ this function uses the data in the module upwindfluxes_data /!\
  !
  ! INPUT :
  !	method = 0 center fluxes, 1 central-upwind fluxes
  !
  SUBROUTINE rk4source()

    IMPLICIT NONE

    INTEGER :: i, j, l
    DOUBLE PRECISION, DIMENSION(2,Nx,Np) :: A0, k1, k2, k3
    DOUBLE PRECISION, DIMENSION(Nx,Np) :: u_ini, ku1, ku2, ku3


    IF (DEBUG==1) THEN
       write(snt,"(i3.3)") nt
       open(unit=12,file='DEBUG/rk4source_'//trim(snt)//'.log',status="replace")
       WRITE (12,*) "###############################"
       WRITE (12,*) "### RK4Source ..."
       WRITE (12,*) "###############################"
       WRITE (12,*) "Nx, Np =",Nx, Np
       WRITE (12,*) "###############################"
       WRITE (12,*) "### K1 ..."
       WRITE (12,*) "###############################"
       write(snt,"(i3.3)") nt
       snt=trim(snt)//'-K1'
    END IF

    !print*,"ATTENTION : f=0"

    !PRINT *,FluxHori(1,1,:)

    ! Computation of k1
    CALL Res()

    A0=A
    k1=DA + S
    A= A0 + dt*k1
    u_ini=u
    ku1= DU + S2
    u(1:Nx,1:Np) = u_ini(1:Nx,1:Np) + dt*ku1(1:Nx,1:Np)
    u(Nx+1,1:Np)=u(1,1:Np)

    !CALL FOMEGA()

    ! Computation of k2
    IF (DEBUG==1) THEN
       WRITE (12,*) "-----------------"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "DA(1,i,j)=",DA(1,i,j)
             WRITE (12,*) "DA(2,i,j)=",DA(2,i,j)
          END DO
       END DO
       WRITE (12,*) "-----------------"
       !DO k=1,N(1)
       !WRITE (12,*) "i=",i
       !WRITE (12,*) "Z(1,i)=",Z(1,i)
       !WRITE (12,*) "Z(2,i)=",Z(2,i)
       !END DO
       WRITE (12,*) "-----------------"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "k1(1,i,j)=",k1(1,i,j)
             WRITE (12,*) "k1(2,i,j)=",k1(2,i,j)
          END DO
       END DO
       write(snt,"(i3.3)") nt
       snt=trim(snt)//'-K2'
       WRITE (12,*) "###############################"
       WRITE (12,*) "### End K1"
       WRITE (12,*) "###############################"
       WRITE (12,*) "###############################"
       WRITE (12,*) "### K2 ..."
       WRITE (12,*) "###############################"
    END IF

    CALL Res()

    DO i=1,Nx
       DO j=1,Np
          k2(1,i,j)=DA(1,i,j) + Sdt2(1,i,j)
          k2(2,i,j)=DA(2,i,j) + Sdt2(2,i,j)
          A(1,i,j) = A0(1,i,j) + dt*k2(1,i,j)/2.d0
          A(2,i,j) = A0(2,i,j) + dt*k2(2,i,j)/2.d0
          ku2(i,j)=DU(i,j) + S2dt2(i,j)
          u(i,j) = u_ini(1,i) + dt*ku2(i,j)/2.d0
       END DO
    END DO
    DO j=1,Np
       u(Nx+1,j)=u(1,j)
    END DO

    !CALL FOMEGA()


    ! Computation of k3
    IF (DEBUG==1) THEN
       WRITE (12,*) "-----------------"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "DA(1,i,j)=",DA(1,i,j)
             WRITE (12,*) "DA(2,i,j)=",DA(2,i,j)
          END DO
       END DO
       WRITE (12,*) "-----------------"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "k2(1,i,j)=",k2(1,i,j)
             WRITE (12,*) "k2(2,i,j)=",k2(2,i,j)
          END DO
       END DO
       write(snt,"(i3.3)") nt
       snt=trim(snt)//'-K3'
       WRITE (12,*) "###############################"
       WRITE (12,*) "### End K2"
       WRITE (12,*) "###############################"
       WRITE (12,*) "###############################"
       WRITE (12,*) "### K3 ..."
       WRITE (12,*) "###############################"
    END IF


    CALL Res()

    DO i=1,Nx
       DO j=1,Np
          k3(1,i,j)= DA(1,i,j) + Sdt2(1,i,j)
          k3(2,i,j)= DA(2,i,j) + Sdt2(2,i,j)
          A(1,i,j) = A0(1,i,j) + dt*k3(1,i,j)
          A(2,i,j) = A0(2,i,j) + dt*k3(2,i,j)
          ku3(i,j)= Du(i,j) + S2dt2(i,j)
          u(i,j) = u_ini(1,i) + dt*ku3(i,j)
       END DO
    END DO
    DO j=1,Np
       u(Nx+1,j)=u(1,j)
    END DO

    !CALL FOMEGA()

    ! Computation of k4
    IF (DEBUG==1) THEN
       WRITE (12,*) "-----------------"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "DA(1,i,j)=",DA(1,i,j)
             WRITE (12,*) "DA(2,i,j)=",DA(2,i,j)
          END DO
       END DO
       WRITE (12,*) "-----------------"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "k3(1,i,j)=",k3(1,i,j)
             WRITE (12,*) "k3(2,i,j)=",k3(2,i,j)
          END DO
       END DO
       write(snt,"(i3.3)") nt
       snt=trim(snt)//'-K4'
       WRITE (12,*) "###############################"
       WRITE (12,*) "### End K3"
       WRITE (12,*) "###############################"
       WRITE (12,*) "###############################"
       WRITE (12,*) "### K4 ..."
       WRITE (12,*) "###############################"
    END IF

    CALL Res()

    ! Results
    DO i=1,Nx
       DO j=1,Np
          A(1,i,j) = A0(1,i,j) + (dt/6.d0)*( k1( 1,i,j) &
               + 2.d0*k2( 1,i,j) + 2.d0*k3( 1,i,j) + DA( 1,i,j) + Sdt( 1,i,j) )
          A(2,i,j) = A0(2,i,j) + (dt/6.d0)*( k1( 2,i,j) &
               + 2.d0*k2( 2,i,j) + 2.d0*k3( 2,i,j) + DA( 2,i,j) + Sdt( 2,i,j) )
          u(i,j) = u_ini(i,j) + (dt/6.d0)*( ku1(i,j) &
               + 2.d0*ku2(i,j) + 2.d0*ku3(i,j) + Du(i,j) + S2dt(i,j) )
       END DO
    END DO
    DO j=1,Np
       u(Nx+1,j)=u(1,j)
    END DO

    !CALL FOMEGA()

    IF (DEBUG==1) THEN
       WRITE (12,*) "-----------------"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "DA(1,i,j)=",DA(1,i,j)
             WRITE (12,*) "DA(2,i,j)=",DA(2,i,j)
          END DO
       END DO
       WRITE (12,*) "-----------------"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "k4(1,i,j)=",DA(1,i,j) + Sdt( 1,i,j)
             WRITE (12,*) "k4(2,i,j)=",DA(2,i,j) + Sdt( 2,i,j)
          END DO
       END DO
       WRITE (12,*) "###############################"
       WRITE (12,*) "### End K4"
       WRITE (12,*) "###############################"
       DO i=1,Nx
          WRITE (12,*) "i=",i
          DO j=1,Np
             WRITE (12,*) "j=",j
             WRITE (12,*) "A(1,i,j)=",A(1,i,j)
             WRITE (12,*) "A(2,i,j)=",A(2,i,j)
          END DO
       END DO
       WRITE (12,*) "###############################"
       WRITE (12,*) "### End RK4Source"
       WRITE (12,*) "###############################"
       CLOSE(12)
    END IF

    !PRINT *,DZ(1,:)
    !PRINT *,Z(1,:)

  END SUBROUTINE rk4source

  !==========

  
END MODULE time_method
