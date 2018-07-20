  !==========

  ! delta(w,q) compute the function gamma
  ! which equals H(-w)H(q-qs)
  ! 
  ! INPUT :
  !	w : double
  !	q : double
  !
  ! OUTPUT :
  !	delta : double
  !
  FUNCTION delta(i,j)

    IMPLICIT NONE

    DOUBLE PRECISION :: delta
    INTEGER, INTENT(in) :: i,j
    DOUBLE PRECISION :: H1, H2

    ! f2py  INTENT(in) :: i, j
    ! f2py  INTENT(out) :: delta

    !IF ( -w(i,j-1)-w(i,j) .GE. 0.d0  ) THEN
    IF ( -w(i,j) .GE. 0.d0  ) THEN
       H1=1
    ELSE
       H1=0
    END IF

    !IF ( (A(2,i,j)-qs(i,j)) .GE. 0.d0 ) THEN
    IF ( (A(2,i,j)-qs(i,j)) .GE. 0.d0 ) THEN
       H2=1
    ELSE
       H2=0
    END IF

    delta=H1*H2
!	IF (delta .GT. 0.d0) THEN
!		write(*,*) i+j
!	END IF
  END FUNCTION delta

  !==========

  ! FUNCTION L(t) compute the latent heat evaopration at time t
  !
  ! INPUT :
  !       i : integer
  ! OUTPUT :
  !       L : double
  FUNCTION L(i,j)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: i,j
    DOUBLE PRECISION :: L

    ! f2py INTENT(in) :: i
    ! f2py INTENT(out) :: L

!    L=2.5008d0*10.d0**(6.d0)-2.3d0*10.d0**(3.d0)*(A(1,i,j)-275.d0)
	!L=0.d0
	
	!! Simplificaton of L
	L=2.5008d0*10.d0**(6.d0)
	
  END FUNCTION L

  !==========

  ! FUNCTION qs(pt,t)
  ! Compute the value of qs at the pressure pt and time t, where
  ! qs is the the specific humidity of saturation
  !
  ! INPUT :
  !    pt : double, pressure
  !    t : double, the time (global value)
  ! OUPUT :
  !    qs : double
  !
  FUNCTION qs(i,j)

    IMPLICIT NONE

    DOUBLE PRECISION :: qs
    INTEGER, INTENT(in) :: i,j

    ! f2py INTENT(in) :: i,j
    ! f2py INTENT(out) :: qs

!	! With 1/pm include
!    !qs=0.622d0*(1/pm(j))*10.d0**( (0.7859d0+0.03477d0*A(1,i,j))/(1.d0+0.00412d0*A(1,i,j)) )
    
!    !Without 1/pm
!	!qs=0.622d0**10.d0**( (0.7859d0+0.03477d0*A(1,i,j))/(1.d0+0.00412d0*A(1,i,j)) )
	
!	! Experimental qs
!	qs=((0.622d0*6.112d0*exp( (17.67d0*(A(1,i,j)-273.15d0))/(A(1,i,j)-29.65d0) ))/(pm(i,j)))!*(1.d0+rqs(i,j))
	!qs_temp(i,j) = qs
	
	!!! Theoretical qs
	qs=0.622d0/pm(i,j)*2.53*10.d0**(8.d0)*exp( -5.43*10**3/A(1,i,j) )
	
!	write(*,*) '1',rqs(i,j)
!	write(*,*) '2',qs
!	!!!!!
!	!qs=1.d0

	END FUNCTION qs

  !==========

  FUNCTION F(i,j)

    IMPLICIT NONE

    DOUBLE PRECISION :: F
    INTEGER, INTENT (in) :: i,j
    DOUBLE PRECISION :: Lt, At, qst

    ! f2py INTENT(int) :: i,j
    ! f2py INTENT(out) :: F

    Lt=L(i,j)
    At=A(1,i,j)
    qst=qs(i,j)

    F = qst*At*( (Lt*R-Cp*Rv*At) / (Cp*Rv*At**2.d0+qst*Lt**2.d0) )
!	write(*,*) 'F',F
!	write(*,*) 'qst',qst
!	write(*,*) 'remains', At*( (Lt*R-Cp*Rv*At) / (Cp*Rv*At**2.d0+qst*Lt**2.d0) )
  END FUNCTION F

  !==========
