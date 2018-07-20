MODULE central_upwind

  ! Central Upwind fluxes

  USE atm_data

  !==========

CONTAINS

	
	INCLUDE "central_upwind_1.f90"

	INCLUDE "central_upwind_2.f90"
	
	INCLUDE "central_upwind_3.f90"
	
	INCLUDE "central_MUSCL.f90"

END MODULE central_upwind
