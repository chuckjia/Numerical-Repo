include ../Makefile.inc

DATA = atm_data.f90
TIME_METHOD = fluxes.f90 time_method.f90
CENTRAL_UPWIND_FLUXES = central_upwind.f90

Files90= $(DATA) $(CENT_FLUXES) $(CENTRAL_UPWIND_FLUXES) $(TIME_METHOD)

libatm:
	$(F2PY) $(F2PYFLAGS) -m $@ -c $(Files90)
	mv $@.so ../

fluxes:
	$(F90) -c $(Files90)
	
CLEAN : 
	rm -rf *.so *.mod *.o *~ ._*
