include Makefile.inc

IMG = IMG/
LIB = lib

# To use this "make all" or "make again"
# creat the fortran lib
# all : lib
all: $(LIB)

again: CLEAN all

lib :
	(cd l_fluxes ; make libatm)

SCP : BIGCLEAN
	scp -r ../$(REP) se-mt7.math.indiana.edu:~/


BIGCLEAN : CLEAN
	rm -rf *.so

CLEAN : run
	rm -rf *.pyc $(EXEC) *.mod *.o *~ img/*png img/*.avi DEBUG/*.log ._* *.log DATA/*.data .DS* ./*/.DS* *.png
	rm -rf ./run_*/divergence_*.png ./run_*/lda_*.png ./run_*/ldax_*.png ./run_*/quiver_*.png ./run_*/*.dat  ./run_*/*.mat ./run_*/*.png
	(cd l_fluxes ; make CLEAN)

run :
	rm -rf run_*/IMG* run_*/DATA run_*/DEBUG run_*/*.pyc run_*/*~ run_*/*.log run_*/fort*
	rm -rf read*/IMG*

VIDEO :
	@echo Creating a video from the images in $(IMG)...
	(cd $(IMG) ; mencoder mf://*.png -mf w=800:h=600:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi)
	mv $(IMG)/output.avi ./
	
SAVE : BIGCLEAN	
	cp -rf ../$(REP) /home/arthbous/Dropbox/Atmosphere/

RESTORE : BIGCLEAN
	cp -prf /home/arthbous/Dropbox/ATM/$(REP) ../
	
TAR : BIGCLEAN
	rm -rf ../$(TAR_NAME).tar.gz
	(cd .. ; tar -cvzf $(TAR_NAME).tar.gz $(REP)/ )
