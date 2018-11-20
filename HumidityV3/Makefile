run:
	@echo "*********************************************************"
	@echo "* Warning: You are about to compile with GCC.           *"
	@echo "* It is advisable to compile with clang for efficiency. *"
	@echo "*********************************************************"
	@echo ">> To compile with clang, use command: make clang. "
	@echo "   If you choose to use the above command, then Clang-6.0"
	@echo "   will be used for compilation. -Ofast will be chosen as"
	@echo "   the optimization level."
	@echo ">> Press Y/y to continue compile with GCC. Press any "
	@echo "   other keys to stop."
	@read -p ">> " key; \
	if [ "$$key" = "y" ] ; \
	then echo ""; echo "[[Compiling with GCC]]" ; echo ""; \
	g++ Driver.cpp -o Output/out -std=c++11 -Ofast ; \
	Output/./out notInEclipse ; \
	else echo "Operation aborted."; \
	fi 

clang:
	clang Driver.cpp -o Output/out -lm -lstdc++ -std=c++11 -Ofast
	Output/./out notInEclipse

cleanall:
	rm -rf Output MovieFrames
	mkdir Output MovieFrames

all:
	run
