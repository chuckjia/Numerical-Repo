VERSION = 2.1
TAR_NAME = atm2d-v$(VERSION)
REP = atm2d-v$(VERSION)

F90= gfortran
F90FLAGS = -g
FC = gfortran
FCFLAGS = -O3 -pg
CC = gcc
CFLAGS = -Wall
LDFLAGS = -lm -Wall
F2PY = f2py2.7
F2PYFLAGS =  --fcompiler=$(F90) #-DF2PY_REPORT_ON_ARRAY_COPY=0
#/usr/lib/libg2c.so.0
