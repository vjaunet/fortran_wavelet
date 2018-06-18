#nom de l'executable
EXEC = cwvttr

#construction des modules
MODULES=$(shell grep -l 'module' *.f90)
MOD=$(MODULES:.f90=.mod)

#source : rech de .f90
SRC=$(wildcard *.f90)

#construction des objets
OBJ=$(SRC:.f90=.o)

#-----------------------------------------------------------
#compilation :

IPATH=-I/usr/include -I/data2/jaunetv/UTILITIES/fortran_lib/MOD
LPATH=-L/usr/lib/x86_64-linux-gnu -L/data2/jaunetv/UTILITIES/fortran_lib/LIB

CC= gfortran
CC+= -O2

LIBS= $(LPATH) -lfftw3 -lpressdata -lspectral

#option de vectorisation et parallelisation:
OPT_para = #
OpenMp  =  #

#option de debugage :
Debug= -Wall -fcheck=all -g

ALL:$(MOD) $(EXEC)

debug: CC += $(Debug)
debug: ALL

parallel: CC += $(OPT_para) $(OpenMp)
parallel: ALL

$(EXEC):$(OBJ)
	$(CC) -o $@ $^ $(LIBS)

%.mod:%.f90
	$(CC) $(IPATH) -c $^ $(LIBS)

%.o: %.f90
	$(CC) $(IPATH) -o $@ -c $^ $(LIBS)

clean:
	rm -rf *.o *~ *.mod
	rm -rf $(EXEC)
