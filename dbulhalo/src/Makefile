OBJS = main.o prof.o ran.o gasdev.o tree.o treef.o kernel.o dmhaloprof.o extf.o diskdeng.o vgrid.o disksgen.o setgrid.o seteos.o hydstat.o diskmetg.o ndenggrid.o expdgen.o gasdgen.o vgasgrid.o para.o setcool.o cool.o dmhaloprof2.o dmhdengen.o dmhvgrid.o dmhpgen.o hpbulext.o buldengen.o bulvgrid.o bulpgen.o exp2dmr.o hgext.o hgdengen.o hgvgrid.o hgpgen.o hgsetang.o hgtemprof.o
EXEC = pgen
#OPTIONS = -O -frecord-marker=4
OPTIONS = -g -mcmodel=medium
#MLIBS = -lm -mmacosx-version-min=10.7
### for MPI
# kepler 
MPI_DIR=/cluster/mpich2_3.0
# mac
MPI_DIR = /usr/local
#
FORT = ${MPI_DIR}/bin/mpif90
MPI_INC = ${MPI_DIR}/include
MPI_LIB = ${MPI_DIR}/lib 
FFLAGS = ${OPTIONS} -I${MPI_INC}
LIBS = ${MLIBS} -L${MPI_LIB} ${OPTIONS} 
### cosmos
#FORT=ifort
#CPP=cpp
#OPTIONS=-g -O3 -xHost -align -ansi-alias -mcmodel=medium -traceback
#FFLAGS = ${OPTIONS}
#LIBS=-lmpi
.SUFFIXES : .o .f

${EXEC} : ${OBJS} define.f
	${FORT} -o ${EXEC} ${OBJS} ${LIBS}
	\mv ${EXEC} ..

.f.o : define.f
	${FORT} ${FFLAGS} -c $<

clean:
	rm *.o
