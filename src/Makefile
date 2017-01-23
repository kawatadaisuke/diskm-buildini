FORT = gfortran
OBJS = main.o rot.o ran.o gasdev.o
EXEC = inigen
MLIBS = -lm
OPTIONS = -g -frecord-marker=4
FFLAGS = ${OPTIONS} 
LIBS = ${MLIBS} ${OPTIONS}
.SUFFIXES : .o .f

${EXEC} : ${OBJS}
	${FORT} -o ${EXEC} ${OBJS} ${LIBS}
	\mv ${EXEC} ..

.f.o :
	${FORT} ${FFLAGS} -c $<

clean:
	rm *.o
