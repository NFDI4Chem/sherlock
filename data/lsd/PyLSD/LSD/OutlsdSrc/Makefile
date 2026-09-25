CC=gcc
OUT=outlsd
OBJS1=outlsd.o
OBJS2=smiles.o jacobi.o matutil.o
OBJS=${OBJS1} ${OBJS2}

${OUT} : ${OBJS}
	${CC} ${OBJS} -lm -o ${OUT}

outlsd.o : outlsd.c jacobi.h matutil.h defs.h 
	${CC} -c outlsd.c

smiles.o : smiles.c 
	${CC} -c smiles.c

jacobi.o : jacobi.c jacobi.h
	${CC} -c jacobi.c

matutil.o : matutil.c matutil.h
	${CC} -c matutil.c

clean :
	rm -f *.o outlsd outlsd.exe *~

