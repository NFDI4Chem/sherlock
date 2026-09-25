CC=gcc
OUT=mol2ab
${OUT} : mol2ab.o 
	${CC} mol2ab.o -lm -o ${OUT}

genpos.o : mol2ab.c
	${CC} -c mol2ab.c

clean :
	rm -f *.o mol2ab mol2ab.exe *~
