CC=gcc
OUT=genpos
${OUT} : genpos.o 
	${CC} genpos.o -lm -o ${OUT}

genpos.o : genpos.c defs.h
	${CC} -c genpos.c

clean :
	rm -f *.o genpos genpos.exe *~
