CC = gcc
AR = ar
RL = ranlib
all :
	cd LsdSrc ; make "CC = ${CC}" "AR = ${AR}" "RL = ${RL}" ; cp lsd ..
	cd OutlsdSrc ; make "CC = ${CC}" ; cp outlsd ..
	cd GenposSrc ; make "CC = ${CC}" ; cp genpos ..
	cd Mol2abSrc ; make "CC = ${CC}" ; cp mol2ab ..

clean :
	cd LsdSrc ; make clean
	cd OutlsdSrc ; make clean
	cd GenposSrc ; make clean
	cd Mol2abSrc ; make clean
	rm -f lsd outlsd genpos mol2ab

cleaner :
	make clean
	rm -f Makefile m_edit solve *~

test :
	cp Data/pinene .
	./solve pinene
