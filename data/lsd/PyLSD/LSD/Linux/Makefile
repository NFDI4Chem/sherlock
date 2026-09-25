all :
	cd LsdSrc ; make ; cp lsd ..
	cd OutlsdSrc ; make ; cp outlsd ..
	cd GenposSrc ; make ; cp genpos ..
	cd Mol2abSrc ; make ; cp mol2ab ..

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
