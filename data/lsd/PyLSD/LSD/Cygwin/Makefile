all :
	cd LsdSrc ; make ; cp lsd.exe ..
	cd OutlsdSrc ; make ; cp outlsd.exe ..
	cd GenposSrc ; make ; cp genpos.exe ..
	cd Mol2abSrc ; make ; cp mol2ab.exe ..

clean :
	cd LsdSrc ; make clean
	cd OutlsdSrc ; make clean
	cd GenposSrc ; make clean
	cd Mol2abSrc ; make clean
	rm -f lsd.exe outlsd.exe genpos.exe mol2ab.exe *~

cleaner :
	make clean
	rm -f Makefile m_edit solve *~

test :
	cp Data/pinene .
	./solve pinene

