CC=gcc
AR=ar
RL=ranlib
OUT=lsd
OBJS=lsd.o lecture.o start.o corr.o doubles.o display.o sub.o equiv.o heap.o remcorr.o chklsd.o distance.o dyntext.o lvm.o connexe.o isom.o
LIB_DIRS=-L.
LIB_FILES=-luniciter -linchizer -lm

$(OUT) : $(OBJS) libuniciter.a libinchizer.a
	$(CC) $(OBJS) -o $(OUT) $(LIB_DIRS) $(LIB_FILES)

libuniciter.a :
	(cd Uniciter ; make "CC=$(CC)" "AR=$(AR)" "RL=$(RL)" ; cp libuniciter.a ..)

libinchizer.a :
	(cd Inchizer ; make "CC=$(CC)" "AR=$(AR)" "RL=$(RL)" ; cp libinchizer.a ..)

$(OBJS) : defs.h

.c.o : $<
	$(CC) -c $<

clean :
	(cd Uniciter ; make clean)
	(cd Inchizer ; make clean)
	rm -f *.o *.a lsd lsd.exe *~ *.sol solncounter* stoplsd*

