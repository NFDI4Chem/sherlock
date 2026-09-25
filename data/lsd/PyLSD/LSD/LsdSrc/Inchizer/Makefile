CC=gcc
AR=ar
RL=ranlib

INCHI_LIB_OBJS = ichican2.o       ichicano.o \
                 ichicans.o       ichiisot.o \
                 ichilnct.o       ichimak2.o \
                 ichimake.o       ichimap1.o \
                 ichimap2.o       ichimap4.o \
                 ichinorm.o       ichiparm.o \
                 ichiprt1.o       ichiprt2.o \
                 ichiprt3.o       ichiqueu.o \
                 ichiring.o       ichisort.o \
                 ichister.o       ichitaut.o \
                 ichi_bns.o       inchi_dll.o \
                 ichiread.o       ichirvr1.o \
                 ichirvr2.o       ichirvr3.o \
                 ichirvr4.o       ichirvr5.o \
                 ichirvr6.o       ichirvr7.o \
                 inchi_dll_main.o runichi.o  \
                 strutil.o        util.o
INCHI_PROG_OBJS = e_inchi_atom.o e_util.o
INCHIZER_OBJ = inchizer.o
OBJS = $(INCHI_LIB_OBJS) $(INCHI_PROG_OBJS) $(INCHIZER_OBJ)

libinchizer.a : $(OBJS)
	rm -f libinchizer.a
	$(AR) qv libinchizer.a $(OBJS)
	$(RL) libinchizer.a

.c.o :
	$(CC) -c $<

demo : demo.o libinchizer.a
	$(CC) demo.o -o demo -L. -linchizer -lm

test : demo
	./demo

clean :
	rm -f *.o *.a *~ demo
