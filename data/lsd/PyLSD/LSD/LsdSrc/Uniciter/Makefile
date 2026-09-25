CC=gcc
AR=ar
RL=ranlib

libuniciter.a : uniciter.o avl.o
	rm -f libuniciter.a
	$(AR) qv libuniciter.a uniciter.o avl.o
	$(RL) libuniciter.a
	
demo : libuniciter.a demo.o
	$(CC) demo.o -o demo -L. -luniciter

.c.o :
	$(CC) -c $<

test : demo
	demo ccc bbb bbb aaa ccc

clean :
	rm -f *.o *.a *~ demo demo.exe
