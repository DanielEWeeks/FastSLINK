CFLAGS= -O -I/usr/local/include
SFLAG=
LDFLAGS=-lm 
CC=gcc

all:		slink unknown isim msim lsim	

unknown:	unknown.c
		$(CC) $(CFLAGS) -o unknown unknown.c $(LDFLAGS)

isim:		isim.c
		$(CC) $(CFLAGS) -o isim isim.c $(LDFLAGS)

msim:		msim.c
		$(CC) $(CFLAGS) -o msim msim.c $(LDFLAGS)

lsim:		lsim.c
		$(CC) $(CFLAGS) -o lsim lsim.c $(LDFLAGS)

slink:		slgetvect.o slinputcode.o sloldnuclear.o\
		    slink.o slautomodified.o slsexmodified.o iostuff.o
		$(CC) $(CFLAGS) -o slink \
		    slgetvect.o slinputcode.o sloldnuclear.o \
		    slink.o slautomodified.o slsexmodified.o iostuff.o \
		    $(LDFLAGS)

sloslink:	slgetvect.o slinputcode.o sloldnuclear.o iostuff.o\
		    sloslink.o slsloautomodified.o slslosexmodified.o iostuff.o
		$(CC) $(CFLAGS) $(SFLAG) -o slink \
		    slgetvect.o slinputcode.o sloldnuclear.o iostuff.o \
		    sloslink.o slsloautomodified.o slslosexmodified.o \
		    $(LDFLAGS)

sloslink.o:	slink.c        commondefs.h sldefs.h  moddefs.h
		$(CC) $(CFLAGS) $(SFLAG) -c -o sloslink.o \
                    slink.c

iostuff.o:	iostuff.c
slink.o:	slink.c      commondefs.h sldefs.h  moddefs.h
slgetvect.o:	slgetvect.c    commondefs.h sldefs.h
slinputcode.o:	slinputcode.c  commondefs.h sldefs.h
sloldnuclear.o:	sloldnuclear.c commondefs.h sldefs.h
slautomodified.o:	slautomodified.c  commondefs.h moddefs.h sldefs.h
slsexmodified.o:	slsexmodified.c   commondefs.h moddefs.h sldefs.h
slsloautomodified.o:	slsloautomodified.c commondefs.h moddefs.h sldefs.h
slslosexmodified.o:	slslosexmodified.c commondefs.h moddefs.h  sldefs.h

clean:
	rm *.o
	rm slink unknown isim msim lsim






