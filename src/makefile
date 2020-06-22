CC= g++ 
CFLAGS= -Wall -c
LFLAGS= -Wall -o

EXEC_DIR= ../bin

MK_RAD_PRO_OBJS= mkRadialProfile.o eos.o binfile.o exception2.o
ALL= mkRadPro 

all : $(ALL) 

mkRadPro : $(MK_RAD_PRO_OBJS) 
	$(CC) $(LFLAGS) $(EXEC_DIR)/mkRadPro $(MK_RAD_PRO_OBJS)  

mkRadialProfile.o : mkRadialProfile.cpp eos.h binfile.h exception2.h paths.h
	$(CC) $(CFLAGS) mkRadialProfile.cpp 

eos.o : eos.cpp eos.h exception2.h
	$(CC) $(CFLAGS) eos.cpp 

binfile.o : binfile.cpp binfile.h
	$(CC) $(CFLAGS) binfile.cpp 

exception2.o : exception2.cpp exception2.h
	$(CC) $(CFLAGS) exception2.cpp

clean : 
	rm -f *.o