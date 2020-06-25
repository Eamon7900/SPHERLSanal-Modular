#Makefile for SPHERLSanal Modular

#make install : build every target and remove build directory 
#make all : build every target and leave object files in build directory
#make clean : remove build directory
#make uninstall remove the build directory and bin directory

#Set source, build, and include paths
BIN=bin/
OBJ=build/
SRC=src/
INC=include/
VPATH= $(SRC):$(INC)

#Flags for compile and link commands
CC= g++ 
CFLAGS= -w -c -I./$(INC) 
LFLAGS= -w -o 

MK_RAD_PRO_OBJS= $(OBJ)mkRadialProfile.o $(OBJ)eos.o $(OBJ)binfile.o $(OBJ)exception2.o

TARGS= $(BIN)mkRadPro 

#Phony targets:
.PHONY= all install dirs clean uninstall
all : $(TARGS) 

install : dirs all clean 

dirs : 
	@mkdir $(BIN)
	@mkdir $(OBJ)

clean : 
	@rm -rf $(OBJ)   

uninstall : clean 
	@rm -rf $(BIN)

#Targets and dependancies
$(BIN)mkRadPro : $(MK_RAD_PRO_OBJS) 
	$(CC) $(LFLAGS) $(BIN)mkRadPro $(MK_RAD_PRO_OBJS)  

$(OBJ)mkRadialProfile.o : mkRadialProfile.cpp eos.h binfile.h exception2.h paths.h
	$(CC) $(CFLAGS) $(SRC)mkRadialProfile.cpp -o $@ 

$(OBJ)eos.o : eos.cpp eos.h exception2.h
	$(CC) $(CFLAGS) $(SRC)eos.cpp -o $@ 

$(OBJ)binfile.o : binfile.cpp binfile.h
	$(CC) $(CFLAGS) $(SRC)binfile.cpp -o $@ 

$(OBJ)exception2.o : exception2.cpp exception2.h
	$(CC) $(CFLAGS) $(SRC)exception2.cpp -o $@  

