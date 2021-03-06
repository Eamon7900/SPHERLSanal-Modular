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
CFLAGS= -Wall -c -I./$(INC) 
LFLAGS= -Wall -o 

MK_RAD_PRO_OBJS= $(OBJ)mkRadialProfile.o $(OBJ)eos.o $(OBJ)datafile.o $(OBJ)binfile.o $(OBJ)exception2.o $(OBJ)argparser.o
MK_2D_SLICE_OBJS= $(OBJ)mk2DSlice.o $(OBJ)eos.o $(OBJ)datafile.o $(OBJ)binfile.o $(OBJ)exception2.o $(OBJ)argparser.o

TARGS= $(BIN)mkRadPro $(BIN)mk2DSlice 

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

$(BIN)mk2DSlice : $(MK_2D_SLICE_OBJS) 
	$(CC) $(LFLAGS) $(BIN)mk2DSlice $(MK_2D_SLICE_OBJS)  

$(OBJ)binfile.o : binfile.h binfile.cpp
	$(CC) $(CFLAGS) $(SRC)binfile.cpp -o $@ 

$(OBJ)mk2DSlice.o : mk2DSlice.cpp eos.h binfile.h
	$(CC) $(CFLAGS) $(SRC)mk2DSlice.cpp -o $@ 

$(OBJ)mkRadialProfile.o : mkRadialProfile.cpp eos.h binfile.h exception2.h paths.h argparser.h
	$(CC) $(CFLAGS) $(SRC)mkRadialProfile.cpp -o $@ 

$(OBJ)eos.o : eos.cpp eos.h binfile.h exception2.h
	$(CC) $(CFLAGS) $(SRC)eos.cpp -o $@ 

$(OBJ)datafile.o : datafile.h binfile.h datafile.cpp 
	$(CC) $(CFLAGS) $(SRC)datafile.cpp -o $@ 

$(OBJ)exception2.o : exception2.cpp exception2.h
	$(CC) $(CFLAGS) $(SRC)exception2.cpp -o $@  

$(OBJ)argparser.o : argparser.cpp argparser.h
	$(CC) $(CFLAGS) $(SRC)argparser.cpp -o $@  
