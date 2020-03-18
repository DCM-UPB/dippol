#=================================
#Do we want debug options for gcc?
#=================================
DEBUG='warning'
#YUKI='-D YUKI'
#OMP=-fopenmp
ifeq ($(DEBUG),'warning')
CFLAGS= -Wall -Wextra -Wno-deprecated -Wno-unused-parameter -O2 -I$(INC) -I$(INC_LOOS) -D DEBUG $(YUKI) $(OMP)#Compiler flags for debug run
else 
CFLAGS= -O2 -Wno-deprecated -Wno-unused-parameter -I$(INC) -I$(INC_LOOS) $(YUKI) $(OMP) #Compiler flags with optimization
endif
LDFLAGS= -lm -Wl,-rpath=${LOOS} -L${LOOS} -lloos #Extra libraries (math.h,...)

#g++ -Wno-deprecated -o dippol dip_pol.cpp  main.cpp  mol_atom.cpp  my_math.cpp  treat_inp.cpp  write.cpp tools.cpp -I/home/hossam/code/dippol/include -lm -I /home/hossam/.lib/loos/include -Wl,-rpath=/home/hossam/.lib/loos -L /home/hossam/.lib/loos -lloos

#=========
#Variables
#=========
LOOS= /home/hossam/.lib/loos
#wildcard  = required for wildcard uses in variables
#addprefix = The value of prefix is prepended to the front of each individual name 
#notdir = If the file name contains no slash, it is left unchanged. Otherwise, everything through the last slash is removed from it.
SHELL=bash
FC=g++
EXEC=dippol
SRC_DIR=src
INC=include
INC_LOOS=${LOOS}/include
OBJ_DIR=obj
SAFE_DIR=safe
PWD=$(shell pwd)
SRC= $(wildcard $(addprefix $(SRC_DIR)/, *.cpp)) #List of the sources
OBJ= $(addprefix $(OBJ_DIR)/,$(notdir $(SRC:.cpp=.o))) #List of the object files 

#==========
#C commands
#==========
# @$(command) here it is @$(FC) = silent mode. No line is written if everything is ok.
# $@ = Name of the target (name before the column = yuki, ...)
# $< = Name of the first dependance (name after the column = $(OBJ), ...)
# % = patern substitution

#Creation of the execuatble
$(EXEC): $(OBJ) #If target (main) is older than $(OBJ)
	 $(FC) -o $@ $(OBJ) $(CFLAGS) $(LDFLAGS)


#Creation of the objects
#Produce .o files if .cpp or .h changed
#%.o: %.cpp %.hpp %.h
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC)/*.h $(INC)/*.hpp
	$(FC) -c -o $@ $< $(CFLAGS)

#g++ -c src/write.cpp -o write.o -Iinclude -I/home/hossam/.lib/loos/include
#============
#============
#Rebuild the dependances (even if another file is name clean)
.PHONY: clean


#Remove the intermediate files
clean: 
	rm -rf $(OBJ_DIR)/*.o

#Remove the intermediate files
very_clean: 
	rm -rf $(OBJ_DIR)/*.o $(SRC_DIR)/TAGS $(EXEC)

#Create tag files to navigate with emacs
tag:
	etags $(SRC_DIR)/*.cpp -o $(SRC_DIR)/TAGS

