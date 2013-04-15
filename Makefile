#################################################################
########################   MakeVars   ###########################
#################################################################

CC=g++
 
CC_OPT=-Wall -Wextra -I"./include" -O3

LD_OPT=-lm

MKDIR=mkdir -p ./obj
 
CIBLE=mc_cpp
 
SRC=$(wildcard ./src/*.cpp)
 
OBJ=$(patsubst ./src/%.cpp,./obj/%.o,$(SRC))
 
#################################################################
########################   Makefile   ###########################
#################################################################
 
all:$(CIBLE)
	@echo "Compilation Success"

$(CIBLE):Makefile
 
./obj/%.o:./src/%.cpp 
	@$(MKDIR)
	$(CC) $(CC_OPT) -c $< -o $@ 

$(CIBLE):$(OBJ)
	$(CC) $(dOBJ) $(fOBJ) $(OBJ) -o $@ $(LD_OPT)

clean:
	rm -f $(CIBLE) ./obj/*.o *.dat *.xyz

rs:
	rsync -avz * ~/Dropbox/monte_carlo/mc_cpp/ 
