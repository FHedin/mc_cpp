#
#  mc_cpp : A basic Monte Carlo simulations software.
#  Copyright (C) 2013  Florent Hedin
#  
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#################################################################
########################   MakeVars   ###########################
#################################################################

CC=g++
#CC=icpc
#CC=clang

# useful warnings for experimental/debug builds
# others (really verbose) : 
# -Wshadow -Wsign-conversion -Wzero-as-null-pointer-constant -Wswitch-default 
# -Wdouble-promotion -Wuninitialized -Winit-self -Wfloat-equal -Wpointer-arith 
# -Wcast-qual -Wcast-align -Wconversion -Wenum-compare -Wlogical-op 
# -Wredundant-decls -Winline

WRN_OPT=-Wall -Wextra  

OPT_OPT=-O2
INC_OPT=-I"./include"
CC_OPT=-std=c++11 $(WRN_OPT) $(OPT_OPT) $(INC_OPT)

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
	$(CC) $(CC_OPT) $(dOBJ) $(fOBJ) $(OBJ) -o $@ $(LD_OPT)

clean:
	rm -f $(CIBLE) ./obj/*.o *.dat *.xyz *.log

rs:
	rsync -avz . ~/Dropbox/monte_carlo/mc_cpp
