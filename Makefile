## Created on 2025-07-23 at 18:16:17 CEST by David Gaspard (ORCID 0000-0002-4449-8782) <david.gaspard@espci.fr> under the MIT License.
## Makefile of the Waffle program.

### MUMPS Settings ###
MUMPS_DIR=/usr/lib/MUMPS_5.8.0
MUMPS_INCL=-I$(MUMPS_DIR)/include
MUMPS_LIBS=-L$(MUMPS_DIR)/lib -lzmumps -lmumps_common -lmetis -lpord -lesmumps -lscotch -lmpiseq -lgfortran
## Note that -lgfortran should be loaded in the last position.

### Compilation options ###
CC=g++
CFLAGS=-W -Wall -Wextra -std=c++17 $(MUMPS_INCL) -Isrc -Jbin -fopenmp -O2
LIBS=-lopenblas -lumfpack -lpng $(MUMPS_LIBS)

##########################
## BUILD MAIN PROGRAM
##########################
PROGNAME = waffle
SRCLIST = $(shell find src/ -name "*.cpp" ! -name "Main.cpp")
BINLIST = $(SRCLIST:src/%.cpp=bin/%.o)

all: directories $(PROGNAME)

$(PROGNAME): $(BINLIST) src/Main.cpp
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

bin/%.o: src/%.cpp 
	$(CC) $(CFLAGS) -c $< -o $@

directories:
	mkdir -p bin/

##########################
## BUILD TEST EXECUTABLE
##########################
TESTSRCLIST = $(shell find test/ -name "Test*.cpp")
TESTEXELIST = $(TESTSRCLIST:test/Test%.cpp=%.test)

test: directories $(TESTEXELIST)

%.test: $(BINLIST) test/Test%.cpp
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@

##################################
## CLEAN ALL BUILDS AND TESTS
##################################
clean:
	rm -rfv $(PROGNAME) bin/*.o $(TESTEXELIST)

###### END OF FILE ######
