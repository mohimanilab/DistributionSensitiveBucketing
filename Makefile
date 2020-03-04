CC=g++
CXXFLAGS=-std=c++11
LDFLAGS=-g
OPTIMIZE=-O3
CCSERVER=/projects/mohimanilab/anaconda2/bin/g++

.SILENT:

GEN=../simulated_scripts
SCRIPT=../simulated_scripts

COMPILE_CMD=$(CC) $(CXXFLAGS) $(OPTIMIZE) $(LDFLAGS)
COMPILE_CMD_ORI=$(CC) $(CXXFLAGS) $(OPTIMIZE) $(LDFLAGS)
COMPILE_SERVER=$(CCSERVER) $(CXXFLAGS) $(OPTIMIZE) $(LDFLAGS)

default: main gen 

server: COMPILE_CMD=$(COMPILE_SERVER)
server: main

main: DSBMain.cpp
	echo "Compiling data-dependent method with post-processing..."
	$(COMPILE_CMD) DSBMain.cpp -o DSBMain
	echo "Done."

gen: DataGeneration.cpp
	echo "Compiling simulation data generation..."
	$(COMPILE_CMD) DataGeneration.cpp -o DataGeneration
	echo "Done."

format: DSBMain.cpp DataGeneration.cpp
	clang-format -style=file -i DSBMain.cpp DataGeneration.cpp

clean:
	rm -f DSBMain
	if [ -d DSBMain.dSYM ]; then
		rm -rf DSBMain.dSYM
	fi
	rm -f DataGeneration
	if [ -d DataGeneration.dSYM ]; then
		rm -rf DataGeneration.dSYM
	fi

current:
	echo "Current compiling command is '$(COMPILE_CMD)'"
