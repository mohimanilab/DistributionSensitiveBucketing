CC=g++
CXXFLAGS=-std=c++11
LDFLAGS=-g
OPTIMIZE=-O3
STATIC=-static
CCSERVER=/projects/mohimanilab/anaconda2/bin/g++

.SILENT:

GEN=../simulated_scripts
SCRIPT=../simulated_scripts

COMPILE_CMD=$(CC) $(CXXFLAGS) $(OPTIMIZE) $(LDFLAGS)
COMPILE_CMD_ORI=$(CC) $(CXXFLAGS) $(OPTIMIZE) $(LDFLAGS)
COMPILE_SERVER=$(CCSERVER) $(CXXFLAGS) $(OPTIMIZE) $(LDFLAGS)

# Check system type. If Linux, we make a static/stand-alone package
# Otherwise OSX, then we just do dynamic link
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
		COMPILE_CMD += -D LINUX $(STATIC)
endif
ifeq ($(UNAME_S),Darvin)
		COMPILE_CMD += -D OSX
endif

default: main gen 

#server: COMPILE_CMD=$(COMPILE_SERVER)
#server: main gen

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
		if [ -f DSBMain ]; then echo -e "removing DSBMain"; rm -f DSBMain; fi
		if [ -d DSBMain.dSYM ]; then rm -rf DSBMain.dSYM; fi
		if [ -f DataGeneration ]; then echo -e "removing DataGeneration"; rm -f DataGeneration; fi
		if [ -d DataGeneration.dSYM ]; then rm -rf DataGeneration.dSYM; fi

current:
		echo "Current compiling command is '$(COMPILE_CMD)'"
