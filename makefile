# The compiler to use
CXX = g++

# Compiler flags 
CXXFLAGS = -std=c++11 -fmax-errors=2 
FASTFLAGS = -O3 
DEBUGFLAGS = -g -O0 -Wall
OPTIMPATH = $(HOME)/optim/header_only_version/

# Linker flags
LDFLAGS = -lsfml-graphics -lsfml-window -lsfml-system -I/usr/include/eigen3 -I$(OPTIMPATH)

# The name of your output program
OUTPUT = run

# The name of your source file
SOURCE = main.cpp

all:
	$(CXX) $(CXXFLAGS) $(FASTFLAGS) $(SOURCE) -o $(OUTPUT) $(LDFLAGS)

debug:
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) $(SOURCE) -o $(OUTPUT) $(LDFLAGS)

depend:
	sudo apt-get install libsfml-dev libeigen3-dev -y

clean:
	rm -f $(OUTPUT)

