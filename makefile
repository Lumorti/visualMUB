# The compiler to use
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11

# Linker flags
LDFLAGS = -lsfml-graphics -lsfml-window -lsfml-system

# The name of your output program
OUTPUT = run

# The name of your source file
SOURCE = main.cpp

all:
	$(CXX) $(CXXFLAGS) $(SOURCE) -o $(OUTPUT) $(LDFLAGS)

clean:
	rm -f $(OUTPUT)

