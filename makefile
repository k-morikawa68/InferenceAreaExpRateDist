CXX = g++
CXXFLAGS = -std=c++17 -I/usr/include/eigen3 
OBJS = \
	   main.o\
	   HalfEdgeMesh.o\
	   AreaExpansionRate.o\

PROGRAM = AER

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(PROGRAM).out

clean: 
	rm -f *.o *.out
