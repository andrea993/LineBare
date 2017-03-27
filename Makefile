CXXFLAGS=-g -pthread -O3 -Wall -msse3 -std=c++14
LDFLAGS=-g
CXX=g++

LineBare: main.o Rls.o
	g++ $(LDFLAGS) -o LineBare main.o Rls.o

main: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

Rls: Rls.h Rls.cpp
	$(CXX) $(CXXFLAGS) -c Rls.cpp
	
	
clean:
	-rm *.o
	-rm LineBare

	
