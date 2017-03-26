CCPPFLAGS=-g -pthread -I/sw/include/root 
LDFLAGS=-gPPFLAGS=-g -pthread -I/sw/include/root -O3
LDFLAGS=-g


main: main.o Rls.o
	g++ $(LDFLAGS) -o main main.o Rls.o

main.o: main.cpp
	g++ $(CPPFLAGS) -c main.cpp

Rls.o: Rls.h Rls.cpp
	g++ $(CPPFLAGS) -c Rls.cpp

	
