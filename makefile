SHELL = /bin/sh

# Select compliler
CMPR = g++

CPPFLAGS = 

CPP =
CFLAGS = -std=c++11 -Wall

OBJ = main.o meshblock.o IC2Dtype.o MUSCL2D.o slopelimiter.o celledges.o riemannS.o savedata.o WENO2D.o  
EXEC = main.x

%.o:%.cpp
	$(CMPR) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

# Link into an executable
main: $(OBJ)
	$(CMPR) $(OBJ) -o $(EXEC)

clean:
	rm -f ./*.o ./*.x
