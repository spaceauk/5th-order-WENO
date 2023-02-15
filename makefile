# Select compliler
CMPR = g++

# Create directory
OBJDIR = ./obj/

CPPFLAGS = 

CPP =
CFLAGS = -std=c++17 -Wall

OBJ = main.o meshblock.o IC2Dtype.o MUSCL2D.o slopelimiter.o celledges.o riemannS.o savedata.o WENO2D.o  
EXEC = main.x

OBJS = $(addprefix $(OBJDIR), $(OBJ))

$(addprefix ./obj/, %.o): %.cpp
	$(CMPR) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

# Link into an executable
main: $(OBJ)
	$(CMPR) $(OBJ) -o $(EXEC)

clean:
	rm -f ./obj/*.o 
	rm -f ./*.x
