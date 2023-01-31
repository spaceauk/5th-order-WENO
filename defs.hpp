// Include C++ headers
#include <cmath>
#include <cstdint>

#include <iostream>
using namespace std;

// Declarations or constants
using real = float;
#define eps 1.0E-6
#define pi 3.14159265359


// Forward declarations needed for class files
class meshblock {
	public: 
		int nx, ny, nvar;
		// Fluid Properties
		real gamma;
		real tEnd;
		// Related arrays
		real*** U; // Conservative variables
		real*** Us;
	       	real*** W; // Primitive variables

		real*** dwdx;
		real*** dwdy;
		real*** wxL;
		real*** wxR;
		real*** wyL;
		real*** wyR;
		real*** res;
		// Types of solver & problem
		string IC;
		string limiter;
		string fluxMth;
	
		// Calulated quantities after input
		real dx, dy, dt;
		real speed, c;
		real ca, cax, cay, cfx, cfy, beta; // When there is M-field
		int count;	
		
	void setSize(int no_x, int no_y, int no_v, real lenx, real leny);
	void setParam (real gammaa, real tEnd0);
	void setType (string IC0, string limiter0, string fluxMth0);
	void setBCs();
	void U2W(string loc);
	void Us2W(string loc);
	void W2U();
	void wavespeed(int i, int j);
};


// All functions
#define MAG(x,y,z) (x*x+y*y+z*z)
#define SQR(x) (x*x)
