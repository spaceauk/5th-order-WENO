#include "defs.hpp"

void firstderiV(real delta,real intm[5][3],real derV[3]);

void viscousflux(meshblock &dom,int i, int j,
		real &tau_xx,real &tau_yy,real &tau_xy,real &q_x,real &q_y) {
	real k_suth=dom.gamma/((dom.gamma-1.)*dom.Pr*dom.Re); // Thermal conductivity
	// Temperature for viscosity
	real T=dom.W[i][j][4]/dom.W[i][j][0]; 
	// Sutherland viscosity
	real mu=(pow(T,1.5))*((1+dom.Suth)/(T+dom.Suth));
	// Bulk viscosity (compressible effects)
	real lambda=(-2./3.)*mu;

	// Obtain relevant derivatives (along x)
	real intm[5][3],derVx[3];
	if (i<=2 or i>=dom.nx-3) {
	// 2nd order accurate central difference for first spatial derivative
		derVx[0]=(dom.W[i+1][j][1]-dom.W[i-1][j][1])/(2*dom.dx);
		derVx[1]=(dom.W[i+1][j][2]-dom.W[i-1][j][2])/(2*dom.dx);
		derVx[2]=(dom.W[i+1][j][4]/dom.W[i+1][j][0]-dom.W[i-1][j][4]/dom.W[i-1][j][0])/(2*dom.dx);
	} else {
		for (int k=0; k<5; k++) {
			intm[k][0]=dom.W[i-2+k][j][1];
			intm[k][1]=dom.W[i-2+k][j][2];
			intm[k][2]=dom.W[i-2+k][j][4]/dom.W[i-2+k][j][0];
		}
		firstderiV(dom.dx,intm,derVx);
	}
	real u_x=derVx[0], v_x=derVx[1];
	real T_x=derVx[2];
	// Obtain relevant derivatives (along y)
	real derVy[3];
	if (j<=2 or j>=dom.ny-3) {
		// 2nd order accurate central difference for first spatial derivative
                derVy[0]=(dom.W[i][j+1][1]-dom.W[i][j-1][1])/(2*dom.dx);
                derVy[1]=(dom.W[i][j+1][2]-dom.W[i][j-1][2])/(2*dom.dx);
                derVy[2]=(dom.W[i][j+1][4]/dom.W[i][j+1][0]-dom.W[i][j-1][4]/dom.W[i][j-1][0])/(2*dom.dy);
         } else {
		for (int k=0; k<5; k++) {
  		        intm[k][0]=dom.W[i][j-2+k][1];
      	                intm[k][1]=dom.W[i][j-2+k][2];
              	        intm[k][2]=dom.W[i][j-2+k][4]/dom.W[i][j-2+k][0];
                }
	}
        firstderiV(dom.dy,intm,derVy);
	real u_y=derVy[0], v_y=derVy[1];
	real T_y=derVy[2];

	// Normal stress rate (include bulk viscosity effects)
	tau_xx=(lambda*(u_x+v_y)+2*mu*u_x)/dom.Re;
	tau_yy=(lambda*(u_x+v_y)+2*mu*v_y)/dom.Re;
	// Shear stress rate
	tau_xy=(mu*(u_y+v_x))/dom.Re;
	// Heat transfer rate due to Fourier's laws of conduction
	q_x=k_suth*T_x;
	q_y=k_suth*T_y;
}

// First spatial derivative ()
void firstderiV(real delta,real intm[5][3],real derV[3]) {

	for (int k=0; k<3; k++) {
		// 4th order accurate central difference for first spatial derivative
		derV[k]=((1./12.)*intm[0][k]-(2./3.)*intm[1][k]+(2./3.)*intm[3][k]-(1./12.)*intm[4][k])/delta;
	}
}
