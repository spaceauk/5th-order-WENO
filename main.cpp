#include<fstream>
#include<iostream>
#include<cstdlib>
#include<math.h>
#include<exception> // Call error exception

#include "defs.hpp"

void IC2Dtype(meshblock &dom);
void MUSCL2D(meshblock &dom, int step);
void WENO2D(meshblock &dom, real*** Q);
void savedata(meshblock &dom);
void savearray(meshblock dom, real*** array, string arrname);

int main() {
	meshblock dom; // Domain using arrays
	string IC,limiter,fluxMth;
	real CFL;

	ifstream inputFile("input.txt");
	if (inputFile.is_open()) {
		// Input Domain-related properties
		int nx, ny, nvar;
		inputFile >> nx >> ny >> nvar;
		real lenx, leny;
		inputFile >> lenx >> leny;
		inputFile >> CFL;
		// Must account for 'ghost' cells thus, +2
		dom.setSize(nx+2,ny+2,nvar,lenx,leny);
		cout<<"Domain size is "<<dom.nx-2<<" x "<<dom.ny-2<<" with nvar="<<nvar<<endl;
		cout<<"CFL="<<CFL<<", dx="<<dom.dx<<" & dy="<<dom.dy<<endl;
		// Choose types of solver & problem
		if (nvar==5) {
			cout<<"-----Hydrodynamics solver-----"<<endl;
		} else if (nvar==8) {
			cout<<"-----MHD solver-----"<<endl;
		}
		inputFile >> IC >> limiter >> fluxMth;
		dom.setType(IC,limiter,fluxMth);
		cout<<"Problems chosen: "<<dom.IC<<endl;
		if (dom.fluxMth=="WENO") {dom.limiter="5th-order";}
		cout<<"Limiter selected: "<<dom.limiter<<endl;
		cout<<"Types of Riemann Solver chosen: "<<dom.fluxMth<<endl;
	} else {
		cout<<"Unable to open the file..."<<endl;
		cout<<"***Terminating program***"<<endl;
		throw exception();
	}


	// Set up initial conditions
	IC2Dtype(dom);
	dom.W2U();
	dom.setBCs();

	// Main loop--------------------------------------
	cout<<"Entering main loop..."<<endl;
	int count=0; dom.count=count;
	real lambda1=0, lambda2=0, lambda3=0, lambda4=0, lambda;
	real t=0.0, dt=0.0; 
	while (t<dom.tEnd and count<1000) {
		// RK2 1st step: update Us
		if (dom.fluxMth=="WENO") {
			WENO2D(dom,dom.U);
		} else {
			MUSCL2D(dom,1);
		}
		for (int i = 1; i < dom.nx-1; i++) {
		    for (int j = 1; j < dom.ny-1; j++) {
		        for (int k = 0; k < dom.nvar; k++) {
            			dom.Us[i][j][k] = dom.U[i][j][k] - dt*dom.res[i][j][k];
		        }
		    }
		}
		dom.setBCs();

		// RK2 2nd step: update U
		if (dom.fluxMth=="WENO") {
			WENO2D(dom,dom.Us);
		} else {
			MUSCL2D(dom,2);
		}
		for (int i = 1; i < dom.nx-1; i++) {
		    for (int j = 1; j < dom.ny-1; j++) {
		        for (int k = 0; k < dom.nvar; k++) {
            			dom.U[i][j][k] = 0.5*(dom.U[i][j][k] + dom.Us[i][j][k] - dt*dom.res[i][j][k]);
		        }
		    }
		}		
		dom.setBCs();

		// Calculate time for advancement
		dom.U2W("main");
		for (int i=1; i<dom.nx-1; i++) {
			for (int j=1; j<dom.ny-1; j++) {
				dom.wavespeed(i,j);
				lambda1=max(lambda1,abs(dom.speed+dom.cfx));
				lambda2=max(lambda2,abs(dom.speed-dom.cfx));
				lambda3=max(lambda3,abs(dom.speed+dom.cfy));				
				lambda4=max(lambda4,abs(dom.speed-dom.cfy));
			}
		}
		lambda=max(lambda1,lambda2);
		lambda=max(lambda,lambda3);
		lambda=max(lambda,lambda4);
		t=t+dt; dom.dt=dt;
		cout<<"At t="<<t<<" where dt="<<dt<<" and loop count at "<<count<<".\n"; 
		count=count+1; dom.count=count;
		dt=CFL*min(dom.dx,dom.dy)/lambda;
		if (dt<0 or isnan(dt)) {
			cout<<"Error as invalid dt produced!!!"<<endl;
			throw exception();
		}
	}

	// Output data
	dom.W2U();
	savedata(dom);	

	return 0;
}
