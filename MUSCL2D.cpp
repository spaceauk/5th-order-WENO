#include <vector>
#include <cmath>
#include <exception>
#include "defs.hpp"

void slopelimiter(meshblock &dom);
void celledges(meshblock &dom);
void riemannS(string fluxMth,vector<real> wL,vector<real> wR,real gamma,char direc,vector<real> &flux);
void savearray(meshblock dom, real*** array, string arrname);

void MUSCL2D(meshblock &dom, int step) {
	// Convert U to W
	if (step==1) {
		dom.U2W("MUSCL");
	} else if (step==2) {
		dom.Us2W("MUSCL");
	}
	dom.setBCs();

	// Apply slope limiter
	slopelimiter(dom);	
	celledges(dom);	

	// Constrained transport method	
	// if (dom.nvar==8) {}
	
	// Compute residuals
	vector<real> wL(dom.nvar),wR(dom.nvar),flux(dom.nvar); 
	for (int i=0; i<dom.nx; i++) {
		for (int j=0; j<dom.ny; j++) {
			for (int k=0; k<dom.nvar; k++) {
				dom.res[i][j][k]=0;
			}
		}
	}
	// (a) In x-direction
	for (int i=2; i<dom.nx-1; i++) {
		for (int j=1; j<dom.ny-1; j++) {
			for (int k=0; k<dom.nvar; k++) {
				wL[k]=dom.wxL[i-1][j][k];
				wR[k]=dom.wxR[i][j][k];
			}
			// Compute flux at i+1/2
			riemannS(dom.fluxMth,wL,wR,dom.gamma,'x',flux);
			// Contribution to the residual of cell (i,j)	
			for (int k=0; k<dom.nvar; k++) {
				dom.res[i-1][j][k]=dom.res[i-1][j][k]+flux[k]/dom.dx;
				dom.res[i][j][k]=dom.res[i][j][k]-flux[k]/dom.dx;
			if (isnan(flux[k])) {
				cout<<"At X: i="<<i<<" j="<<j<<" k="<<k<<":"<<flux[k]<<" dwdx="<<dom.dwdx[i][j][k]<<endl;
				throw exception();
			}
			}
		}
	}
	// (b) In y-direction
	for (int i=1; i<dom.nx-1; i++) {
		for (int j=2; j<dom.ny-1; j++) {
			for (int k=0; k<dom.nvar; k++) {
				wL[k]=dom.wyL[i][j-1][k];
				wR[k]=dom.wyR[i][j][k];
			}
			// Compute flux at j+1/2
			riemannS(dom.fluxMth,wL,wR,dom.gamma,'y',flux);
			// Contribution to the residual of cell (i,j)
			for (int k=0; k<dom.nvar; k++) {
				dom.res[i][j-1][k]=dom.res[i][j-1][k]+flux[k]/dom.dy;
				dom.res[i][j][k]=dom.res[i][j][k]-flux[k]/dom.dy;
			if (isnan(flux[k])) {
				cout<<"At Y: i="<<i<<" j="<<j<<" k="<<k<<":"<<flux[k]<<endl;
				throw exception();
			}
			}
		}
	}
	// Set residual for B.C.s
	// (a) Flux contribution of the most NORTH face: j=M-1
	for (int i=1; i<dom.nx-1; i++) {
		for (int k=0; k<dom.nvar; k++) {
			wL[k]=dom.wyL[i][dom.ny-2][k];
			wR[k]=dom.wyR[i][dom.ny-2][k];
		}
		riemannS(dom.fluxMth,wL,wR,dom.gamma,'y',flux);
		for (int k=0; k<dom.nvar; k++) {
			dom.res[i][dom.ny-2][k]=dom.res[i][dom.ny-2][k]+flux[k]/dom.dy;
		}
	}
	// (b) Flux contribution of the most SOUTH face: j=2
	for (int i=1; i<dom.nx-1; i++) {
		for (int k=0; k<dom.nvar; k++) {
			wL[k]=dom.wyL[i][1][k];
			wR[k]=dom.wyR[i][1][k];
		}
		riemannS(dom.fluxMth,wL,wR,dom.gamma,'y',flux);
		for (int k=0; k<dom.nvar; k++) {
			dom.res[i][1][k]=dom.res[i][1][k]-flux[k]/dom.dy;
		}
	}
	// (c) Flux contribution of the most EAST face: i=N-1
	for (int j=1; j<dom.ny-1; j++) {
		for (int k=0; k<dom.nvar; k++) {
			wL[k]=dom.wxL[dom.nx-2][j][k];
			wR[k]=dom.wxR[dom.nx-2][j][k];
		}
		riemannS(dom.fluxMth,wL,wR,dom.gamma,'x',flux);
		for (int k=0; k<dom.nvar; k++) {
			dom.res[dom.nx-2][j][k]=dom.res[dom.nx-2][j][k]+flux[k]/dom.dx;
		}
	}
	// (d) Flux contribution of the most WEST face: i=2
	for (int j=1; j<dom.ny-1; j++) {
		for (int k=0; k<dom.nvar; k++) {
			wL[k]=dom.wxL[1][j][k];
			wR[k]=dom.wxR[1][j][k];
		}
		riemannS(dom.fluxMth,wL,wR,dom.gamma,'x',flux);
		for (int k=0; k<dom.nvar; k++) {
			dom.res[1][j][k]=dom.res[1][j][k]-flux[k]/dom.dx;
		}
	}
}
