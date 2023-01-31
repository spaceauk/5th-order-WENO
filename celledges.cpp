#include "defs.hpp"

void celledges(meshblock &dom) {
	// Obtain cell left and right edges
	for (int i=2; i<dom.nx-1; i++) { // For x
		for (int j=1; j<dom.ny-1; j++) {
			for (int k=0; k<dom.nvar; k++) {
				dom.wxL[i-1][j][k]=dom.W[i-1][j][k]+0.5*dom.dwdx[i-1][j][k];
				dom.wxR[i][j][k]=dom.W[i][j][k]-0.5*dom.dwdx[i][j][k];
			}
		}
	}

	for (int i=1; i<dom.nx-1; i++) { // For y
		for (int j=2; j<dom.ny-1; j++) {
			for (int k=0; k<dom.nvar; k++) {
				dom.wyL[i][j-1][k]=dom.W[i][j-1][k]+0.5*dom.dwdy[i][j-1][k];
				dom.wyR[i][j][k]=dom.W[i][j][k]-0.5*dom.dwdy[i][j][k];
			}
		}
	}

	// Account for B.C.s - must not involve ghost cells
	for (int k=0; k<dom.nvar; k++) {
		for (int i=1; i<dom.nx-1; i++) {
			dom.wyL[i][dom.ny-2][k]=dom.wyR[i][dom.ny-2][k];
			dom.wyR[i][1][k]=dom.wyL[i][1][k];
		}
		for (int j=1; j<dom.ny-1; j++) {
			dom.wxL[dom.nx-2][j][k]=dom.wxR[dom.nx-2][j][k];
			dom.wxR[1][j][k]=dom.wxL[1][j][k];
		}

	}

}	
