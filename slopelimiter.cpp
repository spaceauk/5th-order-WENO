#include<exception>

#include "defs.hpp"

void minmod(real dw1,real dw2,real dw3,real &mm);

void slopelimiter(meshblock &dom) {
	real ddwdx,ddwdy;

	for (int i=1; i<dom.nx-1; i++) {
		for (int j=1; j<dom.ny-1; j++) {
			for (int k=0; k<dom.nvar; k++) {
				if (dom.limiter=="1st") {
					// Godunov 1st-order:
					dom.dwdx[i][j][k]=0;
					dom.dwdy[i][j][k]=0;
				} else if (dom.limiter=="MC") {
					// Monotonised Centralised (MC) limiter
					real dww=2*(dom.W[i][j][k]-dom.W[i-1][j][k]);
				        real dwe=2*(dom.W[i+1][j][k]-dom.W[i][j][k]);
					real dwc=0.5*(dom.W[i+1][j][k]-dom.W[i-1][j][k]);
				        minmod(dww,dwe,dwc,ddwdx);
					dom.dwdx[i][j][k]=ddwdx;
					real dws=2*(dom.W[i][j][k]-dom.W[i][j-1][k]);
					real dwn=2*(dom.W[i][j+1][k]-dom.W[i][j][k]);
					dwc=0.5*(dom.W[i][j+1][k]-dom.W[i][j-1][k]);
					minmod(dws,dwn,dwc,ddwdy);
					dom.dwdy[i][j][k]=ddwdy;
				} else {
					cout<<"Error as incorrect slope limiter chosen!"<<endl;
					throw exception();
				}
			}
		}
	}
}


void minmod(real dw1,real dw2,real dw3,real &mm) {
	int sign[3],s=0;
	real v[3]={dw1,dw2,dw3};
	
	for (int i=0;i<3;i++) {
		if (v[i]<0) {
			sign[i]=-1;
		} else if (v[i]>0) {
			sign[i]=1;
		} else if (v[i]==0) {
			sign[i]=0;
		}
		s=s+sign[i];
	}
	s=s/3.;
	if (abs(s)==1) {
		mm=s*min(min(abs(dw1),abs(dw2)),abs(dw3));
	} else {
		mm=0.;
	}
}
