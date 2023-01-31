#include<vector>
#include<exception>

#include "defs.hpp"

void RUSA(vector<real> wL,vector<real> wR,real gamma,vector<real> &flux);

void riemannS(string fluxMth,vector<real> wL,vector<real> wR,real gamma,char direc,vector<real> &flux) {
	vector<real> wL_int(8),wR_int(8);
	for (int i=0; i<size(wL); i++) {
		wL_int[i]=wL[i]; wR_int[i]=wR[i];
	}
	if (direc=='x') {
	} else if (direc=='y'){
		wL_int[1]=wL[2]; // Swapped velocity direction
		wL_int[2]=wL[1];
		wR_int[1]=wR[2]; 
		wR_int[2]=wR[1];
		if (size(wL)==8) {
			wL_int[5]=wL[6]; // Swapped magnetic direction
			wL_int[6]=wL[5];
			wR_int[5]=wR[6]; 
			wR_int[6]=wR[5];
		}
	} else {
		cout<<"Wrong direction("<<direc<<") selected!!! \n";
		throw exception();
	}

	vector<real> flux_int(size(flux));
	if (fluxMth=="RUSA") {
		RUSA(wL_int,wR_int,gamma,flux_int);
	} else {
		cout<<"Error as no appropriate Riemann solver chosen! \n";
		throw exception();
	}

	for (int i=0; i<size(wL); i++){flux[i]=flux_int[i];}
	if (direc=='x') {
	} else if (direc=='y'){
		flux[1]=flux_int[2]; // Swapped velocity direction
		flux[2]=flux_int[1];
		if (size(wL)==8) {
			flux[5]=flux_int[6]; // Swapped magnetic firection
			flux[6]=flux_int[5];
		}
	}
}

void RUSA(vector<real> wL,vector<real> wR,real gamma,vector<real> &flux){
	// Left states
	real vnL=wL[1];
	real BnL=wL[5];
	real caL=sqrt(MAG(wL[5],wL[6],wL[7])/wL[0]);
	real aL=sqrt(gamma*wL[4]/wL[0]);
	real canL=sqrt(SQR(BnL)/wL[0]);
	real intm0=SQR(aL)+SQR(caL);
	real cfL=sqrt(0.5*(SQR(aL)+SQR(caL))+0.5*sqrt(SQR(intm0)
			-4*SQR(aL)*SQR(canL)));
	real eL=(wL[4]/(wL[0]*(gamma-1)))+0.5*MAG(wL[1],wL[2],wL[3])
			+0.5*MAG(wL[5],wL[6],wL[7])/wL[0];
	eL=eL*wL[0];
	// Specific enthalpy
	real HL=(eL+wL[4]+0.5*MAG(wL[5],wL[6],wL[7]))/wL[0];
	vector<real> qL(size(wL)), FL(size(wL));
	qL=wL;
	qL[1]=wL[0]*wL[1]; qL[2]=wL[0]*wL[2]; qL[3]=wL[0]*wL[3];
	qL[4]=eL;

	// Right states
        real vnR=wR[1];
        real BnR=wR[5];
        real caR=sqrt(MAG(wR[5],wR[6],wR[7])/wR[0]);
        real aR=sqrt(gamma*wR[4]/wR[0]);
        real canR=sqrt(SQR(BnR)/wR[0]);
	intm0=SQR(aR)+SQR(caR);
        real cfR=sqrt(0.5*(SQR(aR)+SQR(caR))+0.5*sqrt(SQR(intm0)
                        -4*SQR(aR)*SQR(canR)));
        real eR=(wR[4]/(wR[0]*(gamma-1)))+0.5*MAG(wR[1],wR[2],wR[3])
                        +0.5*MAG(wR[5],wR[6],wR[7])/wR[0];
        eR=eR*wR[0];
        // Specific enthalpy
        real HR=(eR+wR[4]+0.5*MAG(wR[5],wR[6],wR[7]))/wR[0];
        vector<real> qR(size(wR)), FR(size(wR));
        qR=wR;
        qR[1]=wR[0]*wR[1]; qR[2]=wR[0]*wR[2]; qR[3]=wR[0]*wR[3];
        qR[4]=eR;

	// Left fluxes
	real ptL=wL[4]+0.5*MAG(wL[5],wL[6],wL[7]);
	FL[0]=wL[0]*vnL;
	FL[1]=wL[0]*vnL*wL[1]+ptL-BnL*wL[5];
	FL[2]=wL[0]*vnL*wL[2]-BnL*wL[6];
	FL[3]=wL[0]*vnL*wL[3]-BnL*wL[7];
	FL[4]=wL[0]*vnL*HL-BnL*(wL[1]*wL[5]+wL[2]*wL[6]+wL[3]*wL[7]);
	FL[5]=vnL*wL[5]-BnL*wL[1];
	FL[6]=vnL*wL[6]-BnL*wL[2];
	FL[7]=vnL*wL[7]-BnL*wL[3];
	// Right fluxes
        real ptR=wR[4]+0.5*MAG(wR[5],wR[6],wR[7]);
        FR[0]=wR[0]*vnR;
        FR[1]=wR[0]*vnR*wR[1]+ptR-BnR*wR[5];
        FR[2]=wR[0]*vnR*wR[2]-BnR*wR[6];
        FR[3]=wR[0]*vnR*wR[3]-BnR*wR[7];
        FR[4]=wR[0]*vnR*HR-BnR*(wR[1]*wR[5]+wR[2]*wR[6]+wR[3]*wR[7]);
        FR[5]=vnR*wR[5]-BnR*wR[1];
        FR[6]=vnR*wR[6]-BnR*wR[2];
        FR[7]=vnR*wR[7]-BnR*wR[3];

	// Compute RUSA flux
	real intm1=cfL+abs(wL[1]), intm2=cfR+abs(wR[1]);
	real intm=max(intm1,intm2);
	for (int k=0; k<size(wL); k++) {
		flux[k]=0.5*(FL[k]+FR[k])-0.5*intm*(qR[k]-qL[k]);
	}
}
