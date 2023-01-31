
#include "defs.hpp"

void WENOflux(real*** Q,real gamma,string direc,real dd,int nx,int ny,int nvar,real*** res);

void WENO2D (meshblock &dom, real*** Q) {

	real resx[dom.nx][dom.ny][dom.nvar];

	// Compute residuals x-direction
	WENOflux(Q,dom.gamma,"x",dom.dx,dom.nx,dom.ny,dom.nvar,dom.res);
	for (int i=0; i<dom.nx; i++) {
		for (int j=0; j<dom.ny; j++) {
			for (int k=0; k<dom.nvar; k++) {
				resx[i][j][k]=dom.res[i][j][k];
			}
		}
	}

	// Compute residuals y-direction
	WENOflux(Q,dom.gamma,"y",dom.dy,dom.nx,dom.ny,dom.nvar,dom.res);
	for (int i=0; i<dom.nx; i++) {
		for (int j=0; j<dom.ny; j++) {
			for (int k=0; k<dom.nvar; k++) {
				dom.res[i][j][k]=dom.res[i][j][k]+resx[i][j][k];
			}
		}
	}

}

void WENOflux(real*** Q,real gamma,string direc,real dd,int nx,int ny,int nvar,real*** res){
	real r,u,v,w,Bx=0.,By=0.,Bz=0.,p,pt,H;
	real vn,Bn=0.;
	real normx,normy;
	real a,ca,can,cf; // Calculate wavespeeds
	if (direc=="x"){
		normx=1; normy=0;
	} else if (direc=="y"){
		normx=0; normy=1;
	}

	// Compute flux from conservative vectors
	real flux[nx][ny][nvar];
	real lambda=0;
	for (int i=0; i<nx; i++) {
                for (int j=0; j<ny; j++) {
                	r = Q[i][j][0];
			u = Q[i][j][1]/r;
			v = Q[i][j][2]/r;
			w = Q[i][j][3]/r;
			vn = u*normx + v*normy;
			p=(gamma-1)*(Q[i][j][4]-0.5*r*MAG(u,v,w));
                        if (nvar==8) {
                                Bx=Q[i][j][5];
                                By=Q[i][j][6];
                                Bz=Q[i][j][7];
                                p=p-(gamma-1)*(0.5*MAG(Q[i][j][5],Q[i][j][6],Q[i][j][7]));
				Bn = Bx*normx + By*normy;
			}
			H=(Q[i][j][4]+p+MAG(Bx,By,Bz))/r;
			pt=p+0.5*MAG(Bx,By,Bz);
			flux[i][j][0]=r*vn;
			flux[i][j][1]=r*vn*u + pt*normx;
			flux[i][j][2]=r*vn*v + pt*normy;
			flux[i][j][3]=r*vn*w + pt*0;
                        flux[i][j][4]=r*vn*H;
			if (nvar==8) {
				flux[i][j][5]=(vn*Bx-Bn*u);
				flux[i][j][6]=(vn*By-Bn*v);
				flux[i][j][7]=(vn*Bz-Bn*w);

				flux[i][j][1]=flux[i][j][1]-Bn*Bx;
				flux[i][j][2]=flux[i][j][2]-Bn*By;
				flux[i][j][3]=flux[i][j][3]-Bn*Bz;
				flux[i][j][4]=flux[i][j][4]-Bn*(u*Bx+v*By+w*Bz);
			}
			a = sqrt(gamma*p/r); // sound speed
			ca=sqrt(MAG(Bx,By,Bz)/r); // Alfven waves
			can=sqrt(SQR(Bn)/r);
			real intm=SQR(a)+SQR(ca);
			cf=sqrt(0.5*(SQR(a)+SQR(ca))+0.5*sqrt(SQR(intm)-4*(SQR(a)*SQR(can)))); // Fast magnetosonic waves	
			intm=abs(vn)+abs(cf);
			lambda=max(lambda,intm);
                }
        }

	// For flux splitting, use global Lax-Friedrichs flux splitting
	int ii,jj;
	real FL[nx][ny][nvar], FR[nx][ny][nvar];
	for (int i=0; i<nx; i++) {
		for (int j=0; j<ny; j++) {
			for (int k=0; k<nvar; k++) {
				FR[i][j][k]=0.5*(flux[i][j][k]+lambda*Q[i][j][k]);
				if (direc=="x"){ // Shift
					ii=i+1; jj=j;
					if (ii>=nx) { ii=ii-(nx-1);}
				} else if (direc=="y"){
					ii=i; jj=j+1;
					if (jj>=ny) { jj=jj-(ny-1);}
				}
				FL[i][j][k]=0.5*(flux[ii][jj][k]-lambda*Q[ii][jj][k]);
			}
		}
	}
	
	// Shift fluxes and calculate fluxes on left & right
	int iii,jjj;
	real FluxL[nx][ny][nvar],FluxR[nx][ny][nvar];
	for (int i=0; i<nx; i++) {
                for (int j=0; j<ny; j++) {
                        for (int k=0; k<nvar; k++) {
				// Shift Forward
                                if (direc=="x"){
                                        ii=i+1; jj=j;
                                        if (ii>=nx) { ii=ii-(nx-1);}
			                iii=i+2; jjj=j;
                                        if (iii>=nx) { iii=iii-(nx-1);}
                                } else if (direc=="y"){
                                        ii=i; jj=j+1;
                                        if (jj>=ny) { jj=jj-(ny-1);}
                                        iii=i; jjj=j+2;
                                        if (jjj>=ny) { jjj=jjj-(ny-1);}
                                }
                                real FLp=FL[ii][jj][k];
				real FLpp=FL[iii][jjj][k];
                                real FRp=FR[ii][jj][k];
				real FRpp=FR[iii][jjj][k];
	
				// Shift Backwards
                                if (direc=="x"){
                                        ii=i-1; jj=j;
                                        if (ii<0) { ii=ii+(nx-1);}
			                iii=i-2; jjj=j;
                                        if (iii<0) { iii=iii+(nx-1);}
                                } else if (direc=="y"){
                                        ii=i; jj=j-1;
                                        if (jj<0) { jj=jj+(ny-1);}
                                        iii=i; jjj=j-2;
                                        if (jjj<0) { jjj=jjj+(ny-1);}
                                }
                                real FLm=FL[ii][jj][k];
				real FLmm=FL[iii][jjj][k];
                                real FRm=FR[ii][jj][k];
				real FRmm=FR[iii][jjj][k];

				// Obtain right flux---------------------
				//Polynomials (from eqn (3.14): High-order WENO FVM on Cartesian Grids with AMR)
				real p0n = (2*FRmm-7*FRm+11*FR[i][j][k])/6;
				real p1n = (-FRm+5*FR[i][j][k]+2*FRp)/6;
				real p2n = (2*FR[i][j][k]+5*FRp-FRpp)/6;
				// Smoothness indicators (Beta factors)
				real B0n = (13./12.)*pow((FRmm-2*FRm+FR[i][j][k]),2) + (1./4.)*pow((FRmm-4*FRm+3*FR[i][j][k]),2);
				real B1n = (13./12.)*pow((FRmm-2*FR[i][j][k]+FRp),2) + (1./4.)*pow((FRm-FRp),2);
				real B2n = (13./12.)*pow((FR[i][j][k]-2*FRp+FRpp),2) + (1./4.)*pow((3*FR[i][j][k]-4*FRp+FRpp),2);
				// Constants
				real d0n=1./10.; 
				real d1n=6./10.; 
				real d2n=3./10.; 
				real epsilon=pow(dd,(5-1));
			        // Alpha weights (from eqn (3.17): High-order WENO FVM on Cartesian Grids with AMR)
				real alpha0n = d0n/pow((epsilon+B0n),2);
				real alpha1n = d1n/pow((epsilon+B1n),2);
				real alpha2n = d2n/pow((epsilon+B2n),2);
				real alphasumn = alpha0n + alpha1n + alpha2n;
				// ENO stencils weights (from eqn (3.17): High-order WENO FVM on Cartesian Grids with AMR)
				real w0n = alpha0n/alphasumn;
				real w1n = alpha1n/alphasumn;
				real w2n = alpha2n/alphasumn;
				// Numerical Flux at cell boundary
				FluxR[i][j][k] = w0n*p0n + w1n*p1n +w2n*p2n;

				// Obtain left flux----------------------
				// Polynomials (from eqn (3.14): High-order WENO FVM on Cartesian Grids with AMR)
				real p0p = (-FLmm+5*FLm+2*FL[i][j][k])/6;
				real p1p = (2*FLm+5*FL[i][j][k]-FLp)/6;
				real p2p = (11*FL[i][j][k]-7*FLp+2*FLpp)/6;
				// Smoothness indicators (Beta factors)
				real B0p = (13./12.)*pow((FLmm-2*FLm+FL[i][j][k]),2) + (1./4.)*pow((FLmm-4*FLm+3*FL[i][j][k]),2);
				real B1p = (13./12.)*pow((FLmm-2*FL[i][j][k]+FLp),2) + (1./4.)*pow((FLm-FLp),2);
				real B2p = (13./12.)*pow((FL[i][j][k]-2*FLp+FLpp),2) + (1./4.)*pow((3*FL[i][j][k]-4*FLp+FLpp),2);
				// Constants
				real d0p=3./10.; 
				real d1p=6./10.; 
				real d2p=1./10.; 
				epsilon=pow(dd,(5-1));
				// Alpha weights (from eqn (3.17): High-order WENO FVM on Cartesian Grids with AMR)
				real alpha0p = d0p/pow((epsilon+B0p),2);
				real alpha1p = d1p/pow((epsilon+B1p),2);
				real alpha2p = d2p/pow((epsilon+B2p),2);
				real alphasump = alpha0p + alpha1p + alpha2p;
				// ENO stencils weights (from eqn (3.17): High-order WENO FVM on Cartesian Grids with AMR)
				real w0p = alpha0p/alphasump;
				real w1p = alpha1p/alphasump;
				real w2p = alpha2p/alphasump;
				// Numerical Flux at cell boundary
				FluxL[i][j][k] = w0p*p0p + w1p*p1p +w2p*p2p;
                        }
                }
        }

	// Obtain FVM residual term, df/dx
	for (int i=1; i<nx-1; i++) {
		for (int j=1; j<ny-1; j++) {
			for (int k=0; k<nvar; k++) {
				if (direc=="x"){ // Shift
                                        ii=i-1; jj=j;
                                        if (ii<0) { ii=ii+(nx-1);}
                                } else if (direc=="y"){
                                        ii=i; jj=j-1;
                                        if (jj<0) { jj=jj+(ny-1);}
                                }
				res[i][j][k]=(FluxL[i][j][k]-FluxL[ii][jj][k] + FluxR[i][j][k]-FluxR[ii][jj][k])/dd;
			}
		}
	}

}
