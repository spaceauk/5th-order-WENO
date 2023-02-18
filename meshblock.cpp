#include <fstream>
#include <string>
#include "defs.hpp"
#include <exception>

void meshblock::setSize(int nx, int ny, int nvar, real lenx, real leny) {
	this->nx=nx;  
	this->ny=ny; 
	this->nvar=nvar;

	// Use new to allocate memory to array for conservative variable
	U = new real**[nx]; // [r,ru,rv,rw,E,Bx,By,Bz]
	Us= new real**[nx];		    
	W = new real**[nx]; // [r,u,v,w,p,Bx,By,Bz]
	dwdx = new real**[nx];			   
	dwdy = new real**[nx];			    
	wxL = new real**[nx];			    
	wxR = new real**[nx];			    
	wyL = new real**[nx];			    
	wyR = new real**[nx];			    
	res = new real**[nx];			          		

	for (int i=0; i<nx; i++) {
		U[i] = new real*[ny];
		Us[i]= new real*[ny];
		W[i] = new real*[ny];
		dwdx[i] = new real*[ny];
		dwdy[i] = new real*[ny];
		wxL[i] = new real*[ny];
		wxR[i] = new real*[ny];
		wyL[i] = new real*[ny];
		wyR[i] = new real*[ny];
		res[i] = new real*[ny];

		for (int j =0; j<ny; j++) {
			U[i][j] = new real[nvar];
			Us[i][j]= new real[nvar];
			W[i][j] = new real[nvar];
			dwdx[i][j] = new real[nvar];
			dwdy[i][j] = new real[nvar];
			wxL[i][j] = new real[nvar];
			wxR[i][j] = new real[nvar];
			wyL[i][j] = new real[nvar];
			wyR[i][j] = new real[nvar];
			res[i][j] = new real[nvar];
		}
	}

	// Calculate quantities (note boundary cells)
	dx = lenx/(nx-2);
	dy = leny/(ny-2);
}

void meshblock::setParam(real gamma, real tEnd) {
	this->gamma=gamma;
	this->tEnd=tEnd;
}

void meshblock::setType(string IC, string limiter, string fluxMth) {
	this->IC=IC;
	this->limiter=limiter;
	this->fluxMth=fluxMth;
}


// Apply B.C.s
void meshblock::setBCs() {

	for (int k=0; k<nvar; k++) {				
		for (int i=0; i<nx; i++) {
			Us[i][0][k]=Us[i][1][k];
			Us[i][ny-1][k]=Us[i][ny-2][k];
			U[i][0][k]=U[i][1][k];
			U[i][ny-1][k]=U[i][ny-2][k];
			W[i][0][k]=W[i][1][k];
			W[i][ny-1][k]=W[i][ny-2][k];
		}
		for (int j=0; j<ny; j++) {
			Us[nx-1][j][k]=Us[nx-2][j][k];
			U[nx-1][j][k]=U[nx-2][j][k];
			W[nx-1][j][k]=W[nx-2][j][k];	
			Us[0][j][k]=Us[1][j][k];
			U[0][j][k]=U[1][j][k];
			W[0][j][k]=W[1][j][k];
		}

	}
}

// Conversion between U & W
void meshblock::U2W(string loc) {
	for (int i=0; i<nx; i++) {
		for (int j=0; j<ny; j++) {
			W[i][j][0]=U[i][j][0];
			W[i][j][1]=U[i][j][1]/U[i][j][0];
			W[i][j][2]=U[i][j][2]/U[i][j][0];
			W[i][j][3]=U[i][j][3]/U[i][j][0];
			W[i][j][4]=(gamma-1)*(U[i][j][4]-0.5*W[i][j][0]*
				   MAG(W[i][j][1],W[i][j][2],W[i][j][3]));
			if (nvar==8) {
				W[i][j][5]=U[i][j][5];
				W[i][j][6]=U[i][j][6];
				W[i][j][7]=U[i][j][7];
				W[i][j][4]=W[i][j][4]-(gamma-1)*(0.5*MAG(U[i][j][5],U[i][j][6],U[i][j][7]));
			}
			// Check for negative pressure
			if (W[i][j][4]<0 and i!=0 and i!=nx-1 and j!=0 and j!=ny-1) {
				cout<<"Error as negative pressure! P="<<W[i][j][4]<<" at i="<<i<<" & j="<<j<<". \n";
				cout<<"U2W in function: "<<loc<<" at loop count="<<count<<endl;
				throw exception();
			}
		}
	}
}
void meshblock::Us2W(string loc) {
        for (int i=0; i<nx; i++) {
                for (int j=0; j<ny; j++) {
                        W[i][j][0]=Us[i][j][0];
                        W[i][j][1]=Us[i][j][1]/Us[i][j][0];
                        W[i][j][2]=Us[i][j][2]/Us[i][j][0];
                        W[i][j][3]=Us[i][j][3]/Us[i][j][0];
                        W[i][j][4]=(gamma-1)*(Us[i][j][4]-0.5*W[i][j][0]*
                                   MAG(W[i][j][1],W[i][j][2],W[i][j][3]));
			if (nvar==8) {
	                        W[i][j][5]=Us[i][j][5];
        	                W[i][j][6]=Us[i][j][6];
                	        W[i][j][7]=Us[i][j][7];
				W[i][j][4]=W[i][j][4]-(gamma-1)*(0.5*MAG(Us[i][j][5],Us[i][j][6],Us[i][j][7]));
			}
			// Check for negative pressure
			if (W[i][j][4]<0 and i!=0 and i!=nx-1 and j!=0 and j!=ny-1) {
				cout<<"Error as negative pressure! P="<<W[i][j][4]<<" at i="<<i<<" & j="<<j<<". \n";
				cout<<"Us2W in function: "<<loc<<" at loop count="<<count<<endl;
				throw exception();
			}
                }
        }
}
void meshblock::W2U() {
	for (int i=0; i<nx; i++) {
		for (int j=0; j<ny; j++) {
			U[i][j][0]=W[i][j][0];
			U[i][j][1]=W[i][j][1]*W[i][j][0];
			U[i][j][2]=W[i][j][2]*W[i][j][0];
			U[i][j][3]=W[i][j][3]*W[i][j][0];
			U[i][j][4]=(W[i][j][4]/(gamma-1))+0.5*W[i][j][0]*
				   MAG(W[i][j][1],W[i][j][2],W[i][j][3]);
			if (nvar==8) {
				U[i][j][5]=W[i][j][5];
				U[i][j][6]=W[i][j][6];
				U[i][j][7]=W[i][j][7];
				U[i][j][4]=U[i][j][4]+0.5*MAG(W[i][j][5],W[i][j][6],W[i][j][7]);
			}
		}
	}
}


// Calculate wave speed
void meshblock::wavespeed(int i, int j) {
	speed=sqrt(MAG(W[i][j][1],W[i][j][2],W[i][j][3]));	
	c=sqrt(gamma*W[i][j][4]/W[i][j][0]);
	ca=sqrt((MAG(W[i][j][5],W[i][j][6],W[i][j][7]))/W[i][j][0]);
	cax=sqrt((SQR(W[i][j][5]))/W[i][j][0]);
	cay=sqrt((SQR(W[i][j][6]))/W[i][j][0]);
	real intm=SQR(c)+SQR(ca);		
	cfx=sqrt(0.5*(SQR(c)+SQR(ca))+0.5*sqrt(SQR(intm)
			-4*(SQR(c)*SQR(cax))));			
	cfy=sqrt(0.5*(SQR(c)+SQR(ca))+0.5*sqrt(SQR(intm)
			-4*(SQR(c)*SQR(cay))));	
}


// Delete memory from arrays

