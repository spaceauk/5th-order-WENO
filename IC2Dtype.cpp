#include<exception> // Call error exception
#include<vector>		    

#include "defs.hpp"

void IC2Dtype(meshblock &dom) {
	real tEnd, gamma;
	vector<real> rinit(4),pinit(4),Einit(4);
	vector<real> uinit(4),vinit(4),winit(4);
       	vector<real> Bxinit(4),Byinit(4),Bzinit(4);

	// Initialize the vectors to zero first
	rinit={0}; pinit={0}; Einit={0};
	uinit={0}; vinit={0}; winit={0};
	Bxinit={0}; Byinit={0}; Bzinit={0};

	// Set conditions for specific I.C.s selected
	if (dom.IC=="SST" and dom.nvar==5) {
		cout<<"2D Sod Shock tube config along x."<<endl;
		tEnd=0.2;
		gamma=2.0;
		pinit={0.1,1.0,1.0,0.1};
		rinit={0.125,1.0,1.0,0.125};

	} else if (dom.IC=="SBI" and dom.nvar==5) {
		cout<<"Shock-bubble interaction problem:"<<endl;
		cout<<"When a Mach 6 shock wave in air impacts on a cylindrical Helium bubble. Note that air & Helium are treated as the same ideal gas fluid for simplicity."<<endl;
		if (dom.nx==dom.ny) {
                        cout<<"Error domain size!"<<endl;
                        throw exception();
                }
		tEnd=0.15;
		pinit={1,41.83,1,0};
		rinit={1,5.268,0.138,0};
		uinit={-3,2.752,-3,0};		
		gamma=1.4;

	} else if (dom.IC=="BWx" and dom.nvar==8) {
		cout<<"2D Brio & Wu shocktube config along x."<<endl;
		tEnd=0.2;
		gamma=2.0;
		pinit={0.1,1.0,1.0,0.1};
		rinit={0.125,1.0,1.0,0.125};
		Bxinit={0.75,0.75,0.75,0.75};
		Byinit={-1.0,1.0,1.0,-1.0};

	} else if (dom.IC=="BWy" and dom.nvar==8) {
		cout<<"2D Brio & Wu shocktube config along y."<<endl;
		tEnd=0.2;
		gamma=2.0;
		pinit={0.1,0.1,1.0,1.0};
		rinit={0.125,0.125,1.0,1.0};
		Byinit={0.75,0.75,0.75,0.75};
		Bxinit={-1.0,-1.0,1.0,1.0};

	} else {
		cout<<"No valid initial condition chosen!"<<endl;
		cout<<"Look at IC2Dtype file to choose the correct I.C.s!"<<endl;
		cout<<"***Program is stopped!***"<<endl;
		throw exception();
	}

	// Assign values to arrays (Must account for boundary cells)
	real x, y;
	if (dom.IC=="rotor") {
	} else if (dom.IC=="SBI") {
		// pre-shocked air, post-shocked air, Helium bubble
		for (int i=1; i<dom.nx-1; i++) {
                        x=i*dom.dx;
                        for (int j=1; j<dom.ny-1; j++) {
                                y=j*dom.dy;
				real radius=sqrt(pow(x-0.25,2)+SQR(y));
				if (radius<=0.15) {
					dom.W[i][j][0]=rinit[2];
                                        dom.W[i][j][1]=uinit[2];
                                        dom.W[i][j][2]=vinit[2];
                                        dom.W[i][j][3]=winit[2];
                                        dom.W[i][j][4]=pinit[2];
				} else if (radius>0.15 and x<=0.05) {
					dom.W[i][j][0]=rinit[1];
                                        dom.W[i][j][1]=uinit[1];
                                        dom.W[i][j][2]=vinit[1];
                                        dom.W[i][j][3]=winit[1];
                                        dom.W[i][j][4]=pinit[1];
				} else if (radius>0.15 and x>0.05) {
					dom.W[i][j][0]=rinit[0];
                                        dom.W[i][j][1]=uinit[0];
                                        dom.W[i][j][2]=vinit[0];
                                        dom.W[i][j][3]=winit[0];
                                        dom.W[i][j][4]=pinit[0];
				} else {
					cout<<"Error as out designated domain!!!"<<endl;
					throw exception();
				}
			}
		}
	} else { 
		for (int i=1; i<dom.nx-1; i++) {
			x=i*dom.dx;			
			for (int j=1; j<dom.ny-1; j++) {
				y=j*dom.dy;
				if (x>=0.5 and y>=0.5) {
					dom.W[i][j][0]=rinit[0];
					dom.W[i][j][1]=uinit[0];
					dom.W[i][j][2]=vinit[0];
					dom.W[i][j][3]=winit[0];
					dom.W[i][j][4]=pinit[0];
					if (dom.nvar==8) {
						dom.W[i][j][5]=Bxinit[0];
						dom.W[i][j][6]=Byinit[0];
						dom.W[i][j][7]=Bzinit[0];				 	  }
				} else if (x<0.5 and y>=0.5) {
					dom.W[i][j][0]=rinit[1];
					dom.W[i][j][1]=uinit[1];
					dom.W[i][j][2]=vinit[1];
					dom.W[i][j][3]=winit[1];
					dom.W[i][j][4]=pinit[1];
					if (dom.nvar==8) {
						dom.W[i][j][5]=Bxinit[1];
						dom.W[i][j][6]=Byinit[1];
						dom.W[i][j][7]=Bzinit[1];
					}
				} else if (x<0.5 and y<0.5) {
					dom.W[i][j][0]=rinit[2];
					dom.W[i][j][1]=uinit[2];
					dom.W[i][j][2]=vinit[2];
					dom.W[i][j][3]=winit[2];
					dom.W[i][j][4]=pinit[2];
					if (dom.nvar==8) {
						dom.W[i][j][5]=Bxinit[2];
						dom.W[i][j][6]=Byinit[2];
						dom.W[i][j][7]=Bzinit[2];
					}
				} else if (x>=0.5 and y<0.5) {
					dom.W[i][j][0]=rinit[3];
					dom.W[i][j][1]=uinit[3];
					dom.W[i][j][2]=vinit[3];
					dom.W[i][j][3]=winit[3];
					dom.W[i][j][4]=pinit[3];
					if (dom.nvar==8) {
						dom.W[i][j][5]=Bxinit[3];
						dom.W[i][j][6]=Byinit[3];
						dom.W[i][j][7]=Bzinit[3];	
					}
				}	
			}
		}
		
	}

	// Assign values from initial condition to class parameters
	dom.setParam(gamma,tEnd);
}
