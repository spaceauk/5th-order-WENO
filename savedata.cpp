#include<fstream>  // input & output

#include"defs.hpp"

void savedata(meshblock &dom) {
	real intm=0.0;

	ofstream Wden;
	Wden.open("savefile/density");
	for (int j=0; j<dom.ny; j++) {
		for (int i=0; i<dom.nx; i++) {
			intm=dom.W[i][j][0];
			Wden<<intm<<" ";
		}
		Wden<<endl;
	}
	Wden.close();

	if (dom.IC=="BWy") {
		ofstream Wdeny;
		Wdeny.open("savefile/denBWy");
		for (int j=0; j<dom.ny; j++) {
			intm=dom.W[dom.nx/2][j][0];
			Wdeny<<intm<<" ";
		}
		Wdeny.close();
	}


	ofstream Wu;
	Wu.open("savefile/u_velocity");
	for (int j=0; j<dom.ny; j++) {
		for (int i=0; i<dom.nx; i++) {
			intm=dom.W[i][j][1];
			Wu<<intm<<" ";
		}
		Wu<<endl;
	}
	Wu.close();

	ofstream Wv;
	Wv.open("savefile/v_velocity");
	for (int j=0; j<dom.ny; j++) {
		for (int i=0; i<dom.nx; i++) {
			intm=dom.W[i][j][2];
			Wv<<intm<<" ";
		}
		Wv<<endl;
	}
	Wv.close();

	ofstream Wp;
	Wp.open("savefile/pressure");
	for (int j=0; j<dom.ny; j++) {
		for (int i=0; i<dom.nx; i++) {
			intm=dom.W[i][j][4];
			Wp<<intm<<" ";
		}
		Wp<<endl;
	}
	Wp.close();

	if (dom.nvar>5) {
		cout<<"Saving magnetic field too..."<<endl;
		ofstream WBx;
		WBx.open("savefile/Bx");
		for (int j=0; j<dom.ny; j++) {
			for (int i=0; i<dom.nx; i++) {
				intm=dom.W[i][j][5];
				WBx<<intm<<" ";
			}
			WBx<<endl;
		}
		WBx.close();
		ofstream WBy;
		WBy.open("savefile/By");
		for (int j=0; j<dom.ny; j++) {
			for (int i=0; i<dom.nx; i++) {
				intm=dom.W[i][j][6];
				WBy<<intm<<" ";
			}
			WBy<<endl;
		}
		WBy.close();
	}

}


void savearray(meshblock dom,real*** array, string arrname) {

	string ARRname="savefile/"+arrname;	
	real intm;
	ofstream arr;
	arr.open(ARRname);
	for (int k=0; k<dom.nvar; k++) {
		arr<<"At ivar="<<k<<endl;
		for (int j=0; j<dom.ny; j++) {
			for (int i=0; i<dom.nx; i++) {
				intm=array[i][j][k];
				arr<<intm<<" ";
			}
			arr<<endl;
		}
		arr<<"\n \n";
	}
	arr.close();

}
