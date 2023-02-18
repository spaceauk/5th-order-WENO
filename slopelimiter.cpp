#include<exception>

#include "defs.hpp"


real slopelimiter(string limiter,real r) {
	real phi=0.;
	if (limiter=="1st") {
		phi=0.;
	} else if (limiter=="VanLeer") {
		phi=(abs(r)+r)/(1.+abs(r));
	} else if (limiter=="Superbee") {
		real intm0=2.*r;
		real intm1=MIN(intm0,1.);
		intm0=MIN(r,2.);
		real intm2=MAX(intm0,intm1);
		phi=MAX(0.,intm2);
	} else if (limiter=="MC") {
		real intm0=2.*r,intm1=0.5*(1.+r);
		real intm2=MIN(intm0,intm1);
		intm0=MIN(intm2,2.); 
		phi=MAX(0.,intm0);
	} else if (limiter=="Koren") {
		real intm0=(1.+2.*r)/3.;
		real intm1=MIN(intm0,2.);
		intm0=2.*r;
		real intm2=MIN(intm0,intm1);
		phi=MAX(0.,intm2);
	} else {
		cout<<"Error as invalid slope limiter chosen!"<<endl;
		throw exception();
	}

	return phi;
}


