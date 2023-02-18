// To allow different types of flux splitting

#include "defs.hpp"

// List of possible flux splitting methods
void GLF(real &FL,real &FR,real flux4L,real Q4L,real lambda,real flux4R,real Q4R) {
	// Lax-Friedrichs flux splitting (can be Global or Local)
	FR=0.5*(flux4R+lambda*Q4R);
	FL=0.5*(flux4L-lambda*Q4L);
}
