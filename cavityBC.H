#ifndef cavityBC_H
#define cavityBC_H

#include "LBM.H"

class cavityBC
{
private:

	LBM& lb_;
	const int& nx = lb_.nx;
	const int& ny = lb_.ny;
	const int& q = lb_.q;

public:


	//- Constructor
	cavityBC(LBM& lb);
	
	void update();
	void equilibrium();
	
};
#endif
