#ifndef LBM_H
#define LBM_H

#include <vector> // vector containers
#include <cmath> // mathematical library
#include <iostream> // for the use of 'cout'
#include <fstream> // file streams
#include <sstream> // string streams
#include <cstdlib> // standard library

#include <armadillo>

using namespace arma;

#include "particle.H"

class LBM
{
protected:
	
	//- Density
	mat rho;
	double rho0;
	double uLid;
	double omega;
	
	//- Velocities
	mat u;
	mat v;
	
	//double** cu;
	
	//- Forces
//	double** fx;
//	double** fy;
	
	//- Populations
	cube fOut;
	cube fEq;
	cube fIn;
	
	//- Number of nodes
	int nx;
	int ny;
	
	//- Model parameters 2DQ9
	const int d = 2;
	const int q = 9;
	

	mat cx;
	mat cy;
	mat t;

	//const int opp[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

public:
	
	//- Constructor
	LBM();

	LBM
	(
		const double& rho0,
		const double& uLid,
		const double& omega,
		const int& nx,
		const int& ny
	);
	
	void initialize();
	void computeMacros();
	void updateDirichletBC();
	void equilibrium();
	void updateBoundsToEq();
	void collision();
	void stream();
////	void collideAndSwap();
////	void streamBySwapping();


////	void computeU();
////	double popEq(int iX, int iY, int iF);
//	
	void writeFluidVTK(const double& time);
//	void writeParticleVTK(const particle& part, const double& time);
//	void writeData(const double& time);
//	
//	friend class cavityBC;
};

#endif
