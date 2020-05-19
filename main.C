#include <vector> // vector containers
#include <cmath> // mathematical library
#include <iostream> // for the use of 'cout'
#include <fstream> // file streams
#include <sstream> // string streams
#include <cstdlib> // standard library

#include "LBM.H"
//#include "cavityBC.H"

using namespace std;

int main(int argc, char* argv[])
{
	int time = 0;
	int maxT = 1001;
	int tPlot = 1000;
	
	int nx = 256 + 1;
	int ny = 256 + 1;
		
	const double uLid = 0.05;
	const double rho0 = 1.0;
	const double Re = 10420;
	const double nu = uLid*ny/Re;
	const double omega = 1./(3.*nu + 0.5);
	
	LBM lbm(rho0, uLid, omega, nx, ny);
	lbm.initialize();
	
	while(time <= maxT)
	{
		lbm.computeMacros();
		lbm.updateDirichletBC();
		lbm.equilibrium();
		lbm.updateBoundsToEq();
		lbm.collision();
		lbm.stream();

		if(time % tPlot == 0)
		{
			lbm.writeFluidVTK(time);
		}

		time ++;
		cout<< "Time: " << time << endl;
	}
	std::cout << arma::arma_version::as_string() << std::endl;
	
	return 0;
}
