#include "cavityBC.H"

cavityBC::cavityBC(LBM& lb)
:
lb_(lb)
{}

void cavityBC::update()
{
	//- Top BC
	for(int iX=1; iX<nx-1; ++iX)
	{
		lb_.u[iX][ny-1] = lb_.uLid;
		lb_.v[iX][ny-1] = 0.0;
		lb_.rho[iX][ny-1] = lb_.rho0;
	}
	
	//- Bottom BC
	for(int iX=0; iX<nx; ++iX)
	{
		lb_.u[iX][0] = 0.0;
		lb_.v[iX][0] = 0.0;
		lb_.rho[iX][0] = lb_.rho0;
	}
	
	//- Left BC
	for(int iY=0; iY<ny; ++iY)
	{
		lb_.u[0][iY] = 0.0;
		lb_.v[0][iY] = 0.0;
		lb_.rho[0][iY] = lb_.rho0;
	}
	
	//- Right BC
	for(int iY=0; iY<ny; ++iY)
	{

		lb_.u[nx-1][iY] = 0.0;
		lb_.v[nx-1][iY] = 0.0;
		lb_.rho[nx-1][iY] = lb_.rho0;
	}
}

void cavityBC::equilibrium()
{
	//- Top BC
	for(int iX=1; iX<nx-1; ++iX)
	{
		for(int iF=0; iF<q; ++iF)
		{
			lb_.fIn[iX][ny-1][iF] = lb_.fEq[iX][ny-1][iF];
		}
	}
	
	//- Bottom BC
	for(int iX=0; iX<nx; ++iX)
	{
		for(int iF=0; iF<q; ++iF)
		{
			lb_.fIn[iX][0][iF] = lb_.fEq[iX][0][iF];
		}
	}
	
	//- Left BC
	for(int iY=0; iY<ny; ++iY)
	{
		for(int iF=0; iF<q; ++iF)
		{
			lb_.fIn[0][iY][iF] = lb_.fEq[0][iY][iF];
		}
	}
	
	//- Right BC
	for(int iY=0; iY<ny; ++iY)
	{
		for(int iF=0; iF<q; ++iF)
		{
			lb_.fIn[nx-1][iY][iF] = lb_.fEq[nx-1][iY][iF];
		}
	}
}
