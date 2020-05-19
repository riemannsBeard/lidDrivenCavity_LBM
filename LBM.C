#include <iostream> // for the use of 'cout'
#include <fstream> // file streams
#include <sstream> // string streams
#include <cstdlib> // standard library

#include "LBM.H"

using namespace std;

LBM::LBM
(
	const double& rho0,
	const double& uLid,
	const double& omega,
	const int& nx,
	const int& ny
)
:
	rho0(rho0),
	uLid(uLid),
	omega(omega),
	nx(nx),
	ny(ny)
{
	t =
	{
		4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.
	};
	
	cx = {0, 1, 0, -1, 0, 1, -1, -1, 1};
	cy = {0, 0, 1, 0, -1, 1, 1, -1, -1};
}

void LBM::computeMacros()
{
	for(int i=0; i<ny; i++)
	{	
		for(int j=0; j<nx; j++)
		{
			u(i, j) = 0.0;
			v(i, j) = 0.0;
			rho(i, j) = 0.0;
				
			for(int k=0; k<q; k++)
			{
				rho(i, j) += fIn(k, i, j);
				
				u(i, j) += fIn(k, i, j)*cx(k);
				v(i, j) += fIn(k, i, j)*cy(k);
			}
		}
	}
	u = u/rho;
	v = v/rho;
	
//	cout <<"rho =\n" << rho << endl;
//	cout << "cx = \n" << cx << endl;
//	
//	cout << "u = \n" << u << endl;
//	cout << "v = \n" << v << endl;
}


void LBM::initialize()
{

  	//- Create folders, delete data file
  	//- Make sure that the VTK folders exist.
  	//- Old file data.dat is deleted, if existing.
	int ignore;
	ignore = system("mkdir -p vtkFluid");
	//ignore = system("mkdir -p vtkParticle");
	//ignore = system("rm -f data.dat");
	
	//- Allocate memory
	rho = mat(ny, nx);
	u = mat(ny, nx);
	v = mat(ny, nx);
	
	fOut = cube(q, ny, nx);
	fIn = cube(q, ny, nx);
	fEq = cube(q, ny, nx);

	for(int i=0; i<ny; i++)
	{
		for(int j=0; j<nx; j++)
		{
			for(int k=0; k<q; k++)
			{
				fIn(k, i, j) = t(k);
			}
		}
	}

	//cout << "fIn = \n" << fIn;
}

void LBM::equilibrium()
{
	for(int i=0; i<ny; i++)
	{
		for(int j=0; j<nx; j++)
		{
			double u2 = u(i,j)*u(i,j);
			double v2 = v(i,j)*v(i,j);
			
			for(int k=0; k<q; k++)
			{
				double cu = 3*(cx(k)*u(i,j) + cy(k)*v(i,j));
//				
			 	fEq(k, i, j) = t(k)*rho(i,j)*(1. + cu +
 					0.5*(cu*cu) - 1.5*(u2 + v2));
//
			}
		}
			
 		//cout << "fEq =\n" << fEq << endl;
	}
}

void LBM::updateDirichletBC()
{
	//- Top BC
	u.row(nx-1) = uLid*ones<mat>(1, nx);
	v.row(nx-1) = zeros<mat>(1, nx);
	rho.row(nx-1) = ones<mat>(1, nx);
		
	//- Bottom BC
	u.row(0) = zeros<mat>(1, nx);
	v.row(0) = zeros<mat>(1, nx);
	rho.row(0) = ones<mat>(1, nx);
	
	//- Left BC
	u.col(0) = zeros<mat>(ny, 1);
	v.col(0) = zeros<mat>(ny, 1);
	rho.col(0) = ones<mat>(ny, 1);
	
	//- Right BC
	u.col(nx-1) = zeros<mat>(ny, 1);
	v.col(nx-1) = zeros<mat>(ny, 1);
	rho.col(nx-1) = ones<mat>(ny, 1);
	
//	cout << "rho =\n" << rho << endl;
}

void LBM::updateBoundsToEq()
{
// BOUNDARY CONDITION: RESET BOUNDARY NODES TO EQUILIBRIUM

	//- Top BC
	fIn(span::all, span(ny-1), span(0, nx-1)) = 
		fEq(span::all, span(ny-1), span(0, nx-1));
		
	//- Bottom BC
	fIn(span::all, span(0), span(0, nx-1)) = 
		fEq(span::all, span(0), span(0, nx-1));

	//- Left BC
	fIn(span::all, span(0, ny-1), span(0)) = 
		fEq(span::all, span(0, ny-1), span(0));	
	

	//- Right BC
	fIn(span::all, span(0, ny-1), span(nx-1)) = 
		fEq(span::all, span(0, ny-1), span(nx-1));

}

void LBM::collision()
{
	for(int k=0; k<9; k++)
	{
		fOut(span(k), span::all, span::all) = 
			fIn(span(k), span::all, span::all) - omega*
			(
				fIn(span(k), span::all, span::all) - 
				fEq(span(k), span::all, span::all)
			);
	}
}

void LBM::stream()
{
	for(int i=0; i<ny; i++)
	{	
		for(int j=0; j<nx; j++)
		{
			for(int k=0; k<q; k++)	
			{	
				int h = i - cy(k);
				if(h<0)
				{
					h = ny-1;
				}
				if(h>ny-1)
				{
					h = 0;
				}
			
				int l = j - cx(k);
				if(l<0)
				{
					l = nx-1;
				}
				if(l>ny-1)
				{
					l = 0;
				}
									
				fIn(k, i, j) = fOut(k, h, l);
			}		
		}
	}
}


void LBM::writeFluidVTK(const double& time)
{
  	//- Create filename
  	stringstream outputFilename;
  	outputFilename << "vtkFluid/fluid_t" << time << ".vtk";
  	ofstream outputFile;

  	//- Open file

  	outputFile.open(outputFilename.str().c_str());

  	//- Write VTK header

  	outputFile << "# vtk DataFile Version 3.0\n";
  	outputFile << "fluid_state\n";
  	outputFile << "ASCII\n";
  	outputFile << "DATASET RECTILINEAR_GRID\n";
  	outputFile << "DIMENSIONS " << nx << " " << ny << " 1" << "\n";
  	outputFile << "X_COORDINATES " << nx << " float\n";

  	for(int i=0; i<nx; i++)
  	{
    	outputFile << i << "\t";
  	}

  	outputFile << "\n";
  	outputFile << "Y_COORDINATES " << ny << " float\n";

  	for(int i=0; i<ny; i++)
  	{
    	outputFile << i << "\t";
  	}

  	outputFile << "\n";
  	outputFile << "Z_COORDINATES 1 float\n";
  	outputFile << 0 << "\n";
  	
  	outputFile << "POINT_DATA " << nx*ny << "\n";

	//- Write density difference

  	outputFile << "SCALARS rho double\n";
	outputFile << "LOOKUP_TABLE default\n";

	for(int i=0; i<ny; i++)
  	{
		for(int j=0; j<nx; j++)
		{
			outputFile << rho(i, j) << "\n";
		}
	}
	
  	//- Write velocity vector
  	outputFile << "VECTORS U double\n";

	for(int i=0; i<ny; i++)
  	{
		for(int j=0; j<nx; j++)
		{
			outputFile << u(i, j) << "\t" << v(i, j) << "\t 0\n";
    	}
  	}

  	//- Close file
	outputFile.close();
}

//void LBM::writeParticleVTK(const particle& part, const double& time)
//{
//  	//- Create filename
//  	stringstream outputFilename;
//  	outputFilename << "vtkParticle/particle_t" << time << ".vtk";
//  	ofstream outputFile;

//  	//- Open file
//  	outputFile.open(outputFilename.str().c_str());

//  	//- Write VTK header
//  	outputFile << "# vtk DataFile Version 3.0\n";
//  	outputFile << "particle_state\n";
//  	outputFile << "ASCII\n";
//  	outputFile << "DATASET POLYDATA\n";

//  	//- Write node positions
//  	outputFile << "POINTS " << part.nNodes() << " float\n";
//  	
//  	//- Store number of nodes
//  	const int k = part.nNodes();

//  	for(int n=0; n<k; ++n)
//  	{
//    	outputFile << part.nodeXY(n).x() << " " << part.nodeXY(n).y() << " 0\n";
//  	}

//  	//- Write lines between neighboring nodes
//  	outputFile << "LINES " << k << " " << 3*k << "\n";

//  	for(int n=0; n<k; ++n)
//  	{
//    	outputFile << "2 " << n << " " << (n + 1)%k << "\n";
//  	}

//  	//- Write vertices

//  	outputFile << "VERTICES 1 " << k << "\n";
//  	outputFile << k << " ";

//  	for(int n=0; n<k; ++n)
//  	{
//    	outputFile << n << " ";
//  	}

//  	//- Close file

//  	outputFile.close();
//}



//void LBM::writeData(const double& time)
//{
//	//- Create filename
////	string outputFilename("data.dat");
////  	ofstream outputFile;
//  	
//  	//- Open file
////	outputFile.open(outputFilename.c_str(), fstream::app);

//  	//- Compute quantities
//	//double force_tot_x = 0;
//  	//double force_tot_y = 0;
//  	//double vel_center_x = 0;
//  	//double vel_center_y = 0;

////  	for(int i = 0; i < particle.num_nodes; ++i)
////  	{
////    	force_tot_x += particle.node[i].force_x;
////    	force_tot_y += particle.node[i].force_y;
////    	vel_center_x += particle.node[i].vel_x;
////    	vel_center_y += particle.node[i].vel_y;
////  	}

//  	//- Write data
////  	output_file << time << " "; // time step
////  	output_file << force_tot_x << " "; // drag force
////  	output_file << force_tot_y << " "; // lift force
////  	output_file << particle.center.x << " "; // center position (x-component)
////  	output_file << particle.center.y << " "; // center position (y-component)
////  	output_file << vel_center_x << " "; // center velocity (x-component)
////  	output_file << vel_center_y << "\n"; // center velocity (y-component)
//}

