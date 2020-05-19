#ifndef particle_H
#define particle_H

#include "node.H"

class particle
{
protected:
	
	//- Number of surface nodes
	int nNodes_;
	
	//- Object radius
	double radius_;
	
	//- Center node
	node center_;	

	//- List of nodes
	node* nodePtr_;

public:
	
	//- Constructor
	particle
	(
		int nNodes,
		double radius,
		node center,
		node* nodePtr 
	);
	
	int nNodes() const;
	node nodeXY(int n) const;
	node nodeXY(int n);
};

#endif
