#include "particle.H"
#include "node.H"
#include <cmath>

particle::particle
(
	int nNodes,
	double radius,
	node center,
	node* nodePtr
)
:
	nNodes_(nNodes),
	radius_(radius),
	center_(center),
	nodePtr_(nodePtr)
{
	const double pi = 3.1415926535897;
	
	for(int n=0; n<nNodes_; ++n)
	{
		nodePtr_[n].x_ = center_.x() +
			+ radius_*sin(2.*pi*static_cast<double>(n)/nNodes_);

		nodePtr_[n].y_ = center_.y() +
			+ radius_*cos(2.*pi*static_cast<double>(n)/nNodes_);

		nodePtr_[n].xRef_ = center_.x() +
			+ radius_*sin(2.*pi*static_cast<double>(n)/nNodes_);

		nodePtr_[n].yRef_ = center_.y() +
			+ radius_*cos(2.*pi*static_cast<double>(n)/nNodes_);	
	}
}

int particle::nNodes() const
{
	return nNodes_;
}

node particle::nodeXY(int n) const
{
	return nodePtr_[n];
}

node particle::nodeXY(int n)
{
	return nodePtr_[n];
}


