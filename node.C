#include "node.H"

node::node
(
	double x,
	double y,
	double xRef,
	double yRef,
	double u,
	double v,
	double fx,
	double fy
)
:
	x_(x), y_(y), xRef_(xRef), yRef_(yRef), u_(u), v_(v), fx_(fx), fy_(fy)
{}
