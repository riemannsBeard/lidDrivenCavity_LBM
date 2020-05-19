#ifndef node_H
#define node_H


class node
{
protected:
	
	//- Current position
	double x_;
	double y_;
	
	//- Reference position
	double xRef_;
	double yRef_;
	
	//- Node velocity
	double u_;
	double v_;
	
	//- Node force
	double fx_;
	double fy_;

public:
	
	//- Constructor
	node
	(
		double x,
		double y,
		double xRef,
		double yRef,
		double u,
		double v,
		double fx,
		double fy
	);
	
	inline double x()
	{
		return x_;
	}

	inline double y()
	{
		return y_;
	}

	friend class particle;

};

#endif
