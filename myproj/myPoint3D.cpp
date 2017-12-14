#include "stdafx.h"
#include "mypoint3d.h"
#include "myvector3d.h"
#include <iostream>

#define min(a,b)            (((a) < (b)) ? (a) : (b))
#define max(a,b)            (((a) > (b)) ? (a) : (b))

myPoint3D::myPoint3D() {}

myPoint3D::myPoint3D(double x, double y, double z)
{
    X = x;
    Y = y;
    Z = z; 
}

myPoint3D myPoint3D::operator+(myVector3D & v1)
{
	return myPoint3D(X+v1.dX, Y+v1.dY, Z+v1.dZ);
}

myPoint3D & myPoint3D::operator+=(myPoint3D & p1)
{
	X += p1.X;
	Y += p1.Y;
	Z += p1.Z;
	return *this;
}

myPoint3D & myPoint3D::operator+=(myVector3D & v1)
{
	X += v1.dX;
	Y += v1.dY;
	Z += v1.dZ;
	return *this;
}

myVector3D & myPoint3D::operator-(myPoint3D & p1)
{
	return myVector3D( X-p1.X, Y-p1.Y, Z-p1.Z );
}

double myPoint3D::dist(myPoint3D p1)
{
	double x, y, z;
	x = pow((this->X - p1.X),2.0);
	y = pow((this->Y - p1.Y),2.0);
	z = pow((this->Z - p1.Z),2.0);

	return sqrt(x+y+z);
}

void myPoint3D::rotate(myVector3D & lp, double theta)
{
	myVector3D tmp(X, Y, Z);
	tmp.rotate(lp, theta);
	X = tmp.dX; Y = tmp.dY; Z = tmp.dZ;
}

void myPoint3D::print(char *s)
{
	  std::cout << s << X << ", " << Y << ", " << Z << "\n";
}

double myPoint3D::dist(myPoint3D *p1, myPoint3D *p2)
{
	double l2 = pow(p1-p2, 2.0);
	if(l2 == 0.0) return dist(*p1);

	//the line segment being p1 + t(p2-p1), we have to determine t to find projection 
	double t = max(0, min(1,((this-p1) * (p2-p1))/l2)); 

	//point projected to line segment
	myPoint3D *projection = p1; 
	myVector3D *tmp = new myVector3D;
	*tmp = (*p1-*p2) * t;
	*projection += *tmp;

	//distance from point to projection point
	return dist(this, projection);
}

double myPoint3D::dist(myPoint3D *p1, myPoint3D *p2, myPoint3D *p3)
{
	return 0.0;
}