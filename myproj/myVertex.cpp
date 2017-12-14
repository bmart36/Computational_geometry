#include "StdAfx.h"
#include "myVertex.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myFace.h"

myVertex::myVertex(void)
{
	point = new myPoint3D();
	originof = NULL;
	normal = new myVector3D(0,0,0);
}

myVertex::~myVertex(void)
{
}

void myVertex::computeNormal()
{
	myHalfedge *tmp = this->originof;
	int d = 0;

	do{
		normal->dX += tmp->adjacent_face->normal->dX;
		normal->dY += tmp->adjacent_face->normal->dY;
		normal->dZ += tmp->adjacent_face->normal->dZ;

		tmp = tmp->twin->next;
		d++;
	}while(tmp!=this->originof);

	*(this->normal) = *(normal)*(1.0/d);

	this->normal->normalize();
}