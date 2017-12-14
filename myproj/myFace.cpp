#include "StdAfx.h"
#include "myFace.h"
#include "myvector3d.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <GL/glew.h>

myFace::myFace(void)
{
	adjacent_halfedge = NULL;
	normal = new myVector3D();
}

myFace::~myFace(void)
{
}

void myFace::computeNormal()
{
	myHalfedge *e, *pe, *ne;
	e = this->adjacent_halfedge;
	pe = e->prev;
	ne = e->next;

	myVector3D v1, v2;
	v1 = *pe->source->point - *e->source->point;
	v2 = *ne->source->point - *e->source->point;

	*(this->normal) = v2.crossproduct(v1);
	this->normal->normalize();
}
