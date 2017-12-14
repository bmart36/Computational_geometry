#include "StdAfx.h"
#include "myHalfedge.h"

myHalfedge::myHalfedge(void)
{
	source = NULL; 
	adjacent_face = NULL; 
	next = NULL;  
	prev = NULL;  
	twin = NULL;  
}


myHalfedge::~myHalfedge(void)
{
}
