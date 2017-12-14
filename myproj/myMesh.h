#pragma once
#include "myFace.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <vector>
#include <string>
#include <map>

class myMesh
{
public:
	std::vector<myVertex *> vertices;
	std::vector<myHalfedge *> halfedges;
	std::vector<myFace *> faces;
	std::string name;

	void checkMesh();
	void readFile(std::string filename);
	void computeNormals();
	void normalize();
	void triangulate();
	void inflateMesh(double dist);
	void smoothenMesh(double dist);
	void splitEdge(myHalfedge *e, myVertex *v);
	void splitFace(myFace *f, myVertex *v);
	void subdivisionCatmullClark();

	myPoint3D edge_point(myHalfedge *e);
	myPoint3D face_point(myFace *f);

	myMesh(void);
	~myMesh(void);
};

