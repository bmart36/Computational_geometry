#include "StdAfx.h"
#include "myMesh.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <GL/glew.h>
#include "myvector3d.h"

using namespace std;

myMesh::myMesh(void)
{
	vertices.clear();
	faces.clear();
	halfedges.clear();
}


myMesh::~myMesh(void)
{
}

void myMesh::checkMesh()
{
	for(auto &e:halfedges)
	{
		if(e->twin->twin != e)
			cout << "twin problem" << endl;
		if(e->prev->next != e)
			cout << "prev next problem" << e << endl; //7
		if (e->next == NULL || e->prev == NULL || e->twin == NULL)
			cout << "NULL value" << endl;
		if(e->next->adjacent_face != e->adjacent_face)
			cout << "adjecent face problem" << e << endl; //2
		if(e->twin->next->source != e->source)
			cout << "source problem" << e << endl; //1
	}
	int count = 0;
	for(auto &f:faces)
	{
		myHalfedge *tmp = f->adjacent_halfedge;
		do{
			count++;
			tmp = tmp->next;
		}while(tmp != f->adjacent_halfedge);

		if(count < 3)
			cout << "face has less than 3 edges" << endl;
	}
}

void myMesh::readFile(std::string filename)
{
	string s, t;
	string tmp;
	map<pair<int, int>, myHalfedge *> twin_map;

	ifstream fin(filename);
	if (!fin.is_open()) throw std::runtime_error( "Unable to open file!");

	while ( getline(fin, s) )
	{
		stringstream myline(s);
		myline >> t;
		if (t == "v")
		{
			myline >> tmp;
			double x=stof(tmp.substr(0, tmp.find("/")));

			myline >> tmp;
			double y=stof(tmp.substr(0, tmp.find("/")));

			myline >> tmp;
			double z=stof(tmp.substr(0, tmp.find("/")));

			myVertex *tmpV = new myVertex();
			tmpV->point = new myPoint3D(x, y, z);
			vertices.push_back(tmpV);

		}
		if (t == "f")
		{

			vector<int> v;
			vector<myHalfedge *> e;

			while (myline >> tmp) {
				v.push_back(stoi(tmp.substr(0, tmp.find("/")))-1);
			}

			if(v.size()==0)
				return;

			for(int i=0; i < v.size();i++){
				e.push_back(new myHalfedge);
			}

			myFace *tmpF= new myFace;
			tmpF->adjacent_halfedge = e[0];
			faces.push_back(tmpF);

			for(int i=0; i < v.size();i++){
				vertices[v[i]]->originof = e[i];
				e[i]->source = vertices[v[i]];
				e[i]->adjacent_face = tmpF;
				e[i]->next = e[(i+1)%e.size()];
				e[i]->prev = e[(i-1+e.size())%e.size()];
				halfedges.push_back(e[i]);
				map<pair<int, int>, myHalfedge *>::iterator it= twin_map.find(make_pair(v[(i+1)%e.size()], v[i]));
				if(it == twin_map.end()){
					twin_map[make_pair(v[i],v[(i+1)%e.size()])] = e[i];
				}
				else{
					myHalfedge *tmpH = it->second;
					e[i]->twin = tmpH;
					tmpH->twin = e[i];
				}
			} 
		}
	}
	checkMesh();
}

void myMesh::triangulate()
{
	vector<myFace *> facesT;
	myHalfedge *adjE;
	myHalfedge *tmpE;

	//Iterate through all faces in the structure
	for(int i=0; i<faces.size(); i++){
		int nbEdges = 0;
		adjE = faces[i]->adjacent_halfedge;
		tmpE = adjE;

		//Count amount of halfedges in the face
		do{
			tmpE = tmpE->next;
			nbEdges++;
		}while(tmpE!=adjE);

		//If face hasn't been triangulated
		if(nbEdges > 3){
			facesT.push_back(faces[i]);
		}
	}

	while(facesT.size()!=0){
		myFace *tmpF = facesT.back();
		adjE = tmpF->adjacent_halfedge;
		tmpE = adjE->next;
		int nbEdges = 0;

		//CReation new face
		myFace *newF = new myFace();

		//Create new halfedges
		myHalfedge *newE = new myHalfedge();
		myHalfedge *newT = new myHalfedge();
		newE->source = tmpE->next->source;
		newE->adjacent_face = newF;
		newE->next = adjE;
		newE->prev = tmpE;
		newE->twin = newT;

		//Modify faces adjecent edge
		tmpF->adjacent_halfedge = newT;

		//Modify face's adjacent edge
		newF->adjacent_halfedge = newE;

		//New halfedge's twin
		newT->source = adjE->source;
		newT->adjacent_face = tmpF;
		newT->next = tmpE->next;
		newT->prev = adjE->prev;
		newT->twin = newE;

		//Link prevoius and next to new twin
		tmpE->next->prev = newT;
		adjE->prev->next = newT;

		//Link previous and next to new halfedge and face
		adjE->prev = newE;
		tmpE->next = newE;
		adjE->adjacent_face = newF;
		tmpE->adjacent_face = newF;

		halfedges.push_back(newE);
		faces.push_back(newF);
		halfedges.push_back(newT);
		
		myHalfedge *e = tmpF->adjacent_halfedge;
		do{
			e = e->next;
			nbEdges++;
		}while(e!=tmpF->adjacent_halfedge);

		//If face has been triangulated
		if(nbEdges <= 3){
			facesT.pop_back();
		}
	}		
}

void myMesh::computeNormals()
{
	for (vector<myFace *>::iterator it = faces.begin(); it != faces.end(); it++)
		(*it)->computeNormal();

	for (vector<myVertex *>::iterator it = vertices.begin(); it != vertices.end(); it++)
		(*it)->computeNormal();
}

void myMesh::normalize()
{
	int i;
	int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;

	for (i=0;i<vertices.size();i++) {
		if (vertices[i]->point->X < vertices[tmpxmin]->point->X) tmpxmin = i;
		if (vertices[i]->point->X > vertices[tmpxmax]->point->X) tmpxmax = i;

		if (vertices[i]->point->Y < vertices[tmpymin]->point->Y) tmpymin = i;
		if (vertices[i]->point->Y > vertices[tmpymax]->point->Y) tmpymax = i;

		if (vertices[i]->point->Z < vertices[tmpzmin]->point->Z) tmpzmin = i;
		if (vertices[i]->point->Z > vertices[tmpzmax]->point->Z) tmpzmax = i;
	}

	double xmin = vertices[tmpxmin]->point->X, xmax = vertices[tmpxmax]->point->X, 
		ymin = vertices[tmpymin]->point->Y, ymax = vertices[tmpymax]->point->Y, 
		zmin = vertices[tmpzmin]->point->Z, zmax = vertices[tmpzmax]->point->Z;

	double scale =  (xmax-xmin) > (ymax-ymin) ? (xmax-xmin) : (ymax-ymin);
	scale = scale > (zmax-zmin) ? scale : (zmax-zmin);

	for (i=0;i<vertices.size();i++) {
		vertices[i]->point->X -= (xmax+xmin)/2;
		vertices[i]->point->Y -= (ymax+ymin)/2;
		vertices[i]->point->Z -= (zmax+zmin)/2;

		vertices[i]->point->X /= scale;
		vertices[i]->point->Y /= scale;
		vertices[i]->point->Z /= scale;
	}
}

void myMesh::inflateMesh(double dist)
{
	myVertex *v;
	for(int i = 0; i<vertices.size(); i++){
		v = vertices[i];

		*(v->point) = *(v->point) + (*(v->normal))*dist;
	}
}

myPoint3D myMesh::face_point(myFace *f)
{
	myPoint3D centroid(0,0,0);
	myHalfedge *tmp = f->adjacent_halfedge;
	int nVertices = 0;

	do{
		centroid += *tmp->source->point;
		nVertices++;
		tmp = tmp->next;
	}while(tmp != f->adjacent_halfedge);

	centroid.X = centroid.X/nVertices;
	centroid.Y = centroid.Y/nVertices;
	centroid.Z = centroid.Z/nVertices;

	return centroid;
}

myPoint3D myMesh::edge_point(myHalfedge *e)
{
	myPoint3D centroid(0,0,0);
	myHalfedge *te = e->twin;

	myPoint3D *face_p_e = new myPoint3D;
	myPoint3D *face_p_te = new myPoint3D;
	*face_p_e = face_point(e->adjacent_face);
	*face_p_te = face_point(te->adjacent_face);

	centroid = *e->source->point;
	centroid += *te->source->point;
	centroid += *face_p_e;
	centroid += *face_p_te;

	centroid.X = centroid.X/4;
	centroid.Y = centroid.Y/4;
	centroid.Z = centroid.Z/4;

	return centroid;
}

//void myMesh::smoothenMesh(double dist)
//{
//	myVertex *v;
//	myPoint3D *x;
//
//	for(int i = 0; i<vertices.size(); i++){
//		v = vertices[i];
//
//		x = averageNeighbours(v);
//		*(v->point) = (*(v->point))*(1-dist) + (*x)*dist;
//	}
//}

void myMesh::splitEdge(myHalfedge *e, myVertex *v)
{
	myHalfedge *te = e->twin;
	myHalfedge *en = new myHalfedge();
	myHalfedge *ten = new myHalfedge();
	map<std::pair<int, int>, myHalfedge *>::iterator it;

	v->originof = en;

	en->next = e->next;
	en->prev = e;
	en->twin = ten;
	en->source = v;
	en->adjacent_face = e->adjacent_face;
	en->next->prev = en;

	ten->next = te;
	ten->prev = te->prev;
	ten->twin = en;
	ten->source = e->next->source;
	ten->adjacent_face = te->adjacent_face;
	ten->prev->next = ten;

	e->next = en;

	te->prev = ten;
	te->source = v;

	ten->source->originof = ten;

	e->adjacent_face->adjacent_halfedge = en;

	vertices.push_back(v);
	halfedges.push_back(en);
	halfedges.push_back(ten);
}

int faceVertices(myFace *f)
{
	int nVertices = 0;

	myHalfedge *tmp = f->adjacent_halfedge;
	do{
		nVertices++;
		tmp = tmp->next;
	}while(tmp != f->adjacent_halfedge);

	return nVertices;
}

void myMesh::splitFace(myFace *f, myVertex *v)
{
	myHalfedge *s = f->adjacent_halfedge;
	vector<myHalfedge *> nes;
	vector<myHalfedge *> tnes;
	vector<myFace *> nfs;
	int n_2 = faceVertices(f)/2;
	int i;

	nes.resize(n_2);
	tnes.resize(n_2);
	nfs.resize(n_2);

	for(i=0; i<n_2; i++)
	{
		nes[i] = new myHalfedge();
		halfedges.push_back(nes[i]);
		tnes[i] = new myHalfedge();
		halfedges.push_back(tnes[i]);
	}

	v->originof = nes[0];
	vertices.push_back(v);

	nfs[0] = f;
	for(i=1; i<n_2; i++)
	{
		nfs[i] = new myFace();
		faces.push_back(nfs[i]);
	}

	for(i=0; i<n_2; i++)
	{
		myHalfedge *sCopy = s->next->next;
		int iplus1 = (i+1) % n_2;
		int iminus1 = (i-1+n_2) % n_2;

		nes[i]->next = s;
		nes[i]->prev = tnes[iplus1];
		nes[i]->twin = tnes[i];
		nes[i]->source = v;
		nes[i]->adjacent_face = nfs[i];

		tnes[i]->next = nes[iminus1];
		tnes[i]->prev = s->prev;
		tnes[i]->twin = nes[i];
		tnes[i]->source = s->source;
		tnes[i]->adjacent_face = nfs[iminus1];

		nfs[i]->adjacent_halfedge = nes[i];

		s->prev = nes[i];
		s->adjacent_face = nfs[i];

		s->next->next = tnes[iplus1];
		s->next->adjacent_face = nfs[i];

		s = sCopy;
	}
}

void getV(myVertex *v, map<myFace *, myVertex *> m)
{
	myPoint3D *v1 = v->point;
	myPoint3D *v2;
	myPoint3D q(0,0,0), r(0,0,0);
	myHalfedge *tmp = v->originof;
	myFace *f;
	int count = 0;

	do{
		//face points average Q
		f=tmp->adjacent_face;
		q += *(m.find(f)->second->point);

		//edge midpoints average R
		v2 = tmp->next->source->point;
		r.X += ((v1->X + v2->X)/2);
		r.Y += ((v1->Y + v2->Y)/2);
		r.Z += ((v1->Z + v2->Z)/2);

		count++;
		tmp = tmp->prev->twin;
	}while(tmp != v->originof);

	//Q
	q.X = q.X/count;
	q.Y = q.Y/count;
	q.Z = q.Z/count;

	//2*R
	r.X = (2*r.X)/count;
	r.Y = (2*r.Y)/count;
	r.Z = (2*r.Z)/count;

	//(n-3)*V
	v->point->X = (count-3)*v->point->X;
	v->point->Y = (count-3)*v->point->Y;
	v->point->Z = (count-3)*v->point->Z;

	v->point->X = (q.X + r.X + v->point->X)/count;
	v->point->Y = (q.Y + r.Y + v->point->Y)/count;
	v->point->Z = (q.Z + r.Z + v->point->Z)/count;
}

void myMesh::subdivisionCatmullClark()
{
	myVertex *v;
	myHalfedge *e;
	myPoint3D *centroid;
	map<myFace *, myVertex *> face_c_map;
	map<myHalfedge *, myVertex *> edge_c_map;
	vector<myFace *> copy_faces = faces;
	vector<myVertex *> copy_vertices = vertices;
	map<pair<myVertex*, myVertex *>, myHalfedge *> twin_map;

	for (auto &f:faces)
	{
		if(face_c_map.find(f)==face_c_map.end())
		{
			v = new myVertex();
			*(v->point) = face_point(f);
			face_c_map.insert(make_pair(f, v));
		}

		e = f->adjacent_halfedge;
		do{
			if(edge_c_map.find(e) == edge_c_map.end())
			{
				v = new myVertex();
				*(v->point) = edge_point(e);
				edge_c_map.insert(make_pair(e, v));
				edge_c_map.insert(make_pair(e->twin, v));
				twin_map.insert(make_pair(make_pair(e->source, e->twin->source), e));
			}
			e = e->next;
		}while(e != f->adjacent_halfedge);
	}

	for(auto &v:copy_vertices)
	{
		//we apply the formula to get the new vector v'= (1/n)*Q + ((2/n)*R) + (((n-3)/n)*v)
		getV(v, face_c_map);
	}

	for(auto const &e : twin_map){
		splitEdge(e.second, edge_c_map.find(e.second)->second);
	}

	for(auto &f:face_c_map)
	{
		splitFace(f.first, f.second);
	}

}