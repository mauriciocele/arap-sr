/** \file Mesh.cpp Definitions for the Mesh class.
    \author Liliya Kharevych
    \date March 2006
*/

#include "Mesh.h"
#include <iostream>
#include <fstream>
#include "util.h"
#include <stack>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>
#include <assert.h>

using namespace std;

Mesh::Mesh(void)
{
   numBoundaryVertices = 0;
}

Mesh::~Mesh(void)
{
   clear();
}

void Mesh::clear() {
   VertexIterator vIter = vertexIterator();
   FaceIterator fIter = faceIterator();
   HalfEdgeIterator eIter = halfEdgeIterator();

   while (!vIter.end()) {
      delete vIter.vertex();
      vIter++;
   }
   while (!fIter.end()) {
      delete fIter.face();
      fIter++;
   }
   while (!eIter.end()) {
      delete eIter.half_edge();
      eIter++;
   }
}

Vertex * Mesh::addVertex(const Vector3 & _p) {
   Vertex * v = new Vertex();
   v->p = _p;
   v->ID = (int)vertices.size();
   vertices.push_back(v);
   return v;
}

Face * Mesh::addFace(std::vector<int> faceVerts) {
   Face * f = new Face();
   f->ID = (int)faces.size();

   Edge * firstEdge = NULL;
   Edge * last = NULL;
   Edge * current = NULL;

   unsigned int i;

   for (i = 0; i < faceVerts.size()-1; i++) {
      current = addEdge(faceVerts[i], faceVerts[i+1]);
   
      check_error (!current, "Problem while parsing mesh faces.");

      current->face = f;

      if (last)
	 last->next = current;
      else
	 firstEdge = current;
      last = current;
   }

   current = addEdge (faceVerts[i], faceVerts[0]);
   check_error (!current, "Problem while parsing mesh faces.");
   
   current->face = f;

   last->next = current;
   current->next = firstEdge;

   f->edge = firstEdge;
   faces.push_back(f);
      
   return f;
}

Edge * Mesh::addEdge (int i, int j) {
   Edge eTmp;
   eTmp.vertex = vertices[i];

   Edge eTmpPair;
   eTmpPair.vertex = vertices[j];
   eTmp.pair = &eTmpPair;

   Mesh::EdgeSet::iterator eIn = edges.find(&eTmp);
   Edge * e;

   if (eIn != edges.end()){
      e = *eIn;
      if (e->face != NULL)
	 return NULL;
   }
   else {
      e = new Edge();
      Edge * pair = new Edge();

      e->vertex = vertices[i];
      pair->vertex = vertices[j];

      pair->pair = e;
      e->pair = pair;

      edges.insert(e);
      edges.insert(pair);

      pair->vertex->edge = pair;
   }   
   return e;
}

/** Cuts the mesh along single edge while maintaining
all the connectivity of the mesh 

Returns false if the mesh is boundary or NULL.

Note that after the mesh the cut is done, the sorting
order of edges in the set is not updated. This could be 
done after all the cuts are performed.
*/
bool Mesh::cutAlongEdge(Edge * forward) {
   if (!forward) return false;
   if (forward->isBoundary()) return false;

   Edge * backward = forward->pair;
   
   Edge * backwardPair = new Edge();
   Edge * forwardPair = new Edge();

   backwardPair->pair = backward;
   forwardPair->pair = forward;
   backwardPair->assignData(*backward);
   forwardPair->assignData(*forward);
   
   Vertex * v = forward->vertex;
   Vertex * endV = backward->vertex;

   if (v->isBoundary()) {
      Vertex * newV = addVertex(v->p);
      newV->assignData(*v);
      Vertex::EdgeAroundIterator aroundV = v->iterator(forward);
      do {
	 aroundV++;
	 aroundV.vertex() = newV;
      }
      while (aroundV.edge_out()->pair->face);

      forwardPair->next = aroundV.edge_out()->pair->next;
      aroundV.edge_out()->pair->next = backwardPair;

      backwardPair->vertex = newV;
      newV->edge = backwardPair;
   }
   else {
      forwardPair->next = backwardPair;
      backwardPair->vertex = v;
      v->edge = backwardPair;
   }

   if (endV->isBoundary()) {
      Vertex * newV = addVertex(endV->p);
      newV->assignData(*endV);
      Vertex::EdgeAroundIterator aroundV = endV->iterator(backward);
      do {
	 aroundV++;
	 aroundV.vertex() = newV;
      }
      while (aroundV.edge_out()->pair->face);

      backwardPair->next = aroundV.edge_out()->pair->next;
      aroundV.edge_out()->pair->next = forwardPair;

      forwardPair->vertex = newV;
      newV->edge = forwardPair;
   }
   else {
      backwardPair->next = forwardPair;
      forwardPair->vertex = endV;
      endV->edge = forwardPair;
   }

   forward->pair = forwardPair;
   backward->pair = backwardPair;

   return (addEdge(forwardPair) && addEdge(backwardPair));
}

/** Compute edge lengths based on vertex positions */
void Mesh::computeLengths() {
   int eIndex = this->numVertices();
   for (Mesh::EdgeIterator eIter = edgeIterator(); !eIter.end(); eIter++, eIndex++) {
      Edge * e = eIter.edge();
	  e->ID = eIndex;
	  e->pair->ID = eIndex;
      e->length = (e->vertex->p - e->next->vertex->p).abs();
      e->pair->length = e->length;
   }
}

static double arc_cosine ( double c )

//****************************************************************************80
//
//  Purpose:
//
//    ARC_COSINE computes the arc cosine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ACOS routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//    This routine truncates arguments outside the range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double C, the argument, the cosine of an angle.
//
//    Output, double ARC_COSINE, an angle whose cosine is C.
//
{
  double angle;
  double pi = 3.141592653589793;

  if ( c <= -1.0 )
  {
    angle = pi;
  }
  else if ( 1.0 <= c )
  {
    angle = 0.0;
  }
  else
  {
    angle = acos ( c );
  }
  return angle;
}

void Mesh::computeCirclePackingMetric()
{
	this->computeLengths();
	this->computeNormals();

	map< int, map<int, double> >		radiusPerFace;

	for( Mesh::FaceIterator fIter = this->faceIterator() ; !fIter.end() ; fIter++ )
	{
		Face				*face = fIter.face();
		map<int, double>	radius;

		for( Face::EdgeAroundIterator eIter = face->iterator() ; !eIter.end() ; eIter++ )
		{
			double		lij = eIter.edge()->length;
			double		ljk = eIter.edge()->next->length;
			double		lki = eIter.edge()->next->next->length;

			double		r = ( lki + lij - ljk ) / 2.0;

			map<int, double>::value_type	vertexPair( eIter.vertex()->ID, r );

			radius.insert( vertexPair );
		}

		map< int, map<int, double> >::value_type		facePair(face->ID, radius);

		radiusPerFace.insert( facePair );
	}

	for( Mesh::VertexIterator vIter = this->vertexIterator() ; !vIter.end() ; vIter++ )
	{
		Vertex::EdgeAroundIterator	around = vIter.vertex()->iterator();

		double	radius = 0.0;
		int		vertexId = vIter.vertex()->ID;

		for( ; !around.end(); around++ )
		{
			if(around.edge_out()->face)
			{
				int		faceId = around.edge_out()->face->ID;

				radius = std::max( radiusPerFace[ faceId ][ vertexId ], radius );
			}
		}

		vIter.vertex()->radius = radius;

		assert( vIter.vertex()->radius > 0.0 );
	}

	Mesh::EdgeIterator eIter = this->edgeIterator();
	for( ; !eIter.end() ; eIter++ )
	{
		double Lij = eIter.edge()->length;
		double ri = eIter.edge()->vertex->radius;
		double rj = eIter.edge()->next->vertex->radius;

		double cosPhi = (Lij*Lij - ri*ri - rj*rj)/( 2*ri*rj );
		double phi = arc_cosine( cosPhi );
		eIter.edge()->phiAngle = phi; //std::max( std::min( phi, M_PI / 2.0 ), 0.0 );
		eIter.edge()->pair->phiAngle = eIter.edge()->phiAngle;

		//assert( eIter.edge()->phiAngle <= M_PI / 2.0 );
	}
}

void Mesh::computeNormals()
{
	map< int, Vector3 >					normals;

	for( Mesh::FaceIterator fIter = this->faceIterator() ; !fIter.end() ; fIter++ )
	{
		Face		*face = fIter.face();
		Vector3		u = face->edge->next->vertex->p - face->edge->vertex->p;
		Vector3		v = face->edge->next->next->vertex->p - face->edge->vertex->p;
		Vector3		n = u.crossProduct(v).norm();
		normals.insert( pair<int, Vector3>(face->ID, n) );
	}

	for( Mesh::VertexIterator vIter = this->vertexIterator() ; !vIter.end() ; vIter++ )
	{
		int		valence = 0;
		vIter.vertex()->n = Vector3(0,0,0);

		for( Vertex::EdgeAroundIterator	around = vIter.vertex()->iterator(); !around.end(); around++ )
		{
			if(around.edge_out()->face)
			{
				int		faceId = around.edge_out()->face->ID;
				vIter.vertex()->n += normals[faceId];
				valence++;
			}
		}

		vIter.vertex()->n /= (double)valence;
		vIter.vertex()->n.normalize();
	}
}


/** Called after the mesh is created to link boundary edges */
void Mesh::linkBoundary() {
   HalfEdgeIterator hIter = halfEdgeIterator();

   for (; !hIter.end(); hIter++) {
      if (!hIter.half_edge()->face) 
	 hIter.half_edge()->vertex->edge = hIter.half_edge();
   }

    VertexIterator vIter = vertexIterator();
   Edge *next, *beg;
   while (!vIter.end()) {
      if (vIter.vertex()->isBoundary()) {
	 beg = vIter.vertex()->edge;
	 next = beg;
	 while (beg->pair->face)
	    beg = beg->pair->next;
	 beg = beg->pair;
	 beg->next = next;
	 numBoundaryVertices++;
      }
      vIter++;
   }

}

/** Computes information about the mesh:

Number of boundary loops,
Number of connected components,
Genus of the mesh.

Only one connected compoment with 0 or 1 boundary loops can be parameterize.
Singularities should be assigned in such a way that Gauss-Bonet theorem is satisfied.
*/
void Mesh::computeMeshInfo() {
   cout << "Topological information about the mesh:" << endl;
   // Number of boundary loops
   Mesh::HalfEdgeIterator hIter = halfEdgeIterator();
   for (; !hIter.end(); hIter++) {
      hIter.half_edge()->check = false;
   }
   numBoundaryLoops = 0;
   for (hIter.reset(); !hIter.end(); hIter++) {
      Edge * e = hIter.half_edge();
      if (e->face)
	 e->check = true;
      else if (!e->check) {
	 Edge * beg = NULL;
	 while (e != beg) {
	    if (!beg) beg = e;
	    check_error(!e, "Building the mesh failed, probem with input format.");
	    e->check = true;
	    e = e->next;
	 }
	 numBoundaryLoops++;
      }
   }
   cout << "Mesh has " << numBoundaryLoops << " boundary loops." << endl;
   // Number of connected components
   numConnectedComponents = 0;
   Mesh::FaceIterator fIter = faceIterator();
   for (; !fIter.end(); fIter++) {
      fIter.face()->check = false;
   }
   stack<Edge *> toVisit;
   for (fIter.reset(); !fIter.end(); fIter++) {
      if (!fIter.face()->check) {
	 numConnectedComponents++;
	 toVisit.push(fIter.face()->edge);
	 while (!toVisit.empty()) {
	    Face * fIn = toVisit.top()->face; 
	    toVisit.pop();
	    fIn->check = true;     
	    Face::EdgeAroundIterator iter = fIn->iterator();
	    for (; !iter.end(); iter++) 
	       if (iter.edge()->pair->face && !iter.edge()->pair->face->check)
		  toVisit.push(iter.edge()->pair);
	 }
      }
   }
   cout << "Mesh has " << numConnectedComponents << " connected components." << endl;
   // Genus number
   check_error(numConnectedComponents == 0, "The mesh read is empty.");
   numGenus = 
      (1 - (numVertices() - numEdges() + numFaces() + numBoundaryLoops ) / 2) / numConnectedComponents;
   cout << "Mesh is genus " << numGenus << "." << endl;
}

/** Check if all the vertices in the mesh have at least on edge coming out of them */
bool Mesh::checkVertexConection() {
   Mesh::FaceIterator fIter = faceIterator();
   Mesh::VertexIterator vIter = vertexIterator();
   bool conectedVert = true;

   for (;!vIter.end(); vIter++)
      vIter.vertex()->check = false;

   for (fIter.reset(); !fIter.end(); fIter++) {
      Face::EdgeAroundIterator around = fIter.face()->iterator();
      for (;!around.end();around++)
	 around.vertex()->check = true;
   }
   for (vIter.reset(); !vIter.end(); vIter++) {
      if (!vIter.vertex()->check) {
	 cerr << "Vertex " << vIter.vertex()->ID << " is not connected." << endl;
	 conectedVert = false;
      }
   }

   return conectedVert;
}

/** Manifoldness check: only one disk should be adjusted on any vertex */
bool Mesh::checkManifold() {
   Mesh::HalfEdgeIterator eIter = halfEdgeIterator();
   Mesh::VertexIterator vIter = vertexIterator();
   bool manifold = true;

   for (;!eIter.end(); eIter++)
      eIter.half_edge()->check = false;

   for (vIter.reset(); !vIter.end(); vIter++) {
      Vertex::EdgeAroundIterator around = vIter.vertex()->iterator();
      for (;!around.end();around++)
	 around.edge_out()->check = true;
   }

   for (eIter.reset(); !eIter.end(); eIter++) {
      if (!eIter.half_edge()->check) {
	 cerr << "Mesh is not manifold - more then one disk at vertex " 
	    << eIter.half_edge()->vertex->ID << endl;
	 manifold = false;
	 break;
      }
   }

   return manifold;
}

void Mesh::CenterAndNormalize()
{
	double maxX, maxY, maxZ, minX, minY, minZ;
	maxX = maxY = maxZ = -1e38;
	minX = minY = minZ = 1e38;

	for( VertexIterator vIter = vertexIterator() ; !vIter.end() ; vIter++ )
	{
		if(vIter.vertex()->p.x() > maxX) maxX = vIter.vertex()->p.x();
		if(vIter.vertex()->p.y() > maxY) maxY = vIter.vertex()->p.y();
		if(vIter.vertex()->p.z() > maxZ) maxZ = vIter.vertex()->p.z();
		if(vIter.vertex()->p.x() < minX) minX = vIter.vertex()->p.x();
		if(vIter.vertex()->p.y() < minY) minY = vIter.vertex()->p.y();
		if(vIter.vertex()->p.z() < minZ) minZ = vIter.vertex()->p.z();
	}

	Vector3 min(minX,minY,minZ);
	Vector3 max(maxX,maxY,maxZ);

	Vector3 center = min + (max - min) * 0.5;

	//Vector3 center(0,0,0);
	//for( VertexIterator vIter = vertexIterator() ; !vIter.end() ; vIter++ )
	//{
	//	center += vIter.vertex()->p;
	//}
	//center = center / (double)this->numVertices();

	for( VertexIterator vIter = vertexIterator() ; !vIter.end() ; vIter++ )
	{
		vIter.vertex()->p -= center;
	}

	double diag = (max - min).abs() / 2.0;
	double scale = 1.0 / diag;
	for( VertexIterator vIter = vertexIterator() ; !vIter.end() ; vIter++ )
	{
		vIter.vertex()->p *= scale;
	}
}

/** Compute initial alpha angles in the mesh. 

This function assumes that the mesh is triangular. 
Alpha angels are stored at each edge and correspond
to an angle opposite to this edge.
*/
//void Mesh::computeInitAngles() {
//   Mesh::HalfEdgeIterator eIter = halfEdgeIterator();
//   for (; !eIter.end(); eIter++) {
//      Edge * e = eIter.half_edge();
//      if (e->face) {
//	 double l1 = e->length;
//	 double l2 = e->next->length;
//	 double l3 = e->next->next->length;
//	 e->alphaOposite = acos((l2*l2 + l3*l3 - l1*l1)/(2.0*l2*l3));
//      }
//      else {
//	 e->alphaOposite = 0;
//      }
//   }
//   for (eIter.reset(); !eIter.end(); eIter++) {
//      Edge * e = eIter.half_edge();
//      e->setTheta(M_PI - e->alphaOposite - e->pair->alphaOposite);
//   }
//}

/** Loads mesh from obj file

Standard format for OBJ file is

v double double double

v double double double

f int int int

f int int int

Files with vt tags also can be parsed
*/
void Mesh::readOBJ(const char * obj_file) {
   string front;
   string v = "v", vt = "vt", f = "f";
   Vector3 vert;
   vector<int> verts;
   vector<Vector3> uvVec;
   vector<int> uvs;
   char etc;
   int id;

   ifstream in(obj_file);

   check_error(!in, "Cannot find input obj file.");

   bool hasUV = false;

   while(!in.eof() || in.peek() != EOF) {
      in >> ws; //extract whitespaces
      if (in.eof() || in.peek() == EOF)
	 break;
      if (in.peek() == '#') {
	 in.ignore(300, '\n');
      }
      else {
	 in >> front;
	 if (front == v) {
	    in >> vert.x() >> vert.y() >> vert.z();
	    addVertex(vert);
	 }
	 else if (front == vt){
	    in >> vert.x() >> vert.y();
	    vert.z() = 0;
	    uvVec.push_back(vert);
	    hasUV = true;
	 }
	 else if (front == f) {
	    verts.clear();
	    uvs.clear();
	    while (in >> id) {

	       check_error(id > numVertices(), "Problem with input OBJ file.");

	       verts.push_back(id-1);
	       bool vtNow = true;
	       if (in.peek() == '/'){
    		  in >> etc;
		  if (in.peek() != '/') {
      		     in >> id;
		     check_warn(id > numVertices(), "Texture coordinate index is greater then number of vertices.");
		     if (id < numVertices() && hasUV) {
			uvs.push_back(id-1);
			vtNow = false;
		     }
		  }
	       }
	       if (in.peek() == '/'){
		  int tmp;
   		  in >> etc;
		  in >> tmp;
	       }
	       if (hasUV && vtNow) {
		  uvs.push_back(id-1);
	       }
	    }
	    in.clear(in.rdstate() & ~ios::failbit);
	    Face * f = addFace(verts);

	    if (hasUV && uvs.size() != 0){
	       int k = 0;
	       for (Face::EdgeAroundIterator e = f->iterator(); !e.end(); e++, k++)
		  e.vertex()->uv = uvVec[uvs[k]];
	    }
	 }
	 else {
	    string line;
	    getline(in, line);
	    cout << "Warning, line: " << line << " ignored." << endl;	
	 }
      }
   }

   in.close();

   // Finnish building the mesh, should be called after each parse.
   finishMesh();
}

/* Write mesh in OBJ format to obj_file parameter */
void Mesh::writeOBJ(const char * obj_file) {
   ofstream out(obj_file);
   check_error (!out, "Cannot find output file.");

   Mesh::VertexIterator vItr = vertexIterator();
   for (vItr.reset(); !vItr.end(); vItr++)
      out << "v " << vItr.vertex()->p << endl;

   for (vItr.reset(); !vItr.end(); vItr++)
      out << "vt " << vItr.vertex()->uv.x() << " " << vItr.vertex()->uv.y() << endl;

   Mesh::FaceIterator fItr = faceIterator();
   for (fItr.reset(); !fItr.end(); fItr++) {
      Face::EdgeAroundIterator around = fItr.face()->iterator();
      out << "f";
      for ( ; !around.end(); around++)
	 out << " " << (around.vertex()->ID + 1) << "/" << (around.vertex()->ID + 1);
      out << endl;
   }
}

/* Write mesh in VT format (only text coodrinates) to vt_file parameter */
void Mesh::writeVT(const char * vt_file) {
   ofstream out(vt_file);
   check_error (!out, "Cannot find output file.");

   Mesh::VertexIterator vItr = vertexIterator();
   for (vItr.reset(); !vItr.end(); vItr++)
      out << "vt " << vItr.vertex()->uv.x() << " " << vItr.vertex()->uv.y() << endl;
}

/** Loads mesh from con (connectivity and edge length) file

Format for CON file is

int (number of vertices)

int (number of faces)

f int int int 

double double double

f int int int 

double double double

Where first double in each line if edge length from first vertex to second one, second double is edge length from 
second vertex to the third one, and the third double is an edge length from third vertex to the first one.
*/
//void Mesh::readCON(const char * conn_file) {
//   string front;
//   string v = "v", f = "f";
//   Vector3 vert;
//   vector<int> verts;
//   int id;
//
//   ifstream in(conn_file);
//
//   check_error(!in, "Cannot find input obj file.");
//
//   int NumFaces, NumVertices;
//
//   in >> NumVertices;
//   in >> NumFaces;
//
//   for (int i = 0; i < NumVertices; i++) {
//      addVertex(vert);
//   }
//   while(!in.eof() || in.peek() != EOF) {
//      in >> ws;
//      if (in.eof() || in.peek() == EOF)
//	 break;
//      if (in.peek() == '#') {
//	 in.ignore(300, '\n');
//      }
//      else {
//	 in >> front;
//	 if (front == f) {
//	    verts.clear();
//	    
//	    for (int i = 0; i < 3; i++) {
//	       in >> id;
//	       check_error(id > numVertices(), "Problem with input CON file.");
//	       verts.push_back(id-1);
//	    }
//	    in.clear(in.rdstate() & ~ios::failbit);
//	    Face * f = addFace(verts);
//
//	    for (Face::EdgeAroundIterator e = f->iterator(); !e.end(); e++) {
//	       in >> e.edge()->length;
//	    }
//	 }
//	 else {
//	    string line;
//	    getline(in, line);
//	    cout << "Warning, line: " << line << " ignored." << endl;	
//	 }
//      }
//   }
//
//   in.close();
//   // Finnish building the mesh, should be called after each parse.
//   finishMesh();
//}

/** Write out the result parameterization in 
"connectivity" format. The edge lengths in the result
file correspond to edge lengths of the final parameterization.
*/
//void Mesh::writeCON(const char * conn_file) {
//   ofstream out(conn_file);
//   check_error (!out, "Cannot find output file.");
//
//   out << numVertices() << endl;
//   out << numFaces() << endl;
//
//   Mesh::FaceIterator fItr = faceIterator();
//   for (fItr.reset(); !fItr.end(); fItr++) {
//      Face::EdgeAroundIterator around = fItr.face()->iterator();
//      out << "f";
//      for ( ; !around.end(); around++)
//	 out << " " << (around.vertex()->ID + 1);
//      out << endl;
//
//      for (around.reset() ; !around.end(); around++) 
//	 out << " " << setprecision(8) << around.edge()->length;
//      out << endl;
//   }
//
//   out.close();
//}

/** Reads edge cuts from the input file and cuts the mesh over those edges.

Format of the file is list of the follwing lines:

vert_from_id vert__to_id face

The information in the file is somewhat redundant, but as we allow for non-regular 
meshes, there are might be two half edges with same end points).

After the mesh is cut it should be a set of developable surfaces, 
so that UV coordinates can be assigned properly.

*/
void Mesh::readCUT_EDGES(const char * edge_file) {
   ifstream in(edge_file);
   check_error(!in, "Cannot find input file for edge cuts.");

   vector<Edge *> edges_to_cut;

   int i, j, f;
   while (in >> i) {
      in >> j >> f;

      if ((i < 0) || (i >= numVertices()) || (j < 0) || (j >= numVertices()))
	 cout << "ID of a vertex in the edges to cut file is out of bounds." << endl;
      else {

	 i--; j--; f--;
   
	 Edge eTmp;
         eTmp.vertex = vertices[i];		
	 Edge eTmpPair;
   	 eTmpPair.vertex = vertices[j];
         eTmp.pair = &eTmpPair;
         eTmpPair.pair = &eTmp;

         Edge * e = NULL;
	 Mesh::EdgeSet::iterator eIn = edges.find(&eTmp);
	 if (eIn != edges.end())
	    e = *eIn;
	 while (eIn != edges.end()) {
	    if ((*eIn)->face && (*eIn)->face->ID == f) {
	       e = *eIn;
	       break;
	    }
	    eIn++;
	 }
	 if (e)
	    edges_to_cut.push_back(e);
      }
   }

   while (!edges_to_cut.empty()) {
      Edge * e = edges_to_cut.back();
      if (e) {
	    cutAlongEdge(e);
      }
      edges_to_cut.pop_back();
   }

   EdgeSet newedges(edges.begin(), edges.end());
   edges = newedges;
}

/** Reads cone singularities allowed by the user.

Format of the file is list of follwing lines:

vertex_id min_val_for_singularity max_val_for_singularity

Values for singularities are expected as multiples of Pi.
Minimum and maximum values for singularities can be equal.
*/
void Mesh::readCONE_VERT(const char * vert_file) {
   ifstream in(vert_file);
   check_error(!in, "Cannot find input file for cone singularities.");

   int id = 0;
   double minC = 0, maxC = 0;

   while (in >> id) {
      in >> minC >> maxC;
      if ((id < 0) || (id >= numVertices()))
	 cout << "ID of a vertex in the cone singularities file is out of bounds." << endl;
      else {
	 id--;
         double valenceAlow = (double)vertices[id]->getValence();
	 check_warn ( maxC <=0, "Invalid cone angle. Cone angle read is less or equal to zero.");
	 check_warn ( minC >= valenceAlow, 
	    "Invalid cone angle. Cone angle read is greater or equal to what valence of vertex allows.");
	 vertices[id]->min_cone_angle = minC;
	 vertices[id]->max_cone_angle = maxC;
	 vertices[id]->constrained = true;
      }
   }
}


/** Compute total Gaussian curvature that cone singularities allow
and check if it satisfy Gauss-Bonet theorem */
void Mesh::checkGaussianCurvature() {
   double minConeGC = 0; 
   double maxConeGC = 0;
  
   for (int i = 0; i < numVertices(); i++)   {
      if (vertices[i]->constrained) {
	 if (vertices[i]->isBoundary()) {
      	    maxConeGC += (M_PI - vertices[i]->min_cone_angle*M_PI);
	    minConeGC += (M_PI - vertices[i]->max_cone_angle*M_PI);
	 }
	 else {
      	    maxConeGC += (2*M_PI - vertices[i]->min_cone_angle*M_PI);
	    minConeGC += (2*M_PI - vertices[i]->max_cone_angle*M_PI);
	 }
      }
      else if (vertices[i]->isBoundary()) {
	 maxConeGC +=  M_PI;
	 minConeGC += -M_PI;;
      }
   }

   double eNum1 = minConeGC / (2*M_PI);
   double eNum2 = maxConeGC / (2*M_PI);
   double eNumCurr = (double)(2 - 2*numGenus - numBoundaryLoops)/(double)numConnectedComponents;
   cout << "Checking Gauss-Bonet Theorem:" << endl;
   cout << "Minimum euler number for this set up is " << eNum1 << endl;
   cout << "Maximum euler number for this set up is " << eNum2 << endl;
   cout << "Euler number of the mesh is " << eNumCurr << "\n" << endl;

   check_warn((eNumCurr < eNum1) || (eNumCurr > eNum2), 
      "The mesh cannot be developed with current set of cone singularities.\n");
}

