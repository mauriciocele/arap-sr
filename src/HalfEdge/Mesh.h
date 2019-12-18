/** \file Mesh.h Declarations for the Mesh class.
    \author Liliya Kharevych
    \date March 2006
*/

#ifndef MESH_CM
#define MESH_CM
#include <vector>
#include <set>
#include <iostream>
#include <algorithm>    // std::min

#include "Edge.h"
#include "Face.h"
#include "Vertex.h"

/** 
Implementation of the half edge datastructure.
*/

class Mesh
{
private:
   /** Comparison functions for edge sets

   Lexicographical order by indices of start end end point of the edge
   */
   struct edge_comp {
      edge_comp() {}
      bool operator()( const Edge *e1, const Edge *e2 ) const {
	 int b1 = e1->vertex->ID, t1 = e1->pair->vertex->ID;
	 int b2 = e2->vertex->ID, t2 = e2->pair->vertex->ID; 

	 int min1, min2, max1, max2;
	 min1 = std::min(b1, t1); min2 = std::min(b2, t2);
	 max1 = std::max(b1, t1); max2 = std::max(b2, t2);
	 if (min1 == min2 && max1 == max2){
	    if (b1 < b2) return true;
	    return (b1 == b2 && t1 < t2);
	 }
	 if (min1 < min2)     return true;
	 return (min1 == min2 && max1 < max2);	
      }
   };


   std::vector<Face *> faces;
   std::multiset<Edge*, Mesh::edge_comp> edges;

public:
	std::vector<Vertex *> vertices;

   /** Number of vertices (or also edges) at the boundary of the mesh   */
   int numBoundaryVertices;

   /** Number of boundary loops in the mesh */
   int numBoundaryLoops;
   /** Number of connected components */
   int numConnectedComponents;

   /** Genus of the mesh */
   int numGenus;

   typedef std::multiset<Edge*, Mesh::edge_comp>  EdgeSet;
   typedef std::vector<Face *> FaceSet;
   typedef std::vector<Vertex *> VertexSet;

   Mesh(void);
   ~Mesh(void);
   void clear();

   Vertex * addVertex(const Vector3 & p);
   Face * addFace(std::vector<int> faceVerts);
   Edge * addEdge(int i, int j);
   Edge * addEdge(Edge * e) {
      edges.insert(e);
      return e;
   }

   bool cutAlongEdge(Edge * forward);

   inline Vertex * vertexAt (int i) {return vertices[i];}
   inline Face * faceAt (int i) {return faces[i];}

   inline int numFaces()     { return (int)faces.size();    }
   inline int numVertices()  { return (int)vertices.size(); }
   inline int numEdges()     { return (int)edges.size()/2;  }
   inline int numBoundary()  { return numBoundaryVertices;  }

   void linkBoundary();
   bool checkManifold();
   bool checkVertexConection();
   void checkGaussianCurvature();

   /**
   Computes number of connected components, number of boundary loops and genus of the mesh.

   Uses variation of Euler Formula: V - E + F = 2 (c - g)
   Where V - number of vertices.
   E - number of edges.
   F - number of faces.
   c - number of connected components.
   */
   void computeMeshInfo();
   void CenterAndNormalize();

   void readOBJ(const char * obj_file);
   void readCONE_VERT(const char * vert_file);
   void readCUT_EDGES(const char * edge_file);
   void writeOBJ(const char * obj_file);
   void writeVT(const char * vt_file);

   void computeLengths();
   void computeCirclePackingMetric();
   void computeNormals();

   void finishMesh() {
      linkBoundary();
      checkManifold();
      checkVertexConection();
      std::cout << "*--------------------*" << std::endl;		
      std::cout << "* Faces:    " << numFaces() << std::endl;		
      std::cout << "* Edges:    " << numEdges() << std::endl;		
      std::cout << "* Vertices: " << numVertices() << std::endl;		
      std::cout << "*--------------------*\n" << std::endl;		
   }

   class FaceIterator {
   private:
      FaceSet::iterator fIter;
      FaceSet * facesPtr;
   public:
      FaceIterator() {
	 facesPtr = NULL;
      }

      FaceIterator(FaceSet * _faces) {
	 facesPtr = _faces;
	 fIter = _faces->begin();
      }
      FaceIterator& operator++( void ){
	 fIter++;
	 return *this;
      }
      FaceIterator operator++( int ){
	 FaceIterator r = *this; ++*this; return r;
      }
      FaceIterator& operator--( void ){
	 fIter--;
	 return *this;
      }
      FaceIterator operator--( int ){
	 FaceIterator r = *this; --*this; return r;
      }

      Face * face( void ) const { return *fIter; }
      //Face * & face( void )     { return *fIter; }

      void reset() {fIter = facesPtr->begin(); }
      bool end(void) { return fIter == facesPtr->end();
      ;}
   };

   class VertexIterator {
   private:
      VertexSet::iterator vIter;
      VertexSet * verticesPtr;
   public:
      VertexIterator() {
	 verticesPtr = NULL;
      }

      VertexIterator(VertexSet * _vertices) {
	 vIter = _vertices->begin();
	 verticesPtr = _vertices;
      }
      VertexIterator& operator++( void ){
	 vIter++;
	 return *this;
      }
      VertexIterator operator++( int ){
	 VertexIterator r = *this; ++*this; return r;
      }
      VertexIterator& operator--( void ){
	 vIter--;
	 return *this;
      }
      VertexIterator operator--( int ){
	 VertexIterator r = *this; --*this; return r;
      }

      Vertex * vertex( void ) const { return *vIter; }
      //Vertex * & vertex( void )     { return *vIter; }

      void reset() {vIter = verticesPtr->begin();}
      bool end(void) { return vIter == verticesPtr->end();}
   };

   class HalfEdgeIterator {
   private:
      EdgeSet::iterator eIter;
      EdgeSet * edgesPtr;
   public:
      HalfEdgeIterator() {
	 edgesPtr = NULL;
      }
      HalfEdgeIterator(EdgeSet * _edges) {
	 eIter = _edges->begin();
	 edgesPtr = _edges;
      }
      HalfEdgeIterator& operator++( void ){
	 eIter++;
	 return *this;
      }
      HalfEdgeIterator operator++( int ){
	 HalfEdgeIterator r = *this; ++*this; return r;
      }
      HalfEdgeIterator& operator--( void ){
	 eIter--;
	 return *this;
      }
      HalfEdgeIterator operator--( int ){
	 HalfEdgeIterator r = *this; --*this; return r;
      }

      Edge * half_edge( void ) const { return *eIter; }
      //		Edge * & operator*( void )     { return *eIter; }

      void find (Edge * eTmp) {
	 eIter = edgesPtr->find(eTmp);
      }

      void reset() {eIter = edgesPtr->begin(); }
      bool end (void) { return eIter == edgesPtr->end();}
   };

   class EdgeIterator {
   private:
      EdgeSet::iterator eIter;
      EdgeSet * edgesPtr;

   public:
      EdgeIterator() {
	 edgesPtr = NULL;
      }

      EdgeIterator(EdgeSet * _edges) {
	 eIter = _edges->begin();
	 edgesPtr = _edges;
      }
      EdgeIterator& operator++( void ){
	 eIter++; eIter++;
	 return *this;
      }
      EdgeIterator operator++( int ){
	 EdgeIterator r = *this; ++*this; return r;
      }
      EdgeIterator& operator--( void ){
	 eIter--; eIter--;
	 return *this;
      }
      EdgeIterator operator--( int ){
	 EdgeIterator r = *this; --*this; return r;
      }

      Edge * edge( void ) const { return *eIter; }
      //		Edge * & operator*( void )     { return *eIter; }

      void reset() {eIter = edgesPtr->begin(); }
      bool end(void) { return eIter == edgesPtr->end();}
   };



   FaceIterator faceIterator() {return FaceIterator(&faces);}
   VertexIterator vertexIterator() {return VertexIterator(&vertices);}
   HalfEdgeIterator halfEdgeIterator() {return HalfEdgeIterator(&edges);}
   EdgeIterator edgeIterator() {return EdgeIterator(&edges);}
};
#endif

