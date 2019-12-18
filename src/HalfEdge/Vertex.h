/** \file Vertex.h Declarations and definitions for the Vertex class.
    \author Liliya Kharevych
    \date March 2006
*/
#ifndef VERTEX_CM
#define VERTEX_CM
#include "Edge.h"

/** Vertex for the mesh datastructure.
*/
class Vertex
{
private:
   const Vertex& operator=(const Vertex& rhs){ return *this;}
public:

   // pointers for mesh structure
   Edge * edge;

   // geometry data
   Vector3 p;
   Vector3 n;
   Vector3 uv;

   // Cone angle data (multiple of pi is stored)
   double min_cone_angle;
   double max_cone_angle;
   bool constrained;
   double cone_angle;

   double radius;

   double weight;

   // distinct id
   int ID;
   int patchID;

   // to check various iterations
   bool check;

   Vertex (const Vector3 & _pos = Vector3(0,0,0)): 
   edge(NULL), p(_pos), uv(0,0,0)
   {
      cone_angle = min_cone_angle = max_cone_angle = 2;
      constrained = false;
	  radius = 0.0;
	  weight = 0.0;
      patchID = 0;
   }

   // Assignment (copy only geometrical information)
   void assignData(const Vertex& rhs){
      if (this != &rhs) {
	 p = rhs.p;
	 uv = rhs.uv;
	 cone_angle = rhs.cone_angle;
	 min_cone_angle = rhs.min_cone_angle;
	 max_cone_angle = rhs.max_cone_angle;
      }
   }


   bool isBoundary() {
      return (edge && !edge->face);
   }

  /** The iterator that visits edges or vertices in one ring of the current face in order. */
  class EdgeAroundIterator {
   private:
      Edge * endI;
      Edge * run;

   public:
      EdgeAroundIterator(Edge * e) {
	 endI = NULL;
	 run = e;
      }
      EdgeAroundIterator& operator++( void ){
	 if (!endI) endI = run;
	 run = run->pair->next;
	 return *this;
      }
      EdgeAroundIterator operator++( int ){
	 EdgeAroundIterator r = *this; ++*this; return r;
      }

      Edge * edge_out( void ) const { return run; }
      Edge * & edge_out( void )     { return run; }

      Vertex * vertex( void ) const { return run->vertex; }
      Vertex * & vertex( void )     { return run->vertex; }

      bool end(void) { return endI == run;}
   };

   EdgeAroundIterator iterator() {return EdgeAroundIterator(edge);}
   EdgeAroundIterator iterator(Edge * eFrom) {return EdgeAroundIterator(eFrom);}

   int getValence() {
      Vertex::EdgeAroundIterator iter = iterator();
      int res = 0;
      for (; !iter.end(); iter++)
	 res++;
      return res;
   }
};
#endif

