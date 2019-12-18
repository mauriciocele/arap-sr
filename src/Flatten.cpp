// Flatten.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#ifdef WIN32
#define NOMINMAX 
#include <windows.h>
#endif

#include <GL/gl.h>

#if _WIN64
#include <GL/x64/glut.h>
#else
#include <GL/glut.h>
#endif

#include "GA/c3ga.h"
#include "GA/c3ga_util.h"
#include "GA/gl_util.h"

#include "primitivedraw.h"
#include "gahelper.h"
#include "laplacian.h"
#include "LengthCorrector.h"

#include <memory>

#include <vector>
#include <queue>
#include <map>
#include <numerics.h>
#include "HalfEdge\Mesh.h"

#include "PriorityQueue.h"

#include <geometry.H>
#include "MatrixNxM.h"
#include "Matrix3x3.h"
#include "TMatrix.h"
#include "ICP.h"
#include <Taucs\taucs_interface.h>
#include <Eigen\core>
#include <Eigen\LU>

#include <ppl.h>

#include "Benchs.h"

const char *WINDOW_TITLE = "Interactive 3D Shape Deformation using Conformal Geometric Algebra";

// GLUT state information
int g_viewportWidth = 800;
int g_viewportHeight = 600;

void display();
void reshape(GLint width, GLint height);
void MouseButton(int button, int state, int x, int y);
void MouseMotion(int x, int y);
void KeyboardUpFunc(unsigned char key, int x, int y);
void SpecialFunc(int key, int x, int y);
void SpecialUpFunc(int key, int x, int y);
void Idle();
void DestroyWindow();

//using namespace boost;
using namespace c3ga;
using namespace std;
using namespace numerics;
using namespace Eigen;

class Camera
{
public:
	float		pos[3];
	float		fw[3];
	float		up[3];
	float		translateVel;
	float		rotateVel;

	Camera()
	{
		float		_pos[] = { 0, 0, -2};
		float		_fw[] = { 0, 0, 1 };
		float		_up[] = { 0, 1, 0 };

		translateVel = 0.005;
		rotateVel = 0.005;
		memcpy(pos, _pos, sizeof(float)*3);
		memcpy(fw, _fw, sizeof(float)*3);
		memcpy(up, _up, sizeof(float)*3);
	}

	void glLookAt()
	{
		gluLookAt( pos[0], pos[1], pos[2], fw[0],  fw[1],  fw[2], up[0],  up[1],  up[2] );
	}
};

class VertexDescriptor
{
public:
	vectorE3GA deformedPosition; //deformed mesh position
	normalizedPoint position; //primary backup
	normalizedPoint positionConstrained; //primary backup
	vectorE3GA laplacianCoordinate; //laplacian Coordinate
	vectorE3GA normal; //for rendering (lighting)
	vectorE3GA normalOrig; //primary backup
	rotor M;
	//rotor R;
	Matrix4d JtJ;
	Vector4d b;
	bool isBoundary;

	VertexDescriptor()
	{
	}
};

class Handle
{
public:
	TRversor R;
	TRversor T;
	dualSphere dS;
	bool fixed;
	dualPlane Plane;
	//int extrema;

	Handle(dualSphere dS, bool fixed, dualPlane dP)
	{
		normalizedPoint x = DualSphereCenter(dS);
		T = _TRversor( exp(-0.5 * _vectorE3GA(x)^ni) );
		R = _TRversor(1.0);
		Plane = dP;
		this->dS = inverse(T) * dS * T;
		this->fixed = fixed;
	}

	TRversor GetTRVersor()
	{
		return T * R;
	}
};

//class Extrema
//{
//public:
//	TRversor R;
//	TRversor T;
//	dualSphere dS;
//	int handle;
//
//	Extrema(dualSphere dS)
//	{
//		normalizedPoint x = DualSphereCenter(dS);
//		T = _TRversor( exp(-0.5 * _vectorE3GA(x)^ni) );
//		R = _TRversor(1.0);
//		this->dS = inverse(T) * dS * T;
//	}
//
//	TRversor GetTRVersor()
//	{
//		return T * R;
//	}
//};

Camera g_camera;
Mesh mesh, meshLow;
map<int, int> mapping, rotorClusters;
vectorE3GA g_prevMousePos;
bool g_rotateModel = false;
bool g_rotateModelOutOfPlane = false;
rotor g_modelRotor = _rotor(1.0);
bool g_rotateKeyRotors = false;
bool g_translateKeyRotors = false;
bool g_computeBasis = false;
float g_dragDistance = -1.0f;
int g_dragObject;
std::shared_ptr<SparseMatrix> A;
std::shared_ptr<SparseMatrix> AHi;
std::set<int> constraints;
std::set<int> fconstraints;
int Lc, LcHi;
int systemType = LaplaceBeltrami; //MeanValue; //LaplaceBeltrami
bool g_showSpheres = true;
bool g_showWires = false;
bool g_iterateManyTimes = false;
double *b3 = nullptr;
double *xyz = nullptr;
double *b3Hi = nullptr;
double *xyzHi = nullptr;
double g_meshArea = 0.0;

bool g_automaticAnimation = false;
bool g_convergence = false;
Benchmarks benchs;

vector<std::shared_ptr<VertexDescriptor>> vertexDescriptors;
vector<std::shared_ptr<VertexDescriptor>> vertexDescriptorsLow;
std::vector<std::shared_ptr<Handle>> handles;
//std::vector<boost::shared_ptr<Extrema>> extremas;


void ComputeLaplacianCoordinates(std::shared_ptr<SparseMatrix> A, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptors)
{
	for( int i = 0 ; i < vertexDescriptors.size() ; ++i )
		vertexDescriptors[i]->laplacianCoordinate = _vectorE3GA(0.0, 0.0, 0.0);
	
	auto numRows = A->numRows();

	for( int i = 0; i < numRows ; ++i)
	{
		SparseMatrix::RowIterator aIter = A->iterator(i);
		for( ; !aIter.end() ; aIter++ )
		{
			auto j = aIter.columnIndex();
			vertexDescriptors[i]->laplacianCoordinate += _vectorE3GA(vertexDescriptors[j]->position) * aIter.value();
		}
	}
}

int PreFactor(std::shared_ptr<SparseMatrix> A, std::set<int>& constraints)
{
	int Lc = CreateMatrix(A->numRows(), A->numColumns());

	auto numRows = A->numRows();

	for( int i = 0; i < numRows ; ++i)
	{
		if(constraints.find(i) == constraints.end())
		{
			SparseMatrix::RowIterator aIter = A->iterator(i);
			for (; !aIter.end(); aIter++)
			{
				auto j = aIter.columnIndex();
				SetMatrixEntry(Lc, i, j, (*A)(i, j));
			}
		}
		else
		{
			SetMatrixEntry(Lc, i, i, 1.0);
		}
	}
	FactorATA(Lc);

	return Lc;
}

TRversor RigidBodyMotion(const vector<normalizedPoint>& p, const vector<normalizedPoint>& q)
{
	if(p.size() <= 3)
	{
		//TODO: If 3 points then FourPlanes method here, 
		//if 2 points then Quaternion method here, 
		//if 1 point then just regular GA here
		return _TRversor(1.0);
	}

	vector<normalizedPoint>	x(p.begin(), p.end());
	set<int> A;
	for(int i = 0; i < p.size() ; ++i)
		A.insert(i);

	PriorityQueue<double> priorities(x.size(), false);
	for(auto iter = A.begin(); iter != A.end() ; iter++ )
	{
		int i = *iter;
		priorities.insert(i, abs(_double(_dualPlane(x[i] - q[i])<<_dualPlane(x[i] - q[i]))) );
	}

	mv V = mv(1.0);
	while(A.size() > 0)
	{
		int a = priorities.getTopIndex();
		if(priorities.getTopPriority() < 1e-4)
		{
			priorities.pop();
			A.erase(a);
			continue;
		}
		priorities.pop();
		dualPlane R = _dualPlane(x[a] - q[a]);
		V = R * V;
		A.erase(a);
		for(auto iter = A.begin(); iter != A.end() ; iter++ )
		{
			int i = *iter;
			x[i] = - R * x[i] * inverse(R);
			priorities.changePriority(i, abs(_double(_dualPlane(x[i] - q[i])<<_dualPlane(x[i] - q[i]))) );
		}
	}
	return _TRversor(V);
}

double meshArea(Mesh *mesh)
{
	double area = 0.0;
	for (Mesh::FaceIterator fIter = mesh->faceIterator(); !fIter.end(); fIter++)
	{
		Face		*face = fIter.face();
		Vector3		u = face->edge->next->vertex->p - face->edge->vertex->p;
		Vector3		v = face->edge->next->next->vertex->p - face->edge->vertex->p;
		area += 0.5 * u.crossProduct(v).abs();
	}
	return area;
}

void meshMapping(Mesh *meshLow, Mesh *mesh, map<int, int>& mapping)
{
	for (Mesh::VertexIterator vIterLow = meshLow->vertexIterator(); !vIterLow.end(); vIterLow++)
	{
		double minDist = 1e38;
		int mappedPoint = -1;
		Vector3& pLow = vIterLow.vertex()->p;
		for (Mesh::VertexIterator vIterHi = mesh->vertexIterator(); !vIterHi.end(); vIterHi++)
		{
			Vector3& p = vIterHi.vertex()->p;
			double dist = (pLow - p).abs();
			if (dist < minDist) {
				minDist = dist;
				mappedPoint = vIterHi.vertex()->ID;
			}
		}
		mapping[vIterLow.vertex()->ID] = mappedPoint;
	}

	int n = mesh->numVertices();
	bool *visited = new bool[n];
	queue<int> assigned;
	for (int i = 0; i < n; ++i)
	{
		visited[i] = false;
	}

	for (map<int, int>::iterator iter = mapping.begin(); iter != mapping.end(); iter++){
		visited[iter->second] = true;
		rotorClusters[iter->second] = iter->first;
		assigned.push(iter->second);
	}

	while (!assigned.empty())
	{
		int i = assigned.front();
		assigned.pop();
		visited[i] = true;

		for (Vertex::EdgeAroundIterator edgeAroundIter = mesh->vertexAt(i)->iterator(); !edgeAroundIter.end(); edgeAroundIter++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;
			if (visited[j] == false)
			{
				rotorClusters[j] = rotorClusters[i];
				assigned.push(j);
			}
		}
	}

	delete[] visited;
}

int main(int argc, char* argv[])
{
	InitTaucsInterface();

	mesh.readOBJ("cactus2.obj"); //armadillo-5k-smooth.obj female.obj david.obj rabbit.obj tyra.obj horse.obj cylinder.obj bar.obj planewithpeaks.obj dragon.obj catHead.obj  cactus.obj  bunny.obj  764_hand-olivier-10kf.obj armadillo.obj
	mesh.CenterAndNormalize();
	mesh.computeNormals();

	meshLow.readOBJ("cactus1.obj"); //armadillo-5k-smooth.obj female.obj david.obj rabbit.obj tyra.obj horse.obj cylinder.obj bar.obj planewithpeaks.obj dragon.obj catHead.obj  cactus.obj  bunny.obj  764_hand-olivier-10kf.obj armadillo.obj
	meshLow.CenterAndNormalize();
	meshLow.computeNormals();

	meshMapping(&meshLow, &mesh, mapping);

	// GLUT Window Initialization:
	glutInit (&argc, argv);
	glutInitWindowSize(g_viewportWidth, g_viewportHeight);
	glutInitDisplayMode( GLUT_RGB | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow(WINDOW_TITLE);

	// Register callbacks:
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutMouseFunc(MouseButton);
	glutMotionFunc(MouseMotion);
	glutKeyboardUpFunc(KeyboardUpFunc);
	glutSpecialFunc(SpecialFunc);
	glutSpecialUpFunc(SpecialUpFunc);
	glutIdleFunc(Idle);
	atexit(DestroyWindow);

	InitializeDrawing();

	dualPlane P1 = _dualPlane(c3gaPoint(.0, 0.8,.0) << (_vectorE3GA(0,1,0)*ni));
	dualPlane P2 = _dualPlane(c3gaPoint(.0,-0.8,.0) << (_vectorE3GA(0,-1,0)*ni));

	//cactus.obj
	handles.push_back(std::shared_ptr<Handle>(new Handle(_dualSphere(c3gaPoint(.0, .04, 0.7) - 0.5*SQR(0.15)*ni), false, P1)));
	handles.push_back(std::shared_ptr<Handle>(new Handle(_dualSphere(c3gaPoint(.0, .04, -0.8) - 0.5*SQR(0.15)*ni), true, P2)));

	////cylinder.obj
	//handles.push_back(std::shared_ptr<Handle>(new Handle(_dualSphere(c3gaPoint(.0, 0.9, .0) - 0.5*SQR(0.25)*ni), false, P1)));
	//handles.push_back(std::shared_ptr<Handle>(new Handle(_dualSphere(c3gaPoint(.0, -0.9, .0) - 0.5*SQR(0.25)*ni), true, P2)));

	//Armadillo Pie y Mano
	//handles.push_back(std::shared_ptr<Handle>(new Handle(_dualSphere(c3gaPoint(-.5, 0.45, -.3) - 0.5*SQR(0.15)*ni), false, P1)));
	//handles.push_back(std::shared_ptr<Handle>(new Handle(_dualSphere(c3gaPoint(.2, -0.6, .1) - 0.5*SQR(0.15)*ni), true, P2)));
	//Armadillo Pubis y Cabeza
	//handles.push_back(std::shared_ptr<Handle>(new Handle(_dualSphere(c3gaPoint(.0, 0.4, -.2) - 0.5*SQR(0.15)*ni), false, P1)));
	//handles.push_back(std::shared_ptr<Handle>(new Handle(_dualSphere(c3gaPoint(.0, -0.05, .1) - 0.5*SQR(0.15)*ni), true, P2)));

	//handles[0]->extrema = 0;
	//handles[1]->extrema = 1;
	//extremas[0]->handle = 0;
	//extremas[1]->handle = 1;

	for( Mesh::VertexIterator vIter = meshLow.vertexIterator() ; !vIter.end() ; vIter++ )
	{
		std::shared_ptr<VertexDescriptor> vd(new VertexDescriptor);

		vd->position = c3gaPoint( vIter.vertex()->p.x(), vIter.vertex()->p.y(), vIter.vertex()->p.z() );
		vd->normalOrig = _vectorE3GA( vIter.vertex()->n.x(), vIter.vertex()->n.y(), vIter.vertex()->n.z() );
		vd->isBoundary = vIter.vertex()->isBoundary();
		//vd->R.Identity();
		vd->M = _rotor(1.0);

		vertexDescriptorsLow.push_back(vd);

		for( int i = 0 ; i < handles.size() ;  ++i )
		{
			TRversor TR = handles[i]->GetTRVersor();

			if( _double(vd->position << (TR * handles[i]->dS * inverse(TR))) > 0 ) //inside the sphere
			//if( _double(vd->position << handles[i]->Plane) > 0 ) //positive side of the plane
			{
				if( !handles[i]->fixed )
				{
					constraints.insert(vIter.vertex()->ID);
				}
				else
				{
					fconstraints.insert(vIter.vertex()->ID);
				}
			}
		}
		if(vd->isBoundary)
		{
			fconstraints.insert(vIter.vertex()->ID);
		}
	}

	A = CreateLaplacianMatrix( &meshLow, systemType );
	
	ComputeLaplacianCoordinates(A, vertexDescriptorsLow);

	TRversor R1 = handles[0]->GetTRVersor();
	for(std::set<int>::iterator citer = constraints.begin(); citer != constraints.end() ; citer++)
		vertexDescriptorsLow[*citer]->positionConstrained = normalize(_point(inverse(R1) * vertexDescriptorsLow[*citer]->position * R1));

	TRversor R2 = handles[1]->GetTRVersor();
	for(std::set<int>::iterator citer = fconstraints.begin(); citer != fconstraints.end() ; citer++)
		vertexDescriptorsLow[*citer]->positionConstrained = normalize(_point(inverse(R2) * vertexDescriptorsLow[*citer]->position * R2));

	b3 = new double [A->numRows() * 3];
	xyz = new double[A->numRows() * 3];

	std::set<int> allconstraints = constraints;
	allconstraints.insert(fconstraints.begin(), fconstraints.end());
	Lc = PreFactor(A, allconstraints);

	g_meshArea = meshArea(&meshLow);

	for (Mesh::VertexIterator vIter = mesh.vertexIterator(); !vIter.end(); vIter++)
	{
		std::shared_ptr<VertexDescriptor> vd(new VertexDescriptor);

		vd->position = c3gaPoint(vIter.vertex()->p.x(), vIter.vertex()->p.y(), vIter.vertex()->p.z());
		vd->normalOrig = _vectorE3GA(vIter.vertex()->n.x(), vIter.vertex()->n.y(), vIter.vertex()->n.z());
		vd->isBoundary = vIter.vertex()->isBoundary();
		//vd->R.Identity();
		vd->M = _rotor(1.0);
		vertexDescriptors.push_back(vd);
	}

	AHi = CreateLaplacianMatrix(&mesh, systemType);
	ComputeLaplacianCoordinates(AHi, vertexDescriptors);
	std::set<int> constraintsHi;
	for (map<int, int>::iterator iter = mapping.begin(); iter != mapping.end(); iter++){
		constraintsHi.insert(iter->second);
	}
	b3Hi = new double[AHi->numRows() * 3];
	xyzHi = new double[AHi->numRows() * 3];
	LcHi = PreFactor(AHi, constraintsHi);

	glutMainLoop();

	return 0;
}

void transferRotations(Mesh *mesh, std::shared_ptr<SparseMatrix> A, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptors, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptorsLow)
{
	for (map<int, int>::iterator iter = rotorClusters.begin(); iter != rotorClusters.end(); iter++) {
		vertexDescriptors[iter->first]->M = vertexDescriptorsLow[iter->second]->M;
	}

	concurrency::parallel_for_each(mesh->vertices.begin(), mesh->vertices.end(), [&](Vertex * vertex)
	//for (auto vertex : mesh->vertices)
	{
		int i = vertex->ID;
		vertexDescriptors[i]->laplacianCoordinate = _vectorE3GA(0.0, 0.0, 0.0);
		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex->iterator(); !edgeAroundIter.end(); edgeAroundIter++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;
			double wij = (*A)(i, j);
			rotor &Ri = vertexDescriptors[i]->M;
			rotor &Rj = vertexDescriptors[j]->M;
			vectorE3GA V = _vectorE3GA(vertexDescriptors[j]->position) - _vectorE3GA(vertexDescriptors[i]->position);
			vectorE3GA Vp = _vectorE3GA(0.5 * wij * (Ri * V * (~Ri) + Rj * V * (~Rj)));
			vertexDescriptors[i]->laplacianCoordinate += Vp;
		}
	//}
	});
}

void SolveLinearSystemTaucs(vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptors)
{
	int n = vertexDescriptors.size();

	for( int i = 0 ; i < n ; ++i )
	{
		auto &laplacianCoordinate = vertexDescriptors[i]->laplacianCoordinate;
		b3[i] = laplacianCoordinate.e1();
		b3[i+n] = laplacianCoordinate.e2();
		b3[i+2*n] = laplacianCoordinate.e3();
	}

	TRversor R2 = handles[1]->GetTRVersor();
	TRversor R2inv = inverse(R2);
	for(std::set<int>::iterator citer = fconstraints.begin(); citer != fconstraints.end() ; citer++)
	{
		int i = *citer;
		auto constraint = normalize(_point(R2 * vertexDescriptors[i]->positionConstrained * R2inv));
		b3[i] = constraint.e1();
		b3[i+n] = constraint.e2();
		b3[i+2*n] = constraint.e3();
	}
	TRversor R1 = handles[0]->GetTRVersor();
	TRversor R1inv = inverse(R1);
	for(std::set<int>::iterator citer = constraints.begin(); citer != constraints.end() ; citer++)
	{
		int i = *citer;
		auto constraint = normalize(_point(R1 * vertexDescriptors[i]->positionConstrained * R1inv));
		b3[i] = constraint.e1();
		b3[i+n] = constraint.e2();
		b3[i+2*n] = constraint.e3();
	}

	SolveATA(Lc, b3, xyz, 3);

	for( int i = 0 ; i < n ; ++i )
	{
		vertexDescriptors[i]->deformedPosition = _vectorE3GA(xyz[i], xyz[i+n], xyz[i+2*n]);
		vertexDescriptors[i]->normal = vertexDescriptors[i]->normalOrig;
	}
}

void SolveLinearSystemTaucsHi(vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptors, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptorsLow)
{
	int n = vertexDescriptors.size();

	for (int i = 0; i < n; ++i)
	{
		auto &laplacianCoordinate = vertexDescriptors[i]->laplacianCoordinate;
		b3Hi[i] = laplacianCoordinate.e1();
		b3Hi[i + n] = laplacianCoordinate.e2();
		b3Hi[i + 2 * n] = laplacianCoordinate.e3();
	}

	for (map<int, int>::iterator iter = mapping.begin(); iter != mapping.end(); iter++){
		int i = iter->second;
		auto &constraint = vertexDescriptorsLow[iter->first]->deformedPosition;
		b3Hi[i] = constraint.e1();
		b3Hi[i + n] = constraint.e2();
		b3Hi[i + 2 * n] = constraint.e3();
	}

	SolveATA(LcHi, b3Hi, xyzHi, 3);

	for (int i = 0; i < n; ++i)
	{
		vertexDescriptors[i]->deformedPosition = _vectorE3GA(xyzHi[i], xyzHi[i + n], xyzHi[i + 2 * n]);
		vertexDescriptors[i]->normal = vertexDescriptors[i]->normalOrig;
	}
}

rotor E3GA(const vector<vectorE3GA>& P, const vector<vectorE3GA>& Q, double wij, const vector<rotor>& Rq, const rotor& Rx, const int Nz)
{
	wij = sqrt(wij);

	int N = Nz * 8 + 1;
	VectorXd  F(N);
	Matrix<double, Dynamic, 4>  J(N, 4);
	Matrix4d  JtJ;
	Vector4d  b;
	Vector4d  x;

	// Inicializamos el vector x:
	// Rotor R0 + R1 e1^e2 + R2 e1^e3 + R3 e2^e3
	//double R0 = 1.0, R1 = 0.0, R2 = 0.0, R3 = 0.0;
	double R0 = _double(Rx), R1 = Rx.e1e2(), R2 = -Rx.e3e1(), R3 = Rx.e2e3();

	J.setZero(N, 4);
	int row;// = 0;
	int M = Nz * 4;
	int K = Nz;
	for (int i = 0; i < K; ++i)
	{
		double p1 = P[i].e1(), p2 = P[i].e2(), p3 = P[i].e3();
		double q1 = Q[i].e1(), q2 = Q[i].e2(), q3 = Q[i].e3();

		double F1 = (q3 + p3) * R2 + (q2 + p2) * R1 + (p1 - q1) * R0; // e1
		double F2 = (q3 + p3) * R3 + ((-q1) - p1) * R1 + (p2 - q2) * R0; // e2
		double F3 = ((-q2) - p2) * R3 + ((-q1) - p1) * R2 + (p3 - q3) * R0; // e3
		double F4 = (p1 - q1) * R3 + (q2 - p2) * R2 + (p3 - q3) * R1; // e1 ^ (e2 ^ e3)

		row = i * 4;

		F(row + 0) = -F1;
		F(row + 1) = -F2;
		F(row + 2) = -F3;
		F(row + 3) = -F4;

		// Calculamos las derivadas parciales de F0 con respecto de los parametros R0, R1, R2, R3
		//dF/dR0
		J(row + 0,0) = (p1 - q1); // e1
		J(row + 1,0) = (p2 - q2); // e2
		J(row + 2,0) = (p3 - q3); // e3
		//dF/dR1
		J(row + 0,1) = (q2 + p2); // e1
		J(row + 1,1) = ((-q1) - p1); // e2
		J(row + 3,1) = (p3 - q3); // e1 ^ (e2 ^ e3)
		//dF/dR2		
		J(row + 0,2) = (q3 + p3); // e1
		J(row + 2,2) = ((-q1) - p1); // e3
		J(row + 3,2) = (q2 - p2); // e1 ^ (e2 ^ e3)
		//dF/dR3		
		J(row + 1,3) = (q3 + p3); // e2
		J(row + 2,3) = ((-q2) - p2); // e3
		J(row + 3,3) = (p1 - q1); // e1 ^ (e2 ^ e3)
		//row += 4;

		row += M;

		double Rj0 = _double(Rq[i]), Rj1 = Rq[i].e1e2(), Rj2 = -Rq[i].e3e1(), Rj3 = Rq[i].e2e3();

		F1 = Rj0 - R0; // 1.0
		F2 = Rj1 - R1; // e1 ^ e2
		F3 = Rj2 - R2; // e1 ^ e3
		F4 = Rj3 - R3; // e2 ^ e3

		F(row + 0) = -wij * F1;
		F(row + 1) = -wij * F2;
		F(row + 2) = -wij * F3;
		F(row + 3) = -wij * F4;

		//dF/dR0
		J(row + 0, 0) = wij * -1; // 1.0
		//dF/dR1
		J(row + 1, 1) = wij * -1; // e1 ^ e2
		//dF/dR2
		J(row + 2, 2) = wij * -1; // e1 ^ e3
		//dF/dR3
		J(row + 3, 3) = wij * -1; // e2 ^ e3
	}

	M *= 2;
	F(M) = -((R3 * R3 + R2 * R2 + R1 * R1 + R0 * R0) - 1.0); // 1.0

	// Calculamos las derivadas parciales de Rot con respecto de los parametros R0, R1, R2 y R3
	//dF/dR0
	J(M,0) = 2.0 * R0; // 1.0
	//dF/dR1
	J(M,1) = 2.0 * R1; // 1.0
	//dF/dR2
	J(M,2) = 2.0 * R2; // 1.0
	//dF/dR3
	J(M,3) = 2.0 * R3; // 1.0

	JtJ = J.transpose() * J;
	b = J.transpose() * F;
	x = JtJ.inverse() * b;

	R0 += x(0);
	R1 += x(1);
	R2 += x(2);
	R3 += x(3);

	rotor R = rotor(rotor_scalar_e1e2_e2e3_e3e1, R0, R1, R3, -R2);
	return R;
}

bool matInverseSymmetric4x4(const Matrix4d &r, Matrix4d &q)
{
	double	a = r(0, 0);
	double	e = r(1, 0), f = r(1, 1);
	double	i = r(2, 0), j = r(2, 1), k = r(2, 2);
	double	m = r(3, 0), n = r(3, 1), o = r(3, 2), p = r(3, 3);

	double kpoo = k * p - o * o;
	double jpon = j * p - o * n;
	double jokn = j * o - k * n;
	double ipom = i * p - o * m;
	double iokm = i * o - k * m;
	double injm = i * n - j * m;
	double enfm = e * n - f * m;
	double q00 = (f * kpoo - j * jpon + n * jokn);
	double q10 = (e * kpoo - j * ipom + n * iokm);
	double q20 = (e * jpon - f * ipom + n * injm);
	double q30 = (e * jokn - f * iokm + j * injm);
	double	determinant = (a * q00 - e * q10 + i * q20 - m * q30);

	if (determinant == 0.0)
		return false;

	determinant = 1.0 / determinant;

	q(0, 0) = determinant * q00;
	q(1, 0) = -determinant * q10;
	q(2, 0) = determinant * q20;
	q(3, 0) = -determinant * q30;

	q(0, 1) = q(1, 0);
	q(1, 1) = determinant * (a * kpoo - i * ipom + m * iokm);
	q(2, 1) = -determinant * (a * jpon - e * ipom + m * injm);
	q(3, 1) = determinant * (a * jokn - e * iokm + i * injm);

	q(0, 2) = q(2, 0);
	q(1, 2) = q(2, 1);
	q(2, 2) = determinant * (a * (f * p - n * n) - e * (e * p - n * m) + m * enfm);
	q(3, 2) = -determinant * (a * (f * o - j * n) - e * (e * o - j * m) + i * enfm);

	q(0, 3) = q(3, 0);
	q(1, 3) = q(3, 1);
	q(2, 3) = q(3, 2);
	q(3, 3) = determinant * (a * (f * k - j * j) - e * (e * k - j * i) + i * (e * j - f * i));

	return true;
}

rotor E3GA_Fast3(const vector<vectorE3GA>& P, const vector<vectorE3GA>& Q, double wij, const vector<rotor>& Rq, const rotor& Rx, const int N)
{
	//wij = sqrt(wij);

	Vector4d F;
	Matrix4d JtJ;
	Vector4d b;
	Vector4d x;

	// Inicializamos el vector x:
	// Rotor R0 + R1 e1^e2 + R2 e1^e3 + R3 e2^e3
	//double R0 = 1.0, R1 = 0.0, R2 = 0.0, R3 = 0.0;
	double R0 = _double(Rx), R1 = Rx.e1e2(), R2 = -Rx.e3e1(), R3 = Rx.e2e3();

	Matrix4d _Jt;

	JtJ.setZero();
	_Jt.setZero();
	double wij2 = wij;// *wij;
	//int N = P.size();

	// Computamos F0(x)
	b.setZero();
	for (int i = 0; i < N; ++i)
	{
		double p1 = P[i].e1(), p2 = P[i].e2(), p3 = P[i].e3();
		double q1 = Q[i].e1(), q2 = Q[i].e2(), q3 = Q[i].e3();

		// R p - q R = vector and trivector
		double F1 = (q3 + p3) * R2 + (q2 + p2) * R1 + (p1 - q1) * R0; // e1
		double F2 = (q3 + p3) * R3 + ((-q1) - p1) * R1 + (p2 - q2) * R0; // e2
		double F3 = ((-q2) - p2) * R3 + ((-q1) - p1) * R2 + (p3 - q3) * R0; // e3
		double F4 = (p1 - q1) * R3 + (q2 - p2) * R2 + (p3 - q3) * R1; // e1 ^ (e2 ^ e3)

		F(0) = -F1;
		F(1) = -F2;
		F(2) = -F3;
		F(3) = -F4;

		// Calculamos las derivadas parciales de F0 con respecto de los parametros R0, R1, R2, R3

		//dF/dR0
		_Jt(0, 0) = (p1 - q1); // e1
		_Jt(0, 1) = (p2 - q2); // e2
		_Jt(0, 2) = (p3 - q3); // e3
		//dF/dR1
		_Jt(1, 0) = (q2 + p2); // e1
		_Jt(1, 1) = ((-q1) - p1); // e2
		_Jt(1, 3) = (p3 - q3); // e1 ^ (e2 ^ e3)
		//dF/dR2		
		_Jt(2, 0) = (q3 + p3); // e1
		_Jt(2, 2) = ((-q1) - p1); // e3
		_Jt(2, 3) = (q2 - p2); // e1 ^ (e2 ^ e3)
		//dF/dR3		
		_Jt(3, 1) = (q3 + p3); // e2
		_Jt(3, 2) = ((-q2) - p2); // e3
		_Jt(3, 3) = (p1 - q1); // e1 ^ (e2 ^ e3)

		b.noalias() += _Jt * F;

		JtJ.noalias() += _Jt * _Jt.transpose();


		// Calculamos las derivadas parciales de F1 con respecto de los parametros R0, R1, R2, R3
		JtJ(0, 0) += wij2;
		JtJ(1, 1) += wij2;
		JtJ(2, 2) += wij2;
		JtJ(3, 3) += wij2;

		double Rj0 = _double(Rq[i]), Rj1 = Rq[i].e1e2(), Rj2 = -Rq[i].e3e1(), Rj3 = Rq[i].e2e3();

		F(0) = Rj0 - R0; // 1.0
		F(1) = Rj1 - R1; // e1 ^ e2
		F(2) = Rj2 - R2; // e1 ^ e3
		F(3) = Rj3 - R3; // e2 ^ e3

		b.noalias() += wij2 * F;
		//b(0) += wij2 * F1;
		//b(1) += wij2 * F2;
		//b(2) += wij2 * F3;
		//b(3) += wij2 * F4;
	}

	//double RF = -((R3 * R3 + R2 * R2 + R1 * R1 + R0 * R0) - 1.0); // R ~R - 1.0   

	// compute partial derivatives of RF: we use F column but should use a vector J
	//F(0) = 2.0 * R0;
	//F(1) = 2.0 * R1;
	//F(2) = 2.0 * R2;
	//F(3) = 2.0 * R3;

	//JtJ.noalias() += F * F.transpose();
	//
	//b.noalias() += F * RF;

	x.noalias() = JtJ.inverse() * b;

	R0 += x(0);
	R1 += x(1);
	R2 += x(2);
	R3 += x(3);

	return rotor(rotor_scalar_e1e2_e2e3_e3e1, R0, R1, R3, -R2);
}

rotor E3GA_Fast2(const vector<vectorE3GA>& P, const vector<vectorE3GA>& Q, double wij, const vector<rotor>& Rq, const rotor& Rx, const int N)
{
	//wij = sqrt(wij);

	Vector4d F;
	Matrix4d JtJ;
	Vector4d b;
	Vector4d x;

	// Inicializamos el vector x:
	// Rotor R0 + R1 e1^e2 + R2 e1^e3 + R3 e2^e3
	//double R0 = 1.0, R1 = 0.0, R2 = 0.0, R3 = 0.0;
	double R0 = _double(Rx), R1 = Rx.e1e2(), R2 = -Rx.e3e1(), R3 = Rx.e2e3();

	Matrix4d _Jt;
	Vector4d u;

	JtJ.setZero();
	_Jt.setZero();
	double wij2 = wij;// *wij;

	// Computamos F0(x)
	b.setZero();
	for (int i = 0; i < N; ++i)
	{
		double p1 = P[i].e1(), p2 = P[i].e2(), p3 = P[i].e3();
		double q1 = Q[i].e1(), q2 = Q[i].e2(), q3 = Q[i].e3();

		double F1 = (q3 + p3) * R2 + (q2 + p2) * R1 + (p1 - q1) * R0; // e1
		double F2 = (q3 + p3) * R3 + ((-q1) - p1) * R1 + (p2 - q2) * R0; // e2
		double F3 = ((-q2) - p2) * R3 + ((-q1) - p1) * R2 + (p3 - q3) * R0; // e3
		double F4 = (p1 - q1) * R3 + (q2 - p2) * R2 + (p3 - q3) * R1; // e1 ^ (e2 ^ e3)

		F(0) = -F1;
		F(1) = -F2;
		F(2) = -F3;
		F(3) = -F4;

		// Calculamos las derivadas parciales de F0 con respecto de los parametros R0, R1, R2, R3

		//dF/dR0
		_Jt(0, 0) = (p1 - q1); // e1
		_Jt(0, 1) = (p2 - q2); // e2
		_Jt(0, 2) = (p3 - q3); // e3
		//dF/dR1
		_Jt(1, 0) = (q2 + p2); // e1
		_Jt(1, 1) = ((-q1) - p1); // e2
		_Jt(1, 3) = (p3 - q3); // e1 ^ (e2 ^ e3)
		//dF/dR2		
		_Jt(2, 0) = (q3 + p3); // e1
		_Jt(2, 2) = ((-q1) - p1); // e3
		_Jt(2, 3) = (q2 - p2); // e1 ^ (e2 ^ e3)
		//dF/dR3		
		_Jt(3, 1) = (q3 + p3); // e2
		_Jt(3, 2) = ((-q2) - p2); // e3
		_Jt(3, 3) = (p1 - q1); // e1 ^ (e2 ^ e3)

		b.noalias() += _Jt * F;

		JtJ.noalias() += _Jt * _Jt.transpose();

		// Calculamos las derivadas parciales de F1 con respecto de los parametros R0, R1, R2, R3
		JtJ(0, 0) += wij2;
		JtJ(1, 1) += wij2;
		JtJ(2, 2) += wij2;
		JtJ(3, 3) += wij2;

		double Rj0 = _double(Rq[i]), Rj1 = Rq[i].e1e2(), Rj2 = -Rq[i].e3e1(), Rj3 = Rq[i].e2e3();

		F(0) = Rj0 - R0; // 1.0
		F(1) = Rj1 - R1; // e1 ^ e2
		F(2) = Rj2 - R2; // e1 ^ e3
		F(3) = Rj3 - R3; // e2 ^ e3

		b.noalias() += wij2 * F;
		//b(0) += wij2 * F1;
		//b(1) += wij2 * F2;
		//b(2) += wij2 * F3;
		//b(3) += wij2 * F4;
	}

	matInverseSymmetric4x4(JtJ, JtJ);
	//JtJ = JtJ.inverse();

	double RF = -((R3 * R3 + R2 * R2 + R1 * R1 + R0 * R0) - 1.0); // 1.0

	// compute partial derivatives of RF: we use F column but should use a vector J
	F(0) = 2.0 * R0;
	F(1) = 2.0 * R1;
	F(2) = 2.0 * R2;
	F(3) = 2.0 * R3;

	//Sherman–Morrison formula
	u.noalias() = JtJ * F;

	double k = 1.0 / (1.0 + u(0, 0) * F(0, 0) + u(1, 0) * F(1, 0) + u(2, 0) * F(2, 0) + u(3, 0) * F(3, 0));

	JtJ(0, 0) = JtJ(0, 0) - k * u(0, 0) * u(0, 0);
	JtJ(0, 1) = JtJ(0, 1) - k * u(0, 0) * u(1, 0);
	JtJ(0, 2) = JtJ(0, 2) - k * u(0, 0) * u(2, 0);
	JtJ(0, 3) = JtJ(0, 3) - k * u(0, 0) * u(3, 0);
	JtJ(1, 0) = JtJ(0, 1);
	JtJ(1, 1) = JtJ(1, 1) - k * u(1, 0) * u(1, 0);
	JtJ(1, 2) = JtJ(1, 2) - k * u(1, 0) * u(2, 0);
	JtJ(1, 3) = JtJ(1, 3) - k * u(1, 0) * u(3, 0);
	JtJ(2, 0) = JtJ(0, 2);
	JtJ(2, 1) = JtJ(1, 2);
	JtJ(2, 2) = JtJ(2, 2) - k * u(2, 0) * u(2, 0);
	JtJ(2, 3) = JtJ(2, 3) - k * u(2, 0) * u(3, 0);
	JtJ(3, 0) = JtJ(0, 3);
	JtJ(3, 1) = JtJ(1, 3);
	JtJ(3, 2) = JtJ(2, 3);
	JtJ(3, 3) = JtJ(3, 3) - k * u(3, 0) * u(3, 0);

	b.noalias() += F * RF;

	x.noalias() = JtJ * b;

	R0 += x(0);
	R1 += x(1);
	R2 += x(2);
	R3 += x(3);

	rotor R = rotor(rotor_scalar_e1e2_e2e3_e3e1, R0, R1, R3, -R2);
	return R;
}

rotor E3GA_Fast4(const vector<vectorE3GA>& P, const vector<vectorE3GA>& Q, double wij, const vector<rotor>& Rq, const rotor& Rx, const int N)
{
	//wij = sqrt(wij);

	Vector4d F;
	F(3) = 0; // e1 ^ (e2 ^ e3)
	Matrix4d JtJ;
	Vector4d b;
	Vector4d x;

	// Inicializamos el vector x:
	// Rotor R0 + R1 e1^e2 + R2 e1^e3 + R3 e2^e3
	//const double R0 = 1.0, R1 = 0.0, R2 = 0.0, R3 = 0.0;
	//double R0 = _double(Rx), R1 = Rx.e1e2(), R2 = -Rx.e3e1(), R3 = Rx.e2e3();

	Matrix4d _Jt;
	//Vector4d u;

	JtJ.setZero();
	//_Jt.setZero();
	_Jt(0, 3) = 0;
	_Jt(1, 2) = 0;
	_Jt(2, 1) = 0;
	_Jt(3, 0) = 0;
	const double wij2 = wij;// *wij;

	double q3p3;
	double q2p2;
	double p1q1;
	double q1p1;
	double p2q2;
	double p3q3;

	// Computamos F0(x)
	b.setZero();
	for (int i = 0; i < N; ++i)
	{
		const vectorE3GA& Pi = P[i];
		const vectorE3GA& Qi = Q[i];

		q1p1 = Qi.e1() + Pi.e1();
		q2p2 = Qi.e2() + Pi.e2();
		q3p3 = Qi.e3() + Pi.e3();
		p1q1 = Pi.e1() - Qi.e1();
		p2q2 = Pi.e2() - Qi.e2();
		p3q3 = Pi.e3() - Qi.e3();

		F(0) = -p1q1; // e1
		F(1) = -p2q2; // e2
		F(2) = -p3q3; // e3

		// Calculamos las derivadas parciales de F0 con respecto de los parametros R0, R1, R2, R3

		//dF/dR0
		_Jt(0, 0) = p1q1; // e1
		_Jt(0, 1) = p2q2; // e2
		_Jt(0, 2) = p3q3; // e3
		//dF/dR1
		_Jt(1, 0) = q2p2; // e1
		_Jt(1, 1) = -q1p1; // e2
		_Jt(1, 3) = p3q3; // e1 ^ (e2 ^ e3)
		//dF/dR2		
		_Jt(2, 0) = q3p3; // e1
		_Jt(2, 2) = -q1p1; // e3
		_Jt(2, 3) = -p2q2; // e1 ^ (e2 ^ e3)
		//dF/dR3		
		_Jt(3, 1) = q3p3; // e2
		_Jt(3, 2) = -q2p2; // e3
		_Jt(3, 3) = p1q1; // e1 ^ (e2 ^ e3)

		b.noalias() += _Jt * F;

		JtJ.noalias() += _Jt * _Jt.transpose();

		b(0) += wij2 * (_double(Rq[i]) - 1.0); // 1.0
		b(1) += wij2 * Rq[i].e1e2(); // e1 ^ e2
		b(2) += wij2 * -Rq[i].e3e1(); // e1 ^ e3
		b(3) += wij2 * Rq[i].e2e3(); // e2 ^ e3
	}

	const double nwij2 = N * wij2 + 1e-6;

	JtJ(0, 0) += nwij2;
	JtJ(1, 1) += nwij2;
	JtJ(2, 2) += nwij2;
	JtJ(3, 3) += nwij2;

	b(0) += 1e-6 * (_double(Rx) - 1.0); // 1.0
	b(1) += 1e-6 * Rx.e1e2(); // e1 ^ e2
	b(2) += 1e-6 * -Rx.e3e1(); // e1 ^ e3
	b(3) += 1e-6 * Rx.e2e3(); // e2 ^ e3

	matInverseSymmetric4x4(JtJ, JtJ);

	x.noalias() = JtJ * b;

	return rotor(rotor_scalar_e1e2_e2e3_e3e1, 1.0 + x(0), x(1), x(3), -x(2));
}

void E3GA_Prep(const vector<vectorE3GA>& P, const vector<vectorE3GA>& Q, double wij2, const int N, Matrix4d &JtJ, Vector4d &b)
{
	JtJ.setZero();

	double q3p3;
	double q2p2;
	double p1q1;
	double q1p1;
	double p2q2;
	double p3q3;
	double p1q1_2;
	double p2q2_2;
	double p3q3_2;
	double q2p2_2;
	double q1p1_2;
	double q3p3_2;

	for (int i = 0; i < N; ++i)
	{
		const vectorE3GA& Pi = P[i];
		const vectorE3GA& Qi = Q[i];

		q1p1 = Qi.e1() + Pi.e1();
		q2p2 = Qi.e2() + Pi.e2();
		q3p3 = Qi.e3() + Pi.e3();
		p1q1 = Pi.e1() - Qi.e1();
		p2q2 = Pi.e2() - Qi.e2();
		p3q3 = Pi.e3() - Qi.e3();

		p1q1_2 = p1q1 * p1q1;
		p2q2_2 = p2q2 * p2q2;
		p3q3_2 = p3q3 * p3q3;
		q2p2_2 = q2p2 * q2p2;
		q1p1_2 = q1p1 * q1p1;
		q3p3_2 = q3p3 * q3p3;

		JtJ(0, 0) += p1q1_2 + p2q2_2 + p3q3_2;
		JtJ(0, 1) += p1q1 * q2p2 - p2q2 * q1p1;
		JtJ(0, 2) += p1q1 * q3p3 - p3q3 * q1p1;
		JtJ(0, 3) += p2q2 * q3p3 - p3q3 * q2p2;
		JtJ(1, 1) += q2p2_2 + q1p1_2 + p3q3_2;
		JtJ(1, 2) += q2p2 * q3p3 - p3q3 * p2q2;
		JtJ(1, 3) += -q1p1 * q3p3 + p3q3 * p1q1;
		JtJ(2, 2) += q3p3_2 + q1p1_2 + p2q2_2;
		JtJ(2, 3) += q1p1 * q2p2 - p2q2 * p1q1;
		JtJ(3, 3) += q3p3_2 + q2p2_2 + p1q1_2;
	}

	b(0) = -JtJ(0, 0);
	b(1) = -JtJ(0, 1);
	b(2) = -JtJ(0, 2);
	b(3) = -JtJ(0, 3);

	const double nwij2 = N * wij2 + 1e-6;

	JtJ(0, 0) += nwij2;
	JtJ(1, 1) += nwij2;
	JtJ(2, 2) += nwij2;
	JtJ(3, 3) += nwij2;

	JtJ(1, 0) = JtJ(0, 1);
	JtJ(2, 0) = JtJ(0, 2);
	JtJ(2, 1) = JtJ(1, 2);
	JtJ(3, 0) = JtJ(0, 3);
	JtJ(3, 1) = JtJ(1, 3);
	JtJ(3, 2) = JtJ(2, 3);

	matInverseSymmetric4x4(JtJ, JtJ);
}

rotor E3GA_Fast5(double wij, const vector<rotor>& Rq, const rotor& Rx, const int N, const Matrix4d &JtJ, const Vector4d &b)
{
	Vector4d x;
	Vector4d b2;
	b2.setZero();
	//const double wij2 = wij;// *wij;

	for (int i = 0; i < N; ++i)
	{
		const rotor& Rqi = Rq[i];
		b2(0) += (_double(Rqi) - 1.0); // 1.0
		b2(1) += Rqi.e1e2(); // e1 ^ e2
		b2(2) += -Rqi.e3e1(); // e1 ^ e3
		b2(3) += Rqi.e2e3(); // e2 ^ e3
	}

	b2 *= wij;

	b2(0) += b(0) + 1e-6 * (_double(Rx) - 1.0); // 1.0
	b2(1) += b(1) + 1e-6 * Rx.e1e2(); // e1 ^ e2
	b2(2) += b(2) + 1e-6 * -Rx.e3e1(); // e1 ^ e3
	b2(3) += b(3) + 1e-6 * Rx.e2e3(); // e2 ^ e3

	x.noalias() = JtJ * b2;

	return rotor(rotor_scalar_e1e2_e2e3_e3e1, 1.0 + x(0), x(1), x(3), -x(2));
}

rotor E3GA4(const vector<vectorE3GA>& P, const vector<vectorE3GA>& Q, const rotor& Rx, const int N)
{
	//wij = sqrt(wij);

	Vector4d F;
	F(3) = 0; // e1 ^ (e2 ^ e3)
	Matrix4d JtJ;
	Vector4d b;
	Vector4d x;

	// Inicializamos el vector x:
	// Rotor R0 + R1 e1^e2 + R2 e1^e3 + R3 e2^e3
	//const double R0 = 1.0, R1 = 0.0, R2 = 0.0, R3 = 0.0;
	//double R0 = _double(Rx), R1 = Rx.e1e2(), R2 = -Rx.e3e1(), R3 = Rx.e2e3();

	Matrix4d _Jt;
	//Vector4d u;

	JtJ.setZero();
	//_Jt.setZero();
	_Jt(0, 3) = 0;
	_Jt(1, 2) = 0;
	_Jt(2, 1) = 0;
	_Jt(3, 0) = 0;

	double q3p3;
	double q2p2;
	double p1q1;
	double q1p1;
	double p2q2;
	double p3q3;
	double p1, q1, p2, q2, p3, q3;

	// Computamos F0(x)
	b.setZero();
	for (int i = 0; i < N; ++i)
	{
		p1 = P[i].e1(), p2 = P[i].e2(), p3 = P[i].e3();
		q1 = Q[i].e1(), q2 = Q[i].e2(), q3 = Q[i].e3();

		q3p3 = q3 + p3;
		q2p2 = q2 + p2;
		p1q1 = p1 - q1;
		q1p1 = q1 + p1;
		p2q2 = p2 - q2;
		p3q3 = p3 - q3;

		F(0) = -p1q1; // e1
		F(1) = -p2q2; // e2
		F(2) = -p3q3; // e3

		// Calculamos las derivadas parciales de F0 con respecto de los parametros R0, R1, R2, R3

		//dF/dR0
		_Jt(0, 0) = p1q1; // e1
		_Jt(0, 1) = p2q2; // e2
		_Jt(0, 2) = p3q3; // e3
		//dF/dR1
		_Jt(1, 0) = q2p2; // e1
		_Jt(1, 1) = -q1p1; // e2
		_Jt(1, 3) = p3q3; // e1 ^ (e2 ^ e3)
		//dF/dR2		
		_Jt(2, 0) = q3p3; // e1
		_Jt(2, 2) = -q1p1; // e3
		_Jt(2, 3) = -p2q2; // e1 ^ (e2 ^ e3)
		//dF/dR3		
		_Jt(3, 1) = q3p3; // e2
		_Jt(3, 2) = -q2p2; // e3
		_Jt(3, 3) = p1q1; // e1 ^ (e2 ^ e3)

		b.noalias() += _Jt * F;

		JtJ.noalias() += _Jt * _Jt.transpose();
	}

	const double wij2 = 1e-6;

	b(0) += wij2 * (_double(Rx) - 1.0); // 1.0
	b(1) += wij2 * Rx.e1e2(); // e1 ^ e2
	b(2) += wij2 * -Rx.e3e1(); // e1 ^ e3
	b(3) += wij2 * Rx.e2e3(); // e2 ^ e3

	JtJ(0, 0) += wij2;
	JtJ(1, 1) += wij2;
	JtJ(2, 2) += wij2;
	JtJ(3, 3) += wij2;

	matInverseSymmetric4x4(JtJ, JtJ);

	x.noalias() = JtJ * b;

	return rotor(rotor_scalar_e1e2_e2e3_e3e1, 1.0 + x(0), x(1), x(3), -x(2));
}

rotor E3GA5(const vector<vectorE3GA>& P, const vector<vectorE3GA>& Q, const int N)
{
	//wij = sqrt(wij);

	Vector4d F;
	F(3) = 0; // e1 ^ (e2 ^ e3)
	Matrix4d JtJ;
	Vector4d b;
	Vector4d b2;
	Vector4d x;

	// Inicializamos el vector x:
	// Rotor R0 + R1 e1^e2 + R2 e1^e3 + R3 e2^e3
	//const double R0 = 1.0, R1 = 0.0, R2 = 0.0, R3 = 0.0;
	//double R0 = _double(Rx), R1 = Rx.e1e2(), R2 = -Rx.e3e1(), R3 = Rx.e2e3();
	rotor Rx = rotor(rotor_scalar_e1e2_e2e3_e3e1, 1.0, 0, 0, 0);

	Matrix4d _Jt;
	//Vector4d u;

	//_Jt.setZero();
	_Jt(0, 3) = 0;
	_Jt(1, 2) = 0;
	_Jt(2, 1) = 0;
	_Jt(3, 0) = 0;

	double q3p3;
	double q2p2;
	double p1q1;
	double q1p1;
	double p2q2;
	double p3q3;
	double p1, q1, p2, q2, p3, q3;

	// Computamos F0(x)
	JtJ.setZero();
	b.setZero();
	for (int i = 0; i < N; ++i)
	{
		p1 = P[i].e1(), p2 = P[i].e2(), p3 = P[i].e3();
		q1 = Q[i].e1(), q2 = Q[i].e2(), q3 = Q[i].e3();

		q3p3 = q3 + p3;
		q2p2 = q2 + p2;
		p1q1 = p1 - q1;
		q1p1 = q1 + p1;
		p2q2 = p2 - q2;
		p3q3 = p3 - q3;

		F(0) = -p1q1; // e1
		F(1) = -p2q2; // e2
		F(2) = -p3q3; // e3

		// Calculamos las derivadas parciales de F0 con respecto de los parametros R0, R1, R2, R3

		//dF/dR0
		_Jt(0, 0) = p1q1; // e1
		_Jt(0, 1) = p2q2; // e2
		_Jt(0, 2) = p3q3; // e3
		//dF/dR1
		_Jt(1, 0) = q2p2; // e1
		_Jt(1, 1) = -q1p1; // e2
		_Jt(1, 3) = p3q3; // e1 ^ (e2 ^ e3)
		//dF/dR2		
		_Jt(2, 0) = q3p3; // e1
		_Jt(2, 2) = -q1p1; // e3
		_Jt(2, 3) = -p2q2; // e1 ^ (e2 ^ e3)
		//dF/dR3		
		_Jt(3, 1) = q3p3; // e2
		_Jt(3, 2) = -q2p2; // e3
		_Jt(3, 3) = p1q1; // e1 ^ (e2 ^ e3)

		b.noalias() += _Jt * F;

		JtJ.noalias() += _Jt * _Jt.transpose();
	}

	const double wij2 = 1e-6;
	JtJ(0, 0) += wij2;
	JtJ(1, 1) += wij2;
	JtJ(2, 2) += wij2;
	JtJ(3, 3) += wij2;
	matInverseSymmetric4x4(JtJ, JtJ);

	rotor newRx;
	double err = 0;
	do {
		b2(0) = b(0) + wij2 * (_double(Rx) - 1.0); // 1.0
		b2(1) = b(1) + wij2 * Rx.e1e2(); // e1 ^ e2
		b2(2) = b(2) + wij2 * -Rx.e3e1(); // e1 ^ e3
		b2(3) = b(3) + wij2 * Rx.e2e3(); // e2 ^ e3

		x.noalias() = JtJ * b2;

		newRx = normalize(rotor(rotor_scalar_e1e2_e2e3_e3e1, 1.0 + x(0), x(1), x(3), -x(2)));
		err = _double(c3ga::norm_e(_rotor(Rx - newRx)));
		Rx = newRx;
	} while (err > 1e-4);
	return Rx;
}

rotor E3GA_Fast(const vector<vectorE3GA>& P, const vector<vectorE3GA>& Q, double wij, const vector<rotor>& Rq, const rotor& Rx, const int N)
{
	wij = sqrt(wij);

	Vector4d F;
	vector<Matrix4d, Eigen::aligned_allocator<Matrix4d>> Jt;
	Matrix4d JtJ;
	Matrix4d JtJ2;
	Vector4d b;
	Vector4d x;
	double E;

	// Inicializamos el vector x:
	// Rotor R0 + R1 e1^e2 + R2 e1^e3 + R3 e2^e3
	//double R0 = 1.0, R1 = 0.0, R2 = 0.0, R3 = 0.0;
	double R0 = _double(Rx), R1 = Rx.e1e2(), R2 = -Rx.e3e1(), R3 = Rx.e2e3();

	Matrix4d _Jt;
	Vector4d u;

	JtJ.setZero();
	_Jt.setZero();
	Jt.reserve(P.size());
	double wij2 = wij*wij;
	for (int i = 0; i < N; ++i)
	{
		double p1 = P[i].e1(), p2 = P[i].e2(), p3 = P[i].e3();
		double q1 = Q[i].e1(), q2 = Q[i].e2(), q3 = Q[i].e3();

		// Calculamos las derivadas parciales de F0 con respecto de los parametros R0, R1, R2, R3

		//dF/dR0
		_Jt(0, 0) = (p1 - q1); // e1
		_Jt(0, 1) = (p2 - q2); // e2
		_Jt(0, 2) = (p3 - q3); // e3
		//dF/dR1
		_Jt(1, 0) = (q2 + p2); // e1
		_Jt(1, 1) = ((-q1) - p1); // e2
		_Jt(1, 3) = (p3 - q3); // e1 ^ (e2 ^ e3)
		//dF/dR2		
		_Jt(2, 0) = (q3 + p3); // e1
		_Jt(2, 2) = ((-q1) - p1); // e3
		_Jt(2, 3) = (q2 - p2); // e1 ^ (e2 ^ e3)
		//dF/dR3		
		_Jt(3, 1) = (q3 + p3); // e2
		_Jt(3, 2) = ((-q2) - p2); // e3
		_Jt(3, 3) = (p1 - q1); // e1 ^ (e2 ^ e3)

		Jt.push_back(_Jt);
		JtJ += _Jt * _Jt.transpose();

		// Calculamos las derivadas parciales de F1 con respecto de los parametros R0, R1, R2, R3
		JtJ(0, 0) += wij2;
		JtJ(1, 1) += wij2;
		JtJ(2, 2) += wij2;
		JtJ(3, 3) += wij2;
	}

	matInverseSymmetric4x4(JtJ, JtJ);
	//JtJ = JtJ.inverse();

	//_J.ZeroMatrix();
	double Energy = 0;
	double Error = 1e38;
	int count = 10;
	while (fabs(Error) > 1e-4 && count-- > 0)
	{
		// Computamos F0(x)
		b.setZero();
		E = 0.0;
		for (int i = 0; i < N; ++i)
		{
			double p1 = P[i].e1(), p2 = P[i].e2(), p3 = P[i].e3();
			double q1 = Q[i].e1(), q2 = Q[i].e2(), q3 = Q[i].e3();

			double F1 = (q3 + p3) * R2 + (q2 + p2) * R1 + (p1 - q1) * R0; // e1
			double F2 = (q3 + p3) * R3 + ((-q1) - p1) * R1 + (p2 - q2) * R0; // e2
			double F3 = ((-q2) - p2) * R3 + ((-q1) - p1) * R2 + (p3 - q3) * R0; // e3
			double F4 = (p1 - q1) * R3 + (q2 - p2) * R2 + (p3 - q3) * R1; // e1 ^ (e2 ^ e3)

			F(0) = -F1;
			F(1) = -F2;
			F(2) = -F3;
			F(3) = -F4;

			E += F1*F1 + F2*F2 + F3*F3 + F4*F4;
			b += Jt[i] * F;

			double Rj0 = _double(Rq[i]), Rj1 = Rq[i].e1e2(), Rj2 = -Rq[i].e3e1(), Rj3 = Rq[i].e2e3();

			F1 = Rj0 - R0; // 1.0
			F2 = Rj1 - R1; // e1 ^ e2
			F3 = Rj2 - R2; // e1 ^ e3
			F4 = Rj3 - R3; // e2 ^ e3

			b(0) += wij2 * F1;
			b(1) += wij2 * F2;
			b(2) += wij2 * F3;
			b(3) += wij2 * F4;

			E += wij2 * (F1*F1 + F2*F2 + F3*F3 + F4*F4);
		}

		double RF = -((R3 * R3 + R2 * R2 + R1 * R1 + R0 * R0) - 1.0); // 1.0

		E += RF * RF;

		// compute partial derivatives of RF: we use F column but should use a vector J
		F(0) = 2.0 * R0;
		F(1) = 2.0 * R1;
		F(2) = 2.0 * R2;
		F(3) = 2.0 * R3;

		//Sherman–Morrison formula
		u = JtJ * F;

		double k = 1.0 / (1.0 + u(0, 0) * F(0, 0) + u(1, 0) * F(1, 0) + u(2, 0) * F(2, 0) + u(3, 0) * F(3, 0));

		JtJ2(0, 0) = JtJ(0, 0) - k * u(0, 0) * u(0, 0);
		JtJ2(0, 1) = JtJ(0, 1) - k * u(0, 0) * u(1, 0);
		JtJ2(0, 2) = JtJ(0, 2) - k * u(0, 0) * u(2, 0);
		JtJ2(0, 3) = JtJ(0, 3) - k * u(0, 0) * u(3, 0);
		JtJ2(1, 0) = JtJ2(0, 1);
		JtJ2(1, 1) = JtJ(1, 1) - k * u(1, 0) * u(1, 0);
		JtJ2(1, 2) = JtJ(1, 2) - k * u(1, 0) * u(2, 0);
		JtJ2(1, 3) = JtJ(1, 3) - k * u(1, 0) * u(3, 0);
		JtJ2(2, 0) = JtJ2(0, 2);
		JtJ2(2, 1) = JtJ2(1, 2);
		JtJ2(2, 2) = JtJ(2, 2) - k * u(2, 0) * u(2, 0);
		JtJ2(2, 3) = JtJ(2, 3) - k * u(2, 0) * u(3, 0);
		JtJ2(3, 0) = JtJ2(0, 3);
		JtJ2(3, 1) = JtJ2(1, 3);
		JtJ2(3, 2) = JtJ2(2, 3);
		JtJ2(3, 3) = JtJ(3, 3) - k * u(3, 0) * u(3, 0);

		b += F * RF;

		x = JtJ2 * b;

		R0 += x(0);
		R1 += x(1);
		R2 += x(2);
		R3 += x(3);

		Error = Energy - E;
		Energy = E;
	}

	rotor R = rotor(rotor_scalar_e1e2_e2e3_e3e1, R0, R1, R3, -R2);
	return R;
}

rotor E3GA(const vector<vectorE3GA>& P, const vector<vectorE3GA>& Q, const rotor& Rx, const int N)
{
	Vector4d F;
	Matrix4d Jt;
	Matrix4d JtJ;
	Vector4d b;
	Vector4d x;

	// Inicializamos el vector x:
	// Rotor R0 + R1 e1^e2 + R2 e1^e3 + R3 e2^e3
	//double R0 = 1.0, R1 = 0.0, R2 = 0.0, R3 = 0.0;
	double R0 = _double(Rx), R1 = Rx.e1e2(), R2 = -Rx.e3e1(), R3 = Rx.e2e3();

	// Computamos en Jacobiano (el Jacobiano es sparse pero aqui lo tenemos denso por ser pequeño)
	Jt.setZero();
	JtJ.setZero();
	b.setZero();
	//int N = P.size();
	// Computamos F(x)
	for (int i = 0; i < N; ++i)
	{
		double p1 = P[i].e1(), p2 = P[i].e2(), p3 = P[i].e3();
		double q1 = Q[i].e1(), q2 = Q[i].e2(), q3 = Q[i].e3();

		double F1 = (q3 + p3) * R2 + (q2 + p2) * R1 + (p1 - q1) * R0; // e1
		double F2 = (q3 + p3) * R3 + ((-q1) - p1) * R1 + (p2 - q2) * R0; // e2
		double F3 = ((-q2) - p2) * R3 + ((-q1) - p1) * R2 + (p3 - q3) * R0; // e3
		double F4 = (p1 - q1) * R3 + (q2 - p2) * R2 + (p3 - q3) * R1; // e1 ^ (e2 ^ e3)

		F(0) = -F1;
		F(1) = -F2;
		F(2) = -F3;
		F(3) = -F4;

		// Calculamos las derivadas parciales de F con respecto de los parametros R0, R1, R2, R3
		//dF/dR0
		Jt(0, 0) = (p1 - q1); // e1
		Jt(0, 1) = (p2 - q2); // e2
		Jt(0, 2) = (p3 - q3); // e3
		//dF/dR1
		Jt(1, 0) = (q2 + p2); // e1
		Jt(1, 1) = ((-q1) - p1); // e2
		Jt(1, 3) = (p3 - q3); // e1 ^ (e2 ^ e3)
		//dF/dR2
		Jt(2, 0) = (q3 + p3); // e1
		Jt(2, 2) = ((-q1) - p1); // e3
		Jt(2, 3) = (q2 - p2); // e1 ^ (e2 ^ e3)
		//dF/dR3
		Jt(3, 1) = (q3 + p3); // e2
		Jt(3, 2) = ((-q2) - p2); // e3
		Jt(3, 3) = (p1 - q1); // e1 ^ (e2 ^ e3)

		JtJ += Jt * Jt.transpose();
		b += Jt * F;

	}

	double RF = -((R3 * R3 + R2 * R2 + R1 * R1 + R0 * R0) - 1.0); // 1.0

	// Calculamos las derivadas parciales de Rot con respecto de los parametros R0, R1, R2 y R3
	//dF/dR0
	F(0) = 2.0 * R0; // 1.0
	//dF/dR1
	F(1) = 2.0 * R1; // 1.0
	//dF/dR2
	F(2) = 2.0 * R2; // 1.0
	//dF/dR3
	F(3) = 2.0 * R3; // 1.0

	JtJ += F * F.transpose();
	b += F * RF;

	//JtJ(0, 0) += 1e-6;
	//JtJ(1, 1) += 1e-6;
	//JtJ(2, 2) += 1e-6;
	//JtJ(3, 3) += 1e-6;

	x = JtJ.inverse() * b;

	R0 += x(0);
	R1 += x(1);
	R2 += x(2);
	R3 += x(3);

	rotor R = rotor(rotor_scalar_e1e2_e2e3_e3e1, R0, R1, R3, -R2);
	return R;
}

rotor E3GA2(const vector<vectorE3GA>& P, const vector<vectorE3GA>& Q, double wij, const vector<rotor>& Rq, const int N)
{
	wij = sqrt(wij);

	MatrixNxM F = MatrixNxM(N * 8 + 1, 1);
	MatrixNxM J = MatrixNxM(N * 8 + 1, 4);
	MatrixNxM JtJ = MatrixNxM(4, 4);
	MatrixNxM b = MatrixNxM(4, 1);
	MatrixNxM x = MatrixNxM(4, 1);
	MatrixNxM E = MatrixNxM(1, 1);

	// Inicializamos el vector x:
	// Rotor R0 + R1 e1^e2 + R2 e1^e3 + R3 e2^e3
	double R0 = 1.0, R1 = 0.0, R2 = 0.0, R3 = 0.0;

	double Energy = 0;
	double Error = 1e38;
	int count = 10;
	while (fabs(Error) > 1e-4 && count-- > 0)
	{
		// Computamos F0(x)
		int row = 0;
		for (int i = 0; i < N; ++i)
		{
			double p1 = P[i].e1(), p2 = P[i].e2(), p3 = P[i].e3();
			double q1 = Q[i].e1(), q2 = Q[i].e2(), q3 = Q[i].e3();

			double F1 = (q3 + p3) * R2 + (q2 + p2) * R1 + (p1 - q1) * R0; // e1
			double F2 = (q3 + p3) * R3 + ((-q1) - p1) * R1 + (p2 - q2) * R0; // e2
			double F3 = ((-q2) - p2) * R3 + ((-q1) - p1) * R2 + (p3 - q3) * R0; // e3
			double F4 = (p1 - q1) * R3 + (q2 - p2) * R2 + (p3 - q3) * R1; // e1 ^ (e2 ^ e3)

			F[row + 0][0] = -F1;
			F[row + 1][0] = -F2;
			F[row + 2][0] = -F3;
			F[row + 3][0] = -F4;
			row += 4;
		}

		// Computamos F1(x)
		for (int i = 0; i < N; ++i)
		{
			double Rj0 = _double(Rq[i]), Rj1 = Rq[i].e1e2(), Rj2 = -Rq[i].e3e1(), Rj3 = Rq[i].e2e3();

			double F1 = (Rj3 * R3 + Rj2 * R2 + Rj1 * R1 + Rj0 * R0) - 1.0; // 1.0
			double F2 = ((-(Rj2 * R3)) + Rj3 * R2) - Rj0 * R1 + Rj1 * R0; // e1 ^ e2
			double F3 = Rj1 * R3 - Rj0 * R2 - Rj3 * R1 + Rj2 * R0; // e1 ^ e3
			double F4 = (-(Rj0 * R3)) - Rj1 * R2 + Rj2 * R1 + Rj3 * R0; // e2 ^ e3

			F[row + 0][0] = -wij * F1;
			F[row + 1][0] = -wij * F2;
			F[row + 2][0] = -wij * F3;
			F[row + 3][0] = -wij * F4;
			row += 4;
		}

		F[row][0] = -((R3 * R3 + R2 * R2 + R1 * R1 + R0 * R0) - 1.0); // 1.0

		// Calculamos las derivadas parciales de Rot con respecto de los parametros R0, R1, R2 y R3
		J.ZeroMatrix();
		row = 0;
		for (int i = 0; i < N; ++i)
		{
			double p1 = P[i].e1(), p2 = P[i].e2(), p3 = P[i].e3();
			double q1 = Q[i].e1(), q2 = Q[i].e2(), q3 = Q[i].e3();

			// Calculamos las derivadas parciales de F0 con respecto de los parametros R0, R1, R2, R3

			//dF/dR0
			J[row + 0][0] = (p1 - q1); // e1
			J[row + 1][0] = (p2 - q2); // e2
			J[row + 2][0] = (p3 - q3); // e3
			//dF/dR1
			J[row + 0][1] = (q2 + p2); // e1
			J[row + 1][1] = ((-q1) - p1); // e2
			J[row + 3][1] = (p3 - q3); // e1 ^ (e2 ^ e3)
			//dF/dR2		
			J[row + 0][2] = (q3 + p3); // e1
			J[row + 2][2] = ((-q1) - p1); // e3
			J[row + 3][2] = (q2 - p2); // e1 ^ (e2 ^ e3)
			//dF/dR3		
			J[row + 1][3] = (q3 + p3); // e2
			J[row + 2][3] = ((-q2) - p2); // e3
			J[row + 3][3] = (p1 - q1); // e1 ^ (e2 ^ e3)
			row += 4;
		}

		for (int i = 0; i < N; ++i)
		{
			// Calculamos las derivadas parciales de F1 con respecto de los parametros R0, R1, R2, R3
			double Rj0 = _double(Rq[i]), Rj1 = Rq[i].e1e2(), Rj2 = -Rq[i].e3e1(), Rj3 = Rq[i].e2e3();

			//dF/dR0
			J[row + 0][0] = wij * Rj0; // 1.0
			J[row + 1][0] = wij * Rj1; // e1 ^ e2
			J[row + 2][0] = wij * Rj2; // e1 ^ e3
			J[row + 3][0] = wij * Rj3; // e2 ^ e3
			//dF/dR1
			J[row + 0][1] = wij * Rj1; // 1.0
			J[row + 1][1] = wij * (-Rj0); // e1 ^ e2
			J[row + 2][1] = wij * (-Rj3); // e1 ^ e3
			J[row + 3][1] = wij * Rj2; // e2 ^ e3
			//dF/dR2
			J[row + 0][2] = wij * Rj2; // 1.0
			J[row + 1][2] = wij * Rj3; // e1 ^ e2
			J[row + 2][2] = wij * (-Rj0); // e1 ^ e3
			J[row + 3][2] = wij * (-Rj1); // e2 ^ e3
			//dF/dR3
			J[row + 0][3] = wij * Rj3; // 1.0
			J[row + 1][3] = wij * (-Rj2); // e1 ^ e2
			J[row + 2][3] = wij * Rj1; // e1 ^ e3
			J[row + 3][3] = wij * (-Rj0); // e2 ^ e3

			row += 4;
		}

		//dF/dR0
		J[P.size() * 8][0] = 2.0 * R0; // 1.0
		//dF/dR1
		J[P.size() * 8][1] = 2.0 * R1; // 1.0
		//dF/dR2
		J[P.size() * 8][2] = 2.0 * R2; // 1.0
		//dF/dR3
		J[P.size() * 8][3] = 2.0 * R3; // 1.0

		TransposedMult(J, J, &JtJ);
		TransposedMult(J, F, &b);
		SolveLU(JtJ, b, &x);

		R0 += x[0][0];
		R1 += x[1][0];
		R2 += x[2][0];
		R3 += x[3][0];

		TransposedMult(F, F, &E);
		Error = Energy - E[0][0];
		Energy = E[0][0];
	}

	rotor R = rotor(rotor_scalar_e1e2_e2e3_e3e1, R0, R1, R3, -R2);
	return R;
}

/*
Exp_approx = {
	1 + _P(1) + _P(1)*_P(1)/2 + _P(1)*_P(1)*_P(1)/6
}
Mi = R0*(1) + R1*(e1^e2) + R2*(e1^e3) + R3*(e2^e3);
twist = B1*(e1^e2) + B2*(e1^e3) + B3*(e2^e3);
?R = Exp_approx(-0.5*twist);
?M = R * Mi;
*/
void calculate2(double B1, double B2, double B3, double &M0, double &M1, double &M2, double &M3)
{
	double R0, R1, R2, R3;
	double X0 = M0, X1 = M1, X2 = M2, X3 = M3;
	R0 = (-(0.125 * B3 * B3)) - 0.125 * B2 * B2 - 0.125 * B1 * B1 + 1.0; // 1.0
	R1 = (0.02083333395421505 * B1 * B3 * B3 + 0.02083333395421505 * B1 * B2 * B2 + 0.02083333395421505 * B1 * B1 * B1) - B1 / 2.0; // e1 ^ e2
	R2 = 0.02083333395421505 * B2 * B3 * B3 + 0.02083333395421505 * B2 * B2 * B2 + (0.02083333395421505 * B1 * B1 - 0.5) * B2; // e1 ^ e3
	R3 = 0.02083333395421505 * B3 * B3 * B3 + ((0.02083333395421505 * B2 * B2 + 0.02083333395421505 * B1 * B1) - 0.5) * B3; // e2 ^ e3
	M0 = (-(R3 * X3)) - R2 * X2 - R1 * X1 + R0 * X0; // 1.0
	M1 = (-(R2 * X3)) + R3 * X2 + R0 * X1 + R1 * X0; // e1 ^ e2
	M2 = (R1 * X3 + R0 * X2) - R3 * X1 + R2 * X0; // e1 ^ e3
	M3 = R0 * X3 - R1 * X2 + R2 * X1 + R3 * X0; // e2 ^ e3
}

rotor LieAlgebraGradientDescent(const vector<vectorE3GA>& P, const vector<vectorE3GA>& Q, const rotor& Rx, const int N)
{
	/*
	Constraint equation:

	(R P ~R) - Q = 0

	Approximating R P ~R to first order using taylor series of e^-B/2 we get:

	R P ~R = (e^-B/2) P (e^B/2) ~= (1 - B/2) P (1 + B/2)
	~= P + 0.5 (P B - B P) - (0.5)^2 B P B
	~= P + P x B + O(B)^2

	dropping the squared term we have:

	R P ~R ~= P + P x B

	Replacing in original equation:

	(P + P x B) - Q = 0
	P x B = Q - P

	That's a LINEAR equation in B.

	[F] [B] = [-P + Q]

	where [F] is the 3x3 matrix of the LHS of equation (1) and [-P x Q] is the 3x1 matrix of the RHS of eq (1) and [B] is the 3x1 matrix solution.
	So having three or more pairs of points we can solve for [B] using Linear Least Squares (LLS).

	That first order aproximation can be accumulated ITERATIVELY in order to get the best approximation.
	*/

	Matrix3d U;
	Vector3d Ub;
	Matrix3d A;
	Vector3d b;

	//rotor M = Rx;
	double R0 = _double(Rx), R1 = Rx.e1e2(), R2 = -Rx.e3e1(), R3 = Rx.e2e3();

	//rotor deltaM;
	//vectorE3GA Pi;

	//dualLine_e1e2_e1e3_e2e3_e1ni_e2ni_e3ni
	//static bivectorE3GA basis[] = { _bivectorE3GA(e1^e2), _bivectorE3GA(e1^e3), _bivectorE3GA(e2^e3) };

	A.setZero();
	double p1, p2, p3, X1, X2, X3;

	U.setZero();
	Ub.setZero();

	for (int i = 0; i < N; ++i)
	{
		const vectorE3GA &Pi = P[i];
		p1 = Pi.e1();
		p2 = Pi.e2();
		p3 = Pi.e3();

		//Xi = M * P[i] * (~M);
		X1 = ((p1 * R3 * R3 + (2.0 * p3 * R1 - 2.0 * p2 * R2) * R3) - p1 * R2 * R2 + 2.0 * p3 * R0 * R2) - p1 * R1 * R1 + 2.0 * p2 * R0 * R1 + p1 * R0 * R0; // e1
		X2 = ((-(p2 * R3 * R3)) + (2.0 * p3 * R0 - 2.0 * p1 * R2) * R3 + p2 * R2 * R2) - 2.0 * p3 * R1 * R2 - p2 * R1 * R1 - 2.0 * p1 * R0 * R1 + p2 * R0 * R0; // e2
		X3 = ((-(p3 * R3 * R3)) + (2.0 * p1 * R1 - 2.0 * p2 * R0) * R3) - p3 * R2 * R2 + ((-(2.0 * p2 * R1)) - 2.0 * p1 * R0) * R2 + p3 * R1 * R1 + p3 * R0 * R0; // e3

		//for (int j = 0; j < 3; ++j)
		//{
		//	vectorE3GA result = Commutator(Pi, basis[j]);
		//	A(0, j) = result.e1();  // e1
		//	A(1, j) = result.e2();  // e2
		//	A(2, j) = result.e3();  // e3
		//}
		A(0, 0) = -X2; // e1
		A(1, 0) = X1; // e2
		A(0, 1) = -X3; // e1
		A(2, 1) = X1; // e3
		A(1, 2) = -X3; // e2
		A(2, 2) = X2; // e3


		const vectorE3GA &Qi = Q[i];
		b(0) = Qi.e1() - X1;  // e1
		b(1) = Qi.e2() - X2;  // e2
		b(2) = Qi.e3() - X3;  // e3

		U.noalias() += A.transpose() * A;
		Ub.noalias() += A.transpose() * b;
	}

	b.noalias() = U.inverse() * Ub;

	calculate2(b(0), b(1), b(2), R0, R1, R2, R3);
	//bivectorE3GA B(__bivectorE3GA_coordinates__::bivectorE3GA_e1e2_e2e3_e3e1, b(0), b(2), -b(1));
	//rotor deltaM = exp(-0.5*B);
	//return deltaM * Rx;

	return rotor(__rotor_coordinates__::rotor_scalar_e1e2_e2e3_e3e1, R0, R1, R3, -R2);
}

void UpdateLaplaciansRotationPrep(Mesh *mesh, std::shared_ptr<SparseMatrix> A, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptors)
{
	//Matrix3x3 m;
	//Matrix3x3 U, W, V, Ut;
	//Matrix3x3 M;
	//ICP icp;

	const double alpha = 0.14;
	concurrency::parallel_for_each(mesh->vertices.begin(), mesh->vertices.end(), [&](Vertex * vertex)
	//for (auto vertex : mesh->vertices)
	{
		vector<vectorE3GA> P;
		vector<vectorE3GA> Q;

		P.resize(32);
		Q.resize(32);
		int i = vertex->ID;

		Vector3 &pi = vertex->p;
		vectorE3GA &tpi = vertexDescriptors[i]->deformedPosition;

		int vertexDegree = 0;
		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex->iterator(); vertexDegree < 32 && !edgeAroundIter.end(); edgeAroundIter++, vertexDegree++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;

			vectorE3GA &tpj = vertexDescriptors[j]->deformedPosition;
			Vector3 &pj = mesh->vertexAt(j)->p;
			Vector3 eij = (pj - pi) * (*A)(i, j);
			vectorE3GA teij = tpj - tpi;

			P[vertexDegree] = _vectorE3GA(eij.x(), eij.y(), eij.z());
			Q[vertexDegree] = teij;
		}
		double invVertexDegree = alpha * g_meshArea / (double)vertexDegree;
		//rotor M = E3GA4(P, Q, vertexDescriptors[i]->M, vertexDegree);
		//rotor M = E3GA5(P, Q, vertexDegree);
		E3GA_Prep(P, Q, invVertexDegree, vertexDegree, vertexDescriptors[i]->JtJ, vertexDescriptors[i]->b);
	//}
	});
}

void UpdateLaplaciansRotationExec(Mesh *mesh, std::shared_ptr<SparseMatrix> A, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptors, bool updateLaplacian)
{
	const double alpha = 0.14;
	concurrency::parallel_for_each(mesh->vertices.begin(), mesh->vertices.end(), [&](Vertex * vertex)
	//for (auto vertex : mesh->vertices)
	{
		vector<rotor> Rq;

		Rq.resize(32);
		int i = vertex->ID;

		Vector3 &pi = vertex->p;
		vectorE3GA &tpi = vertexDescriptors[i]->deformedPosition;

		int vertexDegree = 0;
		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex->iterator(); vertexDegree < 32 && !edgeAroundIter.end(); edgeAroundIter++, vertexDegree++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;
			Rq[vertexDegree] = vertexDescriptors[j]->M;
		}
		double invVertexDegree = alpha * g_meshArea / (double)vertexDegree;
		//rotor M = E3GA4(P, Q, vertexDescriptors[i]->M, vertexDegree);
		//rotor M = E3GA5(P, Q, vertexDegree);
		rotor M = E3GA_Fast5(invVertexDegree, Rq, vertexDescriptors[i]->M, vertexDegree, vertexDescriptors[i]->JtJ, vertexDescriptors[i]->b);
		vertexDescriptors[i]->M = normalize(M);
	//}
	});

	//for (int i = 0; i < vertexDescriptors.size(); ++i)
	//{
	//	vertexDescriptors[i]->M = vertexDescriptors[i]->R;
	//}

	if (!updateLaplacian)
	{
		return;
	}

	concurrency::parallel_for_each(mesh->vertices.begin(), mesh->vertices.end(), [&](Vertex * vertex)
	//for (auto vertex : mesh->vertices)
	{
		int i = vertex->ID;
		vertexDescriptors[i]->laplacianCoordinate = _vectorE3GA(0.0, 0.0, 0.0);
		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex->iterator(); !edgeAroundIter.end(); edgeAroundIter++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;
			double wij = (*A)(i, j);
			rotor &Ri = vertexDescriptors[i]->M;
			rotor &Rj = vertexDescriptors[j]->M;
			vectorE3GA V = _vectorE3GA(vertexDescriptors[j]->position) - _vectorE3GA(vertexDescriptors[i]->position);
			vectorE3GA Vp = _vectorE3GA(0.5 * wij * (Ri * V * (~Ri) + Rj * V * (~Rj)));
			vertexDescriptors[i]->laplacianCoordinate += Vp;
		}
	//}
	});
}

void UpdateLaplaciansRotation(Mesh *mesh, std::shared_ptr<SparseMatrix> A, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptors, bool updateLaplacian)
{
	//Matrix3x3 m;
	//Matrix3x3 U, W, V, Ut;
	//Matrix3x3 M;
	//ICP icp;

	vector<vectorE3GA> P;
	vector<vectorE3GA> Q;
	vector<rotor> Rq;

	P.resize(32);
	Q.resize(32);
	Rq.resize(32);

	const double alpha = 0.08;
	//concurrency::parallel_for_each(mesh->vertices.begin(), mesh->vertices.end(), [&](Vertex * vertex)
	for (auto vertex : mesh->vertices)
	{
		int i = vertex->ID;

		Vector3 &pi = vertex->p;
		vectorE3GA &tpi = vertexDescriptors[i]->deformedPosition;

		int vertexDegree = 0;
		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex->iterator(); vertexDegree < 32 && !edgeAroundIter.end(); edgeAroundIter++, vertexDegree++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;

			vectorE3GA &tpj = vertexDescriptors[j]->deformedPosition;
			Vector3 &pj = mesh->vertexAt(j)->p;
			Vector3 eij = (pj - pi) * (*A)(i, j);
			vectorE3GA teij = tpj - tpi;

			P[vertexDegree] = _vectorE3GA(eij.x(), eij.y(), eij.z());
			Q[vertexDegree] = teij;
			Rq[vertexDegree] = vertexDescriptors[j]->M;
		}
		double invVertexDegree = alpha * g_meshArea / (double)vertexDegree;
		//rotor M = E3GA4(P, Q, vertexDescriptors[i]->M, vertexDegree);
		//rotor M = E3GA5(P, Q, vertexDegree);
		rotor M = E3GA_Fast4(P, Q, invVertexDegree, Rq, vertexDescriptors[i]->M, vertexDegree);
		vertexDescriptors[i]->M = normalize(M);
	}//);

	//for (int i = 0; i < vertexDescriptors.size(); ++i)
	//{
	//	vertexDescriptors[i]->M = vertexDescriptors[i]->R;
	//}

	if (!updateLaplacian)
	{
		return;
	}

	//concurrency::parallel_for_each(mesh->vertices.begin(), mesh->vertices.end(), [&](Vertex * vertex)
	for (auto vertex : mesh->vertices)
	{
		int i = vertex->ID;
		vertexDescriptors[i]->laplacianCoordinate = _vectorE3GA(0.0, 0.0, 0.0);
		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex->iterator(); !edgeAroundIter.end(); edgeAroundIter++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;
			double wij = (*A)(i, j);
			rotor &Ri = vertexDescriptors[i]->M;
			rotor &Rj = vertexDescriptors[j]->M;
			vectorE3GA V = _vectorE3GA(vertexDescriptors[j]->position) - _vectorE3GA(vertexDescriptors[i]->position);
			vectorE3GA Vp = _vectorE3GA(0.5 * wij * (Ri * V * (~Ri) + Rj * V * (~Rj)));
			vertexDescriptors[i]->laplacianCoordinate += Vp;
		}
	}//);
}

double ComputeEnergy(Mesh *mesh, std::shared_ptr<SparseMatrix> A, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptors)
{
	double arapEnergy = 0.0;
	const double alpha = 0.08;
	//concurrency::parallel_for_each(mesh->vertices.begin(), mesh->vertices.end(), [&](Vertex * vertex)
	for (auto vertex : mesh->vertices)
	{
		int i = vertex->ID;

		Vector3 &pi = vertex->p;
		vectorE3GA &tpi = vertexDescriptors[i]->deformedPosition;
		auto& Ri = vertexDescriptors[i]->M;

		int vertexDegree = 0;
		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex->iterator(); !edgeAroundIter.end(); edgeAroundIter++)
		{
			vertexDegree++;
		}

		double invVertexDegree = alpha * g_meshArea / (double)vertexDegree;

		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex->iterator(); !edgeAroundIter.end(); edgeAroundIter++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;

			vectorE3GA &tpj = vertexDescriptors[j]->deformedPosition;
			Vector3 &pj = mesh->vertexAt(j)->p;

			Vector3 eij = (pj - pi);
			auto q_ij = tpj - tpi;
			auto p_ij = _vectorE3GA(eij.x(), eij.y(), eij.z());
			auto& Rj = vertexDescriptors[j]->M;

			arapEnergy += _double((*A)(i, j) * c3ga::norm_e(_vectorE3GA(Ri * p_ij * ~Ri) - q_ij));
			arapEnergy += _double(invVertexDegree * c3ga::norm_e(_rotor(Rj - Ri)));
		}
	}//);

	return arapEnergy;
}

vectorE3GA Color( double d )
{
	static vectorE3GA	c0 = _vectorE3GA( 1, 1, 1);
	static vectorE3GA	c1 = _vectorE3GA( 1, 1, 0);
	static vectorE3GA	c2 = _vectorE3GA( 0, 1, 0);
	static vectorE3GA	c3 = _vectorE3GA( 0, 1, 1);
	static vectorE3GA	c4 = _vectorE3GA( 0, 0, 1);

	if( d < 0.25 )
	{
		double alpha = (d - 0.0) / (0.25-0.0);
		return (1.0 - alpha) * c0 + alpha * c1;
	}
	else if( d < 0.5 )
	{
		double alpha = (d - 0.25) / (0.5-0.25);
		return (1.0 - alpha) * c1 + alpha * c2;
	}
	else if( d < 0.75 )
	{
		double alpha = (d - 0.5) / (0.75-0.5);
		return (1.0 - alpha) * c2 + alpha * c3;
	}
	else
	{
		double alpha = (d - 0.75) / (1.0-0.75);
		return (1.0 - alpha) * c3 + alpha * c4;
	}
}

void display()
{
	/*
	 *	matrices
	 */
	glViewport( 0, 0, g_viewportWidth, g_viewportHeight );
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	pickLoadMatrix();
	GLpick::g_frustumFar = 1000.0;
	GLpick::g_frustumNear = .1;
	gluPerspective( 60.0, (double)g_viewportWidth/(double)g_viewportHeight, GLpick::g_frustumNear, GLpick::g_frustumFar );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glShadeModel(GL_SMOOTH);	//gouraud shading
	glClearDepth(1.0f);
	glClearColor( .75f, .75f, .75f, .0f );
	glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );

	/*
	 *	estados
	 */
	glEnable(GL_CULL_FACE);		//face culling
	glCullFace( GL_BACK );
	glFrontFace( GL_CCW );
	glEnable(GL_DEPTH_TEST);	//z-buffer
	glDepthFunc(GL_LEQUAL);

	/*
	 *	iluminacion
	 */
	float		ambient[] = { .3f, .3f, .3f, 1.f };
	float		diffuse[] = { .3f, .3f, .3f, 1.f };
	float		position[] = { .0f, 0.f, -150.f, 1.f };
	float		specular[] = { 1.f, 1.f, 1.f };

	glLightfv( GL_LIGHT0, GL_AMBIENT, ambient );
	glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuse );
	glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 0);
	glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.0125);
	glEnable(  GL_LIGHT0   );
	glEnable(  GL_LIGHTING );
	glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, specular );
	glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, 50.f );

	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glLoadIdentity();

	g_camera.glLookAt();

	glLightfv( GL_LIGHT0, /*GL_POSITION*/GL_SPOT_DIRECTION, position );

	glPushMatrix();

	rotorGLMult(g_modelRotor);

	benchs.Begin();

	if(!g_computeBasis)
	{
		static bool oneTime = true;
		if (oneTime || (g_rotateKeyRotors || g_translateKeyRotors || g_automaticAnimation))
		{
			if(oneTime == true)
			{
				SolveLinearSystemTaucs(vertexDescriptorsLow);
			}

			oneTime = false;

			for(int i = 0 ; i < 3 ; ++i)
			{
				UpdateLaplaciansRotationPrep(&meshLow, A, vertexDescriptorsLow);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, true);
				SolveLinearSystemTaucs(vertexDescriptorsLow);
			}
			transferRotations(&mesh, AHi, vertexDescriptors, vertexDescriptorsLow);
			SolveLinearSystemTaucsHi(vertexDescriptors, vertexDescriptorsLow);
			//UpdateLaplaciansRotationPrep(&mesh, AHi, vertexDescriptors);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, true);
			//SolveLinearSystemTaucsHi(vertexDescriptors, vertexDescriptorsLow);
		}
	}
	if(g_iterateManyTimes)
	{
		if (g_convergence)
		{
			clock_t begin;
			clock_t end;
			double avgTime = 0.0;

			auto arapEnergy = ComputeEnergy(&meshLow, A, vertexDescriptorsLow);
			int count = 0;
			for (; count < 2500; ++count)
			{
				begin = clock();

				UpdateLaplaciansRotationPrep(&meshLow, A, vertexDescriptorsLow);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, true);
				SolveLinearSystemTaucs(vertexDescriptorsLow);

				end = clock();

				double secDeform = (double)(end - begin) / CLOCKS_PER_SEC;

				avgTime += secDeform;

				auto newArapEnergy = ComputeEnergy(&meshLow, A, vertexDescriptorsLow);
				if (fabs(arapEnergy - newArapEnergy) < 1e-4)
				{
					std::cout << "convergence after " << count + 1 << " iterations" << std::endl;
					break;
				}
				arapEnergy = newArapEnergy;
			}
			arapEnergy = ComputeEnergy(&meshLow, A, vertexDescriptorsLow);
			std::cout << arapEnergy << std::endl;
			std::cout << "Average iteration time (sec): " << avgTime / (double)(count + 1) << endl;
			std::cout << "Total accum time (sec): " << avgTime << endl;

			transferRotations(&mesh, AHi, vertexDescriptors, vertexDescriptorsLow);
			SolveLinearSystemTaucsHi(vertexDescriptors, vertexDescriptorsLow);
			//UpdateLaplaciansRotationPrep(&mesh, AHi, vertexDescriptors);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, true);
			//SolveLinearSystemTaucsHi(vertexDescriptors, vertexDescriptorsLow);
		}
		else
		{
			for (int i = 0; i < 100; ++i)
			{
				UpdateLaplaciansRotationPrep(&meshLow, A, vertexDescriptorsLow);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, false);
				UpdateLaplaciansRotationExec(&meshLow, A, vertexDescriptorsLow, true);
				SolveLinearSystemTaucs(vertexDescriptorsLow);
			}
			transferRotations(&mesh, AHi, vertexDescriptors, vertexDescriptorsLow);
			SolveLinearSystemTaucsHi(vertexDescriptors, vertexDescriptorsLow);
			//UpdateLaplaciansRotationPrep(&mesh, AHi, vertexDescriptors);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, false);
			//UpdateLaplaciansRotationExec(&mesh, AHi, vertexDescriptors, true);
			//SolveLinearSystemTaucsHi(vertexDescriptors, vertexDescriptorsLow);
		}
		g_iterateManyTimes = false;
	}

	benchs.End();

	if (GLpick::g_pickActive) glLoadName((GLuint)-1);

	double alpha = 1.0;

	//glEnable (GL_BLEND);
	//glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//alpha = 0.5;

	//Mesh-Faces Rendering
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL /*GL_LINE GL_FILL GL_POINT*/);
	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1., 1.);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable( GL_COLOR_MATERIAL );
	if (GLpick::g_pickActive) glLoadName((GLuint)10);
	for (Mesh::FaceIterator fIter = mesh.faceIterator(); !fIter.end(); fIter++)
	{
		Face		*face = fIter.face();
		int	v1 = face->edge->vertex->ID;
		int	v2 = face->edge->next->vertex->ID;
		int	v3 = face->edge->next->next->vertex->ID;
		vectorE3GA& n1 = vertexDescriptors[v1]->normal;
		vectorE3GA& n2 = vertexDescriptors[v2]->normal;
		vectorE3GA& n3 = vertexDescriptors[v3]->normal;

		vectorE3GA c0 = Color(0);
		vectorE3GA c1 = Color(0);
		vectorE3GA c2 = Color(0);

		glBegin(GL_TRIANGLES);
		glNormal3d(n1.e1(), n1.e2(), n1.e3());
		glColor4d(c0.e1(), c0.e2(), c0.e3(), alpha);
		glVertex3d(vertexDescriptors[v1]->deformedPosition.e1(), vertexDescriptors[v1]->deformedPosition.e2(), vertexDescriptors[v1]->deformedPosition.e3());
		glNormal3d(n2.e1(), n2.e2(), n2.e3());
		glColor4d(c1.e1(), c1.e2(), c1.e3(), alpha);
		glVertex3d(vertexDescriptors[v2]->deformedPosition.e1(), vertexDescriptors[v2]->deformedPosition.e2(), vertexDescriptors[v2]->deformedPosition.e3());
		glNormal3d(n3.e1(), n3.e2(), n3.e3());
		glColor4d(c2.e1(), c2.e2(), c2.e3(), alpha);
		glVertex3d(vertexDescriptors[v3]->deformedPosition.e1(), vertexDescriptors[v3]->deformedPosition.e2(), vertexDescriptors[v3]->deformedPosition.e3());
		glEnd();
	}

	//for (map<int, int>::iterator iter = mapping.begin(); iter != mapping.end(); iter++){
	//	DrawPoint(c3gaPoint(vertexDescriptors[iter->second]->deformedPosition));
	//}

	//for( Mesh::FaceIterator fIter = meshLow.faceIterator() ; !fIter.end() ; fIter++ )
	//{
	//	Face		*face = fIter.face();
	//	int	v1 = face->edge->vertex->ID;
	//	int	v2 = face->edge->next->vertex->ID;
	//	int	v3 = face->edge->next->next->vertex->ID;
	//	vectorE3GA& n1 = vertexDescriptorsLow[v1]->normal;
	//	vectorE3GA& n2 = vertexDescriptorsLow[v2]->normal;
	//	vectorE3GA& n3 = vertexDescriptorsLow[v3]->normal;

	//	vectorE3GA c0 = Color(0);
	//	vectorE3GA c1 = Color(0);
	//	vectorE3GA c2 = Color(0);

	//	glBegin (GL_TRIANGLES);
	//		glNormal3d( n1.e1(), n1.e2(), n1.e3() );
	//		glColor4d( c0.e1(), c0.e2(), c0.e3(), alpha );
	//		glVertex3d(vertexDescriptorsLow[v1]->deformedPosition.e1(), vertexDescriptorsLow[v1]->deformedPosition.e2(), vertexDescriptorsLow[v1]->deformedPosition.e3());
	//		glNormal3d( n2.e1(), n2.e2(), n2.e3() );
	//		glColor4d( c1.e1(), c1.e2(), c1.e3(), alpha );
	//		glVertex3d(vertexDescriptorsLow[v2]->deformedPosition.e1(), vertexDescriptorsLow[v2]->deformedPosition.e2(), vertexDescriptorsLow[v2]->deformedPosition.e3());
	//		glNormal3d( n3.e1(), n3.e2(), n3.e3() );
	//		glColor4d( c2.e1(), c2.e2(), c2.e3(), alpha );
	//		glVertex3d(vertexDescriptorsLow[v3]->deformedPosition.e1(), vertexDescriptorsLow[v3]->deformedPosition.e2(), vertexDescriptorsLow[v3]->deformedPosition.e3());
	//	glEnd();
	//}

	if(g_showWires)
	{
		if (!GLpick::g_pickActive)
		{
			//Mesh-Edges Rendering (superimposed to faces)
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE /*GL_LINE GL_FILL GL_POINT*/);
			glColor4d( .5, .5, .5, alpha );
			glDisable( GL_LIGHTING );
			for (Mesh::FaceIterator fIter = mesh.faceIterator(); !fIter.end(); fIter++)
			{
				Face		*face = fIter.face();
				int	v1 = face->edge->vertex->ID;
				int	v2 = face->edge->next->vertex->ID;
				int	v3 = face->edge->next->next->vertex->ID;

				glBegin(GL_TRIANGLES);
					glVertex3d(vertexDescriptors[v1]->deformedPosition.e1(), vertexDescriptors[v1]->deformedPosition.e2(), vertexDescriptors[v1]->deformedPosition.e3());
					glVertex3d(vertexDescriptors[v2]->deformedPosition.e1(), vertexDescriptors[v2]->deformedPosition.e2(), vertexDescriptors[v2]->deformedPosition.e3());
					glVertex3d(vertexDescriptors[v3]->deformedPosition.e1(), vertexDescriptors[v3]->deformedPosition.e2(), vertexDescriptors[v3]->deformedPosition.e3());
				glEnd();
			}

			//for( Mesh::FaceIterator fIter = meshLow.faceIterator() ; !fIter.end() ; fIter++ )
			//{
			//	Face		*face = fIter.face();
			//	int	v1 = face->edge->vertex->ID;
			//	int	v2 = face->edge->next->vertex->ID;
			//	int	v3 = face->edge->next->next->vertex->ID;

			//	glBegin (GL_TRIANGLES);
			//		glVertex3d(vertexDescriptorsLow[v1]->deformedPosition.e1(), vertexDescriptorsLow[v1]->deformedPosition.e2(), vertexDescriptorsLow[v1]->deformedPosition.e3());
			//		glVertex3d(vertexDescriptorsLow[v2]->deformedPosition.e1(), vertexDescriptorsLow[v2]->deformedPosition.e2(), vertexDescriptorsLow[v2]->deformedPosition.e3());
			//		glVertex3d(vertexDescriptorsLow[v3]->deformedPosition.e1(), vertexDescriptorsLow[v3]->deformedPosition.e2(), vertexDescriptorsLow[v3]->deformedPosition.e3());
			//	glEnd();
			//}
			glEnable( GL_LIGHTING );
		}
	}
	glDisable( GL_COLOR_MATERIAL );
	glDisable(GL_POLYGON_OFFSET_FILL);

	//glDisable (GL_BLEND);

	if(g_showSpheres)
	{
		//Handles rendering
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL /*GL_LINE GL_FILL GL_POINT*/);

		float	turcoise[] = { .0f, .5f, .5f, 0.3f };
		float	red[] = { .5f, .0f, .0f, 0.3f };

		for( int k = 0 ; k < handles.size() ; ++k)
		{
			if (GLpick::g_pickActive) glLoadName((GLuint)k);
			TRversor R = handles[k]->GetTRVersor();
			DrawTransparentDualSphere( _dualSphere( R * handles[k]->dS * inverse(R) ), turcoise );
		}		
	}

	glPopMatrix();

	glutSwapBuffers();
}

void reshape(GLint width, GLint height)
{
	g_viewportWidth = width;
	g_viewportHeight = height;

	// redraw viewport
	glutPostRedisplay();
}

vectorE3GA mousePosToVector(int x, int y) {
	x -= g_viewportWidth / 2;
	y -= g_viewportHeight / 2;
	return _vectorE3GA((float)-x * e1 - (float)y * e2);
}

void MouseButton(int button, int state, int x, int y)
{
	g_rotateModel = false;
	g_rotateKeyRotors = false;
	g_translateKeyRotors = false;

	if (button == GLUT_LEFT_BUTTON)
	{
		g_prevMousePos = mousePosToVector(x, y);

		GLpick::g_pickWinSize = 1;
		g_dragObject = pick(x, g_viewportHeight - y, display, &g_dragDistance);

		if(g_dragObject == -1 || g_dragObject == 10 )
		{
			vectorE3GA mousePos = mousePosToVector(x, y);
			g_rotateModel = true;

			if ((_Float(norm_e(mousePos)) / _Float(norm_e(g_viewportWidth * e1 + g_viewportHeight * e2))) < 0.2)
				g_rotateModelOutOfPlane = true;
			else g_rotateModelOutOfPlane = false;
		}
		else if(g_dragObject >= 0 && g_dragObject < handles.size())
		{
			g_rotateKeyRotors = true;
		}
	}

	if (button == GLUT_RIGHT_BUTTON)
	{
		g_prevMousePos = mousePosToVector(x, y);

		GLpick::g_pickWinSize = 1;
		g_dragObject = pick(x, g_viewportHeight - y, display, &g_dragDistance);

		if(g_dragObject >= 0 && g_dragObject < handles.size())
			g_translateKeyRotors = true;
	}
}

void MouseMotion(int x, int y)
{
	if (g_rotateModel || g_rotateKeyRotors || g_translateKeyRotors )
	{
		// get mouse position, motion
		vectorE3GA mousePos = mousePosToVector(x, y);
		vectorE3GA motion = mousePos - g_prevMousePos;

		if (g_rotateModel)
		{
			// update rotor
			if (g_rotateModelOutOfPlane)
				g_modelRotor = exp(g_camera.rotateVel * (motion ^ e3) ) * g_modelRotor;
			else 
				g_modelRotor = exp(0.00001f * (motion ^ mousePos) ) * g_modelRotor;
		}
		if(g_rotateKeyRotors)
		{
			//rotor R1 =  _rotor( inverse(g_modelRotor) * exp(-g_camera.rotateVel * (motion ^ e3) ) * g_modelRotor);
			rotor R1 =  _rotor( exp(-g_camera.rotateVel * (motion ^ e3) ) );
			if(g_dragObject < handles.size())
			{
				TRversor R = handles[g_dragObject]->R;
				handles[g_dragObject]->R = normalize(_TRversor( R1 * R  ) );
			}
		}

		if(g_translateKeyRotors)
		{
			normalizedTranslator T1 = _normalizedTranslator(inverse(g_modelRotor) * exp( _freeVector(-g_camera.translateVel*motion*ni) ) * g_modelRotor);
			if(g_dragObject < handles.size())
			{
				TRversor T = handles[g_dragObject]->T;
				handles[g_dragObject]->T = normalize(_TRversor( T1 * T ));
			}
		}

		// remember mouse pos for next motion:
		g_prevMousePos = mousePos;

		// redraw viewport
		glutPostRedisplay();
	}
}

void SpecialFunc(int key, int x, int y)
{
	switch(key) {
		case GLUT_KEY_F1 :
			{
				int mod = glutGetModifiers();
				if(mod == GLUT_ACTIVE_CTRL || mod == GLUT_ACTIVE_SHIFT )
				{
				}
			}
			break;
		case GLUT_KEY_UP:
			{
				if(g_rotateKeyRotors)
				{
					handles[g_dragObject]->dS = ChangeDualSphereRadiusSize(handles[g_dragObject]->dS, 0.025);

					// redraw viewport
					glutPostRedisplay();
				}
			}
			break;
		case GLUT_KEY_DOWN:
			{
				if(g_rotateKeyRotors)
				{
					handles[g_dragObject]->dS = ChangeDualSphereRadiusSize(handles[g_dragObject]->dS, -0.025);

					// redraw viewport
					glutPostRedisplay();
				}
			}
			break;
	}
}


void SpecialUpFunc(int key, int x, int y)
{
}

void KeyboardUpFunc(unsigned char key, int x, int y)
{
	if (key == 'c' || key == 'C')
	{
		g_convergence = !g_convergence;
		vectorE3GA motion = _vectorE3GA(-70.0, -60.0, 0);
		normalizedTranslator T1 = _normalizedTranslator(inverse(g_modelRotor) * exp(_freeVector(-g_camera.translateVel*motion*ni)) * g_modelRotor);
		TRversor T = handles[0]->T;
		handles[0]->T = normalize(_TRversor(T1 * T));

		glutPostRedisplay();
	}

	if (key == 'a' || key == 'A')
	{
		g_automaticAnimation = !g_automaticAnimation;
		if (g_automaticAnimation)
		{
			benchs.Start();
		}
		else
		{
			benchs.Stop();
		}

		glutPostRedisplay();
	}

	if(key == 'w' || key == 'W')
	{
		g_showWires = !g_showWires;
		glutPostRedisplay();
	}
	
	if( key == 'h' || key == 'H' )
	{
		g_showSpheres = !g_showSpheres;
		glutPostRedisplay();
	}

	if( key == 'x' || key == 'X' )
	{
		g_iterateManyTimes = true;
		glutPostRedisplay();
	}
}

void Idle()
{
	// redraw viewport
	if (g_automaticAnimation)
	{
		const int minCount = -50;
		const int maxCount = 50;
		static int mycount = 0;
		static int direction = 0;
		vectorE3GA motion = _vectorE3GA(1.0, 0, 0);
		if (mycount == maxCount)
		{
			direction = 1;
		}
		if (mycount == minCount)
		{
			direction = 0;
		}

		if (mycount < maxCount && direction == 0)
		{
			motion = _vectorE3GA(1.0, 0, 0);
			mycount++;
		}
		if (mycount > minCount && direction == 1)
		{
			motion = _vectorE3GA(-1.0, 0, 0);
			mycount--;
		}

		normalizedTranslator T1 = _normalizedTranslator(inverse(g_modelRotor) * exp(_freeVector(-g_camera.translateVel*motion*ni)) * g_modelRotor);
		TRversor T = handles[0]->T;
		handles[0]->T = normalize(_TRversor(T1 * T));

		if (!benchs.IsStarted())
		{
			benchs.Start();
		}

		glutPostRedisplay();
	}
}

void DestroyWindow()
{
	delete[] b3;
	delete[] xyz;
	delete[] b3Hi;
	delete[] xyzHi;
	ReleaseMatrix(Lc);
	ReleaseMatrix(LcHi);
	ReleaseDrawing();
}

