// Flatten.cpp : Defines the entry point for the console application.
//

#ifdef WIN32
#define NOMINMAX
#include <windows.h>
#endif

#if defined (__APPLE__) || defined (OSX)
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#include "GA/c3ga.h"
#include "GA/c3ga_util.h"
#include "GA/gl_util.h"

#include "primitivedraw.h"
#include "gahelper.h"
#include "laplacian.h"

#include <memory>

#include <vector>
#include <queue>
#include <map>
#include "numerics.h"
#include "HalfEdge/Mesh.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

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
	Eigen::Matrix4d JtJ;
	Eigen::Vector4d b;
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
Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solverLow, solverHi;
int systemType = LaplaceBeltrami; //MeanValue; //LaplaceBeltrami
bool g_showSpheres = true;
bool g_showWires = false;
bool g_iterateManyTimes = false;
Eigen::MatrixXd b3Low;
Eigen::MatrixXd xyzLow;
Eigen::MatrixXd b3Hi;
Eigen::MatrixXd xyzHi;
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
		for( ; !aIter.end() ; ++aIter )
		{
			auto j = aIter.columnIndex();
			vertexDescriptors[i]->laplacianCoordinate += _vectorE3GA(vertexDescriptors[j]->position) * aIter.value();
		}
	}
}

bool is_constrained(std::set<int>& constraints, int vertex)
{
	return constraints.find(vertex) != constraints.end();
}

void PreFactor(std::shared_ptr<SparseMatrix> A, std::set<int>& constraints, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>& solver)
{

	Eigen::SparseMatrix<double> Lc = Eigen::SparseMatrix<double>(A->numRows(), A->numColumns());

	auto numRows = A->numRows();
	for (int i = 0; i < numRows; ++i)
	{
		if (!is_constrained(constraints, i))
		{
			SparseMatrix::RowIterator aIter = A->iterator(i);
			for (; !aIter.end(); ++aIter)
			{
				auto j = aIter.columnIndex();
				Lc.insert(i, j) = (*A)(i, j);
			}
		}
		else
		{
			Lc.insert(i, i) = 1.0;
		}
	}

	Lc.makeCompressed();
	solver.compute(Lc);
	if (solver.info() != Eigen::Success) {
		// TODO: error handling
	}

}

double meshArea(Mesh *mesh)
{
	double area = 0.0;
	for(Face& face :  mesh->getFaces())
	{
		Eigen::Vector3d	u = face.edge->next->vertex->p - face.edge->vertex->p;
		Eigen::Vector3d	v = face.edge->next->next->vertex->p - face.edge->vertex->p;
		area += 0.5 * u.cross(v).norm();
	}
	return area;
}

void meshMapping(Mesh *meshLow, Mesh *mesh, map<int, int>& mapping)
{
	for (Vertex& vertexLow : meshLow->getVertices())
	{
		double minDist = 1e38;
		int mappedPoint = -1;
		Eigen::Vector3d& pLow = vertexLow.p;
		for (Vertex& vertexHi : mesh->getVertices())
		{
			Eigen::Vector3d& p = vertexHi.p;
			double dist = (pLow - p).norm();
			if (dist < minDist) {
				minDist = dist;
				mappedPoint = vertexHi.ID;
			}
		}
		mapping[vertexLow.ID] = mappedPoint;
	}

	vector<bool> visited = vector<bool>(mesh->numVertices(), false);
	queue<int> assigned;

	for (auto&& pair : mapping){
		visited[pair.second] = true;
		rotorClusters[pair.second] = pair.first;
		assigned.push(pair.second);
	}

	while (!assigned.empty())
	{
		int i = assigned.front();
		assigned.pop();
		visited[i] = true;

		for (Vertex::EdgeAroundIterator edgeAroundIter = mesh->vertexAt(i).iterator(); !edgeAroundIter.end(); edgeAroundIter++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;
			if (visited[j] == false)
			{
				rotorClusters[j] = rotorClusters[i];
				assigned.push(j);
			}
		}
	}
}

int main(int argc, char* argv[])
{
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

	for( Vertex& vertexLow : meshLow.getVertices())
	{
		std::shared_ptr<VertexDescriptor> vd(new VertexDescriptor);

		vd->position = c3gaPoint(vertexLow.p.x(), vertexLow.p.y(), vertexLow.p.z() );
		vd->normalOrig = _vectorE3GA(vertexLow.n.x(), vertexLow.n.y(), vertexLow.n.z() );
		vd->isBoundary = vertexLow.isBoundary();
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
					constraints.insert(vertexLow.ID);
				}
				else
				{
					fconstraints.insert(vertexLow.ID);
				}
			}
		}
		if(vd->isBoundary)
		{
			fconstraints.insert(vertexLow.ID);
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

	b3Low = Eigen::MatrixXd(A->numRows(), 3);

	std::set<int> allconstraints = constraints;
	allconstraints.insert(fconstraints.begin(), fconstraints.end());
	PreFactor(A, allconstraints, solverLow);

	g_meshArea = meshArea(&meshLow);

	for (Vertex& vertexHi : mesh.getVertices())
	{
		std::shared_ptr<VertexDescriptor> vd(new VertexDescriptor);

		vd->position = c3gaPoint(vertexHi.p.x(), vertexHi.p.y(), vertexHi.p.z());
		vd->normalOrig = _vectorE3GA(vertexHi.n.x(), vertexHi.n.y(), vertexHi.n.z());
		vd->isBoundary = vertexHi.isBoundary();
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
	b3Hi = Eigen::MatrixXd(AHi->numRows(), 3);
	PreFactor(AHi, constraintsHi, solverHi);

	glutMainLoop();

	return 0;
}

void transferRotations(Mesh *mesh, std::shared_ptr<SparseMatrix> A, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptors, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptorsLow)
{
	for (auto&& pair : rotorClusters) {
		vertexDescriptors[pair.first]->M = vertexDescriptorsLow[pair.second]->M;
	}

	concurrency::parallel_for_each(mesh->getVertices().begin(), mesh->getVertices().end(), [&](Vertex& vertex)
	//for (auto vertex : mesh->vertices)
	{
		int i = vertex.ID;
		vertexDescriptors[i]->laplacianCoordinate = _vectorE3GA(0.0, 0.0, 0.0);
		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex.iterator(); !edgeAroundIter.end(); edgeAroundIter++)
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
		b3Low(i,0) = laplacianCoordinate.e1();
		b3Low(i,1) = laplacianCoordinate.e2();
		b3Low(i,2) = laplacianCoordinate.e3();
	}

	TRversor R2 = handles[1]->GetTRVersor();
	TRversor R2inv = inverse(R2);
	for(std::set<int>::iterator citer = fconstraints.begin(); citer != fconstraints.end() ; citer++)
	{
		int i = *citer;
		auto constraint = normalize(_point(R2 * vertexDescriptors[i]->positionConstrained * R2inv));
		b3Low(i, 0) = constraint.e1();
		b3Low(i, 1) = constraint.e2();
		b3Low(i, 2) = constraint.e3();
	}
	TRversor R1 = handles[0]->GetTRVersor();
	TRversor R1inv = inverse(R1);
	for(std::set<int>::iterator citer = constraints.begin(); citer != constraints.end() ; citer++)
	{
		int i = *citer;
		auto constraint = normalize(_point(R1 * vertexDescriptors[i]->positionConstrained * R1inv));
		b3Low(i, 0) = constraint.e1();
		b3Low(i, 1) = constraint.e2();
		b3Low(i, 2) = constraint.e3();
	}

	xyzLow = solverLow.solve(b3Low);

	for( int i = 0 ; i < n ; ++i )
	{
		vertexDescriptors[i]->deformedPosition = _vectorE3GA(xyzLow(i,0), xyzLow(i, 1), xyzLow(i, 2));
		vertexDescriptors[i]->normal = vertexDescriptors[i]->normalOrig;
	}
}

void SolveLinearSystemTaucsHi(vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptors, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptorsLow)
{
	int n = vertexDescriptors.size();

	for (int i = 0; i < n; ++i)
	{
		auto &laplacianCoordinate = vertexDescriptors[i]->laplacianCoordinate;
		b3Hi(i, 0) = laplacianCoordinate.e1();
		b3Hi(i, 1) = laplacianCoordinate.e2();
		b3Hi(i, 2) = laplacianCoordinate.e3();
	}

	for (map<int, int>::iterator iter = mapping.begin(); iter != mapping.end(); iter++){
		int i = iter->second;
		auto &constraint = vertexDescriptorsLow[iter->first]->deformedPosition;
		b3Hi(i, 0) = constraint.e1();
		b3Hi(i, 1) = constraint.e2();
		b3Hi(i, 2) = constraint.e3();
	}

	xyzHi = solverHi.solve(b3Hi);

	for (int i = 0; i < n; ++i)
	{
		vertexDescriptors[i]->deformedPosition = _vectorE3GA(xyzHi(i, 0), xyzHi(i, 1), xyzHi(i, 2));
		vertexDescriptors[i]->normal = vertexDescriptors[i]->normalOrig;
	}
}

void E3GA_Prep(const vector<vectorE3GA>& P, const vector<vectorE3GA>& Q, double wij2, const int N, Eigen::Matrix4d &JtJ, Eigen::Vector4d &b)
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

	JtJ = JtJ.inverse().eval();
}

rotor E3GA_Fast5(double wij, const vector<rotor>& Rq, const rotor& Rx, const int N, const Eigen::Matrix4d &JtJ, const Eigen::Vector4d &b)
{
	Eigen::Vector4d x;
	Eigen::Vector4d b2;
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

void UpdateLaplaciansRotationPrep(Mesh *mesh, std::shared_ptr<SparseMatrix> A, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptors)
{
	//Matrix3x3 m;
	//Matrix3x3 U, W, V, Ut;
	//Matrix3x3 M;
	//ICP icp;

	const double alpha = 0.14;
	concurrency::parallel_for_each(mesh->getVertices().begin(), mesh->getVertices().end(), [&](Vertex& vertex)
	//for (auto vertex : mesh->vertices)
	{
		vector<vectorE3GA> P;
		vector<vectorE3GA> Q;

		P.resize(32);
		Q.resize(32);
		int i = vertex.ID;

		Eigen::Vector3d &pi = vertex.p;
		vectorE3GA &tpi = vertexDescriptors[i]->deformedPosition;

		int vertexDegree = 0;
		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex.iterator(); vertexDegree < 32 && !edgeAroundIter.end(); edgeAroundIter++, vertexDegree++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;

			vectorE3GA &tpj = vertexDescriptors[j]->deformedPosition;
			Eigen::Vector3d &pj = mesh->vertexAt(j).p;
			Eigen::Vector3d eij = (pj - pi) * (*A)(i, j);
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
	concurrency::parallel_for_each(mesh->getVertices().begin(), mesh->getVertices().end(), [&](Vertex& vertex)
	//for (auto vertex : mesh->vertices)
	{
		vector<rotor> Rq;

		Rq.resize(32);
		int i = vertex.ID;

		Eigen::Vector3d &pi = vertex.p;
		vectorE3GA &tpi = vertexDescriptors[i]->deformedPosition;

		int vertexDegree = 0;
		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex.iterator(); vertexDegree < 32 && !edgeAroundIter.end(); edgeAroundIter++, vertexDegree++)
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

	concurrency::parallel_for_each(mesh->getVertices().begin(), mesh->getVertices().end(), [&](Vertex& vertex)
	//for (auto vertex : mesh->vertices)
	{
		int i = vertex.ID;
		vertexDescriptors[i]->laplacianCoordinate = _vectorE3GA(0.0, 0.0, 0.0);
		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex.iterator(); !edgeAroundIter.end(); edgeAroundIter++)
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

double ComputeEnergy(Mesh *mesh, std::shared_ptr<SparseMatrix> A, vector<std::shared_ptr<VertexDescriptor>>& vertexDescriptors)
{
	double arapEnergy = 0.0;
	const double alpha = 0.08;
	//concurrency::parallel_for_each(mesh->vertices.begin(), mesh->vertices.end(), [&](Vertex * vertex)
	for (Vertex& vertex : mesh->getVertices())
	{
		int i = vertex.ID;

		Eigen::Vector3d& pi = vertex.p;
		vectorE3GA &tpi = vertexDescriptors[i]->deformedPosition;
		auto& Ri = vertexDescriptors[i]->M;

		int vertexDegree = 0;
		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex.iterator(); !edgeAroundIter.end(); edgeAroundIter++)
		{
			vertexDegree++;
		}

		double invVertexDegree = alpha * g_meshArea / (double)vertexDegree;

		for (Vertex::EdgeAroundIterator edgeAroundIter = vertex.iterator(); !edgeAroundIter.end(); edgeAroundIter++)
		{
			int j = edgeAroundIter.edge_out()->pair->vertex->ID;

			vectorE3GA &tpj = vertexDescriptors[j]->deformedPosition;
			Eigen::Vector3d &pj = mesh->vertexAt(j).p;

			Eigen::Vector3d eij = (pj - pi);
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
	for (Face& face : mesh.getFaces())
	{
		int	v1 = face.edge->vertex->ID;
		int	v2 = face.edge->next->vertex->ID;
		int	v3 = face.edge->next->next->vertex->ID;
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
			for (Face& face : mesh.getFaces())
			{
				int	v1 = face.edge->vertex->ID;
				int	v2 = face.edge->next->vertex->ID;
				int	v3 = face.edge->next->next->vertex->ID;

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
	ReleaseDrawing();
}

