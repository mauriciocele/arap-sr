#if !defined(__GA_HELPER__H__)
#define __GA_HELPER__H__

#include "GA/c3ga.h"
#include "GA/c3ga_util.h"
#include <geometry.H>

using namespace c3ga;

typedef tangentVector dualCircle;
#define SQR(x)		((x)*(x))
#define MAX(x, y)	(((x) > (y)) ? (x) : (y))
#define MIN(x, y)	(((x) < (y)) ? (x) : (y))


inline dualCircle _dualCircle(const mv &arg1)
{
	return _tangentVector(arg1);
}

inline normalizedPoint _normalizedPoint( const normalizedFlatPoint& p)
{
	return c3gaPoint(p.e1ni(), p.e2ni(), p.e3ni());
}

inline vectorE3GA _vectorE3GA( const normalizedFlatPoint& p)
{
	return vectorE3GA(c3ga::vectorE3GA_e1_e2_e3, p.e1ni(), p.e2ni(), p.e3ni());
}

inline vectorE3GA _vectorE3GA( const double e1, const double e2, const double e3)
{
	return vectorE3GA(c3ga::vectorE3GA_e1_e2_e3, e1, e2, e3);
}

inline normalizedPoint DualSphereCenter( const dualSphere &dS )
{
	return _normalizedPoint(-0.5*(dS * ni * dS) * inverse(norm_r2(ni << dS)));
}

inline double DualSphereRadius( const dualSphere &dS)
{
	return sqrt(abs(_Float(-1*_scalar(dS * gradeInvolution(dS)) * inverse(norm_r2(ni << dS)))));
}

inline normalizedPoint Lerp(const normalizedPoint& a, const normalizedPoint& b, double t)
{
	return _normalizedPoint(-(1 - t)*(b<<ni)*a - t*(a<<ni)*b + t*(1 - t)*(a<<b)*ni);
}

inline vectorE3GA Slerp(const vectorE3GA& a, const vectorE3GA& b, double t)
{
	//double teta = angle_rad_2d(a,vecmat::vector2<double>(0,0), b);
    //double teta = acos(_double(a << b));
	double _a[3] = { a.e1(), a.e2(), a.e3() };
	double _b[3] = { b.e1(), b.e2(), b.e3() };
	double _co[3] = { 0, 0, 0 };
    double teta = angle_rad_3d( _a, _co, _b );
    return (1.0 / sin(teta)) * (a * sin((1.0 - t) * teta) + b * sin(t * teta));
}

inline dualSphere ChangeDualSphereRadiusSize(const dualSphere& ds, float delta)
{
	normalizedPoint c = DualSphereCenter(ds);
	float r = DualSphereRadius(ds);
	r += delta;
	return _dualSphere( c - 0.5 * SQR(r) * ni );
}

inline bivectorE3GA normalize(const bivectorE3GA &arg1)
{
	point p = _point(arg1);
	return _bivectorE3GA(arg1 * (1.0/_double( norm_r(arg1) )));
}

inline normalizedPoint normalize(const point &p)
{
	return _normalizedPoint(p * (1.0/_double(-ni << p)));
}

inline vectorE3GA normalize(const vectorE3GA &arg1)
{
	return _vectorE3GA(arg1 * (1.0/_double(norm_e(arg1))));
}

inline circle normalize(const circle &arg1)
{
	return _circle(arg1 * (1.0/_double(norm_r(arg1))));
}

inline normalizedTranslator normalize(const translator &arg1)
{
	return _normalizedTranslator(arg1 * (1.0/_double(norm_r(arg1))));
}

inline plane normalize(const plane &arg1)
{
	return _plane(arg1 * (1.0/_double(norm_r(arg1))));
}

inline TRversor normalize(const TRversor &arg1)
{
	return _TRversor(arg1 * (1.0/_double(norm_r(arg1))));
}

inline line normalize(const line &arg1)
{
	return _line(arg1 * (1.0/_double(norm_r(arg1))));
}

inline dualSphere normalize(const dualSphere &arg1)
{
	return _dualSphere(arg1 * (1.0/_double(norm_r(arg1))));
}

inline sphere normalize(const sphere &arg1)
{
	return _sphere(arg1 * (1.0/_double(norm_r(arg1))));
}

inline rotor normalize(const rotor &arg1)
{
	return _rotor(arg1 * (1.0/_double(norm_r(arg1))));
}

#endif
