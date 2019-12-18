#if !defined(__LENGTH_CORRECTOR__H__)
#define __LENGTH_CORRECTOR__H__

#include "GA/c3ga.h"
#include "primitivedraw.h"
#include "gahelper.h"

class LengthCorrector
{
public:
	normalizedPoint A, B;
	normalizedPoint oA, oB;
	normalizedPoint Bprime;
	dualSphere poincareSphere;
	vectorE3GA center;
	double radius;
	TRSversor trs;
	double originalCurveLength;

	LengthCorrector(const normalizedPoint& _A, const normalizedPoint& _B)
	{
		oA = _A;
		oB = _B;
		A = _A;
		B = _B;
		Bprime = _B;
		vectorE3GA a = _vectorE3GA(A);
		vectorE3GA b = _vectorE3GA(B);
		center = 0.5 * (a + b);
		radius = 0.5 * _double(norm_e(b - a));
		trs = _TRSversor(exp(_freeVector(-0.5 * (center * ni)) ) * exp(log(radius)*0.5*no*ni));
		poincareSphere = _dualSphere((no - 0.5*ni));
		originalCurveLength = G(0);
		if(abs(originalCurveLength) < 1e-6 )
			originalCurveLength = 2.0;
	}

	normalizedTranslator GetCorrection(const normalizedPoint& _A, const normalizedPoint& _B)
	{
		A = _A;
		B = _B;

		double alpha = rtnewt( -2.0, 2.0, 1e-2, originalCurveLength );
			
		normalizedTranslator T = normalize( _translator( trs * Translator( alpha ) * inverse(trs)) );

		Bprime = T * B * inverse(T);

		return T;
	}

	void SetPositions(const normalizedPoint& _A, const normalizedPoint& _B)
	{
		A = _A;
		B = _B;
		Bprime = B;
	}

	normalizedTranslator GetCenterTranslator()
	{
		return _normalizedTranslator( exp( _freeVector(-0.5* _vectorE3GA(center)^ni) ) );
	}

	void DrawArc()
	{
		dualCircle geodesicAB = _dualCircle(dual(ToDisk(A) ^ ToDisk(Bprime) ^ poincareSphere));

		::DrawArc(_dualCircle(trs * geodesicAB * inverse(trs)), A, Bprime);
	}

	void DrawPoincareSphere()
	{
		float		objectDiffuse[] = { .5f, .5f, .0f, 0.3f };
		DrawTransparentDualSphere( normalize( _dualSphere(trs * poincareSphere * inverse(trs)) ), objectDiffuse );
	}

	void MoveSphere(const normalizedTranslator &T)
	{
		vectorE3GA a = _vectorE3GA(oA);
		vectorE3GA b = _vectorE3GA(oB);

		bivectorE3GA p = _bivectorE3GA(normalize(b - a) << reverse(I3));
		normalizedPoint centerPtf = _normalizedPoint(T * c3gaPoint(center) * reverse(T));
		vectorE3GA motion = _vectorE3GA(centerPtf) - center;
		vectorE3GA v = _vectorE3GA(( motion << p ) * inverse(p));

		center = center + v;
		radius = _double(norm_e(center - a));
		trs = _TRSversor(exp(_freeVector(-0.5 * (center * ni)) ) * exp(log(radius)*0.5*no*ni));
		originalCurveLength = G(0);
	}

	rotor GetInitialRotation()
	{
		vectorE3GA a = _vectorE3GA(ToDisk(oA));
		vectorE3GA b = _vectorE3GA(ToDisk(oB));
		return  _rotor((1.0 - (b*a)) * (1.0 / sqrt( 2.0*(1.0 - _double(b << a) ) )));
	}

	normalizedPoint ToDisk(const normalizedPoint& p)
	{
		return normalize( _point( inverse(trs) * p * trs ) );
	}

	normalizedPoint ToPlane(const point& p)
	{
		return normalize( _point(trs * p * inverse(trs) ) );
	}

	normalizedTranslator Translator( double alpha )
	{
		vectorE3GA v = normalize(_vectorE3GA(B) - _vectorE3GA(A));
		return _normalizedTranslator( 1.0 - 0.5 * alpha * v * ni );
	}

	circle Circle(double alpha)
	{
		normalizedTranslator T = Translator(alpha);

		return normalize(_circle(ToDisk(A) ^ (T * ToDisk(B) * inverse(T)) ^ poincareSphere));
	}

	normalizedPoint CO(double alpha)
	{
		circle c = Circle(alpha);
		dualCircle dc = _dualCircle( dual(c) );
		return _normalizedPoint(-0.5*(dc * ni * dc) * inverse(norm_r2(ni << dc)));
	}

	double Theta(double alpha)
	{
		normalizedTranslator T = Translator(alpha);
		vectorE3GA a = _vectorE3GA( ToDisk(A) );
		vectorE3GA co = _vectorE3GA( CO(alpha) );
		vectorE3GA b = _vectorE3GA( normalize( _point( T * ToDisk(B) * inverse(T) ) ) );

		double _a[3] = { a.e1(), a.e2(), a.e3() };
		double _b[3] = { b.e1(), b.e2(), b.e3() };
		double _co[3] = { co.e1(), co.e2(), co.e3() };
		return angle_rad_3d( _a, _co, _b );
		//return arc_cosine( _double(normalize(a - co) << normalize(b - co)) );
	}

	double Radius( double alpha )
	{
		circle c = Circle(alpha);
		dualCircle dc = _dualCircle( dual(c) );
		return sqrt(abs(_double(-1*_scalar(dc * gradeInvolution(dc)) * inverse(norm_r2(ni << dc)))));
	}

	double G(double alpha)
	{
		return Theta(alpha)*Radius(alpha);
	}

	void F( const double alpha, const double L, double &f, double &df)
	{
		const double epsilon = 1e-5;
		double fe;

		f = abs(L - G(alpha));
		fe = abs(L - G(alpha + epsilon));
		df = (fe - f) / epsilon;
	}

	//newton raphson
	double rtnewt(const double x1, const double x2, const double xacc, const double L)
	{
		const int JMAX = 20;
		int j;
		double df, dx, f, rtn;

		rtn = 0.5 * ( x1 + x2 );
		for ( j = 0 ; j < JMAX ; j++ )
		{
			F(rtn, L, f, df);

			dx = f / df;

			rtn -= dx;

			//if ( (x1 - rtn ) * ( rtn - x2 ) < 0.0)
			//	throw std::exception("Jumped out of brackets in rtnewt");

			if (fabs(dx) < xacc)
				return rtn;
		}
		//throw std::exception("Maximum number of iterations exceeded in rtnewt");
		return 0.0;
	}
};

#endif