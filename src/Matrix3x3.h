#if !defined(__MATRIX3X3__H__)
#define __MATRIX3X3__H__

#include "svd.h"
#include <geometry.H>
#include <string>
#include <functional>
#include "GA/c3ga.h"
#include "GA/c3ga_util.h"
#include "GA/gl_util.h"

class Matrix3x3
{
public:
    double data[3][3];

	inline const double*	operator [] ( int i ) const { return data[ i ]; }
	inline double*			operator [] ( int i ) { return data[ i ]; }

	inline operator double* () { return (double*)data; }
	inline operator const double* () const { return (const double*)data; }

    Matrix3x3()
	{
	}

    // Copy constructor.
	Matrix3x3(const Matrix3x3& matrix3x3)
    {
		memcpy(data, matrix3x3.data, sizeof(double)*9); 
    }

	// Copy Assignment operator.
	Matrix3x3 operator = (const Matrix3x3& matrix3x3)
	{		
		memcpy(data, matrix3x3.data, sizeof(double)*9); 

		return *this;
	}

    void ZeroMatrix()
    {
		memset(data, 0, sizeof(double)*9); 
    }

    void Identity()
    {
		memset(data, 0, sizeof(double)*9); 
        for (int i = 0; i < 3; ++i)
		{
			data[i][i] = 1.0;
		}
    }

    Matrix3x3 CopyMatrix()
    {
        Matrix3x3 c;
		memcpy(c.data, data, sizeof(double)*9); 
        return c;
    }

	Matrix3x3 Transpose()
	{
		Matrix3x3 c;
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				c.data[j][i] = data[i][j];
		return c;
	}

    double Determinant()
    {
		return( data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) -
				data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0]) +
				data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]) );
    }

	void SVD(Matrix3x3 &U, Matrix3x3 &W, Matrix3x3 &V)
	{
		double u[4][4], v[4][4];
		double w[4];

		for (int i = 1; i <= 3; ++i)
			for (int j = 1; j <= 3; ++j)
				u[i][j] = this->data[i-1][j-1];

		svdcmp3x3(u, 3, 3, w, v);

		for (int i = 1; i <= 3; ++i)
			for (int j = 1; j <= 3; ++j)
				U.data[i-1][j-1] = u[i][j];

		for (int i = 1; i <= 3; ++i)
			for (int j = 1; j <= 3; ++j)
				V.data[i - 1][j - 1] = v[i][j];

		W.ZeroMatrix();

		for (int i = 1; i <= 3; ++i)
			W.data[i - 1][i - 1] = w[i];
	}

	c3ga::rotor ToRotor()
	{
		double trace = data[0][0] + data[1][1] + data[2][2] + 1.0f;
		double qw; // scalar coordinate
		double qx; // coordinate for -e2^e3
		double qy; // coordinate for -e3^e1
		double qz; // coordinate for -e1^e2
		if (trace > 0.00001) {
			double s = 0.5 / (double)sqrt(trace);
			qw = 0.25 / s;
			qw = sqrt(trace) * (0.5);
			qx = (data[2][1] - data[1][2]) * s;
			qy = (data[0][2] - data[2][0]) * s;
			qz = (data[1][0] - data[0][1]) * s;
		}
		else {
			if (data[0][0] > data[1][1] && data[0][0] > data[2][2]) {
				double s = 2.0 * (double)sqrt( 1.0 + data[0][0] - data[1][1] - data[2][2]);
				qx = 0.25 * s;
				qy = (data[0][1] + data[1][0]) / s;
				qz = (data[0][2] + data[2][0]) / s;
				qw = (data[1][2] - data[2][1]) / s;
			}
			else if (data[1][1] > data[2][2]) {
				double s = 2.0 * (double)sqrt( 1.0 + data[1][1] - data[0][0] - data[2][2]);
				qx = (data[0][1] + data[1][0]) / s;
				qy = 0.25 * s;
				qz = (data[1][2] + data[2][1]) / s;
				qw = (data[0][2] - data[2][0]) / s;
			}
			else {
				double s = 2.0 * (double)sqrt( 1.0 + data[2][2] - data[0][0] - data[1][1] );
				qx = (data[0][2] + data[2][0]) / s;
				qy = (data[1][2] + data[2][1]) / s;
				qz = 0.25 * s;
				qw = (data[0][1] - data[1][0]) / s;
			}
		}

		double s = (double) sqrt(qw *qw + qx * qx + qy * qy + qz * qz);

		return c3ga::rotor(c3ga::rotor_scalar_e1e2_e2e3_e3e1, qw / s, -qz / s, -qx / s, -qy / s);
	}

};

inline Matrix3x3 operator *(const Matrix3x3& a, const Matrix3x3& b) 
{
    Matrix3x3 c;

	c.ZeroMatrix();

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
			for( int k = 0 ; k < 3 ; ++k )
				c.data[i][j] += a.data[i][k] * b.data[k][j];

	return c;
}

inline vectorE3GA operator *(const Matrix3x3& a, const vectorE3GA& b) 
{
    vectorE3GA c;

	c.m_c[0] = c.m_c[1] = c.m_c[2] = 0.0;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
			c.m_c[i] += a.data[i][j] * b.m_c[j];

	return c;
}

inline Matrix3x3 operator *(const Matrix3x3& a, double b) 
{
    Matrix3x3 c;

    for ( int i = 0 ; i < 3 ; ++i)
		for ( int j = 0 ; j < 3 ; ++j )
			c.data[i][j] = a.data[i][j] * b;

	return c;
}

inline Matrix3x3 operator *(double b, const Matrix3x3& a) 
{
    Matrix3x3 c;

    for ( int i = 0 ; i < 3 ; ++i)
		for ( int j = 0 ; j < 3 ; ++j )
			c.data[i][j] = a.data[i][j] * b;

	return c;
}

inline Matrix3x3 operator +(const Matrix3x3& a, const Matrix3x3& b) 
{
    Matrix3x3 c;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            c.data[i][j] = a.data[i][j] + b.data[i][j];

    return c;
}


#endif
