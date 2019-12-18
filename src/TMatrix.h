#if !defined(__TMATRIX__H__)
#define __TMATRIX__H__

#include "svd.h"
#include <geometry.H>
#include <string>
#include <functional>
#include "GA/c3ga.h"
#include "GA/c3ga_util.h"
#include "GA/gl_util.h"

template <const int rows, const int cols>
class TMatrix
{
private:
	void MapMemory()
	{
		for(int i = 0 ; i < rows ; ++i)
		{
			data[i] = buffer + i * cols;
		}
	}

public:
    double *data[rows];
    double buffer[rows*cols];

	inline const double*	operator [] ( int i ) const { return data[ i ]; }
	inline double*			operator [] (int i) { return data[i]; }

	inline operator double* () { return buffer; }
	inline operator const double* () const { return buffer; }

    TMatrix()
	{
		MapMemory();
	}

    ~TMatrix()
	{
	}

    // Copy constructor.
	TMatrix(const TMatrix<rows,cols>& matrix)
    {
		MapMemory();
		memcpy(this->buffer, matrix.buffer, sizeof(double)*rows*cols);
    }

	// Copy Assignment operator.
	void operator = (const TMatrix<rows,cols>& matrix)
	{		
		memcpy(this->buffer, matrix.buffer, sizeof(double)*rows*cols);
	}

    void ZeroMatrix()
    {
		memset(this->buffer, 0, sizeof(double)*rows*cols);
    }

    void Identity()
    {
		memset(this->buffer, 0, sizeof(double)*rows*cols);
        for (int i = 0; i < rows; ++i)
		{
			if(i < cols)
			{
				data[i][i] = 1.0;
			}
		}
    }

	void FromArray(double *matrix)
	{
		memcpy(this->buffer, matrix, sizeof(double)*rows*cols);
	}

	TMatrix<cols,rows> Transpose() const
	{
		TMatrix<cols,rows> c;
		Transpose(&c);
		return c;
	}

	void Transpose(TMatrix<cols,rows>* c) const
	{
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				c->data[j][i] = data[i][j];
	}

    double Determinant() const
    {
		if( cols == 3 && rows == 3 )
		{
			return( data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) -
					data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0]) +
					data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]) );
		}
		else if( cols == 2 && rows == 2 )
		{
			return data[0][0] * data[1][1] - data[0][1] * data[1][0];
		}

        double a[rows * rows];
        int pivot[rows];

        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < rows; ++j)
            {
                a[i + j * rows] = data[i][j];
            }
        }

        int info = dge_fa(rows, a, pivot);

        if (info != 0)
            throw std::string("The factorization failed");

        double det = dge_det(rows, a, pivot);

		return det;
    }

	template <const int cols2>
	void TransposedMultAccum(const TMatrix<rows,cols2>& b, TMatrix<cols,cols2>* c) const
	{
		for (int i = 0; i < cols; ++i)
			for (int k = 0; k < rows; ++k)
				for (int j = 0; j < cols2; ++j)
					c->data[i][j] += data[k][i] * b.data[k][j];
	}

	template <const int cols2>
	void TransposedMult(const TMatrix<rows,cols2>& b, TMatrix<cols,cols2>* c) const
	{
		c->ZeroMatrix();
		for (int i = 0; i < cols; ++i)
			for( int k = 0 ; k < rows ; ++k )
				for (int j = 0; j < cols2; ++j)
					c->data[i][j] += data[k][i] * b.data[k][j];
	}

	template <const int cols2>
	void Mult(const TMatrix<cols,cols2>& b, TMatrix<rows,cols2>* c) const
	{
		c->ZeroMatrix();
		for (auto i = 0; i != rows; ++i) {
			const auto row_A = i * cols;
			const auto row_C = i * cols2;
			for (auto k = 0; k != cols; ++k) {
				const auto row_B = k * cols2;
				const auto r = buffer[row_A + k];
				for (auto j = 0; j != cols2; ++j) {
					c->buffer[row_C + j] = c->buffer[row_C + j] + r * b.buffer[row_B + j];
				}
			}
		}	
	
		//c->ZeroMatrix();
		//for (int i = 0; i < rows; ++i)
		//	for (int k = 0; k < cols; ++k)
		//		for (int j = 0; j < cols2; ++j)
		//			c->data[i][j] += data[i][k] * b.data[k][j];
	}

	void Mult(double b, TMatrix<rows,cols>* c) const
	{
		for ( int i = 0 ; i < rows*cols ; ++i)
		{
			c->buffer[i] = buffer[i] * b;
		}
	}

	void Add(const TMatrix<rows,cols>& b, TMatrix<rows,cols>* c) const
	{
		for (int i = 0; i < rows*cols; ++i)
		{
			c->buffer[i] = buffer[i] + b.buffer[i];
		}
	}

	void Sub(const TMatrix<rows,cols>& b, TMatrix<rows,cols>* c) const
	{
		for (int i = 0; i < rows*cols; ++i)
		{
			c->buffer[i] = buffer[i] - b.buffer[i];
		}
	}

	void operator *= (double b) 
	{
		for ( int i = 0 ; i < rows*cols ; ++i)
		{
			buffer[i] *= b;
		}
	}

	void operator += (const TMatrix<rows,cols>& b) 
	{
		for (int i = 0; i < rows*cols; ++i)
		{
			buffer[i] += b.buffer[i];
		}
	}

	void operator -= (const TMatrix<rows,cols>& b) 
	{
		for (int i = 0; i < rows*cols; ++i)
		{
			buffer[i] -= b.buffer[i];
		}
	}

	template <const int cols2>
    void Solve(const TMatrix<rows,cols2>& rhs, TMatrix<cols,cols2>* res) const
    {
        const int N = rows;
        const int RHS_NUM = cols2;

        double a[N * (N + RHS_NUM)];

        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                a[i + j * N] = data[i][j];
            }
            for (int j = 0; j < cols2; ++j)
            {
				a[i + (cols + j) * N] = rhs.data[i][j];
            }
        }

        int solution = r8mat_solve(N, RHS_NUM, a);

        if (solution != 0)
            throw std::string("factorization failed. The solutions could not be computed.");

        for (int i = 0; i < cols; ++i)
        {
            for (int j = 0; j < cols2; ++j)
            {
                res->data[i][j] = a[i + (cols + j) * N];
            }
        }
    }

	template <const int cols2>
    void SolveLU(const TMatrix<rows,cols2>& rhs, TMatrix<cols,cols2>* res) const
    {
        const int N = rows;

        double a[N * N];
        double b[N];
        int	pivot[N];

        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
				a[i + j * N] = data[i][j];
            }
        }

        int info = dge_fa(N, a, pivot);

        if (info != 0)
            throw std::string("The factorization failed");

        for (int j = 0; j < cols2; ++j)
        {
            for (int i = 0; i < rows; ++i)
            {
				b[i] = rhs.data[i][j];
            }

            dge_sl(N, a, pivot, b, 0);

            for (int i = 0; i < rows; ++i)
            {
				res->data[i][j] = b[i];
            }
        }
    }

	template <const int cols2>
    void SolveSVD(const TMatrix<rows,cols2>& rhs, TMatrix<cols,cols2>* res) const
	{
		double w[cols + 1];
		double b[rows + 1];
		double x[cols + 1];
		double *u[rows+1], *v[cols+1];
		double mu[(rows+1)*(cols+1)];
		double mv[(cols+1)*(cols+1)];

		for (int i = 0; i <= rows; ++i)
		{
			u[i] = mu + i*(cols+1);
		}
		for (int i = 0; i <= cols; ++i)
		{
			v[i] = mv + i*(cols+1);
		}

		const double illConditionedThreshold = 1e-6;

		for (int i = 1; i <= rows; ++i)
			for (int j = 1; j <= cols; ++j)
				u[i][j] = data[i - 1][j - 1];

		svdcmp(u, rows, cols, w, v);

		for (int i = 1; i <= cols; ++i)
			if (w[i] < illConditionedThreshold)
				w[i] = 0.0;

		for (int j = 0; j < cols2; ++j)
		{
			for (int i = 0; i < rows; ++i)
			{
				b[i + 1] = rhs.data[i][j];
			}

			svbksb(u, w, v, rows, cols, b, x);

			for (int i = 0; i < cols; ++i)
			{
				res->data[i][j] = x[i + 1];
			}
		}
	}
};

template <const int rows, const int cols, const int cols2>
TMatrix<rows,cols2> operator *(const TMatrix<rows,cols>& a, const TMatrix<cols,cols2>& b)
{
	TMatrix<rows,cols2> c;
	a.Mult(b, &c);
	return c;
}

template <const int rows, const int cols, const int cols2>
TMatrix<cols,cols2> operator /(const TMatrix<rows,cols>& a, const TMatrix<rows,cols2>& b)
{
	TMatrix<cols,cols2> c;
	a.Solve(b, &c);
	return c;
}

template <const int rows, const int cols>
TMatrix<rows,cols> operator *(const TMatrix<rows,cols>& a, double b) 
{
	TMatrix<rows,cols> c;
	a.Mult(b, &c);
	return c;
}

template <const int rows, const int cols>
TMatrix<rows,cols> operator *(double b, const TMatrix<rows,cols>& a) 
{
	TMatrix<rows,cols> c;
	a.Mult(b, &c);
	return c;
}

template <const int rows, const int cols>
TMatrix<rows,cols> operator +(const TMatrix<rows,cols>& a, const TMatrix<rows,cols>& b) 
{
	TMatrix<rows,cols> c;
	a.Add(b, &c);
	return c;
}

template <const int rows, const int cols>
TMatrix<rows,cols> operator -(const TMatrix<rows,cols>& a, const TMatrix<rows,cols>& b) 
{
	TMatrix<rows,cols> c;
	a.Sub(b, &c);
	return c;
}

template <const int rows, const int cols, const int cols2>
void TransposedMult(const TMatrix<rows,cols>& a, const TMatrix<rows,cols2>& b, TMatrix<cols,cols2>* c) 
{
	a.TransposedMult(b, c);
}

template <const int rows, const int cols, const int cols2>
void TransposedMultAccum(const TMatrix<rows,cols>& a, const TMatrix<rows,cols2>& b, TMatrix<cols,cols2>* c) 
{
	a.TransposedMultAccum(b, c);
}

template <const int rows, const int cols, const int cols2>
void Mult(const TMatrix<rows,cols>& a, const TMatrix<cols,cols2>& b, TMatrix<rows,cols2>* c) 
{
	a.Mult(b, c);
}

template <const int rows, const int cols>
void Mult(const TMatrix<rows,cols>& a, double b, TMatrix<rows,cols>* c) 
{
	a.Mult(b, c);
}

template <const int rows, const int cols>
void Add(const TMatrix<rows,cols>& a, const TMatrix<rows,cols>& b, TMatrix<rows,cols>* c) 
{
	a.Add(b, c);
}

template <const int rows, const int cols>
void Sub(const TMatrix<rows,cols>& a, const TMatrix<rows,cols>& b, TMatrix<rows,cols>* c) 
{
	a.Sub(b, c);
}

template <const int rows, const int cols, const int cols2>
void Solve(const TMatrix<rows,cols>& a, const TMatrix<rows,cols2>& rhs, TMatrix<cols,cols2>* res) 
{
	a.Solve(rhs, res);
}

template <const int rows, const int cols, const int cols2>
void SolveLU(const TMatrix<rows,cols>& a, const TMatrix<rows,cols2>& rhs, TMatrix<cols,cols2>* res) 
{
	a.SolveLU(rhs, res);
}

template <const int rows, const int cols, const int cols2>
void SolveSVD(const TMatrix<rows,cols>& a, const TMatrix<rows,cols2>& rhs, TMatrix<cols,cols2>* res) 
{
	a.SolveSVD(rhs, res);
}

#endif
