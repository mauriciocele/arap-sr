#if !defined(__MATRIXNXM__H__)
#define __MATRIXNXM__H__

#include "svd.h"
#include <geometry.H>
#include <string>
#include <functional>
#include "GA/c3ga.h"
#include "GA/c3ga_util.h"
#include "GA/gl_util.h"

class MatrixNxM
{
private:
	void FreeMemory()
	{
		if (data != nullptr)
		{
			delete[] data;
			data = nullptr;
		}
		if (buffer != nullptr)
		{
			delete[] buffer;
			buffer = nullptr;
		}
	}

	void AllocateMemory()
	{
		data = new double*[rows];
		buffer = new double[rows*cols];
		for (int i = 0; i < rows; ++i)
		{
			data[i] = buffer + i * cols;
		}
	}

public:
	double **data;
	double *buffer;
	int rows;
	int cols;

	const double*	operator [] (int i) const { return data[i]; }
	double*			operator [] (int i) { return data[i]; }

	operator double* () { return buffer; }
	operator const double* () const { return buffer; }

	MatrixNxM()
	{
		this->rows = 0;
		this->cols = 0;
		this->data = nullptr;
		this->buffer = nullptr;
	}

	MatrixNxM(int rows, int cols)
	{
		this->rows = rows;
		this->cols = cols;
		AllocateMemory();
	}

	// Copy constructor.
	MatrixNxM(const MatrixNxM& matrixNxM)
	{
		this->rows = matrixNxM.rows;
		this->cols = matrixNxM.cols;
		AllocateMemory();

		memcpy(this->buffer, matrixNxM.buffer, sizeof(double)*rows*cols);
	}

	// Move constructor.
	MatrixNxM(MatrixNxM&& matrixNxM) //move cons
	{
		this->rows = matrixNxM.rows;
		this->cols = matrixNxM.cols;
		this->data = matrixNxM.data;
		this->buffer = matrixNxM.buffer;
		matrixNxM.data = nullptr;
		matrixNxM.buffer = nullptr;
	}

	// Copy Assignment operator.
	MatrixNxM& operator = (const MatrixNxM& matrixNxM)
	{
		if (this->rows != matrixNxM.rows || this->cols != matrixNxM.cols)
		{
			FreeMemory();
			this->rows = matrixNxM.rows;
			this->cols = matrixNxM.cols;
			AllocateMemory();
		}

		memcpy(this->buffer, matrixNxM.buffer, sizeof(double)*rows*cols);

		return *this;
	}

	// Move assignment operator.
	MatrixNxM& operator=(MatrixNxM&& matrixNxM)
	{
		if (this != &matrixNxM)
		{
			FreeMemory();
			this->rows = matrixNxM.rows;
			this->cols = matrixNxM.cols;
			this->data = matrixNxM.data;
			this->buffer = matrixNxM.buffer;
			matrixNxM.data = nullptr;
			matrixNxM.buffer = nullptr;
		}
		return *this;
	}

	~MatrixNxM()
	{
		FreeMemory();
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
			if (i < cols)
			{
				data[i][i] = 1.0;
			}
		}
	}

	void FromArray(double *matrix)
	{
		memcpy(this->buffer, matrix, sizeof(double)*rows*cols);
	}

	MatrixNxM CopyMatrix() const
	{
		MatrixNxM c(rows, cols);

		memcpy(c.buffer, this->buffer, sizeof(double)*rows*cols);

		return c;
	}

	MatrixNxM Transpose() const
	{
		MatrixNxM c(cols, rows);

		for (int i = 0; i < rows; ++i)
		for (int j = 0; j < cols; ++j)
			c.data[j][i] = data[i][j];
		return c;
	}

	void operator *= (double b)
	{
		for (int i = 0; i < rows*cols; ++i)
		{
			buffer[i] *= b;
		}
	}

	void operator += (const MatrixNxM& b)
	{
		for (int i = 0; i < rows*cols; ++i)
		{
			buffer[i] += b.buffer[i];
		}
	}

	void operator -= (const MatrixNxM& b)
	{
		for (int i = 0; i < rows*cols; ++i)
		{
			buffer[i] -= b.buffer[i];
		}
	}

	double Determinant() const
	{
		if (cols == 3 && rows == 3)
		{
			return(data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) -
				data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0]) +
				data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]));
		}
		else if (cols == 2 && rows == 2)
		{
			return data[0][0] * data[1][1] - data[0][1] * data[1][0];
		}
		int N = this->rows;

		double *a = new double[N * N];
		int *pivot = new int[N];

		auto IND = [N](int i, int j) { return i + j * N; };

		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < N; ++j)
			{
				a[IND(i, j)] = data[i][j];
			}
		}

		int info = dge_fa(N, a, pivot);

		if (info != 0)
			throw std::string("The factorization failed");

		double det = dge_det(N, a, pivot);

		delete[] a;
		delete[] pivot;

		return det;
	}

	void SVD(MatrixNxM &U, MatrixNxM &W, MatrixNxM &V) const
	{
		MatrixNxM u(rows + 1, cols + 1);
		MatrixNxM v(cols + 1, cols + 1);

		double *w = new double[cols + 1];

		for (int i = 1; i <= this->rows; ++i)
		for (int j = 1; j <= this->cols; ++j)
			u.data[i][j] = this->data[i - 1][j - 1];

		svdcmp(u.data, rows, cols, w, v.data);

		if (U.rows != rows || U.cols != cols)
			U = MatrixNxM(rows, cols);
		if (V.rows != cols || V.cols != cols)
			V = MatrixNxM(cols, cols);
		if (W.rows != cols || W.cols != cols)
			W = MatrixNxM(cols, cols);

		for (int i = 1; i <= this->rows; ++i)
		for (int j = 1; j <= this->cols; ++j)
			U.data[i - 1][j - 1] = u.data[i][j];

		for (int i = 1; i <= this->cols; ++i)
		for (int j = 1; j <= this->cols; ++j)
			V.data[i - 1][j - 1] = v.data[i][j];

		W.ZeroMatrix();

		for (int i = 1; i <= this->cols; ++i)
			W.data[i - 1][i - 1] = w[i];

		delete[] w;
	}

	void SolveSVD(const MatrixNxM& rhs, MatrixNxM* res) const
	{
		MatrixNxM u(rows + 1, cols + 1);
		MatrixNxM v(cols + 1, cols + 1);
		double *w = new double[cols + 1];
		double *b = new double[rows + 1];
		double *x = new double[cols + 1];

		const double illConditionedThreshold = 1e-6;

		for (int i = 1; i <= this->rows; ++i)
		for (int j = 1; j <= this->cols; ++j)
			u.data[i][j] = this->data[i - 1][j - 1];

		svdcmp(u.data, rows, cols, w, v.data);

		for (int i = 1; i <= cols; ++i)
		if (w[i] < illConditionedThreshold)
			w[i] = 0.0;

		for (int j = 0; j < rhs.cols; ++j)
		{
			for (int i = 0; i < rhs.rows; ++i)
			{
				b[i + 1] = rhs.data[i][j];
			}

			svbksb(u.data, w, v.data, rows, cols, b, x);

			for (int i = 0; i < res->rows; ++i)
			{
				res->data[i][j] = x[i + 1];
			}
		}

		delete[] w;
		delete[] b;
		delete[] x;
	}

	MatrixNxM SolveSVD(const MatrixNxM& rhs) const
	{
		MatrixNxM res(rhs.rows, rhs.cols);
		SolveSVD(rhs, &res);
		return res;
	}

	void Solve(const MatrixNxM& rhs, MatrixNxM* res) const
	{
		int N = this->rows;
		int RHS_NUM = rhs.cols;

		double *a = new double[N * (N + RHS_NUM)];

		auto IND = [N](int i, int j) { return i + j * N; };

		for (int i = 0; i < this->rows; ++i)
		{
			for (int j = 0; j < this->cols; ++j)
			{
				a[IND(i, j)] = this->data[i][j];
			}
			for (int j = 0; j < rhs.cols; ++j)
			{
				a[IND(i, this->cols + j)] = rhs.data[i][j];
			}
		}

		int solution = r8mat_solve(N, RHS_NUM, a);

		if (solution != 0)
			throw std::string("factorization failed. The solutions could not be computed.");

		for (int i = 0; i < res->rows; ++i)
		{
			for (int j = 0; j < rhs.cols; ++j)
			{
				res->data[i][j] = a[IND(i, this->cols + j)];
			}
		}

		delete[] a;
	}

	MatrixNxM Solve(const MatrixNxM& rhs) const
	{
		MatrixNxM res(rhs.rows, rhs.cols);
		Solve(rhs, &res);
		return res;
	}

	void SolveLU(const MatrixNxM& rhs, MatrixNxM* res) const
	{
		const int N = this->rows;

		double *a = new double[N * N];
		double *b = new double[N];
		int	*pivot = new int[N];

		auto IND = [N](int i, int j) { return i + j * N; };

		for (int i = 0; i < this->rows; ++i)
		{
			for (int j = 0; j < this->cols; ++j)
			{
				a[IND(i, j)] = this->data[i][j];
			}
		}

		int info = dge_fa(N, a, pivot);

		if (info != 0)
			throw std::string("The factorization failed");

		for (int j = 0; j < rhs.cols; ++j)
		{
			for (int i = 0; i < rhs.rows; ++i)
			{
				b[i] = rhs.data[i][j];
			}

			dge_sl(N, a, pivot, b, 0);

			for (int i = 0; i < res->rows; ++i)
			{
				res->data[i][j] = b[i];
			}
		}

		delete[] a;
		delete[] b;
		delete[] pivot;
	}

	MatrixNxM SolveLU(const MatrixNxM& rhs) const
	{
		MatrixNxM res(rhs.rows, rhs.cols);
		SolveLU(rhs, &res);
		return res;
	}

	c3ga::rotor ToRotor() const
	{
		double trace = data[0][0] + data[1][1] + data[2][2] + 1.0;
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
				double s = 2.0 * (double)sqrt(1.0 + data[0][0] - data[1][1] - data[2][2]);
				qx = 0.25 * s;
				qy = (data[0][1] + data[1][0]) / s;
				qz = (data[0][2] + data[2][0]) / s;
				qw = (data[1][2] - data[2][1]) / s;
			}
			else if (data[1][1] > data[2][2]) {
				double s = 2.0 * (double)sqrt(1.0 + data[1][1] - data[0][0] - data[2][2]);
				qx = (data[0][1] + data[1][0]) / s;
				qy = 0.25 * s;
				qz = (data[1][2] + data[2][1]) / s;
				qw = (data[0][2] - data[2][0]) / s;
			}
			else {
				double s = 2.0 * (double)sqrt(1.0 + data[2][2] - data[0][0] - data[1][1]);
				qx = (data[0][2] + data[2][0]) / s;
				qy = (data[1][2] + data[2][1]) / s;
				qz = 0.25 * s;
				qw = (data[0][1] - data[1][0]) / s;
			}
		}

		double s = (double)sqrt(qw *qw + qx * qx + qy * qy + qz * qz);

		return c3ga::rotor(c3ga::rotor_scalar_e1e2_e2e3_e3e1, qw / s, -qz / s, -qx / s, -qy / s);
	}
};

inline void TransposedMultAccum(const MatrixNxM& a, const MatrixNxM& b, MatrixNxM* c)
{
	for (int i = 0; i < a.cols; ++i)
	for (int j = 0; j < b.cols; ++j)
	for (int k = 0; k < a.rows; ++k)
		c->data[i][j] += a.data[k][i] * b.data[k][j];
}

inline void TransposedMult(const MatrixNxM& a, const MatrixNxM& b, MatrixNxM* c)
{
	c->ZeroMatrix();
	TransposedMultAccum(a, b, c);
}

inline MatrixNxM TransposedMult(const MatrixNxM& a, const MatrixNxM& b)
{
	MatrixNxM c(a.cols, b.cols);
	TransposedMult(a, b, &c);
	return c;
}

inline void MultAccum(const MatrixNxM& a, const MatrixNxM& b, MatrixNxM* c)
{
	for (int i = 0; i < a.rows; ++i)
	for (int j = 0; j < b.cols; ++j)
	for (int k = 0; k < a.cols; ++k)
		c->data[i][j] += a.data[i][k] * b.data[k][j];
}

inline void Mult(const MatrixNxM& a, const MatrixNxM& b, MatrixNxM* c)
{
	c->ZeroMatrix();
	MultAccum(a, b, c);
}

inline void Mult(const MatrixNxM& a, const double b, MatrixNxM* c)
{
	for (int i = 0; i < a.rows; ++i)
	for (int j = 0; j < a.cols; ++j)
		c->data[i][j] = a.data[i][j] * b;
}

inline void Add(const MatrixNxM& a, const MatrixNxM& b, MatrixNxM* c)
{
	for (int i = 0; i < a.rows; ++i)
	for (int j = 0; j < a.cols; ++j)
		c->data[i][j] = a.data[i][j] + b.data[i][j];
}

inline void Sub(const MatrixNxM& a, const MatrixNxM& b, MatrixNxM* c)
{
	for (int i = 0; i < a.rows; ++i)
	for (int j = 0; j < a.cols; ++j)
		c->data[i][j] = a.data[i][j] - b.data[i][j];
}

inline MatrixNxM operator *(const MatrixNxM& a, const MatrixNxM& b)
{
	MatrixNxM c(a.rows, b.cols);
	Mult(a, b, &c);
	return c;
}

inline MatrixNxM operator *(const MatrixNxM& a, double b)
{
	MatrixNxM c(a.rows, a.cols);
	Mult(a, b, &c);
	return c;
}

inline MatrixNxM operator *(double b, const MatrixNxM& a)
{
	return a * b;
}

inline MatrixNxM operator +(const MatrixNxM& a, const MatrixNxM& b)
{
	MatrixNxM c(a.rows, a.cols);
	Add(a, b, &c);
	return c;
}

inline MatrixNxM operator -(const MatrixNxM& a, const MatrixNxM& b)
{
	MatrixNxM c(a.rows, a.cols);
	Sub(a, b, &c);
	return c;
}

inline MatrixNxM operator /(const MatrixNxM& a, const MatrixNxM& b)
{
	return a.SolveLU(b);
}

inline void Solve(const MatrixNxM& a, const MatrixNxM& rhs, MatrixNxM* res)
{
	a.Solve(rhs, res);
}

inline void SolveLU(const MatrixNxM& a, const MatrixNxM& rhs, MatrixNxM* res)
{
	a.SolveLU(rhs, res);
}

inline void SolveSVD(const MatrixNxM& a, const MatrixNxM& rhs, MatrixNxM* res)
{
	a.SolveSVD(rhs, res);
}

#endif
