#if !defined(__ICP__H__)
#define __ICP__H__

#include "MatrixNxM.h"
#include "GA/c3ga.h"
#include "GA/c3ga_util.h"

class ICP
{
public:
	c3ga::rotor CalculateOptimalTransformation(MatrixNxM &covariance);
};


#endif