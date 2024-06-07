#pragma once
#ifndef _Cs
#define _Cs ::custom::
#endif // !_Cs

#include "Identifier.h"

namespace custom {
	
	void print(const Eigen::MatrixXd& mat, int max_row = -1);

	void print(const Eigen::ArrayXd& arr);

}