///////////////////////////////////////////////////////////////////////////////
///
///	\file    GaussQuadrature.h
///	\author  Paul Ullrich
///	\version February 22, 2024
///
///	<remarks>
///		Copyright 2024 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the BSD-3
///		License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _GAUSSQUADRATURE_H_
#define _GAUSSQUADRATURE_H_

#include <vector>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Quadrature nodes and weights for Gaussian quadrature.
///	</summary>
class GaussQuadrature {

public:
	///	<summary>
	///		Return the Gaussian quadrature points and their corresponding
	///		weights for the given number of points.
	///	</summary>
	static void GetPoints(
		int nCount,
		std::vector<double> & dG,
		std::vector<double> & dW
	);

	///	<summary>
	///		Return the Gaussian quadrature points and their corresponding
	///		weights for the given number of points and 1D reference element.
	///	</summary>
	static void GetPoints(
		int nCount,
		double dXi0,
		double dXi1,
		std::vector<double> & dG,
		std::vector<double> & dW
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif
