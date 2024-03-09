///////////////////////////////////////////////////////////////////////////////
///
///	\file    TriangularQuadrature.h
///	\author  Paul Ullrich
///	\version March 8, 2024
///
///	<remarks>
///		Copyright 2024 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _TRIANGULARQUADRATURE_H_
#define _TRIANGULARQUADRATURE_H_

#include <vector>
#include "array2d.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A wrapper for various triangular quadrature rules.
///	</summary>
///	<note>
///		Triangular quadrature to use for integration
///		Dunavant, D.A. "High Degree Efficient Symmetrical Gaussian Quadrature
///		Rules for the Triangle."  J. Numer. Meth. Eng., 21, pp 1129-1148.
///	</note>
class TriangularQuadratureRule {

public:
	///	<summary>
	///		Constructor with order of quadrature rule.
	///	</summary>
	TriangularQuadratureRule(
		int nOrder
	);

	///	<summary>
	///		Get the order of the rule.
	///	</summary>
	int order() const {
		return m_nOrder;
	}

	///	<summary>
	///		Get the number of points in the rule.
	///	</summary>
	size_t size() const {
		return m_dW.size();
	}

	///	<summary>
	///		Get the specified coordinate.
	///	</summary>
	double g(size_t p, size_t i) const {
		return m_dG(p,i);
	}

	///	<summary>
	///		Get the specified weight.
	///	</summary>
	double w(size_t p) const {
		return m_dW[p];
	}

	///	<summary>
	///		Validate the given triangular quadrature rule.
	///	</summary>
	void validate() const;

protected:
	///	<summary>
	///		Order of the triangular quadrature rule.
	///	</summary>
	int m_nOrder;

	///	<summary>
	///		Coordinates of triangular quadrature rule.
	///	</summary>
	tmp::array2d<double> m_dG;

	///	<summary>
	///		Weights of triangular quadrature rule.
	///	</summary>
	std::vector<double> m_dW;
};

///////////////////////////////////////////////////////////////////////////////

#endif

