///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridUtilities.h
///	\author  Paul Ullrich
///	\version March 8, 2024
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

#ifndef _GRIDUTILITIES_H_
#define _GRIDUTILITIES_H_

#include "GridFace.h"
#include "GridNode.h"
#include "TriangularQuadrature.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the maximum chord length of the edges of a Face.
///	</summary>
double MaxChordLength(
	const Face & face,
	const NodeVector & nodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the area of a face using a single triangular quadrature rule.
///	</summary>
double CalculateFaceAreaTriangularQuadrature(
	const Face & face,
	const NodeVector & nodes,
	const TriangularQuadratureRule & triquadrule
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the area of a face using a single triangular quadrature rule,
///		splitting the face into multiple faces when it's large.
///	</summary>
double CalculateFaceAreaTriangularQuadratureSplit(
	const FaceVector & faces,
	const NodeVector & nodes,
	const TriangularQuadratureRule & triquadrule
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the area of a face using two triangular quadrature rules
///		(chosen depending on edge length) and splitting faces that are
///		too large.
///	</summary>
double CalculateFaceAreaTriangularQuadratureMultiOrder(
	const Face & face,
	const NodeVector & nodes,
	const TriangularQuadratureRule & triquadrule_loworder,
	const TriangularQuadratureRule & triquadrule_highorder
);

///////////////////////////////////////////////////////////////////////////////

#endif // _GRIDUTILITIES_H_
