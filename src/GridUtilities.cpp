///////////////////////////////////////////////////////////////////////////////
///
///	\file	GridElements.cpp
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

#include "GridUtilities.h"

///////////////////////////////////////////////////////////////////////////////

double MaxChordLength(
	const Face & face,
	const NodeVector & nodes
) {
	double dMaxChordLength = nodes[face[face.size()-1]].distance_from(nodes[face[0]]);

    for (size_t k = 0; k < face.size()-1; k++){
		double dChordLength = nodes[face[k]].distance_from(nodes[face[k+1]]);
		if (dChordLength > dMaxChordLength) {
			dMaxChordLength = dChordLength;
		}
	}
	
	return dMaxChordLength;
}

///////////////////////////////////////////////////////////////////////////////

double CalculateFaceAreaTriangularQuadrature(
	const Face & face,
	const NodeVector & nodes,
	const TriangularQuadratureRule & triquadrule
) {
	double dFaceArea = 0.0;

	if (face.size() < 3) {
		_EXCEPTIONT("CalculateFaceAreaTriangularQuadrature requires that Face have at least 3 Nodes");
	}

	const size_t sTriangles = face.size() - 2;

	for (size_t k = 0; k < sTriangles; k++) {

		const Node & node1 = nodes[face[0]];
		const Node & node2 = nodes[face[k+1]];
		const Node & node3 = nodes[face[k+2]];

		double dtriprod = node1.dot(node2.cross(node3));

		for (size_t q = 0; q < triquadrule.size(); q++) {

			Node nodeQ =
				  node1 * triquadrule.g(q,2)
				+ node2 * triquadrule.g(q,0)
				+ node3 * triquadrule.g(q,1);

			double magQ = nodeQ.magnitude();

			double magQ3 = magQ * magQ * magQ;

			dFaceArea += 0.5 * triquadrule.w(q) * dtriprod / magQ3;
		}
	}

	return dFaceArea;
}

///////////////////////////////////////////////////////////////////////////////

double CalculateFaceAreaTriangularQuadratureSplit(
	const FaceVector & faces,
	const NodeVector & nodes,
	const TriangularQuadratureRule & triquadrule
) {
	double dEdgeLengthTolerance = 0.0;
	double dTotalFaceArea = 0.0;

	if (triquadrule.order() >= 8){
		dEdgeLengthTolerance = 0.05;
	} else {
		dEdgeLengthTolerance = 0.003;
	}

    for (size_t i = 0; i < faces.size(); i++){

		if (faces[i].size() < 3) {
			_EXCEPTIONT("CalculateFaceAreaTriQuadratureSplit requires that Face have at least 3 Nodes");
		}

		const size_t nEdges = faces[i].size();
		const size_t nTriangles = nEdges - 2;

		const double dMaxChordLength = MaxChordLength(faces[i], nodes);

		if (dMaxChordLength > dEdgeLengthTolerance) {

			FaceVector facesSplit(4, Face(3));
			NodeVector nodesSplit;
			Node node;

			for (int j = 0; j < nEdges; j++){
				nodesSplit.push_back(nodes[faces[i][j]]);
			}

           	// Insert elements and points
			size_t index = nEdges;
			Node nodeSplit = (nodesSplit[0] + nodesSplit[1]) / 2.0;
			nodeSplit.normalize_in_place();
			nodesSplit.push_back(nodeSplit);

			for (size_t j = 1; j < nEdges-1; j++) {
				Face & face1 = facesSplit[0];
				Face & face2 = facesSplit[1];
				Face & face3 = facesSplit[2];
				Face & face4 = facesSplit[3];

				face1[0] = 0;       face1[1] = index;   face1[2] = index+2;
				face2[0] = index+2; face2[1] = index;   face2[2] = index+1;
				face3[0] = index+1; face3[1] = index;   face3[2] = j;
				face4[0] = index+2; face4[1] = index+1; face4[2] = j+1;

				index++;
				Node nodeSplitJ = (nodesSplit[j] + nodesSplit[j+1]) / 2.0;
				nodeSplitJ.normalize_in_place();
				nodesSplit.push_back(nodeSplitJ);

				index++;
				Node nodeSplit0 = (nodesSplit[0] + nodesSplit[j+1]) / 2.0;
				nodeSplit0.normalize_in_place();
				nodesSplit.push_back(nodeSplit0);		
			}

			dTotalFaceArea += CalculateFaceAreaTriangularQuadratureSplit(facesSplit, nodesSplit, triquadrule);

		} else {
			dTotalFaceArea += CalculateFaceAreaTriangularQuadrature(faces[i], nodes, triquadrule);
		}
	}

	return dTotalFaceArea;
}

///////////////////////////////////////////////////////////////////////////////

double CalculateFaceAreaTriangularQuadratureMultiOrder(
	const Face & face,
	const NodeVector & nodes,
	const TriangularQuadratureRule & triquadrule_loworder,
	const TriangularQuadratureRule & triquadrule_highorder
) {
	double dMaxChordLength = MaxChordLength(face, nodes);

	if (dMaxChordLength < 0.004) {
		return CalculateFaceAreaTriangularQuadrature(face, nodes, triquadrule_loworder);

	} else if (dMaxChordLength < 0.09) {
		return CalculateFaceAreaTriangularQuadrature(face, nodes, triquadrule_highorder);

	} else {
		FaceVector faces;
		faces.push_back(face);
		return CalculateFaceAreaTriangularQuadratureSplit(faces, nodes, triquadrule_highorder);
	}
}

///////////////////////////////////////////////////////////////////////////////
