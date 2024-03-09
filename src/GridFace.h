///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridFace.h
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

#ifndef _GRIDFACE_H_
#define _GRIDFACE_H_

///////////////////////////////////////////////////////////////////////////////

#include "GridNode.h"
#include "GridEdge.h"
#include "Exception.h"

#include <utility>
#include <vector>
#include <map>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A face on a mesh.
///	</summary>
class Face : public NodeIndexVector {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Face(
		size_t node_count = 0
	) :
		NodeIndexVector(node_count, InvalidNode)
	{ }

public:
	///	<summary>
	///		Remove repeated nodes.
	///	</summary>
	void remove_repeated_nodes() {
		if (size() <= 1) {
			return;
		}

		NodeIndexVector new_nodeix;
		for (size_t i = 0; i < size()-1; i++) {
			if ((*this)[i] != (*this)[i+1]) {
				new_nodeix.push_back((*this)[i]);
			}
		}
		if ((*this)[size()-1] != (*this)[0]) {
			new_nodeix.push_back((*this)[size()-1]);
		}

		assign(new_nodeix.begin(), new_nodeix.end());
	}
/*
public:
	///	<summary>
	///		Possible locations of nodes.
	///	</summary>
	enum NodeLocation {
		NodeLocation_Undefined = (-1),
		NodeLocation_Exterior = 0,
		NodeLocation_Default = NodeLocation_Exterior,
		NodeLocation_Interior = 1,
		NodeLocation_Edge = 2,
		NodeLocation_Corner = 3
	};

	///	<summary>
	///		Determine the Edge index corresponding to the given Edge.  If the
	///		Edge is not found an Exception is thrown.
	///	</summary>
	int GetEdgeIndex(const Edge & edge) const;

	///	<summary>
	///		Remove zero Edges (Edges with repeated Node indices)
	///	</summary>
	void RemoveZeroEdges(); /// TODO: REPLACED WITH remove_repeated_nodes
*/
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A vector of Faces.
///	</summary>
typedef std::vector<Face> FaceVector;

///	<summary>
///		A face index.
///	</summary>
typedef int FaceIndex;

///	<summary>
///		An index indicating this Face is invalid.
///	</summary>
static const FaceIndex InvalidFace = (-1);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A pair of face indices, typically on opposite sides of an Edge.
///	</summary>
class FacePair : public std::pair<FaceIndex, FaceIndex> {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	FacePair(
		FaceIndex a_first = InvalidFace,
		FaceIndex a_second = InvalidFace
	) :
		std::pair<FaceIndex, FaceIndex>(a_first, a_second)
	{ }

	///	<summary>
	///		Add a face index to this FacePair.
	///	</summary>
	void add_face(int ixFace) {
		if (first == InvalidFace) {
			first = ixFace;
		} else if (second == InvalidFace) {
			second = ixFace;
		} else {
			_EXCEPTIONT("FacePair already has a full set of Faces");
		}
	}

	///	<summary>
	///		Does this FacePair have a complete set of Faces?
	///	</summary>
	bool is_complete() const {
		return ((first != InvalidFace) && (second != InvalidFace));
	}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A single Edge and the Faces on either side of the Edge.
///	</summary>
typedef std::pair<Edge, FacePair> EdgeMapPair;

///	<summary>
///		A map between Edges and the Faces that are on either side of the Edge.
///	</summary>
typedef std::map<Edge, FacePair> EdgeMap;

///////////////////////////////////////////////////////////////////////////////

#endif

