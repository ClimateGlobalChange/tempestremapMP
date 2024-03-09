///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridEdge.h
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

#ifndef _GRIDEDGE_H_
#define _GRIDEDGE_H_

///////////////////////////////////////////////////////////////////////////////

#include <utility>
#include <vector>
#include <set>
#include <cstdio>

#include "GridNode.h"

///////////////////////////////////////////////////////////////////////////////

enum EdgeType {
	EdgeType_GreatCircleArc = 0,
	EdgeType_Default = EdgeType_GreatCircleArc,
	EdgeType_ConstantLatitude = 1
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An edge connects two nodes.
///	</summary>
class Edge : public std::pair<NodeIndex, NodeIndex> {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Edge(
		NodeIndex nodeix_first = InvalidNode,
		NodeIndex nodeix_second = InvalidNode
	) {
		first = nodeix_first;
		second = nodeix_second;
	}

	///	<summary>
	///		Copy constructor.
	///	</summary>
	Edge(const Edge & edge) {
		first = edge.first;
		second = edge.second;
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	Edge & operator=(const Edge & edge) {
		first = edge.first;
		second = edge.second;
		return (*this);
	}

	///	<summary>
	///		Flip the Edge.
	///	</summary>
	void flip() {
		std::swap(first, second);
	}

	///	<summary>
	///		Order this edge in place, so first <= second.
	///	</summary>
	void order_in_place() {
		if (first > second) {
			flip();
		}
	}

	///	<summary>
	///		Return an equivalent Edge which is ordered with first <= second.
	///	</summary>
	Edge ordered() const {
		if (first < second) {
			return (*this);
		} else {
			return Edge(second, first);
		}
	}

	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator<(const Edge & edge) const {
		Edge edgeThisOrdered = ordered();
		Edge edgeOtherOrdered = edge.ordered();

		if (edgeThisOrdered.first < edgeOtherOrdered.first) {
			return true;
		} else if (edgeThisOrdered.first > edgeOtherOrdered.first) {
			return false;
		} else if (edgeThisOrdered.second < edgeOtherOrdered.second) {
			return true;
		} else {
			return false;
		}
	}

	///	<summary>
	///		Equality operator.
	///	</summary>
	bool operator==(const Edge & edge) const {
		if ((first == edge.first) && (second == edge.second)) {
			return true;
		}
		if ((first == edge.second) && (second == edge.first)) {
			return true;
		}
		return false;
	}

	///	<summary>
	///		Inequality operator.
	///	</summary>
	bool operator!=(const Edge & edge) const {
		return !((*this) == edge);
	}

	///	<summary>
	///		Return the node that is shared between segments.
	///	</summary>
	NodeIndex CommonNode(
		const Edge & edge
	) const {
		if ((first == edge.first) || (first == edge.second)) {
			return first;
		} else if ((second == edge.first) || (second == edge.second)) {
			return second;
		} else {
			return InvalidNode;
		}
	}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A vector of Edges.
///	</summary>
typedef std::vector<Edge> EdgeVector;

///	<summary>
///		A set of Edges.
///	</summary>
typedef std::set<Edge> EdgeSet;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An Edge that contains many interior nodes.
///	</summary>
class MultiEdge : public std::vector<NodeIndex> {

public:
	///	<summary>
	///		Flip the edge.
	///	</summary>
	MultiEdge flipped() const {
		MultiEdge edgeFlip;
		for (int i = size()-1; i >= 0; i--) {
			edgeFlip.push_back((*this)[i]);
		}
		return edgeFlip;
	}

};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A vector of MultiEdges.
///	</summary>
typedef std::vector<MultiEdge> MultiEdgeVector;

///////////////////////////////////////////////////////////////////////////////

#endif

