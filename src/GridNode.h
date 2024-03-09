///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridNode.h
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

#ifndef _GRIDNODE_H_
#define _GRIDNODE_H_

#include "Exception.h"
#include "CoordTransforms.h"

#include <vector>
#include <cstdio>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A coordinate in latitude-longitude space.
///	</summary>
struct LatLon {
	double lat;
	double lon;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A single point in 3D Cartesian geometry.
///	</summary>
class Node {

public:
	///	<summary>
	///		Cartesian x coordinate for this Node.
	///	</summary>
	double x;

	///	<summary>
	///		Cartesian y coordinate for this Node.
	///	</summary>
	double y;

	///	<summary>
	///		Cartesian z coordinate for this Node.
	///	</summary>
	double z;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	Node() :
		x(0.0),
		y(0.0),
		z(0.0)
	{ }

	///	<summary>
	///		Constructor.
	///	</summary>
	Node(
		double a_x,
		double a_y,
		double a_z
	) :
		x(a_x),
		y(a_y),
		z(a_z)
	{ }

	///	<summary>
	///		Copy constructor.
	///	</summary>
	Node(const Node & node) {
		x = node.x;
		y = node.y;
		z = node.z;
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	const Node & operator=(const Node & node) {
		x = node.x;
		y = node.y;
		z = node.z;

		return (*this);
	}

	///	<summary>
	///		Difference between two nodes.
	///	</summary>
	Node operator-(const Node & node) const {
		Node nodeDiff;
		nodeDiff.x = x - node.x;
		nodeDiff.y = y - node.y;
		nodeDiff.z = z - node.z;

		return nodeDiff;
	}

	///	<summary>
	///		Sum of two nodes.
	///	</summary>
	Node operator+(const Node & node) const {
		Node nodeSum;
		nodeSum.x = x + node.x;
		nodeSum.y = y + node.y;
		nodeSum.z = z + node.z;

		return nodeSum;
	}

	///	<summary>
	///		Multiply all coordiantes in node by constant.
	///	</summary>
	Node operator*(const double c) const {
		Node nodeScaled;
		nodeScaled.x = x * c;
		nodeScaled.y = y * c;
		nodeScaled.z = z * c;
		return nodeScaled;
	}

	///	<summary>
	///		Divide all coordinates in node by constant.
	///	</summary>
	Node operator/(const double c) const {
		Node nodeScaled;
		nodeScaled.x = x / c;
		nodeScaled.y = y / c;
		nodeScaled.z = z / c;
		return nodeScaled;
	}

	///	<summary>
	///		Normalize this node to have unit magnitude.
	///	</summary>
	void normalize_in_place() {
		double mag = magnitude();
		x /= mag;
		y /= mag;
		z /= mag;
	}

	///	<summary>
	///		Get a unit magnitude node pointing in the direction of this node.
	///	</summary>
	Node normalized() const {
		Node node(*this);
		node.normalize_in_place();
		return node;
	}

	///	<summary>
	///		Node magnitude.
	///	</summary>
	double magnitude() const {
		return sqrt(x * x + y * y + z * z);
	}

	///	<summary>
	///		Dot product with another node.
	///	</summary>
	double dot(const Node & node) const {
		return x * node.x + y * node.y + z * node.z;
	}

	///	<summary>
	///		Cross product with another node.
	///	</summary>
	Node cross(const Node & node) const {
		Node nodeCross;
		nodeCross.x = y * node.z - z * node.y;
		nodeCross.y = z * node.x - x * node.z;
		nodeCross.z = x * node.y - y * node.x;
		return nodeCross;
	}

	///	<summary>
	///		Cartesian distance from another node.
	///	</summary>
	double distance_from(Node node) const {
		node.x -= x;
		node.y -= y;
		node.z -= z;
		return node.magnitude();
	}

	///	<summary>
	///		Great circle distance (in radians) from another node.
	///		Assumes this Node and argument Node are unit magnitude.
	///	</summary>
	double gc_distance_from(const Node & node) const {
		double coord_distance = distance_from(node);
		if (coord_distance >= 1.0) {
			return M_PI;
		} else {
			return 2.0 * asin(0.5 * coord_distance);
		}
	}

public:
	///	<summary>
	///		Set from a LatLon value in radians.
	///	</summary>
	void from_latlon_rad(const LatLon & latlonrad) {
		x = cos(latlonrad.lon) * cos(latlonrad.lat);
		y = sin(latlonrad.lon) * cos(latlonrad.lat);
		z = sin(latlonrad.lat);
	}

	///	<summary>
	///		Set from a LatLon value in radians.
	///	</summary>
	void from_latlon_deg(const LatLon & latlondeg) {
		double lonrad = DegToRad(latlondeg.lon);
		double latrad = DegToRad(latlondeg.lat);

		x = cos(lonrad) * cos(latrad);
		y = sin(lonrad) * cos(latrad);
		z = sin(latrad);
	}

	///	<summary>
	///		Latitude-longitude coordinates of node (in radians).
	///		Longitude is returned in the interval [0,2pi).
	///		Latitude is returned in the interval [-pi/2,pi/2].
	///	</summary>
	LatLon to_latlon_rad() const {
		LatLon latlon;
		Node node = normalized();

		if (fabs(node.z) < 1.0 - ReferenceTolerance) {
			latlon.lon = atan2(node.y, node.x);
			latlon.lat = asin(node.z);

			if (latlon.lon < 0.0) {
				latlon.lon += 2.0 * M_PI;
			}

		} else if (node.z > 0.0) {
			latlon.lon = 0.0;
			latlon.lat = 0.5 * M_PI;

		} else {
			latlon.lon = 0.0;
			latlon.lat = - 0.5 * M_PI;
		}

		return latlon;
	}

	///	<summary>
	///		Latitude-longitude coordinates of node (in degrees).
	///		Longitude is returned in the interval [0,360.0).
	///		Latitude is returned in the interval [-90,90].
	///	</summary>
	LatLon to_latlon_deg() const {
		LatLon latlon;
		Node node = normalized();

		if (fabs(node.z) < 1.0 - ReferenceTolerance) {
			latlon.lon = RadToDeg(atan2(node.y, node.x));
			latlon.lat = RadToDeg(asin(node.z));

			if (latlon.lon < 0.0) {
				latlon.lon += 360.0;
			}

		} else if (node.z > 0.0) {
			latlon.lon =  0.0;
			latlon.lat = 90.0;

		} else {
			latlon.lon =   0.0;
			latlon.lat = -90.0;
		}

		return latlon;
	}

public:
	///	<summary>
	///		Comparator operator to impose a discrete ordering on Nodes. Note
	///		that one should be careful using this operator as it doesn't
	///		respect floating point tolerance.
	///	</summary>
	bool operator< (const Node & node) const {
		if (x < node.x) {
			return true;
		} else if (x > node.x) {
			return false;
		}

		if (y < node.y) {
			return true;
		} else if (y > node.y) {
			return false;
		}

		if (z < node.z) {
			return true;
		} else if (z > node.z) {
			return false;
		}

		return false;
	}

public:
	///	<summary>
	///		Return a std::string representation of this Node.
	///	</summary>
	std::string to_string() const {
		char szBuffer[80];
		snprintf(szBuffer, 80, "%1.15e %1.15e %1.15e", x, y, z);
		return std::string(szBuffer);
	}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A vector of Nodes.
///	</summary>
typedef std::vector<Node> NodeVector;

///	<summary>
///		A node index.
///	</summary>
typedef int NodeIndex;

///	<summary>
///		A vector for the storage of Node indices.
///	</summary>
typedef std::vector<NodeIndex> NodeIndexVector;

///	<summary>
///		An index indicating this Node is invalid.
///	</summary>
static const NodeIndex InvalidNode = (-1);

///////////////////////////////////////////////////////////////////////////////

#endif

