///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridElements.h
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

#ifndef _GRIDELEMENTS_H_
#define _GRIDELEMENTS_H_

///////////////////////////////////////////////////////////////////////////////

#include "GridNode.h"
#include "GridEdge.h"
#include "GridFace.h"
#include "Exception.h"
#include "netcdfcpp.h"

#include <vector>
#include <set>
#include <map>
#include <string>
#include <cmath>
#include <cassert>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		All Face indices that contain the given Node.
///	</summary>
typedef std::vector< std::set<FaceIndex> > ReverseNodeArray;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A geometric mesh, consisting of Nodes and Faces that connect those Nodes.
///	</summary>
class Mesh {

public:
	///	<summary>
	///		Flag indicating this Mesh contains concave Faces.
	///	</summary>
	bool contains_concave_faces;

public:
	///	<summary>
	///		Vector of Nodes for this mesh.
	///	</summary>
	NodeVector nodes;

	///	<summary>
	///		Vector of Faces for this mesh.
	///	<summary>
	FaceVector faces;

public:
	///	<summary>
	///		Filename for this mesh.
	///	</summary>
	std::string filename;

	///	<summary>
	///		Vector of Face areas.
	///	</summary>
	std::vector<double> face_areas;

	///	<summary>
	///		EdgeMap for this mesh.
	///	</summary>
	EdgeMap edge_map;

	///	<summary>
	///		ReverseNodeArray for this mesh.
	///	</summary>
	ReverseNodeArray reverse_node_array;

public:
	///	<summary>
	///		Vector of source mesh Face indices (assumed to be in
	///		increasing order).
	///	</summary>
	std::vector<FaceIndex> source_face_ix;

	///	<summary>
	///		Vector of target mesh Face indices.
	///	</summary>
	std::vector<FaceIndex> target_face_ix;

	///	<summary>
	///		Vector storing mask variable for this mesh.
	///	</summary>
	std::vector<int> face_mask;

	///	<summary>
	///		Indices of the original Faces for this mesh (for use when
	///		the original mesh has been subdivided).
	///	</summary>
	std::vector<FaceIndex> multi_face_ix;

	///	<summary>
	///		Grid dimensions.
	///	</sumamry>
	std::vector<int> grid_dim_sizes;

	///	<summary>
	///		Grid dimension names.
	///	</sumamry>
	std::vector<std::string> grid_dim_names;

public:
	///	<summary>
	///		Constructor with input mesh filename.
	///	</summary>
	Mesh(
		const std::string & a_filename = "",
		bool a_contains_concave_faces = false
	)  :
		contains_concave_faces(a_contains_concave_faces)
	{
		if (a_filename != "") {
			read(a_filename);
		}
	}

public:
	///	<summary>
	///		Clear the contents of the mesh.
	///	</summary>
	void clear();

	///	<summary>
	///		Calculate total Mesh area.
	///	</summary>
	double calculate_total_area();

	///	<summary>
	///		Calculate Face areas.
	///	</summary>
	void calculate_face_areas();

	///	<summary>
	///		Calculate Face areas from an overlap mesh.
	///	</summary>
	void assign_face_areas_from_overlap_mesh_source_ix(
		const Mesh & meshOverlap
	);

	///	<summary>
	///		Calculate Face areas from an overlap mesh.
	///	</summary>
	void assign_face_areas_from_overlap_mesh_target_ix(
		const Mesh & meshOverlap
	);

	///	<summary>
	///		Build the EdgeMap from the NodeVector and FaceVector.
	///	</summary>
	void build_edge_map();

	///	<summary>
	///		Check if mesh has holes. Note that build_edge_map() needs to be
	///		called before this function can be called.
	///	</summary>
	bool has_holes() const;

	///	<summary>
	///		Build the ReverseNodeArray from the NodeVector and FaceVector.
	///	</summary>
	void build_reverse_node_array();

	///	<summary>
	///		Remove coincident nodes from the Mesh and adjust indices in Faces.
	///	</summary>
	void remove_coincident_nodes(
		bool fVerbose = true
	);

public:
	///	<summary>
	///		Swap the source and target meshes. This further requires reordering
	///		the Faces in source_face_ix and target_face_ix.
	///	</summary>
	void exchange_source_and_target_mesh();

	///	<summary>
	///		Write the mesh to a NetCDF file in Exodus format.
	///	</summary>
	void write_as_exodus(
		const std::string & strFilename,
		NcFile::FileFormat eFileFormat = NcFile::Netcdf4
	) const;

	///	<summary>
	///		Write the mesh to a NetCDF file in SCRIP format.
	///	</summary>
	void write_as_scrip(
		const std::string & strFilename,
		NcFile::FileFormat eFileFormat = NcFile::Netcdf4
	) const;

	///	<summary>
	///		Write the mesh to a NetCDF file in UGRID format.
	///	</summary>
	void write_as_ugrid(
		const std::string & strFilename,
		NcFile::FileFormat eFileFormat = NcFile::Netcdf4
	) const;

	///	<summary>
	///		Read the mesh from a NetCDF file.
	///	</summary>
	void read(const std::string & strFilename);

	///	<summary>
	///		Remove zero edges from all Faces.
	///	</summary>
	void remove_zero_edges();

	///	<summary>
	///		Validate the Mesh.
	///	</summary>
	void validate() const;
};

///////////////////////////////////////////////////////////////////////////////
/*
///	<summary>
///		Location data returned from FindFaceFromNode()
///		Generate a PathSegmentVector describing the path around the face
///		ixCurrentFirstFace.
///	</summary>
struct FindFaceStruct {
	
	///	<summary>
	///		A vector of face indices indicating possible Faces.
	///	</summary>
	std::vector<int> vecFaceIndices;

	///	<summary>
	///		A vector of locations on each Face.  If loc is NodeLocation_Corner,
	///		this corresponds to the associated corner of the Face.  If loc
	///		is NodeLocation_Edge, this corresponds to the associated Edge of
	///		the Face.  If loc is NodeLocation_Interior, this value is
	///		undefined.
	///	</summary>
	std::vector<int> vecFaceLocations;

	///	<summary>
	///		The NodeLocation where this Node lies.
	///	</summary>
	Face::NodeLocation loc;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if an edge is positively oriented
///		(aligned with increasing longitude).
///	</summary>
bool IsPositivelyOrientedEdge(
	const Node & nodeBegin,
	const Node & nodeEnd
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get the local direction vector along the surface of the sphere
///		for the given edge.
///	</summary>
void GetLocalDirection(
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Node & nodeRef,
	const Edge::Type edgetype,
	Node & nodeDir
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get the local direction vector along the surface of the sphere
///		for the given edge.
///	</summary>
void GetLocalDirection(
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	Node & nodeDir
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		For all Nodes on meshSecond that are "nearby" a Node on meshFirst,
///		set the Node equal to the meshFirst Node.
///	</summary>
void EqualizeCoincidentNodes(
	const Mesh & meshFirst,
	Mesh & meshSecond
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Build the mapping function for nodes on meshSecond which are
///		coincident with nodes on meshFirst.
///	</summary>
///	<returns>
///		The number of coincident nodes on meshSecond.
///	</returns>
int BuildCoincidentNodeVector(
	const Mesh & meshFirst,
	const Mesh & meshSecond,
	std::vector<int> & vecSecondToFirstCoincident
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the Jacobian (infinitesmal local area element) for a
///		spherical triangle.  Variables dA and dB specify the coordinates
///		within the spherical triangle:
///		  (dA,dB)=(0,0) corresponds to node1
///		  (dA,dB)=(1,0) corresponds to node2
///		  (dA,dB)=(0,1)=(1,1) corresponds to node3
///	</summary>
double CalculateSphericalTriangleJacobian(
	const Node & node1,
	const Node & node2,
	const Node & node3,
	double dA,
	double dB,
	Node * pnode = NULL
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the Jacobian (infinitesmal local area element) for a
///		spherical triangle.  Variables dA and dB specify the barycentric
///		coordinates within the spherical triangle:
///		  (dA,dB)=(0,0) corresponds to node3
///		  (dA,dB)=(1,0) corresponds to node2
///		  (dA,dB)=(0,1) corresponds to node1
///	</summary>
double CalculateSphericalTriangleJacobianBarycentric(
	const Node & node1,
	const Node & node2,
	const Node & node3,
	double dA,
	double dB,
	Node * pnode = NULL
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the area of a Face using quadrature.
///	</summary>
double CalculateFaceAreaQuadratureMethod(
	const Face & face,
	const NodeVector & nodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the area of a Face using Karney's method (may be poorly
///		conditioned at higher resolutions).
///	</summary>
double CalculateFaceAreaKarneysMethod(
	const Face & face,
	const NodeVector & nodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Average the nodes of a face.
///	</summary>
Node AverageFaceNodes(
	const Face & face,
	const NodeVector & nodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Check if the specified Face is concave.
///	</summary>
bool IsFaceConcave(
	const Face & face,
	const NodeVector & nodes
);
*/
///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the area of a single Face, detecting concave Faces.
///	</summary>
double CalculateFaceArea_Concave(
	const Face & face,
	const NodeVector & nodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the area of a single Face.
///	</summary>
double CalculateFaceArea(
	const Face & face,
	const NodeVector & nodes
);
/*
/// <summary>
///     Calculate triangle area, quadrature.
/// </summary>
double CalculateTriangleAreaQuadratureMethod(
	Node &node1,
	Node &node2,
	Node &node3
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Convert a concave Face into a Convex face.  Based on routine by
///		Mark Bayazit (https://mpen.ca/406/bayazit).
///	</summary>
///	<returns>
///		true if the Face is convex and has been removed from the mesh.faces
///		vector.
///	</returns>
bool ConvexifyFace(
	Mesh & meshin,
	Mesh & meshout,
	int iFace,
	bool fRemoveConcaveFaces,
	bool fVerbose = false
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Convert all concave Faces into the Mesh into Concave faces via
///		subdivision.
///	</summary>
void ConvexifyMesh(
	Mesh & mesh,
	bool fVerbose = false
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Convert concave Mesh meshin to convex Mesh meshout by dividing
///		Faces and populating the MultiFaceMap.
///	</summary>
void ConvexifyMesh(
	Mesh & meshin,
	Mesh & meshout,
	bool fVerbose = false
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find a node within the specified quadrilateral.
///	</summary>
inline Node InterpolateQuadrilateralNode(
	const Node & node0,
	const Node & node1,
	const Node & node2,
	const Node & node3,
	double dA,
	double dB
) {
	Node nodeRef;

	nodeRef.x =
		  (1.0 - dA) * (1.0 - dB) * node0.x
		+        dA  * (1.0 - dB) * node1.x
		+        dA  *        dB  * node2.x
		+ (1.0 - dA) *        dB  * node3.x;

	nodeRef.y =
		  (1.0 - dA) * (1.0 - dB) * node0.y
		+        dA  * (1.0 - dB) * node1.y
		+        dA  *        dB  * node2.y
		+ (1.0 - dA) *        dB  * node3.y;

	nodeRef.z =
		  (1.0 - dA) * (1.0 - dB) * node0.z
		+        dA  * (1.0 - dB) * node1.z
		+        dA  *        dB  * node2.z
		+ (1.0 - dA) *        dB  * node3.z;

	double dMag = sqrt(
		  nodeRef.x * nodeRef.x
		+ nodeRef.y * nodeRef.y
		+ nodeRef.z * nodeRef.z);

	nodeRef.x /= dMag;
	nodeRef.y /= dMag;
	nodeRef.z /= dMag;

	return nodeRef;
}

///////////////////////////////////////////////////////////////////////////////

void Dual(
	Mesh & mesh
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Construct the local dual face around a given node
///	</summary>
void ConstructLocalDualFace(
	const Mesh & mesh,
	NodeVector & meshCenters,
	int & iNodeX,
	Face & faceLocalDual,
	NodeVector & nodesFaceLocal
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the maximum edge length of an element.
///	</summary>
double MaxEdgeLength(
	const Face & face,
	const NodeVector & nodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculates the area of an element; to be used in adaptive element area calculation.
///	</summary>
double CalculateFaceAreaTriQuadrature(
	const Face & face,
	const NodeVector & nodes,
	int nOrder
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Adaptive area calculation that refines an element based on the maximum edge length.
///	</summary>
double CalculateFaceAreaTriQuadratureSplit(
	const FaceVector & faces,
	const NodeVector & nodes,
	int & nOrder
);
*/
///////////////////////////////////////////////////////////////////////////////

#endif

