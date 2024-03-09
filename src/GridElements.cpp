///////////////////////////////////////////////////////////////////////////////
///
///	\file	GridElements.cpp
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

#include "Defines.h"
#include "Announce.h"
#include "GridElements.h"
#include "GridUtilities.h"
#include "STLStringHelper.h"
#include "netcdfcpp.h"

#include "triangle.h"
#include "kdtree.h"
#include "array2d.h"

///////////////////////////////////////////////////////////////////////////////
/// Mesh
///////////////////////////////////////////////////////////////////////////////

void Mesh::clear() {
	nodes.clear();
	faces.clear();

	filename = "";

	source_face_ix.clear();
	target_face_ix.clear();

	face_areas.clear();
	face_mask.clear();
	edge_map.clear();
	reverse_node_array.clear();
	multi_face_ix.clear();
	grid_dim_sizes.clear();
	grid_dim_names.clear();
}

///////////////////////////////////////////////////////////////////////////////

double Mesh::calculate_total_area() {

	// Make sure face areas have been calculated
	if (face_areas.size() != faces.size()) {
		calculate_face_areas();
	}

	// Zero area
	if (faces.size() == 0) {
		return 0.0;
	}

	// Calculate accumulated area carefully
	static const int Jump = 10;

	std::vector<double> accumulated_areas(face_areas);

	for (;;) {
		if (accumulated_areas.size() == 1) {
			break;
		}
		for (int i = 0; i <= (accumulated_areas.size()-1) / Jump; i++) {
			int ixRef = Jump * i;
			accumulated_areas[i] = accumulated_areas[ixRef];
			for (int j = 1; j < Jump; j++) {
				if (ixRef + j >= accumulated_areas.size()) {
					break;
				}
				accumulated_areas[i] += accumulated_areas[ixRef + j];
			}
		}
		accumulated_areas.resize((accumulated_areas.size()-1) / Jump + 1);
	}

	return accumulated_areas[0];
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::calculate_face_areas() {

	// If no faces, then zero area
	if (faces.size() == 0) {
		return;
	}

	// Calculate areas
	int nCountSmallFaces = 0;
	face_areas.resize(faces.size());

	// Calculate the area of each Face (for sometimes concave faces)
	if (contains_concave_faces) {
		for (int i = 0; i < faces.size(); i++) {
			face_areas[i] = CalculateFaceArea_Concave(faces[i], nodes);
			if (face_areas[i] < ReferenceTolerance) {
				nCountSmallFaces++;
			}
		}

	// Calculate the area of each Face (for convex faces)
	} else {
		TriangularQuadratureRule triquadrule_loworder(4);
		TriangularQuadratureRule triquadrule_highorder(8);

		for (int i = 0; i < faces.size(); i++) {
			face_areas[i] =
				CalculateFaceAreaTriangularQuadratureMultiOrder(
					faces[i],
					nodes,
					triquadrule_loworder,
					triquadrule_highorder);

			if (face_areas[i] < ReferenceTolerance) {
				nCountSmallFaces++;
			}
		}
	}

	if (nCountSmallFaces != 0) {
		Announce("WARNING: %i small faces found", nCountSmallFaces);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::assign_face_areas_from_overlap_mesh_source_ix(
	const Mesh & meshOverlap
) {
	if (meshOverlap.face_areas.size() != meshOverlap.faces.size()) {
		_EXCEPTIONT("meshOverlap face_areas is not defined");
	}
	if (meshOverlap.source_face_ix.size() != meshOverlap.faces.size()) {
		_EXCEPTIONT("meshOverlap source_face_ix is not defined");
	}

	// Set all Face areas to zero
	face_areas.resize(faces.size());
	std::fill(face_areas.begin(), face_areas.end(), 0.0);

	// Loop over all Faces in meshOverlap and use to update area
	for (int i = 0; i < meshOverlap.faces.size(); i++) {
		int ixSourceFace = meshOverlap.source_face_ix[i];

		if ((ixSourceFace < 0) || (ixSourceFace >= face_areas.size())) {
			_EXCEPTIONT("Overlap Mesh source_face_ix contains invalid Face index");
		}

		face_areas[ixSourceFace] += meshOverlap.face_areas[i];
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::assign_face_areas_from_overlap_mesh_target_ix(
	const Mesh & meshOverlap
) {
	if (meshOverlap.face_areas.size() != meshOverlap.faces.size()) {
		_EXCEPTIONT("meshOverlap face_areas is not defined");
	}
	if (meshOverlap.target_face_ix.size() != meshOverlap.faces.size()) {
		_EXCEPTIONT("meshOverlap target_face_ix is not defined");
	}

	// Set all Face areas to zero
	face_areas.resize(faces.size());
	std::fill(face_areas.begin(), face_areas.end(), 0.0);

	// Loop over all Faces in meshOverlap and use to update area
	for (int i = 0; i < meshOverlap.faces.size(); i++) {
		int ixTargetFace = meshOverlap.target_face_ix[i];

		if ((ixTargetFace < 0) || (ixTargetFace >= face_areas.size())) {
			_EXCEPTIONT("Overlap Mesh target_face_ix contains invalid Face index");
		}

		face_areas[ixTargetFace] += meshOverlap.face_areas[i];
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::build_edge_map() {
	edge_map.clear();

	for (size_t f = 0; f < faces.size(); f++) {
		const Face & face = faces[f];
		size_t nFaceNodes = face.size();

		for (size_t k = 0; k < nFaceNodes; k++) {
			if (face[k] == face[(k+1)%nFaceNodes]) {
				continue;
			}

			Edge edge(face[k], face[(k+1)%nFaceNodes]);

			auto itEdgeMapPair = edge_map.find(edge);
			if (itEdgeMapPair == edge_map.end()) {
				edge_map.insert(EdgeMapPair(edge, FacePair(f,InvalidFace)));

			} else if (itEdgeMapPair->second.second == InvalidFace) {
				itEdgeMapPair->second = f;

			} else {
				_EXCEPTIONT("FacePair in EdgeMap already has a full set of Faces");
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

bool Mesh::has_holes() const {
	if (edge_map.size() == 0) {
		_EXCEPTIONT("EdgeMap must be initialized before calling has_holes()");
	}

	// Check for any elements of edge_map that don't have
	// a complete set of Faces
	for (auto it = edge_map.begin(); it != edge_map.end(); it++) {
		if (!it->second.is_complete()) {
			return true;
		}
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::build_reverse_node_array() {
	reverse_node_array.clear();
	reverse_node_array.resize(nodes.size());
	
	for (FaceIndex f = 0; f < faces.size(); f++) {
		const Face & face = faces[f];
		for (const NodeIndex ix: face) {
			reverse_node_array[ix].insert(f);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::exchange_source_and_target_mesh() {

	// Verify all vectors are the same size
	if (faces.size() != source_face_ix.size()) {
		_EXCEPTIONT("Mismatch between size of faces and source_face_ix");
	}
	if (faces.size() != target_face_ix.size()) {
		_EXCEPTIONT("Mismatch between size of faces and target_face_ix");
	}

	// Reorder vectors
	FaceVector facesOld = faces;

	std::vector<FaceIndex> source_face_ixOld = source_face_ix;

	// Reordering map
	std::multimap<FaceIndex,size_t> multimapReorder;
	for (size_t i = 0; i < target_face_ix.size(); i++) {
		multimapReorder.insert(std::pair<FaceIndex,FaceIndex>(target_face_ix[i], i));
	}

	// Apply reordering
	faces.clear();
	source_face_ix.clear();
	target_face_ix.clear();

	for (auto itReorder = multimapReorder.begin(); itReorder != multimapReorder.end(); itReorder++) {
		faces.push_back(facesOld[itReorder->second]);
		source_face_ix.push_back(itReorder->first);
		target_face_ix.push_back(source_face_ixOld[itReorder->second]);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::remove_coincident_nodes(
	bool fVerbose
) {

	// Use kdtree to find shortest distance to other nodes
	kdtree * kdt = kd_create(3);
	if (kdt == nullptr) {
		_EXCEPTIONT("Error calling kd_create(3)");
	}

	std::vector<NodeIndex> vecNodeIndex;
	std::vector<NodeIndex> vecUniques;

	vecNodeIndex.reserve(nodes.size());
	vecUniques.reserve(nodes.size());

	kd_insert3(kdt, nodes[0].x, nodes[0].y, nodes[0].z, (void*)(0));
	vecNodeIndex.push_back(0);
	vecUniques.push_back(0);

	for (size_t k = 1; k < nodes.size(); k++) {
		const Node & node = nodes[k];

		kdres * kdresNearest = kd_nearest3(kdt, node.x, node.y, node.z);
		if (kdresNearest == NULL) {
			_EXCEPTIONT("kd_nearest3() failed");
		}
		Node nodeNearest;
		size_t ixNodeNearest = (size_t)(kd_res_item3(kdresNearest, &(nodeNearest.x), &(nodeNearest.y), &(nodeNearest.z)));
		kd_res_free(kdresNearest);

		if (node.distance_from(nodeNearest) < ReferenceTolerance) {
			vecNodeIndex.push_back((NodeIndex)ixNodeNearest);
		} else {
			vecNodeIndex.push_back((NodeIndex)(vecUniques.size()));
			vecUniques.push_back((NodeIndex)k);
			kd_insert3(kdt, node.x, node.y, node.z, (void*)(k));
		}
	}

	kd_free(kdt);

	// Number of uniques 
	if (vecUniques.size() == nodes.size()) {
		return;
	}

	if (fVerbose) {
		Announce("%i duplicate nodes detected", nodes.size() - vecUniques.size());
	}

	// Remove duplicates from nodes vector
	{
		NodeVector nodesOld = nodes;
		nodes.resize(vecUniques.size());
		for (size_t i = 0; i < vecUniques.size(); i++) {
			nodes[i] = nodesOld[vecUniques[i]];
		}
	}

	// Adjust node indices in Faces
	for (Face & face : faces) {
		for (NodeIndex & ixNode : face) {
			ixNode = vecNodeIndex[ixNode];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::write_as_exodus(
	const std::string & strFile,
	NcFile::FileFormat eFileFormat
) const {
	const int ParamFour = 4;
	const int ParamLenString = 33;

	// Temporarily change error reporting
	NcError error_temp(NcError::verbose_fatal);

	// Determine block sizes
	std::vector<int> vecBlockSizes;
	std::vector<int> vecBlockSizeFaces;
	{
		std::map<int, int> mapBlockSizes;
		std::map<int, int>::iterator iterBlockSize;
		int iBlock;
		char szBuffer[ParamLenString];

		for (int i = 0; i < faces.size(); i++) {
			iterBlockSize = mapBlockSizes.find(faces[i].size());

			if (iterBlockSize == mapBlockSizes.end()) {
				mapBlockSizes.insert(
					std::pair<int,int>(faces[i].size(), 1));
			} else {
				(iterBlockSize->second)++;
			}
		}

		vecBlockSizes.resize(mapBlockSizes.size());
		vecBlockSizeFaces.resize(mapBlockSizes.size());

		AnnounceStartBlock("Nodes per element");
		iterBlockSize = mapBlockSizes.begin();
		iBlock = 1;
		for (; iterBlockSize != mapBlockSizes.end(); iterBlockSize++) {
			vecBlockSizes[iBlock-1] = iterBlockSize->first;
			vecBlockSizeFaces[iBlock-1] = iterBlockSize->second;

			Announce("Block %i (%i nodes): %i",
				iBlock, vecBlockSizes[iBlock-1], vecBlockSizeFaces[iBlock-1]);

			iBlock++;
		}
		AnnounceEndBlock(NULL);
	}

	// Output to a NetCDF Exodus file
	NcFile ncOut(strFile.c_str(), NcFile::Replace, NULL, 0, eFileFormat);
	if (!ncOut.is_valid()) {
		_EXCEPTION1("Unable to open grid file \"%s\" for writing",
			strFile.c_str());
	}

	// Auxiliary Exodus dimensions
	NcDim * dimLenString = ncOut.add_dim("len_string", ParamLenString);
	NcDim * dimLenLine = ncOut.add_dim("len_line", 81);
	NcDim * dimFour = ncOut.add_dim("four", ParamFour);
	NcDim * dimTime = ncOut.add_dim("time_step");
	NcDim * dimDimension = ncOut.add_dim("num_dim", 3);

	// Number of nodes
	int nNodeCount = nodes.size();
	NcDim * dimNodes = ncOut.add_dim("num_nodes", nNodeCount);

	// Number of elements
	int nFaceCount = faces.size();
	NcDim * dimElements = ncOut.add_dim("num_elem", nFaceCount);

	// Other dimensions
	NcDim * dimNumQARec = ncOut.add_dim("num_qa_rec", 1);

	// Global attributes
	ncOut.add_att("api_version", 5.00f);
	ncOut.add_att("version", 5.00f);
	ncOut.add_att("floating_point_word_size", 8);
	ncOut.add_att("file_size", 0);

	// Current time
	{
		time_t t = time(0);
		struct tm * timestruct = localtime(&t);

		char szDate[ParamLenString];
		char szTime[ParamLenString];

		strftime(szDate, sizeof(szDate), "%m/%d/%Y", timestruct);
		strftime(szTime, sizeof(szTime), "%X", timestruct);

		char szTitle[128];
		snprintf(szTitle, 128, "tempest(%s) %s: %s", strFile.c_str(), szDate, szTime);
		ncOut.add_att("title", szTitle);

		// Time_whole (unused)
		NcVar * varTimeWhole = ncOut.add_var("time_whole", ncDouble, dimTime);
		if (varTimeWhole == NULL) {
			_EXCEPTIONT("Error creating variable \"time_whole\"");
		}

		// QA records
		char szQARecord[ParamFour][ParamLenString]
			= {"Tempest", "14.0", "", ""};

		strcpy(szQARecord[2], szDate);
		strcpy(szQARecord[3], szTime);

		NcVar * varQARecords =
			ncOut.add_var("qa_records", ncChar,
				dimNumQARec, dimFour, dimLenString);

		if (varQARecords == NULL) {
			_EXCEPTIONT("Error creating variable \"qa_records\"");
		}

		varQARecords->set_cur(0, 0, 0);
		varQARecords->put(&(szQARecord[0][0]), 1, 4, ParamLenString);
	}

	// Coordinate names
	{
		char szCoordNames[3][ParamLenString] = {"x", "y", "z"};

		NcVar * varCoordNames =
			ncOut.add_var("coor_names", ncChar, dimDimension, dimLenString);

		if (varCoordNames == NULL) {
			_EXCEPTIONT("Error creating variable \"coor_names\"");
		}

		varCoordNames->set_cur(0, 0, 0);
		varCoordNames->put(&(szCoordNames[0][0]), 3, ParamLenString);
	}

	// Element blocks
	NcDim * dimNumElementBlocks =
		ncOut.add_dim("num_el_blk", vecBlockSizes.size());

	if (dimNumElementBlocks == NULL) {
		_EXCEPTIONT("Error creating dimension \"num_el_blk\"");
	}

	std::vector<NcDim *> vecElementBlockDim;
	std::vector<NcDim *> vecNodesPerElementDim;
	std::vector<NcDim *> vecAttBlockDim;

	vecElementBlockDim.resize(vecBlockSizes.size());
	vecNodesPerElementDim.resize(vecBlockSizes.size());
	vecAttBlockDim.resize(vecBlockSizes.size());

	char szBuffer[ParamLenString];
	for (int n = 0; n < vecBlockSizes.size(); n++) {
		snprintf(szBuffer, ParamLenString, "num_el_in_blk%i", n+1);
		vecElementBlockDim[n] =
			ncOut.add_dim(szBuffer, vecBlockSizeFaces[n]);

		if (vecElementBlockDim[n] == NULL) {
			_EXCEPTION1("Error creating dimension \"%s\"", szBuffer);
		}

		snprintf(szBuffer, ParamLenString, "num_nod_per_el%i", n+1);
		vecNodesPerElementDim[n] =
			ncOut.add_dim(szBuffer, vecBlockSizes[n]);

		if (vecNodesPerElementDim[n] == NULL) {
			_EXCEPTION1("Error creating dimension \"%s\"", szBuffer);
		}

		snprintf(szBuffer, ParamLenString, "num_att_in_blk%i", n+1);
		vecAttBlockDim[n] =
			ncOut.add_dim(szBuffer, 1);

		if (vecAttBlockDim[n] == NULL) {
			_EXCEPTION1("Error creating dimension \"%s\"", szBuffer);
		}
	}

	// Element block names
	{
		NcVar * varElementBlockNames =
			ncOut.add_var("eb_names", ncChar,
				dimNumElementBlocks, dimLenString);

		if (varElementBlockNames == NULL) {
			_EXCEPTIONT("Error creating dimension \"eb_names\"");
		}
	}

	// Element block status and property
	{
		std::vector<int> vecStatus;
		std::vector<int> vecProp;
		vecStatus.resize(vecBlockSizes.size());
		vecProp.resize(vecBlockSizes.size());

		for (int n = 0; n < vecBlockSizes.size(); n++) {
			vecStatus[n] = 1;
			vecProp[n] = n+1;
		}

		NcVar * varElementBlockStatus =
			ncOut.add_var("eb_status", ncInt, dimNumElementBlocks);

		if (varElementBlockStatus == NULL) {
			_EXCEPTIONT("Error creating variable \"eb_status\"");
		}

		varElementBlockStatus->put(&(vecStatus[0]), vecBlockSizes.size());

		NcVar * varElementProperty =
			ncOut.add_var("eb_prop1", ncInt, dimNumElementBlocks);

		if (varElementProperty == NULL) {
			_EXCEPTIONT("Error creating variable \"eb_prop1\"");
		}

		varElementProperty->put(&(vecProp[0]), vecBlockSizes.size());
		varElementProperty->add_att("name", "ID");
	}

	// Attributes
	{
		for (int n = 0; n < vecBlockSizes.size(); n++) {
			std::vector<double> dAttrib;
			dAttrib.resize(vecBlockSizeFaces[n]);
			for (int i = 0; i < vecBlockSizeFaces[n]; i++) {
				dAttrib[i] = 1.0;
			}

			char szAttribName[ParamLenString];
			snprintf(szAttribName, ParamLenString, "attrib%i", n+1);

			NcVar * varAttrib =
				ncOut.add_var(
					szAttribName, ncDouble,
					vecElementBlockDim[n],
					vecAttBlockDim[n]);

			if (varAttrib == NULL) {
				_EXCEPTION1("Error creating variable \"%s\"", szAttribName);
			}

			varAttrib->set_cur((long)0, (long)0);
			varAttrib->put(&(dAttrib[0]), vecBlockSizeFaces[n], 1);
		}
	}

	// Face-specific variables
	{
		// Face nodes (1-indexed)
		std::vector<NcVar*> vecConnectVar;
		vecConnectVar.resize(vecBlockSizes.size());

		std::vector< tmp::array2d<int> > vecConnect;
		vecConnect.resize(vecBlockSizes.size());

		// Number of elements added to each connectivity array
		std::vector<int> vecConnectCount;
		vecConnectCount.resize(vecBlockSizes.size());

		// Global ids
		std::vector<NcVar*> vecGlobalIdVar;
		vecGlobalIdVar.resize(vecBlockSizes.size());

		std::vector< std::vector<int> > vecGlobalId;
		vecGlobalId.resize(vecBlockSizes.size());

		// Edge types
		std::vector<NcVar*> vecEdgeTypeVar;
		vecEdgeTypeVar.resize(vecBlockSizes.size());

		std::vector< tmp::array2d<int> > vecEdgeType;
		vecEdgeType.resize(vecBlockSizes.size());

		// Parent on source mesh
		std::vector<NcVar*> vecFaceParentAVar;
		vecFaceParentAVar.resize(vecBlockSizes.size());

		std::vector< std::vector<int> > vecFaceParentA;
		vecFaceParentA.resize(vecBlockSizes.size());

		// Parent on target mesh
		std::vector<NcVar*> vecFaceParentBVar;
		vecFaceParentBVar.resize(vecBlockSizes.size());

		std::vector< std::vector<int> > vecFaceParentB;
		vecFaceParentB.resize(vecBlockSizes.size());

		// Initialize block-local storage arrays and create output variables
		for (int n = 0; n < vecBlockSizes.size(); n++) {
			vecConnect[n].resize(vecBlockSizeFaces[n], vecBlockSizes[n], 0);
			vecGlobalId[n].resize(vecBlockSizeFaces[n], 0);
			vecEdgeType[n].resize(vecBlockSizeFaces[n], vecBlockSizes[n], 0);
	
			char szConnectVarName[ParamLenString];
			snprintf(szConnectVarName, ParamLenString, "connect%i", n+1);
			vecConnectVar[n] =
				ncOut.add_var(
					szConnectVarName, ncInt,
					vecElementBlockDim[n],
					vecNodesPerElementDim[n]);

			if (vecConnectVar[n] == NULL) {
				_EXCEPTION1("Error creating variable \"%s\"",
					szConnectVarName);
			}

			char szConnectAttrib[ParamLenString];
			snprintf(szConnectAttrib, ParamLenString, "SHELL%i", vecBlockSizes[n]);
			vecConnectVar[n]->add_att("elem_type", szConnectAttrib);

			char szGlobalIdVarName[ParamLenString];
			snprintf(szGlobalIdVarName, ParamLenString, "global_id%i", n+1);
			vecGlobalIdVar[n] = 
				ncOut.add_var(
					szGlobalIdVarName, ncInt,
					vecElementBlockDim[n]);

			if (vecGlobalIdVar[n] == NULL) {
				_EXCEPTION1("Error creating variable \"%s\"",
					szGlobalIdVarName);
			}

			char szEdgeTypeVarName[ParamLenString];
			snprintf(szEdgeTypeVarName, ParamLenString, "edge_type%i", n+1);
			vecEdgeTypeVar[n] =
				ncOut.add_var(
					szEdgeTypeVarName, ncInt,
					vecElementBlockDim[n],
					vecNodesPerElementDim[n]);

			if (vecEdgeTypeVar[n] == NULL) {
				_EXCEPTION1("Error creating variable \"%s\"",
					szEdgeTypeVarName);
			}

			if (source_face_ix.size() != 0) {
				vecFaceParentA[n].resize(vecBlockSizeFaces[n], 0);

				char szParentAVarName[ParamLenString];
				snprintf(szParentAVarName, ParamLenString, "el_parent_a%i", n+1);
				vecFaceParentAVar[n] =
					ncOut.add_var(
						szParentAVarName, ncInt,
						vecElementBlockDim[n]);

				if (vecFaceParentAVar[n] == NULL) {
					_EXCEPTION1("Error creating variable \"%s\"",
						szParentAVarName);
				}
			}

			if (target_face_ix.size() != 0) {
				vecFaceParentB[n].resize(vecBlockSizeFaces[n], 0);

				char szParentBVarName[ParamLenString];
				snprintf(szParentBVarName, ParamLenString, "el_parent_b%i", n+1);
				vecFaceParentBVar[n] =
					ncOut.add_var(
						szParentBVarName, ncInt,
						vecElementBlockDim[n]);

				if (vecFaceParentBVar[n] == NULL) {
					_EXCEPTION1("Error creating variable \"%s\"",
						szParentBVarName);
				}
			}
		}

		// Rebuild global data structures in local block arrays
		for (int i = 0; i < nFaceCount; i++) {
			int iBlock = 0;
			for (; iBlock < vecBlockSizes.size(); iBlock++) {
				if (vecBlockSizes[iBlock] == faces[i].size()) {
					break;
				}
			}
			if (iBlock == vecBlockSizes.size()) {
				_EXCEPTIONT("Logic error");
			}

			int iLocal = vecConnectCount[iBlock];
			for (int k = 0; k < faces[i].size(); k++) {
				vecConnect[iBlock](iLocal,k) = faces[i][k] + 1;

				//vecEdgeType[iBlock][iLocal][k] =
				//	static_cast<int>(faces[i][k].type);
			}

			vecGlobalId[iBlock][iLocal] = i + 1;

			if (source_face_ix.size() != 0) {
				vecFaceParentA[iBlock][iLocal] = source_face_ix[i] + 1;
			}
			if (target_face_ix.size() != 0) {
				vecFaceParentB[iBlock][iLocal] = target_face_ix[i] + 1;
			}

			vecConnectCount[iBlock]++;
		}

		// Write data to NetCDF file
		for (int n = 0; n < vecBlockSizes.size(); n++) {
			vecConnectVar[n]->set_cur(0, 0);
			vecConnectVar[n]->put(
				vecConnect[n],
				vecBlockSizeFaces[n],
				vecBlockSizes[n]);

			vecGlobalIdVar[n]->set_cur((long)0);
			vecGlobalIdVar[n]->put(
				&(vecGlobalId[n][0]),
				vecBlockSizeFaces[n]);

			//vecEdgeTypeVar[n]->set_cur(0, 0);
			//vecEdgeTypeVar[n]->put(
			//	&(vecEdgeType[n][0][0]),
			//	vecBlockSizeFaces[n],
			//	vecBlockSizes[n]);

			if (source_face_ix.size() != 0) {
				vecFaceParentAVar[n]->set_cur((long)0);
				vecFaceParentAVar[n]->put(
					&(vecFaceParentA[n][0]),
					vecBlockSizeFaces[n]);
			}

			if (target_face_ix.size() != 0) {
				vecFaceParentBVar[n]->set_cur((long)0);
				vecFaceParentBVar[n]->put(
					&(vecFaceParentB[n][0]),
					vecBlockSizeFaces[n]);
			}
		}
	}

	// Node list
	{
		NcVar * varNodes =
			ncOut.add_var("coord", ncDouble, dimDimension, dimNodes);

		if (varNodes == NULL) {
			_EXCEPTIONT("Error creating variable \"coord\"");
		}

		std::vector<double> dCoord(nNodeCount);

		for (int i = 0; i < nNodeCount; i++) {
			dCoord[i] = static_cast<double>(nodes[i].x);
		}
		varNodes->set_cur(0, 0);
		varNodes->put(&(dCoord[0]), 1, nNodeCount);
		for (int i = 0; i < nNodeCount; i++) {
			dCoord[i] = static_cast<double>(nodes[i].y);
		}
		varNodes->set_cur(1, 0);
		varNodes->put(&(dCoord[0]), 1, nNodeCount);
		for (int i = 0; i < nNodeCount; i++) {
			dCoord[i] = static_cast<double>(nodes[i].z);
		}
		varNodes->set_cur(2, 0);
		varNodes->put(&(dCoord[0]), 1, nNodeCount);
	}

	// Grid dimensions
	if (grid_dim_sizes.size() == 2) {
		_ASSERT(grid_dim_names.size() == 2);

		ncOut.add_att("rectilinear", "true");
		ncOut.add_att("rectilinear_dim0_size", grid_dim_sizes[0]);
		ncOut.add_att("rectilinear_dim1_size", grid_dim_sizes[1]);
		ncOut.add_att("rectilinear_dim0_name", grid_dim_names[0].c_str());
		ncOut.add_att("rectilinear_dim1_name", grid_dim_names[1].c_str());
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::write_as_scrip(
	const std::string & strFile,
	NcFile::FileFormat eFileFormat
) const {
	const int ParamLenString = 33;

	// Temporarily change error reporting
	NcError error_temp(NcError::verbose_fatal);

	// Determine block sizes
	std::vector<int> vecBlockSizes;
	std::vector<int> vecBlockSizeFaces;
	{
		std::map<int, int> mapBlockSizes;
		std::map<int, int>::iterator iterBlockSize;
		int iBlock;
		char szBuffer[ParamLenString];

		for (int i = 0; i < faces.size(); i++) {
			iterBlockSize = mapBlockSizes.find(faces[i].size());

			if (iterBlockSize == mapBlockSizes.end()) {
				mapBlockSizes.insert(
					std::pair<int,int>(faces[i].size(), 1));
			} else {
				(iterBlockSize->second)++;
			}
		}

		vecBlockSizes.resize(mapBlockSizes.size());
		vecBlockSizeFaces.resize(mapBlockSizes.size());

		AnnounceStartBlock("Nodes per element");
		iterBlockSize = mapBlockSizes.begin();
		iBlock = 1;
		for (; iterBlockSize != mapBlockSizes.end(); iterBlockSize++) {
			vecBlockSizes[iBlock-1] = iterBlockSize->first;
			vecBlockSizeFaces[iBlock-1] = iterBlockSize->second;

			Announce("Block %i (%i nodes): %i",
				iBlock, vecBlockSizes[iBlock-1], vecBlockSizeFaces[iBlock-1]);

			iBlock++;
		}
		AnnounceEndBlock(NULL);
	}

	// Output to a NetCDF SCRIP file
	NcFile ncOut(strFile.c_str(), NcFile::Replace, NULL, 0, eFileFormat);
	if (!ncOut.is_valid()) {
		_EXCEPTION1("Unable to open grid file \"%s\" for writing",
			strFile.c_str());
	}

	// Find max number of corners oer all faces
	int nFaceCount = faces.size();
	int nCornersMax = 0;
	for (int i = 0; i < nFaceCount; i++) {
		nCornersMax = std::max( nCornersMax, (int)(faces[i].size()) );
	}

	// SCRIP dimensions
	NcDim * dimGridSize   = ncOut.add_dim("grid_size", nFaceCount);
	NcDim * dimGridCorner = ncOut.add_dim("grid_corners", nCornersMax);

	NcDim * dimGridRank;
	if (grid_dim_sizes.size() <= 1) {
		dimGridRank = ncOut.add_dim("grid_rank", 1);
	} else {
		dimGridRank = ncOut.add_dim("grid_rank", grid_dim_sizes.size());
	}

	// Global attributes
	ncOut.add_att("api_version", 5.00f);
	ncOut.add_att("version", 5.00f);
	ncOut.add_att("floating_point_word_size", 8);
	ncOut.add_att("file_size", 0);

	// Grid Area
	{
		NcVar * varArea = ncOut.add_var("grid_area", ncDouble, dimGridSize);
		if (varArea == NULL) {
			_EXCEPTIONT("Error creating variable \"grid_area\"");
		}
		std::vector<double> area(nFaceCount);
		for (int i = 0; i < nFaceCount; i++) {
			area[i] = static_cast<double>( CalculateFaceArea(faces[i], nodes) );
		}
		varArea->set_cur((long)0);
		varArea->put(&(area[0]), nFaceCount);
		varArea->add_att("units", "radians^2");
	}

	// Grid center and corner coordinates
	{
		NcVar * varCenterLat = ncOut.add_var("grid_center_lat", ncDouble, dimGridSize);
		NcVar * varCenterLon = ncOut.add_var("grid_center_lon", ncDouble, dimGridSize);
		NcVar * varCornerLat = ncOut.add_var("grid_corner_lat", ncDouble, dimGridSize, dimGridCorner);
		NcVar * varCornerLon = ncOut.add_var("grid_corner_lon", ncDouble, dimGridSize, dimGridCorner);
		if (varCenterLat == NULL) {
			_EXCEPTIONT("Error creating variable \"grid_center_lat\"");
		}
		if (varCenterLon == NULL) {
			_EXCEPTIONT("Error creating variable \"grid_center_lon\"");
		}
		if (varCornerLat == NULL) {
			_EXCEPTIONT("Error creating variable \"grid_corner_lat\"");
		}
		if (varCornerLon == NULL) {
			_EXCEPTIONT("Error creating variable \"grid_corner_lon\"");
		}

		std::vector<double> centerLat(nFaceCount);
		std::vector<double> centerLon(nFaceCount);
		tmp::array2d<double> cornerLat(nFaceCount, nCornersMax, 0.0);
		tmp::array2d<double> cornerLon(nFaceCount, nCornersMax, 0.0);
		Node corner(0,0,0);
		for (int i = 0; i < nFaceCount; i++) {
			Node center(0,0,0);

			// int nCorners = faces[i].size()+1;
			int nCorners = faces[i].size();
			for (int j = 0; j < nCorners; j++) {
				corner = nodes[ faces[i][j] ];
				XYZtoRLL_Deg(
					corner.x, corner.y, corner.z,
					cornerLon(i,j),
					cornerLat(i,j));
				center = center + corner;
			}

			center = center / nCorners;

			center.normalize_in_place();

			XYZtoRLL_Deg(
				center.x, center.y, center.z,
				centerLon[i],
				centerLat[i]);

			// Adjust corner logitudes
			double lonDiff;
			for (int j = 0; j < nCorners; j++) {

				// First check for polar point
				if ((cornerLat(i,j)==90.) || (cornerLat(i,j)==(-90.))) {
					cornerLon(i,j) = centerLon[i];
				}

				// Next check for corners that wrap around prime meridian
				lonDiff = centerLon[i] - cornerLon(i,j);
				if (lonDiff > 180.0) {
					cornerLon(i,j) = cornerLon(i,j) + (double)360.0;
				}

				if (lonDiff < -180.0) {
					cornerLon(i,j) = cornerLon(i,j) - (double)360.0;
				}
			}
			// Make sure the padded coordinate data is same as last vertex
			for (int j = nCorners; j < nCornersMax; j++) {
				cornerLon(i,j) = cornerLon(i,nCorners-1);
				cornerLat(i,j) = cornerLat(i,nCorners-1);
			}
		}

		varCenterLat->add_att("_FillValue", 9.96920996838687e+36 );
		varCenterLon->add_att("_FillValue", 9.96920996838687e+36 );
		varCornerLat->add_att("_FillValue", 9.96920996838687e+36 );
		varCornerLon->add_att("_FillValue", 9.96920996838687e+36 );

		varCenterLat->set_cur((long)0);
		varCenterLat->put(&(centerLat[0]), nFaceCount);
		varCenterLat->add_att("units", "degrees");

		varCenterLon->set_cur((long)0);
		varCenterLon->put(&(centerLon[0]), nFaceCount);
		varCenterLon->add_att("units", "degrees");

		for (int i=0; i<nFaceCount; i++) {
			varCornerLat->set_cur(i,0);
			varCornerLat->put(&(cornerLat(i,0)), 1, nCornersMax);
			varCornerLon->set_cur(i,0);
			varCornerLon->put(&(cornerLon(i,0)), 1, nCornersMax);
		}
		varCornerLat->add_att("units", "degrees");
		varCornerLon->add_att("units", "degrees");
	}

	// Grid mask (set to 1)
	{
		NcVar * varMask = ncOut.add_var("grid_imask", ncInt, dimGridSize);
		if (varMask == NULL) {
			_EXCEPTIONT("Error creating variable \"grid_imask\"");
		}
		std::vector<int> mask(nFaceCount, 1);
		varMask->set_cur((long)0);
		varMask->put(&(mask[0]), nFaceCount);
	}

	// Grid dims
	{
		NcVar * varDims = ncOut.add_var("grid_dims", ncInt, dimGridRank);
		if (varDims == NULL) {
			_EXCEPTIONT("Error creating variable \"grid_dims\"");
		}
		if (grid_dim_sizes.size() <= 1) {
			int nSize = faces.size();
			varDims->set_cur((long)0);
			varDims->put(&nSize, 1);
		} else {
			varDims->set_cur((long)0);
			int nGridDimSize[2];
			nGridDimSize[0] = grid_dim_sizes[1];
			nGridDimSize[1] = grid_dim_sizes[0];
			varDims->put(&(nGridDimSize[0]), (long)grid_dim_sizes.size());
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::write_as_ugrid(
	const std::string & strFile,
	NcFile::FileFormat eFileFormat
) const {
	const int ParamLenString = 33;

	// Temporarily change error reporting
	NcError error_temp(NcError::verbose_fatal);

	// Output to a NetCDF SCRIP file
	NcFile ncOut(strFile.c_str(), NcFile::Replace, NULL, 0, eFileFormat);
	if (!ncOut.is_valid()) {
		_EXCEPTION1("Unable to open grid file \"%s\" for writing",
			strFile.c_str());
	}

	// Find max number of corners oer all faces
	int nCornersMax = 0;
	for (int i = 0; i < faces.size(); i++) {
		nCornersMax = std::max( nCornersMax, (int)(faces[i].size()) );
	}

	// Generate face nodes array
	tmp::array2d<int> nFaceNodes(faces.size(), nCornersMax);
	for (int i = 0; i < faces.size(); i++) {
		for (int j = 0; j < faces[i].size(); j++) {
			nFaceNodes(i,j) = faces[i][j];
		}
		for (int j = faces[i].size(); j < nCornersMax; j++) {
			nFaceNodes(i,j) = -1;
		}
	}

	// Number of nodes
	NcDim * dimNodes = ncOut.add_dim("nMesh2_node", nodes.size());
	NcDim * dimFaces = ncOut.add_dim("nMesh2_face", faces.size());
	NcDim * dimMaxNodesPerFace = ncOut.add_dim("nMaxMesh2_face_nodes", nCornersMax);

	// Mesh topology
	NcVar * varMesh2 = ncOut.add_var("Mesh2", ncInt);
	varMesh2->add_att("cf_role", "mesh_topology");
	varMesh2->add_att("long_name", "Topology data of 2D unstructured mesh");
	varMesh2->add_att("topology_dimension", 2);
	varMesh2->add_att("node_coordinates", "Mesh2_node_x Mesh2_node_y");
	varMesh2->add_att("node_dimension", "nMesh2_node");
	varMesh2->add_att("face_node_connectivity", "Mesh2_face_nodes");
	varMesh2->add_att("face_dimension", "nMesh2_face");

	// Face nodes
	NcVar * varFaceNodes = ncOut.add_var("Mesh2_face_nodes", ncInt, dimFaces, dimMaxNodesPerFace);
	varFaceNodes->add_att("cf_role", "face_node_connectivity");
	varFaceNodes->add_att("_FillValue", -1);
	varFaceNodes->add_att("start_index", 0);
	varFaceNodes->put(&(nFaceNodes(0,0)), faces.size(), nCornersMax);

	// Mesh node coordinates
	std::vector<double> dNodeLat(nodes.size());
	std::vector<double> dNodeLon(nodes.size());

	for (int i = 0; i < nodes.size(); i++) {
		LatLon nodelatlon = nodes[i].to_latlon_deg();
		dNodeLat[i] = nodelatlon.lat;
		dNodeLon[i] = nodelatlon.lon;
	}

	NcVar * varNodeX = ncOut.add_var("Mesh2_node_x", ncDouble, dimNodes);
	varNodeX->add_att("standard_name", "longitude");
	varNodeX->add_att("long_name", "longitude of 2D mesh nodes");
	varNodeX->add_att("units", "degrees_east");
	varNodeX->put(&(dNodeLon[0]), nodes.size());

	NcVar * varNodeY = ncOut.add_var("Mesh2_node_y", ncDouble, dimNodes);
	varNodeY->add_att("standard_name", "latitude");
	varNodeY->add_att("long_name", "latitude of 2D mesh nodes");
	varNodeY->add_att("units", "degrees_north");
	varNodeY->put(&(dNodeLat[0]), nodes.size());

	// Mask
	if (face_mask.size() == faces.size()) {
		varMesh2->add_att("face_mask", "Mesh2_face_mask");

		NcVar * varIMask = ncOut.add_var("Mesh2_face_mask", ncInt, dimFaces);
		varIMask->add_att("standard_name", "mask");
		varIMask->add_att("long_name", "integer mask of faces");
		varIMask->add_att("units", "none");
		varIMask->put(&(face_mask[0]), faces.size());
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::read(const std::string & strFile) {

	const int ParamFour = 4;
	const int ParamLenString = 33;

	// Store the file name
	filename = strFile;

	// Open the NetCDF file
	if (strFile == "") {
		_EXCEPTIONT("No grid file specified for reading");
	}
	NcFile ncFile(strFile.c_str(), NcFile::ReadOnly);
	if (!ncFile.is_valid()) {
		_EXCEPTION1("Unable to open grid file \"%s\" for reading",
			strFile.c_str());
	}

	// Check for dimension names "grid_size", "grid_rank" and "grid_corners"
	int iSCRIPFormat = 0;
	for (int i = 0; i < ncFile.num_dims(); i++) {
		NcDim * dim = ncFile.get_dim(i);
		std::string strDimName = dim->name();
		if (strDimName == "grid_size") {
			iSCRIPFormat++;
		}
		if (strDimName == "grid_corners") {
			iSCRIPFormat++;
		}
		if (strDimName == "grid_rank") {
			iSCRIPFormat++;
		}
	}

	// Input from a NetCDF SCRIP file
	if (iSCRIPFormat == 3) {
		Announce("SCRIP Format File detected");

		NcDim * dimGridSize = ncFile.get_dim("grid_size");
		NcDim * dimGridCorners = ncFile.get_dim("grid_corners");

		// Get the grid corners
		NcVar * varGridCornerLat = ncFile.get_var("grid_corner_lat");
		NcVar * varGridCornerLon = ncFile.get_var("grid_corner_lon");

		if (varGridCornerLat == NULL) {
			_EXCEPTION1("SCRIP Grid file \"%s\" is missing variable "
					"\"grid_corner_lat\"", strFile.c_str());
		}
		if (varGridCornerLon == NULL) {
			_EXCEPTION1("SCRIP Grid file \"%s\" is missing variable "
					"\"grid_corner_lon\"", strFile.c_str());
		}

		int nGridSize = static_cast<int>(dimGridSize->size());
		int nGridCorners = static_cast<int>(dimGridCorners->size());

		tmp::array2d<double> dCornerLat(nGridSize, nGridCorners);
		tmp::array2d<double> dCornerLon(nGridSize, nGridCorners);

		varGridCornerLat->set_cur(0, 0);
		varGridCornerLat->get(dCornerLat, nGridSize, nGridCorners);

		varGridCornerLon->set_cur(0, 0);
		varGridCornerLon->get(dCornerLon, nGridSize, nGridCorners);

		faces.resize(nGridSize);
		nodes.resize(nGridSize * nGridCorners);

		// Check for units attribute; if "degrees" then convert to radians
		bool fConvertLonToRadians = false;
		NcAtt * attGridCornerLonUnits = varGridCornerLon->get_att("units");

		if (attGridCornerLonUnits != NULL) {
			std::string strLonUnits = attGridCornerLonUnits->as_string(0);
			STLStringHelper::ToLower(strLonUnits);
			if (strLonUnits == "degrees") {
				fConvertLonToRadians = true;
			}
		}

		bool fConvertLatToRadians = false;
		NcAtt * attGridCornerLatUnits = varGridCornerLat->get_att("units");

		if (attGridCornerLatUnits != NULL) {
			std::string strLatUnits = attGridCornerLatUnits->as_string(0);
			STLStringHelper::ToLower(strLatUnits);
			if (strLatUnits == "degrees") {
				fConvertLatToRadians = true;
			}
		}

		// Load mask variable
		NcVar * varMask = ncFile.get_var("grid_imask");
		if (varMask != NULL) {
			if (varMask->num_dims() != 1) {
				_EXCEPTIONT("Unknown format of variable \"grid_imask\": "
					"More than one dimension");
			}
			if (varMask->get_dim(0)->size() != nGridSize) {
				_EXCEPTIONT("Unknown format of variable \"grid_imask\": "
					"Incorrect first dimension size");
			}

			face_mask.resize(nGridSize);
			varMask->get(&(face_mask[0]), nGridSize);
		}

		// Current global node index
		int ixNode = 0;

		for (int i = 0; i < nGridSize; i++) {

			// Create a new Face
			Face faceNew(nGridCorners);
			for (int j = 0; j < nGridCorners; j++) {
				faceNew[j] = ixNode + j;
			}
			faces[i] = faceNew;

			// Insert Face corners into node table
			for (int j = 0; j < nGridCorners; j++) {
				LatLon latlon;
				latlon.lon = dCornerLon(i,j);
				latlon.lat = dCornerLat(i,j);

				if (fConvertLonToRadians) {
					latlon.lon = DegToRad(latlon.lon);
				}
				if (fConvertLatToRadians) {
					latlon.lat = DegToRad(latlon.lat);
				}

				nodes[ixNode].from_latlon_rad(latlon);

				ixNode++;
			}
		}

		// grid_rank
		NcDim * dimGridRank = ncFile.get_dim("grid_rank");
		if (dimGridRank != NULL) {
			if (dimGridRank->size() == 2) {
				NcVar * varGridDims = ncFile.get_var("grid_dims");
				if (varGridDims == NULL) {
					_EXCEPTION1("SCRIP grid file \"%s\" has grid_rank 2 but does not contain variable grid_dims",
						strFile.c_str());
				}
				if (varGridDims->num_dims() != 1) {
					_EXCEPTION1("SCRIP grid file \"%s\" variable grid_dims has more than one dimension",
						strFile.c_str());
				}
				if (varGridDims->get_dim(0)->size() != dimGridRank->size()) {
					_EXCEPTION1("SCRIP grid file \"%s\" variable grid_dims size does not match grid_rank",
						strFile.c_str());
				}

				grid_dim_sizes.resize(2);
				grid_dim_names.resize(2);

				int nGridDims[2];
				varGridDims->get(&(nGridDims[0]), 2);
				grid_dim_sizes[0] = nGridDims[1];
				grid_dim_sizes[1] = nGridDims[0];
				grid_dim_names[0] = "griddim1";
				grid_dim_names[1] = "griddim0";

				if (grid_dim_sizes[0] * grid_dim_sizes[1] != faces.size()) {
					_EXCEPTION4("SCRIP grid file \"%s\" grid_dims (%i,%i) do not agree with grid_size (%i)",
						strFile.c_str(), grid_dim_sizes[0], grid_dim_sizes[1], faces.size());
				}
			}
		}

		// SCRIP does not reference a node table, so we must remove
		// coincident nodes.
		remove_coincident_nodes();

		// Output size
		Announce("Mesh size: Nodes [%i] Elements [%i]",
			nodes.size(), faces.size());

	// Input from a NetCDF Exodus file
	} else {

		// Get version number
		NcAtt * attVersion = ncFile.get_att("version");
		if (attVersion == NULL) {
			_EXCEPTION1("Exodus Grid file \"%s\" is missing attribute "
					"\"version\"", strFile.c_str());
		}
		if (attVersion->type() != ncFloat) {
			_EXCEPTIONT("Exodus Grid type is not of type float");
		}
		float flVersion = attVersion->as_float(0);

		// Number of nodes
		NcDim * dimNodes = ncFile.get_dim("num_nodes");
		if (dimNodes == NULL) {
			_EXCEPTION1("Exodus Grid file \"%s\" is missing dimension "
					"\"num_nodes\"", strFile.c_str());
		}
		int nNodeCount = dimNodes->size();

		// Determine number of blocks
		NcDim * dimElementBlocks = ncFile.get_dim("num_el_blk");
		if (dimElementBlocks == NULL) {
			_EXCEPTION1("Exodus Grid file \"%s\" is missing dimension "
					"\"num_el_blk\"", strFile.c_str());
		}
		int nElementBlocks = dimElementBlocks->size();

		// Total number of elements
		int nTotalElementCount = 0;
		NcDim * dimElements = ncFile.get_dim("num_elem");
		if (dimElements == NULL) {
			for (int ib = 1; ib <= nElementBlocks; ib++) {
				std::string numelblk = std::string("num_el_in_blk" + std::to_string(ib));
				NcDim * dimElementBlockElems = ncFile.get_dim(numelblk.c_str());
				if (dimElementBlockElems != NULL) {
					nTotalElementCount += dimElementBlockElems->size();
				}
			}
			if (nTotalElementCount == 0) {
				_EXCEPTION1("Exodus Grid file \"%s\" is missing dimension "
					"\"num_elem\"", strFile.c_str());
			}

		} else {
			nTotalElementCount = dimElements->size();
		}

		// Output size
		Announce("Mesh size: Nodes [%i] Elements [%i]",
			nNodeCount, nTotalElementCount);

		// Allocate faces
		faces.resize(nTotalElementCount);

		// Loop over all blocks
		for (int n = 0; n < nElementBlocks; n++) {

			// Determine number of nodes per element in this block
			char szNodesPerElement[ParamLenString];
			snprintf(szNodesPerElement, ParamLenString, "num_nod_per_el%i", n+1);
			NcDim * dimNodesPerElement = ncFile.get_dim(szNodesPerElement);
			if (dimNodesPerElement == NULL) {
				_EXCEPTION2("Exodus Grid file \"%s\" is missing dimension "
					"\"%s\"", strFile.c_str(), szNodesPerElement);
			}
			int nNodesPerElement = dimNodesPerElement->size();

			// Number of elements in block
			char szElementsInBlock[ParamLenString];
			snprintf(szElementsInBlock, ParamLenString, "num_el_in_blk%i", n+1);

			NcDim * dimBlockElements = ncFile.get_dim(szElementsInBlock);
			if (dimBlockElements == NULL) {
				_EXCEPTION2("Exodus Grid file \"%s\" is missing dimension "
						"\"%s\"", strFile.c_str(), szElementsInBlock);
			}
			int nFaceCount = dimBlockElements->size();

			// Variables for each face
			tmp::array2d<int> iConnect(nFaceCount, nNodesPerElement, 0);
			std::vector<int> iGlobalId(nFaceCount, 0);

			std::vector<int> iParentA(nFaceCount, InvalidFace);
			std::vector<int> iParentB(nFaceCount, InvalidFace);

			// Load in nodes for all elements in this block
			char szConnect[ParamLenString];
			snprintf(szConnect, ParamLenString, "connect%i", n+1);

			NcVar * varConnect = ncFile.get_var(szConnect);
			if (varConnect == NULL) {
				_EXCEPTION2("Exodus Grid file \"%s\" is missing variable "
						"\"%s\"", strFile.c_str(), szConnect);
			}

			varConnect->set_cur(0, 0);
			varConnect->get(
				iConnect,
				nFaceCount,
				nNodesPerElement);

			// Earlier version didn't have global_id
			if (flVersion == 4.98f) {
				for (int i = 0; i < nFaceCount; i++) {
					iGlobalId[i] = i + 1;
				}

			// Load in global id for all elements in this block
			} else {
				char szGlobalId[ParamLenString];
				snprintf(szGlobalId, ParamLenString, "global_id%i", n+1);

				NcVar * varGlobalId = ncFile.get_var(szGlobalId);
				if (varGlobalId == NULL) {
					_EXCEPTION2("Exodus Grid file \"%s\" is missing variable "
							"\"%s\"", strFile.c_str(), szGlobalId);
				}

				varGlobalId->set_cur((long)0);
				varGlobalId->get(
					&(iGlobalId[0]),
					nFaceCount);
			}

			// Load in edge type for all elements in this block
			char szEdgeType[ParamLenString];
			if (flVersion == 4.98f) {
				snprintf(szEdgeType, ParamLenString, "edge_type");
			} else {
				snprintf(szEdgeType, ParamLenString, "edge_type%i", n+1);
			}

			NcVar * varEdgeType = ncFile.get_var(szEdgeType);
			if (varEdgeType != NULL) {
				//varEdgeType->set_cur(0, 0);
				//varEdgeType->get(
				//	&(iEdgeType[0][0]),
				//	nFaceCount,
				//	nNodesPerElement);
			}

			// Load in parent from A grid for all elements in this block
			char szParentA[ParamLenString];
			if (flVersion == 4.98f) {
				snprintf(szParentA, ParamLenString, "face_source_1");
			} else {
				snprintf(szParentA, ParamLenString, "el_parent_a%i", n+1);
			}

			NcVar * varParentA = ncFile.get_var(szParentA);
			if ((varParentA == NULL) && (source_face_ix.size() != 0)) {
				_EXCEPTION2("Exodus Grid file \"%s\" is missing variable "
						"\"%s\"", strFile.c_str(), szParentA);

			} else if (varParentA != NULL) {
				if (source_face_ix.size() == 0) {
					source_face_ix.resize(nTotalElementCount);
				}

				varParentA->set_cur((long)0);
				varParentA->get(
					&(iParentA[0]),
					nFaceCount);
			}

			// Load in parent from A grid for all elements in this block
			char szParentB[ParamLenString];
			if (flVersion == 4.98f) {
				snprintf(szParentB, ParamLenString, "face_source_2");
			} else {
				snprintf(szParentB, ParamLenString, "el_parent_b%i", n+1);
			}

			NcVar * varParentB = ncFile.get_var(szParentB);
			if ((varParentB == NULL) && (target_face_ix.size() != 0)) {
				_EXCEPTION2("Exodus Grid file \"%s\" is missing variable "
						"\"%s\"", strFile.c_str(), szParentB);

			} else if (varParentB != NULL) {
				if (target_face_ix.size() == 0) {
					target_face_ix.resize(nTotalElementCount);
				}

				varParentB->set_cur((long)0);
				varParentB->get(
					&(iParentB[0]),
					nFaceCount);
			}

			// Put local data into global structures
			for (int i = 0; i < nFaceCount; i++) {
				if (iGlobalId[i] - 1 >= nTotalElementCount) {
					_EXCEPTION2("global_id %i out of range [1,%i]",
						iGlobalId[i], nTotalElementCount);
				}
				faces[iGlobalId[i]-1] = Face(nNodesPerElement);
				for (int k = 0; k < nNodesPerElement; k++) {
					faces[iGlobalId[i] - 1][k] = iConnect(i,k) - 1;
				}

				if (source_face_ix.size() != 0) {
					source_face_ix[iGlobalId[i] - 1] = iParentA[i] - 1;
				}

				if (target_face_ix.size() != 0) {
					target_face_ix[iGlobalId[i] - 1] = iParentB[i] - 1;
				}
			}
		}

		// Earlier version had incorrect parent indexing
		if (flVersion == 4.98f) {
			if (source_face_ix.size() != 0) {
				for (int i = 0; i < nTotalElementCount; i++) {
					source_face_ix[i]++;
				}
			}
			if (target_face_ix.size() != 0) {
				for (int i = 0; i < nTotalElementCount; i++) {
					target_face_ix[i]++;
				}
			}
		}

		// Load in node array
		{
			nodes.resize(nNodeCount);

			NcVar * varNodes = ncFile.get_var("coord");
			if (varNodes == NULL) {
				_EXCEPTION1("Exodus Grid file \"%s\" is missing variable "
						"\"coord\"", strFile.c_str());
			}

			tmp::array2d<double> dNodeCoords(3, nNodeCount);

			// Load in node array
			varNodes->set_cur(0, 0);
			varNodes->get(dNodeCoords, 3, nNodeCount);

			for (int i = 0; i < nNodeCount; i++) {
				nodes[i].x = static_cast<Real>(dNodeCoords(0,i));
				nodes[i].y = static_cast<Real>(dNodeCoords(1,i));
				nodes[i].z = static_cast<Real>(dNodeCoords(2,i));
			}
		}

		// Load in dimension information
		NcAtt * attRectilinear = ncFile.get_att("rectilinear");
		if (attRectilinear != NULL) {
			NcAtt * attDim0Size = ncFile.get_att("rectilinear_dim0_size");
			if (attDim0Size == NULL) {
				_EXCEPTION1("Exodus Grid file \"%s\" has rectilinear attribute set, "
					"but is missing attribute \"rectilinear_dim0_size\"", strFile.c_str());
			}

			NcAtt * attDim1Size = ncFile.get_att("rectilinear_dim1_size");
			if (attDim1Size == NULL) {
				_EXCEPTION1("Exodus Grid file \"%s\" has rectilinear attribute set, "
					"but is missing attribute \"rectilinear_dim1_size\"", strFile.c_str());
			}

			NcAtt * attDim0Name = ncFile.get_att("rectilinear_dim0_name");
			if (attDim0Size == NULL) {
				_EXCEPTION1("Exodus Grid file \"%s\" has rectilinear attribute set, "
					"but is missing attribute \"rectilinear_dim0_name\"", strFile.c_str());
			}

			NcAtt * attDim1Name = ncFile.get_att("rectilinear_dim1_name");
			if (attDim0Size == NULL) {
				_EXCEPTION1("Exodus grid file \"%s\" has rectilinear attribute set, "
					"but is missing attribute \"rectilinear_dim1_name\"", strFile.c_str());
			}

			grid_dim_sizes.resize(2);
			grid_dim_names.resize(2);

			grid_dim_sizes[0] = attDim0Size->as_int(0);
			grid_dim_sizes[1] = attDim1Size->as_int(0);
			grid_dim_names[0] = attDim0Name->as_string(0);
			grid_dim_names[1] = attDim1Name->as_string(0);

			if (grid_dim_sizes[0] * grid_dim_sizes[1] != faces.size()) {
				_EXCEPTION4("Exodus grid file \"%s\" grid_dims (%i,%i) do not agree with grid_size (%i)",
					strFile.c_str(), grid_dim_sizes[0], grid_dim_sizes[1], faces.size());
			}
		}

		// Remove coincident nodes.
		remove_coincident_nodes();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::remove_zero_edges() {
	for (Face & face : faces) {
		face.remove_repeated_nodes();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::validate() const {

	// Valid that Nodes have magnitude 1
	for (size_t i = 0; i < nodes.size(); i++) {
		const Node & node = nodes[i];
		double dMag = node.magnitude();

		if (fabs(dMag - 1.0) > ReferenceTolerance) {
			_EXCEPTION5("Mesh validation failed: "
				"Node[%i] of non-unit magnitude detected (%1.10e, %1.10e, %1.10e) = %1.10e",
				i, node.x, node.y, node.z, dMag);
		}
	}

	// Validate that edges are oriented counter-clockwise
	for (int i = 0; i < faces.size(); i++) {
		const Face & face = faces[i];

		const int nEdges = face.size();
/*
		for (int j = 0; j < nEdges; j++) {

			// Check for zero edges
			for(;;) {
				if (face.edges[j][0] == face.edges[j][1]) {
					j++;
				} else {
					break;
				}
				if (j == nEdges) {
					break;
				}
			}

			if (j == nEdges) {
				break;
			}

			// Find the next non-zero edge
			int jNext = (j + 1) % nEdges;

			for(;;) {
				if (face.edges[jNext][0] == face.edges[jNext][1]) {
					jNext++;
				} else {
					break;
				}
				if (jNext == nEdges) {
					jNext = 0;
				}
				if (jNext == ((j + 1) % nEdges)) {
					_EXCEPTIONT("Mesh validation failed: "
						"No edge information on Face");
				}
			}

			// Get edges
			const Edge & edge0 = face.edges[j];
			const Edge & edge1 = face.edges[(j + 1) % nEdges];

			if (edge0[1] != edge1[0]) {
				_EXCEPTIONT("Mesh validation failed: Edge cyclicity error");
			}

			const Node & node0 = nodes[edge0[0]];
			const Node & node1 = nodes[edge0[1]];
			const Node & node2 = nodes[edge1[1]];

			// Vectors along edges
			Node nodeD1 = node0 - node1;
			Node nodeD2 = node2 - node1;

			// Compute cross-product
			Node nodeCross(CrossProduct(nodeD1, nodeD2));

			// Dot cross product with radial vector
			Real dDot = DotProduct(node1, nodeCross);
*/
/*
			if (dDot > 0.0) {
				printf("\nError detected (orientation):\n");
				printf("  Face %i, Edge %i, Orientation %1.5e\n",
					i, j, dDot);

				printf("  (x,y,z):\n");
				printf("	n0: %1.5e %1.5e %1.5e\n", node0.x, node0.y, node0.z);
				printf("	n1: %1.5e %1.5e %1.5e\n", node1.x, node1.y, node1.z);
				printf("	n2: %1.5e %1.5e %1.5e\n", node2.x, node2.y, node2.z);

				Real dR0 = sqrt(
					node0.x * node0.x + node0.y * node0.y + node0.z * node0.z);
				Real dLat0 = asin(node0.z / dR0);
				Real dLon0 = atan2(node0.y, node0.x);

				Real dR1 = sqrt(
					node1.x * node1.x + node1.y * node1.y + node1.z * node1.z);
				Real dLat1 = asin(node1.z / dR1);
				Real dLon1 = atan2(node1.y, node1.x);

				Real dR2 = sqrt(
					node2.x * node2.x + node2.y * node2.y + node2.z * node2.z);
				Real dLat2 = asin(node2.z / dR2);
				Real dLon2 = atan2(node2.y, node2.x);

				printf("  (lambda, phi):\n");
				printf("	n0: %1.5e %1.5e\n", dLon0, dLat0);
				printf("	n1: %1.5e %1.5e\n", dLon1, dLat1);
				printf("	n2: %1.5e %1.5e\n", dLon2, dLat2);

				printf("  X-Product:\n");
				printf("	%1.5e %1.5e %1.5e\n",
					nodeCross.x, nodeCross.y, nodeCross.z);

				_EXCEPTIONT(
					"Mesh validation failed: Clockwise or concave face detected");
			}
		}
*/
	}
}

///////////////////////////////////////////////////////////////////////////////
// General purpose functions
///////////////////////////////////////////////////////////////////////////////
/*
bool IsPositivelyOrientedEdge(
	const Node & nodeBegin,
	const Node & nodeEnd
) {
	const Real Tolerance = ReferenceTolerance;

	if ((fabs(nodeBegin.x - nodeEnd.x) < Tolerance) &&
		(fabs(nodeBegin.y - nodeEnd.y) < Tolerance) &&
		(fabs(nodeBegin.z - nodeEnd.z) < Tolerance)
	) {
		_EXCEPTIONT("Latitude line of zero length");
	}

	// Both nodes in positive y half-plane
	if ((nodeBegin.y >= 0.0) && (nodeEnd.y >= 0.0)) {
		if (nodeEnd.x < nodeBegin.x) {
			return true;
		} else {
			return false;
		}

	// Both nodes in negative y half-plane
	} else if ((nodeBegin.y <= 0.0) && (nodeEnd.y <= 0.0)) {
		if (nodeEnd.x > nodeBegin.x) {
			return true;
		} else {
			return false;
		}

	// Both nodes in positive x half-plane
	} else if ((nodeBegin.x >= 0.0) && (nodeEnd.x >= 0.0)) {
		if (nodeEnd.y > nodeBegin.y) {
			return true;
		} else {
			return false;
		}

	// Both nodes in negative x half-plane
	} else if ((nodeBegin.x <= 0.0) && (nodeEnd.x <= 0.0)) {
		if (nodeEnd.y < nodeBegin.y) {
			return true;
		} else {
			return false;
		}

	// Arc length too large
	} else {
		_EXCEPTIONT("Arc length too large to determine orientation.");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GetLocalDirection(
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Node & nodeRef,
	const Edge::Type edgetype,
	Node & nodeDir
) {

	// Direction along a great circle arc
	if (edgetype == Edge::Type_GreatCircleArc) {

		// Cartesian direction
		nodeDir = nodeEnd - nodeBegin;

		// Project onto surface of the sphere
		Real dDotDirBegin = DotProduct(nodeDir, nodeRef);
		Real dNormNodeBegin = DotProduct(nodeRef, nodeRef);

		nodeDir.x -= dDotDirBegin / dNormNodeBegin * nodeRef.x;
		nodeDir.y -= dDotDirBegin / dNormNodeBegin * nodeRef.y;
		nodeDir.z -= dDotDirBegin / dNormNodeBegin * nodeRef.z;

	// Direction along a line of constant latitude
	} else if (edgetype == Edge::Type_ConstantLatitude) {
		nodeDir.z = 0.0;

		if (IsPositivelyOrientedEdge(nodeBegin, nodeEnd)) {
			nodeDir.x = - nodeBegin.y;
			nodeDir.y = + nodeBegin.x;
		} else {
			nodeDir.x = + nodeBegin.y;
			nodeDir.y = - nodeBegin.x;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GetLocalDirection(
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	Node & nodeDir
) {

	// Direction along a great circle arc
	if (edgetype == Edge::Type_GreatCircleArc) {

		// Cartesian direction
		nodeDir = nodeEnd - nodeBegin;

		// Project onto surface of the sphere
		Real dDotDirBegin   = DotProduct(nodeDir, nodeBegin);
		Real dNormNodeBegin = DotProduct(nodeBegin, nodeBegin);

		nodeDir.x -= dDotDirBegin / dNormNodeBegin * nodeBegin.x;
		nodeDir.y -= dDotDirBegin / dNormNodeBegin * nodeBegin.y;
		nodeDir.z -= dDotDirBegin / dNormNodeBegin * nodeBegin.z;

	// Direction along a line of constant latitude
	} else if (edgetype == Edge::Type_ConstantLatitude) {
		nodeDir.z = 0.0;

		if (IsPositivelyOrientedEdge(nodeBegin, nodeEnd)) {
			nodeDir.x = - nodeBegin.y;
			nodeDir.y = + nodeBegin.x;
		} else {
			nodeDir.x = + nodeBegin.y;
			nodeDir.y = - nodeBegin.x;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

Real CalculateTriangleAreaQuadratureMethod(Node & node1, Node & node2,
		Node & node3) {

	const int nOrder = 6;

	DataArray1D<double> dG;
	DataArray1D<double> dW;
	GaussQuadrature::GetPoints(nOrder, 0.0, 1.0, dG, dW);

	double dArea = 0.0;

	// Calculate area at quadrature node
	for (int p = 0; p < dW.GetRows(); p++) {
		for (int q = 0; q < dW.GetRows(); q++) {

			double dA = dG[p];
			double dB = dG[q];

			Node dF(
					(1.0 - dB) * ((1.0 - dA) * node1.x + dA * node2.x)
							+ dB * node3.x,
					(1.0 - dB) * ((1.0 - dA) * node1.y + dA * node2.y)
							+ dB * node3.y,
					(1.0 - dB) * ((1.0 - dA) * node1.z + dA * node2.z)
							+ dB * node3.z);

			Node dDaF((1.0 - dB) * (node2.x - node1.x),
					(1.0 - dB) * (node2.y - node1.y),
					(1.0 - dB) * (node2.z - node1.z));

			Node dDbF(-(1.0 - dA) * node1.x - dA * node2.x + node3.x,
					-(1.0 - dA) * node1.y - dA * node2.y + node3.y,
					-(1.0 - dA) * node1.z - dA * node2.z + node3.z);

			double dR = sqrt(dF.x * dF.x + dF.y * dF.y + dF.z * dF.z);

			Node dDaG(
					dDaF.x * (dF.y * dF.y + dF.z * dF.z)
							- dF.x * (dDaF.y * dF.y + dDaF.z * dF.z),
					dDaF.y * (dF.x * dF.x + dF.z * dF.z)
							- dF.y * (dDaF.x * dF.x + dDaF.z * dF.z),
					dDaF.z * (dF.x * dF.x + dF.y * dF.y)
							- dF.z * (dDaF.x * dF.x + dDaF.y * dF.y));

			Node dDbG(
					dDbF.x * (dF.y * dF.y + dF.z * dF.z)
							- dF.x * (dDbF.y * dF.y + dDbF.z * dF.z),
					dDbF.y * (dF.x * dF.x + dF.z * dF.z)
							- dF.y * (dDbF.x * dF.x + dDbF.z * dF.z),
					dDbF.z * (dF.x * dF.x + dF.y * dF.y)
							- dF.z * (dDbF.x * dF.x + dDbF.y * dF.y));

			double dDenomTerm = 1.0 / (dR * dR * dR);

			dDaG.x *= dDenomTerm;
			dDaG.y *= dDenomTerm;
			dDaG.z *= dDenomTerm;

			dDbG.x *= dDenomTerm;
			dDbG.y *= dDenomTerm;
			dDbG.z *= dDenomTerm;

			// Cross product gives local Jacobian
			Node nodeCross = CrossProduct(dDaG, dDbG);

			double dJacobian = sqrt(
					nodeCross.x * nodeCross.x + nodeCross.y * nodeCross.y
							+ nodeCross.z * nodeCross.z);

			dArea += dW[p] * dW[q] * dJacobian;
		}
	}
	return dArea;
}

///////////////////////////////////////////////////////////////////////////////

Real CalculateSphericalTriangleJacobian(
	const Node & node1,
	const Node & node2,
	const Node & node3,
	double dA,
	double dB,
	Node * pnode
) {

	Node dF(
		(1.0 - dB) * ((1.0 - dA) * node1.x + dA * node2.x) + dB * node3.x,
		(1.0 - dB) * ((1.0 - dA) * node1.y + dA * node2.y) + dB * node3.y,
		(1.0 - dB) * ((1.0 - dA) * node1.z + dA * node2.z) + dB * node3.z);

	Node dDaF(
		(1.0 - dB) * (node2.x - node1.x),
		(1.0 - dB) * (node2.y - node1.y),
		(1.0 - dB) * (node2.z - node1.z));

	Node dDbF(
		- (1.0 - dA) * node1.x - dA * node2.x + node3.x,
		- (1.0 - dA) * node1.y - dA * node2.y + node3.y,
		- (1.0 - dA) * node1.z - dA * node2.z + node3.z);

	double dInvR = 1.0 / sqrt(dF.x * dF.x + dF.y * dF.y + dF.z * dF.z);

	if (pnode != NULL) {
		pnode->x = dF.x * dInvR;
		pnode->y = dF.y * dInvR;
		pnode->z = dF.z * dInvR;
	}

	Node dDaG(
		dDaF.x * (dF.y * dF.y + dF.z * dF.z)
			- dF.x * (dDaF.y * dF.y + dDaF.z * dF.z),
		dDaF.y * (dF.x * dF.x + dF.z * dF.z)
			- dF.y * (dDaF.x * dF.x + dDaF.z * dF.z),
		dDaF.z * (dF.x * dF.x + dF.y * dF.y)
			- dF.z * (dDaF.x * dF.x + dDaF.y * dF.y));

	Node dDbG(
		dDbF.x * (dF.y * dF.y + dF.z * dF.z)
			- dF.x * (dDbF.y * dF.y + dDbF.z * dF.z),
		dDbF.y * (dF.x * dF.x + dF.z * dF.z)
			- dF.y * (dDbF.x * dF.x + dDbF.z * dF.z),
		dDbF.z * (dF.x * dF.x + dF.y * dF.y)
			- dF.z * (dDbF.x * dF.x + dDbF.y * dF.y));

	double dDenomTerm = dInvR * dInvR * dInvR;

	dDaG.x *= dDenomTerm;
	dDaG.y *= dDenomTerm;
	dDaG.z *= dDenomTerm;

	dDbG.x *= dDenomTerm;
	dDbG.y *= dDenomTerm;
	dDbG.z *= dDenomTerm;

	// Cross product gives local Jacobian
	Node nodeCross = CrossProduct(dDaG, dDbG);

	double dJacobian = sqrt(
		  nodeCross.x * nodeCross.x
		+ nodeCross.y * nodeCross.y
		+ nodeCross.z * nodeCross.z);

	return dJacobian;
}

///////////////////////////////////////////////////////////////////////////////

Real CalculateSphericalTriangleJacobianBarycentric(
	const Node & node1,
	const Node & node2,
	const Node & node3,
	double dA,
	double dB,
	Node * pnode
) {
	double dC = 1.0 - dA - dB;

	Node dF(
		dA * node1.x + dB * node2.x + dC * node3.x,
		dA * node1.y + dB * node2.y + dC * node3.y,
		dA * node1.z + dB * node2.z + dC * node3.z);

	Node dDaF(
		node1.x - node3.x,
		node1.y - node3.y,
		node1.z - node3.z);

	Node dDbF(
		node2.x - node3.x,
		node2.y - node3.y,
		node2.z - node3.z);

	double dInvR = 1.0 / sqrt(dF.x * dF.x + dF.y * dF.y + dF.z * dF.z);

	if (pnode != NULL) {
		pnode->x = dF.x * dInvR;
		pnode->y = dF.y * dInvR;
		pnode->z = dF.z * dInvR;
	}

	Node dDaG(
		dDaF.x * (dF.y * dF.y + dF.z * dF.z)
			- dF.x * (dDaF.y * dF.y + dDaF.z * dF.z),
		dDaF.y * (dF.x * dF.x + dF.z * dF.z)
			- dF.y * (dDaF.x * dF.x + dDaF.z * dF.z),
		dDaF.z * (dF.x * dF.x + dF.y * dF.y)
			- dF.z * (dDaF.x * dF.x + dDaF.y * dF.y));

	Node dDbG(
		dDbF.x * (dF.y * dF.y + dF.z * dF.z)
			- dF.x * (dDbF.y * dF.y + dDbF.z * dF.z),
		dDbF.y * (dF.x * dF.x + dF.z * dF.z)
			- dF.y * (dDbF.x * dF.x + dDbF.z * dF.z),
		dDbF.z * (dF.x * dF.x + dF.y * dF.y)
			- dF.z * (dDbF.x * dF.x + dDbF.y * dF.y));

	double dDenomTerm = dInvR * dInvR * dInvR;

	dDaG.x *= dDenomTerm;
	dDaG.y *= dDenomTerm;
	dDaG.z *= dDenomTerm;

	dDbG.x *= dDenomTerm;
	dDbG.y *= dDenomTerm;
	dDbG.z *= dDenomTerm;

	// Cross product gives local Jacobian
	Node nodeCross = CrossProduct(dDaG, dDbG);

	double dJacobian = sqrt(
		  nodeCross.x * nodeCross.x
		+ nodeCross.y * nodeCross.y
		+ nodeCross.z * nodeCross.z);

	return 0.5 * dJacobian;
}

///////////////////////////////////////////////////////////////////////////////

Real CalculateFaceAreaQuadratureMethod(
	const Face & face,
	const NodeVector & nodes
) {
	int nTriangles = face.size() - 2;

	const int nOrder = 6;

	DataArray1D<double> dG;
	DataArray1D<double> dW;
	GaussQuadrature::GetPoints(nOrder, 0.0, 1.0, dG, dW);

	double dFaceArea = 0.0;

	// Loop over all sub-triangles of this Face
	for (int j = 0; j < nTriangles; j++) {

		// Calculate the area of the modified Face
		Node node1 = nodes[face[0]];
		Node node2 = nodes[face[j+1]];
		Node node3 = nodes[face[j+2]];

		// Calculate area at quadrature node
		for (int p = 0; p < dW.GetRows(); p++) {
		for (int q = 0; q < dW.GetRows(); q++) {

			double dA = dG[p];
			double dB = dG[q];

			double dJacobian =
				CalculateSphericalTriangleJacobian(
					node1, node2, node3,
					dA, dB);

			dFaceArea += dW[p] * dW[q] * dJacobian;
		}
		}
	}

	return dFaceArea;
}

///////////////////////////////////////////////////////////////////////////////

Real CalculateFaceAreaTriQuadratureMethod(
	const Face & face,
	const NodeVector & nodes
) {

	int nOrder1 = 4;
	int nOrder2 = 8;

	double dFaceArea = 0.0;

	double h = MaxEdgeLength(face,nodes);

	if(h < 0.004){

		dFaceArea = CalculateFaceAreaTriQuadrature(face,nodes,nOrder1);

	}
	else if(h < 0.09){

		dFaceArea = CalculateFaceAreaTriQuadrature(face,nodes,nOrder2);

	}
	else{

		FaceVector faces;
		faces.push_back(face);
		dFaceArea = CalculateFaceAreaTriQuadratureSplit(faces,nodes,nOrder2);

	}

	return dFaceArea;
}

///////////////////////////////////////////////////////////////////////////////

Real CalculateFaceAreaTriQuadrature(
	const Face & face,
	const NodeVector & nodes,
	int nOrder
) {

	double dFaceArea = 0.0;

	TriangularQuadratureRule TriQuadPoints(nOrder);

	DataArray2D<double> dG1 = TriQuadPoints.GetG();
	
	DataArray1D<double> dW1 = TriQuadPoints.GetW();

	int nTriangles = face.size() - 2;

	for (int j = 0; j < nTriangles; j++) {

		Node node1 = nodes[face[0]];
		Node node2 = nodes[face[j+1]];
		Node node3 = nodes[face[j+2]];

		Node nodeCross = CrossProduct(node2,node3);

		double tri_pro = node1.x*nodeCross.x + node1.y*nodeCross.y + node1.z*nodeCross.z;

		for (int q =0; q < dW1.GetRows(); q++){
						
			Node pnts_q = node1*dG1[q][0] + node2*dG1[q][1] + node3*dG1[q][2];

			double nrm_q = pnts_q.Magnitude();

			Node sph_q = pnts_q/nrm_q;

			double nrm_q3 = pow(nrm_q,3);
			
			dFaceArea += 0.5*(tri_pro/nrm_q3)*dW1[q];
		}

	}

	return dFaceArea;


}

///////////////////////////////////////////////////////////////////////////////

Real CalculateFaceAreaTriQuadratureSplit(
	const FaceVector & faces,
	const NodeVector & nodes,
	int & nOrder
) {

	double tol = 0.0;
	double dFaceArea = 0.0;

	if(nOrder >= 8){
		tol = 0.05;
	}
	
	else{
		tol = 0.003;
	}

	double dArea1 = 0.0;

	int nf = faces.size();

	for (int i = 0; i < nf; i++){

		int nv_surf = faces[i].size();
		int nTriangles = faces[i].size() - 2;
		double h = MaxEdgeLength(faces[i],nodes);

		if(h > tol){

			FaceVector surf_fid;
			NodeVector pnts_vor;
			Node node;

			for (int j = 0; j < nv_surf; j++){
				
				pnts_vor.push_back(nodes[faces[i][j]]);
				
			}

		   	// insert elements and points
			int index = nv_surf;
			node = ((pnts_vor[0] + pnts_vor[1])/2.0)/(((pnts_vor[0] + pnts_vor[1])/2.0).Magnitude());
			pnts_vor.push_back(node);

			for (int j = 1; j < nv_surf-1; j++){

				// insert elements
				DataArray2D<int> surf_index;
				surf_index.Allocate(4,3);

				Face face1(3);
				Face face2(3);
				Face face3(3);
				Face face4(3);

				face1.SetNode(0,0), face1.SetNode(1,index), face1.SetNode(2,index+2);
				face2.SetNode(0,index+2), face2.SetNode(1,index), face2.SetNode(2,index+1);
				face3.SetNode(0,index+1), face3.SetNode(1,index), face3.SetNode(2,j);
				face4.SetNode(0,index+2), face4.SetNode(1,index+1), face4.SetNode(2,j+1);

				surf_fid.push_back(face1);
				surf_fid.push_back(face2);
				surf_fid.push_back(face3);
				surf_fid.push_back(face4);

				index += 1;
				node = ((pnts_vor[j] + pnts_vor[j+1])/2.0)/(((pnts_vor[j] + pnts_vor[j+1])/2.0).Magnitude());
				pnts_vor.push_back(node);
				index += 1;
				node = ((pnts_vor[0] + pnts_vor[j+1])/2.0)/(((pnts_vor[0] + pnts_vor[j+1])/2.0).Magnitude());
				pnts_vor.push_back(node);		

			}

			double dArea2 = CalculateFaceAreaTriQuadratureSplit(surf_fid,pnts_vor,nOrder);
			dArea1 += dArea2;

		}
		
		else{

			double dArea2 = CalculateFaceAreaTriQuadrature(faces[i],nodes,nOrder);
			dArea1 += dArea2;
			
		}
		
	}

	return dArea1;
}

///////////////////////////////////////////////////////////////////////////////

double MaxEdgeLength(
	const Face & face,
	const NodeVector & nodes
) {
	int nv_surf = face.size();

	double h = 0.0;

	for (int i = 0; i < nv_surf-1; i++){
		
		Node node_diff = nodes[face[i]]-nodes[face[i+1]];
		h = fmax(h,node_diff.Magnitude());
		
	}
	
	Node node_diff = nodes[face[0]]-nodes[face[nv_surf-1]];
	h = fmax(h,node_diff.Magnitude());
	
	return h;
}

///////////////////////////////////////////////////////////////////////////////

Node AverageFaceNodes(
	const Face & face,
	const NodeVector & nodes
) {
	Node nodeAverage;
	for (size_t k = 0; k < face.size(); k++) {
		const Node & nodeVertex = nodes[face[k]];
		nodeAverage.x += nodeVertex.x;
		nodeAverage.y += nodeVertex.y;
		nodeAverage.z += nodeVertex.z;
	}

	double dMag = nodeAverage.Magnitude();

	nodeAverage.x /= dMag;
	nodeAverage.y /= dMag;
	nodeAverage.z /= dMag;

	return nodeAverage;
}

///////////////////////////////////////////////////////////////////////////////

bool IsFaceConcave(
	const Face & face,
	const NodeVector & nodes
) {
	const int nEdges = face.size();

	MeshUtilitiesFuzzy meshutils;

	bool fHasReflexNodes = false;

	// Search for reflex nodes on Face
	for (int i = 0; i < nEdges; i++) {
		
		int ixLast = (i + nEdges - 1) % nEdges;
		int ixCurr = i;
		int ixNext = (i + 1) % nEdges;

		const Node & nodeLast = nodes[face[ixLast]];
		const Node & nodeCurr = nodes[face[ixCurr]];
		const Node & nodeNext = nodes[face[ixNext]];

		int iSide = meshutils.FindNodeEdgeSide(
			nodeLast,
			nodeCurr,
			Edge::Type_GreatCircleArc,
			nodeNext);

		if (iSide != (-1)) {
			continue;
		}

		return true;
	}

	return false;
}
*/
///////////////////////////////////////////////////////////////////////////////

Real CalculateFaceArea_Concave(
	const Face & face,
	const NodeVector & nodes
) {
/*
	Real dArea = 0.0;

	if (IsFaceConcave(face, nodes)) {

		// Generate a new Mesh only including this Face
		Mesh mesh;
		for (int j = 0; j < face.size(); j++) {
			mesh.nodes.push_back(nodes[face[j]]);
		}
		Face faceTemp(face.size());
		for (int j = 0; j < faceTemp.size(); j++) {
			faceTemp.SetNode(j, j);
		}
		mesh.faces.push_back(faceTemp);

		// Convexify this Face
		Mesh meshout;
		ConvexifyFace(mesh, meshout, 0, false, true);

		if (meshout.faces.size() == 0) {
			_EXCEPTIONT("Call to ConvexifyFace() failed; no convex mesh generated");
		}
		return meshout.CalculateFaceAreas(false);
	}
*/
	return CalculateFaceArea(face, nodes);
}

///////////////////////////////////////////////////////////////////////////////

Real CalculateFaceArea(
	const Face & face,
	const NodeVector & nodes
) {
	return 0.0;
	//return CalculateFaceAreaQuadratureMethod(face, nodes);
}

///////////////////////////////////////////////////////////////////////////////
/*
bool ConvexifyFace(
	Mesh & mesh,
	Mesh & meshout,
	int iFace,
	bool fRemoveConcaveFaces,
	bool fVerbose
) {
	Face & face = mesh.faces[iFace];
	const int nNodes = face.size()-1;
	if(fVerbose) {
		Announce("ConvexifyFace via Triangle package");
		Announce("iFace=%i	nNodes: %i", iFace, nNodes);
	}

	// get center of this face and the local up vector, Z
	Node center(0,0,0);
	for (int i=0; i<nNodes; ++i) center = center + mesh.nodes[face[i]];
	center = center / nNodes;
	Node localZ = center.Normalized();

	// construct tangent space X and Y unit vectors
	Node node0 = (mesh.nodes[face[0]] - center).Normalized();
	Node localY = CrossProduct(localZ,node0);
	Node localX = CrossProduct(localY,localZ);

	// get orthographic projection of nodes onto the tangent plane
	NodeVector planarNodes;
	for (int i=0; i<nNodes; ++i) {
		Node & node3D = mesh.nodes[face[i]];
		Node node2D( DotProduct(node3D,localX), DotProduct(node3D,localY), 0);
		planarNodes.push_back(node2D);
	}

	// fill in triangleio data structures
	struct triangulateio in, out, vorout;

	// initialize data structure for input planar straight-line graph (PSLG)
	in.numberofpoints		   = nNodes;
	in.numberofpointattributes  = 0;
	in.numberofsegments		 = nNodes;
	in.numberofholes			= 0;
	in.numberofregions		  = 0;
	in.pointlist				= (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
	in.segmentlist			  = (int  *) malloc(in.numberofsegments * 2 * sizeof(int));;
	in.pointattributelist	   = (REAL *) NULL;
	in.pointmarkerlist		  = (int  *) NULL;
	in.trianglelist			 = (int  *) NULL;
	in.triangleattributelist	= (REAL *) NULL;
	in.neighborlist			 = (int  *) NULL;
	in.segmentmarkerlist		= (int  *) NULL;
	in.edgelist				 = (int  *) NULL;
	in.edgemarkerlist		   = (int  *) NULL;

	// initialize data structure for output triangulation
	out.pointlist			   = (REAL *) NULL;
	out.pointattributelist	  = (REAL *) NULL;
	out.pointmarkerlist		 = (int  *) NULL;
	out.trianglelist			= (int  *) NULL;
	out.triangleattributelist   = (REAL *) NULL;
	out.neighborlist			= (int  *) NULL;
	out.segmentlist			 = (int  *) NULL;
	out.segmentmarkerlist	   = (int  *) NULL;
	out.edgelist				= (int  *) NULL;
	out.edgemarkerlist		  = (int  *) NULL;

	// initialize data structure for output Voronoi diagram (unused)
	vorout.pointlist			= (REAL *) NULL;
	vorout.pointattributelist   = (REAL *) NULL;
	vorout.edgelist			 = (int  *) NULL;
	vorout.normlist			 = (REAL *) NULL;

	// fill in 2d point list
	for(int i=0; i<nNodes; ++i) {
		const Node & n = planarNodes[i];
		in.pointlist[i*2+0] = n.x;
		in.pointlist[i*2+1] = n.y;
	}

	// fill in segment list
	for(int i=0; i<nNodes; ++i) {
		in.segmentlist[i*2+0] = i;
		in.segmentlist[i*2+1] = (i+1)%nNodes;
	}

	// set options for triangulate function call:
	// p   -> triangulate area in the boundary (PSLG)
	// q5  -> set min triangle angle to 5dg
	// j   -> jettison unused nodes
	// z   -> number nodes starting from zero
	// Y   -> no new nodes on the boundary (so it remains conforming)
	// Q,V -> quiet or verbose output

	if (fVerbose) {
		char options[256] ="pq5jzYV";
		AnnounceBanner();
		triangulate(options, &in, &out, &vorout);
		AnnounceBanner();
	}
	else {
		char options[256] ="pq5jzYQ";
		triangulate(options, &in, &out, &vorout);
	}

	// project new planar nodes onto the unit sphere
	NodeVector newNodes;
	for (int i=0; i<out.numberofpoints; ++i) {
		Node n(out.pointlist[2*i], out.pointlist[2*i+1],0.0);
		Real z = sqrt(1.0 - n.x*n.x - n.y*n.y);
		Node node3D = localZ * z + (localX * n.x) + (localY * n.y);
		node3D = node3D.Normalized();
		newNodes.push_back(node3D);
	}

	// Clean up
	free(in.pointlist);
	free(in.segmentlist);

	// delete concave face from the mesh
	if (fRemoveConcaveFaces) {
		_EXCEPTION();
		//mesh.faces.erase(mesh.faces.begin() + iFace);
	}

	// append new nodes to end of node vector
	int size = meshout.nodes.size();
	meshout.nodes.insert(meshout.nodes.end(), newNodes.begin(), newNodes.end());

	// add new triangles to the mesh
	for (int i=0; i<out.numberoftriangles; ++i) {
		Face newFace(3);
		newFace.SetNode(0, size+out.trianglelist[i*3+0]);
		newFace.SetNode(1, size+out.trianglelist[i*3+1]);
		newFace.SetNode(2, size+out.trianglelist[i*3+2]);
		meshout.faces.push_back(newFace);
	}

	// clean up duplicate nodes created on the boundary
	meshout.RemoveCoincidentNodes();

	return true;
}

///////////////////////////////////////////////////////////////////////////////

bool ConvexifyFaceBayazit(
	Mesh & mesh,
	Mesh & meshout,
	int iFace,
	bool fRemoveConcaveFaces,
	bool fVerbose
) {
	if ((iFace < 0) || (iFace > mesh.faces.size())) {
		_EXCEPTIONT("Face index out of range");
	}

	Face & face = mesh.faces[iFace];

	const int nEdges = face.size();

	MeshUtilitiesFuzzy meshutils;

	bool fHasReflexNodes = false;

	// Search for reflex nodes on Face
	for (int i = 0; i < nEdges; i++) {
		
		int ixLast = (i + nEdges - 1) % nEdges;
		int ixCurr = i;
		int ixNext = (i + 1) % nEdges;

		const Node & nodeLast = mesh.nodes[face[ixLast]];
		const Node & nodeCurr = mesh.nodes[face[ixCurr]];
		const Node & nodeNext = mesh.nodes[face[ixNext]];

		int iSide = meshutils.FindNodeEdgeSide(
			nodeLast,
			nodeCurr,
			Edge::Type_GreatCircleArc,
			nodeNext);

		if (iSide != (-1)) {
			continue;
		}

		if (fVerbose) {
			Announce("Reflex node found: %i", face[ixCurr]);
		}

		// Reflex node found; divide mesh at this node
		int ixDividingNode = (-1);
		double dMinDist = (-1.0);
		for (int j = 0; j < nEdges; j++) {
			if ((j == ixLast) || (j == ixCurr) || (j == ixNext)) {
				continue;
			}

			const Node & nodeCandidate = mesh.nodes[face[j]];

			Node node2;
			node2.x = nodeCurr.x + (nodeNext.x - nodeLast.x);
			node2.y = nodeCurr.y + (nodeNext.y - nodeLast.y);
			node2.z = nodeCurr.z + (nodeNext.z - nodeLast.z);

			double dMag2 = node2.Magnitude();
			node2.x /= dMag2;
			node2.y /= dMag2;
			node2.z /= dMag2;

			// Check that this Node is in the range of the reflex node
			const int iSide0 = meshutils.FindNodeEdgeSide(
				nodeLast,
				nodeCurr,
				Edge::Type_GreatCircleArc,
				nodeCandidate);

			const int iSide1 = meshutils.FindNodeEdgeSide(
				nodeCurr,
				nodeNext,
				Edge::Type_GreatCircleArc,
				nodeCandidate);

			if (iSide0 == (-1)) {
				continue;
			}

			if (iSide1 == (-1)) {
				continue;
			}

			// Check that this Node is of minimum distance
			Node nodeDelta = nodeCandidate - nodeCurr;

			double dDist = nodeDelta.Magnitude();

			if ((dMinDist < 0.0) || (dDist < dMinDist)) {
				ixDividingNode = j;
				dMinDist = dDist;
			}
		}

		// No dividing node found -- add a Steiner vertex
		if (ixDividingNode == (-1)) {

			if (fVerbose) {
				Announce("No dividing node found -- adding a Steiner vertex");
			}

			// Find a node that bisects the two angles
			Node nodeBisect;
			nodeBisect.x = 3.0 * nodeCurr.x - nodeLast.x - nodeNext.x;
			nodeBisect.y = 3.0 * nodeCurr.y - nodeLast.y - nodeNext.y;
			nodeBisect.z = 3.0 * nodeCurr.z - nodeLast.z - nodeNext.z;

			double dBisectMag = nodeBisect.Magnitude();
			nodeBisect.x /= dBisectMag;
			nodeBisect.y /= dBisectMag;
			nodeBisect.z /= dBisectMag;

			int ixIntersectEdge;
			Node nodeClosestIntersect;
			dMinDist = (-1.0);

			for (int j = 0; j < nEdges; j++) {
				if ((j == ixCurr) || (j == ixLast)) {
					continue;
				}

				std::vector<Node> vecIntersections;

				bool fCoincident = meshutils.CalculateEdgeIntersectionsSemiClip(
					mesh.nodes[face[j]],
					mesh.nodes[face[(j+1)%nEdges]],
					Edge::Type_GreatCircleArc,
					nodeCurr,
					nodeBisect,
					Edge::Type_GreatCircleArc,
					vecIntersections);

				if (fCoincident) {
					_EXCEPTIONT("Coincident lines detected");
				}

				if (vecIntersections.size() > 1) {
					_EXCEPTIONT("Logic error");

				} else if (vecIntersections.size() == 1) {
					Node nodeDelta = vecIntersections[0] - nodeCurr;
					double dDist = nodeDelta.Magnitude();

					if ((dMinDist == -1.0) || (dDist < dMinDist)) {
						dMinDist = dDist;
						ixIntersectEdge = j;
						nodeClosestIntersect = vecIntersections[0];
					}
				}
			}

			Edge edgeIntersect = face.edges[ixIntersectEdge];

			if (dMinDist == -1.0) {
				_EXCEPTIONT("Logic error: No intersecting lines found");
			}

			int ixNewIntersectNode = mesh.nodes.size();
			meshout.nodes.push_back(nodeClosestIntersect);

			int iFaceSize1 = 1;
			for (int k = (i+1)%nEdges; face[k] != edgeIntersect[0]; k = (k+1)%nEdges) {
				iFaceSize1++;
			}

			Face faceNew1(iFaceSize1 + 2);
			Face faceNew2(nEdges - iFaceSize1 + 1);

			for (int k = 0; k < iFaceSize1 + 1; k++) {
				faceNew1.SetNode(k, face[(i+k)%nEdges]);
				//printf("%i\n", (i+k)%nEdges);
			}
			faceNew1.SetNode(iFaceSize1 + 1, ixNewIntersectNode);
			for (int k = 0; k < nEdges - iFaceSize1; k++) {
				faceNew2.SetNode(k, face[(ixIntersectEdge+k+1)%nEdges]);
				//printf("%i\n", (ixIntersectEdge+k+1)%nEdges);
			}
			faceNew2.SetNode(nEdges - iFaceSize1, ixNewIntersectNode);

			meshout.faces.push_back(faceNew1);
			meshout.faces.push_back(faceNew2);

			if (fRemoveConcaveFaces) {
				mesh.faces.erase(mesh.faces.begin() + iFace);
			}

			int nFaces = meshout.faces.size();
			ConvexifyFace(meshout, meshout, nFaces-1, true, fVerbose);
			ConvexifyFace(meshout, meshout, nFaces-2, true, fVerbose);

		// Divide the mesh at the reflex node
		} else {
			if (fVerbose) {
				Announce("Dividing node found %i", face[ixDividingNode]);
			}

			int iFaceSize1 = 0;
			for (int k = (i+1)%nEdges; k != ixDividingNode; k = (k+1)%nEdges) {
				iFaceSize1++;
			}

			Face faceNew1(iFaceSize1 + 2);
			Face faceNew2(nEdges - iFaceSize1);

			for (int k = 0; k < iFaceSize1 + 2; k++) {
				faceNew1.SetNode(k, face[(i+k)%nEdges]);
				//printf("%i\n", (i+k)%nEdges);
			}
			for (int k = 0; k < nEdges - iFaceSize1; k++) {
				faceNew2.SetNode(k, face[(ixDividingNode+k)%nEdges]);
				//printf("%i\n", (ixDividingNode+k)%nEdges);
			}

			meshout.faces.push_back(faceNew1);
			meshout.faces.push_back(faceNew2);

			if (fRemoveConcaveFaces) {
				mesh.faces.erase(mesh.faces.begin() + iFace);
			}

			int nFaces = meshout.faces.size();
			ConvexifyFace(meshout, meshout, nFaces-1, true, fVerbose);
			ConvexifyFace(meshout, meshout, nFaces-2, true, fVerbose);
		}

		fHasReflexNodes = true;
		break;
	}

	return fHasReflexNodes;
}

///////////////////////////////////////////////////////////////////////////////

void ConvexifyMesh(
	Mesh & mesh,
	Mesh & meshout,
	bool fVerbose
) {
	char szBuffer[256];

	// Copy all nodes to output mesh
	meshout.nodes.clear();

	// Remove all Faces from output mesh
	meshout.faces.clear();

	// Clear the MultiFaceMap
	meshout.vecMultiFaceMap.clear();

	// Loop through all Faces in the input Mesh
	int nFaces = mesh.faces.size();
	for (int f = 0; f < nFaces; f++) {
		if (fVerbose) {
			snprintf(szBuffer, 256, "Face %i", f);
			AnnounceStartBlock(szBuffer);
		}

		// Adjust current Face index
		int nMeshSize = meshout.faces.size();

		bool fConcaveFace = ConvexifyFace(mesh, meshout, f, false, fVerbose);

		if (fConcaveFace) {
			int nAddedFaces = meshout.faces.size() - nMeshSize;
			for (int i = 0; i < nAddedFaces; i++) {
				meshout.vecMultiFaceMap.push_back(f);
			}

		} else {
			meshout.faces.push_back(mesh.faces[f]);
			meshout.vecMultiFaceMap.push_back(f);
		}

		if (fVerbose) {
			AnnounceEndBlock("Done");
		}
	}

	if (meshout.vecMultiFaceMap.size() != meshout.faces.size()) {
		_EXCEPTIONT("Logic error");
	}
}

///////////////////////////////////////////////////////////////////////////////

void Dual(
	Mesh & mesh
) {
	
	// Generate ReverseNodeArray
	mesh.ConstructReverseNodeArray();

	// Backup Nodes and Faces
	NodeVector nodesOld = mesh.nodes;
	FaceVector facesOld = mesh.faces;

	mesh.nodes.clear();
	mesh.faces.clear();
		
	// Generate new Node array
	for (int i = 0; i < facesOld.size(); i++) {
		Node node;
		for (int j = 0; j < facesOld[i].size(); j++) {
			node.x += nodesOld[facesOld[i][j]].x;
			node.y += nodesOld[facesOld[i][j]].y;
			node.z += nodesOld[facesOld[i][j]].z;
		}

		node.x /= static_cast<double>(facesOld[i].size());
		node.y /= static_cast<double>(facesOld[i].size());
		node.z /= static_cast<double>(facesOld[i].size());

		double dMag = node.Magnitude();

		node.x /= dMag;
		node.y /= dMag;
		node.z /= dMag;

		mesh.nodes.push_back(node);
	}

	// Generate new Face array
	for (int i = 0; i < nodesOld.size(); i++) {
		const int nEdges = mesh.revnodearray[i].size();

		Face face(mesh.revnodearray[i].size());
		Face faceTemp(mesh.revnodearray[i].size());

		int ixNode = 0;
		std::set<int>::const_iterator iter = mesh.revnodearray[i].begin();
		for (; iter != mesh.revnodearray[i].end(); iter++) {
			faceTemp.SetNode(ixNode, *iter);
			ixNode++;
		}

		// Reorient Faces
		Node nodeCentral = nodesOld[i];
		Node node0 = mesh.nodes[faceTemp[0]] - nodeCentral;

		Node nodeCross = CrossProduct(mesh.nodes[faceTemp[0]], nodeCentral);
		
		double dNode0Mag = node0.Magnitude();

		// Determine the angles about the central Node of each Face Node
		std::vector<double> dAngles;
		dAngles.resize(faceTemp.size());
		dAngles[0] = 0.0;

		for (int j = 1; j < nEdges; j++) {
			Node nodeDiff = mesh.nodes[faceTemp[j]] - nodeCentral;
			double dNodeDiffMag = nodeDiff.Magnitude();

			double dSide = DotProduct(nodeCross, nodeDiff);

			double dDotNorm =
				DotProduct(node0, nodeDiff) / (dNode0Mag * dNodeDiffMag);

			double dAngle;
			if (dDotNorm > 1.0) {
				dDotNorm = 1.0;
			}

			dAngles[j] = acos(dDotNorm);

			if (dSide > 0.0) {
				dAngles[j] = - dAngles[j] + 2.0 * M_PI;
			}
		}

		// Orient each Face by putting Nodes in order of increasing angle
		double dCurrentAngle = 0.0;
		face.SetNode(0, faceTemp[0]);
		for (int j = 1; j < nEdges; j++) {
			int ixNextNode = 1;
			double dNextAngle = 2.0 * M_PI;

			for (int k = 1; k < nEdges; k++) {
				if ((dAngles[k] > dCurrentAngle) && (dAngles[k] < dNextAngle)) {
					ixNextNode = k;
					dNextAngle = dAngles[k];
				}
			}

			face.SetNode(j, faceTemp[ixNextNode]);
			dCurrentAngle = dNextAngle;
		}

		mesh.faces.push_back(face);
	}
}

///////////////////////////////////////////////////////////////////////////////

void ConstructLocalDualFace(
	const Mesh & mesh,
	NodeVector & meshCenters,
	int & iNodeX,
	Face & faceLocalDual,
	NodeVector & nodesFaceLocal
) {
		
	std::set<int>::const_iterator iterRevNode = mesh.revnodearray[iNodeX].begin();
	
	for (; iterRevNode != mesh.revnodearray[iNodeX].end(); iterRevNode++) {
		
		Node node;
		
		for (int j = 0; j < mesh.faces[*iterRevNode].size(); j++){
		
			node.x += mesh.nodes[mesh.faces[*iterRevNode][j]].x;
			node.y += mesh.nodes[mesh.faces[*iterRevNode][j]].y;
			node.z += mesh.nodes[mesh.faces[*iterRevNode][j]].z;
			
		}
		
		node.x /= static_cast<double>(mesh.faces[*iterRevNode].size());
		node.y /= static_cast<double>(mesh.faces[*iterRevNode].size());
		node.z /= static_cast<double>(mesh.faces[*iterRevNode].size());

		double dMag = node.Magnitude();

		node.x /= dMag;
		node.y /= dMag;
		node.z /= dMag;

		nodesFaceLocal.push_back(node);
		
		
	}

	const int nEdges = mesh.revnodearray[iNodeX].size();
	
	
	//Face face(mesh.revnodearray[iNodeX].size());
	Face faceTemp(mesh.revnodearray[iNodeX].size());

	int ixNode = 0;
	std::set<int>::const_iterator iter = mesh.revnodearray[iNodeX].begin();
	for (; iter != mesh.revnodearray[iNodeX].end(); iter++) {
		faceTemp.SetNode(ixNode, *iter);
		ixNode++;
	}

	// Reorient Faces
	Node nodeCentral = mesh.nodes[iNodeX];
	
	Node node0 = meshCenters[faceTemp[0]] - nodeCentral;

	Node nodeCross = CrossProduct(meshCenters[faceTemp[0]], nodeCentral);
	
	double dNode0Mag = node0.Magnitude();

	// Determine the angles about the central Node of each Face Node
	std::vector<double> dAngles;
	dAngles.resize(faceTemp.size());
	dAngles[0] = 0.0;

	for (int j = 1; j < nEdges; j++) {
		Node nodeDiff = meshCenters[faceTemp[j]] - nodeCentral;
		double dNodeDiffMag = nodeDiff.Magnitude();

		double dSide = DotProduct(nodeCross, nodeDiff);

		double dDotNorm =
			DotProduct(node0, nodeDiff) / (dNode0Mag * dNodeDiffMag);

		double dAngle;
		if (dDotNorm > 1.0) {
			dDotNorm = 1.0;
		}

		dAngles[j] = acos(dDotNorm);

		if (dSide > 0.0) {
			dAngles[j] = - dAngles[j] + 2.0 * M_PI;
		}
	}

	// Orient each Face by putting Nodes in order of increasing angle
	double dCurrentAngle = 0.0;
	faceLocalDual.SetNode(0, faceTemp[0]);
	for (int j = 1; j < nEdges; j++) {
		int ixNextNode = 1;
		double dNextAngle = 2.0 * M_PI;

		for (int k = 1; k < nEdges; k++) {
			if ((dAngles[k] > dCurrentAngle) && (dAngles[k] < dNextAngle)) {
				ixNextNode = k;
				dNextAngle = dAngles[k];
			}
		}

		faceLocalDual.SetNode(j, faceTemp[ixNextNode]);
		dCurrentAngle = dNextAngle;
	}
}
*/
///////////////////////////////////////////////////////////////////////////////
