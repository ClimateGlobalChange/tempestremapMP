///////////////////////////////////////////////////////////////////////////////
///
///	\file	GenerateOverlapMeshExe.cpp
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

#include "Announce.h"
#include "GridElements.h"
#include "TriangularQuadrature.h"

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error ( NcError::silent_nonfatal );

	try
	{
		AnnounceStartBlock("Loading mesh");
		Mesh mesh("outCSne30.g");
		AnnounceEndBlock("Done");

		AnnounceStartBlock("Calculating grid area");
		mesh.calculate_face_areas();
		AnnounceEndBlock("Done");

		Announce("Total Area: %1.15e", mesh.calculate_total_area());

	}
	catch (Exception & e)
	{
		Announce ( e.ToString().c_str() );
		return ( 0 );
	}
	catch ( ... )
	{
		Announce("Unhandled exception");
		return ( 0 );
	}
}

///////////////////////////////////////////////////////////////////////////////
