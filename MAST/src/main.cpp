//
//  main.cpp
//  MAST
//
//  Created by Manav Bhatia on 9/23/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"

// driver functions
int fluid_driver (libMesh::LibMeshInit& init, GetPot& infile,
                  int argc, char* const argv[]);
int potential_fluid_driver (libMesh::LibMeshInit& init, GetPot& infile,
                            int argc, char* const argv[]);
int structural_driver (libMesh::LibMeshInit& init, GetPot& infile,
                       int argc, char* const argv[]);
int modal_structural_driver (libMesh::LibMeshInit& init, GetPot& infile,
                             int argc, char* const argv[]);
int flutter_driver (libMesh::LibMeshInit& init, GetPot& infile,
                    int argc, char* const argv[]);
int optimization_driver (libMesh::LibMeshInit& init, GetPot& infile,
                       int argc, char* const argv[]);

using namespace libMesh;

int main (int argc, char* const argv[])
{
    libMesh::LibMeshInit init(argc, argv);
    
    // get the input file
    std::string nm = command_line_value("-i", std::string("system_input.in"));
    GetPot infile(nm);
    int rval = 0;
    
    std::string type = infile("analysis_type", "");
    if (type == "fluid")
        rval = fluid_driver(init, infile, argc, argv);
    else if (type == "compressible_potential_fluid")
        rval = potential_fluid_driver(init, infile, argc, argv);

    else if (type == "structures")
        rval = structural_driver(init, infile, argc, argv);
    else if (type == "optimization")
        rval = optimization_driver(init, infile, argc, argv);
    else if (type == "flutter")
        rval = flutter_driver(init, infile, argc, argv);
    else {
        libMesh::out
        << "Invalid analysis type: "
        << type << std::endl
        << "Stopping..." << std::endl;
        libmesh_stop();
    }
    
    return rval;
}

