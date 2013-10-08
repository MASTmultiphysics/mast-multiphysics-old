//
//  main.cpp
//  RealSolver
//
//  Created by Manav Bhatia on 9/23/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/getpot.h"

// driver functions
int fluid_driver (LibMeshInit& init, GetPot& infile,
                  int argc, char* const argv[]);
int potential_fluid_driver (LibMeshInit& init, GetPot& infile,
                            int argc, char* const argv[]);
int static_structural_driver (LibMeshInit& init, GetPot& infile,
                              int argc, char* const argv[]);
int modal_structural_driver (LibMeshInit& init, GetPot& infile,
                             int argc, char* const argv[]);
int flutter_driver (LibMeshInit& init, GetPot& infile,
                    int argc, char* const argv[]);

using namespace libMesh;

int main (int argc, char* const argv[])
{
    LibMeshInit init(argc, argv);
    
    //sleep(100);
    
    // get the input file
    std::string nm = command_line_value("-i", std::string("system_input.in"));
    GetPot infile(nm);
    int rval = 0;
    
    std::string type = infile("analysis_type", "");
    if (type == "fluid")
        rval = fluid_driver(init, infile, argc, argv);
    else if (type == "compressible_potential_fluid")
        rval = potential_fluid_driver(init, infile, argc, argv);
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    else if (type == "structures_static")
        rval = static_structural_driver(init, infile, argc, argv);
    else if (type == "structures_modal")
        rval = modal_structural_driver(init, infile, argc, argv);
#endif
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    else if (type == "flutter")
        rval = flutter_driver(init, infile, argc, argv);
#endif
    else {
        libMesh::out
        << "Invalid analysis type: "
        << type << std::endl
        << "Stopping..." << std::endl;
        libmesh_stop();
    }
    
    return rval;
}