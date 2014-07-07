/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

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
int flutter_convergence_driver (libMesh::LibMeshInit& init, GetPot& infile,
                                int argc, char* const argv[]);



int main (int argc, char* const argv[])
{
    libMesh::LibMeshInit init(argc, argv);
    
    // get the input file
    std::string nm = libMesh::command_line_value("-i", std::string("system_input.in"));
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
    else if (type == "flutter_convergence")
        rval = flutter_convergence_driver(init, infile, argc, argv);
    else {
        libMesh::out
        << "Invalid analysis type: "
        << type << std::endl
        << "Stopping..." << std::endl;
        libmesh_stop();
    }
    
    return rval;
}

