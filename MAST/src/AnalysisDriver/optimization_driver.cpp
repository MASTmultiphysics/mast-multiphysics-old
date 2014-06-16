//
//  optimization_driver.cpp
//  MAST
//
//  Created by Manav Bhatia on 12/20/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//


// C++ includes
#include <iomanip>

// MAST includes
#include "StructuralElems/structural_system_assembly.h"
#include "PropertyCards/element_property_card_base.h"
#include "PropertyCards/element_property_card_3D.h"
#include "PropertyCards/element_property_card_2D.h"
#include "PropertyCards/element_property_card_1D.h"
#include "PropertyCards/material_property_card_base.h"
#include "BoundaryConditions/boundary_condition.h"
#include "BoundaryConditions/displacement_boundary_condition.h"
#include "Mesh/mesh_initializer.h"
#include "Mesh/panel_mesh.h"
#include "Mesh/stiffened_panel.h"
#include "Mesh/nastran_io.h"
#include "AnalysisDriver/topology_optimization.h"
//#include "AnalysisDriver/sizing_optimization.h"
#include "Optimization/gcmma_optimization_interface.h"
#include "beam_postbuckling_sizing.h"

// libmesh includes
#include "libmesh/getpot.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/euler_solver.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/eigen_time_solver.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/xdr_io.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/zero_function.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/exodusII_io.h"




// The main program.
int optimization_driver (libMesh::LibMeshInit& init, GetPot& infile,
                       int argc, char* const argv[])
{
    std::ofstream output;
    output.open("optimization_output.txt", std::ofstream::out);
    
    MAST::GCMMAOptimizationInterface gcmma;
    
    // create and attach topology optimization object
    // MAST::TopologyOptimization func_eval(init, infile, output);
    
    // create and attach sizing optimization object
    MAST::SizingOptimization func_eval(init, infile, output);

    // attach and optimize
    gcmma.attach_function_evaluation_object(func_eval);
    gcmma.optimize();
    
    output.close();
    
    // All done.
    return 0;
}


