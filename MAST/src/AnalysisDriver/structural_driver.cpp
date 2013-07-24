//
//  run_euler.cpp
//  FESystem
//
//  Created by Manav Bhatia on 3/6/13.
//
//

// C++ includes
#include <iomanip>

// FESystem includes
#include "StructuralElems/structural_system_base.h"

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
#include "libmesh/string_to_enum.h"



// The main program.
int main (int argc, char* const argv[])
{
    // Initialize libMesh.
    LibMeshInit init (argc, argv);

    SerialMesh mesh(init.comm());
    mesh.set_mesh_dimension(2);

    MeshTools::Generation::build_square (mesh,
                                         10,
                                         10,
                                         0., 1.,
                                         0., 1.,
                                         TRI3);

    mesh.prepare_for_use();
    
    // Print information about the mesh to the screen.
    mesh.print_info();
    
    
    // Create an equation systems object.
    EquationSystems equation_systems (mesh);
    
    // Declare the system "EulerSystem"
    StructuralSystemBase & system =
    equation_systems.add_system<StructuralSystemBase> ("StructuralSystem");
    
    system.time_solver = AutoPtr<TimeSolver>(new SteadySolver(system));
    
    equation_systems.init ();

    system.time_solver->diff_solver()->quiet = false;
    system.time_solver->diff_solver()->verbose = true;
    system.time_solver->diff_solver()->relative_residual_tolerance =  1.0e-8;
    system.time_solver->diff_solver()->absolute_residual_tolerance =  1.0e-8;
    
    
    // Print information about the system to the screen.
    equation_systems.print_info();
    
    system.solve();
        
    // We write the file in the ExodusII format.
    Nemesis_IO(mesh).write_equation_systems("out.exo",
                                            equation_systems);
    // All done.
    return 0;
}

