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
#include "libmesh/eigen_time_solver.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/condensed_eigen_system.h"



// The main program.
int main_static (int argc, char* const argv[])
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
    
    // Declare the system
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


void assemble_matrices(EquationSystems& es,
                       const std::string& system_name);

void get_dirichlet_dofs(EquationSystems& es,
                        const std::string& system_name,
                        std::set<unsigned int>& dirichlet_dof_ids);


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
    
    CondensedEigenSystem & eigen_system =
    equation_systems.add_system<CondensedEigenSystem> ("Eigensystem");
    
    GetPot infile("structural.in");
    
    eigen_system.add_variable ( "ux", FIRST, LAGRANGE);
    eigen_system.add_variable ( "uy", FIRST, LAGRANGE);
    eigen_system.add_variable ( "uz", FIRST, LAGRANGE);
    eigen_system.add_variable ( "tx", FIRST, LAGRANGE);
    eigen_system.add_variable ( "ty", FIRST, LAGRANGE);
    eigen_system.add_variable ( "tz", FIRST, LAGRANGE);
    
    // Give the system a pointer to the matrix assembly
    // function defined below.
    eigen_system.attach_assemble_function (assemble_matrices);
    
    // Set the type of the problem, here we deal with
    // a generalized Hermitian problem.
    eigen_system.set_eigenproblem_type(GHEP);
    
    // Order the eigenvalues "smallest first"
    eigen_system.eigen_solver->set_position_of_spectrum(LARGEST_MAGNITUDE);

    // Set the number of requested eigenpairs \p n_evals and the number
    // of basis vectors used in the solution algorithm.
    equation_systems.parameters.set<unsigned int>("eigenpairs")    = 10;
    equation_systems.parameters.set<unsigned int>("basis vectors") = 10*3;
    
    // Initialize the data structures for the equation system.
    equation_systems.init();
    
    // Prints information about the system to the screen.
    equation_systems.print_info();
    
    // Pass the Dirichlet dof IDs to the CondensedEigenSystem
    std::set<unsigned int> dirichlet_dof_ids;
    get_dirichlet_dofs(equation_systems, "Eigensystem", dirichlet_dof_ids);
    eigen_system.initialize_condensed_dofs(dirichlet_dof_ids);
    
    // Solve the system "Eigensystem".
    eigen_system.solve();

    // Get the number of converged eigen pairs.
    unsigned int nconv = eigen_system.get_n_converged();

    for (unsigned int i=0; i<nconv; i++)
    {
        std::pair<Real, Real> val = eigen_system.get_eigenpair(i);
        
        std::cout << std::setw(5) << i
        << std::setw(10) << val.first
        << " + i  "
        << std::setw(10) << val.second << std::endl;
        
        std::ostringstream file_name;
        
        // We write the file in the ExodusII format.
        file_name << "out_"
        << std::setw(3)
        << std::setfill('0')
        << std::right
        << i
        << ".exo";
        
        // We write the file in the ExodusII format.
        Nemesis_IO(mesh).write_equation_systems(file_name.str(),
                                                equation_systems);
    }
    // All done.
    return 0;
}

