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
#include "libmesh/xdr_io.h"

#ifndef LIBMESH_USE_COMPLEX_NUMBERS

// The main program.
int static_structural_driver (LibMeshInit& init, GetPot& infile,
                              int argc, char* const argv[])
{
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


void assemble_plate_matrices(EquationSystems& es,
                             const std::string& system_name);
void assemble_beam_matrices(EquationSystems& es,
                            const std::string& system_name);

void get_plate_dirichlet_dofs(EquationSystems& es,
                              const std::string& system_name,
                              std::set<unsigned int>& dirichlet_dof_ids);

void get_beam_dirichlet_dofs(EquationSystems& es,
                             const std::string& system_name,
                             std::set<unsigned int>& dirichlet_dof_ids);


int modal_structural_driver (LibMeshInit& init, GetPot& infile,
                             int argc, char* const argv[])
{
    SerialMesh mesh(init.comm());

    const unsigned int dim = infile("dimension", 0);

    const Real pi = acos(-1.),
    x_length= infile("x_length", 10.),
    y_length= infile("y_length", 10.),
    z_length= infile("z_length", 10.),
    t_by_c =  infile("t_by_c", 0.05),
    chord =   infile("chord", 1.0),
    span =    infile("span", 1.0),
    thickness = 0.5*t_by_c*chord,
    x0=x_length*0.5-chord*0.5, x1=x0+chord, y0=y_length*0.5-span*0.5, y1=y0+span ;


    mesh.set_mesh_dimension(dim-1);
    
    switch (dim-1) {
        case 1:
            MeshTools::Generation::build_line (mesh,
                                               40,
                                               x0, x0+chord,
                                               EDGE2);
            break;

        case 2:
            MeshTools::Generation::build_square (mesh,
                                                 40,
                                                 40,
                                                 x0, x0+chord,
                                                 y0, y0+chord,
                                                 TRI3);
            break;
            
        default:
            libmesh_error();
            break;
    }



    mesh.prepare_for_use();
    
    // Print information about the mesh to the screen.
    mesh.print_info();
    
    
    // Create an equation systems object.
    EquationSystems equation_systems (mesh);
    
    CondensedEigenSystem & eigen_system =
    equation_systems.add_system<CondensedEigenSystem> ("Eigensystem");
    
    eigen_system.add_variable ( "ux", FIRST, LAGRANGE);
    eigen_system.add_variable ( "uy", FIRST, LAGRANGE);
    eigen_system.add_variable ( "uz", FIRST, LAGRANGE);
    eigen_system.add_variable ( "tx", FIRST, LAGRANGE);
    eigen_system.add_variable ( "ty", FIRST, LAGRANGE);
    eigen_system.add_variable ( "tz", FIRST, LAGRANGE);
    
    equation_systems.parameters.set<bool>("if_exchange_AB_matrices") = true;
    
    // Give the system a pointer to the matrix assembly
    // function defined below.
    
    // Set the type of the problem, here we deal with
    // a generalized Hermitian problem.
    eigen_system.set_eigenproblem_type(GHEP);
    
    // Order the eigenvalues "smallest first"
    eigen_system.eigen_solver->set_position_of_spectrum(LARGEST_MAGNITUDE);

    // Set the number of requested eigenpairs \p n_evals and the number
    // of basis vectors used in the solution algorithm.
    const unsigned int n_eig_request = 10;
    equation_systems.parameters.set<unsigned int>("eigenpairs")    = n_eig_request;
    equation_systems.parameters.set<unsigned int>("basis vectors") = n_eig_request*3;
    
    // Initialize the data structures for the equation system.
    equation_systems.init();
    
    // Prints information about the system to the screen.
    equation_systems.print_info();
    
    // Pass the Dirichlet dof IDs to the CondensedEigenSystem
    std::set<unsigned int> dirichlet_dof_ids;
    
    switch (dim-1) {
        case 1:
            eigen_system.attach_assemble_function (assemble_beam_matrices);
            get_beam_dirichlet_dofs(equation_systems, "Eigensystem", dirichlet_dof_ids);
            break;
            
        case 2:
            eigen_system.attach_assemble_function (assemble_plate_matrices);
            get_plate_dirichlet_dofs(equation_systems, "Eigensystem", dirichlet_dof_ids);
            break;

        default:
            libmesh_error();
            break;
    }

    eigen_system.initialize_condensed_dofs(dirichlet_dof_ids);
    
    // Solve the system "Eigensystem".
    eigen_system.solve();

    // Get the number of converged eigen pairs.
    unsigned int nconv = std::min(eigen_system.get_n_converged(),
                                  n_eig_request);

    std::ofstream file;
    file.open("modal_data.in", std::ofstream::out);

    file << "n_eig = " << nconv << std::endl;
    
    for (unsigned int i=0; i<nconv; i++)
    {
        std::pair<Real, Real> val = eigen_system.get_eigenpair(i);
        
        // also add the solution as an independent vector, which will have to
        // be read in
        std::ostringstream vec_name;
        vec_name << "mode_" << i;
        NumericVector<Real>& vec = eigen_system.add_vector(vec_name.str());
        vec = *eigen_system.solution;
        if (equation_systems.parameters.get<bool>("if_exchange_AB_matrices"))
        {
            file << "eig_"  << i << " = " << 1./val.first << std::endl;
            vec.scale(1./sqrt(val.first));
            vec.close();
        }
        else
            file << "eig_"  << i << " = " << val.first << std::endl;

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
    
    file.close();
    
    XdrIO xdr(mesh, true);
    xdr.write("saved_structural_mesh.xdr");
    equation_systems.write("saved_structural_solution.xdr",
                           libMeshEnums::ENCODE,
                           (EquationSystems::WRITE_DATA |
                            EquationSystems::WRITE_ADDITIONAL_DATA));

    // All done.
    return 0;
}

#endif // LIBMESH_USE_COMPLEX_NUMBERS
