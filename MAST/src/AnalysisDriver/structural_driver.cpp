//
//  run_euler.cpp
//  FESystem
//
//  Created by Manav Bhatia on 3/6/13.
//
//

// C++ includes
#include <iomanip>

// MAST includes
#include "StructuralElems/structural_system_assembly.h"
#include "PropertyCards/element_property_card_base.h"
#include "PropertyCards/element_property_card_3D.h"
#include "PropertyCards/element_property_card_2D.h"
#include "PropertyCards/material_property_card_base.h"

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

#ifndef LIBMESH_USE_COMPLEX_NUMBERS


// The main program.
int structural_driver (LibMeshInit& init, GetPot& infile,
                       int argc, char* const argv[])
{
    SerialMesh mesh(init.comm());
    mesh.set_mesh_dimension(2);

    MeshTools::Generation::build_square (mesh, 10, 10, 0., 1., 0., 1., QUAD9);
    //MeshTools::Generation::build_cube (mesh, 5, 5, 5, 0., 1., 0., 1., 0., 1., HEX8);

    mesh.prepare_for_use();
    
    // Print information about the mesh to the screen.
    mesh.print_info();
    
    
    // Create an equation systems object.
    EquationSystems equation_systems (mesh);
    equation_systems.parameters.set<GetPot*>("input_file") = &infile;
    
    // Declare the system
//    NonlinearImplicitSystem & system =
//    equation_systems.add_system<NonlinearImplicitSystem> ("StructuralSystem");
    CondensedEigenSystem & system =
    equation_systems.add_system<CondensedEigenSystem> ("StructuralSystem");
    
    unsigned int o = infile("fe_order", 1);
    std::string fe_family = infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    system.add_variable ( "ux", static_cast<Order>(o), fefamily);
    system.add_variable ( "uy", static_cast<Order>(o), fefamily);
    system.add_variable ( "uz", static_cast<Order>(o), fefamily);
    system.add_variable ( "tx", static_cast<Order>(o), fefamily);
    system.add_variable ( "ty", static_cast<Order>(o), fefamily);
    system.add_variable ( "tz", static_cast<Order>(o), fefamily);

    //system.time_solver = AutoPtr<TimeSolver>(new SteadySolver(system));
    
    // Set the type of the problem, here we deal with
    // a generalized Hermitian problem.
    system.set_eigenproblem_type(GHEP);
    system.extra_quadrature_order = 0;
    
    // Order the eigenvalues "smallest first"
    system.eigen_solver->set_position_of_spectrum(LARGEST_MAGNITUDE);

    
    equation_systems.init ();

//    system.time_solver->diff_solver()->quiet = false;
//    system.time_solver->diff_solver()->verbose = true;
//    system.time_solver->diff_solver()->relative_residual_tolerance =  1.0e-8;
//    system.time_solver->diff_solver()->absolute_residual_tolerance =  1.0e-8;
    
    // Print information about the system to the screen.
    equation_systems.print_info();
    
    MAST::IsotropicMaterialPropertyCard mat;
    MAST::ElementPropertyCard3D prop3d;
    MAST::Solid2DSectionElementPropertyCard prop2d;
    
    MAST::FunctionValue<Real>& E = mat.add<Real>("E", MAST::CONSTANT_FUNCTION),
    &nu = mat.add<Real>("nu", MAST::CONSTANT_FUNCTION),
    &rho = mat.add<Real>("rho", MAST::CONSTANT_FUNCTION),
    &kappa = mat.add<Real>("kappa", MAST::CONSTANT_FUNCTION),
    &h =  prop2d.add<Real>("h", MAST::CONSTANT_FUNCTION);
    E  = 72.e9;
    nu = 0.33;
    rho = 2700.;
    kappa = 5./6.;
    h  = 0.002;
    
    prop3d.set_material(mat);
    prop2d.set_material(mat);
    prop2d.set_diagonal_mass_matrix(false);
    
//    MAST::StructuralSystemAssembly structural_assembly(system,
//                                                       MAST::STATIC,
//                                                       infile);
    MAST::StructuralSystemAssembly structural_assembly(system,
                                                       MAST::MODAL,
                                                       infile);
    structural_assembly.set_property_for_all_elems(prop2d);
    
    // Set the number of requested eigenpairs \p n_evals and the number
    // of basis vectors used in the solution algorithm.
    const unsigned int n_eig_request = 10;
    equation_systems.parameters.set<bool>("if_exchange_AB_matrices") = true;
    equation_systems.parameters.set<unsigned int>("eigenpairs")    = n_eig_request;
    equation_systems.parameters.set<unsigned int>("basis vectors") = n_eig_request*3;

    // Pass the Dirichlet dof IDs to the CondensedEigenSystem
    std::set<unsigned int> dirichlet_dof_ids;
    
    system.attach_assemble_object(structural_assembly);
    structural_assembly.get_dirichlet_dofs(dirichlet_dof_ids);
    system.initialize_condensed_dofs(dirichlet_dof_ids);

    system.solve();

    if (false) {
        
        // We write the file in the ExodusII format.
        Nemesis_IO(mesh).write_equation_systems("out.exo",
                                                equation_systems);
    }
    else {
        
        // Get the number of converged eigen pairs.
        unsigned int nconv = std::min(system.get_n_converged(),
                                      n_eig_request);
        
        std::ofstream file;
        file.open("modal_data.in", std::ofstream::out);
        
        file << "n_eig = " << nconv << std::endl;
        
        for (unsigned int i=0; i<nconv; i++)
        {
            std::pair<Real, Real> val = system.get_eigenpair(i);
            
            // also add the solution as an independent vector, which will have to
            // be read in
            std::ostringstream vec_name;
            vec_name << "mode_" << i;
            NumericVector<Real>& vec = system.add_vector(vec_name.str());
            vec = *system.solution;
            std::complex<Real> eigval;
            if (equation_systems.parameters.get<bool>("if_exchange_AB_matrices"))
            {
                file << "eig_"  << i << " = " << 1./val.first << std::endl;
                vec.scale(1./sqrt(val.first));
                vec.close();
                
                // now write the eigenvalues
                eigval = std::complex<Real>(val.first, val.second);
                eigval = 1./eigval;

                std::cout << std::setw(5) << i
                << std::setw(10) << eigval.real()
                << " + i  "
                << std::setw(10) << eigval.imag() << std::endl;

            }
            else {
                file << "eig_"  << i << " = " << val.first << std::endl;
                
                std::cout << std::setw(5) << i
                << std::setw(10) << val.first
                << " + i  "
                << std::setw(10) << val.second << std::endl;
            }
            
            
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
    }
    
    // now write the data to an output file
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
