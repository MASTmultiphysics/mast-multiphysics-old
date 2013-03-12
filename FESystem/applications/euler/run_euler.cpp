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
#include "euler/assembleEuler.h"
#include "euler/frequency_domain_linearized_euler.h"

// libmesh includes
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/euler_solver.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/error_vector.h"
#include "libmesh/error_estimator.h"
#include "libmesh/uniform_refinement_estimator.h"
#include "libmesh/kelly_error_estimator.h"





// The main program.
int main (int argc, char* const argv[])
{
    // Initialize libMesh.
    LibMeshInit init (argc, argv);
    
    // Parse the input file
    GetPot infile("system_input.in");
    
    // Read in parameters from the input file
    const Real global_tolerance          = infile("global_tolerance", 0.);
    const unsigned int nelem_target      = infile("n_elements", 400);
    const bool transient                 = infile("transient", true);
    const Real deltat                    = infile("deltat", 0.005);
    unsigned int n_timesteps             = infile("n_timesteps", 20);
    const unsigned int write_interval    = infile("write_interval", 5);
    const unsigned int coarsegridsize    = infile("coarsegridsize", 15);
    const unsigned int coarserefinements = infile("coarserefinements", 0);
    const unsigned int max_adaptivesteps = infile("max_adaptivesteps", 0);
    const unsigned int dim               = infile("dimension", 2);
    const bool if_panel_mesh             = infile("use_panel_mesh", true);
    
    // Skip higher-dimensional examples on a lower-dimensional libMesh build
    libmesh_example_assert(dim <= LIBMESH_DIM, "2D/3D support");
    
    // We have only defined 2 and 3 dimensional problems
    libmesh_assert (dim == 2 || dim == 3);
    
    // Create a mesh.
    Mesh mesh;
    
    // And an object to refine it
    MeshRefinement mesh_refinement(mesh);
    mesh_refinement.coarsen_by_parents() = true;
    mesh_refinement.absolute_global_tolerance() = global_tolerance;
    mesh_refinement.nelem_target() = nelem_target;
    mesh_refinement.refine_fraction() = 0.15;
    mesh_refinement.coarsen_fraction() = 0.15;
    mesh_refinement.coarsen_threshold() = 0.1;
    

    if (if_panel_mesh)
    {
        const Real pi = acos(-1.),
        x_length= infile("x_length", 10.),
        y_length= infile("y_length", 10.),
        z_length= infile("y_length", 10.),
        t_by_c =  infile("t_by_c", 0.05),
        chord =   infile("chord", 1.0),
        span =    infile("span", 1.0),
        thickness = 0.5*t_by_c*chord,
        x0=x_length*0.5-chord*0.5, x1=x0+chord, y0=y_length*0.5-span*0.5, y1=y0+span ;
        
        // Use the MeshTools::Generation mesh generator to create a uniform
        // grid on the square [-1,1]^D.  We instruct the mesh generator
        // to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
        // elements in 3D.  Building these higher-order elements allows
        // us to use higher-order approximation, as in example 3.
        if (dim == 2)
            MeshTools::Generation::build_square (mesh,
                                                 coarsegridsize,
                                                 coarsegridsize,
                                                 0., x_length,
                                                 0., y_length,
                                                 QUAD4);
        else if (dim == 3)
            MeshTools::Generation::build_cube (mesh,
                                               coarsegridsize,
                                               coarsegridsize,
                                               coarsegridsize,
                                               0., x_length,
                                               0., y_length,
                                               0., z_length,
                                               HEX8);
        
        
        mesh_refinement.uniformly_refine(coarserefinements);
                
        //march over all the elmeents and tag the sides that all lie on the panel suface
        MeshBase::element_iterator e_it = mesh.elements_begin();
        const MeshBase::element_iterator e_end = mesh.elements_end();
        
        for ( ; e_it != e_end; e_it++)
        {
            for (unsigned int i_side=0; i_side<(*e_it)->n_sides(); i_side++)
            {
                bool side_on_panel = true;
                AutoPtr<Elem> side_elem ((*e_it)->side(i_side).release());
                for (unsigned int i_node=0; i_node<side_elem->n_nodes(); i_node++)
                {
                    const Node& n = *(side_elem->get_node(i_node));
                    if (dim == 2)
                        if ((n(1)>0.) || (n(0) < x0) || (n(0)>x1))
                        {
                            side_on_panel = false;
                            break;
                        }
                    if (dim == 3)
                        if ((n(2)>0.) || (n(0) < x0) || (n(0)>x1) || (n(1) < y0) || (n(1)>y1))
                        {
                            side_on_panel = false;
                            break;
                        }
                }
                if (side_on_panel)
                    mesh.boundary_info->add_side(*e_it, i_side, 10);
            }
        }
        
        
        MeshBase::node_iterator   n_it  = mesh.nodes_begin();
        const Mesh::node_iterator n_end = mesh.nodes_end();
        
        Real x_val, y_val, z_val;
        
        for (; n_it != n_end; n_it++)
        {
            Node& n =  **n_it;
            
            if (dim == 2)
                if ((n(0) >= x0) && (n(0) <= x1))
                {
                    x_val = n(0);
                    y_val = n(1);
                    
                    n(1) += thickness*(1.0-y_val/y_length)*sin(pi*(x_val-x0)/chord);
                }
            
            if (dim == 3)
                if ((n(0) >= x0) && (n(0) <= x1) &&
                    (n(1) >= y0) && (n(1) <= y1))
                {
                    x_val = n(0);
                    y_val = n(1);
                    z_val = n(2);
                    
                    n(2) += thickness*(1.0-z_val/z_length)*sin(pi*(x_val-x0)/chord)*sin(pi*(y_val-y0)/span);
                }
        }
    }
    else
    {
        mesh.set_mesh_dimension(dim);
        const std::string gmsh_input_file = infile("gmsh_input", std::string("mesh.msh"));
        GmshIO gmsh_io(mesh);
        gmsh_io.read(gmsh_input_file);
        mesh.prepare_for_use();
        mesh_refinement.uniformly_refine(coarserefinements);
        
        // Print information about the mesh to the screen.
    }
    

    // Print information about the mesh to the screen.
    mesh.print_info();

    // Create an equation systems object.
    EquationSystems equation_systems (mesh);
    
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    // Declare the system "EulerSystem"
    EulerSystem & system =
    equation_systems.add_system<EulerSystem> ("EulerSystem");
    
    system.attach_init_function (init_euler_variables);

    system.time_solver =
    AutoPtr<TimeSolver>(new EulerSolver(system));
#else
    // Declare the system "EulerSystem"
    FrequencyDomainLinearizedEuler & system =
    equation_systems.add_system<FrequencyDomainLinearizedEuler> ("FrequencyDomainLinearizedEuler");

    system.time_solver =
    AutoPtr<TimeSolver>(new SteadySolver(system));
    libmesh_assert_equal_to (n_timesteps, 1);
#endif
    FluidPostProcessSystem& fluid_post =
    equation_systems.add_system<FluidPostProcessSystem> ("FluidPostProcessSystem");
    
    // Initialize the system
    equation_systems.init ();
    
    system.print_residual_norms = infile("print_residual_norms", false);
    system.print_residuals = infile("print_residuals", false);
    system.print_jacobian_norms = infile("print_jacobian_norms", false);
    system.print_jacobians = infile("print_jacobians", false);
    
    // Set the time stepping options
    system.deltat = deltat;
    
    // And the nonlinear solver options
    NewtonSolver &solver = dynamic_cast<NewtonSolver&>(*(system.time_solver->diff_solver().get()));
    solver.quiet = infile("solver_quiet", true);
    solver.verbose = !solver.quiet;
    solver.max_nonlinear_iterations =
    infile("max_nonlinear_iterations", 15);
    solver.relative_step_tolerance =
    infile("relative_step_tolerance", 1.e-3);
    solver.relative_residual_tolerance =
    infile("relative_residual_tolerance", 0.0);
    solver.absolute_residual_tolerance =
    infile("absolute_residual_tolerance", 0.0);
    solver.continue_after_backtrack_failure =
    infile("continue_after_backtrack_failure", false);
    solver.continue_after_max_iterations =
    infile("continue_after_max_iterations", false);
    solver.require_residual_reduction =
    infile("require_residual_reduction", true);
    
    // And the linear solver options
    solver.max_linear_iterations =
    infile("max_linear_iterations", 50000);
    solver.initial_linear_tolerance =
    infile("initial_linear_tolerance", 1.e-3);
    //solver.brent_line_search = false;
    
    // Print information about the system to the screen.
    equation_systems.print_info();
    
    
    // Now we begin the timestep loop to compute the time-accurate
    // solution of the equations.
    for (unsigned int t_step=0; t_step != n_timesteps; ++t_step)
    {
        // A pretty update message
        std::cout << "\n\nSolving time step " << t_step << ", time = "
        << system.time << std::endl;
        
        // Adaptively solve the timestep
        unsigned int a_step = 0;
        for (; a_step != max_adaptivesteps; ++a_step)
        {
            system.solve();
            
            fluid_post.postprocess();
            
            ErrorVector error;
            
            AutoPtr<ErrorEstimator> error_estimator;
            
            // To solve to a tolerance in this problem we
            // need a better estimator than Kelly
            if (global_tolerance != 0.)
            {
                // We can't adapt to both a tolerance and a mesh
                // size at once
                libmesh_assert_equal_to (nelem_target, 0);
                
                UniformRefinementEstimator *u =
                new UniformRefinementEstimator;
                
                // The lid-driven cavity problem isn't in H1, so
                // lets estimate L2 error
                u->error_norm = L2;
                
                error_estimator.reset(u);
            }
            else
            {
                // If we aren't adapting to a tolerance we need a
                // target mesh size
                libmesh_assert_greater (nelem_target, 0);
                
                // Kelly is a lousy estimator to use for a problem
                // not in H1 - if we were doing more than a few
                // timesteps we'd need to turn off or limit the
                // maximum level of our adaptivity eventually
                error_estimator.reset(new KellyErrorEstimator);
            }
            
            // Calculate error based on u and v (and w?) but not p
            std::vector<Real> weights(dim+2,0.0);  // all set to 1.0
            weights[0] = 1.0;
            // Keep the same default norm type.
            std::vector<FEMNormType>
            norms(1, error_estimator->error_norm.type(0));
            error_estimator->error_norm = SystemNorm(norms, weights);
            
            error_estimator->estimate_error(system, error);
            
            // Print out status at each adaptive step.
            Real global_error = error.l2_norm();
            std::cout << "Adaptive step " << a_step << ": " << std::endl;
            if (global_tolerance != 0.)
                std::cout << "Global_error = " << global_error
                << std::endl;
            if (global_tolerance != 0.)
                std::cout << "Worst element error = " << error.maximum()
                << ", mean = " << error.mean() << std::endl;
            
            if (global_tolerance != 0.)
            {
                // If we've reached our desired tolerance, we
                // don't need any more adaptive steps
                if (global_error < global_tolerance)
                    break;
                mesh_refinement.flag_elements_by_error_tolerance(error);
            }
            else
            {
                // If flag_elements_by_nelem_target returns true, this
                // should be our last adaptive step.
                if (mesh_refinement.flag_elements_by_nelem_target(error))
                {
                    mesh_refinement.refine_and_coarsen_elements();
                    equation_systems.reinit();
                    a_step = max_adaptivesteps;
                    break;
                }
            }
            
            // Carry out the adaptive mesh refinement/coarsening
            mesh_refinement.refine_and_coarsen_elements();
            equation_systems.reinit();
            
            std::cout << "Refined mesh to "
            << mesh.n_active_elem()
            << " active elements and "
            << equation_systems.n_active_dofs()
            << " active dofs." << std::endl;
        }
        // Do one last solve if necessary
        if (a_step == max_adaptivesteps)
        {
            system.solve();
            
            fluid_post.postprocess();
        }
        
        // Advance to the next timestep in a transient problem
        system.time_solver->advance_timestep();
        
#ifdef LIBMESH_HAVE_EXODUS_API
        // Write out this timestep if we're requested to
        if ((t_step+1)%write_interval == 0)
        {
            std::ostringstream file_name;
            
            // We write the file in the ExodusII format.
            file_name << "out_"
            << std::setw(3)
            << std::setfill('0')
            << std::right
            << t_step+1
            << ".pvtu";
            
            VTKIO(mesh).write_equation_systems(file_name.str(),
                                               equation_systems);
//            ExodusII_IO(mesh).write_timestep(file_name.str(),
//                                             equation_systems,
//                                             1, /* This number indicates how many time steps
//                                                 are being written to the file */
//                                             system.time);
        }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
    }
    
    // All done.
    return 0;
}
