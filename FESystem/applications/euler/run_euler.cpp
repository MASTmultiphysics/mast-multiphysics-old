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
#include "libmesh/serial_mesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/ensight_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/xdr_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/euler_solver.h"
#include "libmesh/twostep_time_solver.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/error_vector.h"
#include "libmesh/error_estimator.h"
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/uniform_refinement_estimator.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/string_to_enum.h"




void distributePoints(const unsigned int n_divs, const std::vector<double>& div_locations, const std::vector<unsigned int>& n_subdivs_in_div, const std::vector<double>& relative_mesh_size_at_div, std::vector<double>& points)
{
    libmesh_assert(div_locations.size() == n_divs+1);
    libmesh_assert(relative_mesh_size_at_div.size() == n_divs+1);
    libmesh_assert(n_subdivs_in_div.size() == n_divs);
    
    // calculate total number of points
    unsigned int n_total_points = 1;
    for (unsigned int i=0; i<n_divs; i++)
        n_total_points += n_subdivs_in_div[i];
    
    // resize the points vector and set the first and last points of each division
    points.resize(n_total_points);
    
    unsigned int n=1;
    points[0] = div_locations[0];
    for (unsigned int i=0; i<n_divs; i++)
    {
        n += n_subdivs_in_div[i];
        points[n-1] = div_locations[i+1];
    }
    
    n=1;
    double dx=0.0, growth_factor = 0.0;
    // now calculate the base mesh size, and calculate the nodal points
    for (unsigned int i=0; i<n_divs; i++)
    {
        growth_factor = pow(relative_mesh_size_at_div[i+1]/relative_mesh_size_at_div[i], 1.0/(n_subdivs_in_div[i]-1.0));
        if (fabs(growth_factor-1.0)>1.0e-10)
            dx = (div_locations[i+1]-div_locations[i]) * (1.0-growth_factor)/(1.0-pow(growth_factor, n_subdivs_in_div[i]));
        else
        {
            growth_factor = 1.0;
            dx = (div_locations[i+1]-div_locations[i]) / n_subdivs_in_div[i];
        }
        
        for (unsigned int n_pt=1; n_pt<n_subdivs_in_div[i]; n_pt++)
        {
            points[n+n_pt-1] = points[n+n_pt-2] + dx;
            dx *= growth_factor;
        }
        n += n_subdivs_in_div[i];
    }
}




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
    const Real deltat                    = infile("deltat", 0.005);
    unsigned int n_timesteps             = infile("n_timesteps", 20);
    unsigned int n_timesteps_delta_const = infile("n_timesteps_const_delta", 150);
    const unsigned int write_interval    = infile("write_interval", 5);
    const unsigned int max_adaptivesteps = infile("max_adaptivesteps", 0);
    const unsigned int amr_interval      = infile("amr_interval", 1);
    const Real amr_time_shrink_factor    = infile("amr_time_shrink_factor", 0.25);
    const unsigned int n_uniform_refine  = infile("n_uniform_refine", 0);
    const unsigned int dim               = infile("dimension", 2);
    const bool if_panel_mesh             = infile("use_panel_mesh", true);
    
    // Skip higher-dimensional examples on a lower-dimensional libMesh build
    libmesh_example_assert(dim <= LIBMESH_DIM, "2D/3D support");
    
    // We have only defined 2 and 3 dimensional problems
    libmesh_assert (dim == 2 || dim == 3);
    
    // Create a mesh.
    SerialMesh mesh;
   // ParallelMesh mesh;    

    // And an object to refine it
    MeshRefinement mesh_refinement(mesh);
    mesh_refinement.coarsen_by_parents() = true;
    mesh_refinement.absolute_global_tolerance() = global_tolerance;
    mesh_refinement.nelem_target() = nelem_target;
    mesh_refinement.refine_fraction() = infile("refine_fraction",0.80);
    mesh_refinement.coarsen_fraction() = infile("coarsen_fraction",0.20);
    mesh_refinement.coarsen_threshold() = infile("coarsen_threshold",0.30);
    mesh_refinement.max_h_level() = infile("max_h_level",5);
    std::string strategy = infile("refine_strategy", std::string("error_fraction")),
    error_norm = infile("error_norm", std::string("kelly")),
    elem_type = infile("elem_type", std::string("QUAD4"));
    
    Real mesh_dx, mesh_dy, mesh_dz;
    
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    
    if (if_panel_mesh)
    {
        // first calculate the distributed points
        std::vector<double> div_locations, relative_mesh_size_at_div, x_points, y_points, z_points;
        std::vector<unsigned int> n_subdivs_in_div;
        
        const Real pi = acos(-1.),
        x_length= infile("x_length", 10.),
        y_length= infile("y_length", 10.),
        z_length= infile("z_length", 10.),
        t_by_c =  infile("t_by_c", 0.05),
        chord =   infile("chord", 1.0),
        span =    infile("span", 1.0),
        dx_chordwise_inlet =  infile("dx_chordwise_inlet", 1.0), // these are relative sizes
        dx_chordwise_le =     infile("dx_chordwise_le", 1.0),
        dx_chordwise_te =     infile("dx_chordwise_te", 1.0),
        dx_chordwise_outlet = infile("dx_chordwise_outlet", 1.0),
        dx_spanwise_inlet =   infile("dx_spanwise_inlet", 1.0),
        dx_spanwise_le =      infile("dx_spanwise_le", 1.0),
        dx_spanwise_te =      infile("dx_spanwise_te", 1.0),
        dx_spanwise_outlet =  infile("dx_spanwise_outlet", 1.0),
        dx_vert_surface =     infile("dx_vert_surface", 1.0),
        dx_vert_inf =         infile("dx_vert_inf", 1.0),
        thickness = 0.5*t_by_c*chord,
        x0=x_length*0.5-chord*0.5, x1=x0+chord, y0=y_length*0.5-span*0.5, y1=y0+span ;
        
        const unsigned int
        n_chordwise_le_divs    = infile("n_chordwise_le_divs", 10),
        n_chordwise_panel_divs = infile("n_chordwise_panel_divs", 10),
        n_chordwise_te_divs    = infile("n_chordwise_te_divs", 10),
        n_spanwise_le_divs     = infile("n_spanwise_le_divs", 10),
        n_spanwise_panel_divs  = infile("n_spanwise_panel_divs", 10),
        n_spanwise_te_divs     = infile("n_spanwise_te_divs", 10),
        n_vertical_divs        = infile("n_vertical_divs", 10);
        
        // Use the MeshTools::Generation mesh generator to create a uniform
        // grid on the square [-1,1]^D.  We instruct the mesh generator
        // to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
        // elements in 3D.  Building these higher-order elements allows
        // us to use higher-order approximation, as in example 3.
        {   double vals[] = {0., (x_length-chord)/2., (x_length+chord)/2. , x_length};
            div_locations.assign(vals, vals+4); }
        {   unsigned int vals[] = {n_chordwise_le_divs, n_chordwise_panel_divs, n_chordwise_te_divs};
            n_subdivs_in_div.assign(vals, vals+3); }
        {   double vals[] = {dx_chordwise_inlet, dx_chordwise_le, dx_chordwise_te, dx_chordwise_outlet};
            relative_mesh_size_at_div.assign(vals, vals+4); }
        distributePoints(3, div_locations, n_subdivs_in_div, relative_mesh_size_at_div, x_points);
        mesh_dx = 1./(1.*(n_chordwise_le_divs+n_chordwise_panel_divs+n_chordwise_te_divs));
        
        if (dim == 2)
        {
            {   double vals[] = {0., y_length};
                div_locations.assign(vals, vals+2); }
            {   unsigned int vals[] = {n_vertical_divs};
                n_subdivs_in_div.assign(vals, vals+1); }
            {   double vals[] = {dx_vert_surface, dx_vert_inf};
                relative_mesh_size_at_div.assign(vals, vals+2); }
            distributePoints(1, div_locations, n_subdivs_in_div, relative_mesh_size_at_div, y_points);
            
            mesh_dy = 1./(1.*n_vertical_divs);
            
            MeshTools::Generation::build_square (mesh,
                                                 n_chordwise_le_divs+n_chordwise_panel_divs+n_chordwise_te_divs,
                                                 n_vertical_divs,
                                                 0., 1.,
                                                 0., 1.,
                                                 Utility::string_to_enum<ElemType>(elem_type));
            
            MeshBase::node_iterator   n_it  = mesh.nodes_begin();
            const Mesh::node_iterator n_end = mesh.nodes_end();
            Real x_val, y_val;
            unsigned int x_id, y_id;
            
            for (; n_it != n_end; n_it++)
            {
                Node& n =  **n_it;
                
                x_val = n(0);
                y_val = n(1);
                
                x_id = floor(x_val/mesh_dx);
                // find the correct id
                bool found_id = false;
                unsigned int correct_id = 0;
                for (unsigned int i=0; i<3; i++)
                    if (fabs((x_id+i)*mesh_dx-x_val) <= 1.0e-8)
                    {
                        found_id = true;
                        correct_id = x_id+i;
                        break;
                    }
                libmesh_assert(found_id);
                x_id = correct_id;
                
                found_id = false;
                correct_id = 0;
                y_id = floor(y_val/mesh_dy);
                for (unsigned int i=0; i<3; i++)
                    if (fabs((y_id+i)*mesh_dy-y_val) <= 1.0e-8)
                    {
                        found_id = true;
                        correct_id = y_id+i;
                        break;
                    }
                libmesh_assert(found_id);
                y_id = correct_id;
                
                n(0) = x_points[x_id];
                n(1) = y_points[y_id];
            }
            
        }
        
        else if (dim == 3)
        {
            {   double vals[] = {0., (y_length-span)/2., (y_length+span)/2. , y_length};
                div_locations.assign(vals, vals+4); }
            {   unsigned int vals[] = {n_spanwise_le_divs, n_spanwise_panel_divs, n_spanwise_te_divs};
                n_subdivs_in_div.assign(vals, vals+3); }
            {   double vals[] = {dx_spanwise_inlet, dx_spanwise_le, dx_spanwise_te, dx_spanwise_outlet};
                relative_mesh_size_at_div.assign(vals, vals+4); }
            distributePoints(3, div_locations, n_subdivs_in_div, relative_mesh_size_at_div, y_points);
            
            
            {   double vals[] = {0., z_length};
                div_locations.assign(vals, vals+2); }
            {   unsigned int vals[] = {n_vertical_divs};
                n_subdivs_in_div.assign(vals, vals+1); }
            {   double vals[] = {dx_vert_surface, dx_vert_inf};
                relative_mesh_size_at_div.assign(vals, vals+2); }
            distributePoints(1, div_locations, n_subdivs_in_div, relative_mesh_size_at_div, z_points);
            
            mesh_dy = 1./(1.*(n_spanwise_le_divs+n_spanwise_panel_divs+n_spanwise_te_divs));
            mesh_dz = 1./(1.*n_vertical_divs);
            
            MeshTools::Generation::build_cube (mesh,
                                               n_chordwise_le_divs+n_chordwise_panel_divs+n_chordwise_te_divs,
                                               n_spanwise_le_divs+n_spanwise_panel_divs+n_spanwise_te_divs,
                                               n_vertical_divs,
                                               0., 1.,
                                               0., 1.,
                                               0., 1.,
                                               Utility::string_to_enum<ElemType>(elem_type));
            
            MeshBase::node_iterator   n_it  = mesh.nodes_begin();
            const Mesh::node_iterator n_end = mesh.nodes_end();
            Real x_val, y_val, z_val;
            unsigned int x_id, y_id, z_id;
            
            for (; n_it != n_end; n_it++)
            {
                Node& n =  **n_it;
                x_val = n(0);
                y_val = n(1);
                z_val = n(2);
                
                x_id = floor(x_val/mesh_dx);
                y_id = floor(y_val/mesh_dy);
                z_id = floor(z_val/mesh_dz);
                
                // find the correct id
                bool found_id = false;
                unsigned int correct_id = 0;
                for (unsigned int i=0; i<3; i++)
                    if (fabs((x_id+i)*mesh_dx-x_val) <= 1.0e-8)
                    {
                        found_id = true;
                        correct_id = x_id+i;
                        break;
                    }
                libmesh_assert(found_id);
                x_id = correct_id;
                
                found_id = false;
                correct_id = 0;
                for (unsigned int i=0; i<3; i++)
                    if (fabs((y_id+i)*mesh_dy-y_val) <= 1.0e-8)
                    {
                        found_id = true;
                        correct_id = y_id+i;
                        break;
                    }
                libmesh_assert(found_id);
                y_id = correct_id;
                
                found_id = false;
                correct_id = 0;
                for (unsigned int i=0; i<3; i++)
                    if (fabs((z_id+i)*mesh_dz-z_val) <= 1.0e-8)
                    {
                        found_id = true;
                        correct_id = z_id+i;
                        break;
                    }
                libmesh_assert(found_id);
                z_id = correct_id;
                
                n(0) = x_points[x_id];
                n(1) = y_points[y_id];
                n(2) = z_points[z_id];
            }
            
        }
        
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
    }
#else
    
    mesh.read("saved_mesh.xdr");
    
#endif // LIBMESH_USE_COMPLEX_NUMBERS
    
    // uniformly refine the mesh
    for (unsigned int i=0; i<n_uniform_refine; i++)
        mesh_refinement.uniformly_refine();
    
    // Print information about the mesh to the screen.
    mesh.print_info();
    
    
    // Create an equation systems object.
    EquationSystems equation_systems (mesh);
    
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    // Declare the system "EulerSystem"
    EulerSystem & system =
    equation_systems.add_system<EulerSystem> ("EulerSystem");
    
    system.attach_init_function (init_euler_variables);
    
    FluidPostProcessSystem& fluid_post =
    equation_systems.add_system<FluidPostProcessSystem> ("FluidPostProcessSystem");
    System& delta_val_system =
    equation_systems.add_system<System> ("DeltaValSystem");
    delta_val_system.add_variable("delta", FEType(CONSTANT, MONOMIAL));
    
    
    TwostepTimeSolver *timesolver = new TwostepTimeSolver(system);
    
    //    timesolver->target_tolerance = infile("timesolver_tolerance",);
    //    timesolver->upper_tolerance  = infile("timesolver_upper_tolerance",);
    //    timesolver->component_norm   = SystemNorm(infile("timesolver_norm",));
    timesolver->quiet = infile("timesolver_solver_quiet", true);
    timesolver->min_deltat = infile("timesolver_min_deltat", 1.0e-3);
    timesolver->global_tolerance = infile("timesolver_global_tolerance", false);
    timesolver->max_growth       = infile("timesolver_maxgrowth", 1.5);
    timesolver->max_deltat       = infile("timesolver_max_deltat", 5.0e2);
    timesolver->core_time_solver = AutoPtr<EulerSolver>(new EulerSolver(system));
    system.time_solver = AutoPtr<UnsteadySolver>(timesolver);
    
    //system.time_solver = AutoPtr<UnsteadySolver>(new EulerSolver(system));
    
    equation_systems.init ();
    
#else
    // Declare the system "EulerSystem"
    FrequencyDomainLinearizedEuler & system =
    equation_systems.add_system<FrequencyDomainLinearizedEuler> ("FrequencyDomainLinearizedEuler");
    
    FluidPostProcessSystem& fluid_post =
    equation_systems.add_system<FluidPostProcessSystem> ("DeltaFluidPostProcessSystem");
    
    system.time_solver =
    AutoPtr<TimeSolver>(new SteadySolver(system));
    libmesh_assert_equal_to (n_timesteps, 1);
    
    equation_systems.read<Real>("saved_solution.xdr", libMeshEnums::DECODE);
    VTKIO(mesh).write_equation_systems("steady_solution.pvtu", equation_systems);
    
#endif
    
    
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
    solver.brent_line_search = true;
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
    
    // Print information about the system to the screen.
    equation_systems.print_info();
    
    
    // Now we begin the timestep loop to compute the time-accurate
    // solution of the equations.
    for (unsigned int t_step=0; t_step != n_timesteps; ++t_step)
    {
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
        if (t_step > n_timesteps_delta_const)
            system.if_use_stored_delta = true;
#endif
        // A pretty update message
        std::cout << "\n\nSolving time step " << t_step << ", time = "
        << system.time << std::endl;
        
        // Adaptively solve the timestep
        unsigned int a_step = 0;
        if ((t_step+1)%amr_interval == 0)
            for (; a_step != max_adaptivesteps; ++a_step)
            {
                system.solve();
                system.print_integrated_lift_drag(std::cout);
                system.integrated_force->zero();
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
                delta_val_system.solution->close();
                delta_val_system.update();
#endif
                fluid_post.postprocess();
                
                ErrorVector error;
                
                AutoPtr<ErrorEstimator> error_estimator;
                
                // To solve to a tolerance in this problem we
                // need a better estimator than Kelly
                if (error_norm == "uniform")
                {
                    // We can't adapt to both a tolerance and a mesh
                    // size at once
                    libmesh_assert_greater (global_tolerance, 0);
                    libmesh_assert_equal_to (nelem_target, 0);
                    
                    UniformRefinementEstimator *u = new UniformRefinementEstimator;
                    u->error_norm = L2;
                    error_estimator.reset(u);
                }
                else if (error_norm == "kelly")
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
                else if (error_norm == "patch")
                {
                    error_estimator.reset(new PatchRecoveryErrorEstimator);
                }
                else
                    libmesh_assert(false);
                
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
                
                if (strategy == "error_fraction")
                    mesh_refinement.flag_elements_by_error_fraction(error);
                else if (strategy == "error_tolerance")
                {
                    // If we've reached our desired tolerance, we
                    // don't need any more adaptive steps
                    if (global_error < global_tolerance)
                        break;
                    mesh_refinement.flag_elements_by_error_tolerance(error);
                }
                else if (strategy == "nelem_target")
                {
                    if (mesh_refinement.flag_elements_by_nelem_target(error))
                    {
                        mesh_refinement.refine_and_coarsen_elements();
                        equation_systems.reinit();
                        a_step = max_adaptivesteps;
                        break;
                    }
                }
                else if (strategy == "elem_fraction")
                    mesh_refinement.flag_elements_by_elem_fraction(error);
                else if (strategy == "mean_stddev")
                    mesh_refinement.flag_elements_by_mean_stddev(error);
                else
                    libmesh_assert(false);
                
                // Carry out the adaptive mesh refinement/coarsening
                mesh_refinement.refine_and_coarsen_elements();
                equation_systems.reinit();
                
                // reduce the time step size by a factor
                system.deltat *= amr_time_shrink_factor;
                
                std::cout << "Refined mesh to "
                << mesh.n_active_elem()
                << " active elements and "
                << equation_systems.n_active_dofs()
                << " active dofs." << std::endl;
            }
        // Do one last solve if necessary
       if ((a_step == max_adaptivesteps) || ((t_step+1)%amr_interval != 0))
        {
            system.solve();
            system.print_integrated_lift_drag(std::cout);
            system.integrated_force->zero();
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
            delta_val_system.solution->close();
            delta_val_system.update();
#endif
            fluid_post.postprocess();
        }
        
        // Advance to the next timestep in a transient problem
        system.time_solver->advance_timestep();
        
        // Write out this timestep if we're requested to
        if ((t_step+1)%write_interval == 0)
        {
            std::ostringstream file_name, b_file_name;
            
            // We write the file in the ExodusII format.
            file_name << "out_"
            << std::setw(3)
            << std::setfill('0')
            << std::right
            << t_step+1
            << ".pvtu";
            
            VTKIO(mesh).write_equation_systems(file_name.str(),
                                               equation_systems);

            b_file_name << "b_out_"
            << std::setw(3)
            << std::setfill('0')
            << std::right
            << t_step+1
            << ".pvtu";
            
            std::set<unsigned int> bc_ids; bc_ids.insert(0);
            VTKIO(mesh, true).write_equation_systems(b_file_name.str(),
                                               equation_systems);
            
            //            ExodusII_IO(mesh).write_timestep(file_name.str(),
            //                                             equation_systems,
            //                                             1, /* This number indicates how many time steps
            //                                                 are being written to the file */
            //                                             system.time);
        }
    }
    
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    XdrIO xdr(mesh, true);
    xdr.write("saved_mesh.xdr");
    equation_systems.write("saved_solution.xdr", libMeshEnums::ENCODE);
#endif
    
    // All done.
    return 0;
}
