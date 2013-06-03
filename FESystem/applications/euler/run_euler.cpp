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
#include "solvers/time_solvers/residual_based_adaptive_time_solver.h"

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
#include "libmesh/nemesis_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/euler_solver.h"
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
    const Real terminate_tolerance       = infile("pseudo_time_terminate_tolerance", 1.0e-5);
    unsigned int n_timesteps             = infile("n_timesteps", 1);
    const unsigned int write_interval    = infile("write_interval", 5);
    const bool if_use_amr                = infile("if_use_amr", false);
    const unsigned int max_adaptivesteps = infile("max_adaptivesteps", 0);
    const Real amr_threshold             = infile("amr_threshold", 1.0e1);
    const Real amr_time_shrink_factor    = infile("amr_time_shrink_factor", 0.25);
    const unsigned int n_uniform_refine  = infile("n_uniform_refine", 0);
    const unsigned int dim               = infile("dimension", 2);
    const bool if_panel_mesh             = infile("use_panel_mesh", true);
    
    // Skip higher-dimensional examples on a lower-dimensional libMesh build
    libmesh_example_assert(dim <= LIBMESH_DIM, "2D/3D support");
    
    // We have only defined 2 and 3 dimensional problems
    libmesh_assert (dim == 2 || dim == 3);
    
    // Create a mesh.
    SerialMesh mesh(init.comm());
    //ParallelMesh mesh(init.comm());

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
            
            // the mesh created by the MeshTool has uniform points. Iterate over all nodes
            // and find the corresponding point in the non-uniform distributed points. Use
            // that to assign the point location.
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

            // the mesh created by the MeshTool has uniform points. Iterate over all nodes
            // and find the corresponding point in the non-uniform distributed points. Use
            // that to assign the point location.
            
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
                bool side_on_panel = false, side_on_slip_wall=false;
                AutoPtr<Elem> side_elem ((*e_it)->side(i_side).release());
                for (unsigned int i_node=0; i_node<side_elem->n_nodes(); i_node++)
                {
                    const Node& n = *(side_elem->get_node(i_node));
                    if (dim == 2)
                    {
                        if ((n(1)==0.) && (n(0) >= x0))
                            side_on_panel = true;
                        else if ((n(1)==0.) && (n(0) > 0.) && (n(0) < x_length))
                            side_on_slip_wall = true;
                        else
                        {
                            side_on_panel = false;
                            side_on_slip_wall = false;
                        }
                    }
                    if (dim == 3)
                        if ((n(2)>0.) || (n(0) < x0) || (n(0)>x1) || (n(1) < y0) || (n(1)>y1))
                        {
                            side_on_panel = false;
                            break;
                        }
                }
                if (side_on_panel)
                    mesh.boundary_info->add_side(*e_it, 0, 10);
                if (side_on_slip_wall)
                    mesh.boundary_info->add_side(*e_it, 0, 11);
            }
        }
        
        
        MeshBase::node_iterator   n_it  = mesh.nodes_begin();
        const Mesh::node_iterator n_end = mesh.nodes_end();
        
        Real x_val, y_val, z_val;
        
        for (; n_it != n_end; n_it++)
        {
            Node& n =  **n_it;
            
            if (dim == 2)
            {
                // this is for sine bump
                if ((n(0) >= x0) && (n(0) <= x1))
                {
                    x_val = n(0);
                    y_val = n(1);
                    
                    n(1) += thickness*(1.0-y_val/y_length)*sin(pi*(x_val-x0)/chord);
                }

//                // this is for gaussian bump
//                x_val = n(0);
//                y_val = n(1);
//                
//                n(1) += (1.0-y_val/y_length) * 0.0625 * exp(-25.*pow(x_val-x_length*0.5,2.0));

                
//            // this is for ringleb problem
//            if (dim == 2)
//            {
//                x_val = n(0);
//                if (x_val <= 0.5)
//                    x_val *= 2.;
//                else
//                    x_val = (0.5-x_val)*2.+1.;
//                
//                double exp_val = 0.2, x_break=0.7;
//                if (x_val <= x_break)
//                    x_val = x_val/x_break*pow(x_break, exp_val);
//                else
//                    x_val = pow(x_val,exp_val);
//                
//                y_val = n(1);
//                double kval, qval, aval, gammaval =1.4, rhoval, pval, Jval;
//                kval = y_val * 0.8 + 0.7; // linear variation from 0.7 to 1.5
//                qval = x_val * (kval-0.5) + 0.5; // linear variation from 0.5 to kval
//                aval = sqrt(1.-0.5*(gammaval-1.)*qval*qval);
//                rhoval = pow(aval, 2./(gammaval-1.));
//                pval = pow(aval, 2.*gammaval/(gammaval-1.))/gammaval;
//                Jval = 1./aval + pow(aval,-3.)/3. + pow(aval,-5.)/5. - 0.5*log((1.+aval)/(1.-aval));
//                
//                n(1) = sqrt(1.-pow(qval/kval,2.))/kval/rhoval/qval;
//                if (n(0) > 0.5)
//                    n(1) *= -1.;
//                n(0) = (2./kval/kval - 1./qval/qval)/2./rhoval - Jval/2.;
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
        // ExodusII_IO gmsh_io(mesh);
        // Nemesis_IO gmsh_io(mesh);
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
    
    
    ResidualBaseAdaptiveTimeSolver *timesolver = new ResidualBaseAdaptiveTimeSolver(system);
    EulerSolver *core_time_solver = new EulerSolver(system);
    
    timesolver->quiet              = infile("timesolver_solver_quiet", true);
    timesolver->growth_exponent    = infile("timesolver_growth_exponent", 1.2);
    timesolver->n_iters_per_update = infile("timesolver_update_n_iters", 10);
    timesolver->min_deltat         = infile("timesolver_min_deltat", 1.0e-3);
    timesolver->max_growth         = infile("timesolver_maxgrowth", 4.0);
    timesolver->min_growth         = infile("timesolver_mingrowth", 0.25);
    timesolver->max_deltat         = infile("timesolver_max_deltat", 5.0e2);

    core_time_solver->theta        = infile("timesolver_theta", 1.0);
    
    timesolver->core_time_solver = AutoPtr<UnsteadySolver>(core_time_solver);
    system.time_solver = AutoPtr<UnsteadySolver>(timesolver);
    
    system.dc_recalculate_tolerance = infile("dc_recalculate_tolerance", 10.e-8);
    
    equation_systems.init ();
    
#else
    // Declare the system "EulerSystem"
    FrequencyDomainLinearizedEuler & system =
    equation_systems.add_system<FrequencyDomainLinearizedEuler> ("FrequencyDomainLinearizedEuler");
    
    FrequencyDomainFluidPostProcessSystem& fluid_post =
    equation_systems.add_system<FrequencyDomainFluidPostProcessSystem> ("DeltaFluidPostProcessSystem");
    
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
    //system.verify_analytic_jacobians = 1.0e-3;
    
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
    bool continue_iterations = true;
    unsigned int t_step=0, amr_steps = max_adaptivesteps;
    Real sol_norm = 1.0e10;
    if (!if_use_amr) amr_steps = 0;
    
    while (continue_iterations)
    {
        // A pretty update message
        std::cout << "\n\nSolving time step " << t_step << ", time = "
        << system.time << std::endl;
        
        // Adaptively solve the timestep
        unsigned int a_step = 0;
        if (if_use_amr && (amr_steps > 0) && (sol_norm < amr_threshold))
        {
            system.solve();
            system.print_integrated_lift_drag(std::cout);
            system.integrated_force->zero();
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
            delta_val_system.solution->close();
            delta_val_system.update();
#endif
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
                libmesh_assert_greater (nelem_target, 0);
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
            
            // decrement the amr counter
            amr_steps--;

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
            // tell the solver to recalculte the dc coeffs post refinement
            system.if_use_stored_dc_coeff = false;
#endif
        }

        // Do one last solve before the time step increment
        system.solve();
        system.print_integrated_lift_drag(std::cout);
        system.integrated_force->zero();
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
        system.system_comm->sum(system.entropy_error);
        system.system_comm->sum(system.total_volume);
        std::cout << system.entropy_error << " , " << system.total_volume << " , "
        << sqrt(system.entropy_error/system.total_volume) << " , " << std::endl;
        system.entropy_error = 0.; system.total_volume = 0.;
        system.evaluate_recalculate_dc_flag();
        delta_val_system.solution->close();
        delta_val_system.update();
#endif

        // Advance to the next timestep in a transient problem
        system.time_solver->advance_timestep();
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
        sol_norm = timesolver->_xdot_linf_approx;
#endif

        // check for termination criteria
        t_step++;

        libMesh::out << " Convergence monitor: L-infty norm: " << sol_norm << std::endl;
        if (((t_step >= n_timesteps) || (sol_norm < terminate_tolerance)) &&
            (amr_steps == 0))
        {
            libMesh::out << "\n === Terminating pseudo-time iterations ===" << std::endl;
            continue_iterations = false;
        }
        
        // Write out this timestep if we're requested to
        if (((t_step+1)%write_interval == 0) || !continue_iterations)
        {
            fluid_post.postprocess();

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
