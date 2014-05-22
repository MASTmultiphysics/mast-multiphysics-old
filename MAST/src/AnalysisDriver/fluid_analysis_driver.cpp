//
//  run_euler.cpp
//  MAST
//
//  Created by Manav Bhatia on 3/6/13.
//
//

// C++ includes
#include <iomanip>

// MAST includes
#include "FluidElems/fluid_system.h"
#include "FluidElems/shock_tube_fluid_elem.h"
#include "FluidElems/frequency_domain_linearized_fluid_system.h"
#include "Solvers/residual_based_adaptive_time_solver.h"
#include "FluidElems/aerodynamic_qoi.h"
#include "FluidElems/fluid_newton_solver.h"
#include "BoundaryConditions/rigid_surface_motion.h"
#include "BoundaryConditions/flexible_surface_motion.h"
#include "BoundaryConditions/function_surface_motion.h"
#include "Mesh/panel_mesh.h"
#include "Mesh/gaussian_bump_mesh.h"
#include "Mesh/ringleb_mesh.h"


// libmesh includes
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/xdr_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/euler2_solver.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/error_vector.h"
#include "libmesh/error_estimator.h"
#include "libmesh/patch_recovery_error_estimator.h"
#include "libmesh/uniform_refinement_estimator.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/condensed_eigen_system.h"


#include "Numerics/fem_operator_matrix.h"

int main_fem_operator (int argc, char* const argv[])
{
    libMesh::DenseVector<libMesh::Real> vec1, vec2, shp;
    libMesh::DenseMatrix<libMesh::Real> mat1, mat2, bmat1, bmat2, tmp;
    
    vec1.resize(12); shp.resize(4);
    
    for (unsigned int i=0; i<shp.size(); i++)
        shp(i) = i+1;

    for (unsigned int i=0; i<vec1.size(); i++)
        vec1(i) = i+100;
    
    shp.print(std::cout);
    vec1.print(std::cout);
    mat1.print(std::cout);

    FEMOperatorMatrix b1, b2;
    b1.reinit(3, shp); b2.reinit(3, shp);
    bmat1.resize(3, shp.size()*3); bmat2.resize(3, shp.size()*3);
    for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<shp.size(); j++) {
            bmat1(i, i*shp.size()+j) = shp(j);
            bmat2(i, i*shp.size()+j) = shp(j);
        }

    
    vec2.resize(3);
    b1.vector_mult(vec2, vec1);
    std::cout << "FEMOperator: " << std::endl; vec2.print(std::cout);
    bmat1.vector_mult(vec2, vec1);
    std::cout << "Matrix:      " << std::endl; vec2.print(std::cout);
    
    vec1.resize(3);
    for (unsigned int i=0; i<vec1.size(); i++)
        vec1(i) = 100+i;
    vec2.resize(12);
    b1.vector_mult_transpose(vec2, vec1);
    std::cout << "FEMOperator: " << std::endl; vec2.print(std::cout);
    bmat1.vector_mult_transpose(vec2, vec1);
    std::cout << "Matrix:      " << std::endl; vec2.print(std::cout);
    

    mat1.resize(5,3);
    mat2.resize(5, 12);
    for (unsigned int i=0; i<5; i++)
        for (unsigned int j=0; j<3; j++)
            mat1(i,j) = (i+1)*(j+1);
    std::cout << "multiply matrix" << std::endl;  mat1.print();
    b1.left_multiply(mat2, mat1);
    std::cout << "FEMOperator: " << std::endl;  mat2.print();
    tmp = bmat1;
    tmp.left_multiply(mat1);
    std::cout << "Matrix:      " << std::endl;  tmp.print();


    mat1.resize(5,12);
    mat2.resize(5, 3);
    for (unsigned int i=0; i<5; i++)
        for (unsigned int j=0; j<12; j++)
            mat1(i,j) = (i+1)*(j+1);
    std::cout << "multiply matrix" << std::endl;  mat1.print();
    b1.left_multiply_transpose(mat2, mat1);
    std::cout << "FEMOperator: " << std::endl;  mat2.print();
    bmat1.get_transpose(tmp);
    tmp.left_multiply(mat1);
    std::cout << "Matrix:      " << std::endl;  tmp.print();
    
    
    mat1.resize(12,5);
    mat2.resize(3, 5);
    for (unsigned int i=0; i<12; i++)
        for (unsigned int j=0; j<5; j++)
            mat1(i,j) = (i+1)*(j+1);
    std::cout << "multiply matrix" << std::endl;  mat1.print();
    b1.right_multiply(mat2, mat1);
    std::cout << "FEMOperator: " << std::endl;  mat2.print();
    tmp = bmat1;
    tmp.right_multiply(mat1);
    std::cout << "Matrix:      " << std::endl;  tmp.print();
    
    
    mat1.resize(3, 5);
    mat2.resize(12, 5);
    for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<5; j++)
            mat1(i,j) = (i+1)*(j+1);
    std::cout << "multiply matrix" << std::endl;  mat1.print();
    b1.right_multiply_transpose(mat2, mat1);
    std::cout << "FEMOperator: " << std::endl;  mat2.print();
    bmat1.get_transpose(tmp);
    tmp.right_multiply(mat1);
    std::cout << "Matrix:      " << std::endl;  tmp.print();


    mat2.resize(12, 12);
    b1.right_multiply_transpose(mat2, b1);
    std::cout << "FEMOperator: " << std::endl;  mat2.print();
    bmat1.get_transpose(tmp);
    tmp.right_multiply(bmat1);
    std::cout << "Matrix:      " << std::endl;  tmp.print();
    
}

#include <unistd.h>

// The main program.
int fluid_driver (libMesh::LibMeshInit& init, GetPot& infile,
                  int argc, char* const argv[])
{
    // Read in parameters from the input file
    const libMesh::Real global_tolerance          = infile("global_tolerance", 0.);
    const unsigned int nelem_target      = infile("n_elements", 400);
    const libMesh::Real deltat                    = infile("deltat", 0.005);
    const libMesh::Real terminate_tolerance       = infile("pseudo_time_terminate_tolerance", 1.0e-5);
    unsigned int n_timesteps             = infile("n_timesteps", 1);
    const unsigned int write_interval    = infile("write_interval", 5);
    const bool if_use_amr                = infile("if_use_amr", false);
    const unsigned int max_adaptivesteps = infile("max_adaptivesteps", 0);
    const libMesh::Real amr_threshold             = infile("amr_threshold", 1.0e1);
    const libMesh::Real amr_time_shrink_factor    = infile("amr_time_shrink_factor", 0.25);
    const unsigned int n_uniform_refine  = infile("n_uniform_refine", 0);
    const unsigned int dim               = infile("dimension", 2);
    const bool if_panel_mesh             = infile("use_panel_mesh", true);
    
    // Create a mesh.
    //SerialMesh mesh(init.comm());
    ParallelMesh mesh(init.comm());
    
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
    error_norm = infile("error_norm", std::string("kelly"));
    
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    
    if (if_panel_mesh)
    {
        const unsigned int nx_divs = infile("nx_divs",0),
        ny_divs = infile("ny_divs",0),
        nz_divs = infile("nz_divs",0);
        const libMesh::Real t_by_c =  infile("t_by_c", 0.0);
        ElemType elem_type =
        Utility::string_to_enum<ElemType>(infile("elem_type", "QUAD4"));

        std::vector<libMesh::Real> x_div_loc(nx_divs+1), x_relative_dx(nx_divs+1),
        y_div_loc(ny_divs+1), y_relative_dx(ny_divs+1),
        z_div_loc(nz_divs+1), z_relative_dx(nz_divs+1);
        std::vector<unsigned int> x_divs(nx_divs), y_divs(ny_divs), z_divs(nz_divs);
        std::auto_ptr<MeshInitializer::CoordinateDivisions>
        x_coord_divs (new MeshInitializer::CoordinateDivisions),
        y_coord_divs (new MeshInitializer::CoordinateDivisions),
        z_coord_divs (new MeshInitializer::CoordinateDivisions);
        std::vector<MeshInitializer::CoordinateDivisions*> divs(dim);
        
        // now read in the values: x-coord
        if (nx_divs > 0)
        {
            for (unsigned int i_div=0; i_div<nx_divs+1; i_div++)
            {
                x_div_loc[i_div]     = infile("x_div_loc", 0., i_div);
                x_relative_dx[i_div] = infile( "x_rel_dx", 0., i_div);
                if (i_div < nx_divs) //  this is only till nx_divs
                    x_divs[i_div] = infile( "x_div_nelem", 0, i_div);
            }
            divs[0] = x_coord_divs.get();
            x_coord_divs->init(nx_divs, x_div_loc, x_relative_dx, x_divs);
        }

        // now read in the values: y-coord
        if ((dim > 1) && (ny_divs > 0))
        {
            for (unsigned int i_div=0; i_div<ny_divs+1; i_div++)
            {
                y_div_loc[i_div]     = infile("y_div_loc", 0., i_div);
                y_relative_dx[i_div] = infile( "y_rel_dx", 0., i_div);
                if (i_div < ny_divs) //  this is only till ny_divs
                    y_divs[i_div] = infile( "y_div_nelem", 0, i_div);
            }
            divs[1] = y_coord_divs.get();
            y_coord_divs->init(ny_divs, y_div_loc, y_relative_dx, y_divs);
        }

        // now read in the values: z-coord
        if ((dim == 3) && (nz_divs > 0))
        {
            for (unsigned int i_div=0; i_div<nz_divs+1; i_div++)
            {
                z_div_loc[i_div]     = infile("z_div_loc", 0., i_div);
                z_relative_dx[i_div] = infile( "z_rel_dx", 0., i_div);
                if (i_div < nz_divs) //  this is only till nz_divs
                    z_divs[i_div] = infile( "z_div_nelem", 0, i_div);
            }
            divs[2] = z_coord_divs.get();
            z_coord_divs->init(nz_divs, z_div_loc, z_relative_dx, z_divs);
        }
        
        const std::string mesh_type = infile("mesh_type", std::string(""));
        if (mesh_type == "panel")
        {
            const bool if_cos_bump = infile("if_cos_bump", false);
            const unsigned int n_max_bumps_x = infile("n_max_bumps_x", 1),
            n_max_bumps_y = infile("n_max_bumps_x", 1),
            panel_bc_id = infile("panel_bc_id", 10),
            symmetry_bc_id = infile("symmetry_bc_id", 11);
            
            if (dim == 1)
                MeshInitializer().init(divs, mesh, elem_type);
            else if (dim == 2)
                PanelMesh2D().init(t_by_c, if_cos_bump, n_max_bumps_x,
                                   panel_bc_id, symmetry_bc_id,
                                    divs, mesh, elem_type);
            else if (dim == 3)
                PanelMesh3D().init(t_by_c, if_cos_bump,
                                   n_max_bumps_x, n_max_bumps_y,
                                   panel_bc_id, symmetry_bc_id,
                                   divs, mesh, elem_type);
            else
                libmesh_error();
        }
        else if (mesh_type == "gaussian")
        {
            if (dim == 2)
                GaussianBumpMesh2D().init(t_by_c, divs, mesh, elem_type);
            else if (dim == 3)
                GaussianBumpMesh3D().init(t_by_c, divs, mesh, elem_type);
            else
                libmesh_error();
        }
        else if (mesh_type == "ringleb")
            RinglebMesh().init(nx_divs, ny_divs, mesh, elem_type);
        else
            libmesh_error();
    }
    else
    {
        mesh.set_mesh_dimension(dim);
        const std::string
        mesh_input_file = infile("mesh_input", std::string("")),
        mesh_type = infile("mesh_type", std::string(""));
        
        std::auto_ptr<MeshInput<UnstructuredMesh> > mesh_io;
        if (mesh_type == "gmsh")
            GmshIO(mesh).read(mesh_input_file);
        else if (mesh_type == "exodus")
            ExodusII_IO(mesh).read(mesh_input_file);
        else if (mesh_type == "exodus_parallel")
            ExodusII_IO(mesh).read_parallel(mesh_input_file);
        else if (mesh_type == "nemesis")
            Nemesis_IO(mesh).read(mesh_input_file);
        else
            libmesh_error();
        
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
    libMesh::EquationSystems equation_systems (mesh);
    equation_systems.parameters.set<GetPot*>("input_file") = &infile;
    
    // set data for flight condition
    FlightCondition flight_cond;
    for (unsigned int i=0; i<3; i++)
    {
        flight_cond.body_roll_axis(i)     = infile(    "body_roll_axis", 0., i);
        flight_cond.body_pitch_axis(i)    = infile(   "body_pitch_axis", 0., i);
        flight_cond.body_yaw_axis(i)      = infile(     "body_yaw_axis", 0., i);
        flight_cond.body_euler_angles(i)  = infile( "body_euler_angles", 0., i);
        flight_cond.body_angular_rates(i) = infile("body_angular_rates", 0., i);
    }
    flight_cond.ref_chord       = infile("ref_c",   1.);
    flight_cond.altitude        = infile( "alt",    0.);
    flight_cond.mach            = infile("mach",    .5);
    flight_cond.gas_property.cp = infile(  "cp", 1003.);
    flight_cond.gas_property.cv = infile(  "cv",  716.);
    flight_cond.gas_property.T  = infile("temp",  300.);
    flight_cond.gas_property.rho= infile( "rho",  1.05);
    
    flight_cond.init();
    
    
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    // Declare the system for fluid analysis
    FluidSystem & system =
    equation_systems.add_system<FluidSystem> ("FluidSystem");
    system.flight_condition = &flight_cond;
    
    system.attach_init_function (init_euler_variables);
    
    FluidPostProcessSystem& fluid_post =
    equation_systems.add_system<FluidPostProcessSystem> ("FluidPostProcessSystem");
    fluid_post.flight_condition = &flight_cond;
    
    ResidualBaseAdaptiveTimeSolver *timesolver = new ResidualBaseAdaptiveTimeSolver(system);
    libMesh::Euler2Solver *core_time_solver = new libMesh::Euler2Solver(system);
    
    timesolver->quiet              = infile("timesolver_solver_quiet", true);
    timesolver->growth_exponent    = infile("timesolver_growth_exponent", 1.2);
    timesolver->n_iters_per_update = infile("timesolver_update_n_iters", 10);
    timesolver->min_deltat         = infile("timesolver_min_deltat", 1.0e-3);
    timesolver->max_growth         = infile("timesolver_maxgrowth", 4.0);
    timesolver->min_growth         = infile("timesolver_mingrowth", 0.25);
    timesolver->max_deltat         = infile("timesolver_max_deltat", 5.0e2);

    core_time_solver->theta        = infile("timesolver_theta", 1.0);
    
    timesolver->core_time_solver = AutoPtr<UnsteadySolver>(core_time_solver);
    timesolver->diff_solver().reset(new NewtonSolver(system));
    system.time_solver = AutoPtr<UnsteadySolver>(timesolver);

    system.dc_recalculate_tolerance = infile("dc_recalculate_tolerance", 10.e-8);
    
    equation_systems.init ();
    
    AerodynamicQoI aero_qoi(infile);
    aero_qoi.flight_condition = &flight_cond;
    
    system.attach_qoi(&aero_qoi);
    
//    std::auto_ptr<SurfaceMotionFunction> surface_motion (new SurfaceMotionFunction);
//    //system.surface_motion = surface_motion.get();

//    // initialize the surface motion definition
//    std::auto_ptr<RigidSurfaceMotion> surface_motion(new RigidSurfaceMotion);
//    system.surface_motion = surface_motion.get();
//    
//    surface_motion->pitch_amplitude = 0.01745;
//    surface_motion->pitch_phase = 0.;
//    surface_motion->plunge_amplitude = 0.;
//    
//    surface_motion->pitch_axis(2) = 1.;
//    surface_motion->hinge_location(0) = 0.;
//    surface_motion->init(0., pi/2);

    
#else
    // Declare the system fluid system
    FrequencyDomainLinearizedFluidSystem & system =
    equation_systems.add_system<FrequencyDomainLinearizedFluidSystem>
    ("FrequencyDomainLinearizedFluidSystem");
    system.flight_condition = &flight_cond;

    FrequencyDomainFluidPostProcessSystem& fluid_post =
    equation_systems.add_system<FrequencyDomainFluidPostProcessSystem>
    ("DeltaFluidPostProcessSystem");
    fluid_post.flight_condition = &flight_cond;
    
    system.time_solver =
    AutoPtr<TimeSolver>(new SteadySolver(system));
    libmesh_assert_equal_to (n_timesteps, 1);
    
    equation_systems.read<libMesh::Real>("saved_solution.xdr", libMesh::DECODE,
                                (libMesh::EquationSystems::READ_HEADER |
                                 libMesh::EquationSystems::READ_DATA));
    
    // now initilaize the nonlinear solution
    system.localize_fluid_solution();
    
    // initialize the surface motion definition
    std::auto_ptr<MAST::RigidSurfaceMotion> surface_motion(new MAST::RigidSurfaceMotion);
    system.perturbed_surface_motion = surface_motion.get();
    
    surface_motion->pitch_amplitude = infile("pitch_ampl",0.);
    surface_motion->pitch_phase = infile("pitch_phase",0.);
    surface_motion->plunge_amplitude = infile("plunge_ampl",0.);

    for (unsigned int i=0; i<3; i++)
        surface_motion->pitch_axis(i) = infile("pitch_axis", 0., i);
    
    for (unsigned int i=0; i<3; i++)
        surface_motion->hinge_location(i) = infile("hinge", 0., i);
    
    for (unsigned int i=0; i<3; i++)
        surface_motion->plunge_vector(i) = infile("plunge_vec", 0., i);
    
    libMesh::Real frequency = infile("frequency",0.);
    surface_motion->init(frequency, 0.);
    equation_systems.parameters.set<bool>("if_reduced_freq") =
    infile("if_reduced_freq", false);
    
#endif
    
    
    system.print_residual_norms = infile("print_residual_norms", false);
    system.print_residuals = infile("print_residuals", false);
    system.print_jacobian_norms = infile("print_jacobian_norms", false);
    system.print_jacobians = infile("print_jacobians", false);
    //system.verify_analytic_jacobians = 1.0e-3;
    system.extra_quadrature_order = infile("extra_quadrature_order", 0);
    
    // Set the time stepping options
    system.deltat = deltat;
    
    // And the nonlinear solver options
    NewtonSolver &solver = dynamic_cast<NewtonSolver&>
    (*(system.time_solver->diff_solver().get()));
    solver.quiet = infile("solver_quiet", true);
    solver.verbose = !solver.quiet;
    solver.brent_line_search = false;
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
    unsigned int t_step=0, amr_steps = max_adaptivesteps, a_step = 0;
    libMesh::Real sol_norm = 1.0e10;
    if (!if_use_amr) amr_steps = 0;
    
    while (continue_iterations)
    {
        // A pretty update message
        std::cout << "\n\nSolving time step " << t_step << ", time = "
        << system.time << std::endl;
        
        
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
        // before the first iteration the dc vector needs to be localized
        system.evaluate_recalculate_dc_flag();
#endif

        // Do one last solve before the time step increment
        system.solve();
        //system.assemble_qoi(); // calculate the quantities of interest
        //system.postprocess(); // set the qois to the post-process variables
        
        // Advance to the next timestep in a transient problem
        //system.time_solver->advance_timestep();
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
        sol_norm = timesolver->_x_dot_norm_old;
#endif

        
        // Adaptively solve the timestep
        if (if_use_amr && (amr_steps > 0) && (sol_norm < amr_threshold))
        {
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
            std::vector<libMesh::Real> weights(dim+2,0.0);  // all set to 1.0
            weights[0] = 1.0;
            // Keep the same default norm type.
            std::vector<FEMNormType>
            norms(1, error_estimator->error_norm.type(0));
            error_estimator->error_norm = SystemNorm(norms, weights);
            
            error_estimator->estimate_error(system, error);
            
            // Print out status at each adaptive step.
            libMesh::Real global_error = error.l2_norm();
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
            
            MeshBase::const_element_iterator       el     = mesh.local_elements_begin();
            const MeshBase::const_element_iterator end_el = mesh.local_elements_end();
            for ( ; el != end_el; el++)
                (*el)->print_info();

            
            // reduce the time step size by a factor
            system.deltat *= amr_time_shrink_factor;
            
            std::cout << "Refined mesh to "
            << mesh.n_active_elem()
            << " active elements and "
            << equation_systems.n_active_dofs()
            << " active dofs." << std::endl;
            
            // decrement the amr counter
            amr_steps--;
            a_step++;

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
            // tell the solver to recalculte the dc coeffs post refinement
            system.if_use_stored_dc_coeff = false;
#endif
        }

        // check for termination criteria
        t_step++;

        if (((t_step >= n_timesteps) || (sol_norm < terminate_tolerance)) &&
            (amr_steps == 0))
        {
            libMesh::out << "\n === Terminating pseudo-time iterations ===" << std::endl;
            continue_iterations = false;
        }
        
        // Write out this timestep if we're requested to
        if ((t_step%write_interval == 0) || !continue_iterations)
        {
            fluid_post.postprocess();
            system.assemble_qoi(); // calculate the quantities of interest
            system.postprocess(); // set the qois to the post-process variables

            std::ostringstream file_name, b_file_name;
            
            // We write the file in the ExodusII format.
            file_name << "out_"
            << std::setw(3)
            << std::setfill('0')
            << std::right
            << t_step
            << ".exo";
            
            ExodusII_IO(mesh).write_equation_systems(file_name.str(),
                                                     equation_systems);
            
            
            // output of data along a line
//            std::ostringstream out_name;
//            out_name << "out_"
//            << std::setw(3)
//            << std::setfill('0')
//            << std::right
//            << t_step
//            << ".txt";
//            std::ofstream file;
//            file.open(out_name.str().c_str(), std::ofstream::out);
//            
//            std::vector<unsigned int> v(1); v[0] = system.variable_number("rho");
//            MeshFunction m(equation_systems, *system.solution, system.get_dof_map(), v);
//            m.init();
//            unsigned int ndivs=10000; libMesh::Real dx=5.0/(ndivs*1.);
//            libMesh::Point p; libMesh::DenseVector<libMesh::Real> vals; vals.resize(4);
//            p(0) = 0.; p(1) = .026;
//            while (p(0) < 5.) {
//                m(p, 0., vals);
//                file << p(0) << "  " << vals(0) << std::endl;
//                p(0) += dx;
//            }
//            file.close();

        }
        
    }
    
#ifndef LIBMESH_USE_COMPLEX_NUMBERS
    MeshSerializer serializer(mesh);
    XdrIO xdr(mesh, true);
    xdr.set_write_parallel(false);
    xdr.write("saved_mesh.xdr");
    equation_systems.write("saved_solution.xdr", libMesh::ENCODE,
                           (libMesh::EquationSystems::WRITE_SERIAL_FILES |
                            libMesh::EquationSystems::WRITE_DATA));
#endif
    
    // All done.
    return 0;
}
