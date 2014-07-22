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

// C++ includes
#include <ctime>
#include <sstream>
#include <iomanip>

// MAST includes
#include "Mesh/gaussian_bump_mesh.h"
#include "Flight/flight_condition.h"
#include "FluidElems/fluid_system.h"
#include "FluidElems/linearized_fluid_system.h"
#include "Solvers/residual_based_adaptive_time_solver.h"
#include "FluidElems/aerodynamic_qoi.h"

// libmesh includes
#include "libmesh/getpot.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/dof_map.h"
#include "libmesh/const_function.h"
#include "libmesh/zero_function.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/mesh_function.h"
#include "libmesh/newton_solver.h"
#include "libmesh/euler2_solver.h"
#include "libmesh/exodusII_io.h"





void
gaussian_bump_analysis(libMesh::LibMeshInit& init,
                       GetPot& fluid_infile,
                       const unsigned int p_order,
                       const unsigned int n_panel_divs,
                       const unsigned int n_farfield_divs,
                       const std::string& nm,
                       unsigned int& n_dofs,
                       Real& entropy_error) {
    
    // initialize the fluid data structures
    
    // first the flight condition
    FlightCondition flight_cond;
    for (unsigned int i=0; i<3; i++)
    {
        flight_cond.body_roll_axis(i)     = fluid_infile(    "body_roll_axis", 0., i);
        flight_cond.body_pitch_axis(i)    = fluid_infile(   "body_pitch_axis", 0., i);
        flight_cond.body_yaw_axis(i)      = fluid_infile(     "body_yaw_axis", 0., i);
        flight_cond.body_euler_angles(i)  = fluid_infile( "body_euler_angles", 0., i);
        flight_cond.body_angular_rates(i) = fluid_infile("body_angular_rates", 0., i);
    }
    flight_cond.ref_chord       = fluid_infile("ref_c",   1.);
    flight_cond.altitude        = fluid_infile( "alt",    0.);
    flight_cond.mach            = fluid_infile("mach",    .5);
    flight_cond.gas_property.cp = fluid_infile(  "cp", 1003.);
    flight_cond.gas_property.cv = fluid_infile(  "cv",  716.);
    flight_cond.gas_property.T  = fluid_infile("temp",  300.);
    flight_cond.gas_property.rho= fluid_infile( "rho",  1.05);
    
    flight_cond.init();
    
    // next the fluid mesh
    libMesh::SerialMesh fluid_mesh(init.comm());
    
    std::vector<Real> x_div_loc, x_relative_dx, y_div_loc, y_relative_dx;
    std::vector<unsigned int> x_divs, y_divs;
    std::auto_ptr<MeshInitializer::CoordinateDivisions>
    x_coord_divs (new MeshInitializer::CoordinateDivisions),
    y_coord_divs (new MeshInitializer::CoordinateDivisions);
    std::vector<MeshInitializer::CoordinateDivisions*> divs;
    
    const unsigned int
    fluid_dim     = fluid_infile("dimension",0),
    fluid_nx_divs = fluid_infile("nx_divs",0),
    fluid_ny_divs = fluid_infile("ny_divs",0),
    fluid_nz_divs = fluid_infile("nz_divs",0);
    divs.resize(fluid_dim);
    libMesh::ElemType fluid_elem_type =
    libMesh::Utility::string_to_enum<libMesh::ElemType>(fluid_infile("elem_type", "QUAD4"));
    
    x_div_loc.resize(fluid_nx_divs+1), x_relative_dx.resize(fluid_nx_divs+1),
    y_div_loc.resize(fluid_ny_divs+1), y_relative_dx.resize(fluid_ny_divs+1);
    x_divs.resize(fluid_nx_divs), y_divs.resize(fluid_ny_divs);
    std::auto_ptr<MeshInitializer::CoordinateDivisions>
    fluid_x_coord_divs (new MeshInitializer::CoordinateDivisions),
    fluid_y_coord_divs (new MeshInitializer::CoordinateDivisions);
    std::vector<MeshInitializer::CoordinateDivisions*> fluid_divs(fluid_dim);
    
    // now read in the values: x-coord
    for (unsigned int i_div=0; i_div<fluid_nx_divs+1; i_div++) {
        x_div_loc[i_div]     = fluid_infile("x_div_loc", 0., i_div);
        x_relative_dx[i_div] = fluid_infile( "x_rel_dx", 0., i_div);
    }
    x_divs[0] = n_farfield_divs;
    x_divs[1] = n_panel_divs;
    x_divs[2] = n_farfield_divs;
    
    divs[0] = x_coord_divs.get();
    x_coord_divs->init(fluid_nx_divs, x_div_loc, x_relative_dx, x_divs);
    
    // initialize the y-axis element divisions
    for (unsigned int i_div=0; i_div<fluid_ny_divs+1; i_div++) {
        y_div_loc[i_div]     = fluid_infile("y_div_loc", 0., i_div);
        y_relative_dx[i_div] = fluid_infile( "y_rel_dx", 0., i_div);
    }
    y_divs[0] = n_farfield_divs;
    divs[1] = y_coord_divs.get();
    y_coord_divs->init(fluid_ny_divs, y_div_loc, y_relative_dx, y_divs);
    
    
    const Real t_by_c =  fluid_infile("t_by_c", 0.0);

    // create the mesh
    GaussianBumpMesh2D().init(t_by_c,
                              divs,
                              fluid_mesh,
                              fluid_elem_type);
    
    // Print information about the mesh to the screen.
    fluid_mesh.print_info();
    
    // Create an equation systems object.
    libMesh::EquationSystems fluid_eq_systems(fluid_mesh);
    fluid_eq_systems.parameters.set<GetPot*>("input_file") = &fluid_infile;
    fluid_eq_systems.parameters.set<unsigned int>("p_order") = p_order;
    
    AerodynamicQoI aero_qoi(fluid_infile);
    aero_qoi.flight_condition = &flight_cond;
    
    
    // Declare the system
    FluidSystem& fluid_system =
    fluid_eq_systems.add_system<FluidSystem>("FluidSystem");
    
    FluidPostProcessSystem& fluid_post =
    fluid_eq_systems.add_system<FluidPostProcessSystem> ("FluidPostProcessSystem");
    fluid_post.flight_condition = &flight_cond;

    
    fluid_system.attach_init_function (init_euler_variables);

    fluid_system.attach_qoi(&aero_qoi);

    fluid_system.flight_condition = &flight_cond;
    
    fluid_system.attach_init_function(init_euler_variables);

    ResidualBaseAdaptiveTimeSolver* timesolver =
    new ResidualBaseAdaptiveTimeSolver(fluid_system);
    libMesh::Euler2Solver *core_time_solver = new libMesh::Euler2Solver(fluid_system);
    
    timesolver->quiet              = fluid_infile("timesolver_solver_quiet", true);
    timesolver->growth_exponent    = fluid_infile("timesolver_growth_exponent", 1.2);
    timesolver->n_iters_per_update = fluid_infile("timesolver_update_n_iters", 10);
    timesolver->min_deltat         = fluid_infile("timesolver_min_deltat", 1.0e-3);
    timesolver->max_growth         = fluid_infile("timesolver_maxgrowth", 4.0);
    timesolver->min_growth         = fluid_infile("timesolver_mingrowth", 0.25);
    timesolver->max_deltat         = fluid_infile("timesolver_max_deltat", 5.0e2);
    
    core_time_solver->theta        = fluid_infile("timesolver_theta", 1.0);
    
    timesolver->core_time_solver = libMesh::AutoPtr<libMesh::UnsteadySolver>(core_time_solver);
    timesolver->diff_solver().reset(new libMesh::NewtonSolver(fluid_system));
    fluid_system.time_solver = libMesh::AutoPtr<libMesh::UnsteadySolver>(timesolver);
    
    // now initilaize the nonlinear solution
    fluid_eq_systems.parameters.set<bool>("if_reduced_freq") =
    fluid_infile("if_reduced_freq", false);
    
    fluid_eq_systems.init();
    
    
    libMesh::NewtonSolver &solver = dynamic_cast<libMesh::NewtonSolver&>
    (*timesolver->diff_solver());
    solver.quiet = fluid_infile("solver_quiet", true);
    solver.verbose = !solver.quiet;
    solver.brent_line_search = false;
    solver.max_nonlinear_iterations =
    fluid_infile("max_nonlinear_iterations", 15);
    solver.relative_step_tolerance =
    fluid_infile("relative_step_tolerance", 1.e-3);
    solver.relative_residual_tolerance =
    fluid_infile("relative_residual_tolerance", 0.0);
    solver.absolute_residual_tolerance =
    fluid_infile("absolute_residual_tolerance", 0.0);
    solver.continue_after_backtrack_failure =
    fluid_infile("continue_after_backtrack_failure", false);
    solver.continue_after_max_iterations =
    fluid_infile("continue_after_max_iterations", false);
    solver.require_residual_reduction =
    fluid_infile("require_residual_reduction", true);
    
    // And the linear solver options
    solver.max_linear_iterations =
    fluid_infile("max_linear_iterations", 50000);
    solver.initial_linear_tolerance =
    fluid_infile("initial_linear_tolerance", 1.e-3);
    
    
    // the frequency domain
    //fluid_system_freq.localize_fluid_solution();
    
    // Print information about the system to the screen.
    fluid_eq_systems.print_info();

    
    // solve the system of equations till convergence
    bool continue_iterations = true;
    unsigned int
    t_step = 0,
    n_timesteps                    = fluid_infile("n_timesteps", 1);
    Real sol_norm = 0.,
    terminate_tolerance = fluid_infile("pseudo_time_terminate_tolerance", 1.0e-5);
    
    while (continue_iterations) {
        
        fluid_system.evaluate_recalculate_dc_flag();
        
        fluid_system.solve();
        
        fluid_system.assemble_qoi(); // calculate the quantities of interest
        fluid_system.postprocess(); // set the qois to the post-process variables
        fluid_post.postprocess();
        sol_norm = timesolver->_x_dot_norm_old;
        
        if ((t_step >= n_timesteps) ||
            (sol_norm < terminate_tolerance)) {
            
            libMesh::out << "\n === Terminating pseudo-time iterations ===" << std::endl;
            continue_iterations = false;
        }
        libMesh::ExodusII_IO(fluid_mesh).write_equation_systems("gaussian_bump.exo",
                                                                fluid_eq_systems);

    }
    
    // write the final solution to an output file
    libMesh::ExodusII_IO(fluid_mesh).write_equation_systems("gaussian_bump.exo",
                                                            fluid_eq_systems);
    
    // set the return values
    n_dofs        = fluid_system.n_dofs();
    entropy_error = sqrt(fluid_system.qoi[3]/fluid_system.qoi[2]);
}

int
main(int argc, char* const argv[]) {
    
    libMesh::LibMeshInit init(argc, argv);
    
    GetPot fluid_infile("system_input.in");
    unsigned int
    n_panel_divs = 4,
    n_farfield_divs = 4,
    n_increments = 5,
    n_dofs;
    Real entropy_error=0.;
    
    std::ofstream output;
    output.open("convergence_output.txt", std::ofstream::out);
    
    
    for (unsigned int p_order=1; p_order<5; p_order++) {

        // print the header
        output
        << std::setw(5) << "O"
        << std::setw(5) << "Msh"
        << std::setw(15) << "NDofs"
        << std::setw(35) << "entropy error"
        << std::setw(35) << "time" << std::endl;

        for (unsigned int i=0; i<n_increments; i++) {
            std::cout
            << "**************************************************************************" << std::endl
            << "             Analysis for p = " << p_order << "  mesh = " << i << std::endl
            << "**************************************************************************" << std::endl;
            
            std::clock_t start;
            Real duration;
            
            start = std::clock();
            
            std::ostringstream oss;
            oss << "p_" << p_order << "_m_" << i << "_";
            gaussian_bump_analysis(init,
                                   fluid_infile,
                                   p_order,
                                   n_panel_divs*(i+1),//pow(2,i),
                                   n_farfield_divs*(i+1),//pow(2,i),
                                   oss.str(),
                                   n_dofs,
                                   entropy_error);
            
            duration = ( std::clock() - start ) / (Real) CLOCKS_PER_SEC;
            
            output
            << std::setw(5) << p_order
            << std::setw(5) << i
            << std::setw(15) << n_dofs
            << std::setw(35) << std::setprecision(15) << entropy_error
            << std::setw(35) << std::setprecision(15) << duration << std::endl;
            
            std::cout << std::endl << std::endl;;
        }
        
        output << std::endl << std::endl;
    }
    
    return 0;
}

