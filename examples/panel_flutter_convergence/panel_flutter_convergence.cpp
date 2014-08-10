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


// MAST includes
#include "Optimization/optimization_interface.h"
#include "StructuralElems/structural_system_assembly.h"
#include "PropertyCards/material_property_card_base.h"
#include "PropertyCards/solid_1d_section_element_property_card.h"
#include "PropertyCards/solid_2d_section_element_property_card.h"
#include "PropertyCards/multilayer_1d_section_element_property_card.h"
#include "PropertyCards/multilayer_2d_section_element_property_card.h"
#include "BoundaryConditions/temperature.h"
#include "BoundaryConditions/displacement_boundary_condition.h"
#include "Mesh/mesh_initializer.h"
#include "PropertyCards/isotropic_material_property_card.h"
#include "FluidElems/frequency_domain_linearized_fluid_system.h"
#include "Flight/flight_condition.h"
#include "FluidElems/fluid_system.h"
#include "FluidElems/frequency_domain_linearized_fluid_system.h"
#include "Aeroelasticity/pk_flutter_solver.h"
#include "Aeroelasticity/ug_flutter_solver.h"
#include "Aeroelasticity/noniterative_ug_flutter_solver.h"
#include "Aeroelasticity/coupled_fluid_structure_system.h"
#include "Mesh/panel_mesh.h"


// libmesh includes
#include "libmesh/getpot.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/eigen_solver.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/const_function.h"
#include "libmesh/zero_function.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/sensitivity_data.h"
#include "libmesh/condensed_eigen_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/mesh_function.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/euler2_solver.h"






void
panel_flutter_analysis(libMesh::LibMeshInit& init,
                       GetPot& str_infile,
                       GetPot& fluid_infile,
                       const unsigned int p_order,
                       const unsigned int n_panel_divs,
                       const unsigned int n_farfield_divs,
                       const std::string& nm,
                       unsigned int& n_dofs,
                       Real& flutter_V,
                       Real& flutter_g,
                       Real& flutter_omega ) {
    
    libMesh::SerialMesh beam_mesh(init.comm());
    
    const unsigned int
    dim     = 1,
    nx_divs = 1;
    libMesh::ElemType elem_type =
    libMesh::Utility::string_to_enum<libMesh::ElemType>("EDGE3");
    
    std::vector<Real> x_div_loc(nx_divs+1), x_relative_dx(nx_divs+1),
    y_div_loc, y_relative_dx;
    std::vector<unsigned int> x_divs(nx_divs), y_divs;
    std::auto_ptr<MeshInitializer::CoordinateDivisions>
    x_coord_divs (new MeshInitializer::CoordinateDivisions),
    y_coord_divs (new MeshInitializer::CoordinateDivisions);
    std::vector<MeshInitializer::CoordinateDivisions*> divs(dim);
    
    // now read in the values: x-coord
    x_div_loc[0] = 2.;  // panel LE
    x_div_loc[1] = 4.;  // panel TE
    
    x_relative_dx[0] = 1.;
    x_relative_dx[1] = 1.;
    
    x_divs[0] = n_panel_divs;
    
    divs[0] = x_coord_divs.get();
    x_coord_divs->init(nx_divs, x_div_loc, x_relative_dx, x_divs);
    
    // design data
    const unsigned int n_eig = 4;
    
    // now initialize the mesh
    MeshInitializer().init(divs, beam_mesh, elem_type);
    
    // Print information about the mesh to the screen.
    beam_mesh.print_info();
    
    // Create an equation systems object.
    libMesh::EquationSystems str_eq_systems(beam_mesh);
    str_eq_systems.parameters.set<GetPot*>("input_file") = &str_infile;
    
    // Declare the system
    libMesh::CondensedEigenSystem& eigen_system =
    str_eq_systems.add_system<libMesh::CondensedEigenSystem> ("EigenStructuralSystem");
    
    std::string fe_family = str_infile("fe_family", std::string("LAGRANGE"));
    libMesh::FEFamily fefamily = libMesh::Utility::string_to_enum<libMesh::FEFamily>(fe_family);
    libMesh::Order order = static_cast<libMesh::Order>(p_order);
    
    std::map<std::string, unsigned int> var_id;
    var_id["ux"] = eigen_system.add_variable ( "ux", order, fefamily);
    var_id["uy"] = eigen_system.add_variable ( "uy", order, fefamily);
    var_id["uz"] = eigen_system.add_variable ( "uz", order, fefamily);
    var_id["tx"] = eigen_system.add_variable ( "tx", order, fefamily);
    var_id["ty"] = eigen_system.add_variable ( "ty", order, fefamily);
    var_id["tz"] = eigen_system.add_variable ( "tz", order, fefamily);
    
    MAST::StructuralSystemAssembly eigen_structural_assembly(eigen_system,
                                                             MAST::MODAL,
                                                             str_infile);
    
    eigen_system.attach_assemble_object(eigen_structural_assembly);
    eigen_system.attach_eigenproblem_sensitivity_assemble_object(eigen_structural_assembly);
    
    // apply the boundary conditions
    libMesh::ZeroFunction<Real> zero_function;
    // Pass the Dirichlet dof IDs to the libMesh::CondensedEigenSystem
    std::set<libMesh::boundary_id_type> dirichlet_boundary;
    // read and initialize the boundary conditions
    std::map<libMesh::boundary_id_type, std::vector<unsigned int> > boundary_constraint_map;
    unsigned int n_bc, b_id;
    // first read the boundaries for ux constraint
    
    for (std::map<std::string, unsigned int>::iterator it = var_id.begin();
         it != var_id.end(); it++) {
        
        std::string nm = "n_" + it->first + "_bc"; // name for # bcs
        n_bc = str_infile(nm, 0);
        
        nm = it->first + "_bc";  // name for bc id vector
        
        for (unsigned int i=0; i<n_bc; i++) {
            b_id = str_infile(nm, 0, i); // ith component of the bc id vector
            
            if (!boundary_constraint_map.count(b_id)) // add vector if it does not exist
                boundary_constraint_map[b_id] = std::vector<unsigned int>(0);
            
            boundary_constraint_map[b_id].push_back(var_id[it->first]);
        }
    }
    
    // now iterate over each boundary and create the boudnary condition object
    unsigned int counter=0;
    std::vector<MAST::DisplacementDirichletBoundaryCondition*> bc(boundary_constraint_map.size());
    
    for (std::map<libMesh::boundary_id_type, std::vector<unsigned int> >::iterator
         it = boundary_constraint_map.begin();
         it != boundary_constraint_map.end(); it++) {
        bc[counter] = new MAST::DisplacementDirichletBoundaryCondition;
        bc[counter]->init(it->first, it->second);
        eigen_structural_assembly.add_side_load(it->first, *bc[counter]);
        counter++;
    }
    
    eigen_system.set_eigenproblem_type(libMesh::GHEP);
    eigen_system.eigen_solver->set_position_of_spectrum(libMesh::LARGEST_MAGNITUDE);
    
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
    
    const unsigned int fluid_dim     = fluid_infile("dimension",0),
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
    for (unsigned int i_div=0; i_div<fluid_nx_divs+1; i_div++)
    {
        x_div_loc[i_div]     = fluid_infile("x_div_loc", 0., i_div);
        x_relative_dx[i_div] = fluid_infile( "x_rel_dx", 0., i_div);
    }
    x_divs[0] = n_farfield_divs;
    x_divs[1] = n_panel_divs;
    x_divs[2] = n_farfield_divs;
    
    divs[0] = x_coord_divs.get();
    x_coord_divs->init(fluid_nx_divs, x_div_loc, x_relative_dx, x_divs);

    
    for (unsigned int i_div=0; i_div<fluid_ny_divs+1; i_div++)
    {
        y_div_loc[i_div]     = fluid_infile("y_div_loc", 0., i_div);
        y_relative_dx[i_div] = fluid_infile( "y_rel_dx", 0., i_div);
    }
    y_divs[0] = n_farfield_divs;
    divs[1] = y_coord_divs.get();
    y_coord_divs->init(fluid_ny_divs, y_div_loc, y_relative_dx, y_divs);
    
    
    const bool if_cos_bump = fluid_infile("if_cos_bump", false);
    const unsigned int n_max_bumps_x = fluid_infile("n_max_bumps_x", 1),
    n_max_bumps_y = fluid_infile("n_max_bumps_x", 1),
    panel_bc_id = fluid_infile("panel_bc_id", 10),
    symmetry_bc_id = fluid_infile("symmetry_bc_id", 11);
    const Real t_by_c =  fluid_infile("t_by_c", 0.0);
    
    PanelMesh2D().init(t_by_c, if_cos_bump, n_max_bumps_x,
                       panel_bc_id, symmetry_bc_id,
                       divs, fluid_mesh, fluid_elem_type);
    
    // Print information about the mesh to the screen.
    fluid_mesh.print_info();
    
    // Create an equation systems object.
    libMesh::EquationSystems fluid_eq_systems(fluid_mesh);
    fluid_eq_systems.parameters.set<GetPot*>("input_file") = &fluid_infile;
    fluid_eq_systems.parameters.set<unsigned int>("p_order") = p_order;
    
    // Declare the system
    FluidSystem& fluid_system_nonlin =
    fluid_eq_systems.add_system<FluidSystem>("FluidSystem");
    FrequencyDomainLinearizedFluidSystem& fluid_system_freq =
    fluid_eq_systems.add_system<FrequencyDomainLinearizedFluidSystem>
    ("FrequencyDomainLinearizedFluidSystem");
    
    fluid_system_nonlin.flight_condition = &flight_cond;
    fluid_system_freq.flight_condition = &flight_cond;
    
    fluid_system_nonlin.attach_init_function(init_euler_variables);
    
    fluid_system_nonlin.time_solver =
    libMesh::AutoPtr<libMesh::TimeSolver>(new libMesh::Euler2Solver(fluid_system_nonlin));
    fluid_system_freq.time_solver =
    libMesh::AutoPtr<libMesh::TimeSolver>(new libMesh::SteadySolver(fluid_system_freq));
    fluid_system_freq.time_solver->quiet = false;
    
    // now initilaize the nonlinear solution
    fluid_system_freq.extra_quadrature_order =
    fluid_infile("extra_quadrature_order", 0);
    fluid_eq_systems.parameters.set<bool>("if_reduced_freq") =
    fluid_infile("if_reduced_freq", false);
    
    str_eq_systems.init ();
    fluid_eq_systems.init();
    
    
    libMesh::NewtonSolver &solver = dynamic_cast<libMesh::NewtonSolver&>
    (*(fluid_system_freq.time_solver->diff_solver()));
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
    fluid_system_freq.localize_fluid_solution();
    
    
    // Print information about the system to the screen.
    str_eq_systems.print_info();
    fluid_eq_systems.print_info();
    
    
    // now initialize the flutter solver data structures
    FEMStructuralModel structural_model(eigen_structural_assembly);
    CFDAerodynamicModel aero_model(fluid_system_nonlin,
                                   fluid_system_freq);
    CoupledFluidStructureSystem coupled_system(aero_model,
                                               structural_model);
    coupled_system.nm = nm;
    
    // create the solvers
    MAST::PKFlutterSolver flutter_solver;
    std::string flutter_output_nm = nm + "_flutter_output.txt";
    if (!init.comm().rank())
        flutter_solver.set_output_file(flutter_output_nm);
    flutter_solver.aero_structural_model   = &coupled_system;
    flutter_solver.flight_condition        = &flight_cond;
    flutter_solver.k_red_range.first     = fluid_infile("ug_lower_k", 0.0);
    flutter_solver.k_red_range.second    = fluid_infile("ug_upper_k", 0.35);
    flutter_solver.n_k_red_divs          = fluid_infile("ug_k_divs", 10);
    flutter_solver.v_ref_range.first     = fluid_infile("ug_lower_V", 1.e2);
    flutter_solver.v_ref_range.second    = fluid_infile("ug_upper_V", 3.e2);
    flutter_solver.n_v_ref_divs          = fluid_infile("ug_V_divs", 1);
    
    // Pass the Dirichlet dof IDs to the libMesh::CondensedEigenSystem
    std::set<unsigned int> dirichlet_dof_ids;
    str_eq_systems.parameters.set<bool>("if_exchange_AB_matrices") = true;
    str_eq_systems.parameters.set<unsigned int>("eigenpairs")    = n_eig;
    str_eq_systems.parameters.set<unsigned int>("basis vectors") = n_eig*3;
    str_eq_systems.parameters.set<unsigned int>("nonlinear solver maximum iterations") =
    str_infile("max_nonlinear_iterations", 5);
    eigen_structural_assembly.get_dirichlet_dofs(dirichlet_dof_ids);
    eigen_system.initialize_condensed_dofs(dirichlet_dof_ids);
    
    
    
    // element and material properties
    MAST::ConstantFunction<Real>
    E("E", str_infile("youngs_modulus", 72.e9)),
    nu("nu", str_infile("poisson_ratio", 0.33)),
    rho("rho", str_infile("material_density", 2700.)),
    kappa("kappa", str_infile("shear_corr_factor", 5./6.)),
    alpha("alpha", str_infile("expansion_coefficient", 2.31e-5)),
    h_y("hy", str_infile("thickness", 0.002)),
    h_z("hz", str_infile("width", 0.002)),
    offset_h_y("hy_offset", 0.),
    offset_h_z("hz_offset", 0.);
    
    MAST::IsotropicMaterialPropertyCard mat(0);
    
    // add the properties to the cards
    mat.add(E);
    mat.add(nu);
    mat.add(alpha);
    mat.add(rho);
    mat.add(kappa);
    
    // create values for beam
    MAST::Solid1DSectionElementPropertyCard p(2);
    p.add(h_y); // thickness
    p.add(h_z); // width
    p.add(offset_h_y); // thickness offset
    p.add(offset_h_z); // width offset
    p.y_vector()(1) = 1.; // x-vector along x, y along y
    
    p.set_material(mat);

    eigen_structural_assembly.set_property_for_subdomain(0, p);
    
    // now add the vectors for structural modes, and init the fem structural model
    for (unsigned int i=0; i<n_eig; i++) {
        std::ostringstream vec_nm;
        vec_nm << "mode_" << i ;
        eigen_system.add_vector(vec_nm.str());
    }
    structural_model.eigen_vals.resize(n_eig);
    structural_model.init();
    
    
    libMesh::Point pt; // dummy point object
    
    // now solve the system
    eigen_system.solve();
    
    // the number of converged eigenpairs could be different from the number asked
    unsigned int n_required = std::min(n_eig, eigen_system.get_n_converged());
    
    std::pair<Real, Real> val;
    Complex eigval;
    
    structural_model.eigen_vals.resize(n_eig);
    
    if (str_eq_systems.parameters.get<bool>("if_exchange_AB_matrices"))
        // the total number of eigenvalues is n_eig, but only n_required are usable
        for (unsigned int i=0; i<n_required; i++) {
            val = eigen_system.get_eigenpair(i);
            eigval = Complex(val.first, val.second);
            eigval = 1./eigval;
            
            libMesh::out
            //<< std::setw(35) << std::fixed << std::setprecision(15)
            << eigval.real()
            << std::endl;
            
            // copy modal data into the structural model for flutter analysis
            structural_model.eigen_vals(i) = eigval.real();

            // get the mode
            std::ostringstream vec_nm;
            vec_nm << "mode_" << i;
            libMesh::NumericVector<Real>& vec = eigen_system.get_vector(vec_nm.str());
            eigen_system.solution->scale(sqrt(eigval.real()));
            vec = *(eigen_system.solution);
            vec.close();
            
            // output the mode to a file
            std::ostringstream oss;
            oss << nm <<  "_mode_" <<  i << ".exo";
            libMesh::ExodusII_IO(beam_mesh).write_equation_systems(oss.str(),
                                                                   str_eq_systems);
        }
    else
        libmesh_error(); // should not get here.
    
    // flutter solution
    flutter_solver.scan_for_roots();
    
    if (!init.comm().rank()) {
        flutter_solver.print_sorted_roots();
        flutter_solver.print_crossover_points();
    }
    
    std::pair<bool, const MAST::FlutterRootBase*> root =
    flutter_solver.find_critical_root(1.e-8,10);
    
    if (!init.comm().rank())
        flutter_solver.print_sorted_roots();
    
    n_dofs = fluid_system_freq.n_dofs();
    // get the values of n_dofs and flutter solution for return
    if (root.first) {
        flutter_V = root.second->V;
        flutter_g = root.second->g;
        flutter_omega = root.second->omega;
    }
    else {
        flutter_V     = -1.;
        flutter_g     = -1.;
        flutter_omega = -1.;
    }
    
    // delete the bc
    for (unsigned int i=0; i<bc.size(); i++)
        delete bc[i];
}

int
main(int argc, char* const argv[]) {
    
    libMesh::LibMeshInit init(argc, argv);
    
    GetPot str_infile("structural.in"),
    fluid_infile("system_input.in");
    
    unsigned int
    n_panel_divs = 3,
    n_farfield_divs = 4,
    n_increments = 5,
    n_dofs;
    Real flutter_V, flutter_g, flutter_omega;
    
    std::ofstream output;
    output.open("convergence_output.txt", std::ofstream::out);

    output
    << std::setw(5) << "O"
    << std::setw(5) << "Msh"
    << std::setw(15) << "NDofs"
    << std::setw(35) << "g"
    << std::setw(35) << "omega"
    << std::setw(35) << "V"
    << std::setw(35) << "time" << std::endl;
    
    
    for (unsigned int p_order=1; p_order<3; p_order++) {
        for (unsigned int i=1; i<3; i++) {
            libMesh::out
            << "**************************************************************************" << std::endl
            << "             Analysis for p = " << p_order << "  mesh = " << i << std::endl
            << "**************************************************************************" << std::endl;
            
            std::clock_t start;
            Real duration;
            
            start = std::clock();
            
            std::ostringstream oss;
            oss << "p_" << p_order << "_m_" << i << "_";
            panel_flutter_analysis(init,
                                   str_infile,
                                   fluid_infile,
                                   p_order,
                                   n_panel_divs*pow(2,i),
                                   n_farfield_divs*pow(2,i),
                                   oss.str(),
                                   n_dofs,
                                   flutter_V,
                                   flutter_g,
                                   flutter_omega);

            duration = ( std::clock() - start ) / (Real) CLOCKS_PER_SEC;
            
            output
            << std::setw(5) << p_order
            << std::setw(5) << i
            << std::setw(15) << n_dofs
            << std::setw(35) << std::setprecision(15) << flutter_g
            << std::setw(35) << std::setprecision(15) << flutter_omega
            << std::setw(35) << std::setprecision(15) << flutter_V
            << std::setw(35) << std::setprecision(15) << duration << std::endl;

            libMesh::out << std::endl << std::endl;;
        }
        
        output << std::endl << std::endl;
    }
    
    return 0;
}

