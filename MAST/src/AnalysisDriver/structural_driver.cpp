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
#include "StructuralElems/structural_system_assembly.h"
#include "PropertyCards/element_property_card_base.h"
#include "PropertyCards/element_property_card_3D.h"
#include "PropertyCards/element_property_card_2D.h"
#include "PropertyCards/element_property_card_1D.h"
#include "PropertyCards/material_property_card_base.h"
#include "Optimization/gcmma_optimization_interface.h"
#include "BoundaryConditions/boundary_condition.h"
#include "BoundaryConditions/displacement_boundary_condition.h"
#include "Optimization/topology_optimization.h"
#include "Mesh/mesh_initializer.h"
#include "Mesh/panel_mesh.h"
#include "Mesh/stiffened_panel.h"
#include "Mesh/nastran_io.h"

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
#include "libmesh/parameter_vector.h"
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/zero_function.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/exodusII_io.h"

#ifndef LIBMESH_USE_COMPLEX_NUMBERS


// The main program.
int structural_driver (LibMeshInit& init, GetPot& infile,
                       int argc, char* const argv[])
{
//    MAST::TopologyOptimization topology(init, infile);
//    MAST::GCMMAOptimizationInterface gcmma;
//    gcmma.attach_function_evaluation_object(topology);
//    gcmma.optimize();
//    return 0;
    
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
        ElemType elem_type =
        Utility::string_to_enum<ElemType>(infile("elem_type", "QUAD4"));
        
        std::vector<Real> x_div_loc(nx_divs+1), x_relative_dx(nx_divs+1),
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
            mesh.set_mesh_dimension(dim);
            if (dim == 1)
                MeshInitializer().init(divs, mesh, elem_type);
            else if (dim == 2)
                PanelMesh2D().init(0., false, 0, 100, 101, // dummy values not needed in structures
                                   divs, mesh, elem_type);
            else if (dim == 3)
                PanelMesh3D().init(0., false, 0, 0, 100, 101, // dummy values not needed in structures
                                   divs, mesh, elem_type);
            else
                libmesh_error();
        }
        else if (mesh_type == "stiffened_panel") {
            MAST::StiffenedPanel().init(divs, mesh, elem_type);
        }
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
    EquationSystems equation_systems (mesh);
    equation_systems.parameters.set<GetPot*>("input_file") = &infile;


    // Declare the system
    //NonlinearImplicitSystem & system =
    //equation_systems.add_system<NonlinearImplicitSystem> ("StructuralSystem");
    CondensedEigenSystem* eigen_system = NULL;
    CondensedEigenSystem & system =
    equation_systems.add_system<CondensedEigenSystem> ("StructuralSystem");
    eigen_system = dynamic_cast<CondensedEigenSystem*>(&system);
    
    
    unsigned int o = infile("fe_order", 1);
    std::string fe_family = infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    std::map<std::string, unsigned int> var_id;
    var_id["ux"] = system.add_variable ( "ux", static_cast<Order>(o), fefamily);
    var_id["uy"] = system.add_variable ( "uy", static_cast<Order>(o), fefamily);
    var_id["uz"] = system.add_variable ( "uz", static_cast<Order>(o), fefamily);
    var_id["tx"] = system.add_variable ( "tx", static_cast<Order>(o), fefamily);
    var_id["ty"] = system.add_variable ( "ty", static_cast<Order>(o), fefamily);
    var_id["tz"] = system.add_variable ( "tz", static_cast<Order>(o), fefamily);
    
    MAST::StructuralSystemAssembly structural_assembly(system,
                                                       MAST::MODAL,
                                                       infile);
    
    // Set the type of the problem, here we deal with
    // a generalized Hermitian problem.
    system.extra_quadrature_order = infile("extra_quadrature_order", 0);
    
    
    ConstFunction<Real> press(1.e2);
    MAST::BoundaryCondition bc(MAST::SURFACE_PRESSURE);
    bc.set_function(press);
    std::set<subdomain_id_type> ids;
    mesh.subdomain_ids(ids);
    //structural_assembly.add_side_load(2, bc);
    
    
    system.attach_assemble_object(structural_assembly);
    system.attach_sensitivity_assemble_object(structural_assembly);
    
    // Set the number of requested eigenpairs \p n_evals and the number
    // of basis vectors used in the solution algorithm.
    const unsigned int n_eig_request = 10;
    std::vector<Real> sens;
    // Pass the Dirichlet dof IDs to the CondensedEigenSystem

    std::set<boundary_id_type> dirichlet_boundary;
    // read and initialize the boundary conditions
    std::map<boundary_id_type, std::vector<unsigned int> > boundary_constraint_map;
    unsigned int n_bc, b_id;
    // first read the boundaries for ux constraint

    for (std::map<std::string, unsigned int>::iterator it = var_id.begin();
         it != var_id.end(); it++) {
        
        std::string nm = "n_" + it->first + "_bc"; // name for # bcs
        n_bc = infile(nm, 0);
        
        nm = it->first + "_bc";  // name for bc id vector
        
        for (unsigned int i=0; i<n_bc; i++) {
            b_id = infile(nm, 0, i); // ith component of the bc id vector
            
            if (!boundary_constraint_map.count(b_id)) // add vector if it does not exist
                boundary_constraint_map[b_id] = std::vector<unsigned int>(0);
            
            boundary_constraint_map[b_id].push_back(var_id[it->first]);
        }
    }
    
    std::vector<MAST::BoundaryCondition*> dirichlet_boundary_conditions;
    
    // now iterate over each boundary and create the boudnary condition object
    for (std::map<boundary_id_type, std::vector<unsigned int> >::iterator
         it = boundary_constraint_map.begin();
         it != boundary_constraint_map.end(); it++) {
        MAST::DisplacementBoundaryCondition* bc = new MAST::DisplacementBoundaryCondition;
        bc->init(it->first, it->second);
        dirichlet_boundary_conditions.push_back(bc);
        structural_assembly.add_side_load(it->first, *bc);
    }
    
    if (eigen_system) {
        
        eigen_system->set_eigenproblem_type(GHEP);
        eigen_system->eigen_solver->set_position_of_spectrum(LARGEST_MAGNITUDE);
    }
    
    equation_systems.init ();
    equation_systems.print_info();
    
    // apply the boundary conditions to the eigenproblem if necessary
    if (eigen_system) {        
        std::set<unsigned int> dirichlet_dof_ids;
        equation_systems.parameters.set<bool>("if_exchange_AB_matrices") = true;
        equation_systems.parameters.set<unsigned int>("eigenpairs")    = n_eig_request;
        equation_systems.parameters.set<unsigned int>("basis vectors") = n_eig_request*3;
        
        eigen_system->attach_eigenproblem_sensitivity_assemble_object(structural_assembly);
        structural_assembly.get_dirichlet_dofs(dirichlet_dof_ids);
        eigen_system->initialize_condensed_dofs(dirichlet_dof_ids);
    }
    

    
//    system.time_solver->diff_solver()->quiet = false;
//    system.time_solver->diff_solver()->verbose = true;
//    system.time_solver->diff_solver()->relative_residual_tolerance =  1.0e-8;
//    system.time_solver->diff_solver()->absolute_residual_tolerance =  1.0e-8;
    
    // Print information about the system to the screen.
    
    MAST::IsotropicMaterialPropertyCard mat;
    MAST::ElementPropertyCard3D prop3d;
    MAST::Solid2DSectionElementPropertyCard prop2d;
    MAST::Solid1DSectionElementPropertyCard prop1d;
    
    MAST::FunctionValue<Real>& E = mat.add<Real>("E", MAST::CONSTANT_FUNCTION),
    &nu = mat.add<Real>("nu", MAST::CONSTANT_FUNCTION),
    &rho = mat.add<Real>("rho", MAST::CONSTANT_FUNCTION),
    &kappa = mat.add<Real>("kappa", MAST::CONSTANT_FUNCTION),
    &h =  prop2d.add<Real>("h", MAST::CONSTANT_FUNCTION),
    &th =  prop1d.add<Real>("h", MAST::CONSTANT_FUNCTION),
    &b =  prop1d.add<Real>("b", MAST::CONSTANT_FUNCTION);

    
    E  = infile("youngs_modulus", 72.e9);
    nu = infile("poisson_ratio", 0.33);
    rho =infile("material_density", 2700.);
    kappa = infile("shear_corr_factor", 5./6.);
    h  = infile("thickness", 0.002);
    th  = infile("thickness", 0.002);
    b  = infile("width", 0.002);
    
    DenseVector<Real> prestress; prestress.resize(6);
    prestress(0) = -100.;
    
    prop3d.set_material(mat);
    prop2d.set_material(mat);
    prop2d.set_diagonal_mass_matrix(false);
    prop1d.set_material(mat);
    prop1d.set_diagonal_mass_matrix(false);
    //prop2d.prestress(prestress);
    //prop1d.prestress(prestress);

    //prop2d.set_strain(MAST::VON_KARMAN_STRAIN);
    //prop1d.set_strain(MAST::VON_KARMAN_STRAIN);
    ParameterVector parameters; parameters.resize(1);
    parameters[0] = h.ptr(); // set thickness as a modifiable parameter
    
    
    structural_assembly.set_property_for_all_elems(prop2d);
    structural_assembly.add_parameter(h.ptr(), &h);
    
    MAST::NastranIO(structural_assembly).write("nast.txt");
    return 0;
    
    system.solve();
//    if (eigen_system)
//        eigen_system->sensitivity_solve(parameters, sens);
//    else
//        system.sensitivity_solve(parameters);
    

    if (!eigen_system) {
        
        // We write the file in the ExodusII format.
        Nemesis_IO(mesh).write_equation_systems("out.exo",
                                                equation_systems);
    }
    else {
        
        // Get the number of converged eigen pairs.
        unsigned int nconv = std::min(eigen_system->get_n_converged(),
                                      n_eig_request);
        
        std::ofstream file;
        file.open("modal_data.in", std::ofstream::out);
        
        file << "n_eig = " << nconv << std::endl;
        
        for (unsigned int i=0; i<nconv; i++)
        {
            std::ostringstream file_name;
            
            // We write the file in the ExodusII format.
            file_name << "out_"
            << std::setw(3)
            << std::setfill('0')
            << std::right
            << i
            << ".exo";
            
            // now write the eigenvlaues
            std::pair<Real, Real> val = eigen_system->get_eigenpair(i);

            // We write the file in the ExodusII format.
            Nemesis_IO(mesh).write_equation_systems(file_name.str(),
                                                    equation_systems);
            

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

                std::cout << std::setw(35) << std::fixed << std::setprecision(15) << eigval.real();

//                std::cout << std::setw(5) << i
//                << std::setw(10) << eigval.real()
//                << " + i  "
//                << std::setw(10) << eigval.imag();
                
//                // write the sensitivity
//                Complex eig_sens = -sens[i]/pow(Complex(val.first, val.second),2);
//                if (sens.size())
//                    std::cout << "  deig/dp = "
//                    << std::setw(10) << eig_sens.real()
//                    << " + i  "
//                    << std::setw(10) << eig_sens.imag();
//                
//                std::cout << std::endl;
                
            }
            else {
                file << "eig_"  << i << " = " << val.first << std::endl;
                
                std::cout << std::setw(5) << i
                << std::setw(10) << val.first
                << " + i  "
                << std::setw(10) << val.second;
                
//                // write the sensitivity
//                if (sens.size())
//                    std::cout << "  deig/dp = " << sens[i];

                std::cout<< std::endl;
            }
        }
        std::cout<< std::endl;
        file.close();
    }

    // be sure to delete the boundary condition objects
    for (unsigned int i=0; i<dirichlet_boundary_conditions.size(); i++)
        delete dirichlet_boundary_conditions[i];
    
    // now write the data to an output file
    XdrIO xdr(mesh, true);
    xdr.write("saved_structural_mesh.xdr");
    equation_systems.write("saved_structural_solution.xdr",
                           libMeshEnums::ENCODE,
                           (EquationSystems::WRITE_SERIAL_FILES |
                            EquationSystems::WRITE_DATA |
                            EquationSystems::WRITE_ADDITIONAL_DATA));
    
    // All done.
    return 0;
}

#endif // LIBMESH_USE_COMPLEX_NUMBERS
