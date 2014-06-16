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
#include "PropertyCards/isotropic_material_property_card.h"
#include "PropertyCards/element_property_card_3D.h"
#include "PropertyCards/multilayer_1d_section_element_property_card.h"
#include "PropertyCards/solid_1d_section_element_property_card.h"
#include "PropertyCards/solid_2d_section_element_property_card.h"
#include "BoundaryConditions/boundary_condition.h"
#include "BoundaryConditions/temperature.h"
#include "BoundaryConditions/displacement_boundary_condition.h"
#include "Mesh/mesh_initializer.h"
#include "Mesh/panel_mesh.h"
#include "Mesh/stiffened_panel.h"
#include "Mesh/nastran_io.h"
#include "Numerics/constant_function.h"



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




// The main program.
int structural_driver (libMesh::LibMeshInit& init, GetPot& infile,
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
    const std::string output_mesh_file   = infile("output_mesh_file", "saved_structural_mesh.xdr"),
    output_solution_file   = infile("output_solution_file", "saved_structural_solution.xdr");
    
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
    

    if (if_panel_mesh)
    {
        
        const unsigned int nx_divs = infile("nx_divs",0),
        ny_divs = infile("ny_divs",0),
        nz_divs = infile("nz_divs",0);
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
            bool beam_stiff = infile("beam_stiffeners", false);
            MAST::StiffenedPanel().init(divs, mesh, elem_type, beam_stiff);
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
        else if (mesh_type == "xdr")
            mesh.read(mesh_input_file);
        else
            libmesh_error();
        
        mesh.prepare_for_use();
    }
    
    

    
    // uniformly refine the mesh
    for (unsigned int i=0; i<n_uniform_refine; i++)
        mesh_refinement.uniformly_refine();
    
    // Print information about the mesh to the screen.
    mesh.print_info();
    
    // Create an equation systems object.
    libMesh::EquationSystems equation_systems (mesh);
    equation_systems.parameters.set<GetPot*>("input_file") = &infile;


    // Declare the system
    NonlinearImplicitSystem & static_system =
    equation_systems.add_system<NonlinearImplicitSystem> ("StaticStructuralSystem");
    CondensedEigenSystem & eigen_system =
    equation_systems.add_system<CondensedEigenSystem> ("StructuralSystem");
    
    
    unsigned int o = infile("fe_order", 1);
    std::string fe_family = infile("fe_family", std::string("LAGRANGE"));
    FEFamily fefamily = Utility::string_to_enum<FEFamily>(fe_family);
    
    std::map<std::string, unsigned int> var_id;
    var_id["ux"] = static_system.add_variable ( "ux", static_cast<Order>(o), fefamily);
    var_id["uy"] = static_system.add_variable ( "uy", static_cast<Order>(o), fefamily);
    var_id["uz"] = static_system.add_variable ( "uz", static_cast<Order>(o), fefamily);
    var_id["tx"] = static_system.add_variable ( "tx", static_cast<Order>(o), fefamily);
    var_id["ty"] = static_system.add_variable ( "ty", static_cast<Order>(o), fefamily);
    var_id["tz"] = static_system.add_variable ( "tz", static_cast<Order>(o), fefamily);

    eigen_system.add_variable ( "ux", static_cast<Order>(o), fefamily);
    eigen_system.add_variable ( "uy", static_cast<Order>(o), fefamily);
    eigen_system.add_variable ( "uz", static_cast<Order>(o), fefamily);
    eigen_system.add_variable ( "tx", static_cast<Order>(o), fefamily);
    eigen_system.add_variable ( "ty", static_cast<Order>(o), fefamily);
    eigen_system.add_variable ( "tz", static_cast<Order>(o), fefamily);

    MAST::StructuralSystemAssembly
    static_structural_assembly(static_system,
                               MAST::STATIC,
                               infile),
    eigen_structural_assembly(eigen_system,
                              MAST::MODAL,
                              infile);
    
    // Set the type of the problem, here we deal with
    // a generalized Hermitian problem.
    static_system.extra_quadrature_order = infile("extra_quadrature_order", 0);
    eigen_system.extra_quadrature_order = infile("extra_quadrature_order", 0);
    
    
    MAST::ConstantFunction<libMesh::Real> press("pressure", 1.e2);
    MAST::BoundaryCondition bc(MAST::SURFACE_PRESSURE);
    bc.set_function(press);
    //static_structural_assembly.add_volume_load(0, bc);
    //eigen_structural_assembly.add_volume_load(0, bc);
    MAST::ConstantFunction<libMesh::Real> temp("temp", 10.), ref_temp("ref_temp", 0.);
    MAST::Temperature temp_bc;
    temp_bc.set_function(temp);
    temp_bc.set_reference_temperature_function(ref_temp);
    static_structural_assembly.add_volume_load(0, temp_bc);
    eigen_structural_assembly.add_volume_load(0, temp_bc);
    
    static_system.attach_assemble_object(static_structural_assembly);
    eigen_system.attach_assemble_object(eigen_structural_assembly);
    
    // Pass the Dirichlet dof IDs to the CondensedEigenSystem
    std::set<libMesh::boundary_id_type> dirichlet_boundary;
    // read and initialize the boundary conditions
    std::map<libMesh::boundary_id_type, std::vector<unsigned int> > boundary_constraint_map;
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
    for (std::map<libMesh::boundary_id_type, std::vector<unsigned int> >::iterator
         it = boundary_constraint_map.begin();
         it != boundary_constraint_map.end(); it++) {
        MAST::DisplacementDirichletBoundaryCondition* bc =
        new MAST::DisplacementDirichletBoundaryCondition;
        bc->init(it->first, it->second);
        dirichlet_boundary_conditions.push_back(bc);
        static_structural_assembly.add_side_load(it->first, *bc);
        eigen_structural_assembly.add_side_load(it->first, *bc);
    }
    
    // this needs to be done before equation system init
    eigen_system.set_eigenproblem_type(GHEP);
    eigen_system.eigen_solver->set_position_of_spectrum(LARGEST_MAGNITUDE);
    
    equation_systems.init ();
    equation_systems.print_info();
    
    // apply the boundary conditions to the eigenproblem if necessary
    const unsigned int n_eig_request = 10;
    std::set<unsigned int> dirichlet_dof_ids;
    equation_systems.parameters.set<bool>("if_exchange_AB_matrices") = true;
    equation_systems.parameters.set<unsigned int>("eigenpairs")    = n_eig_request;
    equation_systems.parameters.set<unsigned int>("basis vectors") = n_eig_request*3;
    equation_systems.parameters.set<unsigned int>("nonlinear solver maximum iterations") =
    infile("max_nonlinear_iterations", 5);
    eigen_structural_assembly.get_dirichlet_dofs(dirichlet_dof_ids);
    eigen_system.initialize_condensed_dofs(dirichlet_dof_ids);

    MAST::IsotropicMaterialPropertyCard mat(0);
    MAST::ElementPropertyCard3D prop3d(0);
    MAST::Solid2DSectionElementPropertyCard prop2d(1), prop2d_stiff(2);
    MAST::Solid1DSectionElementPropertyCard prop1d(3);
    
    // create the scalar function values
    MAST::ConstantFunction<libMesh::Real> E("E", infile("youngs_modulus", 72.e9)),
    nu("nu", infile("poisson_ratio", 0.33)),
    alpha("alpha", infile("expansion_coefficient", 2.31e-5)),
    rho("rho", infile("material_density", 2700.)),
    kappa("kappa", infile("shear_corr_factor", 5./6.)),
    h("h", infile("thickness", 0.002)),
    h_stiff("h", infile("thickness", 0.002)),
    hy("hy", 1.00001*infile("thickness", 0.002)),
    hz("hz", infile("width", 0.002)),
    h_off("off", 0.), // plate offset
    off_hy("hy_offset", 0.5*infile("thickness", 0.002)),
    off_hz("hz_offset", 0.);
    
    // add the properties to the cards
    mat.add(E);
    mat.add(nu);
    mat.add(rho);
    mat.add(kappa);
    mat.add(alpha);
    
    prop2d.add(h);
    prop2d.add(h_off);
    prop2d_stiff.add(h_stiff);
    prop2d_stiff.add(h_off);
    prop1d.add(hy);
    prop1d.add(hz);
    prop1d.add(off_hy);
    prop1d.add(off_hz);
    
    
    DenseRealMatrix prestress; prestress.resize(3,3);
    prestress(0,0) = -1.31345e6/.0044*1.0e-3;
    MAST::ConstantFunction<DenseRealMatrix > prestress_func("prestress", prestress);
    
    prop3d.set_material(mat);
    prop2d.set_material(mat); prop2d_stiff.set_material(mat);
    prop2d.set_diagonal_mass_matrix(false); prop2d_stiff.set_diagonal_mass_matrix(false);
    prop1d.set_material(mat);
    prop1d.set_diagonal_mass_matrix(false);
    prop1d.y_vector()(1) = 1.;
    //prop2d.add(prestress_func); // no prestress for stiffener
    prop1d.add(prestress_func);

    
    // multilayer material
    MAST::Multilayer1DSectionElementPropertyCard multi_prop1d(4);
    MAST::Solid1DSectionElementPropertyCard bottom(5), middle(6), top(7);
    MAST::IsotropicMaterialPropertyCard void_material(8);
    MAST::ConstantFunction<libMesh::Real> zeroE("E", 0.),
    hz_bottom("hz", 0.01),
    hz_middle("hz", 0.02),
    hz_top("hz", 0.01);
    void_material.add(zeroE);
    void_material.add(nu);
    void_material.add(kappa);
    // set material
    bottom.set_material(mat);
    middle.set_material(void_material);
    top.set_material(mat);
    // set thickness and width
    // offset_hz is going to be set by the multilayer section card
    bottom.add(hy);
    bottom.add(hz_bottom);
    bottom.add(off_hy);
    middle.add(hy);
    middle.add(hz_middle);
    middle.add(off_hy);
    top.add(hy);
    top.add(hz_top);
    top.add(off_hy);
    // set von karman strain
    multi_prop1d.set_strain(MAST::VON_KARMAN_STRAIN);
    multi_prop1d.y_vector()(1) = 1.;
    std::vector<MAST::Solid1DSectionElementPropertyCard*> layers(3);
    layers[0] = &bottom;
    layers[1] = &middle;
    layers[2] = &top;
    multi_prop1d.set_layers(-1., layers); // offset wrt bottom layer
    
    //prop2d.set_strain(MAST::VON_KARMAN_STRAIN); prop2d_stiff.set_strain(MAST::VON_KARMAN_STRAIN);
    prop1d.set_strain(MAST::VON_KARMAN_STRAIN);
    
    if (dim == 1) {
        static_structural_assembly.set_property_for_subdomain(0, prop1d);
        eigen_structural_assembly.set_property_for_subdomain(0, prop1d);
    }
    else {
        static_structural_assembly.set_property_for_subdomain(0, prop2d);
        eigen_structural_assembly.set_property_for_subdomain(0, prop2d);
    }
    
    const std::string mesh_type = infile("mesh_type", std::string(""));
    if (mesh_type == "stiffened_panel") {
        bool beam_stiff = infile("beam_stiffeners", false);
        unsigned int n_stiff = (infile("nx_divs",0)-1) + (infile("ny_divs",0)-1);
        if (!beam_stiff) {
            // stiffeners using shell elements
            for (unsigned int i=1; i<n_stiff+1; i++) {
                static_structural_assembly.set_property_for_subdomain(i, prop2d_stiff);
                eigen_structural_assembly.set_property_for_subdomain(i, prop2d_stiff);
                static_structural_assembly.add_volume_load(i, temp_bc);
                eigen_structural_assembly.add_volume_load(i, temp_bc);
            }
        }
        else {
            // stiffeners using beam elements with offsets
            off_hz = 0.5*infile("width", 0.002);
            for (unsigned int i=1; i<n_stiff+1; i++) {
                static_structural_assembly.set_property_for_subdomain(i, prop1d);
                eigen_structural_assembly.set_property_for_subdomain(i, prop1d);
                static_structural_assembly.add_volume_load(i, temp_bc);
                eigen_structural_assembly.add_volume_load(i, temp_bc);
            }

//            // iterate over all elements in the region with void and set the
//            // domain id and material property
//            static_structural_assembly.set_property_for_subdomain(n_stiff+1, multi_prop1d);
//            eigen_structural_assembly.set_property_for_subdomain(n_stiff+1, multi_prop1d);
//            MeshBase::const_element_iterator
//            el_it = mesh.elements_begin(),
//            el_end = mesh.elements_end();
//            for ( ; el_it != el_end; el_it++ ) {
//                libMesh::Elem* elem = *el_it;
//                if (elem->dim() != 1)
//                    continue;
//                libMesh::Point p = elem->centroid();
//                if (p(0) >= .1 && p(0) <= .2)
//                    elem->subdomain_id() = n_stiff+1;
//            }
        }
    }

    ParameterVector params;
    params.resize(1); params[0] = hy.ptr();
    static_system.solve();
    static_system.attach_sensitivity_assemble_object(static_structural_assembly);
    static_structural_assembly.add_parameter(hy);
    static_system.sensitivity_solve(params);
    static_system.solution->print();
    static_system.get_sensitivity_solution().print();
    //return 0;
    
    eigen_structural_assembly.set_static_solution_system(&static_system);
    //prop2d.set_strain(MAST::VON_KARMAN_STRAIN); prop2d_stiff.set_strain(MAST::VON_KARMAN_STRAIN);
    //prop1d.set_strain(MAST::VON_KARMAN_STRAIN);
    //eigen_system.solve();
    std::vector<libMesh::Real> sens;
    eigen_structural_assembly.add_parameter(hy);
    eigen_system.attach_eigenproblem_sensitivity_assemble_object(eigen_structural_assembly);
    //eigen_system.sensitivity_solve(params, sens);
    
    std::vector<libMesh::Real> stress;
    //static_structural_assembly.calculate_max_elem_stress(*static_system.solution,
    //                                                     stress, NULL);
    
    // We write the file in the ExodusII format.
    Nemesis_IO(mesh).write_equation_systems("out.exo",
                                            equation_systems);
    // Get the number of converged eigen pairs.
    unsigned int nconv = std::min(eigen_system.get_n_converged(),
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
        std::pair<Real, Real> val = eigen_system.get_eigenpair(i);
        
        // We write the file in the ExodusII format.
        Nemesis_IO(mesh).write_equation_systems(file_name.str(),
                                                equation_systems);
        
        
        // also add the solution as an independent vector, which will have to
        // be read in
        std::ostringstream vec_name;
        vec_name << "mode_" << i;
        libMesh::NumericVector<libMesh::Real>& vec = eigen_system.add_vector(vec_name.str());
        vec = *eigen_system.solution;
        std::complex<libMesh::Real> eigval;
        std::streamsize prec = std::cout.precision();
        if (equation_systems.parameters.get<bool>("if_exchange_AB_matrices"))
        {
            file << "eig_"  << i << " = "
            << std::setw(35) << std::setprecision(15) << 1./val.first
            << std::setprecision(prec) << std::endl;
            vec.scale(1./sqrt(val.first));
            vec.close();
            
            // now write the eigenvalues
            eigval = std::complex<libMesh::Real>(val.first, val.second);
            eigval = 1./eigval;
            std::cout << std::setw(35) << std::fixed << std::setprecision(15) << eigval.real();
            
            if (sens.size() > 0) {
                sens[i] *= -1./pow(val.first,2);
                std::cout << std::setw(35) << std::fixed << std::setprecision(15) << sens[i];
            }
            
            std::cout << std::endl;

        }
        else {
            file << "eig_"  << i << " = "
            << std::setw(35) << std::setprecision(15) << val.first
            << std::setprecision(prec) << std::endl;
            
            std::cout << std::setw(5) << i
            << std::setw(10) << val.first
            << " + i  "
            << std::setw(10) << val.second;
            
            std::cout<< std::endl;
        }
    }
    std::cout<< std::endl;
    file.close();
    
    // now write the data to an output file
    MAST::NastranIO(static_structural_assembly).write("nast.txt");
    XdrIO xdr(mesh, true);
    xdr.write(output_mesh_file);
    equation_systems.write(output_solution_file,
                           libMesh::ENCODE,
                           (libMesh::EquationSystems::WRITE_SERIAL_FILES |
                            libMesh::EquationSystems::WRITE_DATA |
                            libMesh::EquationSystems::WRITE_ADDITIONAL_DATA));

    // be sure to delete the boundary condition objects
    for (unsigned int i=0; i<dirichlet_boundary_conditions.size(); i++)
        delete dirichlet_boundary_conditions[i];
    
    // All done.
    return 0;
}


