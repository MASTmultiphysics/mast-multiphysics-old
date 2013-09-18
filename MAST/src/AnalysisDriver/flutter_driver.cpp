
// C++ includes
#include <string>


// MAST includes
#include "Aeroelasticity/ug_flutter_solver.h"
#include "Flight/flight_condition.h"
#include "Aeroelasticity/coupled_fluid_structure_system.h"
#include "FluidElems/frequency_domain_linearized_fluid_system.h"


// libMesh includes
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/xdr_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/euler_solver.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/condensed_eigen_system.h"


#ifdef LIBMESH_USE_COMPLEX_NUMBERS

int main (int argc, char* const argv[])
{
    // Initialize libMesh.
    LibMeshInit init (argc, argv);

    GetPot infile("system_input.in");

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

    
    // *****************************************
    // initialize the fluid system
    // *****************************************
    ParallelMesh fluid_mesh(init.comm());
    fluid_mesh.read("saved_mesh.xdr");
    EquationSystems fluid_equation_systems (fluid_mesh);
    FrequencyDomainLinearizedFluidSystem & linearized_fluid_system =
    fluid_equation_systems.add_system<FrequencyDomainLinearizedFluidSystem>
    ("FrequencyDomainLinearizedFluidSystem");
    linearized_fluid_system.flight_condition = &flight_cond;
    linearized_fluid_system.time_solver =
    AutoPtr<TimeSolver>(new SteadySolver(linearized_fluid_system));
    linearized_fluid_system.time_solver->quiet = false;
    // read the fluid system from the saved file
    fluid_equation_systems.read<Real>("saved_solution.xdr",
                                      libMeshEnums::DECODE);
    // now initilaize the nonlinear solution
    linearized_fluid_system.localize_fluid_solution();
    linearized_fluid_system.extra_quadrature_order =
    infile("extra_quadrature_order", 0);
    fluid_equation_systems.parameters.set<bool>("if_reduced_freq") =
    infile("if_reduced_freq", false);
    
    // print the information
    fluid_mesh.print_info();
    fluid_equation_systems.print_info();
    
    NewtonSolver &solver = dynamic_cast<NewtonSolver&>
    (*(linearized_fluid_system.time_solver->diff_solver()));
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

    // *****************************************
    // now initialize the structural system
    // *****************************************
    ParallelMesh structural_mesh(init.comm());
    structural_mesh.read("saved_structural_mesh.xdr");
    
    EquationSystems structural_equation_systems (structural_mesh);
    CondensedEigenSystem & structural_system =
    structural_equation_systems.add_system<CondensedEigenSystem> ("Eigensystem");
    structural_equation_systems.read<Real>("saved_structural_solution.xdr",
                                           libMeshEnums::DECODE,
                                           (EquationSystems::READ_HEADER |
                                            EquationSystems::READ_DATA |
                                            EquationSystems::READ_ADDITIONAL_DATA));
    
    // Prints information about the system to the screen.
    structural_mesh.print_info();
    structural_equation_systems.print_info();
    
    
    // read the eigenvalues
    FEMStructuralModel structural_model(structural_system);
    GetPot modal_data("modal_data.in");
    unsigned int neig = modal_data("n_eig", 0);
    structural_model.eigen_vals.resize(neig);
    for (unsigned int i=0; i<neig; i++)
    {
        std::ostringstream oss;
        oss << "eig_" << i;
        structural_model.eigen_vals(i) = modal_data(oss.str(), 0.);
    }
    structural_model.init();
    
    // attach the fluid and structural systems to the models
    System& nonlinear_fluid_system =
    fluid_equation_systems.get_system<System>("EulerSystem");
    CFDAerodynamicModel aero_model(nonlinear_fluid_system,
                                   linearized_fluid_system);
    CoupledFluidStructureSystem
    coupled_system(aero_model, structural_model);

    // create the solvers
    std::string nm  = infile("flutter_output", "flutter_output.txt");
    UGFlutterSolver flutter_solver;
    if (!init.comm().rank())
        flutter_solver.set_output_file(nm);
    flutter_solver.aero_structural_model = &coupled_system;
    flutter_solver.flight_condition      = &flight_cond;
    flutter_solver.k_ref_range.first     = infile("ug_lower_k", 0.0);
    flutter_solver.k_ref_range.second    = infile("ug_upper_k", 0.35);
    flutter_solver.n_k_divs              = infile("ug_k_divs", 10);
    flutter_solver.scan_for_roots();
    if (!init.comm().rank())
        flutter_solver.print_crossover_points();
    flutter_solver.find_critical_root();
    if (!init.comm().rank())
        flutter_solver.print_sorted_roots();
}

#endif // LIBMESH_USE_COMPLEX_NUMBERS

