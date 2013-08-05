
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
void assemble_matrices(EquationSystems& es,
                       const std::string& system_name);

void get_dirichlet_dofs(EquationSystems& es,
                        const std::string& system_name,
                        std::set<unsigned int>& dirichlet_dof_ids);


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
    FrequencyDomainLinearizedFluidSystem & fluid_system =
    fluid_equation_systems.add_system<FrequencyDomainLinearizedFluidSystem>
    ("FrequencyDomainLinearizedFluidSystem");
    fluid_system.flight_condition = &flight_cond;
    fluid_system.time_solver =
    AutoPtr<TimeSolver>(new SteadySolver(fluid_system));
    // read the fluid system from the saved file
    fluid_equation_systems.read<Real>("saved_solution.xdr",
                                      libMeshEnums::DECODE);
    // now initilaize the nonlinear solution
    fluid_system.localize_fluid_solution();
    fluid_system.extra_quadrature_order = 2;
    
    // print the information
    fluid_mesh.print_info();
    fluid_equation_systems.print_info();
    

    // *****************************************
    // now initialize the structural system
    // *****************************************
    SerialMesh structural_mesh(init.comm());
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

    // create the solvers
    UGFlutterSolver flutter_solver;
    
    // attach the fluid and structural systems to the models
    CFDAerodynamicModel aero_model(fluid_system);
    CoupledFluidStructureSystem
    coupled_system(aero_model, structural_model);
    
    ComplexMatrixX a;
    coupled_system.get_aero_operator_matrix(0.1, a);
}

#endif // LIBMESH_USE_COMPLEX_NUMBERS
