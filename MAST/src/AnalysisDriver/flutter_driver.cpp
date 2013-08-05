
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

    // *****************************************
    // initialize the fluid system
    // *****************************************
    ParallelMesh fluid_mesh(init.comm());
    fluid_mesh.read("saved_mesh.xdr");
    EquationSystems fluid_equation_systems (fluid_mesh);
    FrequencyDomainLinearizedFluidSystem & fluid_system =
    fluid_equation_systems.add_system<FrequencyDomainLinearizedFluidSystem>
    ("FrequencyDomainLinearizedFluidSystem");
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
    structural_mesh.set_mesh_dimension(1);
    MeshTools::Generation::build_square (structural_mesh,
                                         10,
                                         0., 1.,
                                         EDGE2);
    
    structural_mesh.prepare_for_use();
    
    EquationSystems structural_equation_systems (structural_mesh);
    CondensedEigenSystem & structural_system =
    structural_equation_systems.add_system<CondensedEigenSystem> ("Eigensystem");
    
    structural_system.add_variable ( "ux", FIRST, LAGRANGE);
    structural_system.add_variable ( "uy", FIRST, LAGRANGE);
    structural_system.add_variable ( "uz", FIRST, LAGRANGE);
    structural_system.add_variable ( "tx", FIRST, LAGRANGE);
    structural_system.add_variable ( "ty", FIRST, LAGRANGE);
    structural_system.add_variable ( "tz", FIRST, LAGRANGE);

    // Give the system a pointer to the matrix assembly
    // function defined below.
    structural_system.attach_assemble_function (assemble_matrices);
    
    // Set the type of the problem, here we deal with
    // a generalized Hermitian problem.
    structural_system.set_eigenproblem_type(GHEP);
    
    // Order the eigenvalues "smallest first"
    structural_system.eigen_solver->set_position_of_spectrum(LARGEST_MAGNITUDE);
    
    // Set the number of requested eigenpairs \p n_evals and the number
    // of basis vectors used in the solution algorithm.
    structural_equation_systems.parameters.set<unsigned int>("eigenpairs")    = 10;
    structural_equation_systems.parameters.set<unsigned int>("basis vectors") = 10*3;
    
    // Initialize the data structures for the equation system.
    structural_equation_systems.init();
    
    // Prints information about the system to the screen.
    structural_equation_systems.print_info();
    
    // Pass the Dirichlet dof IDs to the CondensedEigenSystem
    std::set<unsigned int> dirichlet_dof_ids;
    get_dirichlet_dofs(structural_equation_systems,
                       "Eigensystem", dirichlet_dof_ids);
    structural_system.initialize_condensed_dofs(dirichlet_dof_ids);
    
    // Solve the system "Eigensystem".
    structural_system.solve();

    // Get the number of converged eigen pairs.
    unsigned int nconv = structural_system.get_n_converged();

    // create the solvers
    UGFlutterSolver flutter_solver;
    FlightCondition flight_cond;
    
    // attach the fluid and structural systems to the models
    CFDAerodynamicModel aero_model(fluid_system);
    FEMStructuralModel structural_model(structural_system);
    CoupledFluidStructureSystem
    coupled_system(aero_model, structural_model);
    
    ComplexMatrixX a;
    coupled_system.get_aero_operator_matrix(0.1, a);
}

#endif // LIBMESH_USE_COMPLEX_NUMBERS
