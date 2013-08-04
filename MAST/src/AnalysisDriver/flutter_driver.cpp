
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

#ifdef LIBMESH_USE_COMPLEX_NUMBERS

int main_flutter (int argc, char* const argv[])
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
    StructuralSystemBase & structural_system =
    structural_equation_systems.add_system<StructuralSystemBase>
    ("StructuralSystem");
    
    structural_system.time_solver =
    AutoPtr<TimeSolver>(new SteadySolver(structural_system));
    structural_equation_systems.init ();
    
    // print information for structural system
    structural_mesh.print_info();
    structural_equation_systems.print_info();
    
    
    
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
