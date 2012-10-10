//
//  TestingIncludes.h
//  FESystem
//
//  Created by Manav Bhatia on 9/26/12.
//
//

#ifndef FESystem_TestingIncludes_h
#define FESystem_TestingIncludes_h

// C++ includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory.h>
#include <cassert>

// FESystem includes
#include "Base/FESystemBase.h"
#include "Utils/AutoPtrTable.h"
// Plotting
#include "Plotting/PLPlot.h"
// geometry
#include "Geom/Point.h"
#include "Geom/RectangularCoordinateSystem.h"
// mesh
#include "Mesh/MeshBase.h"
#include "Mesh/ElementType.h"
#include "Mesh/Node.h"
#include "Mesh/EdgeElemBase.h"
#include "Mesh/Edge2.h"
#include "Mesh/Edge3.h"
#include "Mesh/Quad4.h"
#include "Mesh/Quad9.h"
#include "Mesh/Quad8.h"
#include "Mesh/Tri3.h"
#include "Mesh/Tri6.h"
#include "Mesh/Tri7.h"
// finite elements
#include "FiniteElems/FELagrange.h"
#include "Functions/FunctionMappingBase.h"
// quadrature
#include "Quadrature/TrapezoidQuadrature.h"
// degrees of freedom
#include "Base/DegreeOfFreedomMap.h"
#include "Base/DegreeOfFreedomUnit.h"
// numerics
#include "Numerics/DenseMatrix.h"
#include "Numerics/SparseMatrix.h"
#include "Numerics/SparsityPattern.h"
// solver
#include "Solvers/LinearSolvers/LapackLinearSolver.h"
#include "Solvers/LinearSolvers/LUFactorizationLinearSolver.h"
#include "Solvers/LinearSolvers/PseudoTimeSteppingLinearSolver.h"
#include "Solvers/LinearSolvers/QRFactorizationLinearSolver.h"
#include "Solvers/EigenSolvers/QRMethodLinearEigensolver.h"
#include "Solvers/EigenSolvers/LapackLinearEigenSolver.h"
#include "Solvers/EigenSolvers/ArpackLinearEigenSolver.h"
#include "Solvers/TransientSolvers/NewmarkTransientSolver.h"
#include "Solvers/TransientSolvers/ExplicitRungeKuttaTransientSolver.h"
#include "Solvers/NonlinearSolvers/NewtonIterationNonlinearSolver.h"
// output processors
#include "OutputProcessors/GmshOutputProcessor.h"
#include "OutputProcessors/VtkOutputProcessor.h"
#include "OutputProcessors/TecplotOutputProcessor.h"
// Structural elements
#include "Disciplines/Structure/ReissnerMindlinPlate.h"
#include "Disciplines/Structure/DKTPlate.h"
#include "Disciplines/Structure/EulerBernoulliBeam.h"
#include "Disciplines/Structure/TimoshenkoBeam.h"
#include "Disciplines/Structure/ExtensionBar.h"
#include "Disciplines/Structure/VonKarmanStrain1D.h"
#include "Disciplines/Structure/VonKarmanStrain2D.h"
#include "Disciplines/Structure/Membrane.h"
#include "Disciplines/Structure/PistonTheory1D.h"
#include "Disciplines/Structure/PistonTheory2D.h"
// Fluid Elements
#include "Disciplines/Fluid/FluidElementBase.h"

enum MeshType{RIGHT_DIAGONAL, LEFT_DIAGONAL, CROSS, INVALID_MESH};

void createLineMesh(FESystem::Mesh::ElementType elem_type, FESystem::Mesh::MeshBase& mesh, FESystem::Geometry::Point& origin,
                    FESystemUInt nx, FESystemDouble x_length, FESystemUInt& n_elem_nodes, MeshType m_type, FESystemBoolean local_cs_same_as_global);

void createPlaneMesh(FESystem::Mesh::ElementType elem_type, FESystem::Mesh::MeshBase& mesh, FESystem::Geometry::Point& origin,
                       FESystemUInt nx, FESystemUInt ny, FESystemDouble x_length, FESystemDouble y_length, FESystemUInt& n_elem_nodes,
                       MeshType m_type, FESystemBoolean local_cs_same_as_global);

void calculateBeamStructuralMatrices(FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                                     const FESystem::Mesh::MeshBase& mesh, FESystem::Numerics::VectorBase<FESystemDouble>& global_sol,
                                     FESystem::Numerics::VectorBase<FESystemDouble>& internal_force,
                                     FESystem::Numerics::VectorBase<FESystemDouble>& external_force,
                                     FESystem::Numerics::MatrixBase<FESystemDouble>& global_stiffness_mat,
                                     FESystem::Numerics::VectorBase<FESystemDouble>& global_mass_vec);

void staticAnalysis(FESystemUInt dim, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map, const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                    const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Numerics::MatrixBase<FESystemDouble>& stiff_mat,
                    const FESystem::Numerics::VectorBase<FESystemDouble>& rhs, FESystem::Numerics::VectorBase<FESystemDouble>& sol);

void nonlinearSolution(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const
                       FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                       const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                       const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs,
                       FESystem::Numerics::MatrixBase<FESystemDouble>& stiff_mat,
                       FESystem::Numerics::VectorBase<FESystemDouble>& rhs, FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                       void (*calculateStructuralMatrices)(FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes,
                                                           const FESystem::Base::DegreeOfFreedomMap& dof_map,
                                                           const FESystem::Mesh::MeshBase& mesh,
                                                           FESystem::Numerics::VectorBase<FESystemDouble>& global_sol,
                                                           FESystem::Numerics::VectorBase<FESystemDouble>& internal_force,
                                                           FESystem::Numerics::VectorBase<FESystemDouble>& external_force,
                                                           FESystem::Numerics::MatrixBase<FESystemDouble>& global_stiffness_mat,
                                                           FESystem::Numerics::VectorBase<FESystemDouble>& global_mass_vec));

void modalAnalysis(FESystemUInt dim, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map, const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                   const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Numerics::MatrixBase<FESystemDouble>& stiff_mat,
                   const FESystem::Numerics::VectorBase<FESystemDouble>& mass_vec, const FESystemUInt n_modes, std::vector<FESystemUInt>& sorted_ids, FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals,
                   FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec);

void transientAnalysis(FESystemUInt dim, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map, const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                       const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Numerics::MatrixBase<FESystemDouble>& stiff_mat,
                       const FESystem::Numerics::VectorBase<FESystemDouble>& mass_vec, const std::vector<FESystemUInt>& sorted_ids, const FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals,
                       const FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec);

void transientFlutterSolution(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes,
                              const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map, const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                              const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Numerics::MatrixBase<FESystemDouble>& stiff_mat,
                              const FESystem::Numerics::VectorBase<FESystemDouble>& mass_vec, const FESystemUInt n_modes, std::vector<FESystemUInt>& sorted_ids, FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals,
                              FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec,
                              void (*calculatePistonTheoryMatrices)(FESystemDouble mach, FESystemDouble rho, FESystemDouble gamma, FESystemDouble a_inf, FESystemDouble u_inf,
                                                                    FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                                                                    const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Mesh::MeshBase& mesh,
                                                                    FESystem::Numerics::VectorBase<FESystemDouble>& global_sol, FESystem::Numerics::VectorBase<FESystemDouble>& global_vel,
                                                                    FESystem::Numerics::VectorBase<FESystemDouble>& force, FESystem::Numerics::VectorBase<FESystemDouble>& generalized_force,
                                                                    FESystem::Numerics::MatrixBase<FESystemDouble>& aero_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_stiffness_mat,
                                                                    FESystem::Numerics::MatrixBase<FESystemDouble>& aero_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_damp_mat,
                                                                    const FESystemUInt n_modes, std::vector<FESystemUInt>& sorted_ids, FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals, FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec));


void calculatePlatePistonTheoryMatrices(FESystemDouble mach, FESystemDouble rho, FESystemDouble gamma, FESystemDouble a_inf, FESystemDouble u_inf,
                                        FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                                        const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Mesh::MeshBase& mesh,
                                        FESystem::Numerics::VectorBase<FESystemDouble>& global_sol, FESystem::Numerics::VectorBase<FESystemDouble>& global_vel,
                                        FESystem::Numerics::VectorBase<FESystemDouble>& force, FESystem::Numerics::VectorBase<FESystemDouble>& generalized_force,
                                        FESystem::Numerics::MatrixBase<FESystemDouble>& aero_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_stiffness_mat,
                                        FESystem::Numerics::MatrixBase<FESystemDouble>& aero_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_damp_mat,
                                        const FESystemUInt n_modes, std::vector<FESystemUInt>& sorted_ids, FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals, FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec);


void calculateBeamPistonTheoryMatrices(FESystemDouble mach, FESystemDouble rho, FESystemDouble gamma, FESystemDouble a_inf, FESystemDouble u_inf,
                                       FESystemBoolean if_nonlinear, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                                       const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs, const FESystem::Mesh::MeshBase& mesh,
                                       FESystem::Numerics::VectorBase<FESystemDouble>& global_sol, FESystem::Numerics::VectorBase<FESystemDouble>& global_vel,
                                       FESystem::Numerics::VectorBase<FESystemDouble>& force, FESystem::Numerics::VectorBase<FESystemDouble>& generalized_force,
                                       FESystem::Numerics::MatrixBase<FESystemDouble>& aero_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_stiffness_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_stiffness_mat,
                                       FESystem::Numerics::MatrixBase<FESystemDouble>& aero_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& reduced_damp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& generalized_damp_mat,
                                       const FESystemUInt n_modes, std::vector<FESystemUInt>& sorted_ids, FESystem::Numerics::VectorBase<FESystemDouble>& eig_vals, FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vec);

int test_ode_integration_second_order(int argc, char * const argv[]);

int test_ode_integration_first_order(int argc, char * const argv[]);

int plate_analysis_driver(int argc, char * const argv[]);

int beam_analysis_driver(int argc, char * const argv[]);

int euler_analysis_driver(int argc, char * const argv[]);

#endif
