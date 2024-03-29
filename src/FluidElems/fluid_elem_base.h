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

#ifndef __MAST__euler_elem_base__
#define __MAST__euler_elem_base__

// C++ includes
#include <ostream>
#include <map>

// MAST includes
#include "Numerics/fem_operator_matrix.h"
#include "Flight/flight_condition.h"
#include "BoundaryConditions/boundary_surface_motion.h"

// libMesh includes
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/fem_context.h"
#include "libmesh/system.h"
#include "libmesh/analytic_function.h"
#include "libmesh/parallel.h"




// Forward declerations
class GetPot;

enum FluidPrimitiveVars
{ RHO_PRIM, VEL1, VEL2, VEL3, TEMP };
enum FluidConservativeVars
{ RHO_CONS, RHOVEL1, RHOVEL2, RHOVEL3, ETOT };
enum FluidBoundaryConditionType
{ NO_SLIP_WALL, SYMMETRY_WALL, SLIP_WALL, FAR_FIELD, EXHAUST, ISOTHERMAL, ADIABATIC };


class PrimitiveSolution
{
public:
    PrimitiveSolution();
    
    void zero();
    
    void init(const unsigned int dim,
              const DenseRealVector& conservative_sol,
              const Real cp_val,
              const Real cv_val,
              bool if_viscous);
    
    void print(std::ostream& out) const;
    
    Real c_pressure(const Real p0, const Real q0) const;
    
    void get_uvec(DenseRealVector& u) const;
    
    DenseRealVector primitive_sol;
    
    unsigned int dimension;
    
    Real cp, cv;
    
    Real rho, u1, u2, u3, T, p, a, e_tot, k, entropy, mach;
    
    // viscous quantities
    Real Pr, k_thermal, mu, lambda;
};


template <typename ValType>
class SmallPerturbationPrimitiveSolution
{
public:
    SmallPerturbationPrimitiveSolution();
    
    void zero();
    void init(const PrimitiveSolution& sol, const libMesh::DenseVector<ValType>& delta_sol);
    void print(std::ostream& out) const;
    
    ValType c_pressure(const Real q0) const;
    void get_duvec(libMesh::DenseVector<ValType>& du) const;
    
    libMesh::DenseVector<ValType> perturb_primitive_sol;
    ValType drho, du1, du2, du3, dT, dp, da, de_tot, dk, dentropy, dmach;
    
    const PrimitiveSolution* primitive_sol;
};






// The Navier-Stokes system class.
// libMesh::FEMSystem, TimeSolver and  libMesh::NewtonSolver will handle most tasks,
// but we must specify element residuals
class FluidElemBase
{
public:
    // Constructor
    FluidElemBase(GetPot& infile):
    _if_viscous(false), _if_full_linearization(false),
    _if_update_stabilization_per_quadrature_point(true),
    _include_pressure_switch(false),
    surface_motion(NULL),
    flight_condition(NULL),
    dim(0),
    _dissipation_scaling(1.),
    _infile(infile)
    { }

    virtual ~FluidElemBase();

    
    void init_data();

    
    void get_infinity_vars( DenseRealVector& vars_inf ) const;
    
    /*!
     *    This defines the surface motion for use with the nonlinear 
     *    fluid solver. This can be used to define either a time-dependent
     *    motion, or a steady-state motion.
     */
    MAST::SurfaceMotionBase* surface_motion;
    
    FlightCondition* flight_condition;

    unsigned int dim;
    
    bool if_viscous() const {
        return _if_viscous;
    }

    void calculate_dxidX (const std::vector<unsigned int>& vars,
                          const unsigned int qp, libMesh::FEMContext& c,
                          DenseRealMatrix& dxi_dX,
                          DenseRealMatrix& dX_dxi);
    
    
    void update_solution_at_quadrature_point
    ( const std::vector<unsigned int>& vars, const unsigned int qp,
     libMesh::FEMContext& c, const bool if_elem_domain,
     const DenseRealVector& elem_solution, DenseRealVector& conservative_sol,
     PrimitiveSolution& primitive_sol, FEMOperatorMatrix& B_mat,
     std::vector<FEMOperatorMatrix>& dB_mat);
    
    
    
    void calculate_advection_flux(const unsigned int calculate_dim,
                                  const PrimitiveSolution& sol,
                                  DenseRealVector& flux);

    void calculate_diffusion_flux(const unsigned int calculate_dim,
                                  const PrimitiveSolution& sol,
                                  const DenseRealMatrix& stress_tensor,
                                  const DenseRealVector& temp_gradient,
                                  DenseRealVector& flux);

    /*!
     *    calculates and returns the stress tensor in \p stress_tensor.
     */
    void calculate_diffusion_tensors(const DenseRealVector& elem_sol,
                                     const std::vector<FEMOperatorMatrix>& dB_mat,
                                     const DenseRealMatrix& dprim_dcons,
                                     const PrimitiveSolution& psol,
                                     DenseRealMatrix& stress_tensor,
                                     DenseRealVector& temp_gradient);

    void calculate_conservative_variable_jacobian(const PrimitiveSolution& sol,
                                                  DenseRealMatrix& dcons_dprim,
                                                  DenseRealMatrix& dprim_dcons);

    void calculate_advection_flux_jacobian(const unsigned int calculate_dim,
                                           const PrimitiveSolution& sol,
                                           DenseRealMatrix& mat);

    void calculate_advection_flux_jacobian_rho_derivative(const unsigned int calculate_dim,
                                                          const PrimitiveSolution& sol,
                                                          DenseRealMatrix& mat);

    void calculate_advection_flux_jacobian_u1_derivative(const unsigned int calculate_dim,
                                                         const PrimitiveSolution& sol,
                                                         DenseRealMatrix& mat);

    void calculate_advection_flux_jacobian_u2_derivative(const unsigned int calculate_dim,
                                                         const PrimitiveSolution& sol,
                                                         DenseRealMatrix& mat);

    void calculate_advection_flux_jacobian_u3_derivative(const unsigned int calculate_dim,
                                                         const PrimitiveSolution& sol,
                                                         DenseRealMatrix& mat);

    void calculate_advection_flux_jacobian_T_derivative(const unsigned int calculate_dim,
                                                        const PrimitiveSolution& sol,
                                                        DenseRealMatrix& mat);

    
    void calculate_diffusion_flux_jacobian(const unsigned int flux_dim,
                                           const unsigned int deriv_dim,
                                           const PrimitiveSolution& sol,
                                           DenseRealMatrix& mat);

    void calculate_advection_flux_jacobian_sensitivity_for_conservative_variable
    (const unsigned int calculate_dim,
     const PrimitiveSolution& sol,
     std::vector<DenseRealMatrix >& mat);
    
    
    void calculate_advection_flux_jacobian_sensitivity_for_primitive_variable
    (const unsigned int calculate_dim, const unsigned int primitive_var,
     const PrimitiveSolution& sol, DenseRealMatrix& mat);
    
    
    void calculate_advection_left_eigenvector_and_inverse_for_normal
    (const PrimitiveSolution& sol, const libMesh::Point& normal,
     DenseRealMatrix& eig_vals, DenseRealMatrix& l_eig_mat,
     DenseRealMatrix& l_eig_mat_inv_tr);

    
    void calculate_advection_left_eigenvector_and_inverse_rho_derivative_for_normal
    (const PrimitiveSolution& sol, const libMesh::Point& normal,
     DenseRealMatrix& eig_vals, DenseRealMatrix& l_eig_mat,
     DenseRealMatrix& l_eig_mat_inv_tr);

    
    void calculate_advection_left_eigenvector_and_inverse_u1_derivative_for_normal
    (const PrimitiveSolution& sol, const libMesh::Point& normal,
     DenseRealMatrix& eig_vals, DenseRealMatrix& l_eig_mat,
     DenseRealMatrix& l_eig_mat_inv_tr);

    
    void calculate_advection_left_eigenvector_and_inverse_u2_derivative_for_normal
    (const PrimitiveSolution& sol, const libMesh::Point& normal,
     DenseRealMatrix& eig_vals, DenseRealMatrix& l_eig_mat,
     DenseRealMatrix& l_eig_mat_inv_tr);

    
    void calculate_advection_left_eigenvector_and_inverse_u3_derivative_for_normal
    (const PrimitiveSolution& sol, const libMesh::Point& normal,
     DenseRealMatrix& eig_vals, DenseRealMatrix& l_eig_mat,
     DenseRealMatrix& l_eig_mat_inv_tr);

    
    void calculate_advection_left_eigenvector_and_inverse_T_derivative_for_normal
    (const PrimitiveSolution& sol, const libMesh::Point& normal,
     DenseRealMatrix& eig_vals, DenseRealMatrix& l_eig_mat,
     DenseRealMatrix& l_eig_mat_inv_tr);

    
    void calculate_entropy_variable_jacobian(const PrimitiveSolution& sol,
                                             DenseRealMatrix& dUdV,
                                             DenseRealMatrix& dVdU);
    
    
    void calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
    (const PrimitiveSolution& sol,
     const Real ui_ni, const libMesh::Point& nvec,
     const DenseRealVector& dnvec, DenseRealMatrix& mat);
    
    
    void calculate_advection_flux_spatial_derivative
    (const unsigned int i, DenseRealVector* flux,
     DenseRealMatrix* dflux_dU);
    
    
    
    void calculate_diffusive_flux_jacobian(unsigned int div_coord,
                                           unsigned int flux_coord,
                                           DenseRealMatrix& mat);
    
    
    bool calculate_barth_tau_matrix(const std::vector<unsigned int>& vars,
                                    const unsigned int qp, libMesh::FEMContext& c,
                                    const PrimitiveSolution& sol,
                                    DenseRealMatrix& tau,
                                    std::vector<DenseRealMatrix >& tau_sens);
    
    bool calculate_aliabadi_tau_matrix(const std::vector<unsigned int>& vars,
                                       const unsigned int qp, libMesh::FEMContext& c,
                                       const PrimitiveSolution& sol,
                                       DenseRealMatrix& tau,
                                       std::vector<DenseRealMatrix >& tau_sens);
    
    void calculate_hartmann_discontinuity_operator
    (const std::vector<unsigned int>& vars,
     const unsigned int qp,
     libMesh::FEMContext& c,
     const PrimitiveSolution& sol,
     const DenseRealVector& elem_solution,
     const std::vector<FEMOperatorMatrix>& dB_mat,
     const DenseRealMatrix& Ai_Bi_advection,
     libMesh::DenseVector<Real>& discontinuity_val);

    
    void calculate_aliabadi_discontinuity_operator
    (const std::vector<unsigned int>& vars,
     const unsigned int qp,
     libMesh::FEMContext& c,
     const PrimitiveSolution& sol,
     const DenseRealVector& elem_solution,
     const std::vector<FEMOperatorMatrix>& dB_mat,
     const DenseRealMatrix& Ai_Bi_advection,
     libMesh::DenseVector<Real>& discontinuity_val);


    template <typename ValType>
    void calculate_small_disturbance_aliabadi_discontinuity_operator
    (const std::vector<unsigned int>& vars,
     const unsigned int qp,
     libMesh::FEMContext& c,
     const PrimitiveSolution& sol,
     const SmallPerturbationPrimitiveSolution<ValType>& dsol,
     const DenseRealVector& elem_solution,
     const std::vector<FEMOperatorMatrix>& dB_mat,
     const DenseRealMatrix& Ai_Bi_advection,
     libMesh::DenseVector<Real>& discontinuity_val);
    
    
    void calculate_differential_operator_matrix
    (const std::vector<unsigned int>& vars, const unsigned int qp,
     libMesh::FEMContext& c, const DenseRealVector& elem_solution,
     const PrimitiveSolution& sol, const FEMOperatorMatrix& B_mat,
     const std::vector<FEMOperatorMatrix>& dB_mat,
     const std::vector<DenseRealMatrix >& Ai_advection,
     const DenseRealMatrix& Ai_Bi_advection,
     const std::vector<std::vector<DenseRealMatrix > >& Ai_sens,
     DenseRealMatrix& LS_operator, DenseRealMatrix& LS_sens);
    
protected:
    
    std::vector<FluidPrimitiveVars> _active_primitive_vars;

    std::vector<FluidConservativeVars> _active_conservative_vars;
    
    bool _if_viscous, _if_full_linearization,
    _if_update_stabilization_per_quadrature_point,
    _include_pressure_switch;
    
    Real _dissipation_scaling;
    
    GetPot& _infile;
    
    std::multimap<unsigned int, FluidBoundaryConditionType> _boundary_condition;

};



#endif /* defined(__MAST__euler_elem_base__) */
