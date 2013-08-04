//
//  euler_elem_base.h
//  MAST
//
//  Created by Manav Bhatia on 3/6/13.
//
//

#ifndef __FESystem__euler_elem_base__
#define __FESystem__euler_elem_base__

// C++ includes
#include <ostream>
#include <map>

// MAST includes
#include "Numerics/fem_operator_matrix.h"
#include "Flight/flight_condition.h"
#include "FluidElems/surface_motion.h"

// libMesh includes
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/fem_context.h"
#include "libmesh/system.h"
#include "libmesh/analytic_function.h"
#include "libmesh/parallel.h"


using namespace libMesh;

enum FluidPrimitiveVars
{ RHO_PRIM, VEL1, VEL2, VEL3, TEMP };
enum FluidConservativeVars
{ RHO_CONS, RHOVEL1, RHOVEL2, RHOVEL3, ETOT };
enum FluidBoundaryConditionType
{ NO_SLIP_WALL, SYMMETRY_WALL, SLIP_WALL, FAR_FIELD, ISOTHERMAL, ADIABATIC };


class PrimitiveSolution
{
public:
    PrimitiveSolution();
    
    void zero();
    void init(const unsigned int dim, const DenseVector<Real>& conservative_sol, const Real cp_val, const Real cv_val,
              bool if_viscous);
    void print(std::ostream& out) const;
    
    Real c_pressure(const Real p0, const Real q0) const;
    void get_uvec(DenseVector<Real>& u) const;
    
    DenseVector<Real> primitive_sol;
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
    void init(const PrimitiveSolution& sol, const DenseVector<ValType>& delta_sol);
    void print(std::ostream& out) const;
    
    ValType c_pressure(const Real q0) const;
    void get_duvec(DenseVector<ValType>& du) const;
    
    DenseVector<ValType> perturb_primitive_sol;
    ValType drho, du1, du2, du3, dT, dp, da, de_tot, dk, dentropy, dmach;
    
    const PrimitiveSolution* primitive_sol;
};






// The Navier-Stokes system class.
// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class EulerElemBase
{
public:
    // Constructor
    EulerElemBase():
    _if_viscous(false), _if_full_linearization(false),
    _if_update_stabilization_per_quadrature_point(true)
    {}

    virtual ~EulerElemBase();

    
    void init_data();

    
    void get_infinity_vars( DenseVector<Real>& vars_inf ) const;
    
    /*!
     *    This defines the surface motion for use with the nonlinear 
     *    fluid solver. This can be used to define either a time-dependent
     *    motion, or a steady-state motion.
     */
    AutoPtr<SurfaceMotionBase> surface_motion;
    
    FlightCondition* flight_condition;

protected:
    
    
    void calculate_dxidX (const std::vector<unsigned int>& vars,
                          const unsigned int qp, FEMContext& c,
                          DenseMatrix<Real>& dxi_dX,
                          DenseMatrix<Real>& dX_dxi);
    
    
    void update_solution_at_quadrature_point
    ( const std::vector<unsigned int>& vars, const unsigned int qp,
     FEMContext& c, const bool if_elem_domain,
     const DenseVector<Real>& elem_solution, DenseVector<Real>& conservative_sol,
     PrimitiveSolution& primitive_sol, FEMOperatorMatrix& B_mat,
     std::vector<FEMOperatorMatrix>& dB_mat);
    
    
    
    void calculate_advection_flux(const unsigned int calculate_dim,
                                  const PrimitiveSolution& sol,
                                  DenseVector<Real>& flux);

    void calculate_diffusion_flux(const unsigned int calculate_dim,
                                  const PrimitiveSolution& sol,
                                  const DenseMatrix<Real>& stress_tensor,
                                  const DenseVector<Real>& temp_gradient,
                                  DenseVector<Real>& flux);

    /*!
     *    calculates and returns the stress tensor in \p stress_tensor.
     */
    void calculate_diffusion_tensors(const DenseVector<Real>& elem_sol,
                                     const std::vector<FEMOperatorMatrix>& dB_mat,
                                     const DenseMatrix<Real>& dprim_dcons,
                                     const PrimitiveSolution& psol,
                                     DenseMatrix<Real>& stress_tensor,
                                     DenseVector<Real>& temp_gradient);

    void calculate_conservative_variable_jacobian(const PrimitiveSolution& sol,
                                                  DenseMatrix<Real>& dcons_dprim,
                                                  DenseMatrix<Real>& dprim_dcons);

    void calculate_advection_flux_jacobian(const unsigned int calculate_dim,
                                           const PrimitiveSolution& sol,
                                           DenseMatrix<Real>& mat);

    void calculate_diffusion_flux_jacobian(const unsigned int flux_dim,
                                           const unsigned int deriv_dim,
                                           const PrimitiveSolution& sol,
                                           DenseMatrix<Real>& mat);

    void calculate_advection_flux_jacobian_sensitivity_for_conservative_variable
    (const unsigned int calculate_dim,
     const PrimitiveSolution& sol,
     std::vector<DenseMatrix<Real> >& mat);
    
    
    void calculate_advection_flux_jacobian_sensitivity_for_primitive_variable
    (const unsigned int calculate_dim, const unsigned int primitive_var,
     const PrimitiveSolution& sol, DenseMatrix<Real>& mat);
    

    void calculate_advection_left_eigenvector_and_inverse_for_normal
    (const PrimitiveSolution& sol, const Point& normal,
     DenseMatrix<Real>& eig_vals, DenseMatrix<Real>& l_eig_mat,
     DenseMatrix<Real>& l_eig_mat_inv_tr);
    
    void calculate_entropy_variable_jacobian(const PrimitiveSolution& sol,
                                             DenseMatrix<Real>& dUdV,
                                             DenseMatrix<Real>& dVdU);
    
    
    void calculate_advection_flux_jacobian_for_moving_solid_wall_boundary
    (const PrimitiveSolution& sol,
     const Real ui_ni, const Point& nvec,
     const DenseVector<Real>& dnvec, DenseMatrix<Real>& mat);
    
    
    void calculate_advection_flux_spatial_derivative
    (const unsigned int i, DenseVector<Real>* flux,
     DenseMatrix<Real>* dflux_dU);
    
    
    
    void calculate_diffusive_flux_jacobian(unsigned int div_coord,
                                           unsigned int flux_coord,
                                           DenseMatrix<Real>& mat);
    
    
    bool calculate_barth_tau_matrix(const std::vector<unsigned int>& vars,
                                    const unsigned int qp, FEMContext& c,
                                    const PrimitiveSolution& sol,
                                    DenseMatrix<Real>& tau,
                                    std::vector<DenseMatrix<Real> >& tau_sens);
    
    bool calculate_aliabadi_tau_matrix(const std::vector<unsigned int>& vars,
                                       const unsigned int qp, FEMContext& c,
                                       const PrimitiveSolution& sol,
                                       DenseMatrix<Real>& tau,
                                       std::vector<DenseMatrix<Real> >& tau_sens);
    
    
    void calculate_differential_operator_matrix
    (const std::vector<unsigned int>& vars, const unsigned int qp,
     FEMContext& c, const DenseVector<Real>& elem_solution,
     const PrimitiveSolution& sol, const FEMOperatorMatrix& B_mat,
     const std::vector<FEMOperatorMatrix>& dB_mat,
     const std::vector<DenseMatrix<Real> >& Ai_advection,
     const DenseMatrix<Real>& Ai_Bi_advection,
     const DenseMatrix<Real>& A_inv_entropy,
     const std::vector<std::vector<DenseMatrix<Real> > >& Ai_sens,
     DenseMatrix<Real>& LS_operator, DenseMatrix<Real>& LS_sens,
     Real& discontinuity_val);
    
    
    unsigned int dim;
    
    std::vector<FluidPrimitiveVars> _active_primitive_vars;

    std::vector<FluidConservativeVars> _active_conservative_vars;
    
    bool _if_viscous, _if_full_linearization,
    _if_update_stabilization_per_quadrature_point;
;
    
    std::multimap<unsigned int, FluidBoundaryConditionType> _boundary_condition;

};



#endif /* defined(__FESystem__euler_elem_base__) */
