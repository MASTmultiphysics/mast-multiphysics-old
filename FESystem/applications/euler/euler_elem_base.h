//
//  euler_elem_base.h
//  FESystem
//
//  Created by Manav Bhatia on 3/6/13.
//
//

#ifndef __FESystem__euler_elem_base__
#define __FESystem__euler_elem_base__

// C++ includes
#include <ostream>


// libMesh includes
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/fem_context.h"
#include "libmesh/system.h"
#include "libmesh/analytic_function.h"
#include "libmesh/parallel.h"


using namespace libMesh;

enum FluidPrimitiveVars {RHO_PRIM, VEL1, VEL2, VEL3, TEMP };
enum FluidConservativeVars {RHO_CONS, RHOVEL1, RHOVEL2, RHOVEL3, ETOT };


class PrimitiveSolution
{
public:
    PrimitiveSolution();
    
    void zero();
    void init(const unsigned int dim, const DenseVector<Real>& conservative_sol, const Real cp_val, const Real cv_val);
    void print(std::ostream& out) const;
    
    Real c_pressure(const Real p0, const Real q0) const;
    
    DenseVector<Real> primitive_sol;
    unsigned int dimension;
    Real cp, cv;
    Real rho, u1, u2, u3, T, p, a, e_tot, k, entropy, mach;
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
    
    DenseVector<ValType> perturb_primitive_sol;
    ValType drho, du1, du2, du3, dT, dp, da, de_tot, dk, dentropy, dmach;
    
    const PrimitiveSolution* primitive_sol;
};




class FluidPostProcessSystem : public System
{
public:
    // Constructor
    FluidPostProcessSystem(EquationSystems& es,
                const std::string& name_in,
                const unsigned int number_in)
    : System(es, name_in, number_in)
    {}
    
    virtual void init_data();
    
    virtual void postprocess();

    unsigned int u, v, w, T, s, p, cp, a, M;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    unsigned int rho_im, u_im, v_im, w_im, T_im, s_im, p_im, cp_im, a_im, M_im;
#endif // LIBMESH_USE_COMPLEX_NUMBERS

};


// The Navier-Stokes system class.
// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class EulerElemBase
{
public:
    // Constructor
    EulerElemBase():
    system_comm(NULL),
    aoa(0.), rho_inf(0.), mach_inf(0.), temp_inf(0.), cp(0.), cv(0.), R(0.), gamma(0.), a_inf(0.), u1_inf(0.), u2_inf(0.), u3_inf(0.), q0_inf(0.), p_inf(0.)
    {}

    virtual ~EulerElemBase();

    
    void init_data();

    
    void get_infinity_vars( DenseVector<Real>& vars_inf ) const;
    

    void print_integrated_lift_drag(std::ostream& o);
    

    DenseVector<Number>* integrated_force;
    
    Real entropy_error, total_volume;

protected:
    
    
    void calculate_dxidX (const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c,
                          DenseMatrix<Real>& dxi_dX);
    
    
    void update_solution_at_quadrature_point( const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c, 
                                             const bool if_elem_domain, const DenseVector<Real>& elem_solution,
                                             DenseVector<Real>& conservative_sol, PrimitiveSolution& primitive_sol,
                                             DenseMatrix<Real>& B_mat);
    
    
    void update_jacobian_at_quadrature_point( const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c, PrimitiveSolution& primitive_sol,
                                             std::vector<DenseMatrix<Real> >& dB_mat,
                                             std::vector<DenseMatrix<Real> >& Ai_advection, DenseMatrix<Real>& A_entropy, DenseMatrix<Real>& A_inv_entropy );
    
    
    
    void calculate_advection_flux(const unsigned int calculate_dim,
                                  const PrimitiveSolution& sol,
                                  DenseVector<Real>& flux);
    
    void calculate_conservative_variable_jacobian(const PrimitiveSolution& sol,
                                                  DenseMatrix<Real>& dcons_dprim,
                                                  DenseMatrix<Real>& dprim_dcons);

    void calculate_advection_flux_jacobian(const unsigned int calculate_dim,
                                           const PrimitiveSolution& sol,
                                           DenseMatrix<Real>& mat);
    
    void calculate_advection_flux_jacobian_sensitivity_for_conservative_variable(const unsigned int calculate_dim,
                                                                                 const PrimitiveSolution& sol,
                                                                                 std::vector<DenseMatrix<Real> >& mat);
    
    
    void calculate_advection_flux_jacobian_sensitivity_for_primitive_variable(const unsigned int calculate_dim,
                                                                              const unsigned int primitive_var,
                                                                              const PrimitiveSolution& sol,
                                                                              DenseMatrix<Real>& mat);

    void calculate_advection_left_eigenvector_and_inverse_for_normal(const PrimitiveSolution& sol,
                                                                     const Point& normal, DenseMatrix<Real>& eig_vals,
                                                                     DenseMatrix<Real>& l_eig_mat, DenseMatrix<Real>& l_eig_mat_inv_tr);
    
    void calculate_entropy_variable_jacobian(const PrimitiveSolution& sol,
                                             DenseMatrix<Real>& dUdV, DenseMatrix<Real>& dVdU);
    
    
    void calculate_advection_flux_jacobian_for_moving_solid_wall_boundary(const PrimitiveSolution& sol,
                                                                          const Real xi_ni, const Point& nvec, DenseMatrix<Real>& mat);
    
    
    void calculate_advection_flux_spatial_derivative(const unsigned int i, DenseVector<Real>* flux, DenseMatrix<Real>* dflux_dU);
    
    
    
    void calculate_diffusive_flux_jacobian(unsigned int div_coord, unsigned int flux_coord, DenseMatrix<Real>& mat);
    
    
    void calculate_barth_tau_matrix(const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c,
                                    const PrimitiveSolution& sol,
                                    DenseMatrix<Real>& tau, std::vector<DenseMatrix<Real> >& tau_sens);
    
    void calculate_aliabadi_tau_matrix(const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c,
                                       const PrimitiveSolution& sol,
                                       DenseMatrix<Real>& tau, std::vector<DenseMatrix<Real> >& tau_sens);
    
    
    void calculate_differential_operator_matrix(const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c, const DenseVector<Real>& elem_solution, const PrimitiveSolution& sol,
                                                const DenseMatrix<Real>& B_mat, const std::vector<DenseMatrix<Real> >& dB_mat,
                                                const std::vector<DenseMatrix<Real> >& Ai_advection, const DenseMatrix<Real>& Ai_Bi_advection,
                                                const DenseMatrix<Real>& A_inv_entropy, const std::vector<std::vector<DenseMatrix<Real> > >& Ai_sens,
                                                DenseMatrix<Real>& LS_operator, DenseMatrix<Real>& LS_sens, Real& discontinuity_val);
    
    
    unsigned int dim;
    
    std::vector<FluidPrimitiveVars> _active_primitive_vars;

    std::vector<FluidConservativeVars> _active_conservative_vars;
    
    
    const Parallel::Communicator* system_comm;
    
public:

        
    Real aoa, rho_inf, mach_inf, temp_inf, cp, cv, R, gamma, a_inf, u1_inf, u2_inf, u3_inf, q0_inf, p_inf;
    
};



#endif /* defined(__FESystem__euler_elem_base__) */
