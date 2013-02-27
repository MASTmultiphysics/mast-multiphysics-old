//
//  assembleEuler.h
//  FESystem
//
//  Created by Manav Bhatia on 2/21/13.
//
//

#ifndef __FESystem__assembleEuler__
#define __FESystem__assembleEuler__

// DiffSystem framework files
#include "libmesh/fem_system.h"

using namespace libMesh;

class PrimitiveSolution
{
public:
    PrimitiveSolution();
    
    void zero();
    void init(const unsigned int dim, const DenseVector<Real>& conservative_sol, const Real cp, const Real cv);
    void print(std::ostream& out);
    
    Real c_pressure(const Real p0, const Real q0);
    
    DenseVector<Real> primitive_sol;
    Real rho, u1, u2, u3, T, p, a, e_tot, k, entropy, mach;
};

// The Navier-Stokes system class.
// FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
// but we must specify element residuals
class EulerSystem : public FEMSystem
{
public:
    // Constructor
    EulerSystem(EquationSystems& es,
                const std::string& name_in,
                const unsigned int number_in)
    : FEMSystem(es, name_in, number_in),
    aoa(0.), rho_inf(0.), mach_inf(0.), temp_inf(0.), cp(0.), cv(0.), R(0.), gamma(0.), a_inf(0.), u1_inf(0.), u2_inf(0.), u3_inf(0.), q0_inf(0.), p_inf(0.)
    {}
    
    void init_data();
    
    virtual void init_context(DiffContext &context);
    
    virtual void postprocess();

    virtual bool element_time_derivative (bool request_jacobian,
                                          DiffContext &context);
    
    virtual bool side_constraint (bool request_jacobian,
                                     DiffContext &context);
    
    virtual bool mass_residual (bool request_jacobian,
                                DiffContext& context);
    
    
protected:
    
    void get_infinity_vars( DenseVector<Real>& vars_inf );
    
    
    void calculate_dxidX (const unsigned int qp, FEMContext& c,
                          DenseMatrix<Real>& dxi_dX);
    

    void update_solution_at_quadrature_point(const unsigned int qp, FEMContext& c, const bool if_elem_time_derivative,
                                             const bool if_elem_domain,
                                              DenseVector<Real>& conservative_sol, PrimitiveSolution& primitive_sol,
                                              DenseMatrix<Real>& B_mat);
    
    
    void update_jacobian_at_quadrature_point(const unsigned int qp, FEMContext& c, PrimitiveSolution& primitive_sol,
                                             std::vector<DenseMatrix<Real> >& dB_mat,
                                             std::vector<DenseMatrix<Real> >& Ai_advection, DenseMatrix<Real>& A_entropy, DenseMatrix<Real>& A_inv_entropy );
    
    
    
    void calculate_advection_flux(const unsigned int calculate_dim,
                                  const PrimitiveSolution& sol,
                                  DenseVector<Real>& flux);

    void calculate_advection_flux_jacobian(const unsigned int calculate_dim,
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

    
    void calculate_artificial_diffusion_operator(const unsigned int qp, FEMContext& c,
                                                 const PrimitiveSolution& sol,
                                                 DenseMatrix<Real>& streamline_operator);
    
    
    void calculate_differential_dperator_matrix(const unsigned int qp, FEMContext& c, const bool if_elem_time_derivative, const PrimitiveSolution& sol,
                                                const DenseMatrix<Real>& B_mat, const std::vector<DenseMatrix<Real> >& dB_mat,
                                                const std::vector<DenseMatrix<Real> >& Ai_advection, const DenseMatrix<Real>& Ai_Bi_advection,
                                                const DenseMatrix<Real>& A_inv_entropy,
                                                DenseMatrix<Real>& LS_operator, Real& discontinuity_val);
    
    unsigned int dim;
    
    Real aoa, rho_inf, mach_inf, temp_inf, cp, cv, R, gamma, a_inf, u1_inf, u2_inf, u3_inf, q0_inf, p_inf;
    
    std::vector<unsigned int> vars;
};

#endif /* defined(__FESystem__assembleEuler__) */
