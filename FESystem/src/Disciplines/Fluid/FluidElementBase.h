//
//  FluidElementBase.h
//  FESystem
//
//  Created by Manav Bhatia on 9/19/12.
//
//

#ifndef __FESystem__FluidElementBase__
#define __FESystem__FluidElementBase__

// C++ includes
#include <vector>

// FESystem includes
#include "Base/FESystemTypes.h"



namespace FESystem
{
    // Forward declerations
    namespace Numerics {template<typename ValType> class VectorBase;}
    namespace Numerics {template<typename ValType> class MatrixBase;}
    namespace Quadrature {class QuadratureBase;}
    namespace FiniteElement {class FiniteElementBase;}
    namespace Mesh {class ElemBase;}
    namespace Geometry {class Point;}

    
    namespace Fluid
    {
        class FluidElementBase
        {
        public:
            FluidElementBase();
            
            ~FluidElementBase();
            
            
            virtual void clear();
            
            
            void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                            FESystemDouble dt_val, FESystemDouble cp_val, FESystemDouble cv_val, const FESystem::Numerics::VectorBase<FESystemDouble>& sol, const FESystem::Numerics::VectorBase<FESystemDouble>& vel,
                            FESystemBoolean if_update_dc, std::vector<FESystemDouble>& dc_vals_at_q_pts);
            
            void calculateMixedBoundaryCondition(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, const FESystem::Numerics::VectorBase<FESystemDouble>& U_vec,
                                                 FESystem::Numerics::VectorBase<FESystemDouble>& bc_vec);
            
            void calculateTangentMatrixForMixedBoundaryCondition(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, const FESystem::Numerics::VectorBase<FESystemDouble>& U_vec,
                                                                 FESystem::Numerics::MatrixBase<FESystemDouble>& bc_mat);
            
            void calculateFluxBoundaryCondition(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, const FESystem::Numerics::VectorBase<FESystemDouble>& mass_flux,
                                                const FESystem::Numerics::MatrixBase<FESystemDouble>& momentum_flux_tensor, const FESystem::Numerics::VectorBase<FESystemDouble>& energy_flux,
                                                FESystem::Numerics::VectorBase<FESystemDouble>& bc_vec);
            
            void calculateSolidWallFluxBoundaryCondition(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, FESystem::Numerics::VectorBase<FESystemDouble>& bc_vec);
            
            void calculateTangentMatrixForSolidWallFluxBoundaryCondition(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, FESystem::Numerics::MatrixBase<FESystemDouble>& jac);

            void calculateFluxBoundaryConditionUsingLocalSolution(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, FESystem::Numerics::VectorBase<FESystemDouble>& bc_vec);
            
            void calculateTangentMatrixForFluxBoundaryConditionUsingLocalSolution(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, FESystem::Numerics::MatrixBase<FESystemDouble>& jac);
            
            void calculateResidualVector(FESystem::Numerics::VectorBase<FESystemDouble>& res);
            
            void calculateTangentMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& dres_dx, FESystem::Numerics::MatrixBase<FESystemDouble>& dres_dxdot);
            
            void calculatePrimitiveVariableValues(const FESystem::Numerics::VectorBase<FESystemDouble>& conservative_sol, const FESystemDouble p0, FESystem::Numerics::VectorBase<FESystemDouble>& primitive_sol, FESystemDouble& press, FESystemDouble& entropy, FESystemDouble& mach, FESystemDouble& cp);
            
        protected:
            
            void calculateAdvectionFlux(const FESystemUInt i, FESystem::Numerics::VectorBase<FESystemDouble>& flux);
            
            void calculateAdvectionFluxSpatialDerivative(const FESystemUInt i, FESystem::Numerics::VectorBase<FESystemDouble>* flux, FESystem::Numerics::MatrixBase<FESystemDouble>* dflux_dU);
            
            void calculateEntropyVariableJacobian(FESystem::Numerics::MatrixBase<FESystemDouble>& dUdV, FESystem::Numerics::MatrixBase<FESystemDouble>& dVdU);

            void calculatePressureJacobianOnSolidWall(FESystem::Numerics::VectorBase<FESystemDouble>& dpdU);

            void calculateAdvectionFluxJacobian(FESystemUInt div_coord, FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            void calculateAdvectionLeftEigenvectorAndInverseForNormal(const FESystem::Numerics::VectorBase<FESystemDouble>& normal, FESystem::Numerics::MatrixBase<FESystemDouble>& eig_vals,
                                                                      FESystem::Numerics::MatrixBase<FESystemDouble>& l_eig_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& l_eig_mat_inv);
            
            void calculateDiffusiveFluxJacobian(FESystemUInt div_coord, FESystemUInt flux_coord, FESystem::Numerics::MatrixBase<FESystemDouble>& mat);

            void calculateArtificialDiffusionOperator(FESystem::Numerics::MatrixBase<FESystemDouble>& streamline_operator);
            
            void calculateDifferentialOperatorMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& LS_operator, FESystemDouble& discontinuity_val, FESystem::Numerics::VectorBase<FESystemDouble>& discont_operator_sens);

            void calculateOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat, FESystemBoolean if_strain, FESystemUInt deriv_dim);
            
            void calculateOperatorMatrixForBoundary(const FESystemUInt b_id, const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat, FESystemBoolean if_strain, FESystemUInt deriv_dim);
            
            void updateVariablesAtQuadraturePoint(const FESystem::Geometry::Point& pt);

            void updateVariablesAtQuadraturePointForBoundary(const FESystemUInt b_id, const FESystem::Geometry::Point& pt);

            void updateVariablesForInterpolationOperator(const FESystem::Numerics::MatrixBase<FESystemDouble>& Bmat);
            
            FESystemDouble estimateJacobianSpectralRadius(const FESystem::Numerics::MatrixBase<FESystemDouble>& mat);

            FESystemBoolean if_initialized;
            
            const FESystem::Mesh::ElemBase* geometric_elem;
            
            const FESystem::Quadrature::QuadratureBase* quadrature;
            
            const FESystem::FiniteElement::FiniteElementBase* finite_element;
            
            FESystemBoolean if_include_diffusion_flux;
            
            FESystem::Numerics::VectorBase<FESystemDouble>* interpolated_sol;

            FESystem::Numerics::VectorBase<FESystemDouble>* interpolated_vel;

            const FESystem::Numerics::VectorBase<FESystemDouble>* solution;

            const FESystem::Numerics::VectorBase<FESystemDouble>* velocity;

            FESystem::Numerics::MatrixBase<FESystemDouble>  *dX_dxi, *dxi_dX, *B_mat, *A_entropy, *A_inv_entropy, *Ai_Bi_advection;
            
            FESystem::Numerics::VectorBase<FESystemDouble> *N_vec;

            std::vector<FESystem::Numerics::VectorBase<FESystemDouble>*> N_vec_dx;

            std::vector<FESystem::Numerics::MatrixBase<FESystemDouble>*> B_mat_dxi;
            
            std::vector<FESystem::Numerics::MatrixBase<FESystemDouble>*> Ai_advection;

            // shape function Jacobian
            FESystemDouble jac;
            
            // Fluid properties: given as user input
            FESystemDouble dt, cp, cv, gamma, R, s0, p0, T0;
            
            // Fluid variables
            FESystemDouble rho, p, T, u1, u2, u3, e_tot, k, e, a;
            
            // discontinuity capturing operator value at quadrature points
            std::vector<FESystemDouble>* discontinuity_capturing_value;
            
            // flag to specifiy whether to use the provided discontinuity capturing values, or to calculate and update them
            FESystemBoolean if_update_discont_values;
        };
    }
}



#endif /* defined(__FESystem__FluidElementBase__) */
