//
//  FluidElementBase.h
//  FESystem
//
//  Created by Manav Bhatia on 9/19/12.
//
//

#ifndef __FESystem__FluidElementBase__
#define __FESystem__FluidElementBase__


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
            
            
            void clear();
            
            
            void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                            FESystemDouble cp_val, FESystemDouble cv_val);
            
            void calculateFluxBoundaryCondition(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, const FESystem::Numerics::VectorBase<FESystemDouble>& mass_flux,
                                                const FESystem::Numerics::MatrixBase<FESystemDouble>& momentum_flux_tensor, const FESystem::Numerics::VectorBase<FESystemDouble>& energy_flux,
                                                FESystem::Numerics::VectorBase<FESystemDouble>& bc_vec);
            
            void calculateResidualVector(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, const FESystem::Numerics::VectorBase<FESystemDouble>& vel,
                                         FESystem::Numerics::VectorBase<FESystemDouble>& res);
            
            void calculateTangentMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, const FESystem::Numerics::VectorBase<FESystemDouble>& vel,
                                        FESystem::Numerics::MatrixBase<FESystemDouble>& dres_dx, FESystem::Numerics::MatrixBase<FESystemDouble>& dres_dxdot);
            
        protected:
            
            void calculateAdvectionFlux(const FESystemUInt i, FESystem::Numerics::VectorBase<FESystemDouble>& flux);
            
            void calculateAdvectionFluxSpatialDerivative(const FESystemUInt i, const FESystem::Numerics::MatrixBase<FESystemDouble>& dBmat_dx, FESystem::Numerics::VectorBase<FESystemDouble>* flux, FESystem::Numerics::MatrixBase<FESystemDouble>* dflux_dU);
            
            void updateVariablesAtQuadraturePoint(const FESystem::Numerics::MatrixBase<FESystemDouble>& Bmat);
            
            void calculateConservationVariableJacobian(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            void calculateAdvectionFluxJacobian(FESystemUInt div_coord, FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            void calculateDiffusiveFluxJacobian(FESystemUInt div_coord, FESystemUInt flux_coord, FESystem::Numerics::MatrixBase<FESystemDouble>& mat);

            void calculateArtificialDiffusionOperator(FESystem::Numerics::MatrixBase<FESystemDouble>& streamline_operator, FESystem::Numerics::MatrixBase<FESystemDouble>& discontinuity_operator);
            
            void calculateDifferentialOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& mat);

            void calculateOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat, FESystemBoolean if_strain, FESystemUInt deriv_dim);
            
            
            FESystemBoolean if_initialized;
            
            const FESystem::Mesh::ElemBase* geometric_elem;
            
            const FESystem::Quadrature::QuadratureBase* quadrature;
            
            const FESystem::FiniteElement::FiniteElementBase* finite_element;
            
            FESystemBoolean if_include_diffusion_flux;
            
            FESystem::Numerics::VectorBase<FESystemDouble>* interpolated_sol;
            
            const FESystem::Numerics::VectorBase<FESystemDouble>* solution;
            
            const FESystem::Numerics::VectorBase<FESystemDouble>* velocity;
           
            // Fluid properties: given as user input
            FESystemDouble cp, cv, gamma, R, s0, p0, T0;
            
            // Fluid variables
            FESystemDouble mu, rho, p, T, u1, u2, u3, e_tot;
            
            // Derived properties
            FESystemDouble s, v, a, d, h, e, alpha_p, beta_T, gamma_bar, k, e_1, e_1_bar;
            
            // Derivatives
            FESystemDouble dp_drho_T, dp_dT_rho, de_drho_T, de_dT_rho, drho_dp_T, drho_dT_p, de_dp_T, de_dT_p;
            
            // Other intermediate constants
            FESystemDouble e_rho_1, e_rho_2, e_rho_3, e_rho_4, e_p_1, e_p_2, e_p_3, e_p_4, s1, s2, e_c_1, e_c_2, e_c_3, e_c_4;
        };
    }
}



#endif /* defined(__FESystem__FluidElementBase__) */
