//
//  TimoshenkoBeam.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_TimoshenkoBeam_h
#define FESystem_TimoshenkoBeam_h

// FESystem includes
#include  "Disciplines/Structure/LinearBeamElementBase.h"


namespace FESystem
{
    namespace Geometry {class Point;}
    
    namespace Structures
    {
        class TimoshenkoBeam: public FESystem::Structures::LinearBeamElementBase
        {
        public:
            TimoshenkoBeam();
            
            virtual ~TimoshenkoBeam();
            
            virtual void clear();
            
            
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_bend,
                                    const FESystem::Quadrature::QuadratureBase& q_shear, FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble I_tr, FESystemDouble I_ch,
                                    FESystemDouble A, FESystemBoolean if_lateral);
            
            virtual void calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            
            virtual void getStressTensor(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                         FESystem::Numerics::MatrixBase<FESystemDouble>& mat) ;
            
            
        protected:
            
            void getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& bend_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& shear_mat);
            
            void calculateBendingOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);
            
            void calculateShearOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);
            
            const FESystem::Quadrature::QuadratureBase*  quadrature_bending;
            
            const FESystem::Quadrature::QuadratureBase*  quadrature_shear;
            
            FESystemDouble kappa;
        };
    }
}


#endif
