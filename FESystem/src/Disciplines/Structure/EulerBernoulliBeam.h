//
//  EulerBernoulliBeam.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_EulerBernoulliBeam_h
#define FESystem_EulerBernoulliBeam_h


// FESystem includes
#include  "Disciplines/Structure/LinearBeamElementBase.h"


namespace FESystem
{
    namespace Geometry {class Point;}
    
    namespace Structures
    {
        class EulerBernoulliBeam: public FESystem::Structures::LinearBeamElementBase
        {
        public:
            EulerBernoulliBeam();
            
            virtual ~EulerBernoulliBeam();
            
            virtual void clear();
            
            
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble I_tr, FESystemDouble I_ch, FESystemDouble A);
            
            virtual void calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            
            virtual void getStressTensor(const FESystem::Geometry::Point& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                         FESystem::Numerics::MatrixBase<FESystemDouble>& mat) ;
            
            
        protected:
            
            void getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& bend_mat);
            
            void calculateBendingOperatorMatrix(const FESystem::Geometry::Point& pt, FESystemDouble length, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);
        };
    }
}


#endif
