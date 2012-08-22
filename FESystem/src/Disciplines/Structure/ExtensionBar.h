//
//  ExtensionBar.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_Extension_Bar_h
#define FESystem_Extension_Bar_h

// FESystem includes
#include "Disciplines/Structure/Structural1DElementBase.h"


namespace FESystem
{
    namespace Geometry {class Point;}
    
    namespace Structures
    {
        class ExtensionBar: public FESystem::Structures::Structural1DElementBase
        {
        public:
            ExtensionBar();
            
            virtual ~ExtensionBar();
            
            virtual void clear();
            
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble area);

            virtual void calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            virtual void calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
        protected:
          
            void calculateOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat, FESystemBoolean if_strain);

            void getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            void getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            FESystemDouble area_val;
        };
    }
}


#endif // FESystem_Extension_Bar_h
