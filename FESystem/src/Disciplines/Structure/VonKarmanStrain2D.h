//
//  VonKarmanStrain2D.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_VonKarmanStrain2D_h
#define FESystem_VonKarmanStrain2D_h

// FESystem includes
#include "Disciplines/Structure/Structural2DElementBase.h"


namespace FESystem
{
    // Forward declerations
    namespace Geometry {class Point;}    
    
    namespace Structures
    {
        // Forward declerations
        class Membrane;
        
        class VonKarmanStrain2D: public FESystem::Structures::Structural2DElementBase
        {
        public:
            VonKarmanStrain2D();
            
            virtual ~VonKarmanStrain2D();
            
            virtual void clear();
            
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystem::Structures::Membrane& membrane);
            
            virtual void calculateTangentStiffnessMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& membrane_sol, FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
        protected:
            
            void getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            void calculateOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);
          
            FESystem::Structures::Membrane* membrane_elem;
            
        };
    }
}


#endif

