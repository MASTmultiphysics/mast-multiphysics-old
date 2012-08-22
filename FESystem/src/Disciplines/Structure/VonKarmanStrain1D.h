//
//  VonKarmanStrain1D.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_VonKarmanStrain1D_h
#define FESystem_VonKarmanStrain1D_h

// FESystem includes
#include "Disciplines/Structure/Structural1DElementBase.h"


namespace FESystem
{
    // Forward declerations
    namespace Geometry {class Point;}
    
    namespace Structures
    {
        // Forward declerations
        class ExtensionBar;
        class LinearBeamElementBase;
        
        class VonKarmanStrain1D: public FESystem::Structures::Structural1DElementBase
        {
        public:
            VonKarmanStrain1D();
            
            virtual ~VonKarmanStrain1D();
            
            virtual void clear();
            
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystem::Structures::ExtensionBar& bar, FESystem::Structures::LinearBeamElementBase& beam);

            virtual void calculateTangentStiffnessMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& bar_sol, FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
                        
        protected:
            
            void calculateOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);
            
            FESystem::Structures::ExtensionBar* bar_elem;

            FESystem::Structures::LinearBeamElementBase* beam_elem;
        
        };
    }
}


#endif
