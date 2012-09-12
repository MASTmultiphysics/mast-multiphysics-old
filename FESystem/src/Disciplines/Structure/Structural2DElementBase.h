//
//  Structural2DElementBase.h
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

#ifndef FESystem_Structural2DElementBase_h
#define FESystem_Structural2DElementBase_h


// FESystem includes
#include "Disciplines/Structure/StructuralElementBase.h"


namespace FESystem
{
    namespace Structures
    {
        
        class Structural2DElementBase: public FESystem::Structures::StructuralElementBase
        {
        public:
            Structural2DElementBase();
            
            virtual ~Structural2DElementBase();

            virtual void clear();
            
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble th);

            FESystemDouble getThickness() const;
            
            
        protected:
            
            FESystemDouble th_val;
        };
    }
}



#endif
