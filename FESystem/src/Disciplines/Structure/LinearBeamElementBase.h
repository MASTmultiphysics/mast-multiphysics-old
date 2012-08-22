//
//  LinearBeamElementBase.h
//  FESystem
//
//  Created by Manav Bhatia on 8/18/12.
//
//

#ifndef FESystem_LinearBeamElementBase_h
#define FESystem_LinearBeamElementBase_h

// FESystem includes
#include "Disciplines/Structure/Structural1DElementBase.h"


namespace FESystem
{
    namespace Structures
    {
        class LinearBeamElementBase: public FESystem::Structures::Structural1DElementBase
        {
        public:
            
            LinearBeamElementBase():
            if_include_chordwise_stiffness(false),
            I_tr_val(0.0),
            I_ch_val(0.0)
            { }
            
            ~LinearBeamElementBase() {}
            
            FESystemBoolean ifIncludeChordwiseStiffness() {return this->if_include_chordwise_stiffness;}
            
        protected:
            
            virtual void clear()
            {
                this->if_include_chordwise_stiffness = false;
                this->I_tr_val = 0.0;
                this->I_ch_val = 0.0;
            }
            
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble I_tr, FESystemDouble I_ch, FESystemBoolean if_chordwise)
            {
                FESystem::Structures::StructuralElementBase::initialize(elem,fe,q_rule,E,nu,rho);
                this->if_include_chordwise_stiffness = if_include_chordwise_stiffness;
                this->I_tr_val = I_tr;
                this->I_ch_val = I_ch;
            }
            
            FESystemBoolean if_include_chordwise_stiffness;
            
            FESystemDouble I_tr_val;
            
            FESystemDouble I_ch_val;
        };
    
    }
}


#endif
