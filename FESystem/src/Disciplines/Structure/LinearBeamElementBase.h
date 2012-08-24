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
            
            LinearBeamElementBase();
            
            ~LinearBeamElementBase();
            
            FESystemBoolean ifIncludeChordwiseStiffness();

            virtual FESystemUInt getNElemDofs() const;
            
            /*!
             *    Returns the indices for the matrix entries that this element will contribute to. For instance, a bar element will only have the extensional stiffness
             *    matrix, and the indices will correspond to the location of the stiffness terms corresponding to the u-dofs.
             */
            virtual void getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec);

        protected:
            
            virtual void clear();
            
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble I_tr, FESystemDouble I_ch, FESystemBoolean if_chordwise);
            
            FESystemBoolean if_include_chordwise_stiffness;
            
            FESystemDouble I_tr_val;
            
            FESystemDouble I_ch_val;
        };
    
    }
}


#endif
