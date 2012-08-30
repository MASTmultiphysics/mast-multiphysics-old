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
    
    // Forward declerations
    namespace  Geometry {class Point;}
    
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
            
            virtual void calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            virtual void calculateDiagonalMassMatrix(FESystem::Numerics::VectorBase<FESystemDouble>& vec);
            
            virtual void calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat) = 0;

        protected:
            
            virtual void clear();
            
            virtual void initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                    FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble I_tr, FESystemDouble I_ch,  FESystemDouble A, FESystemBoolean if_lateral);
            
            void getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            void calculateInertiaOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);
            
            FESystemBoolean if_include_lateral_stiffness;
            
            FESystemDouble I_tr_val;
            
            FESystemDouble I_ch_val;
        };
    
    }
}


#endif
