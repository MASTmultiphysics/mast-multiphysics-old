//
//  LinearPlateElementBase.h
//  FESystem
//
//  Created by Manav Bhatia on 8/23/12.
//
//

#ifndef __FESystem__LinearPlateElementBase__
#define __FESystem__LinearPlateElementBase__

// FESystem includes
#include "Disciplines/Structure/Structural2DElementBase.h"


namespace FESystem
{
    namespace Structures
    {
        class LinearPlateElementBase: public FESystem::Structures::Structural2DElementBase
        {
        public:
            LinearPlateElementBase();
            
            virtual ~LinearPlateElementBase();
            
            virtual FESystemUInt getNElemDofs() const;
            
            /*!
             *    Returns the indices for the matrix entries that this element will contribute to. For instance, a bar element will only have the extensional stiffness
             *    matrix, and the indices will correspond to the location of the stiffness terms corresponding to the u-dofs.
             */
            virtual void getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec);
            
            virtual void transformMatrixToGlobalSystem(const FESystem::Numerics::MatrixBase<FESystemDouble>& elem_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& global_mat);

            virtual void transformVectorToGlobalSystem(const FESystem::Numerics::VectorBase<FESystemDouble> &elem_vec, FESystem::Numerics::VectorBase<FESystemDouble> &global_vec);

        protected:
            
        };
    }
}



#endif /* defined(__FESystem__LinearPlateElementBase__) */
