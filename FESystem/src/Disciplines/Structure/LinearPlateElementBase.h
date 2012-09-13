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
    // Forward declerations
    namespace Geometry {class Point;}
    
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
            
            void calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
            void calculateDiagonalMassMatrix(FESystem::Numerics::VectorBase<FESystemDouble>& vec);
            
            virtual void calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat) = 0;

            virtual void calculateDistributedLoad(FESystemDouble p_val, FESystem::Numerics::VectorBase<FESystemDouble>& vec);
            
            
        protected:
            
            void calculateInertiaOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat);

            void getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
            
        };
    }
}



#endif /* defined(__FESystem__LinearPlateElementBase__) */
