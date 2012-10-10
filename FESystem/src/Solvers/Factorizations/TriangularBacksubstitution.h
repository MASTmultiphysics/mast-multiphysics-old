//
//  TriangularBacksubstitution.h
//  FESystem
//
//  Created by Manav Bhatia on 4/6/11.
//  Copyright 2011 . All rights reserved.
//

#ifndef __fesystem_triangular_backsubstitution_h__
#define __fesystem_triangular_backsubstitution_h__

namespace FESystem
{
    // Forward declerations
    namespace Numerics {template <typename ValType> class MatrixBase;}
    namespace Numerics {template <typename ValType> class VectorBase;}

    namespace FactorizationSolvers
    {
        
        /*!
         *    triangular matrix type enumeration
         */
        enum TriangularType
        {
            UPPER_TRIANGULAR,
            LOWER_TRIANGULAR,
            INVALID_TRIANGULAR_PART
        };
        
        /*!
         *    This class provides the methods to perform triangular backsubstitution, i.e. solve \f$ T x = b \f$
         *    where T is a triangular matrix (upper or lower), and b is a given RHS to solve for x. The 
         *    method setTriangularMatrixType is used to tell the solver about the nature of the matrix, whether it is
         *    upper or lower triangular. 
         */
        template <typename ValType> 
        class TriangularBacksubstitution
        {
        public:
            /*!
             *   Constructor
             */
            TriangularBacksubstitution();
            
            virtual ~TriangularBacksubstitution();
            
            /*!
             *   Clears the data structure
             */
            virtual void clear();
            
            void setTriangularMatrixType(FESystem::FactorizationSolvers::TriangularType t);

            FESystem::FactorizationSolvers::TriangularType getTriangularMatrixType();
            
            void setMatrix(const FESystem::Numerics::MatrixBase<ValType>& m);

            const FESystem::Numerics::MatrixBase<ValType>& getMatrix() const;

            void backSubstitute(const FESystem::Numerics::VectorBase<ValType>& v1,
                                FESystem::Numerics::VectorBase<ValType>& res);

            void backSubstitute(const FESystem::Numerics::MatrixBase<ValType>& m1,
                                FESystem::Numerics::MatrixBase<ValType>& res);
            
        protected:
            
            FESystem::FactorizationSolvers::TriangularType tri_type;
                        
            const FESystem::Numerics::MatrixBase<ValType>* mat;            
        };
    }
}


#endif // __fesystem_triangular_backsubstitution_h__
