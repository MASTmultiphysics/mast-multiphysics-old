//
//  LocalMatrixBase.h
//  FESystem
//
//  Created by Manav Bhatia on 6/13/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_local_matrix_base_h__
#define __fesystem_local_matrix_base_h__


// C++ includes
#include <set>


// FESystem includes
#include "Numerics/MatrixBase.h"


namespace FESystem
{
    namespace Numerics
    {
        
        template <typename ValType>
        class LocalMatrixBase: public FESystem::Numerics::MatrixBase<ValType>
        {
        public:
            
            LocalMatrixBase();
            
            virtual ~LocalMatrixBase();

            virtual void clear();
            
            const std::pair<FESystemUInt, FESystemUInt> getSize() const;

            void zero();
            
            void scale(const ValType& t);
            
            ValType getMinVal() const;

            ValType getMaxVal() const;

            ValType getRayleighQuotient(const FESystem::Numerics::VectorBase<ValType>& v) const;
            
            void normalizeRowBasis();
            
            void normalizeColumnBasis();
            
            typename RealOperationType(ValType) getL1Norm() const;
            
            typename RealOperationType(ValType) getFrobeniusNorm() const;
            
            ValType getDeterminant() const;
            
            ValType getCofactor(FESystemUInt first_column, const std::set<FESystemUInt>& rows) const;
            
            void getInverse(FESystem::Numerics::MatrixBase<ValType>& mat) const;

            ValType getRayleighQuotient(const FESystem::Numerics::VectorBase<ValType>& v);

            virtual const ValType* getMatrixValues() const;
            
            
            virtual ValType* getMatrixValues();
            
            virtual void write(std::ostream& out) const;
            
            
            virtual void writeDetailed(std::ostream& out) const;
            
        protected:

            void resizeValPtr(FESystemUInt n);

            void resizeValPtr(ValType* v, FESystemUInt n);

            /*!
             *   All values are stored locally for the local matrices
             */
            ValType* mat_vals;
         
            /*!
             *   Dimension of this vector
             */
            FESystemUInt n_vals;
            
            /*!
             *   True if this vector owns the pointer to the vector vec_vals, and created it. Used to track whether or not this needs to be deleted.
             */ 
            FESystemBoolean if_owns_pointer;

            /*!
             *   Indicates the size of this matrix
             */
            std::pair<FESystemUInt, FESystemUInt> dims;

        };
        
        // specializations of the output routines
        template <> void  FESystem::Numerics::LocalMatrixBase<FESystemFloat>::write(std::ostream &out) const;
        template <> void  FESystem::Numerics::LocalMatrixBase<FESystemDouble>::write(std::ostream &out) const;
        template <> void  FESystem::Numerics::LocalMatrixBase<FESystemComplexFloat>::write(std::ostream &out) const;
        template <> void  FESystem::Numerics::LocalMatrixBase<FESystemComplexDouble>::write(std::ostream &out) const;
    }
}



#endif // __fesystem_local_matrix_base_h__
