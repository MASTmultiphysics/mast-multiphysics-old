
/*
 *  LocalVector.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_local_vector_h__
#define __fesystem_local_vector_h__

// C++ includes
#include <vector>

// FESystem includes
#include "Numerics/VectorBase.h"

namespace FESystem
{
	namespace Numerics
	{
		template <typename ValType>
		class LocalVector: public FESystem::Numerics::VectorBase<ValType>
		{
		public: 
            /*!
             *    Default constructor. The user can use the resize() method to initialize the 
             *    vector dimension.
             */
			LocalVector();
			
            /*!
             *    This constructor uses the dimension to initialize the vector. 
             */
			LocalVector(FESystemUInt n);

			virtual ~LocalVector();
            
            virtual FESystem::Numerics::VectorType getType() const;
            
			virtual FESystemUInt getSize() const;
            
            /*!
             *   This initializes the vector and the vector will create an array of the specified size
             */
            virtual void resize(FESystemUInt i); 

            /*!
             *   The vector is associated to the provided array pointer of dimension i. The vector is not zeroed out in this case. 
             */
            virtual void resize(FESystemUInt i, ValType* vals); 

            /*!
             *   clears the data structure
             */
            void clear();
            
            virtual void copyVector(const VectorBase<ValType>& t);

            
            virtual void copyVectorVals(const VectorBase<ValType>& t);

            
            virtual void setSubVectorVals(FESystemUInt row1, FESystemUInt row2, 
                                          FESystemUInt v_row1, FESystemUInt v_row2, 
                                          const FESystem::Numerics::VectorBase<ValType>& v);
            
            virtual void addSubVectorVals(FESystemUInt row1, FESystemUInt row2, 
                                          FESystemUInt v_row1, FESystemUInt v_row2, 
                                          const ValType f, const FESystem::Numerics::VectorBase<ValType>& v);

            virtual void getSubVectorVals(FESystemUInt row1, FESystemUInt row2, 
                                          FESystemUInt v_row1, FESystemUInt v_row2, 
                                          FESystem::Numerics::VectorBase<ValType>& v) const;

            /*!
             *    copies the vecetor values from the given vector into the rows given in \p rows. It is assumed that the row indices are given in 
             *    ascending order. 
             */
            virtual void setSubVectorValsFromIndices(const std::vector<FESystemUInt>& rows, const FESystem::Numerics::VectorBase<ValType>& vec);
			
            /*!
             *    copies the vecetor values from the rows of this vector given in \p rows into the given vector. It is assumed that the row indices are given in 
             *    ascending order. 
             */
            virtual void getSubVectorValsFromIndices(const std::vector<FESystemUInt>& rows, FESystem::Numerics::VectorBase<ValType>& vec) const;

            virtual const ValType& getVal(const FESystemUInt i) const; 			
			
            virtual void zero();
            
            /*!
             *   Returns the sum of all values in the vector
             */
            virtual ValType getSum() const;

            /*!
             *   Returns the minimum value
             */
            virtual ValType getMinVal() const;
            
            /*!
             *   Returns the maximum value
             */
            virtual ValType getMaxVal() const;

            /*!
             *     Sets the \p i^th element of this vector to \p val
             */
			virtual void setVal(const FESystemUInt i, const ValType& val); 
            
            /*!
             *     Adds the \p i^th element of this vector to \p val
             */
			virtual void addVal(const FESystemUInt i, const ValType& val); 

            /*!
             *     Adds the elements of \p val to the elements of this vector identified by the indices in \p indices.
             */
			virtual void addVal(const std::vector<FESystemUInt>& indices, const FESystem::Numerics::VectorBase<ValType>& val);

            /*!
             *     Sets all elements of this vector to \p val
             */
			virtual void setAllVals(const ValType& val); 

			virtual void scale(const ValType& t);

            virtual void complexConjugate();
            
            /*!
             *    scales elements n1 to n2 with factor t.
             */
            virtual void scaleSubVector(FESystemUInt n1, FESystemUInt n2, const ValType& t);
			
            virtual void scaleToUnitLength();

            virtual void initializeToRandomUnitVector(ValType low, ValType up);

            virtual typename RealOperationType(ValType) getL1Norm() const;

            virtual typename RealOperationType(ValType) getL2Norm() const;

            virtual void elementwiseMultiply(const FESystem::Numerics::VectorBase<ValType>& t);
            
			virtual void elementwiseMultiply(const FESystem::Numerics::VectorBase<ValType>& t, FESystem::Numerics::VectorBase<ValType>& res);
            
            virtual ValType dotProduct(const FESystem::Numerics::VectorBase<ValType>& t) const;
            
            /*!
             *   calculates \f$ v2 = v0 \times v1 \f$ where v0 is this vector, and v1 and v2 are the given vectors in the input arguments
             */
            virtual void crossProduct(const FESystem::Numerics::VectorBase<ValType>& v1, FESystem::Numerics::VectorBase<ValType>& v2) const;
            
			virtual void add (ValType f, const VectorBase<ValType>& t);
			
            virtual const ValType* getVectorValues() const;
            
            virtual ValType* getVectorValues();
            
            virtual void write(std::ostream& out) const;
            
		protected:

            /*!
             *   pointer to array of values for this vector
             */
			ValType* vec_vals;
            
            /*!
             *   Dimension of this vector
             */
            FESystemUInt n_vals;
            
            /*!
             *   True if this vector owns the pointer to the vector vec_vals, and created it. Used to track whether or not this needs to be deleted.
             */ 
            FESystemBoolean if_owns_pointer;
            
		};
        
        
        // Specialization of complex conjugate method
        template <> void FESystem::Numerics::LocalVector<FESystemDouble>::complexConjugate();
        template <> void FESystem::Numerics::LocalVector<FESystemFloat>::complexConjugate();
        template <> void FESystem::Numerics::LocalVector<FESystemComplexDouble>::complexConjugate();
        template <> void FESystem::Numerics::LocalVector<FESystemComplexFloat>::complexConjugate();

        
        DeclareException0(VectorNotInitialized, 
                          << "Vector Not Initialized Before Usage\n");
        
        DeclareException2(VectorSizeMismatch, 
                          FESystemUInt, FESystemUInt, 
                          << "Vector Size Mismatch\n"
                          << "Vector 1 Dimension: " << Arg1 << "\n"
                          << "Vector 2 Dimension: " << Arg2); 
        
        DeclareException2(ElementLocationExceedesVectorSize, 
                          FESystemUInt, FESystemUInt, 
                          << "Element Location Exceedes Vector Size\n"
                          << "Matrix Dimension: " << Arg1 << "\n"
                          << "Element Location: " << Arg2 );
        
        
	}
}


#endif // __fesystem_local_vector_h__
