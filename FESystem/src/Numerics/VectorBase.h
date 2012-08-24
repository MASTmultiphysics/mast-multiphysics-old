
/*
 *  VectorBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/27/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_vector_base_h__
#define __fesystem_vector_base_h__

// C++ includes
#include <memory>
#include <vector>

// FESystem includes
#include "Base/FESystemTypes.h"
#include "Base/FESystemExceptions.h"


namespace FESystem
{
	namespace Numerics
	{
        enum VectorType
        {LOCAL_VECTOR};
        
        
		template <typename ValType>
		class VectorBase
		{
		public:
            
            virtual ~VectorBase(){}
            
            virtual FESystem::Numerics::VectorType getType() const = 0;
            
			virtual FESystemUInt getSize() const = 0;
			
            virtual void resize(FESystemUInt i) = 0; 
            
            virtual void copyVector (const VectorBase<ValType>& t)=0;
            
            virtual void copyVectorVals(const VectorBase<ValType>& t)=0;
            
            virtual void addSubVectorVals(FESystemUInt row1, FESystemUInt row2, 
                                          FESystemUInt v_row1, FESystemUInt v_row2, 
                                          const ValType f, const FESystem::Numerics::VectorBase<ValType>& v) = 0;

            virtual void setSubVectorVals(FESystemUInt row1, FESystemUInt row2, 
                                          FESystemUInt v_row1, FESystemUInt v_row2, 
                                          const FESystem::Numerics::VectorBase<ValType>& v) = 0;
            
            virtual void getSubVectorVals(FESystemUInt row1, FESystemUInt row2, 
                                          FESystemUInt v_row1, FESystemUInt v_row2, 
                                          FESystem::Numerics::VectorBase<ValType>& v) const = 0;

            /*!
             *   copies the values from the given vector to the location given in \p rows. It is assumed that the 
             *   indices are given in ascending order. 
             */
            virtual void setSubVectorValsFromIndices(const std::vector<FESystemUInt>& rows, const FESystem::Numerics::VectorBase<ValType>& vec) = 0;

            /*!
             *    copies the vecetor values from the rows of this vector given in \p rows into the given vector. It is assumed that the row indices are given in 
             *    ascending order. 
             */
            virtual void getSubVectorValsFromIndices(const std::vector<FESystemUInt>& rows, FESystem::Numerics::VectorBase<ValType>& vec) const =0;

            virtual const ValType& getVal(const FESystemUInt i) const = 0; 			
			
            virtual void zero() = 0;

            /*!
             *   Returns the sum of all values in the vector
             */
            virtual ValType getSum() const=0;

            /*!
             *   Returns the minimum value
             */
            virtual ValType getMinVal() const=0;
            
            /*!
             *   Returns the maximum value
             */
            virtual ValType getMaxVal() const=0;

            /*!
             *     Sets the \p i^th element of this vector to \p val
             */
			virtual void setVal(const FESystemUInt i, const ValType& val) = 0; 

            /*!
             *     Adds the \p i^th element of this vector to \p val
             */
			virtual void addVal(const FESystemUInt i, const ValType& val) = 0; 

            /*!
             *     Adds the elements of \p val to the elements of this vector identified by the indices in \p indices.
             */
			virtual void addVal(const std::vector<FESystemUInt>& indices, const FESystem::Numerics::VectorBase<ValType>& val) = 0;

            /*!
             *     Sets all elements of this vector to \p val
             */
			virtual void setAllVals(const ValType& val) = 0; 

            /*!
             *    scales all values with the factor t.
             */
			virtual  void scale(const ValType& t)=0;

            /*!
             *    scales elements n1 to n2 with factor t.
             */
            virtual  void scaleSubVector(FESystemUInt n1, FESystemUInt n2, const ValType& t)=0;
            
			virtual  void scaleToUnitLength()=0;

            virtual void initializeToRandomUnitVector(ValType low, ValType up)=0;
            
            virtual typename RealOperationType(ValType) getL1Norm() const=0;
            
            virtual typename RealOperationType(ValType) getL2Norm() const=0;

			virtual  void 
            elementwiseMultiply(const FESystem::Numerics::VectorBase<ValType>& t)=0;
            
			virtual  void 
            elementwiseMultiply(const FESystem::Numerics::VectorBase<ValType>& t,
                                FESystem::Numerics::VectorBase<ValType>& res)=0;
            
            virtual  ValType
            dotProduct(const FESystem::Numerics::VectorBase<ValType>& t) const=0;

            /*!
             *   calculates \f$ v2 = v0 \times v1 \f$ where v0 is this vector, and v1 and v2 are the given vectors in the input arguments
             */
            virtual  void
            crossProduct(const FESystem::Numerics::VectorBase<ValType>& v1, FESystem::Numerics::VectorBase<ValType>& v2) const=0;
                        
			virtual void add (ValType f, const VectorBase<ValType>& t)=0;
			
            virtual const ValType* getVectorValues() const = 0;
            
            virtual ValType* getVectorValues() = 0;
            
            virtual void write(std::ostream& out) const = 0;
            			
		protected:
			
			
		};
        
        template <typename ValType>
        std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > VectorCreate(FESystem::Numerics::VectorType vec_type);
        
	}
    
}


#endif // __fesystem_vector_base_h__
