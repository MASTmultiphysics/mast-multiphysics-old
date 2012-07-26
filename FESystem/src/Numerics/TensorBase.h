
/*
 *  TensorBase.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/15/09.
 *  Copyright 2009 . All rights reserved.
 *
 */


#ifndef __fesystem_tensor_base_h__
#define __fesystem_tensor_base_h__

// C++ includes
#include <vector>


// FEsystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{
	namespace Numerics
	{
		// Forward declerations
		class TensorIndexBase;
		
		template <typename ValType>
		class TensorBase 
		{
		public: 
			
			TensorBase(const FESystemUInt d, 
					   const FESystemUInt r);
			
			virtual ~TensorBase() = 0;
			
			FESystemUInt getRank() const;
			
			FESystemUInt getDimension() const;
			
			const ValType& getVal(const FESystem::Numerics::TensorIndexBase& index) const; 

			void setVal(const FESystem::Numerics::TensorIndexBase& index, 
						const ValType& val) const; 
			

			const  TensorBase<ValType> operator* (const ValType& t);

			const  TensorBase<ValType> operator* (const TensorBase<ValType>& t);

			const  TensorBase<ValType> operator+ (const TensorBase<ValType>& t);

			const  TensorBase<ValType> operator- (const TensorBase<ValType>& t);

			const  TensorBase<ValType> & operator*= (const ValType& t);

			const  TensorBase<ValType> & operator*= (const TensorBase<ValType>& t);

			const  TensorBase<ValType> & operator+= (const TensorBase<ValType>& t);

			const  TensorBase<ValType> & operator-= (const TensorBase<ValType>& t);
			
			template <typename ValType2>
			const  TensorBase< typename ReturnType(ValType, ValType2) > 
			operator* (const ValType2 t);

			template <typename ValType2>
			const  TensorBase< typename ReturnType(ValType, ValType2) > 
			operator* (const TensorBase<ValType2> t);

			template <typename ValType2>
			const  TensorBase< typename ReturnType(ValType, ValType2) > 
			operator+ (const TensorBase<ValType2> t);

			template <typename ValType2>
			const  TensorBase< typename ReturnType(ValType, ValType2) > 
			operator- (const TensorBase<ValType2> t);
						
		protected:
			
			const FESystemUInt dim;
			
			const FESystemUInt rank;
			
			std::vector<ValType> values;
			
		};
	}
		
}



#endif // __fesystem_tensor_base_h__
