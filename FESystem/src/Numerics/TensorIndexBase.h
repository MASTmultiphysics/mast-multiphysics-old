
/*
 *  TensorIndex.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/15/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_tensor_index_base_h__
#define __fesystem_tensor_index_base_h__

// C++ includes 
#include <vector>

// FEsystem includes
#include "Base/FESystemTypes.h"


namespace FESystem
{
	namespace Numerics
	{
		class TensorIndexBase
		{
		public:
			TensorIndexBase(const FESystemUInt d,
							const FESystemUInt r);
			
			virtual ~TensorIndexBase() = 0;
			
			void reset(); 
			
			FESystemUInt getRank() const;
						
			FESystemUInt getDimension() const;
			
			void setIndex(const FESystemUInt r, const FESystemUInt d);

			void setIndex(const std::vector<FESystemUInt>& i);

			const std::vector<FESystemUInt>&  getIndex() const;
			
			FESystemUInt getIndex(const FESystemUInt r) const; 
			
		protected:
			
			const FESystemUInt dim;
			
			const FESystemUInt rank;
			
			std::vector<FESystemUInt> index;
		};
	}
}

#endif // __fesystem_tensor_index_base_h__
