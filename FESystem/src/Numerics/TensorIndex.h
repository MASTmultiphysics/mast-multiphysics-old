
/*
 *  TensorIndex.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/15/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_tensor_index_h__
#define __fesystem_tensor_index_h__

// FEsystem include
#include "Numerics/TensorIndexBase.h"


namespace FESystem
{
	namespace Numerics
	{
		template <FESystemUInt dim, FESystemUInt rank>
		class TensorIndex : public FESystem::Numerics::TensorIndexBase
		{
		public:
			TensorIndex(): FESystem::Numerics::TensorIndexBase(dim, rank){}
			
			virtual ~TensorIndex(){}
			
		protected:
			
		};
	}
}


#endif // __fesystem_tensor_index_h__

