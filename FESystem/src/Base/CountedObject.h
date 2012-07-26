
/*
 *  CountedObject.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/4/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_counted_object_h__
#define __fesystem_counted_object_h__

// FESystem includes
#include "Base/Counter.h"


namespace FESystem
{
	namespace Base
	{
		template<typename T>
		class CountedObject : public FESystem::Base::Counter 
		{
		public:
			virtual ~CountedObject();
			
		protected:
			CountedObject();
		};
		
	}
	
}


#endif // __fesystem_counted_object_h__

