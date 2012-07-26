
/*
 *  Communicator.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/4/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_communicator_base_h__
#define __fesystem_communicator_base_h__

// FESystem includes
#include "Base/FESystemTypes.h"

// mpi includes
//#include "mpi.h"


namespace FESystem
{
	namespace Parallel
	{
		class CommunicatorBase
		{
		public:
			CommunicatorBase();
			
			virtual ~CommunicatorBase() = 0;
			
			FESystemInt getSize();
			
			FESystemInt getRank();
						
		protected:
			
		};
	}
}


#endif // __fesystem_communicator_base_h__

