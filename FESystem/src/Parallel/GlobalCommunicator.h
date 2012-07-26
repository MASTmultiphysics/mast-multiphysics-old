
/*
 *  GlobalCommunicator.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/4/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_global_communicator_h__
#define __fesystem_global_communicator_h__

// FESystem includes
#include "Parallel/CommunicatorBase.h"


namespace FESystem 
{
	namespace Base
	{
		class FESystemLibrary;
	}
	
	namespace Parallel 
	{
		class GlobalCommunicator : public FESystem::Parallel::CommunicatorBase 
		{
		public:
			
			virtual ~GlobalCommunicator();
			
		private:

			GlobalCommunicator();
			
			//MPI_Comm  mpi_comm_world;
			
			friend class FESystem::Base::FESystemLibrary;
		};
	}
}


#endif // __fesystem_global_communicator_h__

