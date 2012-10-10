
/*
 *  Counter.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/4/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_counter_h__
#define __fesystem_counter_h__

// C++ includes
#include <map>
#include <string>


// FESystem includes
#include "Base/FESystemTypes.h"



namespace FESystem
{
	namespace Base
	{
		// Forward declerations
		class FESystemLibrary;
		
		class Counter
		{
		public:
			
			virtual ~Counter();
			
			void incrementCounter(const std::string& name); 
			
			void decrementCounter(const std::string& name); 
			
		private:
			
			Counter(const std::string& class_name);
			
			static std::map<std::string, std::pair<FESystemInt, FESystemInt> >  counter_map;
			
			friend class FESystem::Base::FESystemLibrary;
			
		};
		
	}
	
}


#endif // __fesystem_counter_h__

