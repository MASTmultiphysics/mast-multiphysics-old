
/*
 *  EnumName.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/11/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_enum_name_h__
#define __fesystem_enum_name_h__

// C++ includes
#include <string>

// FESystem includes
#include "Base/FESystemBase.h"

namespace FESystem 
{
		
	namespace Utility 
	{
		std::string libraryEnumToName(const FESystem::Base::LibraryEnumerations& lib);
		
		FESystem::Base::LibraryEnumerations libraryNameToEnum(const std::string& lib);
	}
}


#endif  // __fesystem_enum_name_h__

