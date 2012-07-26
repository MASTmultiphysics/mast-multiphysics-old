
/*
 *  EnumName.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/11/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// C++ include
#include <iostream>

// FESystem includes
#include "Utils/EnumName.h"


std::string 
FESystem::Utility::libraryEnumToName(const FESystem::Base::LibraryEnumerations& lib)
{
	std::string name;
	
	switch (lib) 
	{
		case FESystem::Base::FESystem:
			name = "FESystem";
			break;

		default:
			FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, lib);
			break;
	}
	
	return name;
}


FESystem::Base::LibraryEnumerations 
FESystem::Utility::libraryNameToEnum(const std::string& lib)
{

}



