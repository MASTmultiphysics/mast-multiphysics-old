
/*
 *  CommandLineArguments.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/11/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_command_line_arguments_h__
#define __fesystem_command_line_arguments_h__

// C++ includes
#include <string>
#include <vector>
#include <sstream>


// FESystem includes
#include "Base/FESystemTypes.h"

namespace FESystem
{
	namespace Utility
	{
		class CommandLineArguments
		{
		public:
			
			CommandLineArguments();
			
			~CommandLineArguments();
			
			void appendArgument(std::string& str);

			const std::string& getArguments() const;
			
			void printArguments(std::ostream& output) const;
			
			const std::pair<FESystemUInt, char* const*> getArgumentsPair() const;
			
			FESystemUInt getNumArguments() const;
			
		protected:
			
			std::string combined_string;
			
			std::vector< std::string > argument_list;

			std::vector<char*> argument_pointer;
		};
	}
}

#endif // __fesystem_command_line_arguments_h__

