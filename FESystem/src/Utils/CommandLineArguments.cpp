
/*
 *  CommandLineArguments.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/11/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

// BOOST includes
#include <boost/algorithm/string.hpp>


// FESystem includes
#include "CommandLineArguments.h"


FESystem::Utility::CommandLineArguments::CommandLineArguments()
{
	
}


FESystem::Utility::CommandLineArguments::~CommandLineArguments()
{
	
}


void
FESystem::Utility::CommandLineArguments::appendArgument(std::string& str)
{
	// append this to the combined string
	this->combined_string += " ";
	this->combined_string += str;

	boost::algorithm::trim(this->combined_string);
		
	// now split this and overwrite the previous split result
	boost::algorithm::split(this->argument_list, this->combined_string, 
							boost::algorithm::is_space(), 
							boost::algorithm::token_compress_on);
	
	// also update the list of pointers
	FESystemUInt num = this->argument_list.size();
	this->argument_pointer.resize(num);
	for (FESystemUInt i=0 ; i < num; i++) 
		this->argument_pointer[i] = &(this->argument_list[i][0]);
}



const std::string&
FESystem::Utility::CommandLineArguments::getArguments() const
{
	return this->combined_string;
}



void
FESystem::Utility::CommandLineArguments::printArguments(std::ostream& output) const
{
	std::vector<std::string>::const_iterator it, end;
	it = this->argument_list.begin();
	end = this->argument_list.end();
	for ( ; it != end; it++) 
		output << *it << " ";
	output << std::endl;
}



const std::pair<FESystemUInt, char* const*> 
FESystem::Utility::CommandLineArguments::getArgumentsPair() const
{
	return std::make_pair<FESystemUInt, char* const*> (this->argument_pointer.size(),
													   &(this->argument_pointer[0]));
}



FESystemUInt
FESystem::Utility::CommandLineArguments::getNumArguments() const
{
	return this->argument_pointer.size();
}


