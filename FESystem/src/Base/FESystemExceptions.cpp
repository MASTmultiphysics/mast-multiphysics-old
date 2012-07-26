
/*
 *  FESystemExceptions.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/4/09.
 *  Copyright 2009 . All rights reserved.
 *
 */


// FESystem includes
#include "Base/FESystemExceptions.h"


FESystem::Exception::ExceptionBase::ExceptionBase():
std::exception(),
date(),
modification_date(),
compilation_date(),
file(),
line(),
function()
{}


FESystem::Exception::ExceptionBase::~ExceptionBase() throw()
{}
	

void
FESystem::Exception::ExceptionBase::printInfo(std::ostream& output) const
{
	
	output << this->getMessage() << std::endl
	<< "Tested Condition: " << this->condition << std::endl
	<< "In function: " << this->function << std::endl
	<< "Line: " << this->line << "  " << std::endl
	<< "File: " << this->file << std::endl
	<< "Compiled on : " << this->compilation_date << std::endl
	<< "Modified on : " << this->modification_date << std::endl;
}
	


void 
FESystem::Exception::ExceptionBase::setCondName(const std::string& cond)
{
	this->condition = cond;
}
	

const char*
FESystem::Exception::ExceptionBase::what() const throw ()
{
	std::ostringstream msg;
	msg << this->line << std::endl;
	return msg.str().c_str();
}


void 
FESystem::Exception::ExceptionBase::setExceptionData(std::string d, std::string t1, 
													 std::string t2, std::string f, 
													 FESystemUInt l, std::string func)
{
	this->date = d;
	this->modification_date = t1;
	this->compilation_date = t2;
	this->file = f;
	this->line = l;
	this->function = func;
}


