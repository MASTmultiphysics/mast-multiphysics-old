
/*
 *  FESystemExceptions.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/4/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_fesystem_exceptions_h__
#define __fesystem_fesystem_exceptions_h__

// FESystem includes
#include "Base/FESystemTypes.h"

// C++ includes
#include <exception>
#include <string>
#include <iostream>
#include <sstream>


namespace FESystem
{	
	namespace Exception
	{
		class ExceptionBase : public std::exception
		{
		public:
			ExceptionBase();

			virtual ~ExceptionBase() throw();
			
			void printInfo(std::ostream& output) const;
			
			virtual std::string getMessage() const = 0;
			
			virtual const char* what() const throw();
			
			void setCondName(const std::string& cond);
			
			void setExceptionData(std::string d, std::string t1, 
								  std::string t2, std::string f, 
								  FESystemUInt l, std::string func);
		protected:
			
			
			std::string date;
			
			std::string modification_date;
			
			std::string compilation_date;
			
			std::string file;
			
			FESystemUInt line;
			
			std::string function;
			
			std::string condition;
		};

	}
	
}


#define DeclareException0(name, message)                       \
class name: public FESystem::Exception::ExceptionBase          \
{                                                              \
public:                                                        \
name(): ExceptionBase() {} \
virtual ~name() throw() {}                                     \
virtual std::string getMessage() const                         \
{  std::ostringstream msg;                                     \
   msg message;                                                \
   return msg.str();                                           \
}                                                              \
}


#define DeclareException1(name, Type1, message)                \
class name: public FESystem::Exception::ExceptionBase          \
{                                                              \
public:                                                        \
name(const Type1& A1):                                         \
ExceptionBase(), Arg1(A1) {} \
virtual ~name() throw() {}                                     \
virtual std::string getMessage() const                         \
{  std::ostringstream msg;                                     \
   msg message;                                                \
   return msg.str();                                           \
}                                                              \
protected:                                                     \
Type1 Arg1;                                                    \
}



#define DeclareException2(name, Type1, Type2, message)         \
class name: public FESystem::Exception::ExceptionBase          \
{                                                              \
public:                                                        \
name(const Type1& A1, const Type2& A2):                        \
ExceptionBase(), Arg1(A1), Arg2(A2) {} \
virtual ~name() throw() {}                                     \
virtual std::string getMessage() const                         \
{  std::ostringstream msg;                                     \
msg message;                                                   \
return msg.str();                                              \
}                                                              \
protected:                                                     \
Type1 Arg1;                                                    \
Type2 Arg2;                                                    \
}



#define DeclareException3(name, Type1, Type2, Type3, message) \
class name: public FESystem::Exception::ExceptionBase          \
{                                                              \
public:                                                        \
name(const Type1& A1, const Type2& A2, const Type3& A3):                        \
ExceptionBase(), Arg1(A1), Arg2(A2), Arg3(A3) {} \
virtual ~name() throw() {}                                     \
virtual std::string getMessage() const                         \
{  std::ostringstream msg;                                     \
msg message;                                                   \
return msg.str();                                              \
}                                                              \
protected:                                                     \
Type1 Arg1;                                                    \
Type2 Arg2;                                                    \
Type3 Arg3;                                                    \
}


#define DeclareException4(name, Type1, Type2, Type3, Type4, message) \
class name: public FESystem::Exception::ExceptionBase          \
{                                                              \
public:                                                        \
name(const Type1& A1, const Type2& A2, const Type3& A3, const Type4& A4):\
ExceptionBase(), Arg1(A1), Arg2(A2), Arg3(A3), Arg4(A4) {} \
virtual ~name() throw() {}                                     \
virtual std::string getMessage() const                         \
{  std::ostringstream msg;                                     \
msg message;                                                   \
return msg.str();                                              \
}                                                              \
protected:                                                     \
Type1 Arg1;                                                    \
Type2 Arg2;                                                    \
Type3 Arg3;                                                    \
Type4 Arg4;                                                    \
}





#define FESystemStrictAssert0(cond, exception_name)            \
if (!(cond))                                                   \
{                                                              \
exception_name exc;                                            \
exc.setExceptionData(__DATE__,__TIMESTAMP__,__TIME__,__FILE__,__LINE__,__PRETTY_FUNCTION__); \
exc.setCondName(#cond);                                        \
exc.printInfo(std::cout);                                      \
throw exc;                                                     \
}


#define FESystemStrictAssert1(cond, exception_name, A1)        \
if (!(cond))                                                   \
{                                                              \
exception_name exc(A1);                                        \
exc.setExceptionData(__DATE__,__TIMESTAMP__,__TIME__,__FILE__,__LINE__,__PRETTY_FUNCTION__); \
exc.setCondName(#cond);                                        \
exc.printInfo(std::cout);                                      \
throw exc;                                                     \
}


#define FESystemStrictAssert2(cond, exception_name, A1, A2)    \
if (!(cond))                                                   \
{                                                              \
exception_name exc(A1, A2);                                    \
exc.setExceptionData(__DATE__,__TIMESTAMP__,__TIME__,__FILE__,__LINE__,__PRETTY_FUNCTION__); \
exc.setCondName(#cond);                                        \
exc.printInfo(std::cout);                                      \
throw exc;                                                     \
}


#define FESystemStrictAssert3(cond, exception_name, A1, A2, A3)    \
if (!(cond))                                                   \
{                                                              \
exception_name exc(A1, A2, A3);                                \
exc.setExceptionData(__DATE__,__TIMESTAMP__,__TIME__,__FILE__,__LINE__,__PRETTY_FUNCTION__); \
exc.setCondName(#cond);                                        \
exc.printInfo(std::cout);                                      \
throw exc;                                                     \
}


#define FESystemStrictAssert4(cond, exception_name, A1, A2, A3, A4)    \
if (!(cond))                                                   \
{                                                              \
exception_name exc(A1, A2, A3, A4);                            \
exc.setExceptionData(__DATE__,__TIMESTAMP__,__TIME__,__FILE__,__LINE__,__PRETTY_FUNCTION__); \
exc.setCondName(#cond);                                        \
exc.printInfo(std::cout);                                      \
throw exc;                                                     \
}


#ifdef FESYSTEM_DEBUG                                                   
#define FESystemAssert0(cond, exception_name)                  \
FESystemStrictAssert0(cond, exception_name)                    
#else
#define FESystemAssert0(cond, exception_name) {}
#endif


#ifdef FESYSTEM_DEBUG                                                   
#define FESystemAssert1(cond, exception_name, A1)              \
FESystemStrictAssert1(cond, exception_name, A1)                    
#else
#define FESystemAssert1(cond, exception_name, A1) {}
#endif


#ifdef FESYSTEM_DEBUG                                                   
#define FESystemAssert2(cond, exception_name, A1, A2)          \
FESystemStrictAssert2(cond, exception_name, A1, A2)                    
#else
#define FESystemAssert2(cond, exception_name, A1, A2) {}
#endif


#ifdef FESYSTEM_DEBUG                                                   
#define FESystemAssert3(cond, exception_name, A1, A2, A3)          \
FESystemStrictAssert3(cond, exception_name, A1, A2, A3)                    
#else
#define FESystemAssert3(cond, exception_name, A1, A2, A3) {}
#endif


#ifdef FESYSTEM_DEBUG                                                   
#define FESystemAssert4(cond, exception_name, A1, A2, A3, A4)          \
FESystemStrictAssert4(cond, exception_name, A1, A2, A3, A4)                    
#else
#define FESystemAssert4(cond, exception_name, A1, A2, A3, A4) {}
#endif


namespace FESystem
{
	namespace Exception
	{
		DeclareException0(NULLQuantity, 
						  << "NULL Pointer Found.");

		DeclareException0(InvalidFunctionCall, 
						  << "This function does not support any functionality. Use a different function instead.");

		DeclareException0(EnumNotHandled, 
						  << "Enumeration Not Handled.");

        DeclareException0(InvalidValue, 
						  << "Specified Value Is Invalid.");

        DeclareException0(InvalidState, 
						  << "Object In Invalid State Before Operation.");

		DeclareException1(EnumerationNotHandled, 
						  FESystemUInt, 
						  << "Enumeration Not Handled : " << Arg1);

        DeclareException1(InvalidID, 
						  FESystemUInt, 
						  << "Invalid ID Given : " << Arg1);

        DeclareException2(InvalidTag,
						  std::string,
						  std::string,
						  << "Invalid tag found. Found: " << Arg1
						  << "  Expected: " << Arg2 );

		DeclareException2(IndexOutOfBound,
						  FESystemUInt,
						  FESystemUInt,
						  << "Index out of bound. Found: " << Arg1 
						  << ". Should be less than: " << Arg2 );

		DeclareException2(DimensionsDoNotMatch, 
						  FESystemUInt,
						  FESystemUInt,
						  << "Dimension does not match. Found: " << Arg1 
						  << ". Should be: " << Arg2 );

        DeclareException2(PositiveDifferenceNeeded, 
						  FESystemUInt,
						  FESystemUInt,
						  << "Positive difference needed between two values. Val1: " << Arg1 
						  << "  should be less than Val2: " << Arg2 );


	}
}


#endif // __fesystem_fesystem_exceptions_h__


