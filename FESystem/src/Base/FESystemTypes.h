
/*
 *  FESystemTypes.h
 *  FESystem
 *
 *  Created by Manav Bhatia on 11/16/09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#ifndef __fesystem_fesystem_types_h__
#define __fesystem_fesystem_types_h__

// C++ include
#include <float.h>
#include <complex>

// FESystem includes
#include  "Base/FESystemConfig.h"

// constants
#define MACHINE_EPSILON 1.0e-14

// define number types here
typedef unsigned int  FESystemUInt;
typedef int  FESystemInt;
typedef double FESystemDouble;
typedef float FESystemFloat;
typedef bool FESystemBoolean;
typedef std::complex<double> FESystemComplexDouble;
typedef std::complex<float> FESystemComplexFloat;


namespace FESystem
{
	namespace Base
	{
        template <typename ValTypeArg, typename ValTypeRet> ValTypeRet magnitude(const ValTypeArg& v);
        template <> FESystemFloat inline magnitude<FESystemFloat, FESystemFloat>(const FESystemFloat& v) {return fabs(v);}
        template <> FESystemDouble inline magnitude<FESystemDouble, FESystemDouble>(const FESystemDouble& v) {return fabs(v);}
        template <> FESystemFloat inline magnitude<FESystemComplexFloat, FESystemFloat>(const FESystemComplexFloat& v) {return abs(v);}
        template <> FESystemDouble inline magnitude<FESystemComplexDouble, FESystemDouble>(const FESystemComplexDouble& v) {return abs(v);}
        
        template <typename ValTypeArg, typename ValTypeRet> ValTypeRet real(const ValTypeArg& v);
        template <> FESystemFloat inline real(const FESystemFloat& v) {return v;}
        template <> FESystemDouble inline real(const FESystemDouble& v) {return v;}
        template <> FESystemFloat inline real(const FESystemComplexFloat& v) {return real(v);}
        template <> FESystemDouble inline real(const FESystemComplexDouble& v) {return real(v);}
        
        template <typename ValTypeArg, typename ValTypeRet> ValTypeRet imag(const ValTypeArg& v);
        template <> FESystemFloat inline imag(const FESystemFloat& v) {return 0.0;}
        template <> FESystemDouble inline imag(const FESystemDouble& v) {return 0.0;}
        template <> FESystemFloat inline imag(const FESystemComplexFloat& v) {return imag(v);}
        template <> FESystemDouble inline imag(const FESystemComplexDouble& v) {return imag(v);}

        template <typename ValTypeArg, typename ValTypeRet> ValTypeRet comparisonValue(const ValTypeArg& v);
        template <> FESystemFloat inline comparisonValue(const FESystemFloat& v) {return v;}
        template <> FESystemDouble inline comparisonValue(const FESystemDouble& v) {return v;}
        template <> FESystemFloat inline comparisonValue(const FESystemComplexFloat& v) {return abs(v);}
        template <> FESystemDouble inline comparisonValue(const FESystemComplexDouble& v) {return abs(v);}

        template <typename ValType> ValType getMachineMin();
        template <typename ValType> ValType getMachineMax();
        
        template <> FESystemFloat inline getMachineMin() {return FLT_MIN;}
        template <> FESystemDouble inline getMachineMin() {return DBL_MIN;}

        template <> FESystemFloat inline getMachineMax() {return FLT_MAX;}
        template <> FESystemDouble inline getMachineMax() {return DBL_MAX;}

		template <typename Type1, typename Type2>
		struct MultiTypeOperation{typedef void return_type;};
        
		template <typename Type1>
		struct RealValueType{typedef void return_type;};

		template <typename Type1>
		struct ComplexValueType{typedef void return_type;};
    }
    
}


#define TypeComparion(Type1, Type2, Type3)                                      \
namespace FESystem                                                              \
{                                                                               \
namespace Base                                                               \
{                                                                            \
template < >                                                              \
struct MultiTypeOperation<Type1, Type2>                                  \
{typedef Type3 return_type;};                                            \
\
template < >                                                              \
struct MultiTypeOperation<Type2, Type1>                                  \
{typedef Type3 return_type;};                                            \
}                                                                            \
}


#define RealValueReturnType(Type1, Type2)                              \
namespace FESystem                                                              \
{                                                                               \
namespace Base                                                               \
{                                                                            \
template < >                                                              \
struct RealValueType<Type1>                                        \
{typedef Type2 return_type;};                                            \
}                                                                            \
}

#define ComplexValueReturnType(Type1, Type2)                              \
namespace FESystem                                                              \
{                                                                               \
namespace Base                                                               \
{                                                                            \
template < >                                                              \
struct ComplexValueType<Type1>                                        \
{typedef Type2 return_type;};                                            \
}                                                                            \
}


TypeComparion(FESystemFloat, FESystemDouble, FESystemDouble);
TypeComparion(FESystemFloat, FESystemComplexFloat, FESystemComplexFloat);
TypeComparion(FESystemFloat, FESystemComplexDouble, FESystemComplexDouble);
TypeComparion(FESystemDouble, FESystemComplexFloat, FESystemComplexDouble);
TypeComparion(FESystemDouble, FESystemComplexDouble, FESystemComplexDouble);


RealValueReturnType(FESystemFloat, FESystemFloat);
RealValueReturnType(FESystemDouble, FESystemDouble);
RealValueReturnType(FESystemComplexFloat, FESystemFloat);
RealValueReturnType(FESystemComplexDouble, FESystemDouble);

ComplexValueReturnType(FESystemFloat, FESystemComplexFloat);
ComplexValueReturnType(FESystemDouble, FESystemComplexDouble);
ComplexValueReturnType(FESystemComplexFloat, FESystemComplexFloat);
ComplexValueReturnType(FESystemComplexDouble, FESystemComplexDouble);


#define ReturnType(Type1, Type2) \
FESystem::Base::MultiTypeOperation<Type1, Type2>::return_type

#define RealOperationType(Type1) \
FESystem::Base::RealValueType<Type1>::return_type

#define ComplexOperationType(Type1) \
FESystem::Base::ComplexValueType<Type1>::return_type


#endif // __fesystem_fesystem_types_h__


