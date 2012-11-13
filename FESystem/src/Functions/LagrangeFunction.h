//
//  LagrangeFunction.h
//  FESystem
//
//  Created by Manav Bhatia on 3/23/12.
//  Copyright (c) 2012. All rights reserved.
//

#ifndef __fesystem_lagrange_function_h__
#define __fesystem_lagrange_function_h__


// FESystem includes
#include "Functions/FunctionBase.h"
#include "Base/FESystemExceptions.h"
#include "Numerics/LocalVector.h"


namespace FESystem
{
    namespace Functions
    {
        /*!
         *   This class derives from the base class FunctionBase, and implements the Lagrange function 
         *   with abcissa and ordinate types as ValTypes. 
         */
        template <typename ValType>
        class LagrangeFunction: public FESystem::Numerics::FunctionBase<ValType, ValType>
        {
        public:
            /*!
             *    default constructor
             */
            LagrangeFunction();
            
            /*!
             *    destructor
             */
            virtual ~LagrangeFunction();
            
            /*!
             *   clears the data structure for this object.
             */ 
            virtual void clear();
 
            /*!
             *   This initializes the polyonomial using the given \p vals
             */ 
            void reinit(const FESystem::Numerics::VectorBase<ValType>& vals);
            
            /*!
             *   Returns the number of shape functions that will be evaluated in this Lagrange functions. This is 
             *   equal to the number of discrete points. 
             */ 
            FESystemUInt getNPoints() const;            
            
            
            /*!
             *   dummy function for LagrangeFunction since the evaluation needs the point for which the 
             *   function needs to be evaluated. Use the overloaded function instead
             */
            virtual ValType getFunctionValue(const ValType& abcissa) const;

            
            /*!
             *   @returns value of the function at the given abcissa 
             *   @param abcissa value at which the function should be evaluated
             *   @param point id of the point specified in vals for which this shape function will be evaluated
             */
            ValType getFunctionValue(const ValType& abcissa, const FESystemUInt& point) const;

            
            /*! 
             *   dummy function for LagrangeFunction since the evaluation needs the point for which the 
             *   function needs to be evaluated. Use the overloaded function instead
             */
            virtual ValType getFunctionDerivative(const ValType& abcissa, const FESystemUInt& order) const;
            
            
            /*!
             *   @returns derivative of the function at the given abcissa 
             *   @param abcissa value at which the function should be evaluated 
             *   @param order derivative  
             */
            ValType getFunctionDerivative(const ValType& abcissa, const FESystemUInt& point, const FESystemUInt& order) const;
            
        protected:
            
            /*!
             *   This calculates the order^th derivative of the Lagrangian (sub)polynomial where the numerator product terms are 
             *   given in \p vals
             */
            ValType calculateDerivative(const std::vector<ValType>& vals, const FESystemUInt& order) const;
            
            /*!
             *    stores the initialization state for the function
             */
            FESystemBoolean if_initialized;
            
            /*!
             *    The abcissa points used to calculate the Lagrange polynomial
             */ 
            FESystem::Numerics::VectorBase<ValType>* abcissa_vals;
            
            /*!
             *    The denominator for each term
             */ 
            FESystem::Numerics::VectorBase<ValType>* denominators;

        };

        /*!
         *   Checks whether any of the denominators in the Lagrange function are zero, which 
         *   implies coincident points. 
         */
        DeclareException0(CoincidentPoints, 
						  << "Singularity found in Lagrange function. Coincident points given.");

    }
}



template <typename ValType>
inline 
FESystem::Functions::LagrangeFunction<ValType>::LagrangeFunction():
FESystem::Numerics::FunctionBase<ValType, ValType>(),
if_initialized(false)
{
    this->abcissa_vals = new FESystem::Numerics::LocalVector<ValType>;
    this->denominators = new FESystem::Numerics::LocalVector<ValType>;
}


template <typename ValType>
inline 
FESystem::Functions::LagrangeFunction<ValType>::~LagrangeFunction()
{
    delete this->abcissa_vals;
    delete this->denominators;
}


template <typename ValType>
inline 
void 
FESystem::Functions::LagrangeFunction<ValType>::reinit(const FESystem::Numerics::VectorBase<ValType>& vals)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemUInt n = vals.getSize();
    this->abcissa_vals->resize(n);
    this->abcissa_vals->copyVectorVals(vals);

    // calculate the denominators for each term
    this->denominators->resize(n);
    this->denominators->setAllVals(1.0);
    // now iterate over all terms and calculate the products
    for (FESystemUInt i=0; i<n; i++)
        for (FESystemUInt j=0; j<n; j++)
            if (i != j)
                this->denominators->setVal(i, this->denominators->getVal(i)*(vals.getVal(j)-vals.getVal(i)));
    
    // now make sure that none of the values are zero
    for (FESystemUInt i=0; i<n; i++)
        if (FESystem::Base::magnitude<ValType, typename RealOperationType(ValType)>(this->denominators->getVal(i)) < MACHINE_EPSILON)
            FESystemAssert0(false, FESystem::Functions::CoincidentPoints);
    
    this->if_initialized = true;
}


template <typename ValType>
inline 
void
FESystem::Functions::LagrangeFunction<ValType>::clear() 
{
    this->abcissa_vals->zero();
    this->denominators->zero();
    
    this->if_initialized = false;
}


template <typename ValType>
inline 
FESystemUInt 
FESystem::Functions::LagrangeFunction<ValType>::getNPoints() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    return this->abcissa_vals->getSize();
}


template <typename ValType>
inline 
ValType
FESystem::Functions::LagrangeFunction<ValType>::getFunctionValue(const ValType& abcissa) const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}


template <typename ValType>
inline 
ValType
FESystem::Functions::LagrangeFunction<ValType>::getFunctionValue(const ValType& abcissa, const FESystemUInt& point) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    FESystemUInt n = this->abcissa_vals->getSize();
    
    // make sure that the point number is valid
    FESystemAssert2(point < n, FESystem::Exception::IndexOutOfBound, n, point);
    
    // calculate the numerator terms
    std::vector<ValType> sub_vals;
    if (sub_vals.size() != (n-1))
        sub_vals.resize(n-1);
    FESystemUInt n_sub = 0;
    for (FESystemUInt i=0; i<n; i++)
        if (i != point)
            sub_vals[n_sub++] = (this->abcissa_vals->getVal(i) - abcissa);

    ValType val = this->calculateDerivative(sub_vals,0); // calculate the zeroth order derivative
    
    return val/this->denominators->getVal(point);
}


template <typename ValType>
inline 
ValType
FESystem::Functions::LagrangeFunction<ValType>::getFunctionDerivative(const ValType& abcissa, const FESystemUInt& order) const
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}



template <typename ValType>
inline 
ValType
FESystem::Functions::LagrangeFunction<ValType>::getFunctionDerivative(const ValType& abcissa, const FESystemUInt& point, const FESystemUInt& order) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    FESystemUInt n = this->abcissa_vals->getSize();
    
    // make sure that the point number is valid
    FESystemAssert2(point < n, FESystem::Exception::IndexOutOfBound,n, point);
    
    // calculate the numerator terms
    std::vector<ValType> sub_vals;
    if (sub_vals.size() != n-1)
        sub_vals.resize(n-1);
    
    FESystemUInt n_sub = 0;
    for (FESystemUInt i=0; i<n; i++)
        if (i != point)
            sub_vals[n_sub++] = (this->abcissa_vals->getVal(i) - abcissa);
    
    ValType val = this->calculateDerivative(sub_vals, order); // calculate the zeroth order derivative
    
    return val/this->denominators->getVal(point);
}


template <typename ValType>
inline 
ValType
FESystem::Functions::LagrangeFunction<ValType>::calculateDerivative(const std::vector<ValType>& vals, const FESystemUInt& order) const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    // if the order of differentiation is greater than the number of terms in val, then the differentiation value is equal to 0.0.
    // if the order of differentiation, and the number of terms in val are equal to 1, then the differentiation value is equal to -1.0. 
    // Else, subdivide the number of terms into summation of N products, where N = vals.size()
    FESystemUInt n = vals.size(), n_sub=0;
    ValType dval = 1.0;

    if (order == 0) // return the product of the terms
        for (FESystemUInt i=0; i<n; i++)
            dval *= vals[i];
    else if (order > n)
        dval = 0.0;
    else if ((order == n) && (n==1)) 
        dval = -1.0;
    else 
    {
        dval = 0.0;
        // this will store the terms for each differential term, and there are n=vals.size() of those.
        std::vector<ValType> sub_vals;
        if (sub_vals.size() != n-1)
            sub_vals.resize(n-1);
        for (FESystemUInt i=0; i<n; i++)
        {
            n_sub = 0;
            for (FESystemUInt j=0; j<n; j++)
                if (i != j)
                    sub_vals[n_sub++] = vals[j];
            
            dval -= this->calculateDerivative(sub_vals, order-1);
        }
    }
    
    return dval;
}


#endif  //__fesystem_lagrange_function_h__

