/*
 *  VectorBase.cpp
 *  FESystem
 *
 *  Created by Manav Bhatia on 1/9/11.
 *  Copyright 2011 . All rights reserved.
 *
 */

// FESystem includes
#include "Numerics/VectorBase.h"
#include "Numerics/LocalVector.h"


template <typename ValType>
std::auto_ptr<FESystem::Numerics::VectorBase<ValType> >
FESystem::Numerics::VectorCreate(FESystem::Numerics::VectorType vec_type)
{
    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > vec;
    
    switch (vec_type) {
        case FESystem::Numerics::LOCAL_VECTOR:
            vec.reset(new FESystem::Numerics::LocalVector<ValType>());
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            break;
    }
    
    return vec;
}


/***************************************************************************************/
// Template instantiations for some generic classes

template std::auto_ptr<FESystem::Numerics::VectorBase<FESystemFloat> > FESystem::Numerics::VectorCreate(FESystem::Numerics::VectorType vec_type);
template std::auto_ptr<FESystem::Numerics::VectorBase<FESystemDouble> > FESystem::Numerics::VectorCreate(FESystem::Numerics::VectorType vec_type);
template std::auto_ptr<FESystem::Numerics::VectorBase<FESystemComplexFloat> > FESystem::Numerics::VectorCreate(FESystem::Numerics::VectorType vec_type);
template std::auto_ptr<FESystem::Numerics::VectorBase<FESystemComplexDouble> > FESystem::Numerics::VectorCreate(FESystem::Numerics::VectorType vec_type);


/***************************************************************************************/


