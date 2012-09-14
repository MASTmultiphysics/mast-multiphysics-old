//
//  EigenSolverBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 3/26/11.
//  Copyright 2011 . All rights reserved.
//

// C++ includes
#include <float.h>
#include <set>

// FESystem includes
#include "Solvers/EigenSolvers/LinearEigenSolverBase.h"
#include "Numerics/MatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::LinearEigenSolverBase():
n_max_iterations(1000),
shift_type(FESystem::EigenSolvers::NO_SHIFT),
shift_value(0),
matrices_are_set(false),
n_converged_eig_vals(0),
convergence_tolerance(MACHINE_EPSILON),
A_mat(NULL),
B_mat(NULL)
{
    
}


template <typename ValType>
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::~LinearEigenSolverBase()
{
    
}


template <typename ValType>
void
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::setEigenProblemType
(FESystem::EigenSolvers::LinearEigenProblemType p)
{
    this->problem_type = p;
}


template <typename ValType>
FESystem::EigenSolvers::LinearEigenProblemType
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getEigenProblemType() const
{
    return this->problem_type;
}


template <typename ValType>
void
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::setEigenSpectrumType
(FESystem::EigenSolvers::EigenSpectrumType s)
{
    this->spectrum = s;
}


template <typename ValType>
FESystem::EigenSolvers::EigenSpectrumType
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getEigenSpectrumType() const
{
    return this->spectrum;
}


template <typename ValType>
void
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::setEigenShiftType
(FESystem::EigenSolvers::EigenShiftType p)
{
    this->shift_type = p;
}


template <typename ValType>
FESystem::EigenSolvers::EigenShiftType
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getEigenShiftType() const
{
    return this->shift_type;
}



template <typename ValType>
void
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::setEigenShiftValue(ValType s)
{
    this->shift_value = s;
}


template <typename ValType>
ValType
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getEigenShiftValue() const
{
    return this->shift_value;
}


template <typename ValType>
void
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::setMaxIterations
(FESystemUInt n)
{
    this->n_max_iterations = n;
}


template <typename ValType>
FESystemUInt
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getMaxIterations() const
{
    return this->n_max_iterations;
}



template <typename ValType>
void
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::setConvergenceTolerance(FESystemDouble& val)
{
    this->convergence_tolerance = val;
}

template <typename ValType>
FESystemDouble
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getConvergenceTolerance() const
{
    return this->convergence_tolerance;
}



template <typename ValType>
void
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::setMatrix
(FESystem::Numerics::MatrixBase<ValType>* lhs,
 FESystem::Numerics::MatrixBase<ValType>* rhs)
{
    switch (this->problem_type)
    {
        case FESystem::EigenSolvers::HERMITIAN:
        case FESystem::EigenSolvers::NONHERMITIAN:
        {
            FESystemAssert0(lhs != NULL, FESystem::Exception::NULLQuantity);
            FESystemAssert0(rhs == NULL, FESystem::EigenSolvers::GeneralizedEigenproblemDoesNotHaveRHSMatrix);
            
            this->A_mat = lhs;
        }
            break;
            
        case FESystem::EigenSolvers::GENERALIZED_HERMITIAN:
        case FESystem::EigenSolvers::GENERALIZED_NONHERMITIAN:
        {
            FESystemAssert0(lhs != NULL, FESystem::Exception::NULLQuantity);
            FESystemAssert0(rhs != NULL, FESystem::Exception::NULLQuantity);
            
            this->A_mat = lhs;
            this->B_mat = rhs;
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    
    std::pair<FESystemUInt, FESystemUInt> s = this->getAMatrix().getSize();
    
    // now initialize the matrix and vector to store the eigenvectors and eigenvalues
    switch (this->problem_type)
    {
        case FESystem::EigenSolvers::HERMITIAN:
        case FESystem::EigenSolvers::GENERALIZED_HERMITIAN:
        {
            this->eig_vec_mat_right.reset(FESystem::Numerics::MatrixCreate<ValType>(FESystem::Numerics::LOCAL_DENSE_MATRIX).release()); // eigenvector storage
            this->eig_vec_mat_left.reset(FESystem::Numerics::MatrixCreate<ValType>(FESystem::Numerics::LOCAL_DENSE_MATRIX).release()); // eigenvector storage
            this->eig_val_vec.reset(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release()); // eigenvalue storage
            
            
            this->eig_vec_mat_right->resize(s.first, s.second);
            this->eig_vec_mat_left->resize(s.first, s.second);
            this->eig_val_vec->resize(s.first);
        }
            break;
            
        case FESystem::EigenSolvers::NONHERMITIAN:
        case FESystem::EigenSolvers::GENERALIZED_NONHERMITIAN:
        {
            this->eig_vec_mat_right_complex.reset(FESystem::Numerics::MatrixCreate<typename ComplexOperationType(ValType)>(FESystem::Numerics::LOCAL_DENSE_MATRIX).release()); // eigenvector storage
            this->eig_vec_mat_left_complex.reset(FESystem::Numerics::MatrixCreate<typename ComplexOperationType(ValType)>(FESystem::Numerics::LOCAL_DENSE_MATRIX).release()); // eigenvector storage
            this->eig_val_vec_complex.reset(FESystem::Numerics::VectorCreate<typename ComplexOperationType(ValType)>(FESystem::Numerics::LOCAL_VECTOR).release()); // eigenvalue storage
            
            
            this->eig_vec_mat_right_complex->resize(s.first, s.second);
            this->eig_vec_mat_left_complex->resize(s.first, s.second);
            this->eig_val_vec_complex->resize(s.first);
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    this->n_converged_eig_vals = 0; // set the counter to 0
    
    this->matrices_are_set = true;
}



template <typename ValType>
FESystem::Numerics::MatrixBase<ValType>&
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getAMatrix()
{
    FESystemAssert0(this->A_mat != NULL, FESystem::Exception::NULLQuantity);
    return *(this->A_mat);
}


template <typename ValType>
FESystem::Numerics::MatrixBase<ValType>&
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getBMatrix()
{
    FESystemAssert0(this->matrices_are_set, FESystem::EigenSolvers::MatrixNotSet);
    FESystemAssert0((this->getEigenProblemType() == FESystem::EigenSolvers::GENERALIZED_HERMITIAN) ||
                    (this->getEigenProblemType() == FESystem::EigenSolvers::GENERALIZED_NONHERMITIAN),
                    FESystem::EigenSolvers::GeneralizedEigenproblemDoesNotHaveRHSMatrix);
    
    FESystemAssert0(this->B_mat != NULL, FESystem::Exception::NULLQuantity);
    
    return *(this->B_mat);
}



template <typename ValType>
const FESystem::Numerics::MatrixBase<ValType>&
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getRightEigenVectorMatrix() const
{
    FESystemAssert0(this->eig_vec_mat_right.get() != NULL, FESystem::Exception::NULLQuantity);
    
    switch (this->getEigenProblemType()){
        case FESystem::EigenSolvers::HERMITIAN:
        case FESystem::EigenSolvers::GENERALIZED_HERMITIAN:
            return *(this->eig_vec_mat_right);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
}



template <typename ValType>
const FESystem::Numerics::MatrixBase<ValType>&
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getLeftEigenVectorMatrix() const
{
    FESystemAssert0(this->eig_vec_mat_left.get() != NULL, FESystem::Exception::NULLQuantity);
    
    switch (this->getEigenProblemType()){
        case FESystem::EigenSolvers::HERMITIAN:
        case FESystem::EigenSolvers::GENERALIZED_HERMITIAN:
            return *(this->eig_vec_mat_left);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
}



template <typename ValType>
const FESystem::Numerics::VectorBase<ValType>&
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getEigenValues() const
{
    FESystemAssert0(this->eig_val_vec.get() != NULL, FESystem::Exception::NULLQuantity);
    
    switch (this->getEigenProblemType()){
        case FESystem::EigenSolvers::HERMITIAN:
        case FESystem::EigenSolvers::GENERALIZED_HERMITIAN:
            return *(this->eig_val_vec);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
}




template <typename ValType>
const FESystem::Numerics::MatrixBase<typename ComplexOperationType(ValType)>&
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getRightComplexEigenVectorMatrix() const
{
    FESystemAssert0(this->eig_vec_mat_right_complex.get() != NULL, FESystem::Exception::NULLQuantity);
    
    switch (this->getEigenProblemType()){
        case FESystem::EigenSolvers::NONHERMITIAN:
        case FESystem::EigenSolvers::GENERALIZED_NONHERMITIAN:
            return *(this->eig_vec_mat_right_complex);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    
}



template <typename ValType>
const FESystem::Numerics::MatrixBase<typename ComplexOperationType(ValType)>&
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getLeftComplexEigenVectorMatrix() const
{
    FESystemAssert0(this->eig_vec_mat_left_complex.get() != NULL, FESystem::Exception::NULLQuantity);
    
    switch (this->getEigenProblemType()){
        case FESystem::EigenSolvers::NONHERMITIAN:
        case FESystem::EigenSolvers::GENERALIZED_NONHERMITIAN:
            return *(this->eig_vec_mat_left_complex);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    
}



template <typename ValType>
const FESystem::Numerics::VectorBase<typename ComplexOperationType(ValType)>&
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::getComplexEigenValues() const
{
    FESystemAssert0(this->eig_val_vec_complex.get() != NULL, FESystem::Exception::NULLQuantity);
    
    switch (this->getEigenProblemType()){
        case FESystem::EigenSolvers::NONHERMITIAN:
        case FESystem::EigenSolvers::GENERALIZED_NONHERMITIAN:
            return *(this->eig_val_vec_complex);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
}


template <typename ValType>
void
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::prepareSortingVector(FESystem::EigenSolvers::EigenValueSortingCriteria sort,
                                                                        std::vector<FESystemUInt>& sorted_ids) const
{
    FESystemAssert0(this->n_converged_eig_vals>0 , FESystem::Exception::InvalidState);
    
    sorted_ids.resize(this->n_converged_eig_vals);
    
    switch (this->getEigenProblemType())
    {
        case FESystem::EigenSolvers::HERMITIAN:
        case FESystem::EigenSolvers::GENERALIZED_HERMITIAN:
        {
            FESystemAssert0(this->eig_val_vec.get() != NULL, FESystem::Exception::NULLQuantity);
            FESystemAssert2(this->eig_val_vec->getSize() == this->n_converged_eig_vals, FESystem::Exception::DimensionsDoNotMatch, this->eig_val_vec->getSize(), this->n_converged_eig_vals);
            
            sorted_ids.resize(this->n_converged_eig_vals);
            std::set<FESystemUInt> sorted_id_set;
            
            switch (sort)
            {
                case FESystem::EigenSolvers::VALUE:
                {
                    FESystemBoolean found;
                    ValType minval, maxval, val=0.0;
                    FESystemUInt minvalid, maxvalid, n_sweeps=0;
                    // for each sweep, look at the ids that have not been sorted yet (i.e., not in sorted_id_set)
                    while (sorted_id_set.size() < this->n_converged_eig_vals)
                    {
                        minval=FESystem::Base::getMachineMax<typename RealOperationType(ValType)>(); maxval=-minval;
                        minvalid=0; maxvalid=0;
                        found = false;
                        
                        for (FESystemUInt i=0; i<this->n_converged_eig_vals; i++)
                            if (sorted_id_set.count(i) == 0)
                            {
                                val = this->eig_val_vec->getVal(i);
                                if (FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(val) < FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(minval))
                                {
                                    minval = val;
                                    minvalid = i;
                                    found = true;
                                }
                                if (FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(val) > FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(maxval))
                                {
                                    maxval = val;
                                    maxvalid = i;
                                    found = true;
                                }
                            }
                        // make sure something was found
                        FESystemAssert0(found, FESystem::Exception::InvalidValue);
                        
                        // now add the ids to the sorted vector
                        if (maxvalid == minvalid) // happens when the last value in the vector remains
                        {
                            sorted_id_set.insert(maxvalid);
                            sorted_ids[n_sweeps] = maxvalid;
                        }
                        else
                        {
                            sorted_ids[n_sweeps] = minvalid;
                            sorted_ids[this->n_converged_eig_vals-1-n_sweeps] = maxvalid;
                            sorted_id_set.insert(maxvalid);
                            sorted_id_set.insert(minvalid);
                        }
                        
                        // increment the sweep counter
                        n_sweeps++;
                    }
                }
                    break;
                    
                default:
                    // only the magnitude sorting is available for Herminitian problems
                    FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
            }
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
}




template <typename ValType>
void
FESystem::EigenSolvers::LinearEigenSolverBase<ValType>::completePowerIterations
(FESystem::Numerics::MatrixBase<ValType>& mat,
 FESystem::Numerics::VectorBase<ValType>& eig_vals,
 FESystem::Numerics::MatrixBase<ValType>& eig_vecs)
{
    FESystemUInt n = mat.getSize().second;
    
    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> >
    vec1(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release()),
    vec2(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release());
    
    vec1->resize(n);
    vec2->resize(n);
    
    // only the first eigenvalue is calculated
    // assuem an initial unit vector
    for (FESystemUInt i=0; i<n; i++)
        vec1->setVal(i, 1.0);
    
    typename RealOperationType(ValType) conv = 1.0e6;
    ValType eig0 = 0.0, eig1 = 0.0;
    FESystemBoolean v1_to_v2 = true;
    
    while (fabs(conv) >= this->getConvergenceTolerance())
    {
        if (v1_to_v2)
        {
            vec1->scaleToUnitLength();
            mat.rightVectorMultiply(*vec1, *vec2);  // v2 = A v1
            eig0 = vec2->dotProduct(*vec1); // (v1' A v1)
            v1_to_v2 = false;
        }
        else
        {
            vec2->scaleToUnitLength();
            mat.rightVectorMultiply(*vec2, *vec1);  // v1 = A v2
            eig0 = vec1->dotProduct(*vec2); // (v2' A v2) / (v2' v2)
            v1_to_v2 = true;
        }
        
        // convergence metric
        conv = FESystem::Base::magnitude<ValType, typename RealOperationType(ValType)>(eig0 - eig1);
        eig1 = eig0; // eig = v2' A v1 / (v2' v2)
    }
    
    eig_vals.setVal(0, eig1); // copy the eigenvalue
    if (v1_to_v2)
        eig_vecs.setColumnVals(0, 0, n-1, *vec2); // copy the eigenvector from v2
    else
        eig_vecs.setColumnVals(0, 0, n-1, *vec1); // copy the eigenvector from v1
}





/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::EigenSolvers::LinearEigenSolverBase);


/***************************************************************************************/

