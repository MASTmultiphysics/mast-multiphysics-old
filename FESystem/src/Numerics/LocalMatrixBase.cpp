//
//  LocalMatrixBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 6/13/12.
//  Copyright (c) 2012. All rights reserved.
//


// C++ includes
#include <iostream>
#include <iomanip>

// FESystem includes
#include "Numerics/LocalMatrixBase.h"
#include "Numerics/VectorBase.h"
#include "Solvers/Factorizations/HouseHolderTriangulation.h"
#include "Solvers/Factorizations/TriangularBacksubstitution.h"
#include "Base/macros.h"


template <typename ValType>
FESystem::Numerics::LocalMatrixBase<ValType>::LocalMatrixBase():
mat_vals(NULL),
n_vals(0),
if_owns_pointer(false)
{

    
} 


template <typename ValType>
FESystem::Numerics::LocalMatrixBase<ValType>::~LocalMatrixBase()
{
    this->clear();
} 




template <typename ValType>
void
FESystem::Numerics::LocalMatrixBase<ValType>::clear()
{
    if (this->if_owns_pointer)
        if (mat_vals != NULL)
            delete[] this->mat_vals;
    this->mat_vals = NULL;
    this->n_vals = 0;
    this->if_owns_pointer = false;
}



template <typename ValType>
void
FESystem::Numerics::LocalMatrixBase<ValType>::resizeValPtr(FESystemUInt n)
{
    if (this->if_owns_pointer)
    {
        if (this->n_vals != n)
        {
            delete[] this->mat_vals;
            this->mat_vals = new ValType[n];
        }
    }
    else // if it does not own the pointer, then reset the pointer to a data created locally
    {
        this->mat_vals = NULL;
        this->mat_vals = new ValType[n];
    }
    
    this->n_vals = n;
    this->if_owns_pointer = true;
    
    this->zero();
}



template <typename ValType>
void
FESystem::Numerics::LocalMatrixBase<ValType>::resizeValPtr(ValType* v, FESystemUInt n)
{
    this->clear();
    this->mat_vals = v;
    this->n_vals = n;
    this->if_owns_pointer = false;
}



template <typename ValType>
const std::pair<FESystemUInt, FESystemUInt>
FESystem::Numerics::LocalMatrixBase<ValType>::getSize() const
{
    return this->dims;
} 


template <typename ValType>
void 
FESystem::Numerics::LocalMatrixBase<ValType>::zero()
{
    for (FESystemUInt i=0; i<this->n_vals; i++)
        this->mat_vals[i] = 0.0;
}




template <typename ValType>
void
FESystem::Numerics::LocalMatrixBase<ValType>::scale(const ValType& t)
{
    for (FESystemUInt i=0; i < this->n_vals; i++)
        this->mat_vals[i] *= t;
}



template <typename ValType>
ValType
FESystem::Numerics::LocalMatrixBase<ValType>::getMinVal() const
{
    ValType v_min;
    typename RealOperationType(ValType) v = FESystem::Base::getMachineMax<typename RealOperationType(ValType)>();
    for (FESystemUInt i=0; i<this->n_vals; i++)
        if ( FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(this->mat_vals[i]) < v)
        {
            v_min = this->mat_vals[i];
            v = FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(v_min);
        }
    
    return v_min;
}



template <typename ValType>
ValType
FESystem::Numerics::LocalMatrixBase<ValType>::getMaxVal() const
{
    ValType v_max;
    typename RealOperationType(ValType) v = FESystem::Base::getMachineMin<typename RealOperationType(ValType)>();
    for (FESystemUInt i=0; i<this->n_vals; i++)
        if ( FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(this->mat_vals[i]) > v)
        {
            v_max = this->mat_vals[i];
            v = FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(v_max);
        }
    
    return v_max;
}



template <typename ValType>
ValType
FESystem::Numerics::LocalMatrixBase<ValType>::getRayleighQuotient(const FESystem::Numerics::VectorBase<ValType>& v)
{   
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    FESystemAssert0(s.first == s.second, FESystem::Numerics::MatrixMustBeSquareForOperation);
    
    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > 
    vec1(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release());
    
    vec1->resize(s.second);    
    
    this->rightVectorMultiply(v, *vec1); // vec1 =  A x
    ValType val = v.dotProduct(*vec1);  // x' vec1
    val /= pow(v.getL2Norm(),2);// (x' A x) / (x' x)
    return val; 
}



template <typename ValType>
void 
FESystem::Numerics::LocalMatrixBase<ValType>::normalizeRowBasis()
{
    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > vec(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR));
    
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    vec->resize(s.second);
    
    // for each row vector, get the values of the row, and scale it with its L2 norm
    for (FESystemUInt i=s.first; i<s.first; i++)
    {
        this->getRowVals(i, 0, s.second, *vec); 
        this->scaleRow(i, vec->getL2Norm());
    }
}



template <typename ValType>
void 
FESystem::Numerics::LocalMatrixBase<ValType>::normalizeColumnBasis()
{
    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > vec(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR));
    
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    vec->resize(s.first);
    
    // for each row vector, get the values of the row, and scale it with its L2 norm
    for (FESystemUInt i=s.second; i<s.second; i++)
    {
        this->getColumnVals(i, 0, s.first, *vec); 
        this->scaleColumn(i, vec->getL2Norm());
    }
}



template <typename ValType>
typename RealOperationType(ValType) 
FESystem::Numerics::LocalMatrixBase<ValType>::getL1Norm() const
{   
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > 
    vec1(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release());
    
    vec1->resize(s.first); 
    FESystemDouble norm = 0.0, val1=0.0;
    
    for (FESystemUInt i=0; i<s.second; i++) {
        this->getColumnVals(i, 0, s.second-1, *vec1);
        val1 = vec1->getL1Norm();
        if (val1 > norm)
            norm = val1;
    }
    
    return norm; 
}



template <typename ValType>
typename RealOperationType(ValType) 
FESystem::Numerics::LocalMatrixBase<ValType>::getFrobeniusNorm() const
{   
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    std::auto_ptr<FESystem::Numerics::VectorBase<ValType> > 
    vec1(FESystem::Numerics::VectorCreate<ValType>(FESystem::Numerics::LOCAL_VECTOR).release());
    
    vec1->resize(s.first); 
    FESystemDouble norm = 0.0;
    
    for (FESystemUInt i=0; i<s.second; i++) {
        this->getColumnVals(i, 0, s.second-1, *vec1);
        norm += pow(vec1->getL2Norm(),2);
    }
    
    norm = sqrt(norm);
    
    return norm; 
}





template <typename ValType>
ValType
FESystem::Numerics::LocalMatrixBase<ValType>::getDeterminant() const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    // make sure that the matrix is square
    FESystemAssert0(s.first == s.second, FESystem::Numerics::MatrixMustBeSquareForOperation);
    
    // specialize for 1, 2 and 3-D matrices and if the dimension is greater than that, use the QR method
    if (s.first == 1)
        return this->getVal(0,0);
    else if (s.first == 2)
        return (this->getVal(0,0)*this->getVal(1,1)-this->getVal(0,1)*this->getVal(1,0));
    else if (s.first == 3)
    {
        ValType v00=this->getVal(0,0),v01=this->getVal(0,1),v02=this->getVal(0,2),
        v10=this->getVal(1,0),v11=this->getVal(1,1),v12=this->getVal(1,2),
        v20=this->getVal(2,0),v21=this->getVal(2,1),v22=this->getVal(2,2),
        det=-v02*v11*v20 + v01*v12*v20 + v02*v10*v21 - v00*v12*v21 - v01*v10*v22 + v00*v11*v22;
        return det;
    }
    else
    {
        std::set<FESystemUInt> included_rows;
        // initialize the ids
        for (FESystemUInt i=0; i<s.second; i++)
            included_rows.insert(i);
        
        return this->getCofactor(0, included_rows);
    }
}




template <typename ValType>
ValType
FESystem::Numerics::LocalMatrixBase<ValType>::getCofactor(FESystemUInt first_column, const std::set<FESystemUInt>& rows) const
{
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize();
    
    // make sure that the submatrix is square
    FESystemAssert0(s.first-first_column == rows.size(), FESystem::Numerics::MatrixMustBeSquareForOperation);
    
    ValType v=0.0, factor;
    
    std::set<FESystemUInt> included_rows(rows);
    
    // get the cofactor value for the submatrix
    if (s.second-first_column == 1)  // if the last column and one value remain
        return this->getVal(*(rows.begin()), first_column); // the row should have one element ,which is obtained by the row.begin()
    else
    {
        // iterate over all included rows 
        std::set<FESystemUInt>::const_iterator it=rows.begin(), end=rows.end();
        FESystemUInt row_num = 0;
        for ( ; it != end; it++)
        {
            included_rows.erase(*it); // remove the row number for that will create the submatrix
            factor = this->getVal(*it, first_column);
            if (FESystem::Base::comparisonValue<ValType, typename RealOperationType(ValType)>(factor) != 0.0)
            {
                // the power is identified as according to the addition of column number (which is always zero),
                // and the row number "in the submatrix" (and not in the global matrix)
                factor *= pow(-1.0,row_num++); 
                factor *= this->getCofactor(first_column+1,included_rows);
                v += factor;  // the submatrix is defined by the starting column and the retained set of rows
            }
            included_rows.insert(*it); // replace the row that was removed earlier
        }
        
        return v;
    }
}



template <typename ValType>
void
FESystem::Numerics::LocalMatrixBase<ValType>::getInverse(FESystem::Numerics::MatrixBase<ValType>& mat) const
{
    // specialize for 1, 2 and 3-D matrices and if the dimension is greater than that, use the QR method
    std::pair<FESystemUInt, FESystemUInt> s = this->getSize(), s_m=mat.getSize();    
    FESystemAssert0(s.first == s.second, FESystem::Numerics::MatrixMustBeSquareForOperation);
    FESystemAssert4(((s.first==s_m.first) && (s.second==s_m.second)), FESystem::Numerics::MatrixSizeMismatch,
                    s.first, s.second, s_m.first, s_m.second);
    
    if (s.first == 1)
    {
        ValType v=1.0; v/=this->getVal(0,0);
        mat.setVal(0,0, v);
    }
    else if (s.first == 2)
    {
        ValType v,det=(this->getVal(0,0)*this->getVal(1,1)-this->getVal(0,1)*this->getVal(1,0));

        v=this->getVal(1,1); v/=det;
        mat.setVal(0,0, v);

        v=this->getVal(0,0); v/=det;
        mat.setVal(1,1, v);

        v=-this->getVal(0,1); v/=det;
        mat.setVal(0,1, v);

        v=-this->getVal(1,0); v/=det;
        mat.setVal(1,0, v);
    }
    else if (s.first == 3)
    {
        ValType v, v00=this->getVal(0,0),v01=this->getVal(0,1),v02=this->getVal(0,2),
        v10=this->getVal(1,0),v11=this->getVal(1,1),v12=this->getVal(1,2),
        v20=this->getVal(2,0),v21=this->getVal(2,1),v22=this->getVal(2,2),
        det=-v02*v11*v20 + v01*v12*v20 + v02*v10*v21 - v00*v12*v21 - v01*v10*v22 + v00*v11*v22;
        
        v=-v12*v21 + v11*v22; v/=det;
        mat.setVal(0,0,v);

        v=v02*v21 - v01*v22; v/=det;
        mat.setVal(0,1,v);

        v=-v02*v11 + v01*v12; v/=det;
        mat.setVal(0,2,v);

        v=v12*v20 - v10*v22; v/=det;
        mat.setVal(1,0,v);

        v=-v02*v20 + v00*v22; v/=det;
        mat.setVal(1,1,v);

        v=v02*v10 - v00*v12; v/=det;
        mat.setVal(1,2,v);

        v=-v11*v20 + v10*v21; v/=det;
        mat.setVal(2,0,v);

        v=v01*v20 - v00*v21; v/=det;
        mat.setVal(2,1,v);

        v=-v01*v10 + v00*v11; v/=det;
        mat.setVal(2,2,v);

    }
    else
    {        
        // calculate the QR Factorization of this matrix, and used that during the inversion operations
        FESystem::FactorizationSolvers::HouseholderTriangulation<ValType> qr_householder;
        qr_householder.setMatrix(this);
        qr_householder.factorize();  // A = QR
        
        FESystem::Numerics::MatrixBase<ValType>& q_mat = qr_householder.getQMatrix();
        FESystem::Numerics::MatrixBase<ValType>& r_mat = qr_householder.getRMatrix();
        std::auto_ptr<FESystem::Numerics::MatrixBase<ValType> > q_mat_tr(FESystem::Numerics::MatrixCreate<ValType>(FESystem::Numerics::LOCAL_DENSE_MATRIX));
        q_mat_tr->resize(q_mat.getSize().first, q_mat.getSize().second);
        q_mat_tr->copyMatrixTranspose(q_mat);
        
        FESystem::FactorizationSolvers::TriangularBacksubstitution<ValType> back_substitute;
        back_substitute.setTriangularMatrixType(FESystem::FactorizationSolvers::UPPER_TRIANGULAR);
        back_substitute.setMatrix(r_mat);    
        back_substitute.backSubstitute(*q_mat_tr,mat); // A^{-1} = R^{-1} Q^T 
    }
}




template <typename ValType>
const ValType*
FESystem::Numerics::LocalMatrixBase<ValType>::getMatrixValues() const
{
    return &(this->mat_vals[0]);
}



template <typename ValType>
ValType*
FESystem::Numerics::LocalMatrixBase<ValType>::getMatrixValues() 
{
    return &(this->mat_vals[0]);
}


template <>
void
FESystem::Numerics::LocalMatrixBase<FESystemFloat>::write(std::ostream& out) const
{
    unsigned int m = this->getSize().first;
    unsigned int n = this->getSize().second;
    FESystemUInt width = 17, precision = 9;
    
    out << "Size: " << m << ",  " << n << std::endl;
    for (unsigned int i=0; i<m; i++)
    {
        for (unsigned int j=0; j < n; j++)
            if (this->getVal(i,j) == 0.0)
                out << std::setw(width) << this->getVal(i,j);
            else if (log10(fabs(this->getVal(i,j))) >= -1)
                out << std::setw(width) << std::setprecision(precision) << this->getVal(i,j);
            else
                out << std::setw(width) << std::scientific << std::setprecision(precision) << this->getVal(i,j);
        out << std::endl;
    }
}


template <>
void
FESystem::Numerics::LocalMatrixBase<FESystemDouble>::write(std::ostream& out) const
{
    unsigned int m = this->getSize().first;
    unsigned int n = this->getSize().second;
    FESystemUInt width = 17, precision = 9;
    
    out << "Size: " << m << ",  " << n << std::endl;
    for (unsigned int i=0; i<m; i++)
    {
        for (unsigned int j=0; j < n; j++)
            if (this->getVal(i,j) == 0.0)
                out << std::setw(width) << this->getVal(i,j);
            else if (log10(fabs(this->getVal(i,j))) >= -1)
                out << std::setw(width) << std::setprecision(precision) << this->getVal(i,j);
            else
                out << std::setw(width) << std::showpoint << std::setprecision(precision) << this->getVal(i,j);
        out << std::endl;
    }
}



template <>
void
FESystem::Numerics::LocalMatrixBase<FESystemComplexFloat>::write(std::ostream& out) const
{
    unsigned int m = this->getSize().first;
    unsigned int n = this->getSize().second;
    FESystemUInt width = 25, precision = 6;
    
    out << "Size: " << m << ",  " << n << std::endl;
    for (unsigned int i=0; i<m; i++)
    {
        for (unsigned int j=0; j < n; j++)
            out << std::setw(width) << std::setprecision(precision) << this->getVal(i,j);
        out << std::endl;
    }
}


template <>
void
FESystem::Numerics::LocalMatrixBase<FESystemComplexDouble>::write(std::ostream& out) const
{
    unsigned int m = this->getSize().first;
    unsigned int n = this->getSize().second;
    FESystemUInt width = 25, precision = 6;
    
    out << "Size: " << m << ",  " << n << std::endl;
    for (unsigned int i=0; i<m; i++)
    {
        for (unsigned int j=0; j < n; j++)
            out << std::setw(width) << std::setprecision(precision) << this->getVal(i,j);
        out << std::endl;
    }
}




template <typename ValType>
void
FESystem::Numerics::LocalMatrixBase<ValType>::writeDetailed(std::ostream& out) const
{
    unsigned int m = this->getSize().first;
    unsigned int n = this->getSize().second;
    
    out << "Dense Matrix" << std::endl;
    out << "Size: " << m << ",  " << n << std::endl;
    for (unsigned int i=0; i<m; i++)
    {
        out << "Row #" << i << "  :  ";
        
        for (unsigned int j=0; j < n; j++)
        {
            //            out.width(15);
            //            out.precision(10);
            out << "   :("<<j<<")->    " << this->getVal(i,j) ;
        }
        out << std::endl;
    }
}

/***************************************************************************************/
// Template instantiations for some generic classes

INSTANTIATE_CLASS_FOR_ALL_DATA_TYPES(FESystem::Numerics::LocalMatrixBase);


/***************************************************************************************/


