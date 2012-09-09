//
//  LinearPlateElementBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/23/12.
//
//


// FESystem includes
#include "Disciplines/Structure/LinearPlateElementBase.h"
#include "Mesh/ElemBase.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "Base/FESystemExceptions.h"
#include "Geom/Point.h"



FESystem::Structures::LinearPlateElementBase::LinearPlateElementBase():
FESystem::Structures::Structural2DElementBase()
{
}


FESystem::Structures::LinearPlateElementBase::~LinearPlateElementBase()
{
    
}



FESystemUInt
FESystem::Structures::LinearPlateElementBase::getNElemDofs() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    return 3*this->geometric_elem->getNNodes();
}


void
FESystem::Structures::LinearPlateElementBase::getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec)
{
    FESystemUInt n = this->geometric_elem->getNNodes();
    vec.resize(3*n);
    
    for (FESystemUInt i=0; i<n; i++) vec[i] = 2*n + i; // w-displacement
    for (FESystemUInt i=0; i<n; i++) vec[i+n] = 3*n + i; // theta-x
    for (FESystemUInt i=0; i<n; i++) vec[i+2*n] = 4*n + i; // theta-y
}

void
FESystem::Structures::LinearPlateElementBase::calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == 3*n) && (s.second == 3*n)), FESystem::Numerics::MatrixSizeMismatch, 3*n, 3*n, s.first, s.second);
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, C_mat, tmp_mat1, tmp_mat2;
    C_mat.resize(3,3); B_mat.resize(3, 3*n); tmp_mat1.resize(3, 3*n), tmp_mat2.resize(3*n, 3*n);
    C_mat.zero(); B_mat.zero(); tmp_mat1.zero(); tmp_mat2.zero();
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    mat.zero();
    this->getMaterialMassMatrix(C_mat);
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        this->calculateInertiaOperatorMatrix(*(q_pts[i]), B_mat);
        
        C_mat.matrixRightMultiply(1.0, B_mat, tmp_mat1);
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
        
        mat.add(q_weight[i]*jac, tmp_mat2);
    }
    
}



void
FESystem::Structures::LinearPlateElementBase::calculateDiagonalMassMatrix(FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    const FESystemUInt n = this->geometric_elem->getNNodes();
    
    FESystemAssert2(vec.getSize() == 3*n, FESystem::Exception::DimensionsDoNotMatch, 3*n, vec.getSize());
    
    vec.zero();
    FESystemDouble wt = this->geometric_elem->getElementSize(*(this->finite_element), *(this->quadrature)) * this->th_val * this->rho_val;
    
    wt /= (1.0*n);
    vec.setAllVals(wt*10e-4);
    for (FESystemUInt i=0; i<n; i++)
        vec.setVal(i, wt);
}



void
FESystem::Structures::LinearPlateElementBase::calculateInertiaOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    FESystemAssert4(((s.first == 3) && (s.second== 3*n)), FESystem::Numerics::MatrixSizeMismatch, 3, 3*n, s.first, s.second);
    
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n);
    B_mat.zero();
    
    Nvec.zero();
    this->finite_element->getShapeFunction(pt, Nvec);
    B_mat.setRowVals(0,   0,   n-1, Nvec); // w
    B_mat.setRowVals(1,   n, 2*n-1, Nvec); // theta_x
    B_mat.setRowVals(2, 2*n, 3*n-1, Nvec); // theta_y
}



void
FESystem::Structures::LinearPlateElementBase::getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == 3) && (s.second== 3)), FESystem::Numerics::MatrixSizeMismatch, 3, 3, s.first, s.second);
    
    mat.setVal(0, 0, this->rho_val * this->th_val);
    mat.setVal(1, 1, 1.0e-12 * this->rho_val * this->th_val);
    mat.setVal(2, 2, 1.0e-12 * this->rho_val * this->th_val);
}



