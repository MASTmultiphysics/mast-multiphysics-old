//
//  Membrane.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

// FESystem includes
#include "Disciplines/Structure/Membrane.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"


FESystem::Structures::Membrane::Membrane():
FESystem::Structures::Structural2DElementBase()
{
}


FESystem::Structures::Membrane::~Membrane()
{
    
}


FESystemUInt
FESystem::Structures::Membrane::getNElemDofs() const
{
    return 2*this->geometric_elem->getNNodes();
}


void
FESystem::Structures::Membrane::getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec)
{
    FESystemUInt n = this->geometric_elem->getNNodes();
    vec.resize(2*n);
    
    for (FESystemUInt i=0; i<2*n; i++) vec[i] = i; // u- and v-displacement
}



void
FESystem::Structures::Membrane::calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->finite_element->getNShapeFunctions();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == 2*n) && (s.second== 2*n)), FESystem::Numerics::MatrixSizeMismatch, 2*n, 2*n, s.first, s.second);
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, C_mat, tmp_mat1, tmp_mat2;
    C_mat.resize(2,2); B_mat.resize(2, 2*n); tmp_mat1.resize(2, 2*n), tmp_mat2.resize(2*n, 2*n);
    C_mat.zero(); B_mat.zero(); tmp_mat1.zero(); tmp_mat2.zero();
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    mat.zero();
    this->getMaterialMassMatrix(C_mat);
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        this->calculateOperatorMatrix(*(q_pts[i]), B_mat, false);
        
        C_mat.matrixRightMultiply(1.0, B_mat, tmp_mat1);
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
        
        mat.add(q_weight[i]*jac, tmp_mat2);
    }
}



void
FESystem::Structures::Membrane::calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->finite_element->getNShapeFunctions();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == 2*n) && (s.second== 2*n)), FESystem::Numerics::MatrixSizeMismatch, 2*n, 2*n, s.first, s.second);
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, C_mat, tmp_mat1, tmp_mat2;
    C_mat.resize(3,3); B_mat.resize(3, 2*n); tmp_mat1.resize(3, 2*n), tmp_mat2.resize(2*n, 2*n);
    C_mat.zero(); B_mat.zero(); tmp_mat1.zero(); tmp_mat2.zero();
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    mat.zero();
    this->getMaterialComplianceMatrix(C_mat);

    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        this->calculateOperatorMatrix(*(q_pts[i]), B_mat, true);
        
        C_mat.matrixRightMultiply(1.0, B_mat, tmp_mat1);
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
        
        mat.add(q_weight[i]*jac, tmp_mat2);
    }
}



void
FESystem::Structures::Membrane::calculateOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat, FESystemBoolean if_strain)
{
    const FESystemUInt n = this->finite_element->getNShapeFunctions();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    FESystemAssert4(((s.first == 3) && (s.second== 2*n)), FESystem::Numerics::MatrixSizeMismatch, 3, 2*n, s.first, s.second);
    
    static std::vector<FESystemUInt> derivatives_x(2),derivatives_y(2);
    derivatives_x[0] = 1; derivatives_x[1] = 0;
    derivatives_y[0] = 0; derivatives_y[1] = 1;
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n); Nvec.zero();
    B_mat.zero();
    
    if (if_strain)
    {
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_x, pt, Nvec);
        B_mat.setRowVals(0, 0, n-1, Nvec); // epsilon x
        B_mat.setRowVals(2, n, 2*n-1, Nvec); // gamma xy
        
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_y, pt, Nvec);
        B_mat.setRowVals(1, n, 2*n-1, Nvec); // epsilon y
        B_mat.setRowVals(2, 0, n-1, Nvec); // gamma xy
    }
    else
    {
        this->finite_element->getShapeFunction(pt, Nvec);
        B_mat.setRowVals(0, 0, n-1, Nvec); // u
        B_mat.setRowVals(1, n, 2*n-1, Nvec); // v
    }
}




void
FESystem::Structures::Membrane::getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == 2) && (s.second== 2)), FESystem::Numerics::MatrixSizeMismatch, 2, 2, s.first, s.second);
    
    mat.setVal(0, 0, this->rho_val * this->th_val);
    mat.setVal(1, 1, this->rho_val * this->th_val);
}



void
FESystem::Structures::Membrane::getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->finite_element->getNShapeFunctions();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == 3) && (s.second== 2*n)), FESystem::Numerics::MatrixSizeMismatch, 3, 2*n, s.first, s.second);

    FESystemDouble val = this->E_val/(1.0-this->nu_val*this->nu_val);
    
    mat.setVal(0, 0, val);
    mat.setVal(0, 1, this->nu_val*val);
    mat.setVal(1, 0, this->nu_val*val);
    mat.setVal(1, 1, val);
    mat.setVal(2, 2, this->G_val);
    mat.scale(this->th_val);
}

