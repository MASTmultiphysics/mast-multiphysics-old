//
//  TorsionBar.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//


// FESystem includes
#include "Disciplines/Structure/TorsionBar.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"



FESystem::Structures::TorsionBar::TorsionBar():
FESystem::Structures::Structural1DElementBase(),
polar_inertia_val(0.0),
J_val(0.0)
{
    
}


FESystem::Structures::TorsionBar::~TorsionBar()
{
    
}



void
FESystem::Structures::TorsionBar::clear()
{
    this->polar_inertia_val = 0.0;
    this->J_val = 0.0;
    FESystem::Structures::StructuralElementBase::clear();
}



FESystemUInt
FESystem::Structures::TorsionBar::getNElemDofs() const
{
    return this->geometric_elem->getNNodes();
}


void
FESystem::Structures::TorsionBar::getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec)
{
    FESystemUInt n = this->geometric_elem->getNNodes();
    vec.resize(n);
    
    for (FESystemUInt i=0; i<n; i++) vec[i] = 3*n + i; // theta-x
}


void
FESystem::Structures::TorsionBar::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                             FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble polar_inertia, FESystemDouble J)
{
    FESystem::Structures::StructuralElementBase::initialize(elem, fe, q_rule, E, nu, rho);
    
    this->polar_inertia_val = polar_inertia;
    this->J_val = J;
}



void
FESystem::Structures::TorsionBar::calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();

    FESystemAssert4(((s.first == n) && (s.second== n)), FESystem::Numerics::MatrixSizeMismatch, n, n, s.first, s.second);

    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, C_mat, tmp_mat1, tmp_mat2;
    C_mat.resize(1,1); B_mat.resize(1, n); tmp_mat1.resize(1, n), tmp_mat2.resize(n, n);
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
FESystem::Structures::TorsionBar::calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == n) && (s.second== n)), FESystem::Numerics::MatrixSizeMismatch, n, n, s.first, s.second);
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, C_mat, tmp_mat1, tmp_mat2;
    C_mat.resize(1,1); B_mat.resize(1, n); tmp_mat1.resize(1, n), tmp_mat2.resize(n, n);
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
FESystem::Structures::TorsionBar::getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    FESystemAssert4(((s.first == 1) && (s.second== 1)), FESystem::Numerics::MatrixSizeMismatch, 1, 1, s.first, s.second);
    
    mat.setVal(0, 0, this->rho_val*this->polar_inertia_val);
}



void
FESystem::Structures::TorsionBar::getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    FESystemAssert4(((s.first == 1) && (s.second== 1)), FESystem::Numerics::MatrixSizeMismatch, 1, 1, s.first, s.second);
    
    mat.setVal(0, 0, this->G_val*this->J_val);
}



void
FESystem::Structures::TorsionBar::calculateOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat, FESystemBoolean if_strain)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    FESystemAssert4(((s.first == 1) && (s.second== n)), FESystem::Numerics::MatrixSizeMismatch, 1, n, s.first, s.second);
    
    static std::vector<FESystemUInt> derivatives(1);
    derivatives[0] = 1;
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n); Nvec.zero();
    B_mat.zero();
    
    if (if_strain)
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives, pt, Nvec);
    else
        this->finite_element->getShapeFunction(pt, Nvec);
    
    B_mat.setRowVals(0, 0, n-1, Nvec);
}


