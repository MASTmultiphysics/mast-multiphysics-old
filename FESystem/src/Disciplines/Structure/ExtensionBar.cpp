//
//  ExtensionBar.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

// FESystem includes
#include "Disciplines/Structure/ExtensionBar.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"


FESystem::Structures::ExtensionBar::ExtensionBar():
FESystem::Structures::Structural1DElementBase(),
area_val(0.0)
{
}


FESystem::Structures::ExtensionBar::~ExtensionBar()
{
    
}


void
FESystem::Structures::ExtensionBar::clear()
{
    this->area_val = 0.0;
    FESystem::Structures::StructuralElementBase::clear();
}



void
FESystem::Structures::ExtensionBar::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                               FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble area)
{
    FESystem::Structures::StructuralElementBase::initialize(elem, fe, q_rule, E, nu, rho);
    
    this->area_val = area;
}



void
FESystem::Structures::ExtensionBar::calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->finite_element->getNShapeFunctions();
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
FESystem::Structures::ExtensionBar::calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->finite_element->getNShapeFunctions();
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
FESystem::Structures::ExtensionBar::getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    FESystemAssert0(((s.first == s.second) && (s.first == 1)), FESystem::Exception::InvalidValue);
    
    mat.setVal(0, 0, this->rho_val*this->area_val);
}



void
FESystem::Structures::ExtensionBar::getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    FESystemAssert4(((s.first == 1) && (s.second== 1)), FESystem::Numerics::MatrixSizeMismatch, 1, 1, s.first, s.second);
    
    mat.setVal(0, 0, this->E_val*this->area_val);
}


void
FESystem::Structures::ExtensionBar::calculateOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat, FESystemBoolean if_strain)
{
    const FESystemUInt n = this->finite_element->getNShapeFunctions();
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


