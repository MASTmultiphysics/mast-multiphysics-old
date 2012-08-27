//
//  VonKarmanStrain2D.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

// FESystem includes
#include "Disciplines/Structure/VonKarmanStrain2D.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"
#include "Disciplines/Structure/Membrane.h"



FESystem::Structures::VonKarmanStrain2D::VonKarmanStrain2D():
FESystem::Structures::Structural2DElementBase(),
membrane_elem(NULL)
{
    
}


FESystem::Structures::VonKarmanStrain2D::~VonKarmanStrain2D()
{
    
}



void
FESystem::Structures::VonKarmanStrain2D::clear()
{
    FESystem::Structures::Structural2DElementBase::clear();
    this->membrane_elem = NULL;
}



void
FESystem::Structures::VonKarmanStrain2D::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                    FESystem::Structures::Membrane& membrane)
{
    FESystem::Structures::StructuralElementBase::initialize(elem, fe, q_rule, 0, 0, 0);
    this->membrane_elem = & membrane;
}



void
FESystem::Structures::VonKarmanStrain2D::calculateTangentStiffnessMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& membrane_sol, FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == n) && (s.second== n)), FESystem::Numerics::MatrixSizeMismatch, n, n, s.first, s.second);
    
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, stress_mat, tmp_mat1, tmp_mat2;
    static FESystem::Numerics::LocalVector<FESystemDouble> pt;
    stress_mat.resize(2,2); B_mat.resize(2, n); tmp_mat1.resize(2, n), tmp_mat2.resize(n, n);
    stress_mat.zero(); B_mat.zero(); tmp_mat1.zero(); tmp_mat2.zero();
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    mat.zero();
    this->membrane_elem->getStressTensor(pt, membrane_sol, stress_mat);
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        this->calculateOperatorMatrix(*(q_pts[i]), B_mat);
        
        stress_mat.matrixRightMultiply(1.0, B_mat, tmp_mat1); // sigma B
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2); // B^T sigma B
        
        mat.add(q_weight[i]*jac, tmp_mat2);
    }
}



void
FESystem::Structures::VonKarmanStrain2D::calculateOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    FESystemAssert4(((s.first == 2) && (s.second== n)), FESystem::Numerics::MatrixSizeMismatch, 2, n, s.first, s.second);
    
    static std::vector<FESystemUInt> derivatives_x(2),derivatives_y(2);
    derivatives_x[0] = 1; derivatives_x[1] = 0;
    derivatives_y[0] = 0; derivatives_y[1] = 1;
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n);
    B_mat.zero();
    
    Nvec.zero();
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_x, pt, Nvec);
    B_mat.setRowVals(0, 0,  n-1, Nvec); // dw/dx
    
    Nvec.zero();
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_y, pt, Nvec);
    B_mat.setRowVals(1, 0,  n-1, Nvec); // dw/dy
}


