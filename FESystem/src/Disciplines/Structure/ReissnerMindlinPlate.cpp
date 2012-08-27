//
//  ReissnerMindlinPlate.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

// FESystem includes
#include "Disciplines/Structure/ReissnerMindlinPlate.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"


FESystem::Structures::ReissnerMindlinPlate::ReissnerMindlinPlate():
FESystem::Structures::LinearPlateElementBase(),
quadrature_bending(NULL),
quadrature_shear(NULL),
kappa(5.0/6.0)
{
}


FESystem::Structures::ReissnerMindlinPlate::~ReissnerMindlinPlate()
{
    
}



void
FESystem::Structures::ReissnerMindlinPlate::clear()
{
    FESystem::Structures::Structural2DElementBase::clear();
    this->quadrature_bending = NULL;
    this->quadrature_shear = NULL;
}





void
FESystem::Structures::ReissnerMindlinPlate::getStressTensor(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                                            FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




void
FESystem::Structures::ReissnerMindlinPlate::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_bend,
                                                       const FESystem::Quadrature::QuadratureBase& q_shear, FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble th)
{
    FESystem::Structures::Structural2DElementBase::initialize(elem, fe, q_bend, E, nu, rho, th);
    this->quadrature_bending = &q_bend;
    this->quadrature_shear = &q_shear;
}



void
FESystem::Structures::ReissnerMindlinPlate::calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
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
FESystem::Structures::ReissnerMindlinPlate::calculateDiagonalMassMatrix(FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
    
    FESystemAssert2(vec.getSize() == 3*n, FESystem::Exception::DimensionsDoNotMatch, 3*n, vec.getSize());
    
    vec.zero();
    FESystemDouble wt = this->geometric_elem->getElementSize(*(this->finite_element), *(this->quadrature_bending)) * this->th_val * this->rho_val;
    
    wt /= (1.0*n);
    vec.setAllVals(wt*10e-4);
    for (FESystemUInt i=0; i<n; i++)
        vec.setVal(i, wt);
}



void
FESystem::Structures::ReissnerMindlinPlate::calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == 3*n) && (s.second== 3*n)), FESystem::Numerics::MatrixSizeMismatch, 3*n, 3*n, s.first, s.second);
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, C_mat_bend, C_mat_shear, tmp_mat1, tmp_mat2;
    C_mat_bend.resize(3,3); C_mat_shear.resize(3,3); B_mat.resize(3, 3*n); tmp_mat1.resize(3, 3*n), tmp_mat2.resize(3*n, 3*n);
    C_mat_bend.zero(); C_mat_shear.zero(); B_mat.zero(); tmp_mat1.zero(); tmp_mat2.zero();
    
    const std::vector<FESystem::Geometry::Point*>& q_pts_bend = this->quadrature_bending->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight_bend = this->quadrature_bending->getQuadraturePointWeights();

    const std::vector<FESystem::Geometry::Point*>& q_pts_shear = this->quadrature_shear->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight_shear = this->quadrature_shear->getQuadraturePointWeights();

    FESystemDouble jac=0.0;
    mat.zero();
    this->getMaterialComplianceMatrix(C_mat_bend, C_mat_shear);
        
    // bending contribution
    for (FESystemUInt i=0; i<q_pts_bend.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts_bend[i]));
        this->calculateBendingOperatorMatrix(*(q_pts_bend[i]), B_mat);
        
        C_mat_bend.matrixRightMultiply(1.0, B_mat, tmp_mat1);
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
        
        mat.add(q_weight_bend[i]*jac, tmp_mat2);
    }

    // shear contribution
    B_mat.zero();
    for (FESystemUInt i=0; i<q_pts_shear.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts_shear[i]));
        this->calculateShearOperatorMatrix(*(q_pts_shear[i]), B_mat);
        
        C_mat_shear.matrixRightMultiply(1.0, B_mat, tmp_mat1);
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
        
        mat.add(q_weight_shear[i]*jac, tmp_mat2);
    }
}




void
FESystem::Structures::ReissnerMindlinPlate::calculateInertiaOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
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
FESystem::Structures::ReissnerMindlinPlate::calculateBendingOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    FESystemAssert4(((s.first == 3) && (s.second== 3*n)), FESystem::Numerics::MatrixSizeMismatch, 3, 3*n, s.first, s.second);
    
    static std::vector<FESystemUInt> derivatives_x(2),derivatives_y(2);
    derivatives_x[0] = 1; derivatives_x[1] = 0;
    derivatives_y[0] = 0; derivatives_y[1] = 1;
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n); 
    B_mat.zero();

    Nvec.zero();
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_x, pt, Nvec);
    B_mat.setRowVals(0, 2*n,  3*n-1, Nvec); // epsilon-x: thetay
    Nvec.scale(-1.0);
    B_mat.setRowVals(2,   n,  2*n-1, Nvec); // gamma-xy : thetax
    
    Nvec.zero();
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_y, pt, Nvec);
    B_mat.setRowVals(2, 2*n,  3*n-1, Nvec); // gamma-xy : thetay
    Nvec.scale(-1.0);
    B_mat.setRowVals(1,   n,  2*n-1, Nvec); // epsilon-y: thetax
}



void
FESystem::Structures::ReissnerMindlinPlate::calculateShearOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();;
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    FESystemAssert4(((s.first == 3) && (s.second== 3*n)), FESystem::Numerics::MatrixSizeMismatch, 3, 3*n, s.first, s.second);
    
    static std::vector<FESystemUInt> derivatives_x(2),derivatives_y(2);
    derivatives_x[0] = 1; derivatives_x[1] = 0;
    derivatives_y[0] = 0; derivatives_y[1] = 1;
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n);
    B_mat.zero();
    
    Nvec.zero();
    this->finite_element->getShapeFunction(pt, Nvec);
    B_mat.setRowVals(0, 2*n,  3*n-1, Nvec); // gamma-xz:  thetay
    Nvec.scale(-1.0);
    B_mat.setRowVals(1,   n,  2*n-1, Nvec); // gamma-yz : thetax

    Nvec.zero();
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_x, pt, Nvec);
    B_mat.setRowVals(0, 0,  n-1, Nvec); // gamma-xz:  w

    Nvec.zero();
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_y, pt, Nvec);
    B_mat.setRowVals(1, 0,  n-1, Nvec); // gamma-yz : w
}



void
FESystem::Structures::ReissnerMindlinPlate::getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == 3) && (s.second== 3)), FESystem::Numerics::MatrixSizeMismatch, 3, 3, s.first, s.second);
    
    mat.setVal(0, 0, this->rho_val * this->th_val);
    mat.setVal(1, 1, 1.0e-12 * this->rho_val * this->th_val);
    mat.setVal(2, 2, 1.0e-12 * this->rho_val * this->th_val);
}



void
FESystem::Structures::ReissnerMindlinPlate::getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& bend_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& shear_mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s_b = bend_mat.getSize();
    const std::pair<FESystemUInt, FESystemUInt> s_s = shear_mat.getSize();
    
    FESystemAssert4(((s_b.first == 3) && (s_b.second== 3)), FESystem::Numerics::MatrixSizeMismatch, 3, 3, s_b.first, s_b.second);
    FESystemAssert4(((s_s.first == 3) && (s_s.second== 3)), FESystem::Numerics::MatrixSizeMismatch, 3, 3, s_s.first, s_s.second);
    
    FESystemDouble val = this->E_val/(1.0-this->nu_val*this->nu_val);
    
    bend_mat.setVal(0, 0, val);
    bend_mat.setVal(0, 1, this->nu_val*val);
    bend_mat.setVal(1, 0, this->nu_val*val);
    bend_mat.setVal(1, 1, val);
    bend_mat.setVal(2, 2, this->G_val);
    bend_mat.scale(pow(this->th_val,3)/12.0);
    
    shear_mat.setVal(0, 0, this->G_val);
    shear_mat.setVal(1, 1, this->G_val);
    shear_mat.scale(this->kappa*this->th_val);

}



