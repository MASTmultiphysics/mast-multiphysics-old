//
//  TimoshenkoBeam.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/26/12.
//
//

// FESystem includes
#include "Disciplines/Structure/TimoshenkoBeam.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"


FESystem::Structures::TimoshenkoBeam::TimoshenkoBeam():
FESystem::Structures::LinearBeamElementBase(),
quadrature_bending(NULL),
quadrature_shear(NULL),
kappa(5.0/6.0)
{
}


FESystem::Structures::TimoshenkoBeam::~TimoshenkoBeam()
{
    
}



void
FESystem::Structures::TimoshenkoBeam::clear()
{
    FESystem::Structures::LinearBeamElementBase::clear();
    this->quadrature_bending = NULL;
    this->quadrature_shear = NULL;
}





void
FESystem::Structures::TimoshenkoBeam::getStressTensor(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                                            FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




void
FESystem::Structures::TimoshenkoBeam::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_bend,
                                                 const FESystem::Quadrature::QuadratureBase& q_shear, FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble I_tr, FESystemDouble I_ch,
                                                 FESystemDouble A, FESystemBoolean if_lateral)
{
    FESystem::Structures::LinearBeamElementBase::initialize(elem, fe, q_bend, E, nu, rho, I_tr, I_ch, A, if_lateral);
    this->quadrature_bending = &q_bend;
    this->quadrature_shear = &q_shear;
}



void
FESystem::Structures::TimoshenkoBeam::calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, C_mat, tmp_mat1, tmp_mat2;

    if (this->if_include_lateral_stiffness)
    {
        FESystemAssert4(((s.first == 4*n) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 4*n, 4*n, s.first, s.second);
        C_mat.resize(4,4); B_mat.resize(4, 4*n); tmp_mat1.resize(4, 4*n), tmp_mat2.resize(4*n, 4*n);
    }
    else
    {
        FESystemAssert4(((s.first == 2*n) && (s.second == 2*n)), FESystem::Numerics::MatrixSizeMismatch, 2*n, 2*n, s.first, s.second);
        C_mat.resize(2,2); B_mat.resize(2, 2*n); tmp_mat1.resize(2, 2*n), tmp_mat2.resize(2*n, 2*n);
    }
    
    
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
FESystem::Structures::TimoshenkoBeam::calculateDiagonalMassMatrix(FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const FESystemUInt s = vec.getSize();
    
    if (this->if_include_lateral_stiffness)
    {
        FESystemAssert2(s == 4*n, FESystem::Exception::DimensionsDoNotMatch, s, 4*n);
    }
    else
    {
        FESystemAssert2(s == 2*n, FESystem::Exception::DimensionsDoNotMatch, s, 2*n);
    }
    
    vec.zero();
    FESystemDouble wt = this->geometric_elem->getElementSize(*(this->finite_element), *(this->quadrature_bending)) * this->area_val * this->rho_val;
    
    wt /= (1.0*n);
    vec.setAllVals(wt*10e-4);

    if (this->if_include_lateral_stiffness)
        for (FESystemUInt i=0; i<2*n; i++)
            vec.setVal(i, wt);
    else
        for (FESystemUInt i=0; i<n; i++)
            vec.setVal(i, wt);

}



void
FESystem::Structures::TimoshenkoBeam::calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, B_mat_shear, C_mat_bend, C_mat_shear, tmp_mat1, tmp_mat2;
    
    if (this->if_include_lateral_stiffness)
    {
        FESystemAssert4(((s.first == 4*n) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 4*n, 4*n, s.first, s.second);
        C_mat_bend.resize(1,1); C_mat_shear.resize(2,2); B_mat.resize(1, 4*n); B_mat_shear.resize(2, 4*n); tmp_mat1.resize(1, 4*n), tmp_mat2.resize(4*n, 4*n);
    }
    else
    {
        FESystemAssert4(((s.first == 2*n) && (s.second == 2*n)), FESystem::Numerics::MatrixSizeMismatch, 2*n, 2*n, s.first, s.second);
        C_mat_bend.resize(1,1); C_mat_shear.resize(1,1); B_mat.resize(1, 2*n); B_mat_shear.resize(1, 2*n); tmp_mat1.resize(1, 2*n), tmp_mat2.resize(2*n, 2*n);
    }
    
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
    
    // transverse shear contribution
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
FESystem::Structures::TimoshenkoBeam::calculateInertiaOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n);
    Nvec.zero(); B_mat.zero();
    this->finite_element->getShapeFunction(pt, Nvec);

    if (this->if_include_lateral_stiffness)
    {
        FESystemAssert4(((s.first == 4) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 1, 4*n, s.first, s.second);

        B_mat.setRowVals(0,   0,   n-1, Nvec); // v
        B_mat.setRowVals(1,   n, 2*n-1, Nvec); // w
        B_mat.setRowVals(2, 2*n, 3*n-1, Nvec); // theta_y
        B_mat.setRowVals(3, 3*n, 4*n-1, Nvec); // theta_z
    }
    else
    {
        FESystemAssert4(((s.first == 4) && (s.second == 2*n)), FESystem::Numerics::MatrixSizeMismatch, 1, 2*n, s.first, s.second);
        B_mat.resize(1, 2*n);
        
        B_mat.setRowVals(0,   0,   n-1, Nvec); // w
        B_mat.setRowVals(1,   n, 2*n-1, Nvec); // theta_y
    }
}



void
FESystem::Structures::TimoshenkoBeam::calculateBendingOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    static std::vector<FESystemUInt> derivatives_x(1);
    derivatives_x[0] = 1;
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n);
    B_mat.zero();
    Nvec.zero();
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_x, pt, Nvec);
    
    if (this->if_include_lateral_stiffness)
    {
        FESystemAssert4(((s.first == 1) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 1, 4*n, s.first, s.second);

        B_mat.setRowVals(0, 2*n,  3*n-1, Nvec); // epsilon-x: thetay
        Nvec.scale(-1.0 * this->I_ch_val/this->I_tr_val);
        B_mat.setRowVals(0, 3*n,  4*n-1, Nvec); // epsilon-x: thetaz
    }
    else
    {
        FESystemAssert4(((s.first == 1) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 1, 4*n, s.first, s.second);
        
        Nvec.scale(-1.0);
        B_mat.setRowVals(0, n,  2*n-1, Nvec); // epsilon-x: thetay
    }
}



void
FESystem::Structures::TimoshenkoBeam::calculateShearOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    static std::vector<FESystemUInt> derivatives_x(2);
    derivatives_x[0] = 1;
    
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n);
    B_mat.zero();

    if (this->if_include_lateral_stiffness)
    {
        FESystemAssert4(((s.first == 2) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 2, 4*n, s.first, s.second);
        
        Nvec.zero();
        this->finite_element->getShapeFunction(pt, Nvec);
        B_mat.setRowVals(0, 2*n,  3*n-1, Nvec); // gamma-xz:  thetay
        Nvec.scale(-1.0);
        B_mat.setRowVals(1, 3*n,  4*n-1, Nvec); // gamma-xy : thetaz
        
        Nvec.zero();
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_x, pt, Nvec);
        B_mat.setRowVals(0,   n,  2*n-1, Nvec); // gamma-xz:  w
        B_mat.setRowVals(1,   0,    n-1, Nvec); // gamma-xy:  v
    }
    else
    {
        FESystemAssert4(((s.first == 1) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 1, 4*n, s.first, s.second);
        
        Nvec.zero();
        this->finite_element->getShapeFunction(pt, Nvec);
        B_mat.setRowVals(0,   n,  2*n-1, Nvec); // gamma-xz:  thetay
        
        Nvec.zero();
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives_x, pt, Nvec);
        B_mat.setRowVals(0,   0,    n-1, Nvec); // gamma-xz:  w
    }
}



void
FESystem::Structures::TimoshenkoBeam::getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    if (this->if_include_lateral_stiffness)
    {
        FESystemAssert4(((s.first == 4) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 4, 4*n, s.first, s.second);
        mat.setVal(0, 0, this->rho_val * this->area_val);
        mat.setVal(1, 1, this->rho_val * this->area_val);
        mat.setVal(3, 3, 1.0e-12 * this->rho_val * this->area_val);
        mat.setVal(4, 4, 1.0e-12 * this->rho_val * this->area_val);
    }
    else
    {
        FESystemAssert4(((s.first == 2) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 2, 4*n, s.first, s.second);
        mat.setVal(0, 0, this->rho_val * this->area_val);
        mat.setVal(2, 2, 1.0e-12 * this->rho_val * this->area_val);
    }
    
}



void
FESystem::Structures::TimoshenkoBeam::getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& bend_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& shear_mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s_b = bend_mat.getSize();
    const std::pair<FESystemUInt, FESystemUInt> s_s = shear_mat.getSize();
    
    if (this->if_include_lateral_stiffness)
    {
        FESystemAssert4(((s_b.first == 1) && (s_b.second== 1)), FESystem::Numerics::MatrixSizeMismatch, 1, 1, s_b.first, s_b.second);
        FESystemAssert4(((s_s.first == 2) && (s_s.second== 2)), FESystem::Numerics::MatrixSizeMismatch, 2, 2, s_s.first, s_s.second);

        shear_mat.setVal(0, 0, this->G_val);
        shear_mat.setVal(1, 1, this->G_val);
    }
    else
    {
        FESystemAssert4(((s_b.first == 1) && (s_b.second== 1)), FESystem::Numerics::MatrixSizeMismatch, 1, 1, s_b.first, s_b.second);
        FESystemAssert4(((s_s.first == 1) && (s_s.second== 1)), FESystem::Numerics::MatrixSizeMismatch, 1, 1, s_s.first, s_s.second);
        
        shear_mat.setVal(0, 0, this->G_val);
    }

    bend_mat.setVal(0, 0, this->E_val * this->I_tr_val);
        
    shear_mat.scale(this->kappa*this->area_val);
}


