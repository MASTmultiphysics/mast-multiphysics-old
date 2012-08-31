//
//  LinearBeamElementBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/23/12.
//
//


// FESystem includes
#include "Disciplines/Structure/LinearBeamElementBase.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"


FESystem::Structures::LinearBeamElementBase::LinearBeamElementBase():
FESystem::Structures::Structural1DElementBase(),
I_tr_val(0.0),
I_ch_val(0.0)
{
    
}
            


FESystem::Structures::LinearBeamElementBase::~LinearBeamElementBase()
{
    
}


FESystemUInt
FESystem::Structures::LinearBeamElementBase::getNElemDofs() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    return 4*this->geometric_elem->getNNodes();
}





void
FESystem::Structures::LinearBeamElementBase::clear()
{
    FESystem::Structures::Structural1DElementBase::clear();
    this->I_tr_val = 0.0;
    this->I_ch_val = 0.0;
}


void
FESystem::Structures::LinearBeamElementBase::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                        FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble I_tr, FESystemDouble I_ch,  FESystemDouble A)
{
    FESystem::Structures::Structural1DElementBase::initialize(elem,fe,q_rule,E,nu,rho,A);
    this->I_tr_val = I_tr;
    this->I_ch_val = I_ch;
}
            



void
FESystem::Structures::LinearBeamElementBase::getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec)
{
    FESystemUInt n = this->geometric_elem->getNNodes();
    
    vec.resize(4*n);
    for (FESystemUInt i=0; i<2*n; i++) vec[i] = n + i; // v-, w-displacement
    for (FESystemUInt i=0; i<2*n; i++) vec[2*n+i] = 4*n + i; // theta-y and theta-z
}



void
FESystem::Structures::LinearBeamElementBase::calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> B_mat, C_mat, tmp_mat1, tmp_mat2;
    
    FESystemAssert4(((s.first == 4*n) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 4*n, 4*n, s.first, s.second);
    C_mat.resize(4,4); B_mat.resize(4, 4*n); tmp_mat1.resize(4, 4*n), tmp_mat2.resize(4*n, 4*n);
    
    
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
FESystem::Structures::LinearBeamElementBase::calculateDiagonalMassMatrix(FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const FESystemUInt s = vec.getSize();
    
    FESystemAssert2(s == 4*n, FESystem::Exception::DimensionsDoNotMatch, s, 4*n);
    
    vec.zero();
    FESystemDouble wt = this->geometric_elem->getElementSize(*(this->finite_element), *(this->quadrature)) * this->area_val * this->rho_val;
    
    wt /= (1.0*n);
    vec.setAllVals(wt*10e-4);
    
    for (FESystemUInt i=0; i<2*n; i++)
        vec.setVal(i, wt);
}




void
FESystem::Structures::LinearBeamElementBase::calculateDistributedLoad(FESystemDouble v_val, FESystemDouble w_val, FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const FESystemUInt s = vec.getSize();
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    static std::vector<FESystemUInt> v_dof_indices, w_dof_indices;
    Nvec.resize(n); v_dof_indices.resize(n); w_dof_indices.resize(n);
    Nvec.zero(); vec.zero();
    
    for (FESystemUInt i=0; i<n; i++)
    {
        v_dof_indices[i] =   i;
        w_dof_indices[i] = n+i;
    }
    
    FESystemAssert2(s == 4*n, FESystem::Exception::DimensionsDoNotMatch, s, 4*n);
    
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));

        this->finite_element->getShapeFunction(*(q_pts[i]), Nvec);
        Nvec.scale(q_weight[i]*jac*v_val);
        vec.addVal(v_dof_indices, Nvec); // add to the v-forces
        
        this->finite_element->getShapeFunction(*(q_pts[i]), Nvec);
        Nvec.scale(q_weight[i]*jac*w_val);
        vec.addVal(w_dof_indices, Nvec); // add to the v-forces
    }
}




void
FESystem::Structures::LinearBeamElementBase::calculateInertiaOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n);
    Nvec.zero(); B_mat.zero();
    this->finite_element->getShapeFunction(pt, Nvec);
    
    FESystemAssert4(((s.first == 4) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 1, 4*n, s.first, s.second);
    
    B_mat.setRowVals(0,   0,   n-1, Nvec); // v
    B_mat.setRowVals(1,   n, 2*n-1, Nvec); // w
    B_mat.setRowVals(2, 2*n, 3*n-1, Nvec); // theta_y
    B_mat.setRowVals(3, 3*n, 4*n-1, Nvec); // theta_z
}



void
FESystem::Structures::LinearBeamElementBase::getMaterialMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == 4) && (s.second == 4*n)), FESystem::Numerics::MatrixSizeMismatch, 4, 4*n, s.first, s.second);
    mat.setVal(0, 0, this->rho_val * this->area_val);
    mat.setVal(1, 1, this->rho_val * this->area_val);
    mat.setVal(3, 3, 1.0e-12 * this->rho_val * this->area_val);
    mat.setVal(4, 4, 1.0e-12 * this->rho_val * this->area_val);
}




