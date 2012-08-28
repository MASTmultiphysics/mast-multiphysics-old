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


FESystem::Structures::LinearBeamElementBase::LinearBeamElementBase():
FESystem::Structures::Structural1DElementBase(),
if_include_lateral_stiffness(false),
I_tr_val(0.0),
I_ch_val(0.0),
area_val(0.0)
{
    
}
            


FESystem::Structures::LinearBeamElementBase::~LinearBeamElementBase()
{
    
}


FESystemUInt
FESystem::Structures::LinearBeamElementBase::getNElemDofs() const
{
    if (this->if_include_lateral_stiffness)
        return 4*this->geometric_elem->getNNodes();
    else
        return 2*this->geometric_elem->getNNodes();
}



FESystemBoolean
FESystem::Structures::LinearBeamElementBase::ifIncludeChordwiseStiffness()
{
    return this->if_include_lateral_stiffness;
}
            




void
FESystem::Structures::LinearBeamElementBase::clear()
{
    FESystem::Structures::StructuralElementBase::clear();
    this->if_include_lateral_stiffness = false;
    this->I_tr_val = 0.0;
    this->I_ch_val = 0.0;
}


void
FESystem::Structures::LinearBeamElementBase::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                        FESystemDouble E, FESystemDouble nu, FESystemDouble rho, FESystemDouble I_tr, FESystemDouble I_ch,  FESystemDouble A, FESystemBoolean if_lateral)
{
    FESystem::Structures::StructuralElementBase::initialize(elem,fe,q_rule,E,nu,rho);
    this->if_include_lateral_stiffness = if_lateral;
    this->I_tr_val = I_tr;
    this->I_ch_val = I_ch;
    this->area_val = A;
}
            



void
FESystem::Structures::LinearBeamElementBase::getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec)
{
    FESystemUInt n = this->geometric_elem->getNNodes();
    
    if (this->if_include_lateral_stiffness)
    {
        vec.resize(4*n);
        for (FESystemUInt i=0; i<4*n; i++) vec[i] = n + i; // v-, w-displacement, theta-x and theta-y
    }
    else
    {
        vec.resize(2*n);
        for (FESystemUInt i=0; i<n; i++) vec[i] = 2*n + i; // w-displacement
        for (FESystemUInt i=0; i<n; i++) vec[i+n] = 4*n + i; // theta-y
    }
    
}




void
FESystem::Structures::LinearBeamElementBase::transformMatrixToGlobalSystem(const FESystem::Numerics::MatrixBase<FESystemDouble> &elem_mat, FESystem::Numerics::MatrixBase<FESystemDouble> &global_mat)
{
    FESystemUInt n = this->geometric_elem->getNNodes(), ndof_elem=this->getNElemDofs();
    const std::pair<FESystemUInt, FESystemUInt> s = elem_mat.getSize(), s_g = global_mat.getSize();
    
    FESystemAssert4((s.first == ndof_elem) && (s.second == ndof_elem), FESystem::Numerics::MatrixSizeMismatch, ndof_elem, ndof_elem, s.first, s.second);
    FESystemAssert4((s_g.first == n*6) && (s_g.second == n*6), FESystem::Numerics::MatrixSizeMismatch, n*6, n*6, s_g.first, s_g.second);
    
    static std::vector<FESystemUInt> indices;
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> T_mat, tmp_mat;
    T_mat.resize(n*6, n*6); T_mat.zero();
    tmp_mat.resize(n*6, n*6); tmp_mat.zero();
    
    this->calculateDeformationTransformationMatrix(T_mat);
    
    this->getActiveElementMatrixIndices(indices);
    
    global_mat.zero();
    global_mat.addVal(indices, indices, elem_mat);
    
    // the transformation matrix as returned by the coordinate system is actually T^T. Hence, instead of
    // T^T K T, the operation performed here is T K T^T
    global_mat.matrixRightMultiplyTranspose(1.0, T_mat, tmp_mat);
    global_mat.zero();
    T_mat.matrixRightMultiply(1.0, tmp_mat, global_mat);
    
    //    global_mat.write(std::cout);
}



void
FESystem::Structures::LinearBeamElementBase::transformVectorToGlobalSystem(const FESystem::Numerics::VectorBase<FESystemDouble> &elem_vec, FESystem::Numerics::VectorBase<FESystemDouble> &global_vec)
{
    FESystemUInt n = this->geometric_elem->getNNodes();
    
    static std::vector<FESystemUInt> indices;
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> T_mat;
    static FESystem::Numerics::LocalVector<FESystemDouble> tmp_vec;
    T_mat.resize(n*6, n*6); T_mat.zero();
    tmp_vec.resize(n*6); tmp_vec.zero();
    
    this->calculateDeformationTransformationMatrix(T_mat);
    
    this->getActiveElementMatrixIndices(indices);
    
    tmp_vec.zero();
    tmp_vec.addVal(indices, elem_vec);
    
    T_mat.rightVectorMultiply(tmp_vec, global_vec);
    
}


