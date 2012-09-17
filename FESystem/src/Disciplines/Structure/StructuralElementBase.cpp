//
//  StructuralElementBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/12/12.
//
//


// FESystem includes
#include "Disciplines/Structure/StructuralElementBase.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Geom/CoordinateSystemBase.h"
#include "Geom/Point.h"
#include "Functions/FunctionMappingBase.h"
#include "Numerics/DenseMatrix.h"


FESystem::Structures::StructuralElementBase::StructuralElementBase():
if_initialized(false),
geometric_elem(NULL),
quadrature(NULL),
finite_element(NULL),
E_val(0.0),
nu_val(0.0),
G_val(0.0),
rho_val(0.0)
{
    
}


FESystem::Structures::StructuralElementBase::~StructuralElementBase()
{
    
}


std::string
FESystem::Structures::StructuralElementBase::getVariableName(FESystem::Structures::StructuralVariable var)
{
    std::string name;
    
    switch (var)
    {
        case FESystem::Structures::U_DISP:
            name = "U_DISP";
            break;
            
        case FESystem::Structures::V_DISP:
            name = "V_DISP";
            break;

        case FESystem::Structures::W_DISP:
            name = "W_DISP";
            break;

        case FESystem::Structures::THETA_X:
            name = "THETA_X";
            break;

        case FESystem::Structures::THETA_Y:
            name = "THETA_Y";
            break;

        case FESystem::Structures::THETA_Z:
            name = "THETA_Z";
            break;
        
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
    
    return name;
}


FESystem::Structures::StructuralVariable
FESystem::Structures::StructuralElementBase::getVariableEnum(const std::string &var)
{
    FESystem::Structures::StructuralVariable name;
    
    if (var == "U_DISP")
        name = FESystem::Structures::U_DISP;
    else if (var == "V_DISP")
        name = FESystem::Structures::V_DISP;
    else if (var == "W_DISP")
        name = FESystem::Structures::W_DISP;
    else if (var == "THETA_X")
        name = FESystem::Structures::THETA_X;
    else if (var == "THETA_Y")
        name = FESystem::Structures::THETA_Y;
    else if (var == "THETA_Z")
        name = FESystem::Structures::THETA_Z;
    else
        FESystemAssert0(false, FESystem::Exception::InvalidValue);
    
    return name;
}


void
FESystem::Structures::StructuralElementBase::clear()
{
    this->geometric_elem = NULL;
    this->quadrature = NULL;
    this->finite_element = NULL;
    this->E_val  = 0.0;
    this->nu_val = 0.0;
    this->G_val = 0.0;
    this->rho_val = 0.0;

    this->if_initialized = false;
}



void
FESystem::Structures::StructuralElementBase::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe,
                                                        const FESystem::Quadrature::QuadratureBase& q_rule, FESystemDouble E, FESystemDouble nu, FESystemDouble rho)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidFunctionCall);
    
    this->if_initialized = true;
    this->geometric_elem = &elem;
    this->quadrature = &q_rule;
    this->finite_element = &fe;

    this->E_val  = E;
    this->nu_val = nu;
    this->G_val = E/2.0/(1.0+nu);
    this->rho_val = rho;
}


void
FESystem::Structures::StructuralElementBase::calculateDeformationTransformationMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    static FESystem::Geometry::Point pt(3);
    
    // valid combinations of displacement include
    FESystemUInt n_elem_nodes = this->geometric_elem->getNNodes();
    
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4((s.first == n_elem_nodes*6) && (s.second == n_elem_nodes*6), FESystem::Numerics::MatrixSizeMismatch, 6*n_elem_nodes, 6*n_elem_nodes, s.first, s.second);
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> T_mat_sub;
    T_mat_sub.resize(3,3);
    
    this->geometric_elem->getLocalPhysicalCoordinateSystem().getFunctionMappingObject().getMappingJacobian(pt, T_mat_sub);
    mat.zero();
    
    for (FESystemUInt i=0; i<n_elem_nodes; i++)
        for (FESystemUInt j=0; j<3; j++)
            for (FESystemUInt k=0; k<3; k++)
            {
                mat.setVal(    j*n_elem_nodes+i,     k*n_elem_nodes+i, T_mat_sub.getVal(j,k)); // displacement
                mat.setVal((j+3)*n_elem_nodes+i, (k+3)*n_elem_nodes+i, T_mat_sub.getVal(j,k)); // rotations
            }
}



void
FESystem::Structures::StructuralElementBase::transformMatrixToGlobalSystem(const FESystem::Numerics::MatrixBase<FESystemDouble> &elem_mat, FESystem::Numerics::MatrixBase<FESystemDouble> &global_mat)
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
}



void
FESystem::Structures::StructuralElementBase::transformVectorToGlobalSystem(const FESystem::Numerics::VectorBase<FESystemDouble> &elem_vec, FESystem::Numerics::VectorBase<FESystemDouble> &global_vec)
{
    FESystemUInt n = this->geometric_elem->getNNodes();
    
    static std::vector<FESystemUInt> indices;
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> T_mat;
    static FESystem::Numerics::LocalVector<FESystemDouble> tmp_vec;
    T_mat.resize(n*6, n*6); T_mat.zero();
    tmp_vec.resize(n*6); tmp_vec.zero();
    
    this->calculateDeformationTransformationMatrix(T_mat);
    
    this->getActiveElementMatrixIndices(indices);
    
    global_vec.zero();
    tmp_vec.zero();
    tmp_vec.addVal(indices, elem_vec);
    
    T_mat.rightVectorMultiply(tmp_vec, global_vec);
}




