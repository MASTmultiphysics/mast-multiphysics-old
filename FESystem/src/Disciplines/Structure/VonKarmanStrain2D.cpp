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
#include "Disciplines/Structure/LinearPlateElementBase.h"


FESystem::Structures::VonKarmanStrain2D::VonKarmanStrain2D():
FESystem::Structures::Structural2DElementBase(),
if_constant_extension_stress(false),
inplane_stress(NULL),
membrane_elem(NULL),
plate_elem(NULL)
{
    this->inplane_stress = new FESystem::Numerics::DenseMatrix<FESystemDouble>;
    this->inplane_stress->resize(2, 2);
}


FESystem::Structures::VonKarmanStrain2D::~VonKarmanStrain2D()
{
    
}



FESystemUInt
FESystem::Structures::VonKarmanStrain2D::getNElemDofs() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    // get the number of degrees of freedom in the membrane element
    if (this->if_constant_extension_stress)
        return this->plate_elem->getNElemDofs();
    else
        return this->plate_elem->getNElemDofs() + this->membrane_elem->getNElemDofs();
}



void
FESystem::Structures::VonKarmanStrain2D::clear()
{
    FESystem::Structures::Structural2DElementBase::clear();
    this->if_constant_extension_stress = false;
    this->inplane_stress->zero();
    this->membrane_elem = NULL;
    this->plate_elem = NULL;
}



void
FESystem::Structures::VonKarmanStrain2D::getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec)
{
    static std::vector<FESystemUInt> plate_vec, membrane_vec;
    
    if (if_constant_extension_stress)
        this->plate_elem->getActiveElementMatrixIndices(vec);
    else
    {
        this->membrane_elem->getActiveElementMatrixIndices(membrane_vec);
        this->plate_elem->getActiveElementMatrixIndices(plate_vec);
        
        vec.resize(plate_vec.size() + membrane_vec.size());
        for (FESystemUInt i=0; i<membrane_vec.size(); i++)
            vec[i] = membrane_vec[i];
        for (FESystemUInt i=0; i<plate_vec.size(); i++)
            vec[i+membrane_vec.size()] = plate_vec[i];
    }
}



void
FESystem::Structures::VonKarmanStrain2D::getStressTensor(const FESystem::Geometry::Point& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                                         FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




void
FESystem::Structures::VonKarmanStrain2D::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                    FESystem::Structures::Membrane& mem, FESystem::Structures::LinearPlateElementBase& plt)
{
    FESystemAssert0(mem.E_val == plt.E_val, FESystem::Exception::InvalidValue);
    FESystemAssert0(mem.getThickness() == plt.getThickness(), FESystem::Exception::InvalidValue);
    FESystemAssert0(mem.rho_val == mem.rho_val, FESystem::Exception::InvalidValue);
    
    FESystem::Structures::Structural2DElementBase::initialize(elem, fe, q_rule, plt.E_val, plt.nu_val, plt.rho_val, plt.getThickness());
    this->if_constant_extension_stress = false;
    this->th_val = plt.getThickness();
    this->membrane_elem = &mem;
    this->plate_elem = &plt;
}



void
FESystem::Structures::VonKarmanStrain2D::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                    FESystem::Numerics::MatrixBase<FESystemDouble>& stress, FESystem::Structures::LinearPlateElementBase& plt)
{
    const std::pair<FESystemUInt, FESystemUInt> s = stress.getSize();
    
    FESystemAssert4((s.first == 2) && (s.second == 2), FESystem::Numerics::MatrixSizeMismatch, 2, 2, s.first, s.first);
    
    FESystem::Structures::Structural2DElementBase::initialize(elem, fe, q_rule, plt.E_val, plt.nu_val, plt.rho_val, plt.getThickness());
    this->if_constant_extension_stress = true;
    this->inplane_stress->copyMatrix(stress);
    this->th_val = plt.getThickness();
    this->membrane_elem = NULL;
    this->plate_elem = &plt;
}



void
FESystem::Structures::VonKarmanStrain2D::calculateInternalForceVector(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    const FESystemUInt n = this->geometric_elem->getNNodes();
    
    FESystemAssert2(sol.getSize() == this->getNElemDofs(), FESystem::Exception::DimensionsDoNotMatch, sol.getSize(), this->getNElemDofs());
    FESystemAssert2(vec.getSize() == this->getNElemDofs(), FESystem::Exception::DimensionsDoNotMatch, vec.getSize(), this->getNElemDofs());
    
    static FESystem::Numerics::LocalVector<FESystemDouble> membrane_dof, w_dof, stress_vec, strain_vec, tmp_vk, tmp_vec1, tmp_vec2, plate_sol, plate_internal_force;
    static FESystem::Numerics::DenseMatrix<FESystemDouble> stress_tensor, B_vk_mat, B_plate_mat, B_membrane_mat, plate_stiff_mat, C_mat;
    static std::vector<FESystemUInt> membrane_dof_indices, plate_dof_indices, w_dof_indices;
    B_membrane_mat.resize(1, 2*n); B_plate_mat.resize(1, n); B_vk_mat.resize(3, 2); stress_tensor.resize(2, 2); C_mat.resize(3, 3);
    stress_vec.resize(3), strain_vec.resize(3); tmp_vec1.resize(3); tmp_vec2.resize(n); tmp_vk.resize(2);
    plate_stiff_mat.resize(this->plate_elem->getNElemDofs(), this->plate_elem->getNElemDofs()); membrane_dof.resize(n); w_dof.resize(n);
    plate_sol.resize(this->plate_elem->getNElemDofs()); plate_internal_force.resize(this->plate_elem->getNElemDofs());
    membrane_dof_indices.resize(2*n); w_dof_indices.resize(n); plate_dof_indices.resize(this->plate_elem->getNElemDofs());
    
    
    
    // get the degree of freedom indices in the element
    if (this->if_constant_extension_stress)
    {
        for (FESystemUInt i=0; i<n; i++)
            w_dof_indices[i] =   i;
        for (FESystemUInt i=0; i<3*n; i++)
            plate_dof_indices[i] = i;
    }
    else
    {
        this->membrane_elem->getActiveElementMatrixIndices(membrane_dof_indices);
        sol.getSubVectorValsFromIndices(membrane_dof_indices, membrane_dof);
        
        for (FESystemUInt i=0; i<n; i++)
            w_dof_indices[i] = 2*n+i;
        
        for (FESystemUInt i=0; i<3*n; i++)
            plate_dof_indices[i] = 2*n+i; // offset by the u,v-dofs
    }
    
    sol.getSubVectorValsFromIndices(w_dof_indices, w_dof);
    sol.getSubVectorValsFromIndices(plate_dof_indices, plate_sol);
    
    
    // get the compliance matrix
    this->getMaterialComplianceMatrix(C_mat);
    
    // calculate the vector values
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble  jac=0.0;
    vec.zero();
    
    // bending contribution
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        
        // get the membrane stress
        if (this->if_constant_extension_stress)
        {
            stress_vec.setVal(0, this->inplane_stress->getVal(0,0)); // sigma_xx
            stress_vec.setVal(1, this->inplane_stress->getVal(0,1)); // sigma_xy
            stress_vec.setVal(2, this->inplane_stress->getVal(1,1)); // sigma_yy
        }
        else
        {
            stress_tensor.zero();
            this->membrane_elem->getStressTensor(*(q_pts[i]), membrane_dof, stress_tensor);
            stress_vec.setVal(0, stress_tensor.getVal(0,0)); // sigma_xx
            stress_vec.setVal(1, stress_tensor.getVal(0,1)); // sigma_xy
            stress_vec.setVal(2, stress_tensor.getVal(1,1)); // sigma_yy
        }
        
        // calculate the vonKarman stress
        this->calculateTransverseDisplacementOperatorMatrix(*(q_pts[i]), B_plate_mat);
        
        // calculate the extension
        tmp_vk.zero();
        B_plate_mat.rightVectorMultiply(w_dof, tmp_vk); // gives [dw/dx, dw/dy]
        strain_vec.setVal(0,       0.5*pow(tmp_vk.getVal(0),2)); // 1/2  (dw/dx)^2
        strain_vec.setVal(0,       0.5*pow(tmp_vk.getVal(1),2)); // 1/2  (dw/dy)^2
        strain_vec.setVal(0, tmp_vk.getVal(0)*tmp_vk.getVal(1)); // (dw/dx)(dw/dy)
        C_mat.rightVectorMultiply(strain_vec, tmp_vec1); // von Karman strain
        stress_vec.add(1.0, tmp_vec1); // membrane + von Karman strain
        stress_vec.scale(this->th_val); // this gives the force per unit length
        
        // contribution from w-displacement
        B_vk_mat.zero();
        B_vk_mat.setVal(0, 0, tmp_vk.getVal(0)); // dw/dx
        B_vk_mat.setVal(1, 1, tmp_vk.getVal(1)); // dw/dy
        B_vk_mat.setVal(2, 0, tmp_vk.getVal(1)); // dw/dy
        B_vk_mat.setVal(2, 1, tmp_vk.getVal(0)); // dw/dx
        
        B_vk_mat.leftVectorMultiply(stress_vec, strain_vec);
        B_plate_mat.leftVectorMultiply(strain_vec, tmp_vec1);
        tmp_vec1.scale(q_weight[i]*jac);
        vec.addVal(w_dof_indices, tmp_vec1); // B_vk^T sigma

        
        // contribution from membrane strain
        if (!this->if_constant_extension_stress)
        {
            tmp_vec1.zero();
            this->membrane_elem->calculateOperatorMatrix(*(q_pts[i]), B_membrane_mat, true);
            B_membrane_mat.leftVectorMultiply(stress_vec, tmp_vec1);
            tmp_vec1.scale(q_weight[i]*jac);
            vec.addVal(membrane_dof_indices, tmp_vec1); // B_m^T sigma
        }
    }
    
    // plate strain contribution
    sol.getSubVectorValsFromIndices(plate_dof_indices, plate_sol);
    this->plate_elem->calculateStiffnessMatrix(plate_stiff_mat);
    plate_stiff_mat.rightVectorMultiply(plate_sol, plate_internal_force);
    vec.addVal(plate_dof_indices, plate_internal_force);
}



void
FESystem::Structures::VonKarmanStrain2D::calculateTangentStiffnessMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert2(sol.getSize() == this->getNElemDofs(), FESystem::Exception::DimensionsDoNotMatch, sol.getSize(), this->getNElemDofs());
    FESystemAssert4((s.first == this->getNElemDofs()) && (s.second == this->getNElemDofs()), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, this->getNElemDofs(), this->getNElemDofs());
    
    
    static FESystem::Numerics::LocalVector<FESystemDouble> membrane_dof, w_dof, stress_vec, strain_vec, tmp_vk, tmp_vec1, tmp_vec2, plate_sol, plate_internal_force;
    static FESystem::Numerics::DenseMatrix<FESystemDouble> stress_tensor, B_vk_mat, B_plate_mat, B_membrane_mat, plate_stiff_mat, C_mat;
    static std::vector<FESystemUInt> membrane_dof_indices, plate_dof_indices, w_dof_indices;
    B_membrane_mat.resize(1, 2*n); B_plate_mat.resize(1, n); B_vk_mat.resize(3, 2); stress_tensor.resize(2, 2); C_mat.resize(3, 3);
    stress_vec.resize(3), strain_vec.resize(3); tmp_vec1.resize(3); tmp_vec2.resize(n); tmp_vk.resize(2);
    plate_stiff_mat.resize(this->plate_elem->getNElemDofs(), this->plate_elem->getNElemDofs()); membrane_dof.resize(n); w_dof.resize(n);
    plate_sol.resize(this->plate_elem->getNElemDofs()); plate_internal_force.resize(this->plate_elem->getNElemDofs());
    membrane_dof_indices.resize(2*n); w_dof_indices.resize(n); plate_dof_indices.resize(this->plate_elem->getNElemDofs());
    
    
    
    // get the degree of freedom indices in the element
    if (this->if_constant_extension_stress)
    {
        for (FESystemUInt i=0; i<n; i++)
            w_dof_indices[i] =   i;
        for (FESystemUInt i=0; i<3*n; i++)
            plate_dof_indices[i] = i;
    }
    else
    {
        this->membrane_elem->getActiveElementMatrixIndices(membrane_dof_indices);
        sol.getSubVectorValsFromIndices(membrane_dof_indices, membrane_dof);
        
        for (FESystemUInt i=0; i<n; i++)
            w_dof_indices[i] = 2*n+i;
        
        for (FESystemUInt i=0; i<3*n; i++)
            plate_dof_indices[i] = 2*n+i; // offset by the u,v-dofs
    }
    
    sol.getSubVectorValsFromIndices(w_dof_indices, w_dof);
    sol.getSubVectorValsFromIndices(plate_dof_indices, plate_sol);
    
    
    // get the compliance matrix
    this->getMaterialComplianceMatrix(C_mat);
    
    // calculate the vector values
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble  jac=0.0;
    mat.zero();
    
    // bending contribution
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        
        // get the membrane stress
        if (this->if_constant_extension_stress)
        {
            stress_vec.setVal(0, this->inplane_stress->getVal(0,0)); // sigma_xx
            stress_vec.setVal(1, this->inplane_stress->getVal(0,1)); // sigma_xy
            stress_vec.setVal(2, this->inplane_stress->getVal(1,1)); // sigma_yy
        }
        else
        {
            stress_tensor.zero();
            this->membrane_elem->getStressTensor(*(q_pts[i]), membrane_dof, stress_tensor);
            stress_vec.setVal(0, stress_tensor.getVal(0,0)); // sigma_xx
            stress_vec.setVal(1, stress_tensor.getVal(0,1)); // sigma_xy
            stress_vec.setVal(2, stress_tensor.getVal(1,1)); // sigma_yy
        }
        
        // calculate the vonKarman stress
        this->calculateTransverseDisplacementOperatorMatrix(*(q_pts[i]), B_plate_mat);
        
        // calculate the extension
        tmp_vk.zero();
        B_plate_mat.rightVectorMultiply(w_dof, tmp_vk); // gives [dw/dx, dw/dy]
        strain_vec.setVal(0,       0.5*pow(tmp_vk.getVal(0),2)); // 1/2  (dw/dx)^2
        strain_vec.setVal(0,       0.5*pow(tmp_vk.getVal(1),2)); // 1/2  (dw/dy)^2
        strain_vec.setVal(0, tmp_vk.getVal(0)*tmp_vk.getVal(1)); // (dw/dx)(dw/dy)
        C_mat.rightVectorMultiply(strain_vec, tmp_vec1); // von Karman strain
        stress_vec.add(1.0, tmp_vec1); // membrane + von Karman strain
        stress_vec.scale(this->th_val); // this gives the force per unit length

        // set values in the stress tensor
        stress_tensor.setVal(0, 0, stress_vec.getVal(0)); // sigma_xx
        stress_tensor.setVal(1, 0, stress_vec.getVal(2)); // sigma_xy
        stress_tensor.setVal(0, 1, stress_vec.getVal(2)); // sigma_xy
        stress_tensor.setVal(1, 1, stress_vec.getVal(1)); // sigma_yy

        // contribution from w-displacement
        B_vk_mat.zero();
        B_vk_mat.setVal(0, 0, tmp_vk.getVal(0)); // dw/dx
        B_vk_mat.setVal(1, 1, tmp_vk.getVal(1)); // dw/dy
        B_vk_mat.setVal(2, 0, tmp_vk.getVal(1)); // dw/dy
        B_vk_mat.setVal(2, 1, tmp_vk.getVal(0)); // dw/dx
        
        // calculate stiffness matrix contribution due to dsigma_vk / dw : both membrane and von Karman stresses will be needed
        // first calculate the matrix dsigma_vk / dw
//        C_mat.matrixRightMultiply(1.0, B_vk_mat, tmp_mat1); // C B_vk
//        tmp_mat1.matrixRightMultiply(1.0, B_plate_mat, dsigma_dw); // C B_vk B_w
//        
//        // now calculate the contribution from von Karman strain
//        B_vk_mat.matrixTransposeRightMultiply(1.0, dsigma_dw, tmp_mat2); // B_vk^T dsigma/dw
//        B_plate_mat.matrixTransposeRightMultiply(1.0, tmp_mat2, tmp_mat3); // B_w^T B_vk^T dsigma/dw
//        mat.addVal(w_dof_indices, w_dof_indices, tmp_mat3); // B_w^T B_vk^T dsigma/dw
//        
//        // von Karman (w,w) coupling B_w^T sigma B_w
//        stress_tensor.matrixRightMultiply(1.0, B_plate_mat, tmp_mat2);
//        B_plate_mat.matrixTransposeRightMultiply(1.0, tmp_mat2, tmp_mat3);
//        tmp_mat3.scale(q_weight[i]*jac*this->th_val);
//        mat.addVal(w_dof_indices, w_dof_indices, tmp_mat3); // B_w^T sigma B_w
//        
//        // contribution from membrane strain
//        if (!this->if_constant_extension_stress)
//        {
//            // and the stiffness contribution from membrane due to dsigma_vk / dw (u,w)
//            B_membrane_mat.matrixTransposeRightMultiply(1.0, dsigma_dw, tmp_mat2); // B_m^T dsigma/dw
//            mat.addVal(membrane_dof_indices, w_dof_indices, tmp_mat2); // B_w^T B_vk^T dsigma/dw
//            
//            // finally, the stiffness coupling of von Karman and membrane dofs (w,u)
//            C_mat.matrixRightMultiply(1.0, B_membrane_mat, tmp_mat1);
//            B_vk_mat.matrixTransposeRightMultiply(1.0, tmp_mat1, tmp_mat2);
//            B_plate_mat.matrixTransposeRightMultiply(1.0, tmp_mat2, tmp_mat3);
//            tmp_mat3.scale(q_weight[i]*jac*this->th_val);
//            mat.addVal(w_dof_indices, membrane_dof_indices, tmp_mat3); // B_w^T B_vk^T C B_m
//         }
    }
//
//    
//    // contribution from membrane stiffness matrix
//    if (!this->if_constant_extension_stress)
//    {
//        this->membrane_elem->calculateStiffnessMatrix(membrane_stiff_mat);
//        mat.addVal(membrane_dof_indices, membrane_dof_indices, membrane_stiff_mat); // B_m^T C B_m
//    }

    // contribution from plate stiffness matrix
    this->plate_elem->calculateStiffnessMatrix(plate_stiff_mat);
    mat.addVal(plate_dof_indices, plate_dof_indices, plate_stiff_mat);
}



void
FESystem::Structures::VonKarmanStrain2D::calculateTransverseDisplacementOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    FESystemAssert4(((s.first == 2) && (s.second== n)), FESystem::Numerics::MatrixSizeMismatch, 2, n, s.first, s.second);
    
    static std::vector<FESystemUInt> derivatives(2);
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n); Nvec.zero();

    B_mat.zero();

    // get dN/dx
    derivatives[0] = 1; derivatives[1] = 0;
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives, pt, Nvec);
    B_mat.setRowVals(0, 0, n-1, Nvec);

    // get dN/dy
    derivatives[0] = 0; derivatives[1] = 1;
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives, pt, Nvec);
    B_mat.setRowVals(1, 0, n-1, Nvec);
}



void
FESystem::Structures::VonKarmanStrain2D::getMaterialComplianceMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4(((s.first == 3) && (s.second== 3)), FESystem::Numerics::MatrixSizeMismatch, 3, 3, s.first, s.second);
    
    FESystemDouble val = this->E_val/(1.0-this->nu_val*this->nu_val);
    
    mat.setVal(0, 0, val);
    mat.setVal(0, 1, this->nu_val*val);
    mat.setVal(1, 0, this->nu_val*val);
    mat.setVal(1, 1, val);
    mat.setVal(2, 2, this->G_val);
    mat.scale(pow(this->th_val,3)/12.0);
}



