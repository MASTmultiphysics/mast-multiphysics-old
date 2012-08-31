//
//  VonKarmanStrain1D.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

// FESystem includes
#include "Disciplines/Structure/VonKarmanStrain1D.h"
#include "Base/FESystemExceptions.h"
#include "Mesh/ElemBase.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"
#include "Disciplines/Structure/ExtensionBar.h"
#include "Disciplines/Structure/LinearBeamElementBase.h"


FESystem::Structures::VonKarmanStrain1D::VonKarmanStrain1D():
FESystem::Structures::Structural1DElementBase(),
if_constant_extension_stress(false),
extension_stress(0.0),
bar_elem(NULL),
beam_elem(NULL)
{
    
}


FESystem::Structures::VonKarmanStrain1D::~VonKarmanStrain1D()
{
    
}



FESystemUInt
FESystem::Structures::VonKarmanStrain1D::getNElemDofs() const
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);

    // get the number of degrees of freedom in the bar element
    if (this->if_constant_extension_stress)
        return this->beam_elem->getNElemDofs();
    else
        return this->beam_elem->getNElemDofs() + this->bar_elem->getNElemDofs();
}



void
FESystem::Structures::VonKarmanStrain1D::clear()
{
    FESystem::Structures::Structural1DElementBase::clear();
    this->if_constant_extension_stress = false;
    this->extension_stress = 0.0;
    this->area_val = 0.0;
    this->bar_elem = NULL;
    this->beam_elem = NULL;
}



void
FESystem::Structures::VonKarmanStrain1D::getActiveElementMatrixIndices(std::vector<FESystemUInt>& vec)
{
    static std::vector<FESystemUInt> beam_vec, bar_vec;

    if (if_constant_extension_stress)
        this->beam_elem->getActiveElementMatrixIndices(vec);
    else
    {
        this->bar_elem->getActiveElementMatrixIndices(bar_vec);
        this->beam_elem->getActiveElementMatrixIndices(beam_vec);

        vec.resize(beam_vec.size() + bar_vec.size());
        for (FESystemUInt i=0; i<bar_vec.size(); i++)
            vec[i] = bar_vec[i];
        for (FESystemUInt i=0; i<beam_vec.size(); i++)
            vec[i+bar_vec.size()] = beam_vec[i];
    }
}



void
FESystem::Structures::VonKarmanStrain1D::getStressTensor(const FESystem::Numerics::VectorBase<FESystemDouble>& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                                      FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(false, FESystem::Exception::InvalidFunctionCall);
}




void
FESystem::Structures::VonKarmanStrain1D::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                    FESystem::Structures::ExtensionBar& bar, FESystem::Structures::LinearBeamElementBase& beam)
{
    FESystemAssert0(bar.E_val == beam.E_val, FESystem::Exception::InvalidValue);
    FESystemAssert0(bar.getArea() == beam.getArea(), FESystem::Exception::InvalidValue);
    FESystemAssert0(bar.rho_val == beam.rho_val, FESystem::Exception::InvalidValue);
    
    FESystem::Structures::Structural1DElementBase::initialize(elem, fe, q_rule, beam.nu_val, beam.nu_val, beam.rho_val, beam.getArea());
    this->if_constant_extension_stress = false;
    this->area_val = beam.getArea();
    this->bar_elem = &bar;
    this->beam_elem = &beam;
}



void
FESystem::Structures::VonKarmanStrain1D::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                                    FESystemDouble stress, FESystem::Structures::LinearBeamElementBase& beam)
{
    FESystem::Structures::Structural1DElementBase::initialize(elem, fe, q_rule, beam.nu_val, beam.nu_val, beam.rho_val, beam.getArea());
    this->if_constant_extension_stress = true;
    this->extension_stress = stress;
    this->area_val = beam.getArea();
    this->bar_elem = NULL;
    this->beam_elem = &beam;
}



void
FESystem::Structures::VonKarmanStrain1D::calculateInternalForceVector(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, FESystem::Numerics::VectorBase<FESystemDouble>& vec)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    const FESystemUInt n = this->geometric_elem->getNNodes();

    FESystemAssert2(sol.getSize() == this->getNElemDofs(), FESystem::Exception::DimensionsDoNotMatch, sol.getSize(), this->getNElemDofs());
    FESystemAssert2(vec.getSize() == this->getNElemDofs(), FESystem::Exception::DimensionsDoNotMatch, vec.getSize(), this->getNElemDofs());

    static FESystem::Numerics::LocalVector<FESystemDouble> u_dofs, v_dof, w_dof, strain, tmp_vec1, tmp_vec2, beam_sol, beam_internal_force;
    static FESystem::Numerics::DenseMatrix<FESystemDouble> stress_tensor, B_beam_mat, B_bar_mat, beam_stiff_mat;
    static std::vector<FESystemUInt> u_dof_indices, beam_dof_indices, v_dof_indices, w_dof_indices;
    B_bar_mat.resize(1, n); B_beam_mat.resize(1, n); strain.resize(1); tmp_vec1.resize(1); tmp_vec2.resize(n);
    beam_stiff_mat.resize(this->beam_elem->getNElemDofs(), this->beam_elem->getNElemDofs()); v_dof.resize(n); w_dof.resize(n); v_dof_indices.resize(n); w_dof_indices.resize(n);
    beam_sol.resize(this->beam_elem->getNElemDofs()); beam_internal_force.resize(this->beam_elem->getNElemDofs());
    beam_dof_indices.resize(4*n);

    
    // get the degree of freedom indices in the element
    if (this->if_constant_extension_stress)
    {
        for (FESystemUInt i=0; i<n; i++)
        {
            v_dof_indices[i] =   i;
            w_dof_indices[i] = n+i;
            
        }
        for (FESystemUInt i=0; i<4*n; i++)
            beam_dof_indices[i] = i;
    }
    else
    {
        this->bar_elem->getActiveElementMatrixIndices(u_dof_indices);
        u_dofs.resize(this->bar_elem->getNElemDofs()); stress_tensor.resize(1, 1);
        sol.getSubVectorValsFromIndices(u_dof_indices, u_dofs);

        for (FESystemUInt i=0; i<n; i++)
        {
            v_dof_indices[i] =   n+i;
            w_dof_indices[i] = 2*n+i;
        }
        for (FESystemUInt i=0; i<4*n; i++)
            beam_dof_indices[i] = n+i; // offset by the u-dofs
    }
    
    sol.getSubVectorValsFromIndices(v_dof_indices, v_dof);
    sol.getSubVectorValsFromIndices(w_dof_indices, w_dof);
    
    
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    // resize the quantities, and get the stress values
    FESystemDouble u_stress = 0.0, v_nonlinear_stress = 0.0, w_nonlinear_stress = 0.0, jac=0.0;
    vec.zero();
    
    // bending contribution
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));

        // get the bar stress
        if (this->if_constant_extension_stress)
            u_stress = this->extension_stress;
        else
        {
            this->bar_elem->getStressTensor(*(q_pts[i]), u_dofs, stress_tensor);
            u_stress = stress_tensor.getVal(0, 0);
        }
        
        // calculate the extension
        this->calculateTransverseDisplacementOperatorMatrix((*(q_pts[i])), B_beam_mat);
        // nonlinear strain due to v and w displacement
        strain.zero();
        B_beam_mat.rightVectorMultiply(v_dof, strain);
        v_nonlinear_stress = 0.5 * this->E_val * pow(strain.getVal(0),2); // strain due to v-displacement
        
        strain.zero();
        B_beam_mat.rightVectorMultiply(w_dof, strain);
        w_nonlinear_stress = 0.5 * this->E_val * pow(strain.getVal(0),2); // strain due to w-displacement
        
        // this is the total stress (linear + nonlinear)
        tmp_vec1.setVal(0, (u_stress + v_nonlinear_stress + w_nonlinear_stress)*this->area_val);
        

        // contribution from bar strain
        if (!this->if_constant_extension_stress)
        {
            tmp_vec2.zero();
            this->bar_elem->calculateOperatorMatrix(*(q_pts[i]), B_bar_mat, true);
            B_bar_mat.leftVectorMultiply(tmp_vec1, tmp_vec2);
            tmp_vec2.scale(q_weight[i]*jac);
            vec.addVal(u_dof_indices, tmp_vec2); // B_m^T sigma
        }
        
        
        // contribution from v-displacement
        tmp_vec2.zero();
        B_beam_mat.leftVectorMultiply(tmp_vec1, tmp_vec2);
        tmp_vec2.scale(q_weight[i]*jac);
        vec.addVal(v_dof_indices, tmp_vec2); // B_v^T sigma
        vec.addVal(w_dof_indices, tmp_vec2); // B_w^T sigma
    }

    // beam strain contribution
    this->beam_elem->calculateStiffnessMatrix(beam_stiff_mat);
    beam_stiff_mat.rightVectorMultiply(beam_sol, beam_internal_force);
    vec.addVal(beam_dof_indices, beam_internal_force);
}



void
FESystem::Structures::VonKarmanStrain1D::calculateTangentStiffnessMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert2(sol.getSize() == this->getNElemDofs(), FESystem::Exception::DimensionsDoNotMatch, sol.getSize(), this->getNElemDofs());
    FESystemAssert4((s.first == this->getNElemDofs()) && (s.second == this->getNElemDofs()), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, this->getNElemDofs(), this->getNElemDofs());
    
    
    static FESystem::Numerics::LocalVector<FESystemDouble> u_dofs, v_dof, w_dof, strain, tmp_vec1, tmp_vec2;
    static FESystem::Numerics::DenseMatrix<FESystemDouble> stress_tensor, B_beam_mat, B_bar_mat, beam_stiff_mat, tmp_mat1, tmp_mat2;
    static std::vector<FESystemUInt> u_dof_indices, beam_dof_indices, v_dof_indices, w_dof_indices;
    B_bar_mat.resize(1, n); B_beam_mat.resize(1, n); strain.resize(1);
    tmp_mat1.resize(1, n); tmp_mat1.resize(n,n); tmp_vec1.resize(n); tmp_vec2.resize(n);
    beam_stiff_mat.resize(this->beam_elem->getNElemDofs(), this->beam_elem->getNElemDofs()); v_dof.resize(n); w_dof.resize(n); v_dof_indices.resize(n); w_dof_indices.resize(n);
    
    // get the degrees of freedom indices in the element
    if (this->if_constant_extension_stress)
    {
        for (FESystemUInt i=0; i<n; i++)
        {
            v_dof_indices[i] =   i;
            w_dof_indices[i] = n+i;
            
        }
        for (FESystemUInt i=0; i<4*n; i++)
            beam_dof_indices[i] = i;
    }
    else
    {
        this->bar_elem->getActiveElementMatrixIndices(u_dof_indices);
        u_dofs.resize(this->bar_elem->getNElemDofs()); stress_tensor.resize(1, 1);
        sol.getSubVectorValsFromIndices(u_dof_indices, u_dofs);

        for (FESystemUInt i=0; i<n; i++)
        {
            v_dof_indices[i] =   n+i;
            w_dof_indices[i] = 2*n+i;
        }
        for (FESystemUInt i=0; i<4*n; i++)
            beam_dof_indices[i] = n+i; // offset by the u-dofs
    }
    
    sol.getSubVectorValsFromIndices(v_dof_indices, v_dof);
    sol.getSubVectorValsFromIndices(w_dof_indices, w_dof);
    

    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    // resize the quantities, and get the stress values
    FESystemDouble u_stress = 0.0, v_nonlinear_stress = 0.0, w_nonlinear_stress = 0.0, dv_dx = 0.0, dw_dx = 0.0, jac=0.0;
    mat.zero();
    
    // numerical quadrature of the tangent stiffness matrix
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        
        // get the bar stress
        if (this->if_constant_extension_stress)
            u_stress = this->extension_stress;
        else
        {
            this->bar_elem->getStressTensor(*(q_pts[i]), u_dofs, stress_tensor);
            u_stress = stress_tensor.getVal(0, 0);
        }
        
        // calculate the extension
        this->calculateTransverseDisplacementOperatorMatrix((*(q_pts[i])), B_beam_mat);
        // nonlinear strain due to v and w displacement
        strain.zero();
        B_beam_mat.rightVectorMultiply(v_dof, strain);
        dv_dx = strain.getVal(0);
        v_nonlinear_stress = 0.5 * this->E_val * pow(strain.getVal(0),2); // strain due to v-displacement
        
        strain.zero();
        B_beam_mat.rightVectorMultiply(w_dof, strain);
        dw_dx = strain.getVal(0);
        w_nonlinear_stress = 0.5 * this->E_val * pow(strain.getVal(0),2); // strain due to w-displacement
        
        // this is the total stress (linear + nonlinear)
        tmp_vec1.setVal(0, (u_stress + v_nonlinear_stress + w_nonlinear_stress)*this->area_val);
        
        
        // contribution from bar strain
        if (!this->if_constant_extension_stress)
        {
            // B_m^T sigma
            tmp_vec2.zero();
            this->bar_elem->calculateOperatorMatrix(*(q_pts[i]), B_bar_mat, true);
            // add sensitivity due to bar dofs: same as bar stiffness matrix
            // B_m^T A E B_m
            tmp_mat2.zero();
            B_bar_mat.matrixTransposeRightMultiply(this->area_val*this->E_val*q_weight[i]*jac, B_bar_mat, tmp_mat2);
            mat.addVal(u_dof_indices, u_dof_indices, tmp_mat2); // add contribution from u
            // add sensitivity due to beam dofs: first v
            // B_m^T A E B_v
            tmp_mat2.zero();
            B_bar_mat.matrixTransposeRightMultiply(dv_dx*this->area_val*this->E_val*q_weight[i]*jac, B_beam_mat, tmp_mat2);
            mat.addVal(u_dof_indices, v_dof_indices, tmp_mat2); // add contribution from v
            // add sensitivity due to beam dofs: now w
            // B_m^T A E B_w
            tmp_mat2.zero();
            B_bar_mat.matrixTransposeRightMultiply(dw_dx*this->area_val*this->E_val*q_weight[i]*jac, B_beam_mat, tmp_mat2);
            mat.addVal(u_dof_indices, w_dof_indices, tmp_mat2); // add contribution from w
            
            // now the dsigma_du term for the beam dofs: first v
            // B_v^T E A dv/dx B_m
            tmp_mat2.zero();
            B_beam_mat.matrixTransposeRightMultiply(dv_dx*this->area_val*this->E_val*q_weight[i]*jac, B_bar_mat, tmp_mat2);
            mat.addVal(v_dof_indices, u_dof_indices, tmp_mat2); // add contribution from u
            
            // the dsigma_du term for the beam dofs: now w
            // B_w^T E A dw/dx B_m
            tmp_mat2.zero();
            B_beam_mat.matrixTransposeRightMultiply(dw_dx*this->area_val*this->E_val*q_weight[i]*jac, B_bar_mat, tmp_mat2);
            mat.addVal(w_dof_indices, u_dof_indices, tmp_mat2); // add contribution from u
        }
        
        
        // contribution from v-displacement
        // B_v^T E A dv/dx B_v
        tmp_mat2.zero();
        B_beam_mat.matrixTransposeRightMultiply(pow(dv_dx,2)*this->area_val*this->E_val*q_weight[i]*jac, B_beam_mat, tmp_mat2);
        mat.addVal(v_dof_indices, v_dof_indices, tmp_mat2); // add contribution from v
        
        // B_w^T E A (dw/dx) (dv/dx) B_v
        tmp_mat2.zero();
        B_beam_mat.matrixTransposeRightMultiply(dw_dx*dv_dx*this->area_val*this->E_val*q_weight[i]*jac, B_beam_mat, tmp_mat2);
        mat.addVal(w_dof_indices, v_dof_indices, tmp_mat2); // add contribution from w

        // B_v^t A sigma B_v
        tmp_mat2.zero();
        B_beam_mat.matrixTransposeRightMultiply((u_stress+v_nonlinear_stress+w_nonlinear_stress)*this->area_val*q_weight[i]*jac, B_beam_mat, tmp_mat2);
        mat.addVal(v_dof_indices, v_dof_indices, tmp_mat2); // add contribution from v

        // contribution from w-displacement
        // B_w^T E A (dw/dx)^2 B_w
        tmp_mat2.zero();
        B_beam_mat.matrixTransposeRightMultiply(pow(dw_dx,2)*this->area_val*this->E_val*q_weight[i]*jac, B_beam_mat, tmp_mat2);
        mat.addVal(w_dof_indices, w_dof_indices, tmp_mat2); // add contribution from w

        // B_v^T E A (dv/dx) (dw/dx) B_w
        tmp_mat2.zero();
        B_beam_mat.matrixTransposeRightMultiply(dw_dx*dv_dx*this->area_val*this->E_val*q_weight[i]*jac, B_beam_mat, tmp_mat2);
        mat.addVal(v_dof_indices, w_dof_indices, tmp_mat2); // add contribution from w

        // B_w^t A sigma B_w
        tmp_mat2.zero();
        B_beam_mat.matrixTransposeRightMultiply((u_stress+v_nonlinear_stress+w_nonlinear_stress)*this->area_val*q_weight[i]*jac, B_beam_mat, tmp_mat2);
        mat.addVal(w_dof_indices, w_dof_indices, tmp_mat2); // add contribution from w
    }

    // beam stiffness matrix contribution
    this->beam_elem->calculateStiffnessMatrix(beam_stiff_mat);
    mat.addVal(beam_dof_indices, beam_dof_indices, beam_stiff_mat);

}



void
FESystem::Structures::VonKarmanStrain1D::calculateTransverseDisplacementOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat)
{
    const FESystemUInt n = this->geometric_elem->getNNodes();
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    FESystemAssert4(((s.first == 1) && (s.second== n)), FESystem::Numerics::MatrixSizeMismatch, 1, n, s.first, s.second);
    
    static std::vector<FESystemUInt> derivatives(1);
    derivatives[0] = 1;
    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n); Nvec.zero();
    B_mat.zero();
    
    this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives, pt, Nvec);
    
    B_mat.setRowVals(0, 0, n-1, Nvec);
}


