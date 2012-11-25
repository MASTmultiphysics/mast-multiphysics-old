//
//  EulerTest.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/26/12.
//
//

// FESystem include
#include "TestingIncludes.h"

// TBB includes
#include "tbb/tbb.h"
#include "tbb/mutex.h"
#include "tbb/spin_mutex.h"

enum AnalysisCase
{
    AIRFOIL_BUMP,
    RAMP,
    HYPERSONIC_CYLINDER
};



const FESystemDouble  rho=1.05, u1=1400.0, temp = 300.0, cp= 1.003e3, cv = 0.716e3, R=cp-cv, p = R*rho*temp, time_step=1.0e-0, final_t=1.0e6;
const FESystemDouble x_length = 2.0, y_length = 0.5, nonlin_tol = 1.0e-6;
const FESystemDouble t_by_c = 0.02, chord = 0.5, thickness = 0.5*t_by_c*chord, x0=x_length/2-chord/2, x1=x0+chord; // airfoilf data
const FESystemDouble rc = 0.5, rx= 1.5, ry = 3.0, theta = 5.0*PI_VAL/12.0; // hypersonic cylinder data
const FESystemDouble x_init = 0.2, ramp_slope = 0.05; // ramp data
const FESystemUInt nx=100, ny=60, dim = 2, max_nonlin_iters = 3, n_vars=4;
const AnalysisCase case_type = AIRFOIL_BUMP;





void modifyMeshForCase(FESystem::Mesh::MeshBase& mesh)
{
    {
        std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
        std::vector<FESystem::Mesh::Node*>::iterator it=nodes.begin(), end=nodes.end();
        
        FESystemDouble x_val, y_val;
        
        switch (case_type)
        {
            case AIRFOIL_BUMP:
            {
                for ( ; it!=end; it++)
                {
                    if (((*it)->getVal(0) >= x0) && ((*it)->getVal(0) <= x1))
                    {
                        x_val = (*it)->getVal(0);
                        y_val = (*it)->getVal(1);
                        
                        y_val += thickness*(1.0-y_val/y_length)*sin(PI_VAL*(x_val-x0)/chord);
                        
                        (*it)->setVal(1, y_val);
                    }
                }
            }
                break;

            case RAMP:
            {
                for ( ; it!=end; it++)
                {
                    x_val = (*it)->getVal(0);
                    y_val = (*it)->getVal(1);
                    
                    if (x_val >= x_init*x_length)
                        y_val += (x_val-x_init*x_length)*(1.0-y_val/y_length)*ramp_slope;
                    
                    (*it)->setVal(1, y_val);
                }
            }
                break;

            case HYPERSONIC_CYLINDER:
            {
                for ( ; it!=end; it++)
                {
                    x_val = (*it)->getVal(0)/x_length;
                    y_val = (*it)->getVal(1)/y_length;
                    
                    (*it)->setVal(0, -(rx-(rx-rc)*x_val)*cos(theta*(2*y_val-1)));
                    (*it)->setVal(1, (ry-(ry-rc)*x_val)*sin(theta*(2*y_val-1)));
                }
            }
                break;

            default:
                FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, case_type);
                break;
        }
        
    }
    
    
    // now update the elements after mesh modification
    {
        std::vector<FESystem::Mesh::ElemBase*> elems = mesh.getElements();
        std::vector<FESystem::Mesh::ElemBase*>::iterator it=elems.begin(), end=elems.end();
        
        for ( ; it!=end; it++)
            (*it)->updateAfterMeshDeformation();
    }
}




void setBoundaryConditionTag(FESystem::Mesh::MeshBase& mesh, std::set<FESystemUInt>& bc_dofs)
{
    std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    switch (case_type)
    {
        case AIRFOIL_BUMP:
        case RAMP:
        {
            for (FESystemUInt i=0; i<nodes.size(); i++)
            {
                if ((nodes[i]->getVal(0) == 0.0)) // left boundary nodes
                {
                    bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // rho
                    bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // rho u1
                    bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(2).global_dof_id[0]); // rho u2
                    bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0]); // rho e
                    const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
                    std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
                    for ( ; it != end; it++)
                        if (!(*it)->checkForTag(0))
                            (*it)->setTag(0);
                }
                
                if (fabs(nodes[i]->getVal(0) - x_length) <= FESystem::Base::getMachineEpsilon<FESystemDouble>()) // right boundary nodes
                {
                    const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
                    std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
                    for ( ; it != end; it++)
                        if (!(*it)->checkForTag(1))
                            (*it)->setTag(1);
                }
                
                if ((nodes[i]->getVal(1) == 0.0)) // lower boundary nodes
                {
                    const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
                    std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
                    for ( ; it != end; it++)
                        if (!(*it)->checkForTag(2))
                            (*it)->setTag(2);
                }
                
                
                if (fabs(nodes[i]->getVal(1) - y_length) <= FESystem::Base::getMachineEpsilon<FESystemDouble>()) // upper boundary nodes
                {
                    const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
                    std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
                    for ( ; it != end; it++)
                        if (!(*it)->checkForTag(3))
                            (*it)->setTag(3);
                }
            }
        }
            break;

        case HYPERSONIC_CYLINDER:
        {
            for (FESystemUInt i=0; i<nodes.size(); i++)
            {
                if ((nodes[i]->getVal(0) == 0.0)) // inlet
                {
                    bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0]); // rho
                    bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0]); // rho u1
                    bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(2).global_dof_id[0]); // rho u2
                    bc_dofs.insert(nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0]); // rho e
                    const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
                    std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
                    for ( ; it != end; it++)
                        if (!(*it)->checkForTag(0))
                            (*it)->setTag(0);
                }
                
                if ((nodes[i]->getVal(0) == x_length)) // solid wall
                {
                    const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
                    std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
                    for ( ; it != end; it++)
                        if (!(*it)->checkForTag(1))
                            (*it)->setTag(1);
                }

                if ((nodes[i]->getVal(1) == 0.0)) // outlet
                {
                    const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
                    std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
                    for ( ; it != end; it++)
                        if (!(*it)->checkForTag(2))
                            (*it)->setTag(2);
                }
                
                if ((nodes[i]->getVal(1) == y_length)) // outlet
                {
                    const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
                    std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
                    for ( ; it != end; it++)
                        if (!(*it)->checkForTag(3))
                            (*it)->setTag(3);
                }
            }
        }
            break;

        default:
            FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, case_type);
            break;
    }
}


void evaluateBoundaryConditionData(const FESystem::Mesh::ElemBase& elem, const FESystem::Quadrature::QuadratureBase& q_boundary, FESystem::Fluid::FluidElementBase& fluid_elem,
                                   FESystem::Numerics::MatrixBase<FESystemDouble>& momentum_flux_tensor, FESystem::Numerics::VectorBase<FESystemDouble>& mass_flux, FESystem::Numerics::VectorBase<FESystemDouble>& energy_flux,
                                   FESystemBoolean if_calculate_jacobian,
                                   FESystem::Numerics::VectorBase<FESystemDouble>& tmp_vec, FESystem::Numerics::MatrixBase<FESystemDouble>& tmp_mat,
                                   FESystem::Numerics::VectorBase<FESystemDouble>& bc_vec, FESystem::Numerics::MatrixBase<FESystemDouble>& elem_mat)
{
    tmp_vec.zero();
    tmp_mat.zero();
    
    bc_vec.zero();
    elem_mat.zero();

    switch (case_type)
    {
        case AIRFOIL_BUMP:
        case RAMP:
        {
            if (elem.checkForTag(0)) // left edge
            {
                fluid_elem.calculateFluxBoundaryCondition(3, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, tmp_vec);
                bc_vec.add(1.0, tmp_vec);
            }
            if (elem.checkForTag(1)) // right edge
            {
                fluid_elem.calculateFluxBoundaryConditionUsingLocalSolution(1, q_boundary, tmp_vec);
                bc_vec.add(1.0, tmp_vec);

                if (if_calculate_jacobian)
                {
                    fluid_elem.calculateTangentMatrixForFluxBoundaryConditionUsingLocalSolution(1, q_boundary, tmp_mat);
                    elem_mat.add(1.0, tmp_mat);
                }
            }
            // set the flux value for the lower and upper boundary
            if (elem.checkForTag(2)) // lower edge
            {
                fluid_elem.calculateSolidWallFluxBoundaryCondition(0, q_boundary, tmp_vec);
                bc_vec.add(1.0, tmp_vec);
                
                if (if_calculate_jacobian)
                {
                    fluid_elem.calculateTangentMatrixForSolidWallFluxBoundaryCondition(0, q_boundary, tmp_mat);
                    elem_mat.add(1.0, tmp_mat);
                }
            }
            if (elem.checkForTag(3)) // upper edge
            {
                fluid_elem.calculateSolidWallFluxBoundaryCondition(2, q_boundary, tmp_vec);
                bc_vec.add(1.0, tmp_vec);

                if (if_calculate_jacobian)
                {
                    fluid_elem.calculateTangentMatrixForSolidWallFluxBoundaryCondition(2, q_boundary, tmp_mat);
                    elem_mat.add(1.0, tmp_mat);
                }
            }
        }
            break;

        case HYPERSONIC_CYLINDER:
        {
            if (elem.checkForTag(0)) // inlet
            {
                fluid_elem.calculateFluxBoundaryCondition(3, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, tmp_vec);
                bc_vec.add(1.0, tmp_vec);
            }
            if (elem.checkForTag(1)) // solid wall
            {
                fluid_elem.calculateSolidWallFluxBoundaryCondition(1, q_boundary, tmp_vec);
                bc_vec.add(1.0, tmp_vec);
                
                if (if_calculate_jacobian)
                {
                    fluid_elem.calculateTangentMatrixForSolidWallFluxBoundaryCondition(1, q_boundary, tmp_mat);
                    elem_mat.add(1.0, tmp_mat);
                }
            }
            // set the flux value for the lower and upper boundary
            if (elem.checkForTag(2)) // outlet
            {
                fluid_elem.calculateFluxBoundaryConditionUsingLocalSolution(0, q_boundary, tmp_vec);
                bc_vec.add(1.0, tmp_vec);
                
                if (if_calculate_jacobian)
                {
                    fluid_elem.calculateTangentMatrixForFluxBoundaryConditionUsingLocalSolution(0, q_boundary, tmp_mat);
                    elem_mat.add(1.0, tmp_mat);
                }
            }
            if (elem.checkForTag(3)) // outlet
            {
                fluid_elem.calculateFluxBoundaryConditionUsingLocalSolution(2, q_boundary, tmp_vec);
                bc_vec.add(1.0, tmp_vec);
                
                if (if_calculate_jacobian)
                {
                    fluid_elem.calculateTangentMatrixForFluxBoundaryConditionUsingLocalSolution(2, q_boundary, tmp_mat);
                    elem_mat.add(1.0, tmp_mat);
                }
            }
        }
            break;

        default:
            FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, case_type);
            break;
    }
}



tbb::mutex assembly_mutex;

class AssembleElementMatrices
{
public:
    AssembleElementMatrices(const std::vector<FESystem::Mesh::ElemBase*>& e, FESystemDouble t_step,
                            FESystemUInt n,
                            const FESystem::Base::DegreeOfFreedomMap& dof,
                            const FESystem::Quadrature::QuadratureBase& q,
                            const FESystem::Quadrature::QuadratureBase& q_b,
                            const FESystem::Numerics::VectorBase<FESystemDouble>& s,
                            const FESystem::Numerics::VectorBase<FESystemDouble>& v,
                            const FESystemBoolean if_jacobian,
                            FESystem::Numerics::VectorBase<FESystemDouble>& r,
                            FESystem::Numerics::MatrixBase<FESystemDouble>& stiff,
                            FESystem::Numerics::MatrixBase<FESystemDouble>& mass):
    elems(e),
    n_elem_dofs(n),
    dt(t_step),
    dof_map(dof),
    q_rule(q),
    q_boundary(q_b),
    sol(s),
    vel(v),
    if_calculate_jacobian(if_jacobian),
    residual(r),
    global_stiffness_mat(stiff),
    global_mass_mat(mass)
    {
        
    }
    
    void operator() (const tbb::blocked_range<FESystemUInt>& r) const
    {
        FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat1, elem_mat2, tmp_mat, bc_mat;
        FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, elem_sol, elem_vel, bc_vec, tmp_vec, tmp_vec2;
        elem_mat1.resize(n_elem_dofs, n_elem_dofs); elem_mat2.resize(n_elem_dofs, n_elem_dofs); elem_vec.resize(n_elem_dofs);
        elem_sol.resize(n_elem_dofs); elem_vel.resize(n_elem_dofs); bc_vec.resize(n_elem_dofs); tmp_vec.resize(n_elem_dofs); tmp_vec2.resize(n_elem_dofs);
        tmp_mat.resize(n_elem_dofs, n_elem_dofs); bc_mat.resize(n_elem_dofs, n_elem_dofs);
        
        FESystem::Numerics::LocalVector<FESystemDouble> mass_flux, energy_flux;
        FESystem::Numerics::DenseMatrix<FESystemDouble> momentum_flux_tensor;
        mass_flux.resize(2); energy_flux.resize(2); momentum_flux_tensor.resize(2, 2);
        switch (case_type)
        {
            case AIRFOIL_BUMP:
            case RAMP:
            case HYPERSONIC_CYLINDER:
            {
                // set the flux value for the left and right boundary
                mass_flux.zero(); mass_flux.setVal(0, rho*u1);
                momentum_flux_tensor.zero(); momentum_flux_tensor.setVal(0, 0, rho*u1*u1+p); momentum_flux_tensor.setVal(1, 1, p);
                energy_flux.zero(); energy_flux.setVal(0, rho*u1*(cv*temp+0.5*u1*u1)+p*u1);
            }
                break;

            default:
                FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, case_type);
                break;
        }

        
        FESystem::FiniteElement::FELagrange fe;
        FESystem::Fluid::FluidElementBase fluid_elem;

        for (FESystemUInt i=r.begin(); i!=r.end(); i++)
        {
            fe.clear();
            fe.reinit(*(elems[i]));
            elem_mat1.zero();
            elem_mat2.zero();
            elem_vec.zero();
            
            dof_map.getFromGlobalVector(*(elems[i]), sol, elem_sol);
            dof_map.getFromGlobalVector(*(elems[i]), vel, elem_vel);
            
            fluid_elem.clear();
            fluid_elem.initialize(*(elems[i]), fe, q_rule, dt, cp, cv, elem_sol, elem_vel);
            
            evaluateBoundaryConditionData(*elems[i], q_boundary, fluid_elem, momentum_flux_tensor, mass_flux, energy_flux, if_calculate_jacobian, tmp_vec, tmp_mat,  bc_vec, bc_mat);

            fluid_elem.calculateResidualVector(elem_vec);
            elem_vec.add(1.0, bc_vec);

            if (if_calculate_jacobian)
            {
                fluid_elem.calculateTangentMatrix(elem_mat1, elem_mat2);
                elem_mat1.add(1.0, bc_mat);
            }
            
            {
                tbb::mutex::scoped_lock my_lock(assembly_mutex);
                dof_map.addToGlobalVector(*(elems[i]), elem_vec, residual);
                if (if_calculate_jacobian)
                {
                    dof_map.addToGlobalMatrix(*(elems[i]), elem_mat1, global_stiffness_mat);
                    dof_map.addToGlobalMatrix(*(elems[i]), elem_mat2, global_mass_mat);
                }
            }
        }
    }
    
    void assembleAll() const
    {
        FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat1, elem_mat2, tmp_mat, bc_mat;
        FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, elem_sol, elem_vel, bc_vec, tmp_vec, tmp_vec2;
        elem_mat1.resize(n_elem_dofs, n_elem_dofs); elem_mat2.resize(n_elem_dofs, n_elem_dofs); elem_vec.resize(n_elem_dofs);
        elem_sol.resize(n_elem_dofs); elem_vel.resize(n_elem_dofs); bc_vec.resize(n_elem_dofs); tmp_vec.resize(n_elem_dofs); tmp_vec2.resize(n_elem_dofs);
        tmp_mat.resize(n_elem_dofs, n_elem_dofs); bc_mat.resize(n_elem_dofs, n_elem_dofs);
        
        FESystem::Numerics::LocalVector<FESystemDouble> mass_flux, energy_flux;
        FESystem::Numerics::DenseMatrix<FESystemDouble> momentum_flux_tensor;
        mass_flux.resize(2); energy_flux.resize(2); momentum_flux_tensor.resize(2, 2);
        switch (case_type)
        {
            case AIRFOIL_BUMP:
            case RAMP:
            case HYPERSONIC_CYLINDER:
            {
                // set the flux value for the left and right boundary
                mass_flux.zero(); mass_flux.setVal(0, rho*u1);
                momentum_flux_tensor.zero(); momentum_flux_tensor.setVal(0, 0, rho*u1*u1+p); momentum_flux_tensor.setVal(1, 1, p);
                energy_flux.zero(); energy_flux.setVal(0, rho*u1*(cv*temp+0.5*u1*u1)+p*u1);
            }
                break;
                
            default:
                FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, case_type);
                break;
        }

        FESystem::FiniteElement::FELagrange fe;
        FESystem::Fluid::FluidElementBase fluid_elem;
        
        for (FESystemUInt i=0; i<elems.size(); i++)
        {
            fe.clear();
            fe.reinit(*(elems[i]));
            elem_mat1.zero();
            elem_mat2.zero();
            elem_vec.zero();
            
            dof_map.getFromGlobalVector(*(elems[i]), sol, elem_sol);
            dof_map.getFromGlobalVector(*(elems[i]), vel, elem_vel);
            
            fluid_elem.clear();
            fluid_elem.initialize(*(elems[i]), fe, q_rule, dt, cp, cv, elem_sol, elem_vel);
            
            evaluateBoundaryConditionData(*elems[i], q_boundary, fluid_elem, momentum_flux_tensor, mass_flux, energy_flux, if_calculate_jacobian, tmp_vec, tmp_mat,  bc_vec, bc_mat);
            
            fluid_elem.calculateResidualVector(elem_vec);
            elem_vec.add(1.0, bc_vec);
            
            if (if_calculate_jacobian)
            {
                fluid_elem.calculateTangentMatrix(elem_mat1, elem_mat2);
                elem_mat1.add(1.0, bc_mat);
            }
            
            {
                dof_map.addToGlobalVector(*(elems[i]), elem_vec, residual);
                if (if_calculate_jacobian)
                {
                    dof_map.addToGlobalMatrix(*(elems[i]), elem_mat1, global_stiffness_mat);
                    dof_map.addToGlobalMatrix(*(elems[i]), elem_mat2, global_mass_mat);
                }
            }
        }
    }

    
protected:
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems;
    const FESystemUInt n_elem_dofs;
    const FESystemDouble dt;
    const FESystem::Base::DegreeOfFreedomMap& dof_map;
    const FESystem::Quadrature::QuadratureBase& q_rule;
    const FESystem::Quadrature::QuadratureBase& q_boundary;
    const FESystem::Numerics::VectorBase<FESystemDouble>& sol;
    const FESystem::Numerics::VectorBase<FESystemDouble>& vel;
    const FESystemBoolean if_calculate_jacobian;
    FESystem::Numerics::VectorBase<FESystemDouble>& residual;
    FESystem::Numerics::MatrixBase<FESystemDouble>& global_stiffness_mat;
    FESystem::Numerics::MatrixBase<FESystemDouble>& global_mass_mat;
};




class EvaluatePrimitiveVariables
{
public:
    EvaluatePrimitiveVariables(const std::vector<FESystem::Mesh::ElemBase*>& e,
                               const std::vector<FESystem::Mesh::Node*>& nd,
                               FESystemDouble t_step, FESystemUInt n,
                               const FESystem::Base::DegreeOfFreedomMap& dof,
                               const FESystem::Quadrature::QuadratureBase& q,
                               const FESystem::Numerics::VectorBase<FESystemDouble>& s,
                               const FESystem::Numerics::VectorBase<FESystemDouble>& v,
                               FESystem::Numerics::VectorBase<FESystemDouble>& pr,
                               FESystem::Numerics::VectorBase<FESystemDouble>& ad):
    elems(e),
    nodes(nd),
    n_elem_dofs(n),
    dt(t_step),
    dof_map(dof),
    q_rule(q),
    sol(s),
    vel(v),
    primitive(pr),
    additional(ad)
    {
        
    }
    
    void operator() (const tbb::blocked_range<FESystemUInt>& r) const
    {
        FESystem::Numerics::LocalVector<FESystemDouble> elem_sol, elem_vel, tmp_vec, tmp_vec2;
        elem_sol.resize(n_elem_dofs); elem_vel.resize(n_elem_dofs); tmp_vec.resize(n_vars); tmp_vec2.resize(n_vars);
        
        FESystem::FiniteElement::FELagrange fe;
        FESystem::Fluid::FluidElementBase fluid_elem;
        
        fe.clear();
        fe.reinit(*(elems[0]));
        dof_map.getFromGlobalVector(*(elems[0]), sol, elem_sol);
        dof_map.getFromGlobalVector(*(elems[0]), vel, elem_vel);
        fluid_elem.clear();
        fluid_elem.initialize(*(elems[0]), fe, q_rule, dt, cp, cv, elem_sol, elem_vel);

        FESystemDouble press, entropy;
        
        for (FESystemUInt i=r.begin(); i!=r.end(); i++)
        {
            for (FESystemUInt j=0; j<n_vars; j++)
                tmp_vec.setVal(j, sol.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]));
            fluid_elem.calculatePrimitiveVariableValues(tmp_vec, tmp_vec2, press, entropy);
            for (FESystemUInt j=0; j<n_vars; j++)
                primitive.setVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0], tmp_vec2.getVal(j));
            additional.setVal(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0], press);
            additional.setVal(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0], entropy);
        }
    }
    
protected:
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems;
    const std::vector<FESystem::Mesh::Node*>& nodes;
    const FESystemUInt n_elem_dofs;
    const FESystemDouble dt;
    const FESystem::Base::DegreeOfFreedomMap& dof_map;
    const FESystem::Quadrature::QuadratureBase& q_rule;
    const FESystem::Numerics::VectorBase<FESystemDouble>& sol;
    const FESystem::Numerics::VectorBase<FESystemDouble>& vel;
    FESystem::Numerics::VectorBase<FESystemDouble>& primitive;
    FESystem::Numerics::VectorBase<FESystemDouble>& additional;
};



void calculateEulerQuantities(FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                              const FESystem::Mesh::MeshBase& mesh,
                              const FESystemDouble dt,
                              const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                              const FESystem::Numerics::VectorBase<FESystemDouble>& vel,
                              const FESystemBoolean if_calculate_jacobian,
                              FESystem::Numerics::VectorBase<FESystemDouble>& residual,
                              FESystem::Numerics::MatrixBase<FESystemDouble>& global_stiffness_mat,
                              FESystem::Numerics::MatrixBase<FESystemDouble>& global_mass_mat,
                              FESystem::Numerics::VectorBase<FESystemDouble>& primitive_sol,
                              FESystem::Numerics::VectorBase<FESystemDouble>& additional_sol)
{
    FESystemUInt n_elem_dofs;
    
    n_elem_dofs = n_vars*n_elem_nodes;
    
    // prepare the quadrature rule and FE for the
    FESystem::Quadrature::TrapezoidQuadrature q_rule, q_boundary;
    
    switch (elem_type)
    {
        case FESystem::Mesh::QUAD4:
        case FESystem::Mesh::TRI3:
            q_rule.init(2, 3);  // bending quadrature is higher than the shear quadrature for reduced integrations
            q_boundary.init(1,3);
            break;
            
        case FESystem::Mesh::QUAD9:
        case FESystem::Mesh::TRI6:
            q_rule.init(2, 5);
            q_boundary.init(1, 5);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    residual.zero(); global_stiffness_mat.zero(); global_mass_mat.zero();
    primitive_sol.zero(); additional_sol.zero();
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    tbb::parallel_for(tbb::blocked_range<FESystemUInt>(0, elems.size()),
                      AssembleElementMatrices(elems, dt, n_elem_dofs, dof_map, q_rule, q_boundary, sol, vel, if_calculate_jacobian, residual, global_stiffness_mat, global_mass_mat));
    
    
//    AssembleElementMatrices a(elems, dt, n_elem_dofs, dof_map, q_rule, q_boundary, sol, vel, residual, global_stiffness_mat, global_mass_mat);
//    a.assembleAll();
    
    tbb::parallel_for(tbb::blocked_range<FESystemUInt>(0, nodes.size()),
                      EvaluatePrimitiveVariables(elems, nodes, dt, n_elem_dofs, dof_map, q_rule, sol, vel, primitive_sol, additional_sol));
}


void testJacobian(FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                  const FESystem::Mesh::MeshBase& mesh,
                  const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                  const FESystem::Numerics::VectorBase<FESystemDouble>& vel)
{
    FESystemUInt n_elem_dofs;
    
    n_elem_dofs = 4*n_elem_nodes;
    
    
    FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat1, elem_mat2, momentum_flux_tensor;
    FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, elem_sol, elem_vel, bc_vec, mass_flux, energy_flux, elem_res, vec, elem_sol_delta;
    elem_mat1.resize(n_elem_dofs, n_elem_dofs); elem_mat2.resize(n_elem_dofs, n_elem_dofs); elem_vec.resize(n_elem_dofs);  elem_sol.resize(n_elem_dofs); elem_vel.resize(n_elem_dofs), bc_vec.resize(n_elem_dofs);
    mass_flux.resize(2); energy_flux.resize(2); momentum_flux_tensor.resize(2, 2); vec.resize(n_elem_dofs), elem_sol_delta.resize(n_elem_dofs); elem_res.resize(n_elem_dofs);
    
    // prepare the quadrature rule and FE for the
    FESystem::Quadrature::TrapezoidQuadrature q_rule, q_boundary;
    FESystem::FiniteElement::FELagrange fe;
    FESystem::Fluid::FluidElementBase fluid_elem;
    
    switch (elem_type)
    {
        case FESystem::Mesh::QUAD4:
        case FESystem::Mesh::TRI3:
            q_rule.init(2, 3);  // bending quadrature is higher than the shear quadrature for reduced integrations
            q_boundary.init(1,3);
            break;
            
        case FESystem::Mesh::QUAD9:
        case FESystem::Mesh::TRI6:
            q_rule.init(2, 5);
            q_boundary.init(1, 5);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
        
    const std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    FESystemUInt elem_num=0;
    
    fe.clear();
    fe.reinit(*(elems[elem_num]));
    elem_mat1.zero();
    elem_mat2.zero();
    elem_vec.zero();
    
    dof_map.getFromGlobalVector(*(elems[elem_num]), sol, elem_sol);
    dof_map.getFromGlobalVector(*(elems[elem_num]), vel, elem_vel);

    // set the flux value for the left and right boundary
    mass_flux.zero(); mass_flux.setVal(0, rho*u1);
    momentum_flux_tensor.zero(); momentum_flux_tensor.setVal(0, 0, rho*u1*u1+p); momentum_flux_tensor.setVal(1, 1, p);
    energy_flux.zero(); energy_flux.setVal(0, rho*u1*(cv*temp+0.5*u1*u1)+p*u1);

    
    // calculate the baseline quantities
    elem_res.zero();
    fluid_elem.clear();
    fluid_elem.initialize(*(elems[elem_num]), fe, q_rule, time_step, cp, cv, elem_sol, elem_vel);
    
    if (elems[elem_num]->checkForTag(0)) // left edge
    {
        fluid_elem.calculateFluxBoundaryCondition(3, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
        elem_res.add(1.0, bc_vec);
    }
    if (elems[elem_num]->checkForTag(1)) // right edge
    {
        fluid_elem.calculateFluxBoundaryCondition(1, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
        elem_res.add(1.0, bc_vec);
    }
    // set the flux value for the lower and upper boundary
    if (elems[elem_num]->checkForTag(2)) // lower edge
    {
        fluid_elem.calculateFluxBoundaryCondition(0, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
        elem_res.add(1.0, bc_vec);
    }
    if (elems[elem_num]->checkForTag(3)) // upper edge
    {
        fluid_elem.calculateFluxBoundaryCondition(2, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
        elem_res.add(1.0, bc_vec);
    }
    
    fluid_elem.calculateResidualVector(elem_vec);
    elem_res.add(1.0, elem_vec);
    
    fluid_elem.calculateTangentMatrix(elem_mat1, elem_mat2); // mat1 is the jacobian to be compared
    
    FESystemDouble delta = 0.01;
    
    // now calculate the perturbed quantity
    for (FESystemUInt i=0; i<n_elem_dofs; i++)
    {
        elem_sol_delta.copyVector(elem_sol);
        elem_sol_delta.setVal(i, elem_sol.getVal(i)+delta);
        vec.zero();
        
        fluid_elem.clear();
        fluid_elem.initialize(*(elems[elem_num]), fe, q_rule, time_step, cp, cv, elem_sol_delta, elem_vel);
        
        if (elems[elem_num]->checkForTag(0)) // left edge
        {
            fluid_elem.calculateFluxBoundaryCondition(3, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
            vec.add(1.0, bc_vec);
        }
        if (elems[elem_num]->checkForTag(1)) // right edge
        {
            fluid_elem.calculateFluxBoundaryCondition(1, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
            vec.add(1.0, bc_vec);
        }
        // set the flux value for the lower and upper boundary
        if (elems[elem_num]->checkForTag(2)) // lower edge
        {
            fluid_elem.calculateFluxBoundaryCondition(0, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
            vec.add(1.0, bc_vec);
        }
        if (elems[elem_num]->checkForTag(3)) // upper edge
        {
            fluid_elem.calculateFluxBoundaryCondition(2, q_boundary, mass_flux, momentum_flux_tensor, energy_flux, bc_vec);
            vec.add(1.0, bc_vec);
        }
        
        fluid_elem.calculateResidualVector(elem_vec);
        vec.add(1.0, elem_vec);
        vec.add(-1.0, elem_res);
        vec.scale(1.0/delta);
        
        elem_mat2.setColumnVals(i, 0, n_elem_dofs-1, vec);
    }
    
    elem_mat1.write(std::cout);
    elem_mat2.write(std::cout);
}


void transientEulerAnalysis(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map)
{
    FESystem::LinearSolvers::LUFactorizationLinearSolver<FESystemDouble> linear_solver;
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    FESystem::Numerics::SparseMatrix<FESystemDouble> stiff_mat, mass;
    FESystem::Numerics::LocalVector<FESystemDouble>  sol, vel, primitive_sol, additional_sol;
    
    sol.resize(dof_map.getNDofs()); vel.resize(dof_map.getNDofs());
    stiff_mat.resize(dof_map.getSparsityPattern());  mass.resize(dof_map.getSparsityPattern());
    primitive_sol.resize(dof_map.getNDofs()); additional_sol.resize(dof_map.getNDofs());
    
    // set the initial condition
    sol.zero();
    for (FESystemUInt i=0; i<nodes.size(); i++)
    {
        sol.setVal(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0], rho);
        sol.setVal(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0], rho * u1);
        sol.setVal(nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0], rho * (temp*cv + u1*u1*0.5));
    }

    // initialize the solver
    FESystem::TransientSolvers::NewmarkTransientSolver<FESystemDouble> transient_solver;
    FESystem::Numerics::SparsityPattern ode_sparsity;
    std::vector<FESystemBoolean> ode_order_include(1); ode_order_include[0] = true;
    std::vector<FESystemDouble> int_constants(1); int_constants[0]=1.0;
    transient_solver.initialize(1, dof_map.getNDofs(), int_constants);
    //transient_solver.enableAdaptiveTimeStepping(4, 1.2, 1.0e-1);
    transient_solver.setConvergenceTolerance(nonlin_tol, max_nonlin_iters);
    transient_solver.setActiveJacobianTerm(ode_order_include);
    transient_solver.setMassMatrix(false, &mass);
    
    
    transient_solver.setInitialTimeData(0, time_step, sol);
    
    transient_solver.setJacobianMatrix(stiff_mat);
    transient_solver.setLinearSolver(linear_solver, false);
    transient_solver.setLinearSolverDataStructureReuse(true);
    
    FESystem::OutputProcessor::VtkOutputProcessor output;
    std::fstream output_file;
    std::vector<FESystemUInt> vars(4); vars[0]=0; vars[1]=1; vars[2]=2; vars[3] = 3; //vars[4] = 4; // write all solutions
    
    FESystemUInt n_skip=0, n_count=0, n_write=0;
    FESystem::TransientSolvers::TransientSolverCallBack call_back;
    while (transient_solver.getCurrentTime()<final_t)
    {
        call_back = transient_solver.incrementTimeStep();
        switch (call_back)
        {
            case FESystem::TransientSolvers::TIME_STEP_CONVERGED:
            {
                if (n_count == n_skip)
                {
                    std::stringstream oss;
                    oss << "sol_" << n_write << ".vtk";
                    output_file.open(oss.str().c_str(),std::fstream::out);
                    output.writeMesh(output_file, mesh, dof_map);
                    output.writeSolution(output_file, "Sol", mesh, dof_map, vars, transient_solver.getCurrentStateVector());
                    output.writeSolution(output_file, "Vel", mesh, dof_map, vars, transient_solver.getCurrentStateVelocityVector());
                    output.writeSolution(output_file, "Primitive", mesh, dof_map, vars, primitive_sol);
                    output.writeSolution(output_file, "Additional", mesh, dof_map, vars, additional_sol);
                    additional_sol.setAllVals(transient_solver.getCurrentTime());
                    output.writeSolution(output_file, "TimeValue", mesh, dof_map, vars, additional_sol);
                    
                    output_file.close();
                    
                    n_write++;
                    n_count=0;
                }
                else
                    n_count++;
            }
                break;
                
            case FESystem::TransientSolvers::EVALUATE_X_DOT:
                // the Jacobian is not updated since it is constant with respect to time
            {
                //testJacobian(elem_type, n_elem_nodes, dof_map, mesh, sol, vel, vel_func);
                calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, transient_solver.getCurrentStepSize(),
                                         transient_solver.getCurrentStateVector(), transient_solver.getCurrentStateVelocityVector(), false,
                                         transient_solver.getCurrentVelocityFunctionVector(), transient_solver.getCurrentJacobianMatrix(),
                                         mass, primitive_sol, additional_sol);
            }
                break;

            case FESystem::TransientSolvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN:
                // the Jacobian is not updated since it is constant with respect to time
            {
                //testJacobian(elem_type, n_elem_nodes, dof_map, mesh, sol, vel, vel_func);
                calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, transient_solver.getCurrentStepSize(),
                                         transient_solver.getCurrentStateVector(), transient_solver.getCurrentStateVelocityVector(), true,
                                         transient_solver.getCurrentVelocityFunctionVector(), transient_solver.getCurrentJacobianMatrix(),
                                         mass, primitive_sol, additional_sol);
            }
                break;
                
            default:
                FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
                break;
        }
    }
}



void nonlinearEulerSolution(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const
                            FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                            const FESystem::Numerics::SparsityPattern& nonbc_sparsity_pattern,
                            const std::map<FESystemUInt, FESystemUInt>& old_to_new_id_map, const std::vector<FESystemUInt>& nonbc_dofs)
{
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    FESystem::Numerics::SparseMatrix<FESystemDouble> stiff_mat, mass;
    FESystem::Numerics::LocalVector<FESystemDouble> residual, sol, vel, primitive_sol, additional_sol;
    
    residual.resize(dof_map.getNDofs()); sol.resize(dof_map.getNDofs()); vel.resize(dof_map.getNDofs());
    primitive_sol.resize(dof_map.getNDofs()); additional_sol.resize(dof_map.getNDofs());
    stiff_mat.resize(dof_map.getSparsityPattern()); mass.resize(dof_map.getSparsityPattern());

    
    FESystem::LinearSolvers::LUFactorizationLinearSolver<FESystemDouble> linear_solver;
    FESystem::NonlinearSolvers::NewtonIterationNonlinearSolver<FESystemDouble> nonlinear_solver;
    
    nonlinear_solver.initialize(stiff_mat, linear_solver);
    nonlinear_solver.setConvergenceLimits(20, 1.0e-10);
    
    // for output
    FESystem::OutputProcessor::VtkOutputProcessor output;
    std::fstream output_file;
    std::vector<FESystemUInt> vars(5); vars[0]=0; vars[1]=1; vars[2]=2; vars[3]=3; vars[4]=4; // write all solutions
    
    
    FESystem::NonlinearSolvers::NonlinearSolverCallBack call_back = nonlinear_solver.getCurrentCallBack();
    
    while (call_back != FESystem::NonlinearSolvers::SOLUTION_CONVERGED)
    {
        switch (call_back)
        {
            case FESystem::NonlinearSolvers::WAITING_TO_START:
                // nothing to be done
                break;
                
            case FESystem::NonlinearSolvers::SET_INITIAL_GUESS:
            {
                // initialize solution to a uniform density and velocity value
                FESystem::Numerics::VectorBase<FESystemDouble>& vec = nonlinear_solver.getCurrentSolution();
                vec.zero();
                for (FESystemUInt i=0; i<nodes.size(); i++)
                {
                    vec.setVal(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0], rho); // rho
                    vec.setVal(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0], rho * u1); // rho u
                    vec.setVal(nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0], rho * (temp*cv + u1*u1*0.5)); // rho * (cv * T + u^2/2)
                }
            }
                break;
                
            case FESystem::NonlinearSolvers::SOLUTION_CONVERGED:
            {
                // write the solution
                std::stringstream oss;
                oss << "sol_" << nonlinear_solver.getCurrentIterationNumber() << ".vtk";
                output_file.open(oss.str().c_str(),std::fstream::out);
                output.writeMesh(output_file, mesh, dof_map);
                output.writeSolution(output_file, "Sol", mesh, dof_map, vars, nonlinear_solver.getCurrentSolution());
                output_file.close();
            }
                break;
                
            case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL:
            {
                calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, 1.0, nonlinear_solver.getCurrentSolution(), vel, false,
                                         nonlinear_solver.getResidualVector(), nonlinear_solver.getJacobianMatrix(), mass,
                                         primitive_sol, additional_sol);
            }
                break;
                
            case FESystem::NonlinearSolvers::EVALUATE_JACOBIAN:
            case FESystem::NonlinearSolvers::EVALUATE_RESIDUAL_AND_JACOBIAN:
            {
                calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, 1.0, nonlinear_solver.getCurrentSolution(), vel, true,
                                         nonlinear_solver.getResidualVector(), nonlinear_solver.getJacobianMatrix(), mass,
                                         primitive_sol, additional_sol);
            }
                break;
                
            default:
                break;
        }
        call_back = nonlinear_solver.incrementSolution();
        if (call_back == FESystem::NonlinearSolvers::MAXIMUM_ITERATIONS_REACHED)
        {
            std::cout << "Maximum iterations reached" << std::endl;
            break;
        }
    }
    
    // copy the last solution and output
    sol.copyVector(nonlinear_solver.getCurrentSolution());
    // write the solution for each node
    for (FESystemUInt i=0; i<nodes.size(); i++)
    {
        std::cout << "Node: " << std::setw(8) << i;
        // write location
        for (FESystemUInt j=0; j<dim; j++)
            std::cout << std::setw(15) << nodes[i]->getVal(j);
        for (FESystemUInt j=0; j<3; j++)
            std::cout << std::setw(15) << sol.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]);
        std::cout << std::endl;
    }
}







int euler_analysis_driver(int argc, char * const argv[])
{
    FESystem::Mesh::ElementType elem_type;
    
    // create the mesh object
    FESystem::Mesh::MeshBase mesh;
    
    // create a nx x ny grid of nodes
    FESystemUInt n_elem_nodes, n_elem_dofs;
    FESystem::Geometry::Point origin(3);
    
    elem_type = FESystem::Mesh::QUAD4;
    createPlaneMesh(elem_type, mesh, origin, nx, ny, x_length, y_length, n_elem_nodes, CROSS, true);
    
    n_elem_dofs = 6*n_elem_nodes;
    
    // now add the degrees of freedom
    FESystem::Base::DegreeOfFreedomMap dof_map(mesh);
    std::string name;
    name = "rho"; dof_map.addVariable(name, 0);
    name = "rhou"; dof_map.addVariable(name, 0);
    name = "rhov"; dof_map.addVariable(name, 0);
    //name = "rhow"; dof_map.addVariable(name, 0);
    name = "rhoe"; dof_map.addVariable(name, 0);
    dof_map.reinit(); // distribute the dofs
    
    // create the finite element and initialize the shape functions
    FESystem::Numerics::SparseMatrix<FESystemDouble> global_stiffness_mat, global_mass_mat;
    FESystem::Numerics::LocalVector<FESystemDouble> rhs, sol;
    
    rhs.resize(dof_map.getNDofs()); sol.resize(dof_map.getNDofs());
    global_mass_mat.resize(dof_map.getSparsityPattern());
    global_stiffness_mat.resize(dof_map.getSparsityPattern());
    
    // apply boundary condition and place a load on the last dof
    std::set<FESystemUInt> bc_dofs;
    setBoundaryConditionTag(mesh, bc_dofs);
    
    // now create the vector of ids that do not have bcs
    std::vector<FESystemUInt> nonbc_dofs;
//    for (FESystemUInt i=0; i<dof_map.getNDofs(); i++)
//        if (!bc_dofs.count(i))
//            nonbc_dofs.push_back(i);
    
    // modify mesh after application of boundary condition
    modifyMeshForCase(mesh);

    
//    // prepare a map of the old to new ID
//    std::vector<FESystemUInt>::const_iterator dof_it=nonbc_dofs.begin(), dof_end=nonbc_dofs.end();
//    std::map<FESystemUInt, FESystemUInt> old_to_new_id_map;
//    FESystemUInt n=0;
//    for ( ; dof_it!=dof_end; dof_it++)
//        old_to_new_id_map.insert(std::map<FESystemUInt, FESystemUInt>::value_type(*dof_it, n++));
    
//    FESystem::Numerics::SparsityPattern nonbc_sparsity_pattern;
//    dof_map.getSparsityPattern().initSparsityPatternForNonConstrainedDOFs(nonbc_dofs, old_to_new_id_map, nonbc_sparsity_pattern);
    
    //nonlinearEulerSolution(dim, elem_type, n_elem_nodes, mesh, dof_map, nonbc_sparsity_pattern, old_to_new_id_map, nonbc_dofs);
    transientEulerAnalysis(dim, elem_type, n_elem_nodes, mesh, dof_map);
    
    //exit(1);
    
    return 0;
}





