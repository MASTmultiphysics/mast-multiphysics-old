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
    NACA_AIRFOIL,
    RAMP,
    HYPERSONIC_CYLINDER
};



const FESystemDouble  rho=1.05, u1=277.832, temp = 300.0, cp= 1.003e3, cv = 0.716e3, R=cp-cv, q0 = 0.5*rho*u1*u1, p = R*rho*temp, time_step=1.0e-4, final_t=1.0e6;
const FESystemDouble x_length = 22.0, y_length = 8.00, nonlin_tol = 1.0e-6, fd_delta = 1.0e-7;
const FESystemDouble t_by_c = 0.12, chord = 1.0, thickness = 0.5*t_by_c*chord, x0=x_length/2-chord/2, x1=x0+chord; // airfoilf data
const FESystemDouble rc = 0.5, rx= 1.5, ry = 3.0, theta = 5.0*PI_VAL/12.0; // hypersonic cylinder data
const FESystemDouble x_init = 0.2, ramp_slope = 0.05; // ramp data
const FESystemUInt dim = 2, max_nonlin_iters = 0, n_vars=4, dc_freeze_iter_num = 120;
const AnalysisCase case_type = NACA_AIRFOIL;
const FESystemBoolean if_fd = false;

std::vector<std::vector<FESystemDouble> > dc_vals;

FESystemUInt x_n_divs, y_n_divs;
std::vector<FESystemDouble> x_div_locations, x_relative_mesh_size_in_div, x_points, y_div_locations, y_relative_mesh_size_in_div, y_points;
std::vector<FESystemUInt> x_n_subdivs_in_div, y_n_subdivs_in_div;


void initMeshParameters()
{
    switch (case_type)
    {
        case AIRFOIL_BUMP:
        case NACA_AIRFOIL:
        {
            x_n_divs = 3;
            x_div_locations.resize(x_n_divs+1);
            x_relative_mesh_size_in_div.resize(x_n_divs+1);
            x_n_subdivs_in_div.resize(x_n_divs);
            
            x_div_locations[0] = 0.0;
            x_div_locations[1] = (x_length-chord)/2.0;
            x_div_locations[2] = (x_length+chord)/2.0;
            x_div_locations[3] = x_length;            
            x_relative_mesh_size_in_div[0] = 100.0;
            x_relative_mesh_size_in_div[1] = 1.0;
            x_relative_mesh_size_in_div[2] = 2.0;
            x_relative_mesh_size_in_div[3] = 100.0;
            
            x_n_subdivs_in_div[0] = 30;
            x_n_subdivs_in_div[1] = 40;
            x_n_subdivs_in_div[2] = 15;

            y_n_divs = 1;
            y_div_locations.resize(y_n_divs+1);
            y_relative_mesh_size_in_div.resize(y_n_divs+1);
            y_n_subdivs_in_div.resize(y_n_divs);
            
            y_div_locations[0] = 0.0;
            y_div_locations[1] = y_length;
            
            y_relative_mesh_size_in_div[0] = 1.0;
            y_relative_mesh_size_in_div[1] = 50.0;
            
            y_n_subdivs_in_div[0] = 30;
        }
            break;

        case RAMP:
        {
            x_n_divs = 2;
            x_div_locations.resize(x_n_divs+1);
            x_relative_mesh_size_in_div.resize(x_n_divs+1);
            x_n_subdivs_in_div.resize(x_n_divs);
            
            x_div_locations[0] = 0.0;
            x_div_locations[1] = x_init;
            x_div_locations[2] = x_length;
            
            x_relative_mesh_size_in_div[0] = 30.0;
            x_relative_mesh_size_in_div[1] = 1.0;
            x_relative_mesh_size_in_div[2] = 1.0;
            
            x_n_subdivs_in_div[0] = 20;
            x_n_subdivs_in_div[1] = 100;
            
            y_n_divs = 1;
            y_div_locations.resize(y_n_divs+1);
            y_relative_mesh_size_in_div.resize(y_n_divs+1);
            y_n_subdivs_in_div.resize(y_n_divs);
            
            y_div_locations[0] = 0.0;
            y_div_locations[1] = y_length;
            
            y_relative_mesh_size_in_div[0] = 1.0;
            y_relative_mesh_size_in_div[1] = 20.0;
            
            y_n_subdivs_in_div[0] = 60;
        }
            break;

        case HYPERSONIC_CYLINDER:
        {
            x_n_divs = 1;
            x_div_locations.resize(x_n_divs+1);
            x_relative_mesh_size_in_div.resize(x_n_divs+1);
            x_n_subdivs_in_div.resize(x_n_divs);
            
            x_div_locations[0] = 0.0;
            x_div_locations[1] = x_length;
            
            x_relative_mesh_size_in_div[0] = 20.0;
            x_relative_mesh_size_in_div[1] = 1.0;
            
            x_n_subdivs_in_div[0] = 60;
            
            y_n_divs = 1;
            y_div_locations.resize(y_n_divs+1);
            y_relative_mesh_size_in_div.resize(y_n_divs+1);
            y_n_subdivs_in_div.resize(y_n_divs);
            
            y_div_locations[0] = 0.0;
            y_div_locations[1] = y_length;
            
            y_relative_mesh_size_in_div[0] = 1.0;
            y_relative_mesh_size_in_div[1] = 1.0;
            
            y_n_subdivs_in_div[0] = 80;
        }
            break;

        default:
            break;
    }
}


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
                
            case NACA_AIRFOIL:
            {
                FESystemDouble eta=0.0;
                
                for ( ; it!=end; it++)
                {
                    if (((*it)->getVal(0) >= x0) && ((*it)->getVal(0) <= x1))
                    {
                        x_val = (*it)->getVal(0);
                        y_val = (*it)->getVal(1);
                        
                        eta = (x_val-x0)/chord;
                        
                        y_val += (1.0-y_val/y_length)*2.0*thickness/0.2*(.2969*sqrt(eta)-.1260*eta-.3516*eta*eta+.2843*pow(eta,3)-.1015*pow(eta,4));
                        
                        (*it)->setVal(1, y_val);
                    }
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
        case NACA_AIRFOIL:
        {
            for (FESystemUInt i=0; i<nodes.size(); i++)
            {
                if ((nodes[i]->getVal(0) == 0.0)) // left boundary nodes
                {
                    const std::set<FESystem::Mesh::ElemBase*>& e_set = nodes[i]->getElementConnectivitySet();
                    std::set<FESystem::Mesh::ElemBase*>::const_iterator it = e_set.begin(), end = e_set.end();
                    for ( ; it != end; it++)
                        if (!(*it)->checkForTag(0))
                            (*it)->setTag(0);
                }
                
                if (fabs(nodes[i]->getVal(0) - x_length) <= 100*FESystem::Base::getMachineEpsilon<FESystemDouble>()) // right boundary nodes
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
                    
                    // set tag 4 for airfoil surface
                    it = e_set.begin(), end = e_set.end();
                    if ((nodes[i]->getVal(0)>x0) && (nodes[i]->getVal(0)<x1))
                        for ( ; it != end; it++)
                            if (!(*it)->checkForTag(4))
                                (*it)->setTag(4);
                }
                
                
                if (fabs(nodes[i]->getVal(1) - y_length) <= 100*FESystem::Base::getMachineEpsilon<FESystemDouble>()) // upper boundary nodes
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
                                   const FESystem::Numerics::VectorBase<FESystemDouble>& sol_inf,
                                   const FESystem::Numerics::MatrixBase<FESystemDouble>& momentum_flux_tensor, const FESystem::Numerics::VectorBase<FESystemDouble>& mass_flux,
                                   const FESystem::Numerics::VectorBase<FESystemDouble>& energy_flux, FESystemBoolean if_calculate_jacobian,
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
        case NACA_AIRFOIL:
        {
            if (elem.checkForTag(0)) // left edge
            {
                fluid_elem.calculateMixedBoundaryCondition(3, q_boundary, sol_inf, tmp_vec);
                bc_vec.add(1.0, tmp_vec);
                
                if (if_calculate_jacobian)
                {
                    fluid_elem.calculateTangentMatrixForMixedBoundaryCondition(3, q_boundary, tmp_mat);
                    elem_mat.add(1.0, tmp_mat);
                }
            }
            if (elem.checkForTag(1)) // right edge
            {
                fluid_elem.calculateMixedBoundaryCondition(1, q_boundary, sol_inf, tmp_vec);
                bc_vec.add(1.0, tmp_vec);

                if (if_calculate_jacobian)
                {
                    fluid_elem.calculateTangentMatrixForMixedBoundaryCondition(1, q_boundary, tmp_mat);
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
                //fluid_elem.calculateSolidWallFluxBoundaryCondition(2, q_boundary, tmp_vec);
                fluid_elem.calculateMixedBoundaryCondition(2, q_boundary, sol_inf, tmp_vec);
                bc_vec.add(1.0, tmp_vec);

                if (if_calculate_jacobian)
                {
                    //fluid_elem.calculateTangentMatrixForSolidWallFluxBoundaryCondition(2, q_boundary, tmp_mat);
                    fluid_elem.calculateTangentMatrixForMixedBoundaryCondition(2, q_boundary, tmp_mat);
                    elem_mat.add(1.0, tmp_mat);
                }
            }
        }
            break;

        case HYPERSONIC_CYLINDER:
        {
            if (elem.checkForTag(0)) // inlet
            {
                fluid_elem.calculateMixedBoundaryCondition(3, q_boundary, sol_inf, tmp_vec);
                bc_vec.add(1.0, tmp_vec);
                
                if (if_calculate_jacobian)
                {
                    fluid_elem.calculateTangentMatrixForMixedBoundaryCondition(3, q_boundary, tmp_mat);
                    elem_mat.add(1.0, tmp_mat);
                }
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
                fluid_elem.calculateMixedBoundaryCondition(0, q_boundary, sol_inf, tmp_vec);
                bc_vec.add(1.0, tmp_vec);
                
                if (if_calculate_jacobian)
                {
                    fluid_elem.calculateTangentMatrixForMixedBoundaryCondition(0, q_boundary, tmp_mat);
                    elem_mat.add(1.0, tmp_mat);
                }
            }
            if (elem.checkForTag(3)) // outlet
            {
                fluid_elem.calculateMixedBoundaryCondition(2, q_boundary, sol_inf, tmp_vec);
                bc_vec.add(1.0, tmp_vec);
                
                if (if_calculate_jacobian)
                {
                    fluid_elem.calculateTangentMatrixForMixedBoundaryCondition(2, q_boundary, tmp_mat);
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



extern tbb::mutex assembly_mutex;

class AssembleElementMatrices
{
public:
    AssembleElementMatrices(const std::vector<FESystem::Mesh::ElemBase*>& e, FESystemDouble t_step,
                            FESystemUInt n_nd, FESystemUInt n_dof,
                            const FESystem::Base::DegreeOfFreedomMap& dof,
                            const FESystem::Quadrature::QuadratureBase& q,
                            const FESystem::Quadrature::QuadratureBase& q_b,
                            const FESystem::Numerics::VectorBase<FESystemDouble>& s,
                            const FESystem::Numerics::VectorBase<FESystemDouble>& v,
                            std::vector<std::vector<FESystemDouble> >& dc,
                            const FESystemBoolean if_update_dc,
                            const FESystemBoolean if_jacobian,
                            FESystem::Numerics::VectorBase<FESystemDouble>& r,
                            FESystem::Numerics::MatrixBase<FESystemDouble>& stiff,
                            FESystem::Numerics::MatrixBase<FESystemDouble>& mass):
    elems(e),
    n_elem_nodes(n_nd),
    n_elem_dofs(n_dof),
    dt(t_step),
    dof_map(dof),
    q_rule(q),
    q_boundary(q_b),
    sol(s),
    vel(v),
    elem_dc_vals(dc),
    if_update_elem_dc_vals(if_update_dc),
    if_calculate_jacobian(if_jacobian),
    residual(r),
    global_stiffness_mat(stiff),
    global_mass_mat(mass)
    {
    }
    
    
    void operator() (const tbb::blocked_range<FESystemUInt>& r) const
    {
        this->assemble(r.begin(), r.end());
    }
    

    void assembleAll() const
    {
        this->assemble(0, elems.size());
    }

    
    void assemble(FESystemUInt e1, FESystemUInt e2) const
    {
        FESystem::Numerics::DenseMatrix<FESystemDouble> elem_mat1, elem_mat2, tmp_mat, tmp_mat2, bc_mat;
        FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, elem_sol, elem_vel, bc_vec, tmp_vec, tmp_vec2, delta_sol, delta_res;
        elem_mat1.resize(n_elem_dofs, n_elem_dofs); elem_mat2.resize(n_elem_dofs, n_elem_dofs); elem_vec.resize(n_elem_dofs); tmp_mat2.resize(n_elem_dofs, n_elem_dofs);
        elem_sol.resize(n_elem_dofs); elem_vel.resize(n_elem_dofs); bc_vec.resize(n_elem_dofs); tmp_vec.resize(n_elem_dofs); tmp_vec2.resize(n_elem_dofs);
        tmp_mat.resize(n_elem_dofs, n_elem_dofs); bc_mat.resize(n_elem_dofs, n_elem_dofs); delta_sol.resize(n_elem_dofs); delta_res.resize(n_elem_dofs);
        
        
        FESystem::Numerics::LocalVector<FESystemDouble> mass_flux, energy_flux, sol_inf;
        FESystem::Numerics::DenseMatrix<FESystemDouble> momentum_flux_tensor;
        
        mass_flux.resize(2); energy_flux.resize(2); momentum_flux_tensor.resize(2, 2); sol_inf.resize(n_elem_dofs);
        switch (case_type)
        {
            case AIRFOIL_BUMP:
            case RAMP:
            case HYPERSONIC_CYLINDER:
            case NACA_AIRFOIL:
            {
                // set the flux value for the left and right boundary
                mass_flux.zero(); mass_flux.setVal(0, rho*u1);
                momentum_flux_tensor.zero(); momentum_flux_tensor.setVal(0, 0, rho*u1*u1+p); momentum_flux_tensor.setVal(1, 1, p);
                energy_flux.zero(); energy_flux.setVal(0, rho*u1*(cv*temp+0.5*u1*u1)+p*u1);
                for (FESystemUInt i=0; i<n_elem_nodes; i++)
                {
                    sol_inf.setVal(i, rho);
                    sol_inf.setVal(n_elem_nodes+i, rho*u1);
                    sol_inf.setVal(n_elem_nodes*2+i, 0.0);
                    sol_inf.setVal(n_elem_nodes*3+i, rho*(cv*temp+0.5*u1*u1));
                }
            }
                break;
                
            default:
                FESystemAssert1(false, FESystem::Exception::EnumerationNotHandled, case_type);
                break;
        }

        
        FESystem::FiniteElement::FELagrange fe;
        FESystem::Fluid::FluidElementBase fluid_elem;
        
        for (FESystemUInt i=e1; i!=e2; i++)
        {
            fe.clear();
            fe.reinit(*(elems[i]));
            elem_mat1.zero();
            elem_mat2.zero();
            elem_vec.zero();
            
            dof_map.getFromGlobalVector(*(elems[i]), sol, elem_sol);
            dof_map.getFromGlobalVector(*(elems[i]), vel, elem_vel);
            
            fluid_elem.clear();
            fluid_elem.initialize(*(elems[i]), fe, q_rule, dt, cp, cv, elem_sol, elem_vel, if_update_elem_dc_vals, elem_dc_vals[i]);
            
            if (if_calculate_jacobian && !if_fd)
                evaluateBoundaryConditionData(*elems[i], q_boundary, fluid_elem, sol_inf, momentum_flux_tensor, mass_flux, energy_flux, true, tmp_vec, tmp_mat, bc_vec, bc_mat);
            else
                evaluateBoundaryConditionData(*elems[i], q_boundary, fluid_elem, sol_inf, momentum_flux_tensor, mass_flux, energy_flux, false, tmp_vec, tmp_mat, bc_vec, bc_mat);
            
            fluid_elem.calculateResidualVector(elem_vec);
            elem_vec.add(1.0, bc_vec);
            
            if (if_calculate_jacobian && !if_fd)
            {
                fluid_elem.calculateTangentMatrix(elem_mat1, elem_mat2);
                elem_mat1.add(1.0, bc_mat);
            }
            else if (if_calculate_jacobian && if_fd)
                this->calculateFiniteDifferenceJacobianForElem(*(elems[i]), sol_inf, momentum_flux_tensor, mass_flux, energy_flux,
                                                               elem_sol, elem_vel, elem_vec, elem_dc_vals[i], elem_mat1, delta_sol, delta_res, tmp_mat, tmp_mat, tmp_vec);
            
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
    

    
    void calculateFiniteDifferenceJacobianForElem(const FESystem::Mesh::ElemBase& elem, const FESystem::Numerics::VectorBase<FESystemDouble>& sol_inf,
                                                  const FESystem::Numerics::MatrixBase<FESystemDouble>& momentum_flux_tensor,
                                                  const FESystem::Numerics::VectorBase<FESystemDouble>& mass_flux, const FESystem::Numerics::VectorBase<FESystemDouble>& energy_flux,
                                                  const FESystem::Numerics::VectorBase<FESystemDouble>& elem_sol, const FESystem::Numerics::VectorBase<FESystemDouble>& elem_vel,
                                                  const FESystem::Numerics::VectorBase<FESystemDouble>& res0, std::vector<FESystemDouble>& elem_dc,
                                                  FESystem::Numerics::MatrixBase<FESystemDouble>& jac, FESystem::Numerics::VectorBase<FESystemDouble>& delta_sol,
                                                  FESystem::Numerics::VectorBase<FESystemDouble>& delta_res, FESystem::Numerics::MatrixBase<FESystemDouble>& tmp_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& tmp_mat2,
                                                  FESystem::Numerics::VectorBase<FESystemDouble>& tmp_vec) const
    {
        FESystem::FiniteElement::FELagrange fe;
        FESystem::Fluid::FluidElementBase fluid_elem;
        FESystemDouble val1 = 0.0, delta_val = 0.0;
        fe.clear(); fe.reinit(elem);
        
        for (FESystemUInt i=0; i<elem_sol.getSize(); i++)
        {
            delta_sol.copyVector(elem_sol);
            val1 = delta_sol.getVal(i);
            if (FESystem::Base::comparisonValue<FESystemDouble, FESystemDouble>(val1) > 1.0e-7)
            {
                delta_val = fd_delta * val1;
                val1 *=  (1.0 + fd_delta);
            }
            else
            {
                delta_val = fd_delta;
                val1 =  fd_delta;
            }
            delta_sol.setVal(i, val1);
            
            // now calculate the residual
            fluid_elem.clear();
            fluid_elem.initialize(elem, fe, q_rule, dt, cp, cv, delta_sol, elem_vel, false, elem_dc); // reuse the specified dc vals
            
            evaluateBoundaryConditionData(elem, q_boundary, fluid_elem, sol_inf, momentum_flux_tensor, mass_flux, energy_flux, false, tmp_vec, tmp_mat, delta_res, tmp_mat2);
            fluid_elem.calculateResidualVector(tmp_vec);
            
            delta_res.add(1.0, tmp_vec); delta_res.add(-1.0, res0); delta_res.scale(1.0/delta_val);
            jac.setColumnVals(i, 0, elem_sol.getSize()-1, delta_res);
        }
    }
    

protected:
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems;
    const FESystemUInt n_elem_nodes, n_elem_dofs;
    const FESystemDouble dt;
    const FESystem::Base::DegreeOfFreedomMap& dof_map;
    const FESystem::Quadrature::QuadratureBase& q_rule;
    const FESystem::Quadrature::QuadratureBase& q_boundary;
    const FESystem::Numerics::VectorBase<FESystemDouble>& sol;
    const FESystem::Numerics::VectorBase<FESystemDouble>& vel;
    const FESystemBoolean if_calculate_jacobian;
    const FESystemBoolean if_update_elem_dc_vals;
    std::vector<std::vector<FESystemDouble> >& elem_dc_vals;
    FESystem::Numerics::VectorBase<FESystemDouble>& residual;
    FESystem::Numerics::MatrixBase<FESystemDouble>& global_stiffness_mat;
    FESystem::Numerics::MatrixBase<FESystemDouble>& global_mass_mat;
};




class EvaluateLinearizedForce
{
public:
    EvaluateLinearizedForce(const std::vector<FESystem::Mesh::ElemBase*>& e, FESystemDouble t_step,
                            FESystemUInt n_nd, FESystemUInt n_dof,
                            const FESystem::Base::DegreeOfFreedomMap& dof,
                            const FESystem::Quadrature::QuadratureBase& q,
                            const FESystem::Quadrature::QuadratureBase& q_b,
                            const FESystem::Numerics::VectorBase<FESystemDouble>& s,
                            const FESystem::Numerics::VectorBase<FESystemDouble>& v,
                            const FESystem::Numerics::VectorBase<FESystemDouble>& str_vel,
                            std::vector<std::vector<FESystemDouble> >& dc,
                            FESystem::Numerics::VectorBase<FESystemDouble>& f):
    elems(e),
    n_elem_nodes(n_nd),
    n_elem_dofs(n_dof),
    dt(t_step),
    dof_map(dof),
    q_rule(q),
    q_boundary(q_b),
    sol(s),
    vel(v),
    structural_vel(str_vel),
    elem_dc_vals(dc),
    force_vec(f)
    {
    }
    
    
    void operator() (const tbb::blocked_range<FESystemUInt>& r) const
    {
        this->assemble(r.begin(), r.end());
    }
    
    
    void assembleAll() const
    {
        this->assemble(0, elems.size());
    }
    
    
    void assemble(FESystemUInt e1, FESystemUInt e2) const
    {
        FESystem::Numerics::LocalVector<FESystemDouble> elem_vec, elem_sol, elem_vel, motion_vec;
        elem_vec.resize(n_elem_dofs); elem_sol.resize(n_elem_dofs); elem_vel.resize(n_elem_dofs); motion_vec.resize(dim*n_elem_nodes);
        
        FESystem::FiniteElement::FELagrange fe;
        FESystem::Fluid::FluidElementBase fluid_elem;
        
        motion_vec.zero();
        for (FESystemUInt i=0; i<n_elem_nodes; i++)
            motion_vec.setVal(n_elem_nodes+i, 1.0);
        
        for (FESystemUInt i=e1; i!=e2; i++)
        {
            fe.clear();
            fe.reinit(*(elems[i]));
            elem_vec.zero();
            
            dof_map.getFromGlobalVector(*(elems[i]), sol, elem_sol);
            dof_map.getFromGlobalVector(*(elems[i]), vel, elem_vel);
            
            fluid_elem.clear();
            fluid_elem.initialize(*(elems[i]), fe, q_rule, dt, cp, cv, elem_sol, elem_vel, false, elem_dc_vals[i]);
            
            
            if (elems[i]->checkForTag(4))
            {
                fluid_elem.calculateLinearizedForcingFunctionForBoundaryMotion(0, q_boundary, motion_vec, elem_vec);
                
                {
                    tbb::mutex::scoped_lock my_lock(assembly_mutex);
                    dof_map.addToGlobalVector(*(elems[i]), elem_vec, force_vec);
                }
            }
        }
    }
    
    
protected:
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems;
    const FESystemUInt n_elem_nodes, n_elem_dofs;
    const FESystemDouble dt;
    const FESystem::Base::DegreeOfFreedomMap& dof_map;
    const FESystem::Quadrature::QuadratureBase& q_rule;
    const FESystem::Quadrature::QuadratureBase& q_boundary;
    const FESystem::Numerics::VectorBase<FESystemDouble>& sol;
    const FESystem::Numerics::VectorBase<FESystemDouble>& vel;
    const FESystem::Numerics::VectorBase<FESystemDouble>& structural_vel;
    std::vector<std::vector<FESystemDouble> >& elem_dc_vals;
    FESystem::Numerics::VectorBase<FESystemDouble>& force_vec;
};



class EvaluateComplexCpPerturbation
{
public:
    EvaluateComplexCpPerturbation(const std::vector<FESystem::Mesh::ElemBase*>& e,
                                  const std::vector<FESystem::Mesh::Node*>& nd,
                                  FESystemDouble t_step,
                                  FESystemUInt n_nd, FESystemUInt n_dof,
                                  const FESystem::Base::DegreeOfFreedomMap& dof,
                                  const FESystem::Quadrature::QuadratureBase& q,
                                  const FESystem::Quadrature::QuadratureBase& q_b,
                                  const FESystem::Numerics::VectorBase<FESystemDouble>& s,
                                  const FESystem::Numerics::VectorBase<FESystemDouble>& v,
                                  std::vector<std::vector<FESystemDouble> >& dc,
                                  const FESystem::Numerics::VectorBase<FESystemComplexDouble>& cmplx_sol,
                                  FESystem::Numerics::VectorBase<FESystemDouble>& additional_vals):
    elems(e),
    nodes(nd),
    n_elem_nodes(n_nd),
    n_elem_dofs(n_dof),
    dt(t_step),
    dof_map(dof),
    q_rule(q),
    q_boundary(q_b),
    sol(s),
    vel(v),
    elem_dc_vals(dc),
    complex_sol(cmplx_sol),
    cp_vals(additional_vals)
    {
    }
    
    void operator() (const tbb::blocked_range<FESystemUInt>& r) const
    {
        this->assemble(r.begin(), r.end());
    }
    
    
    void assembleAll() const
    {
        this->assemble(0, elems.size());
    }
    
    void assemble(const FESystemUInt n1, const FESystemUInt n2) const
    {
        FESystem::Numerics::LocalVector<FESystemDouble> elem_sol, elem_vel, tmp_vec;
        FESystem::Numerics::LocalVector<FESystemComplexDouble> tmp_vec_complex;
        elem_sol.resize(n_elem_dofs); elem_vel.resize(n_elem_dofs); tmp_vec.resize(n_vars); tmp_vec_complex.resize(n_vars);
        
        FESystem::FiniteElement::FELagrange fe;
        FESystem::Fluid::FluidElementBase fluid_elem;
        
        fe.clear();
        fe.reinit(*(elems[0]));
        dof_map.getFromGlobalVector(*(elems[0]), sol, elem_sol);
        dof_map.getFromGlobalVector(*(elems[0]), vel, elem_vel);
        fluid_elem.clear();
        fluid_elem.initialize(*(elems[0]), fe, q_rule, dt, cp, cv, elem_sol, elem_vel, false, elem_dc_vals[0]);
        
        FESystemComplexDouble complex_cp;
        
        for (FESystemUInt i=n1; i!=n2; i++)
        {
            for (FESystemUInt j=0; j<n_vars; j++)
                tmp_vec.setVal(j, sol.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]));
            for (FESystemUInt j=0; j<n_vars; j++)
                tmp_vec_complex.setVal(j, complex_sol.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]));
            
            fluid_elem.calculateComplexCpPerturbation(tmp_vec, tmp_vec_complex, q0, complex_cp);
            
            cp_vals.setVal(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0], FESystem::Base::real<FESystemComplexDouble, FESystemDouble>(complex_cp));
            cp_vals.setVal(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0], FESystem::Base::imag<FESystemComplexDouble, FESystemDouble>(complex_cp));
        }
    }
    
    
protected:
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems;
    const std::vector<FESystem::Mesh::Node*>& nodes;
    const FESystemUInt n_elem_nodes, n_elem_dofs;
    const FESystemDouble dt;
    const FESystem::Base::DegreeOfFreedomMap& dof_map;
    const FESystem::Quadrature::QuadratureBase& q_rule;
    const FESystem::Quadrature::QuadratureBase& q_boundary;
    const FESystem::Numerics::VectorBase<FESystemDouble>& sol;
    const FESystem::Numerics::VectorBase<FESystemDouble>& vel;
    std::vector<std::vector<FESystemDouble> >& elem_dc_vals;
    const FESystem::Numerics::VectorBase<FESystemComplexDouble>& complex_sol;
    FESystem::Numerics::VectorBase<FESystemDouble>& cp_vals;
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
                               std::vector<std::vector<FESystemDouble> >& dc,
                               const FESystemBoolean if_update_dc,
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
    if_update_elem_dc_vals(if_update_dc),
    elem_dc_vals(dc),
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
        fluid_elem.initialize(*(elems[0]), fe, q_rule, dt, cp, cv, elem_sol, elem_vel, false, elem_dc_vals[0]);

        FESystemDouble press, entropy, mach, cp;
        
        for (FESystemUInt i=r.begin(); i!=r.end(); i++)
        {
            for (FESystemUInt j=0; j<n_vars; j++)
                tmp_vec.setVal(j, sol.getVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0]));
            fluid_elem.calculatePrimitiveVariableValues(tmp_vec, q0, p, tmp_vec2, press, entropy, mach, cp);
            for (FESystemUInt j=0; j<n_vars; j++)
                primitive.setVal(nodes[i]->getDegreeOfFreedomUnit(j).global_dof_id[0], tmp_vec2.getVal(j));
            additional.setVal(nodes[i]->getDegreeOfFreedomUnit(0).global_dof_id[0], press);
            additional.setVal(nodes[i]->getDegreeOfFreedomUnit(1).global_dof_id[0], entropy);
            additional.setVal(nodes[i]->getDegreeOfFreedomUnit(2).global_dof_id[0], mach);
            additional.setVal(nodes[i]->getDegreeOfFreedomUnit(3).global_dof_id[0], cp);
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
    const FESystemBoolean if_update_elem_dc_vals;
    std::vector<std::vector<FESystemDouble> >& elem_dc_vals;
    FESystem::Numerics::VectorBase<FESystemDouble>& primitive;
    FESystem::Numerics::VectorBase<FESystemDouble>& additional;
};



void calculateEulerQuantities(FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                              const FESystem::Mesh::MeshBase& mesh,
                              const FESystemDouble dt, const FESystemUInt iter_num,
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
            q_rule.init(2, 1);  // bending quadrature is higher than the shear quadrature for reduced integrations
            q_boundary.init(1,1);
            break;
            
        case FESystem::Mesh::QUAD9:
        case FESystem::Mesh::TRI6:
            q_rule.init(2, 4);
            q_boundary.init(1, 4);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    residual.zero(); global_stiffness_mat.zero(); global_mass_mat.zero();
    primitive_sol.zero(); additional_sol.zero();
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
 
    // prepare the vector for dc values
    if (dc_vals.size() != elems.size())
    {
        FESystemUInt nq = q_rule.getQuadraturePoints().size();
        dc_vals.resize(elems.size());
        for (FESystemUInt i=0; i<dc_vals.size(); i++)
            dc_vals[i].resize(nq);
    }
    FESystemBoolean if_update_dc_vals = true;
    if (iter_num > dc_freeze_iter_num) if_update_dc_vals = false;
    
    tbb::parallel_for(tbb::blocked_range<FESystemUInt>(0, elems.size()),
                      AssembleElementMatrices(elems, dt, n_elem_nodes, n_elem_dofs, dof_map, q_rule, q_boundary, sol, vel, dc_vals, if_update_dc_vals, if_calculate_jacobian, residual, global_stiffness_mat, global_mass_mat));
    
    
//    AssembleElementMatrices a(elems, dt, n_elem_nodes, n_elem_dofs, dof_map, q_rule, q_boundary, sol, vel, dc_vals, if_update_dc_vals, if_calculate_jacobian, residual, global_stiffness_mat, global_mass_mat);
//    a.assembleAll();
    
    tbb::parallel_for(tbb::blocked_range<FESystemUInt>(0, nodes.size()),
                      EvaluatePrimitiveVariables(elems, nodes, dt, n_elem_dofs, dof_map, q_rule, sol, vel, dc_vals,if_update_dc_vals, primitive_sol, additional_sol));
}



void transientEulerAnalysis(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                            FESystem::Numerics::VectorBase<FESystemDouble>& final_sol)
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
    transient_solver.enableAdaptiveTimeStepping(4, 0.4, 1.0e3);
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
    FESystemBoolean if_converged = false;
    
    while (!if_converged)
    {
        call_back = transient_solver.incrementTimeStep();
        
        if (transient_solver.getCurrentTime()>final_t)
            if_converged = true;
        
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
                
                if (transient_solver.getCurrentStateVelocityVector().getLInfNorm() <= 1.0e-10)
                    if_converged = true;
            }
                break;
                
            case FESystem::TransientSolvers::EVALUATE_X_DOT:
                // the Jacobian is not updated since it is constant with respect to time
            {
                //testJacobian(elem_type, n_elem_nodes, dof_map, mesh, sol, vel, vel_func);
                calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, transient_solver.getCurrentStepSize(), transient_solver.getCurrentIterationNumber(),
                                         transient_solver.getCurrentStateVector(), transient_solver.getCurrentStateVelocityVector(), false,
                                         transient_solver.getCurrentVelocityFunctionVector(), transient_solver.getCurrentJacobianMatrix(),
                                         mass, primitive_sol, additional_sol);
            }
                break;

            case FESystem::TransientSolvers::EVALUATE_X_DOT_AND_X_DOT_JACOBIAN:
                // the Jacobian is not updated since it is constant with respect to time
            {
                //testJacobian(elem_type, n_elem_nodes, dof_map, mesh, sol, vel, vel_func);
                calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, transient_solver.getCurrentStepSize(), transient_solver.getCurrentIterationNumber(),
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
    
    // return the final vector to the calling funciton
    final_sol.copyVector(transient_solver.getCurrentStateVector());
}




void frequencyDomainAnalysis(FESystemUInt dim, FESystem::Mesh::ElementType elem_type, FESystemUInt n_elem_nodes, const FESystem::Mesh::MeshBase& mesh, const FESystem::Base::DegreeOfFreedomMap& dof_map,
                             const FESystem::Numerics::VectorBase<FESystemDouble>& final_sol)
{
    FESystem::LinearSolvers::LUFactorizationLinearSolver<FESystemComplexDouble> linear_solver;
    const std::vector<FESystem::Mesh::Node*>& nodes = mesh.getNodes();
    
    FESystem::Numerics::SparseMatrix<FESystemDouble> stiff_mat, mass;
    FESystem::Numerics::LocalVector<FESystemDouble>  vel, f_vec, res, dummy_vec, additional_sol;
    FESystem::Numerics::SparseMatrix<FESystemComplexDouble> complex_mat, tmp_complex_mat;
    FESystem::Numerics::LocalVector<FESystemComplexDouble> complex_f_vec, complex_sol;
    
    vel.resize(dof_map.getNDofs()); f_vec.resize(dof_map.getNDofs()); res.resize(dof_map.getNDofs()); additional_sol.resize(dof_map.getNDofs());
    stiff_mat.resize(dof_map.getSparsityPattern());  mass.resize(dof_map.getSparsityPattern());
    
    complex_f_vec.resize(dof_map.getNDofs()); complex_sol.resize(dof_map.getNDofs());
    complex_mat.resize(dof_map.getSparsityPattern()); tmp_complex_mat.resize(dof_map.getSparsityPattern());
    

    
    FESystemUInt n_elem_dofs;
    
    n_elem_dofs = n_vars*n_elem_nodes;
    
    // prepare the quadrature rule and FE for the
    FESystem::Quadrature::TrapezoidQuadrature q_rule, q_boundary;
    
    switch (elem_type)
    {
        case FESystem::Mesh::QUAD4:
        case FESystem::Mesh::TRI3:
            q_rule.init(2, 1);  // bending quadrature is higher than the shear quadrature for reduced integrations
            q_boundary.init(1,1);
            break;
            
        case FESystem::Mesh::QUAD9:
        case FESystem::Mesh::TRI6:
            q_rule.init(2, 4);
            q_boundary.init(1, 4);
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::EnumNotHandled);
    }
    
    additional_sol.zero();
    
    const std::vector<FESystem::Mesh::ElemBase*>& elems = mesh.getElements();
    
    // calculate the Jacobian and mass matrices
    calculateEulerQuantities(elem_type, n_elem_nodes, dof_map, mesh, 1.0e3, 1000, final_sol, vel, true, res, stiff_mat, mass, additional_sol, additional_sol);
    
    // calculate the force vector
    tbb::parallel_for(tbb::blocked_range<FESystemUInt>(0, elems.size()),
                      EvaluateLinearizedForce(elems, 1.0e3, n_elem_nodes, n_elem_dofs, dof_map, q_rule, q_boundary, final_sol, vel, dummy_vec, dc_vals, f_vec));
    
    FESystemDouble omega = 50.0;
    complex_mat.copyRealMatrix(mass); complex_mat.scale(FESystemComplexDouble(0.0,1.0)*omega);
    tmp_complex_mat.copyRealMatrix(stiff_mat);
    complex_mat.add(FESystemComplexDouble(-1.0,0.0), tmp_complex_mat);
    
    complex_f_vec.copyRealVector(f_vec); complex_f_vec.scale(FESystemComplexDouble(0.0, 1.0)*omega);
    complex_sol.zero();
    
    linear_solver.setSystemMatrix(complex_mat, false);
    linear_solver.solve(complex_f_vec, complex_sol);

    // calculate the complex cp values
    additional_sol.zero();
    tbb::parallel_for(tbb::blocked_range<FESystemUInt>(0, nodes.size()),
                      EvaluateComplexCpPerturbation(elems, nodes, 1.0e3, n_elem_nodes, n_elem_dofs, dof_map, q_rule, q_boundary, final_sol, vel, dc_vals, complex_sol, additional_sol));
    
    FESystem::OutputProcessor::VtkOutputProcessor output;
    std::fstream output_file;
    std::vector<FESystemUInt> vars(4); vars[0]=0; vars[1]=1; vars[2]=2; vars[3] = 3; //vars[4] = 4; // write all solutions
    output_file.open("perturb_sol.vtk",std::fstream::out);
    output.writeMesh(output_file, mesh, dof_map);
    output.writeSolution(output_file, "Sol", mesh, dof_map, vars, final_sol);
    output.writeSolution(output_file, "Perturbed", mesh, dof_map, vars, additional_sol);
    output_file.close();
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
    
    initMeshParameters();
    
    distributePoints(x_n_divs, x_div_locations, x_n_subdivs_in_div, x_relative_mesh_size_in_div, x_points);
    distributePoints(y_n_divs, y_div_locations, y_n_subdivs_in_div, y_relative_mesh_size_in_div, y_points);
    
    createPlaneMesh(elem_type, mesh, origin, x_points, y_points, n_elem_nodes, CROSS, true);
    
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
    
    // modify mesh after application of boundary condition
    modifyMeshForCase(mesh);
    
    FESystem::Numerics::LocalVector<FESystemDouble> sol_vec;
    sol_vec.resize(dof_map.getNDofs());
    
    // transient solution
    transientEulerAnalysis(dim, elem_type, n_elem_nodes, mesh, dof_map, sol_vec);

    // frequency domain solution
    frequencyDomainAnalysis(dim, elem_type, n_elem_nodes, mesh, dof_map, sol_vec);
    
    //exit(1);
    
    return 0;
}





