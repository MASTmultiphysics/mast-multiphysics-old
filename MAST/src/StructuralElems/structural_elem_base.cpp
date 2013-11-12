//
//  structural_elememt_base.cpp
//  MAST
//
//  Created by Manav Bhatia on 10/16/13.
//
//

// libMesh includes
#include "libmesh/quadrature.h"

// MAST includes
#include "StructuralElems/structural_elem_base.h"
#include "PropertyCards/element_property_card_base.h"
#include "Numerics/fem_operator_matrix.h"
#include "StructuralElems/structural_element_3D.h"
#include "StructuralElems/structural_element_2D.h"
//#include "StructuralElems/structural_element_1D.h"
#include "ThermalElems/temperature_function.h"
#include "Numerics/sensitivity_parameters.h"
#include "Base/boundary_condition.h"


MAST::StructuralElementBase::StructuralElementBase(System& sys,
                                                   const Elem& elem,
                                                   const MAST::ElementPropertyCardBase& p):
_system(sys),
_elem(elem),
_property(p),
_temperature(NULL)
{ }



MAST::StructuralElementBase::~StructuralElementBase()
{ }




bool
MAST::StructuralElementBase::damping_force (bool request_jacobian,
                                            DenseVector<Real>& f,
                                            DenseMatrix<Real>& jac)
{
    libmesh_error(); // to be implemented
    
}



bool
MAST::StructuralElementBase::damping_force_sensitivity(bool request_jacobian,
                                                       DenseVector<Real>& f,
                                                       DenseMatrix<Real>& jac)
{
    libmesh_error(); // to be implemented
}



bool
MAST::StructuralElementBase::inertial_force (bool request_jacobian,
                                             DenseVector<Real>& f,
                                             DenseMatrix<Real>& jac)
{
    FEMOperatorMatrix Bmat;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<std::vector<Real> >& phi = _fe->get_phi();
    const unsigned int n_phi = (unsigned int)phi.size(), n1=6, n2=6*n_phi;
    
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2, local_jac;
    DenseVector<Real>  phi_vec, tmp_vec1_n1, tmp_vec2_n2, local_f;
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2); local_jac.resize(n2, n2);
    phi_vec.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n2);
    local_f.resize(n2);
    
    _property.calculate_matrix(_elem,
                               MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX,
                               material_mat);
    
    if (_property.if_diagonal_mass_matrix()) {
        Real vol = 0.;
        const unsigned int nshp = _fe->n_shape_functions();
        for (unsigned int i=0; i<JxW.size(); i++)
            vol += JxW[i];
        vol /= (1.* nshp);
        for (unsigned int i_var=0; i_var<6; i_var++)
            for (unsigned int i=0; i<nshp; i++)
                local_jac(i_var*nshp+i, i_var*nshp+i) =
                vol*material_mat(i_var, i_var);
        local_jac.vector_mult(local_f, local_acceleration);
    }
    else {
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            
            // now set the shape function values
            for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
                phi_vec(i_nd) = phi[i_nd][qp];
            
            Bmat.reinit(_system.n_vars(), phi_vec);
            
            Bmat.left_multiply(tmp_mat1_n1n2, material_mat);
            
            tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, local_acceleration);
            Bmat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
            local_f.add(JxW[qp], tmp_vec2_n2);
            
            if (request_jacobian) {
                
                Bmat.right_multiply_transpose(tmp_mat2_n2n2,
                                              tmp_mat1_n1n2);
                local_jac.add(JxW[qp], tmp_mat2_n2n2);
            }
            
        }
    }
    
    // now transform to the global coorodinate system
    if (_elem.dim() < 3) {
        _transform_to_global_system(local_f, tmp_vec2_n2);
        f.add(1., tmp_vec2_n2);
        if (request_jacobian) {
            _transform_to_global_system(local_jac, tmp_mat2_n2n2);
            jac.add(1., tmp_mat2_n2n2);
        }
    }
    else {
        f.add(1., local_f);
        if (request_jacobian)
            jac.add(1., local_jac);
    }

    return request_jacobian;
}



bool
MAST::StructuralElementBase::inertial_force_sensitivity(bool request_jacobian,
                                                        DenseVector<Real>& f,
                                                        DenseMatrix<Real>& jac)
{
    // this should be true if the function is called
    libmesh_assert(this->sensitivity_params);
    libmesh_assert(!this->sensitivity_params->shape_sensitivity()); // this is not implemented for now
    libmesh_assert(this->sensitivity_params->total_order() == 1); // only first order sensitivity
    
    // check if the material property or the provided exterior
    // values, like temperature, are functions of the sensitivity parameter
    bool calculate = false;
    if (_temperature)
        calculate = calculate || _temperature->depends_on(*(this->sensitivity_params));
    calculate = calculate || _property.depends_on(*(this->sensitivity_params));
    
    // nothing to be calculated if the element does not depend on the
    // sensitivity parameter.
    if (!calculate)
        return false;
    
    FEMOperatorMatrix Bmat;
    
    const std::vector<Real>& JxW = _fe->get_JxW();
    const std::vector<std::vector<Real> >& phi = _fe->get_phi();
    const unsigned int n_phi = (unsigned int)phi.size(), n1=6, n2=6*n_phi;
    
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2, local_jac;
    DenseVector<Real>  phi_vec, tmp_vec1_n1, tmp_vec2_n2, local_f;
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2); local_jac.resize(n2, n2);
    phi_vec.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n2);
    local_f.resize(n2);
    
    _property.calculate_matrix_sensitivity(_elem,
                                           MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX,
                                           material_mat,
                                           *(this->sensitivity_params));
    
    if (_property.if_diagonal_mass_matrix()) {
        Real vol = 0.;
        const unsigned int nshp = _fe->n_shape_functions();
        for (unsigned int i=0; i<JxW.size(); i++)
            vol += JxW[i];
        vol /= (1.* nshp);
        for (unsigned int i_var=0; i_var<6; i_var++)
            for (unsigned int i=0; i<nshp; i++)
                local_jac(i_var*nshp+i, i_var*nshp+i) =
                vol*material_mat(i_var, i_var);
        local_jac.vector_mult(local_f, local_acceleration);
    }
    else {
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            
            // now set the shape function values
            for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
                phi_vec(i_nd) = phi[i_nd][qp];
            
            Bmat.reinit(_system.n_vars(), phi_vec);
            
            Bmat.left_multiply(tmp_mat1_n1n2, material_mat);
            
            tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, local_acceleration);
            Bmat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
            local_f.add(JxW[qp], tmp_vec2_n2);
            
            if (request_jacobian) {
                
                Bmat.right_multiply_transpose(tmp_mat2_n2n2,
                                              tmp_mat1_n1n2);
                local_jac.add(JxW[qp], tmp_mat2_n2n2);
            }
            
        }
    }
    
    // now transform to the global coorodinate system
    if (_elem.dim() < 3) {
        _transform_to_global_system(local_f, tmp_vec2_n2);
        f.add(1., tmp_vec2_n2);
        if (request_jacobian) {
            _transform_to_global_system(local_jac, tmp_mat2_n2n2);
            jac.add(1., tmp_mat2_n2n2);
        }
    }
    else {
        f.add(1., local_f);
        if (request_jacobian)
            jac.add(1., local_jac);
    }
    
    return request_jacobian;
}




bool
MAST::StructuralElementBase::side_external_force(bool request_jacobian,
                                                 DenseVector<Real> &f,
                                                 DenseMatrix<Real> &jac) {
    
    // iterate over the boundary ids given in the provided force map
    std::multimap<unsigned int, MAST::BoundaryCondition*>::const_iterator
    bc_it1, bc_it2;
    
    const BoundaryInfo& binfo = *_system.get_mesh().boundary_info;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    boundary_id_type bc_id;
    bool calculate_jac = false;
    
    for (unsigned short int n=0; n<_elem.n_sides(); n++) {
        bc_id = binfo.boundary_id(&_elem, n);
        if ((bc_id != BoundaryInfo::invalid_id) &&
            _side_bc_map.count(bc_id)) {
            // find the loads on this boundary and evaluate the f and jac
            switch (bc_it1->second->type()) {
                case MAST::SURFACE_PRESSURE:
                    calculate_jac = (calculate_jac ||
                                     surface_pressure_force(request_jacobian,
                                                            f, jac,
                                                            n,
                                                            *bc_it1->second));
                    break;
                    
                default:
                    // not implemented yet
                    libmesh_error();
                    break;
            }
        }
    }
    
    return (request_jacobian && calculate_jac);
}



bool
MAST::StructuralElementBase::side_external_force_sensitivity(bool request_jacobian,
                                                             DenseVector<Real> &f,
                                                             DenseMatrix<Real> &jac) {
    
    libmesh_assert(false);
}



bool
MAST::StructuralElementBase::volume_external_force(bool request_jacobian,
                                                   DenseVector<Real> &f,
                                                   DenseMatrix<Real> &jac) {
    
    libmesh_error();
    
    bool calculate_jac = false; // start with false, and then set true if
                                // any one of the forces provides a Jacobian
    
    if (_temperature)  // thermal load
        calculate_jac = (calculate_jac ||
                         this->thermal_force(request_jacobian, f, jac));
    
    return (request_jacobian && calculate_jac);
}




bool
MAST::StructuralElementBase::volume_external_force_sensitivity(bool request_jacobian,
                                                               DenseVector<Real> &f,
                                                               DenseMatrix<Real> &jac) {
    
    libmesh_error();
    
    bool calculate_jac = false; // start with false, and then set true if
                                // any one of the forces provides a Jacobian
    
    if (_temperature)  // thermal load
        calculate_jac = (calculate_jac ||
                         this->thermal_force(request_jacobian, f, jac));
    
    return (request_jacobian && calculate_jac);
}




void
MAST::StructuralElementBase::_transform_to_global_system(const DenseMatrix<Real>& local_mat,
                                                        DenseMatrix<Real>& global_mat) const {
    libmesh_assert_equal_to( local_mat.m(),  local_mat.n());
    libmesh_assert_equal_to(global_mat.m(), global_mat.n());
    libmesh_assert_equal_to( local_mat.m(), global_mat.m());
    
    const unsigned int n_dofs = _fe->n_shape_functions();
    global_mat.zero();
    DenseMatrix<Real> tmp_mat;
    tmp_mat.resize(6*n_dofs, 6*n_dofs);
    
    const DenseMatrix<Real>& Tmat = _transformation_matrix();
    // now initialize the global T matrix
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++) {
                tmp_mat(j*n_dofs+i, k*n_dofs+i) = Tmat(j,k); // for u,v,w
                tmp_mat((j+3)*n_dofs+i, (k+3)*n_dofs+i) = Tmat(j,k); // for tx,ty,tz
            }
    
    // right multiply with T^T, and left multiply with T.
    global_mat = local_mat;
    global_mat.right_multiply_transpose(tmp_mat);
    global_mat.left_multiply(tmp_mat);
}



void
MAST::StructuralElementBase::_transform_to_local_system(const DenseVector<Real>& global_vec,
                                                       DenseVector<Real>& local_vec) const {
    libmesh_assert_equal_to( local_vec.size(),  global_vec.size());
    
    const unsigned int n_dofs = _fe->n_shape_functions();
    local_vec.zero();
    DenseMatrix<Real> tmp_mat;
    tmp_mat.resize(6*n_dofs, 6*n_dofs);
    
    const DenseMatrix<Real>& Tmat = _transformation_matrix();
    // now initialize the global T matrix
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++) {
                tmp_mat(j*n_dofs+i, k*n_dofs+i) = Tmat(j,k); // for u,v,w
                tmp_mat((j+3)*n_dofs+i, (k+3)*n_dofs+i) = Tmat(j,k); // for tx,ty,tz
            }
    
    // left multiply with T^T
    tmp_mat.vector_mult_transpose(local_vec, global_vec);
}



void
MAST::StructuralElementBase::_transform_to_global_system(const DenseVector<Real>& local_vec,
                                                        DenseVector<Real>& global_vec) const {
    libmesh_assert_equal_to( local_vec.size(),  global_vec.size());
    
    const unsigned int n_dofs = _fe->n_shape_functions();
    global_vec.zero();
    DenseMatrix<Real> tmp_mat;
    tmp_mat.resize(6*n_dofs, 6*n_dofs);
    
    const DenseMatrix<Real>& Tmat = _transformation_matrix();
    // now initialize the global T matrix
    for (unsigned int i=0; i<n_dofs; i++)
        for (unsigned int j=0; j<3; j++)
            for (unsigned int k=0; k<3; k++) {
                tmp_mat(j*n_dofs+i, k*n_dofs+i) = Tmat(j,k); // for u,v,w
                tmp_mat((j+3)*n_dofs+i, (k+3)*n_dofs+i) = Tmat(j,k); // for tx,ty,tz
            }

    // left multiply with T
    tmp_mat.vector_mult(global_vec, local_vec);
}





void
MAST::StructuralElementBase::_init_fe_and_qrule( const Elem& e) {
    
    unsigned int nv = _system.n_vars();

    libmesh_assert (nv);
    FEType fe_type = _system.variable_type(0); // all variables are assumed to be of same type
    

    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _system.variable_type(i));

    // Create an adequate quadrature rule
    _fe.reset(FEBase::build(e.dim(), fe_type).release());
    _qrule.reset(fe_type.default_quadrature_rule
                 (e.dim(),
                  _system.extra_quadrature_order +  // system extra quadrature
                  _property.extra_quadrature_order(e, fe_type)).release()); // elem extra quadrature
    _fe->attach_quadrature_rule(_qrule.get());
    _fe->get_phi();
    _fe->get_JxW();
    _fe->get_dphi();
    
    _fe->reinit(&e);
}


void
MAST::StructuralElementBase::_get_side_fe_and_qrule(unsigned int s,
                                                    std::auto_ptr<FEBase>& fe,
                                                    std::auto_ptr<QBase>& qrule) {
//    unsigned int nv = _system.n_vars();
//    
//    libmesh_assert (nv);
//    FEType fe_type = _system.variable_type(0); // all variables are assumed to be of same type
//    
//    
//    for (unsigned int i=1; i != nv; ++i)
//        libmesh_assert(fe_type == _system.variable_type(i));
//    
//    // Create an adequate quadrature rule
//    qrule.reset(fe_type.default_quadrature_rule
//                (_elem.dim(), _system.extra_quadrature_order).release());
//    fe.reset(FEBase::build(_elem.dim(), fe_type).release());
//    
//    side_qrule = hardest_fe_type.default_quadrature_rule
//    (dim-1, sys.extra_quadrature_order).release();
//    if (dim == 3)
//        edge_qrule = hardest_fe_type.default_quadrature_rule
//        (1, _system.extra_quadrature_order).release();
//    
//    // Next, create finite element objects
//    // Preserving backward compatibility here for now
//    // Should move to the protected/FEAbstract interface
//    _element_fe_var.resize(nv);
//    _side_fe_var.resize(nv);
//    if (dim == 3)
//        _edge_fe_var.resize(nv);
//    
//    for (unsigned int i=0; i != nv; ++i)
//    {
//        FEType fe_type = sys.variable_type(i);
//        
//        if ( _element_fe[fe_type] == NULL )
//        {
//            _element_fe[fe_type] = FEAbstract::build(dim, fe_type).release();
//            _element_fe[fe_type]->attach_quadrature_rule(element_qrule);
//            _side_fe[fe_type] = FEAbstract::build(dim, fe_type).release();
//            _side_fe[fe_type]->attach_quadrature_rule(side_qrule);
//            
//            if (dim == 3)
//            {
//                _edge_fe[fe_type] = FEAbstract::build(dim, fe_type).release();
//                _edge_fe[fe_type]->attach_quadrature_rule(edge_qrule);
//            }
//        }
//        _element_fe_var[i] = _element_fe[fe_type];
//        _side_fe_var[i] = _side_fe[fe_type];
//        if (dim == 3)
//            _edge_fe_var[i] = _edge_fe[fe_type];
//        
//    }
}



std::auto_ptr<MAST::StructuralElementBase>
MAST::build_structural_element(System& sys,
                               const Elem& elem,
                               const MAST::ElementPropertyCardBase& p) {
    
    std::auto_ptr<MAST::StructuralElementBase> e;
    
    switch (elem.dim()) {
        case 3:
            e.reset(new MAST::StructuralElement3D(sys, elem, p));
            break;
            
        case 2:
            e.reset(new MAST::StructuralElement2D(sys, elem, p));
            break;

        case 1:
        default:
            libmesh_error();
            break;
    }
    
    return e;
}


