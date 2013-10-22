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


MAST::StructuralElementBase::StructuralElementBase(System& sys,
                                                   const Elem& elem,
                                                   const MAST::ElementPropertyCardBase& p):
_system(sys),
_elem(elem),
_property(p)
{ }



MAST::StructuralElementBase::~StructuralElementBase()
{ }




bool
MAST::StructuralElementBase::damping_force (bool request_jacobian,
                                            DenseVector<Real>& f,
                                            DenseMatrix<Real>& jac)
{
    FEMOperatorMatrix Bmat;
    FEBase* fe = NULL;
    QBase* qrule = NULL;
    _get_fe_and_qrule(fe, qrule);
    
    const std::vector<Real>& JxW = fe->get_JxW();
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseVector<Real>  phi, tmp_vec1_n1, tmp_vec2_n2;
    
    _property.calculate_matrix(_elem,
                               MAST::SECTION_INTEGRATED_MATERIAL_DAMPING_MATRIX,
                               material_mat);
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        Bmat.reinit(_system.n_vars(), phi);
        
        Bmat.left_multiply(tmp_mat1_n1n2, material_mat);
        
        tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, local_velocity);
        Bmat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        f.add(JxW[qp], tmp_vec2_n2);
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(tmp_mat2_n2n2,
                                          tmp_mat1_n1n2);
            jac.add(JxW[qp], tmp_mat2_n2n2);
        }
        
    }
    
    delete fe;
    delete qrule;
    
    return request_jacobian;
}



bool
MAST::StructuralElementBase::inertial_force (bool request_jacobian,
                                             DenseVector<Real>& f,
                                             DenseMatrix<Real>& jac)
{
    FEMOperatorMatrix Bmat;
    FEBase* fe = NULL;
    QBase* qrule = NULL;
    _get_fe_and_qrule(fe, qrule);
    
    const std::vector<Real>& JxW = fe->get_JxW();
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseVector<Real>  phi, tmp_vec1_n1, tmp_vec2_n2;
    
    _property.calculate_matrix(_elem,
                               MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX,
                               material_mat);
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        Bmat.reinit(_system.n_vars(), phi);
        
        Bmat.left_multiply(tmp_mat1_n1n2, material_mat);
        
        tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, local_acceleration);
        Bmat.vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        f.add(JxW[qp], tmp_vec2_n2);
        
        if (request_jacobian) {
            
            Bmat.right_multiply_transpose(tmp_mat2_n2n2,
                                          tmp_mat1_n1n2);
            jac.add(JxW[qp], tmp_mat2_n2n2);
        }
        
    }
    
    delete fe;
    delete qrule;
    
    return request_jacobian;
}



bool
MAST::StructuralElementBase::side_external_force(bool request_jacobian,
                                                 DenseVector<Real> &f,
                                                 DenseMatrix<Real> &jac) {
    
    libmesh_assert(false);
}



bool
MAST::StructuralElementBase::volume_external_force(bool request_jacobian,
                                                   DenseVector<Real> &f,
                                                   DenseMatrix<Real> &jac) {
    
    libmesh_assert(false);
}


void
MAST::StructuralElementBase::_transform_to_global_system(const DenseMatrix<Real>& local_mat,
                                                        DenseMatrix<Real>& global_mat) {
    
}



void
MAST::StructuralElementBase::_transform_to_local_system(const DenseVector<Real>& global_vec,
                                                       DenseVector<Real>& local_vec) {
    
}



void
MAST::StructuralElementBase::_transform_to_global_system(const DenseVector<Real>& local_vec,
                                                        DenseVector<Real>& global_vec) {
    
}




void
MAST::StructuralElementBase::_get_fe_and_qrule(FEBase* fe, QBase* qrule) {
    
    //    // We need to know which of our variables has the hardest
    //    // shape functions to numerically integrate.
    //
    //    unsigned int nv = sys.n_vars();
    //
    //    libmesh_assert (nv);
    //    FEType hardest_fe_type = sys.variable_type(0);
    //
    //    for (unsigned int i=0; i != nv; ++i)
    //    {
    //        FEType fe_type = sys.variable_type(i);
    //
    //        // FIXME - we don't yet handle mixed finite elements from
    //        // different families which require different quadrature rules
    //        // libmesh_assert_equal_to (fe_type.family, hardest_fe_type.family);
    //
    //        if (fe_type.order > hardest_fe_type.order)
    //            hardest_fe_type = fe_type;
    //    }
    //
    //    // Create an adequate quadrature rule
    //    element_qrule = hardest_fe_type.default_quadrature_rule
    //    (dim, sys.extra_quadrature_order).release();
    //    side_qrule = hardest_fe_type.default_quadrature_rule
    //    (dim-1, sys.extra_quadrature_order).release();
    //    if (dim == 3)
    //        edge_qrule = hardest_fe_type.default_quadrature_rule
    //        (1, sys.extra_quadrature_order).release();
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


void
MAST::StructuralElementBase::_get_side_fe_and_qrule(unsigned int s,
                                                    FEBase* fe, QBase* qrule) {
    
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
            
        case 1:
        case 2:
        default:
            libmesh_error();
            break;
    }
    
    return e;
}


