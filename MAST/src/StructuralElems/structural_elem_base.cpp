//
//  structural_elememt_base.cpp
//  MAST
//
//  Created by Manav Bhatia on 10/16/13.
//
//

// MAST includes
#include "StructuralElems/structural_elem_base.h"
#include "StructuralElems/strain_operator.h"
#include "PropertyCards/element_property_card_base.h"

// libMesh includes
#include "libmesh/quadrature.h"


bool
MAST::StructuralElementBase::internal_force (bool request_jacobian,
                                       DenseVector<Real>& f,
                                       DenseMatrix<Real>& jac)
{
    std::auto_ptr<MAST::StrainOperator>
    strain_operator(MAST::build_strain_operator().release());
 
    const std::vector<Real>& JxW = strain_operator->get_JxW();
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseVector<Real>  tmp_vec1_n1, tmp_vec2_n2;
    
    _property->calculate_matrix(*_elem,
                                MAST::SECTION_INTEGRATED_STRESS_MATERIAL_MATRIX,
                                material_mat);
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        strain_operator->initialize_for_qp(qp);
        strain_operator->left_multiply(tmp_mat1_n1n2, material_mat);
        
        tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, _local_solution);
        strain_operator->vector_mult_transpose(tmp_vec2_n2, tmp_vec1_n1);
        f.add(JxW[qp], tmp_vec2_n2);
        
        if (request_jacobian) {

            strain_operator->right_multiply_transpose(tmp_mat2_n2n2,
                                                      tmp_mat1_n1n2);
            jac.add(JxW[qp], tmp_mat2_n2n2);
        }
        
    }
    
    return request_jacobian;
}



bool
MAST::StructuralElementBase::damping_force (bool request_jacobian,
                                            DenseVector<Real>& f,
                                            DenseMatrix<Real>& jac)
{
    FEMOperatorMatrix Bmat;
    FEBase* fe = NULL;
    QBase* qrule = NULL;
    get_fe_and_qrule(fe, qrule);
    
    const std::vector<Real>& JxW = fe->get_JxW();
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseVector<Real>  phi, tmp_vec1_n1, tmp_vec2_n2;
    
    _property->calculate_matrix(*_elem,
                                MAST::SECTION_INTEGRATED_DAMPING_MATERIAL_MATRIX,
                                material_mat);
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        Bmat.reinit(_system->n_vars(), phi);
        
        Bmat.left_multiply(tmp_mat1_n1n2, material_mat);
        
        tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, _local_velocity);
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
    get_fe_and_qrule(fe, qrule);
    
    const std::vector<Real>& JxW = fe->get_JxW();
    DenseMatrix<Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2;
    DenseVector<Real>  phi, tmp_vec1_n1, tmp_vec2_n2;
    
    _property->calculate_matrix(*_elem,
                                MAST::SECTION_INTEGRATED_INERTIA_MATERIAL_MATRIX,
                                material_mat);
    
    for (unsigned int qp=0; qp<JxW.size(); qp++) {
        
        Bmat.reinit(_system->n_vars(), phi);
        
        Bmat.left_multiply(tmp_mat1_n1n2, material_mat);
        
        tmp_mat1_n1n2.vector_mult(tmp_vec1_n1, _local_acceleration);
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
MAST::StructuralElementBase::get_fe_and_qrule(FEBase* fe, QBase* qrule) {

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

