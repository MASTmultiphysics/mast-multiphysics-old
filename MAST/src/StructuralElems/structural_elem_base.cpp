//
//  structural_elememt_base.cpp
//  MAST
//
//  Created by Manav Bhatia on 10/16/13.
//
//

#include <iomanip>
// libMesh includes
#include "libmesh/quadrature.h"

// MAST includes
#include "StructuralElems/structural_elem_base.h"
#include "PropertyCards/element_property_card_base.h"
#include "Numerics/fem_operator_matrix.h"
#include "Numerics/utility.h"
#include "StructuralElems/structural_element_3D.h"
#include "StructuralElems/structural_element_2D.h"
#include "StructuralElems/structural_element_1D.h"
#include "BoundaryConditions/boundary_condition.h"
#include "BoundaryConditions/small_disturbance_motion.h"


MAST::StructuralElementBase::StructuralElementBase(libMesh::System& sys,
                                                   const libMesh::Elem& elem,
                                                   const MAST::ElementPropertyCardBase& p):
sensitivity_param(nullptr),
follower_forces(false),
_system(sys),
_elem(elem),
_property(p)
{ }



MAST::StructuralElementBase::~StructuralElementBase()
{ }




bool
MAST::StructuralElementBase::damping_force (bool request_jacobian,
                                            libMesh::DenseVector<libMesh::Real>& f,
                                            libMesh::DenseMatrix<libMesh::Real>& jac)
{
    libmesh_error(); // to be implemented
    
}



bool
MAST::StructuralElementBase::damping_force_sensitivity(bool request_jacobian,
                                                       libMesh::DenseVector<libMesh::Real>& f,
                                                       libMesh::DenseMatrix<libMesh::Real>& jac)
{
    libmesh_error(); // to be implemented
}



bool
MAST::StructuralElementBase::inertial_force (bool request_jacobian,
                                             libMesh::DenseVector<libMesh::Real>& f,
                                             libMesh::DenseMatrix<libMesh::Real>& jac)
{
    FEMOperatorMatrix Bmat;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const std::vector<std::vector<libMesh::Real> >& phi = _fe->get_phi();
    const unsigned int n_phi = (unsigned int)phi.size(), n1=6, n2=6*n_phi;
    
    libMesh::DenseMatrix<libMesh::Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2, local_jac;
    libMesh::DenseVector<libMesh::Real>  phi_vec, tmp_vec1_n1, tmp_vec2_n2, local_f;
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2); local_jac.resize(n2, n2);
    phi_vec.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n2);
    local_f.resize(n2);
    
    std::auto_ptr<MAST::FieldFunction<libMesh::DenseMatrix<libMesh::Real> > > mat_inertia
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX,
                            *this).release());

    if (_property.if_diagonal_mass_matrix()) {
        // as an approximation, get matrix at the first quadrature point
        (*mat_inertia)(xyz[0], _system.time, material_mat);
        
        libMesh::Real vol = 0.;
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
        libMesh::Point p;
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            
            this->global_coordinates(xyz[qp], p);
            
            (*mat_inertia)(p, _system.time, material_mat);
            
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
        transform_to_global_system(local_f, tmp_vec2_n2);
        f.add(1., tmp_vec2_n2);
        if (request_jacobian) {
            transform_to_global_system(local_jac, tmp_mat2_n2n2);
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
                                                        libMesh::DenseVector<libMesh::Real>& f,
                                                        libMesh::DenseMatrix<libMesh::Real>& jac)
{
    // this should be true if the function is called
    libmesh_assert(this->sensitivity_param);
    libmesh_assert(!this->sensitivity_param->is_shape_parameter()); // this is not implemented for now
    
    
    // check if the material property or the provided exterior
    // values, like temperature, are functions of the sensitivity parameter
    bool calculate = false;
    calculate = calculate || _property.depends_on(*(this->sensitivity_param));
    
    // nothing to be calculated if the element does not depend on the
    // sensitivity parameter.
    if (!calculate)
        return false;
    
    FEMOperatorMatrix Bmat;
    
    const std::vector<libMesh::Real>& JxW = _fe->get_JxW();
    const std::vector<Point>& xyz = _fe->get_xyz();
    const std::vector<std::vector<libMesh::Real> >& phi = _fe->get_phi();
    const unsigned int n_phi = (unsigned int)phi.size(), n1=6, n2=6*n_phi;
    
    libMesh::DenseMatrix<libMesh::Real> material_mat, tmp_mat1_n1n2, tmp_mat2_n2n2, local_jac;
    libMesh::DenseVector<libMesh::Real>  phi_vec, tmp_vec1_n1, tmp_vec2_n2, local_f;
    tmp_mat1_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2); local_jac.resize(n2, n2);
    phi_vec.resize(n_phi); tmp_vec1_n1.resize(n1); tmp_vec2_n2.resize(n2);
    local_f.resize(n2);
    
    std::auto_ptr<MAST::FieldFunction<libMesh::DenseMatrix<libMesh::Real> > > mat_inertia
    (_property.get_property(MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX,
                            *this).release());
    
    if (_property.if_diagonal_mass_matrix()) {
        
        mat_inertia->total(*this->sensitivity_param,
                           xyz[0], _system.time, material_mat);
        
        libMesh::Real vol = 0.;
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
        libMesh::Point p;
        
        for (unsigned int qp=0; qp<JxW.size(); qp++) {
            
            this->global_coordinates(xyz[qp], p);
            
            mat_inertia->total(*this->sensitivity_param,
                               p, _system.time, material_mat);

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
        transform_to_global_system(local_f, tmp_vec2_n2);
        f.add(1., tmp_vec2_n2);
        if (request_jacobian) {
            transform_to_global_system(local_jac, tmp_mat2_n2n2);
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




template <typename ValType>
bool
MAST::StructuralElementBase::side_external_force(bool request_jacobian,
                                                 libMesh::DenseVector<libMesh::Real> &f,
                                                 libMesh::DenseMatrix<libMesh::Real> &jac,
                                                 std::multimap<libMesh::boundary_id_type, MAST::BoundaryCondition*>& bc) {
    
    // iterate over the boundary ids given in the provided force map
    std::pair<std::multimap<libMesh::boundary_id_type, MAST::BoundaryCondition*>::const_iterator,
    std::multimap<libMesh::boundary_id_type, MAST::BoundaryCondition*>::const_iterator> it;
    
    const BoundaryInfo& binfo = *_system.get_mesh().boundary_info;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    for (unsigned short int n=0; n<_elem.n_sides(); n++) {
        if (!binfo.n_boundary_ids(&_elem, n))
            continue;
        
        std::vector<libMesh::boundary_id_type> bc_ids = binfo.boundary_ids(&_elem, n);
        std::vector<libMesh::boundary_id_type>::const_iterator bc_it = bc_ids.begin();
        for ( ; bc_it != bc_ids.end(); bc_it++) {
            // find the loads on this boundary and evaluate the f and jac
            it =bc.equal_range(*bc_it);
            
            for ( ; it.first != it.second; it.first++) {
                // apply all the types of loading
                switch (it.first->second->type()) {
                    case MAST::SURFACE_PRESSURE:

                        calculate_jac = (calculate_jac ||
                                         surface_pressure_force(request_jacobian,
                                                                f, jac,
                                                                n,
                                                                *it.first->second));
                        break;
                        
                    case MAST::SMALL_DISTURBANCE_MOTION:
                        calculate_jac = (calculate_jac ||
                                         small_disturbance_surface_pressure_force<ValType>(request_jacobian,
                                                                                           f, jac,
                                                                                           n,
                                                                                           *it.first->second));
                        break;

                    case MAST::DISPLACEMENT_DIRICHLET:
                        // nothing to be done here
                        break;
                        
                    default:
                        // not implemented yet
                        libmesh_error();
                        break;
                }
            }
        }
    }
    return (request_jacobian && calculate_jac);
}



template <typename ValType>
bool
MAST::StructuralElementBase::side_external_force_sensitivity(bool request_jacobian,
                                                             libMesh::DenseVector<libMesh::Real> &f,
                                                             libMesh::DenseMatrix<libMesh::Real> &jac,
                                                             std::multimap<libMesh::boundary_id_type, MAST::BoundaryCondition*>& bc) {
    
    // iterate over the boundary ids given in the provided force map
    std::pair<std::multimap<libMesh::boundary_id_type, MAST::BoundaryCondition*>::const_iterator,
    std::multimap<libMesh::boundary_id_type, MAST::BoundaryCondition*>::const_iterator> it;
    
    const BoundaryInfo& binfo = *_system.get_mesh().boundary_info;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    for (unsigned short int n=0; n<_elem.n_sides(); n++) {
        if (!binfo.n_boundary_ids(&_elem, n))
            continue;
        
        std::vector<libMesh::boundary_id_type> bc_ids = binfo.boundary_ids(&_elem, n);
        std::vector<libMesh::boundary_id_type>::const_iterator bc_it = bc_ids.begin();
        for ( ; bc_it != bc_ids.end(); bc_it++) {
            // find the loads on this boundary and evaluate the f and jac
            it =bc.equal_range(*bc_it);
            
            for ( ; it.first != it.second; it.first++) {
                // apply all the types of loading
                switch (it.first->second->type()) {

                    case MAST::SURFACE_PRESSURE:
                        calculate_jac = (calculate_jac ||
                                         surface_pressure_force_sensitivity(request_jacobian,
                                                                            f, jac,
                                                                            n,
                                                                            *it.first->second));
                        break;
                        
                    case MAST::DISPLACEMENT_DIRICHLET:
                        // nothing to be done here
                        break;
                        
                    default:
                        // not implemented yet
                        libmesh_error();
                        break;
                }
            }
        }
    }
    return (request_jacobian && calculate_jac);
}



template <typename ValType>
bool
MAST::StructuralElementBase::volume_external_force(bool request_jacobian,
                                                   libMesh::DenseVector<libMesh::Real> &f,
                                                   libMesh::DenseMatrix<libMesh::Real> &jac,
                                                   std::multimap<libMesh::subdomain_id_type, MAST::BoundaryCondition*>& bc) {
    // iterate over the boundary ids given in the provided force map
    std::pair<std::multimap<libMesh::subdomain_id_type, MAST::BoundaryCondition*>::const_iterator,
    std::multimap<libMesh::subdomain_id_type, MAST::BoundaryCondition*>::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    libMesh::subdomain_id_type sid = _elem.subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {

            case MAST::SURFACE_PRESSURE:
                calculate_jac = (calculate_jac ||
                                 surface_pressure_force(request_jacobian,
                                                        f, jac,
                                                        *it.first->second));
                break;
            case MAST::TEMPERATURE:
                calculate_jac = (calculate_jac ||
                                 thermal_force(request_jacobian,
                                               f, jac,
                                               *it.first->second));
                break;

            case MAST::SMALL_DISTURBANCE_MOTION:
                calculate_jac = (calculate_jac ||
                                 small_disturbance_surface_pressure_force<ValType>(request_jacobian,
                                                                                   f, jac,
                                                                                   *it.first->second));
                break;

            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    
    return (request_jacobian && calculate_jac);
}




template <typename ValType>
bool
MAST::StructuralElementBase::volume_external_force_sensitivity(bool request_jacobian,
                                                               libMesh::DenseVector<libMesh::Real> &f,
                                                               libMesh::DenseMatrix<libMesh::Real> &jac,
                                                               std::multimap<libMesh::subdomain_id_type, MAST::BoundaryCondition*>& bc) {
    // iterate over the boundary ids given in the provided force map
    std::pair<std::multimap<libMesh::subdomain_id_type, MAST::BoundaryCondition*>::const_iterator,
    std::multimap<libMesh::subdomain_id_type, MAST::BoundaryCondition*>::const_iterator> it;
    
    // for each boundary id, check if any of the sides on the element
    // has the associated boundary
    bool calculate_jac = false;
    
    libMesh::subdomain_id_type sid = _elem.subdomain_id();
    // find the loads on this boundary and evaluate the f and jac
    it =bc.equal_range(sid);
    
    for ( ; it.first != it.second; it.first++) {
        // apply all the types of loading
        switch (it.first->second->type()) {

            case MAST::SURFACE_PRESSURE:
                calculate_jac = (calculate_jac ||
                                 surface_pressure_force_sensitivity(request_jacobian,
                                                                    f, jac,
                                                                    *it.first->second));
                break;

            case MAST::TEMPERATURE:
                calculate_jac = (calculate_jac ||
                                 thermal_force_sensitivity(request_jacobian,
                                                           f, jac,
                                                           *it.first->second));
                break;
                
            default:
                // not implemented yet
                libmesh_error();
                break;
        }
    }
    return (request_jacobian && calculate_jac);
}




bool
MAST::StructuralElementBase::surface_pressure_force(bool request_jacobian,
                                                    libMesh::DenseVector<libMesh::Real> &f,
                                                    libMesh::DenseMatrix<libMesh::Real> &jac,
                                                    const unsigned int side,
                                                    MAST::BoundaryCondition &p) {

    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    FEMOperatorMatrix Bmat;
    
    // get the function from this boundary condition
    MAST::FieldFunction<libMesh::Real>& func =
    dynamic_cast<MAST::FieldFunction<libMesh::Real>&>(p.function());
    std::auto_ptr<libMesh::FEBase> fe;
    std::auto_ptr<libMesh::QBase> qrule;
    _get_side_fe_and_qrule(this->get_elem_for_quadrature(), side, fe, qrule);
    
    const std::vector<libMesh::Real> &JxW = fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = fe->get_xyz();
    const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();
    const unsigned int n_phi = (unsigned int)phi.size();
    const unsigned int n1=3, n2=6*n_phi;
    
    // boundary normals
    const std::vector<Point>& face_normals = fe->get_normals();
    libMesh::Real press;
    
    libMesh::DenseVector<libMesh::Real> phi_vec, force, local_f, tmp_vec_n2;
    phi_vec.resize(n_phi); force.resize(2*n1); local_f.resize(n2);
    tmp_vec_n2.resize(n2);
    
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++)
    {
        this->global_coordinates(qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure value
        func(pt, _system.time, press);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = press * face_normals[qp](i_dim);
        
        Bmat.vector_mult_transpose(tmp_vec_n2, force);
        
        local_f.add(-JxW[qp], tmp_vec_n2);
    }
    
    // now transform to the global system and add
    if (_elem.dim() < 3) {
        transform_to_global_system(local_f, tmp_vec_n2);
        f.add(1., tmp_vec_n2);
    }
    else
        f.add(1., local_f);
    

    return (request_jacobian && follower_forces);
}




bool
MAST::StructuralElementBase::surface_pressure_force(bool request_jacobian,
                                                    libMesh::DenseVector<libMesh::Real> &f,
                                                    libMesh::DenseMatrix<libMesh::Real> &jac,
                                                    MAST::BoundaryCondition &p) {

    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    FEMOperatorMatrix Bmat;
    
    // get the function from this boundary condition
    MAST::FieldFunction<libMesh::Real>& func =
    dynamic_cast<MAST::FieldFunction<libMesh::Real>&>(p.function());
    const std::vector<libMesh::Real> &JxW = _fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = _fe->get_xyz();
    const std::vector<std::vector<libMesh::Real> >& phi = _fe->get_phi();
    const unsigned int n_phi = (unsigned int)phi.size();
    const unsigned int n1=3, n2=6*n_phi;
    
    // normal for face integration
    libMesh::Point normal;
    // direction of pressure assumed to be normal (along local z-axis)
    // to the element face for 2D and along local y-axis for 1D element.
    normal(_elem.dim()) = 1.;
    
    libMesh::Real press;
    
    libMesh::DenseVector<libMesh::Real> phi_vec, force, local_f, tmp_vec_n2;
    phi_vec.resize(n_phi); force.resize(2*n1); local_f.resize(n2);
    tmp_vec_n2.resize(n2);
    
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++)
    {
        
        this->global_coordinates(qpoint[qp], pt);

        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure value
        func(pt, _system.time, press);
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = press * normal(i_dim);
        
        Bmat.vector_mult_transpose(tmp_vec_n2, force);
        
        local_f.add(-JxW[qp], tmp_vec_n2);
    }
    
    // now transform to the global system and add
    transform_to_global_system(local_f, tmp_vec_n2);
    f.add(1., tmp_vec_n2);
    

    return (request_jacobian && follower_forces);
}


template <typename ValType>
bool
MAST::StructuralElementBase::small_disturbance_surface_pressure_force(bool request_jacobian,
                                                                      libMesh::DenseVector<libMesh::Real> &f,
                                                                      libMesh::DenseMatrix<libMesh::Real> &jac,
                                                                      const unsigned int side,
                                                                      MAST::BoundaryCondition &p) {
//#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    libmesh_assert_equal_to(p.type(), MAST::SMALL_DISTURBANCE_MOTION);
    
    MAST::SmallDisturbanceMotion& sd_motion = dynamic_cast<MAST::SmallDisturbanceMotion&>(p);
    MAST::SurfaceMotionBase& surf_motion = sd_motion.get_deformation();
    MAST::SmallDisturbanceSurfacePressure& surf_press = sd_motion.get_pressure();
    
    FEMOperatorMatrix Bmat;
    
    // get the function from this boundary condition
    std::auto_ptr<libMesh::FEBase> fe;
    std::auto_ptr<libMesh::QBase> qrule;
    _get_side_fe_and_qrule(this->local_elem().local_elem(), side, fe, qrule);
    
    const std::vector<libMesh::Real> &JxW = fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = fe->get_xyz();
    const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();
    const unsigned int n_phi = (unsigned int)phi.size();
    const unsigned int n1=3, n2=6*n_phi;
    
    // boundary normals
    const std::vector<Point>& face_normals = fe->get_normals();
    
    libMesh::Real press;
    ValType dpress;
    libMesh::DenseVector<libMesh::Real> phi_vec;
    libMesh::DenseVector<ValType> wtrans, utrans, dn_rot, force, local_f, tmp_vec_n2;
    wtrans.resize(3); utrans.resize(3); dn_rot.resize(3);
    phi_vec.resize(n_phi); force.resize(2*n1); local_f.resize(n2);
    tmp_vec_n2.resize(n2);
    
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++)
    {
        this->global_coordinates(qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);

        // get pressure and deformation information
        surf_press.surface_pressure<ValType>(_system.time, pt, press, dpress);
        surf_motion.surface_velocity(_system.time, pt, face_normals[qp],
                                     wtrans, utrans, dn_rot);
        
        //            press = 0.;
        //            dpress = Complex(2./4.*std::real(dn_rot(0)),  2./4./.1*std::imag(utrans(1)));
        //            std::cout << q_point[qp](0)
        //            << "  " << std::real(utrans(1))
        //            << "  " << std::imag(utrans(1))
        //            << "  " << std::real(dn_rot(0))
        //            << "  " << std::imag(dn_rot(0))
        //            << "  " << std::real(press)
        //            << "  " << std::imag(press)
        //            << "  " << std::real(dpress)
        //            << "  " << std::imag(dpress) << std::endl;

        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) =  ( press * dn_rot(i_dim) + // steady pressure
                             dpress * face_normals[qp](i_dim)); // unsteady pressure

        
        Bmat.vector_mult_transpose(tmp_vec_n2, force);
        
        local_f.add(-JxW[qp], tmp_vec_n2);
    }
    
    // now transform to the global system and add
    if (_elem.dim() < 3) {
        transform_to_global_system(local_f, tmp_vec_n2);
        MAST::add_to_assembled_vector(f, tmp_vec_n2);
    }
    else
        MAST::add_to_assembled_vector(f, local_f);
    

    return (request_jacobian && follower_forces);
}




template <typename ValType>
bool
MAST::StructuralElementBase::small_disturbance_surface_pressure_force(bool request_jacobian,
                                                                      libMesh::DenseVector<libMesh::Real> &f,
                                                                      libMesh::DenseMatrix<libMesh::Real> &jac,
                                                                      MAST::BoundaryCondition &p) {
//#ifdef LIBMESH_USE_COMPLEX_NUMBERS
    libmesh_assert(_elem.dim() < 3); // only applicable for lower dimensional elements.
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    libmesh_assert_equal_to(p.type(), MAST::SMALL_DISTURBANCE_MOTION);

    MAST::SmallDisturbanceMotion& sd_motion = dynamic_cast<MAST::SmallDisturbanceMotion&>(p);
    MAST::SurfaceMotionBase& surf_motion = sd_motion.get_deformation();
    MAST::SmallDisturbanceSurfacePressure& surf_press = sd_motion.get_pressure();

    FEMOperatorMatrix Bmat;
    
    // get the function from this boundary condition
    const std::vector<libMesh::Real> &JxW = _fe->get_JxW();
    
    // Physical location of the quadrature points
    const std::vector<Point>& qpoint = _fe->get_xyz();
    const std::vector<std::vector<libMesh::Real> >& phi = _fe->get_phi();
    const unsigned int n_phi = (unsigned int)phi.size();
    const unsigned int n1=3, n2=6*n_phi;
    
    // normal for face integration
    libMesh::Point normal;
    // direction of pressure assumed to be normal (along local z-axis)
    // to the element face for 2D and along local y-axis for 1D element.
    normal(_elem.dim()) = 1.;
    
    libMesh::Real press;
    ValType dpress;
    libMesh::DenseVector<libMesh::Real> phi_vec;
    libMesh::DenseVector<ValType> wtrans, utrans, dn_rot, force, local_f, tmp_vec_n2;
    wtrans.resize(3); utrans.resize(3); dn_rot.resize(3);
    phi_vec.resize(n_phi); force.resize(2*n1); local_f.resize(n2);
    tmp_vec_n2.resize(n2);
    
    libMesh::Point pt;
    
    for (unsigned int qp=0; qp<qpoint.size(); qp++)
    {
        this->global_coordinates(qpoint[qp], pt);
        
        // now set the shape function values
        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            phi_vec(i_nd) = phi[i_nd][qp];
        
        Bmat.reinit(2*n1, phi_vec);
        
        // get pressure and deformation information
        surf_press.surface_pressure<ValType>(_system.time, pt, press, dpress);
        surf_motion.surface_velocity(_system.time, pt, normal,
                                     wtrans, utrans, dn_rot);
//        std::cout << std::setw(15) << pt(0)
//        << std::setw(15) << std::real(press)
//        << std::setw(15) << std::imag(press)
//        << std::setw(15) << std::real(dpress)
//        << std::setw(15) << std::imag(dpress)
//        << std::setw(15) << std::real(utrans(1))
//        << std::setw(15) << std::imag(utrans(1))
//        << std::setw(15) << std::real(dn_rot(0))
//        << std::setw(15) << std::imag(dn_rot(0)) << std::endl;
        
        
        // calculate force
        for (unsigned int i_dim=0; i_dim<n1; i_dim++)
            force(i_dim) = ( press * dn_rot(i_dim) + // steady pressure
                            dpress * normal(i_dim)); // unsteady pressure
        
        Bmat.vector_mult_transpose(tmp_vec_n2, force);
        
        local_f.add(-JxW[qp], tmp_vec_n2);
    }
    
    // now transform to the global system and add
    transform_to_global_system(local_f, tmp_vec_n2);
    MAST::add_to_assembled_vector(f, tmp_vec_n2);
    

    return (request_jacobian && follower_forces);
}



bool
MAST::StructuralElementBase::surface_pressure_force_sensitivity(bool request_jacobian,
                                                                libMesh::DenseVector<libMesh::Real> &f,
                                                                libMesh::DenseMatrix<libMesh::Real> &jac,
                                                                const unsigned int side,
                                                                MAST::BoundaryCondition &p) {
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    // currently not implemented.
    return false;
    
    // this should be true if the function is called
    libmesh_assert(this->sensitivity_param);
    libmesh_assert(!this->sensitivity_param->is_shape_parameter()); // this is not implemented for now
    
    
    // check if the material property or the provided exterior
    // values, like temperature, are functions of the sensitivity parameter
    bool calculate = false;
    calculate = calculate || _property.depends_on(*(this->sensitivity_param));
    
    // nothing to be calculated if the element does not depend on the
    // sensitivity parameter.
    if (!calculate)
        return false;
    
    
    return (request_jacobian && follower_forces);
}




bool
MAST::StructuralElementBase::surface_pressure_force_sensitivity(bool request_jacobian,
                                                                libMesh::DenseVector<libMesh::Real> &f,
                                                                libMesh::DenseMatrix<libMesh::Real> &jac,
                                                                MAST::BoundaryCondition &p) {
    libmesh_assert(!follower_forces); // not implemented yet for follower forces
    
    // currently not implemented.
    return false;
    
    // this should be true if the function is called
    libmesh_assert(this->sensitivity_param);
    libmesh_assert(!this->sensitivity_param->is_shape_parameter()); // this is not implemented for now
    
    
    // check if the material property or the provided exterior
    // values, like temperature, are functions of the sensitivity parameter
    bool calculate = false;
    calculate = calculate || _property.depends_on(*(this->sensitivity_param));
    
    // nothing to be calculated if the element does not depend on the
    // sensitivity parameter.
    if (!calculate)
        return false;
    
    
    return (request_jacobian && follower_forces);
}





template <typename ValType>
void
MAST::StructuralElementBase::transform_to_global_system(const libMesh::DenseMatrix<ValType>& local_mat,
                                                        libMesh::DenseMatrix<ValType>& global_mat) const {
    libmesh_assert_equal_to( local_mat.m(),  local_mat.n());
    libmesh_assert_equal_to(global_mat.m(), global_mat.n());
    libmesh_assert_equal_to( local_mat.m(), global_mat.m());
    
    const unsigned int n_dofs = _fe->n_shape_functions();
    global_mat.zero();
    libMesh::DenseMatrix<ValType> tmp_mat;
    tmp_mat.resize(6*n_dofs, 6*n_dofs);
    
    const libMesh::DenseMatrix<libMesh::Real>& Tmat = _transformation_matrix();
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



template <typename ValType>
void
MAST::StructuralElementBase::transform_to_local_system(const libMesh::DenseVector<ValType>& global_vec,
                                                        libMesh::DenseVector<ValType>& local_vec) const {
    libmesh_assert_equal_to( local_vec.size(),  global_vec.size());
    
    const unsigned int n_dofs = _fe->n_shape_functions();
    local_vec.zero();
    libMesh::DenseMatrix<libMesh::Real> tmp_mat;
    tmp_mat.resize(6*n_dofs, 6*n_dofs);
    
    const libMesh::DenseMatrix<libMesh::Real>& Tmat = _transformation_matrix();
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



template <typename ValType>
void
MAST::StructuralElementBase::transform_to_global_system(const libMesh::DenseVector<ValType>& local_vec,
                                                         libMesh::DenseVector<ValType>& global_vec) const {
    libmesh_assert_equal_to( local_vec.size(),  global_vec.size());
    
    const unsigned int n_dofs = _fe->n_shape_functions();
    global_vec.zero();
    libMesh::DenseMatrix<libMesh::Real> tmp_mat;
    tmp_mat.resize(6*n_dofs, 6*n_dofs);
    
    const libMesh::DenseMatrix<libMesh::Real>& Tmat = _transformation_matrix();
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
MAST::StructuralElementBase::_init_fe_and_qrule( const libMesh::Elem& e) {
    
    unsigned int nv = _system.n_vars();
    
    libmesh_assert (nv);
    FEType fe_type = _system.variable_type(0); // all variables are assumed to be of same type
    
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _system.variable_type(i));
    
    // Create an adequate quadrature rule
    _fe.reset(libMesh::FEBase::build(e.dim(), fe_type).release());
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
MAST::StructuralElementBase::_get_side_fe_and_qrule(const libMesh::Elem& e,
                                                    unsigned int s,
                                                    std::auto_ptr<libMesh::FEBase>& fe,
                                                    std::auto_ptr<libMesh::QBase>& qrule) {
    unsigned int nv = _system.n_vars();
    
    libmesh_assert (nv);
    FEType fe_type = _system.variable_type(0); // all variables are assumed to be of same type
    
    
    for (unsigned int i=1; i != nv; ++i)
        libmesh_assert(fe_type == _system.variable_type(i));
    
    // Create an adequate quadrature rule
    fe.reset(libMesh::FEBase::build(e.dim(), fe_type).release());
    qrule.reset(fe_type.default_quadrature_rule
                (e.dim()-1,
                 _system.extra_quadrature_order).release());  // system extra quadrature
    fe->attach_quadrature_rule(qrule.get());
    fe->get_phi();
    fe->get_JxW();
    
    fe->reinit(&e, s);
}




std::auto_ptr<MAST::StructuralElementBase>
MAST::build_structural_element(libMesh::System& sys,
                               const libMesh::Elem& elem,
                               const MAST::ElementPropertyCardBase& p) {
    
    std::auto_ptr<MAST::StructuralElementBase> e;
    
    switch (elem.dim()) {
        case 1:
            e.reset(new MAST::StructuralElement1D(sys, elem, p));
            break;
            
        case 2:
            e.reset(new MAST::StructuralElement2D(sys, elem, p));
            break;
            
        case 3:
            e.reset(new MAST::StructuralElement3D(sys, elem, p));
            break;

        default:
            libmesh_error();
            break;
    }
    
    return e;
}



// template instantiations
template
void
MAST::StructuralElementBase::transform_to_global_system<libMesh::Real>(const libMesh::DenseMatrix<libMesh::Real>& local_mat,
                                                              libMesh::DenseMatrix<libMesh::Real>& global_mat) const;

template
void
MAST::StructuralElementBase::transform_to_local_system<libMesh::Real>(const libMesh::DenseVector<libMesh::Real>& global_vec,
                                                             libMesh::DenseVector<libMesh::Real>& local_vec) const;

template
void
MAST::StructuralElementBase::transform_to_global_system<libMesh::Real>(const libMesh::DenseVector<libMesh::Real>& local_vec,
                                                              libMesh::DenseVector<libMesh::Real>& global_vec) const;

//#ifdef LIBMESH_USE_COMPLEX_NUMBERS

template
void
MAST::StructuralElementBase::transform_to_global_system<libMesh::Complex>(const libMesh::DenseMatrix<libMesh::Complex>& local_mat,
                                                              libMesh::DenseMatrix<libMesh::Complex>& global_mat) const;

template
void
MAST::StructuralElementBase::transform_to_local_system<libMesh::Complex>(const libMesh::DenseVector<libMesh::Complex>& global_vec,
                                                             libMesh::DenseVector<libMesh::Complex>& local_vec) const;

template
void
MAST::StructuralElementBase::transform_to_global_system<libMesh::Complex>(const libMesh::DenseVector<libMesh::Complex>& local_vec,
                                                              libMesh::DenseVector<libMesh::Complex>& global_vec) const;

template
bool
MAST::StructuralElementBase::side_external_force<libMesh::Real>(bool request_jacobian,
                                                                libMesh::DenseVector<libMesh::Real> &f,
                                                                libMesh::DenseMatrix<libMesh::Real> &jac,
                                                                std::multimap<libMesh::boundary_id_type,
                                                                MAST::BoundaryCondition *> &bc);

template
bool
MAST::StructuralElementBase::side_external_force<libMesh::Complex>(bool request_jacobian,
                                                                   libMesh::DenseVector<libMesh::Real> &f,
                                                                   libMesh::DenseMatrix<libMesh::Real> &jac,
                                                                   std::multimap<libMesh::boundary_id_type,
                                                                   MAST::BoundaryCondition *> &bc);


template
bool
MAST::StructuralElementBase::side_external_force_sensitivity<libMesh::Real>(bool request_jacobian,
                                                                            libMesh::DenseVector<libMesh::Real> &f,
                                                                            libMesh::DenseMatrix<libMesh::Real> &jac,
                                                                            std::multimap<libMesh::boundary_id_type,
                                                                            MAST::BoundaryCondition *> &bc);

template
bool
MAST::StructuralElementBase::side_external_force_sensitivity<libMesh::Complex>(bool request_jacobian,
                                                                               libMesh::DenseVector<libMesh::Real> &f,
                                                                               libMesh::DenseMatrix<libMesh::Real> &jac,
                                                                               std::multimap<libMesh::boundary_id_type,
                                                                               MAST::BoundaryCondition *> &bc);

template
bool
MAST::StructuralElementBase::volume_external_force<libMesh::Real>(bool request_jacobian,
                                                                  libMesh::DenseVector<libMesh::Real>& f,
                                                                  libMesh::DenseMatrix<libMesh::Real>& jac,
                                                                  std::multimap<libMesh::subdomain_id_type, MAST::BoundaryCondition*>& bc);

template
bool
MAST::StructuralElementBase::volume_external_force<libMesh::Complex>(bool request_jacobian,
                                                                     libMesh::DenseVector<libMesh::Real>& f,
                                                                     libMesh::DenseMatrix<libMesh::Real>& jac,
                                                                     std::multimap<libMesh::subdomain_id_type, MAST::BoundaryCondition*>& bc);


template
bool
MAST::StructuralElementBase::volume_external_force_sensitivity<libMesh::Real>(bool request_jacobian,
                                                                              libMesh::DenseVector<libMesh::Real>& f,
                                                                              libMesh::DenseMatrix<libMesh::Real>& jac,
                                                                              std::multimap<libMesh::subdomain_id_type, MAST::BoundaryCondition*>& bc);

template
bool
MAST::StructuralElementBase::volume_external_force_sensitivity<libMesh::Complex>(bool request_jacobian,
                                                                                 libMesh::DenseVector<libMesh::Real>& f,
                                                                                 libMesh::DenseMatrix<libMesh::Real>& jac,
                                                                                 std::multimap<libMesh::subdomain_id_type, MAST::BoundaryCondition*>& bc);



