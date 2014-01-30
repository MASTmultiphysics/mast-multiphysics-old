//
//  solid_1d_section_element_property_card.cpp
//  MAST
//
//  Created by Manav Bhatia on 1/30/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "PropertyCards/solid_1d_section_element_property_card.h"



MAST::Solid1DSectionElementPropertyCard::SectionIntegratedExtensionStiffnessMatrix::
SectionIntegratedExtensionStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat,
                                          MAST::FieldFunction<Real>& h):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedExtensionStiffnessMatrix1D"),
_material_stiffness(mat),
_h(h) {
    _functions.insert(mat.master());
    _functions.insert(h.master());
}





void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedExtensionStiffnessMatrix::operator() (const Point& p,
                                                       const Real t,
                                                       DenseMatrix<Real>& v) const {
    // [C]*h
    Real h;
    _h(p, t, h);
    _material_stiffness(p, t, v);
    v.scale(h);
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedExtensionStiffnessMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                               const Point& p,
                                                               const Real t,
                                                               DenseMatrix<Real>& v) const {
    
    libmesh_assert_equal_to(par.total_order(), 1);
    const MAST::FunctionBase& f = par.get_first_order_derivative_parameter();
    
    v.resize(3,3);
    
    // [C]*dh
    if (_h.is_equal(f))
        _material_stiffness(p, t, v);
    
    // += [dC]*h
    DenseMatrix<Real> m;
    Real h;
    _h(p, t, h);
    _material_stiffness.partial_derivative(par, p, t, m);
    v.add(h, m);
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedExtensionStiffnessMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                             const Point& p,
                                                             const Real t,
                                                             DenseMatrix<Real>& v) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    
    // sensitivity of h
    Real h, dh;
    _h.total_derivative(par, p, t, dh);
    if (fabs(dh) > 0.) {
        MAST::SensitivityParameters dhpar; dhpar.add_parameter(&_h, 1);
        this->partial_derivative(dhpar, p, t, v);
        v.scale(dh);
    }
    
    
    // sensitivity of C
    DenseMatrix<Real> dm;
    _h(p, t, h);
    _material_stiffness.total_derivative(par, p, t, dm);
    if (dm.linfty_norm() > 0.)
        v.add(h, dm);
}






MAST::Solid1DSectionElementPropertyCard::SectionIntegratedExtensionBendingStiffnessMatrix::
SectionIntegratedExtensionBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat,
                                                 MAST::FieldFunction<Real>& h):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedExtensionBendingStiffnessMatrix1D"),
_material_stiffness(mat),
_h(h) {
    _functions.insert(mat.master());
    _functions.insert(h.master());
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedExtensionBendingStiffnessMatrix::operator() (const Point& p,
                                                              const Real t,
                                                              DenseMatrix<Real>& v) const {
    // this is zero for solid sections unless an offset is added
    v.resize(3,3);
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedExtensionBendingStiffnessMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                                      const Point& p,
                                                                      const Real t,
                                                                      DenseMatrix<Real>& v) const {
    
    // this is zero for solid sections unless an offset is added
    v.resize(3,3);
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedExtensionBendingStiffnessMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                                    const Point& p,
                                                                    const Real t,
                                                                    DenseMatrix<Real>& v) const {
    // this is zero for solid sections unless an offset is added
    v.resize(3,3);
}




MAST::Solid1DSectionElementPropertyCard::SectionIntegratedBendingStiffnessMatrix::
SectionIntegratedBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat,
                                        MAST::FieldFunction<Real>& h):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedBendingStiffnessMatrix1D"),
_material_stiffness(mat),
_h(h) {
    _functions.insert(mat.master());
    _functions.insert(h.master());
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedBendingStiffnessMatrix::operator() (const Point& p,
                                                     const Real t,
                                                     DenseMatrix<Real>& v) const {
    // [C]*h
    Real h;
    _h(p, t, h);
    _material_stiffness(p, t, v);
    v.scale(pow(h,3)/12.);
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedBendingStiffnessMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                             const Point& p,
                                                             const Real t,
                                                             DenseMatrix<Real>& v) const {
    
    libmesh_assert_equal_to(par.total_order(), 1);
    const MAST::FunctionBase& f = par.get_first_order_derivative_parameter();
    
    v.resize(3,3);
    Real h;
    _h(p, t, h);
    
    // [C]*dh
    if (_h.is_equal(f)) {
        _material_stiffness(p, t, v);
        v.scale(pow(h,2)/4.);
    }
    
    
    // += [dC]*h
    DenseMatrix<Real> m;
    _material_stiffness.partial_derivative(par, p, t, m);
    v.add(pow(h,3)/12., m);
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedBendingStiffnessMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                           const Point& p,
                                                           const Real t,
                                                           DenseMatrix<Real>& v) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    
    // sensitivity of h
    Real h, dh;
    _h.total_derivative(par, p, t, dh);
    if (fabs(dh) > 0.) {
        MAST::SensitivityParameters dhpar; dhpar.add_parameter(&_h, 1);
        this->partial_derivative(dhpar, p, t, v);
        v.scale(dh);
    }
    
    
    // sensitivity of C
    DenseMatrix<Real> dm;
    _h(p, t, h);
    _material_stiffness.total_derivative(par, p, t, dm);
    if (dm.linfty_norm() > 0.)
        v.add(pow(h,3)/12., dm);
}





MAST::Solid1DSectionElementPropertyCard::SectionIntegratedInertiaMatrix::
SectionIntegratedInertiaMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat,
                               MAST::FieldFunction<Real>& h):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedInertiaMatrix1D"),
_material_inertia(mat),
_h(h) {
    _functions.insert(mat.master());
    _functions.insert(h.master());
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedInertiaMatrix::operator() (const Point& p,
                                            const Real t,
                                            DenseMatrix<Real>& v) const {
    DenseMatrix<Real> m;
    v.resize(6, 6);
    // this only returns the material inertia
    _material_inertia(p, t, m);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            v(i,j) = m(i,j);
        }
        v(i+3,i+3) = m(i,i) * 1.0e-12;
    }
    
    Real h;
    _h(p, t, h);
    v.scale(h);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedInertiaMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                    const Point& p,
                                                    const Real t,
                                                    DenseMatrix<Real>& v) const {
    
    libmesh_assert_equal_to(par.total_order(), 1);
    const MAST::FunctionBase& f = par.get_first_order_derivative_parameter();
    
    DenseMatrix<Real> m;
    v.resize(6, 6);
    
    
    // [C]*dh
    if (_h.is_equal(f))
        _material_inertia(p, t, v);
    
    // += [dC]*h
    Real h;
    _h(p, t, h);
    _material_inertia.partial_derivative(par, p, t, m);
    for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<3; j++)
            v(i,j) += h*m(i,j);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedInertiaMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                  const Point& p,
                                                  const Real t,
                                                  DenseMatrix<Real>& v) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    
    v.resize(6, 6);
    
    Real h, dh;
    _h.total_derivative(par, p, t, dh);
    if (fabs(dh) > 0.) {
        MAST::SensitivityParameters dhpar; dhpar.add_parameter(&_h, 1);
        this->partial_derivative(dhpar, p, t, v);
        v.scale(dh);
    }
    
    // += [dC]*h
    DenseMatrix<Real> m;
    _h(p, t, h);
    _material_inertia.total_derivative(par, p, t, m);
    for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<3; j++)
            v(i,j) += h*m(i,j);
}




MAST::Solid1DSectionElementPropertyCard::SectionIntegratedThermalExpansionMatrix::
SectionIntegratedThermalExpansionMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat_stiff,
                                        MAST::FieldFunction<DenseMatrix<Real> > &mat_expansion,
                                        MAST::FieldFunction<Real>& h):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedThermalExpansionMatrix1D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_h(h) {
    _functions.insert(mat_stiff.master());
    _functions.insert(mat_expansion.master());
}




void
MAST::Solid1DSectionElementPropertyCard::SectionIntegratedThermalExpansionMatrix::operator() (const Point& p,
                                                                                              const Real t,
                                                                                              DenseMatrix<Real>& v) const {
    DenseMatrix<Real> m;
    _material_stiffness(p, t, v);
    _material_expansion(p, t, m);
    v.right_multiply(m);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedThermalExpansionMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                             const Point& p,
                                                             const Real t,
                                                             DenseMatrix<Real>& v) const {
    libmesh_error();
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedThermalExpansionMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                           const Point& p,
                                                           const Real t,
                                                           DenseMatrix<Real>& v) const {
    libmesh_error();
}




MAST::Solid1DSectionElementPropertyCard::SectionIntegratedPrestressMatrix::
SectionIntegratedPrestressMatrix(MAST::FieldFunction<DenseMatrix<Real> > &prestress,
                                 MAST::FieldFunction<Real>& h):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedPrestressMatrix1D"),
_prestress(prestress),
_h(h) {
    _functions.insert(prestress.master());
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedPrestressMatrix::operator() (const Point& p,
                                              const Real t,
                                              DenseMatrix<Real>& v) const {
    libmesh_error();
}





void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedPrestressMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                      const Point& p,
                                                      const Real t,
                                                      DenseMatrix<Real>& v) const {
    libmesh_error();
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedPrestressMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                    const Point& p,
                                                    const Real t,
                                                    DenseMatrix<Real>& v) const {
    libmesh_error();
}



//void
//MAST::Solid1DSectionElementPropertyCard::SectionIntegratedPrestressMatrix::prestress_vector(MAST::ElemenetPropertyMatrixType t,
//                                                                                const DenseMatrix<Real>& T,
//                                                                                DenseVector<Real>& v) const {
//    v.resize(6); // zero, if the stress has not been defined
//    if (_prestress.m() != 0) {
//        for (unsigned int i=0; i<3; i++)
//            v(i) = _prestress(i,i);
//        v(3) = _prestress(0,1);  // tau_xy
//        v(4) = _prestress(1,2);  // tau_yz
//        v(5) = _prestress(0,2);  // tau_xz
//    }
//}






inline Real
MAST::Solid1DSectionElementPropertyCard::value(const std::string& val) const {
    
    Real h_y = this->get<Real>("h_y")(), // section height
    h_z = this->get<Real>("h_z")(),        // section width
    off_h_y = 0., off_h_z = 0.;
    
    if (_properties.count("off_h_y"))
        off_h_y = this->get<Real>("off_h_y")();
    if (_properties.count("off_h_z"))
        off_h_z = this->get<Real>("off_h_z")();
    
    Real Area = h_y*h_z,
    Iyy = h_y*(pow(h_z,3)/12. + pow(off_h_z,2)*h_z),
    Izz = h_z*(pow(h_y,3)/12. + pow(off_h_y,2)*h_y),
    Iyz = (off_h_y*off_h_z*h_y*h_z),
    J=1.;
    
    if (val == "A")
        return Area;
    else if (val == "IYY")
        return Iyy;
    else if (val == "IZZ")
        return Izz;
    else if (val == "IYZ")
        return Iyz;
    else if (val == "J")
        return J;
    
    // should not get here
    libmesh_error();
    return 0.;
}




inline void
MAST::Solid1DSectionElementPropertyCard::calculate_matrix(const libMesh::Elem &elem,
                                                          MAST::ElemenetPropertyMatrixType t,
                                                          DenseMatrix<Real>& m) const
{
    libmesh_assert(_material); // should have been set
    
    switch (elem.dim()) {
            
        case 1: {
            Real h_y = this->get<Real>("h_y")(), // section height
            h_z = this->get<Real>("h_z")(),      // section width
            off_h_y = 0., off_h_z = 0.;          // offset values of mid-plane from longitudinal axis
            
            if (_properties.count("off_h_y"))
                off_h_y = this->get<Real>("off_h_y")();
            if (_properties.count("off_h_z"))
                off_h_z = this->get<Real>("off_h_z")();
            
            Real Area = h_y*h_z,
            Iyy = h_y*(pow(h_z,3)/12. + pow(off_h_z,2)*h_z),   // used for w-bending
            Izz = h_z*(pow(h_y,3)/12. + pow(off_h_y,2)*h_y),   // used for v-bending
            Iyz = (off_h_y*off_h_z*h_y*h_z),
            J=1.;
            
            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX: {
                    _material->calculate_1d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m);
                    m.scale_row(0, Area);
                    m.scale_row(1, J);
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX: {
                    // for solid sections with isotropic material this is zero
                    // m(0, 0) and m(0, 1) identify coupling of extension and bending
                    // m(1, 0) and m(1, 1) identify coupling of torsion and bending
                    Real E = _material->get<Real>("E")();
                    m.resize(2,2);
                    m(0,0) = E*h_y*h_z*off_h_y;  // u <-> v
                    m(0,1) = E*h_y*h_z*off_h_z;  // u <-> w
                    m(1,0) = 0.;           // tx <-> v
                    m(1,1) = 0.;           // tx <-> w
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX: {
                    m.resize(2,2);
                    Real E = _material->get<Real>("E")();
                    m(0,0) = E*Izz;
                    m(0,1) = E*Iyz;
                    m(1,0) = E*Iyz;
                    m(1,1) = E*Iyy;
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX:
                    _material->calculate_1d_matrix(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                   m);
                    m.scale(Area);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
                    _material->calculate_1d_matrix(MAST::MATERIAL_INERTIA_MATRIX,
                                                   m);
                    m.scale(Area);
                    // now scale the rotation dofs with small factors
                    for (unsigned int i=0; i<3; i++) {
                        m(i+3, i+3) *= 1.0e-6;
                    }
                    break;
                    
                default:
                    libmesh_error();
                    break;
            }
        }
            break;
            
        case 2:
        case 3:
        default:
            libmesh_error();
            break;
    }
}



inline void
MAST::Solid1DSectionElementPropertyCard::calculate_matrix_sensitivity(const libMesh::Elem &elem,
                                                                      MAST::ElemenetPropertyMatrixType t,
                                                                      DenseMatrix<Real>& m,
                                                                      const MAST::SensitivityParameters& p) const
{
    libmesh_assert(_material); // should have been set
    
    // only first order sensitivities are calculated at this point
    libmesh_assert_equal_to(p.total_order(), 1);
    const MAST::FunctionBase& f = p.get_first_order_derivative_parameter();
    
    DenseMatrix<Real> dm;
    
    const bool depends = this->depends_on(f);
    
    switch (elem.dim()) {
            
        case 1: {
            Real h_y = this->get<Real>("h_y")(), // section height
            h_z = this->get<Real>("h_z")(),        // section width
            off_h_y = 0., off_h_z = 0.;          // offset values of mid-plane from longitudinal axis
            
            if (_properties.count("off_h_y"))
                off_h_y = this->get<Real>("off_h_y")();
            if (_properties.count("off_h_z"))
                off_h_z = this->get<Real>("off_h_z")();
            
            Real Area = h_y*h_z,
            Iyy = h_y*(pow(h_z,3)/12. + pow(off_h_z,2)*h_z),   // used for w-bending
            Izz = h_z*(pow(h_y,3)/12. + pow(off_h_y,2)*h_y),   // used for v-bending
            Iyz = (off_h_y*off_h_z*h_y*h_z),
            
            dAreadhy = h_z, dAreadhz = h_y,
            
            dIyydhy = (pow(h_z,3)/12. + pow(off_h_z,2)*h_z),
            dIyydhz = h_y*(pow(h_z,2)/4. + pow(off_h_z,2)),
            dIzzdhy = h_z*(pow(h_y,2)/4. + pow(off_h_y,2)),
            dIzzdhz = (pow(h_y,3)/12. + pow(off_h_y,2)*h_y),
            
            dIyzdhy = off_h_y*off_h_z*h_z,
            dIyzdhz = off_h_y*off_h_z*h_y,
            
            J = 1., dJdhy = 0., dJdhz = 0.;
            
            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX: {
                    _material->calculate_1d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m);
                    _material->calculate_1d_matrix_sensitivity(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                               dm, p);
                    dm.scale_row(0, Area);
                    dm.scale_row(1, J);
                    
                    if (depends && f.name() == "h_y") {
                        m.scale_row(0, dAreadhy);
                        m.scale_row(1, dJdhy);
                    }
                    else if (depends && f.name() == "h_z") {
                        m.scale_row(0, dAreadhz);
                        m.scale_row(1, dJdhz);
                    }
                    else
                        m.zero();
                    
                    m.add(1., dm);
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
                    // for solid sections with isotropic material this is zero
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX: {
                    m.resize(2,2); dm.resize(2,2);
                    Real E = _material->get<Real>("E")(),
                    dEdp = _material->get<Real>("E").sensitivity(p);
                    dm(0,0) = dEdp*Iyy;
                    dm(0,1) = dEdp*Iyz;
                    dm(1,0) = dEdp*Iyz;
                    dm(1,1) = dEdp*Izz;
                    
                    if (depends && f.name() == "h_y") {
                        m(0,0) = E*dIyydhy;
                        m(0,1) = E*dIyzdhy;
                        m(1,0) = E*dIyzdhy;
                        m(1,1) = E*dIzzdhy;
                    }
                    else if (depends && f.name() == "h_z") {
                        m(0,0) = E*dIyydhz;
                        m(0,1) = E*dIyzdhz;
                        m(1,0) = E*dIyzdhz;
                        m(1,1) = E*dIzzdhz;
                    }
                    else {
                        m.zero();
                    }
                    
                    m.add(1., dm);
                }
                    break;
                    
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX: {
                    _material->calculate_1d_matrix(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                   m);
                    _material->calculate_1d_matrix_sensitivity(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                               dm, p);
                    if (depends && f.name() == "h_y")
                        m.scale(dAreadhy);
                    else if (depends && f.name() == "h_z")
                        m.scale(dAreadhz);
                    else
                        m.zero();
                    
                    m.add(Area, dm);
                }
                    break;
                    
                default:
                    libmesh_error();
                    break;
            }
        }
            break;
            
        case 2:
        case 3:
        default:
            libmesh_error();
            break;
    }
}



inline void
MAST::Solid1DSectionElementPropertyCard::prestress_vector(MAST::ElemenetPropertyMatrixType t,
                                                          const DenseMatrix<Real>& T,
                                                          DenseVector<Real>& v) const {
    // prestress is independent of cross-section.
    
    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            _prestress_vector(T, v);
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
            // for solid sections with isotropic material this is zero
            v.resize(2);
            break;
            
        default:
            libmesh_error();
            break;
    }
}




inline void
MAST::Solid1DSectionElementPropertyCard::prestress_vector_sensitivity(MAST::ElemenetPropertyMatrixType t,
                                                                      const DenseMatrix<Real>& T,
                                                                      DenseVector<Real>& v,
                                                                      const MAST::SensitivityParameters& p) const {
    // only first order sensitivities are calculated at this point
    libmesh_assert_equal_to(p.total_order(), 1);
    // prestress is independent of cross-section.
    v.resize(2);
}





inline void
MAST::Solid1DSectionElementPropertyCard::prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                                          const DenseMatrix<Real>& T,
                                                          DenseMatrix<Real>& m) const {
    // prestress is independent of cross-section.
    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            _prestress_matrix(T, m);
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
            // for solid sections with isotropic material this is zero
            m.resize(2,2);
            break;
            
        default:
            libmesh_error();
            break;
    }
}



inline void
MAST::Solid1DSectionElementPropertyCard::prestress_matrix_sensitivity(MAST::ElemenetPropertyMatrixType t,
                                                                      const DenseMatrix<Real>& T,
                                                                      DenseMatrix<Real>& m,
                                                                      const MAST::SensitivityParameters& p) const {
    // only first order sensitivities are calculated at this point
    libmesh_assert_equal_to(p.total_order(), 1);
    
    // prestress is independent of cross-section.
    m.resize(2,2);
}


