//
//  solid_2d_section_element_property_card.cpp
//  MAST
//
//  Created by Manav Bhatia on 1/30/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "PropertyCards/solid_2d_section_element_property_card.h"




MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionStiffnessMatrix::
SectionIntegratedExtensionStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat,
                                          MAST::FieldFunction<Real>& h):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedExtensionStiffnessMatrix2D"),
_material_stiffness(mat),
_h(h) {
    _functions.insert(mat.master());
    _functions.insert(h.master());
}



void
MAST::Solid2DSectionElementPropertyCard::
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
MAST::Solid2DSectionElementPropertyCard::
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
MAST::Solid2DSectionElementPropertyCard::
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






MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionBendingStiffnessMatrix::
SectionIntegratedExtensionBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat,
                                                 MAST::FieldFunction<Real>& h):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedExtensionBendingStiffnessMatrix2D"),
_material_stiffness(mat),
_h(h) {
    _functions.insert(mat.master());
    _functions.insert(h.master());
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedExtensionBendingStiffnessMatrix::operator() (const Point& p,
                                                              const Real t,
                                                              DenseMatrix<Real>& v) const {
    // this is zero for solid sections unless an offset is added
    v.resize(3,3);
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedExtensionBendingStiffnessMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                                      const Point& p,
                                                                      const Real t,
                                                                      DenseMatrix<Real>& v) const {
    
    // this is zero for solid sections unless an offset is added
    v.resize(3,3);
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedExtensionBendingStiffnessMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                                    const Point& p,
                                                                    const Real t,
                                                                    DenseMatrix<Real>& v) const {
    // this is zero for solid sections unless an offset is added
    v.resize(3,3);
}




MAST::Solid2DSectionElementPropertyCard::SectionIntegratedBendingStiffnessMatrix::
SectionIntegratedBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat,
                                        MAST::FieldFunction<Real>& h):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedBendingStiffnessMatrix2D"),
_material_stiffness(mat),
_h(h) {
    _functions.insert(mat.master());
    _functions.insert(h.master());
}



void
MAST::Solid2DSectionElementPropertyCard::
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
MAST::Solid2DSectionElementPropertyCard::
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
MAST::Solid2DSectionElementPropertyCard::
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





MAST::Solid2DSectionElementPropertyCard::SectionIntegratedInertiaMatrix::
SectionIntegratedInertiaMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat,
                               MAST::FieldFunction<Real>& h):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedInertiaMatrix2D"),
_material_inertia(mat),
_h(h) {
    _functions.insert(mat.master());
    _functions.insert(h.master());
}




void
MAST::Solid2DSectionElementPropertyCard::
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
MAST::Solid2DSectionElementPropertyCard::
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
MAST::Solid2DSectionElementPropertyCard::
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




MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionMatrix::
SectionIntegratedThermalExpansionMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat_stiff,
                                        MAST::FieldFunction<DenseMatrix<Real> > &mat_expansion,
                                        MAST::FieldFunction<Real>& h):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedThermalExpansionMatrix2D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_h(h) {
    _functions.insert(mat_stiff.master());
    _functions.insert(mat_expansion.master());
}




void
MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionMatrix::operator() (const Point& p,
                                                                                  const Real t,
                                                                                  DenseMatrix<Real>& v) const {
    DenseMatrix<Real> m;
    _material_stiffness(p, t, v);
    _material_expansion(p, t, m);
    v.right_multiply(m);
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedThermalExpansionMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                             const Point& p,
                                                             const Real t,
                                                             DenseMatrix<Real>& v) const {
    libmesh_error();
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedThermalExpansionMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                           const Point& p,
                                                           const Real t,
                                                           DenseMatrix<Real>& v) const {
    libmesh_error();
}




MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressMatrix::
SectionIntegratedPrestressMatrix(MAST::FieldFunction<DenseMatrix<Real> > &prestress,
                                 MAST::FieldFunction<Real>& h):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedPrestressMatrix2D"),
_prestress(prestress),
_h(h) {
    _functions.insert(prestress.master());
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedPrestressMatrix::operator() (const Point& p,
                                              const Real t,
                                              DenseMatrix<Real>& v) const {
    libmesh_error();
}





void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedPrestressMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                      const Point& p,
                                                      const Real t,
                                                      DenseMatrix<Real>& v) const {
    libmesh_error();
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedPrestressMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                    const Point& p,
                                                    const Real t,
                                                    DenseMatrix<Real>& v) const {
    libmesh_error();
}



//void
//MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressMatrix::prestress_vector(MAST::ElemenetPropertyMatrixType t,
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




