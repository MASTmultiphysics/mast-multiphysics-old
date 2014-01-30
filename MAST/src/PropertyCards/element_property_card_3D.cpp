//
//  element_property_card_3D.cpp
//  MAST
//
//  Created by Manav Bhatia on 1/30/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "PropertyCards/element_property_card_3D.h"


MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix::
SectionIntegratedStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedStiffnessMatrix3D"),
_material_stiffness(mat) {
    _functions.insert(mat.master());
}


void
MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix::operator() (const Point& p,
                                                                           const Real t,
                                                                           DenseMatrix<Real>& v) const {
    // this only returns the material stiffness
    _material_stiffness(p, t, v);
}
            

void
MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                                                   const Point& p,
                                                                                   const Real t,
                                                                                   DenseMatrix<Real>& v) const {
    // this only returns the material stiffness
    _material_stiffness.partial_derivative(par, p, t, v);
}
            

void
MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                                                 const Point& p,
                                                                                 const Real t,
                                                                                 DenseMatrix<Real>& v) const {
    // this only returns the material stiffness
    _material_stiffness.total_derivative(par, p, t, v);
}
            



MAST::ElementPropertyCard3D::SectionIntegratedInertiaMatrix::
SectionIntegratedInertiaMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedInertiaMatrix3D"),
_material_inertia(mat) {
    _functions.insert(mat.master());
}
            


void
MAST::ElementPropertyCard3D::SectionIntegratedInertiaMatrix::operator() (const Point& p,
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
}



void
MAST::ElementPropertyCard3D::SectionIntegratedInertiaMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                                                 const Point& p,
                                                                                 const Real t,
                                                                                 DenseMatrix<Real>& v) const {
    DenseMatrix<Real> m;
    v.resize(6, 6);
    // this only returns the material inertia
    // sensitivity of rotary inertia is assumed to be zero
    _material_inertia.partial_derivative(par, p, t, m);
    for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<3; j++)
            v(i,j) = m(i,j);
}
            

void
MAST::ElementPropertyCard3D::SectionIntegratedInertiaMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                                               const Point& p,
                                                                               const Real t,
                                                                               DenseMatrix<Real>& v) const {
    DenseMatrix<Real> m;
    v.resize(6, 6);
    // this only returns the material inertia
    // sensitivity of rotary inertia is assumed to be zero
    _material_inertia.total_derivative(par, p, t, m);
    for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<3; j++)
            v(i,j) = m(i,j);
}




MAST::ElementPropertyCard3D::SectionIntegratedThermalExpansionMatrix::
SectionIntegratedThermalExpansionMatrix(MAST::FieldFunction<DenseMatrix<Real> > &mat_stiff,
                                        MAST::FieldFunction<DenseMatrix<Real> > &mat_expansion):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedThermalExpansionMatrix3D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion) {
    _functions.insert(mat_stiff.master());
    _functions.insert(mat_expansion.master());
}
            



void
MAST::ElementPropertyCard3D::SectionIntegratedThermalExpansionMatrix::operator() (const Point& p,
                                                                                  const Real t,
                                                                                  DenseMatrix<Real>& v) const {
    DenseMatrix<Real> m;
    _material_stiffness(p, t, v);
    _material_expansion(p, t, m);
    v.right_multiply(m);
}




void
MAST::ElementPropertyCard3D::SectionIntegratedThermalExpansionMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                                                          const Point& p,
                                                                                          const Real t,
                                                                                          DenseMatrix<Real>& v) const {
    libmesh_error();
}
            



void
MAST::ElementPropertyCard3D::SectionIntegratedThermalExpansionMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                                                        const Point& p,
                                                                                        const Real t,
                                                                                        DenseMatrix<Real>& v) const {
    libmesh_error();
}




MAST::ElementPropertyCard3D::SectionIntegratedPrestressMatrix::
SectionIntegratedPrestressMatrix(MAST::FieldFunction<DenseMatrix<Real> > &prestress):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedPrestressMatrix3D"),
_prestress(prestress){
    _functions.insert(prestress.master());
}
            



void
MAST::ElementPropertyCard3D::SectionIntegratedPrestressMatrix::operator() (const Point& p,
                                                                           const Real t,
                                                                           DenseMatrix<Real>& v) const {
    libmesh_error();
}





void
MAST::ElementPropertyCard3D::SectionIntegratedPrestressMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                                                   const Point& p,
                                                                                   const Real t,
                                                                                   DenseMatrix<Real>& v) const {
    libmesh_error();
}



void
MAST::ElementPropertyCard3D::SectionIntegratedPrestressMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                                                 const Point& p,
                                                                                 const Real t,
                                                                                 DenseMatrix<Real>& v) const {
    libmesh_error();
}



//void
//MAST::ElementPropertyCard3D::SectionIntegratedPrestressMatrix::prestress_vector(MAST::ElemenetPropertyMatrixType t,
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


