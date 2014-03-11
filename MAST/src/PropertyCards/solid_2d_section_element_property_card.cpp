//
//  solid_2d_section_element_property_card.cpp
//  MAST
//
//  Created by Manav Bhatia on 1/30/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "PropertyCards/solid_2d_section_element_property_card.h"
#include "StructuralElems/bending_structural_elem.h"



MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionStiffnessMatrix::
SectionIntegratedExtensionStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                          MAST::FieldFunction<Real>* h):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedExtensionStiffnessMatrix2D"),
_material_stiffness(mat),
_h(h) {
    _functions.insert(mat);
    _functions.insert(h);
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedExtensionStiffnessMatrix::operator() (const Point& p,
                                                       const Real t,
                                                       DenseMatrix<Real>& m) const {
    // [C]*h
    Real h;
    (*_h)(p, t, h);
    (*_material_stiffness)(p, t, m);
    m.scale(h);
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedExtensionStiffnessMatrix::partial (const MAST::FieldFunctionBase& f,
                                                    const Point& p,
                                                    const Real t,
                                                    DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    Real h, dhdf;
    (*_h)(p, t, h); _h->partial(f, p, t, dhdf);
    (*_material_stiffness)(p, t, m); _material_stiffness->partial(f, p, t, dm);
    
    // [C]*dh
    m.scale(dhdf);
    
    // += [dC]*h
    m.add(h, dm);
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedExtensionStiffnessMatrix::total (const MAST::FieldFunctionBase& f,
                                                  const Point& p,
                                                  const Real t,
                                                  DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    Real h, dhdf;
    (*_h)(p, t, h); _h->total(f, p, t, dhdf);
    (*_material_stiffness)(p, t, m); _material_stiffness->total(f, p, t, dm);
    
    // [C]*dh
    m.scale(dhdf);
    
    // += [dC]*h
    m.add(h, dm);
}






MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionBendingStiffnessMatrix::
SectionIntegratedExtensionBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                                 MAST::FieldFunction<Real>* h,
                                                 MAST::FieldFunction<Real> *off):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedExtensionBendingStiffnessMatrix2D"),
_material_stiffness(mat),
_h(h),
_off(off) {
    _functions.insert(mat);
    _functions.insert(h);
    _functions.insert(off);
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedExtensionBendingStiffnessMatrix::operator() (const Point& p,
                                                              const Real t,
                                                              DenseMatrix<Real>& m) const {
    // [C]*h
    Real h, off;
    (*_h)(p, t, h);
    (*_off)(p, t, off);
    (*_material_stiffness)(p, t, m);
    m.scale(h*off);
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedExtensionBendingStiffnessMatrix::partial (const MAST::FieldFunctionBase& f,
                                                           const Point& p,
                                                           const Real t,
                                                           DenseMatrix<Real>& m) const {
    
    DenseMatrix<Real> dm;
    m.resize(3,3); dm.resize(3, 3);
    Real h, off, dh, doff;

    (*_h)(p, t, h); _h->partial(f, p, t, dh);
    (*_off)(p, t, off); _off->partial(f, p, t, doff);
    (*_material_stiffness)(p, t, m); _material_stiffness->partial(f, p, t, dm);
    m.scale(dh*off + h*doff);
    m.add(h*off, dm);
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedExtensionBendingStiffnessMatrix::total (const MAST::FieldFunctionBase& f,
                                                         const Point& p,
                                                         const Real t,
                                                         DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(3,3); dm.resize(3, 3);
    Real h, off, dh, doff;
    
    (*_h)(p, t, h); _h->total(f, p, t, dh);
    (*_off)(p, t, off); _off->total(f, p, t, doff);
    (*_material_stiffness)(p, t, m); _material_stiffness->total(f, p, t, dm);
    m.scale(dh*off + h*doff);
    m.add(h*off, dm);
}




MAST::Solid2DSectionElementPropertyCard::SectionIntegratedBendingStiffnessMatrix::
SectionIntegratedBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                        MAST::FieldFunction<Real> *h,
                                        MAST::FieldFunction<Real> *off):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedBendingStiffnessMatrix2D"),
_material_stiffness(mat),
_h(h),
_off(off) {
    _functions.insert(mat);
    _functions.insert(h);
    _functions.insert(off);
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedBendingStiffnessMatrix::operator() (const Point& p,
                                                     const Real t,
                                                     DenseMatrix<Real>& m) const {
    // [C]*h
    Real h, off;
    (*_h)(p, t, h);
    (*_off)(p, t, off);
    (*_material_stiffness)(p, t, m);
    m.scale(pow(h,3)/12. + h*pow(off,2));
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedBendingStiffnessMatrix::partial (const MAST::FieldFunctionBase& f,
                                                  const Point& p,
                                                  const Real t,
                                                  DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(3,3); dm.resize(3, 3);
    Real h, dhdf, off, doff;
    (*_h)(p, t, h); _h->partial(f, p, t, dhdf);
    (*_off)(p, t, off); _h->partial(f, p, t, doff);
    (*_material_stiffness)(p, t, m); _material_stiffness->partial(f, p, t, dm);
    
    // [C]*dh
    m.scale(pow(h,2)/4.*dhdf + dhdf*pow(off,2) + h*2.*off*doff);
    
    // += [dC]*h
    m.add(pow(h,3)/12. + h*pow(off, 2), dm);
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedBendingStiffnessMatrix::total (const MAST::FieldFunctionBase& f,
                                                const Point& p,
                                                const Real t,
                                                DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(3,3); dm.resize(3, 3);
    Real h, dhdf, off, doff;
    (*_h)(p, t, h); _h->total(f, p, t, dhdf);
    (*_off)(p, t, off); _h->total(f, p, t, doff);
    (*_material_stiffness)(p, t, m); _material_stiffness->total(f, p, t, dm);
    
    // [C]*dh
    m.scale(pow(h,2)/4.*dhdf + dhdf*pow(off,2) + h*2.*off*doff);
    
    // += [dC]*h
    m.add(pow(h,3)/12. + h*pow(off, 2), dm);
}





MAST::Solid2DSectionElementPropertyCard::SectionIntegratedInertiaMatrix::
SectionIntegratedInertiaMatrix(MAST::FieldFunction<Real> *rho,
                               MAST::FieldFunction<Real> * h,
                               MAST::FieldFunction<Real> *off):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedInertiaMatrix2D"),
_rho(rho),
_h(h),
_off(off) {
    _functions.insert(rho);
    _functions.insert(h);
    _functions.insert(off);
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedInertiaMatrix::operator() (const Point& p,
                                            const Real t,
                                            DenseMatrix<Real>& m) const {
    m.resize(6, 6);
    Real h, rho, off;
    (*_h)(p, t, h);
    (*_off)(p, t, off);
    (*_rho)(p, t, rho);

    for (unsigned int i=0; i<3; i++)
        m(i,i) = h;
    
    m(0,4) = off*h; m(4,0) = m(0,4);       // extension-bending coupling
    m(1,3) = -off*h; m(3,1) = m(1,3);      // extension-bending coupling
    m(3,3) = pow(h,3)/12. + h*pow(off,2);  // rotary inertia
    m(4,4) = pow(h,3)/12. + h*pow(off,2);  // rotary inertia
    m(5,5) = pow(h,3)/12.*1.0e-12; // neglect the rotary inertia wrt theta_z
    
    // reduce the rotation inertia component
    for (unsigned int i=0; i<2; i++)
        m(i+3,i+3) *= 1.0e-16;

    m.scale(rho);
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedInertiaMatrix::partial (const MAST::FieldFunctionBase& f,
                                         const Point& p,
                                         const Real t,
                                         DenseMatrix<Real>& m) const {
    m.resize(6,6);
    Real h, dhdf, rho, drhodf, off, doff;
    (*_h)(p, t, h); _h->partial(f, p, t, dhdf);
    (*_off)(p, t, off); _off->partial(f, p, t, doff);
    (*_rho)(p, t, rho); _rho->partial(f, p, t, drhodf);
    
    for (unsigned int i=0; i<3; i++)
        m(i,i) = drhodf*h + rho*dhdf;
    
    m(0,4) = doff*h+off*dhdf; m(4,0) = m(0,4);        // extension-bending coupling
    m(1,3) = -doff*h-off*dhdf; m(3,1) = m(1,3);      // extension-bending coupling
    m(3,3) = drhodf*pow(h,3)/12.+rho*pow(h,2)/4.*dhdf;  // rotary inertia
    m(4,4) = drhodf*pow(h,3)/12.+rho*pow(h,2)/4.*dhdf;  // rotary inertia
    m(5,5) = (drhodf*pow(h,3)/12.+rho*pow(h,2)/4.*dhdf)*1.0e-12; // neglect the rotary inertia wrt theta_z
    
    // reduce the rotation inertia component
    for (unsigned int i=0; i<2; i++)
        m(i+3,i+3) *= 1.0e-16;
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedInertiaMatrix::total (const MAST::FieldFunctionBase& f,
                                       const Point& p,
                                       const Real t,
                                       DenseMatrix<Real>& m) const {
    m.resize(6,6);
    Real h, dhdf, rho, drhodf, off, doff;
    (*_h)(p, t, h); _h->total(f, p, t, dhdf);
    (*_off)(p, t, off); _off->total(f, p, t, doff);
    (*_rho)(p, t, rho); _rho->total(f, p, t, drhodf);
    
    for (unsigned int i=0; i<3; i++)
        m(i,i) = drhodf*h + rho*dhdf;
    
    m(0,4) = doff*h+off*dhdf; m(4,0) = m(0,4);        // extension-bending coupling
    m(1,3) = -doff*h-off*dhdf; m(3,1) = m(1,3);      // extension-bending coupling
    m(3,3) = drhodf*pow(h,3)/12.+rho*pow(h,2)/4.*dhdf;  // rotary inertia
    m(4,4) = drhodf*pow(h,3)/12.+rho*pow(h,2)/4.*dhdf;  // rotary inertia
    m(5,5) = (drhodf*pow(h,3)/12.+rho*pow(h,2)/4.*dhdf)*1.0e-12; // neglect the rotary inertia wrt theta_z
    
    // reduce the rotation inertia component
    for (unsigned int i=0; i<2; i++)
        m(i+3,i+3) *= 1.0e-16;
}




MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionAMatrix::
SectionIntegratedThermalExpansionAMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat_stiff,
                                        MAST::FieldFunction<DenseMatrix<Real> > *mat_expansion,
                                        MAST::FieldFunction<Real> * h):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedThermalExpansionAMatrix2D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_h(h) {
    _functions.insert(mat_stiff);
    _functions.insert(mat_expansion);
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedThermalExpansionAMatrix::operator() (const Point& p,
                                                     const Real t,
                                                     DenseMatrix<Real>& m) const {
    DenseMatrix<Real> at;
    Real h;
    (*_h)(p, t, h);
    (*_material_stiffness)(p, t, m);
    (*_material_expansion)(p, t, at);
    
    m.right_multiply(at);
    m.scale(h);
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedThermalExpansionAMatrix::partial (const MAST::FieldFunctionBase& f,
                                                  const Point& p,
                                                  const Real t,
                                                  DenseMatrix<Real>& m) const {
    DenseMatrix<Real> m1, at, dm, dat;
    Real h, dh;
    (*_h)(p, t, h); _h->partial(f, p, t, dh);
    (*_material_stiffness)(p, t, m1); _material_stiffness->partial(f, p, t, dm);
    (*_material_expansion)(p, t, at); _material_expansion->partial(f, p, t, dat);
    
    m = m1;
    
    m.right_multiply(at);
    m.scale(dh);
    
    m1.right_multiply(dat);
    dm.right_multiply(at);
    m1.add(1., dm);
    
    m.add(h, m1);
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedThermalExpansionAMatrix::total (const MAST::FieldFunctionBase& f,
                                                 const Point& p,
                                                 const Real t,
                                                 DenseMatrix<Real>& m) const {
    DenseMatrix<Real> m1, at, dm, dat;
    Real h, dh;
    (*_h)(p, t, h); _h->total(f, p, t, dh);
    (*_material_stiffness)(p, t, m1); _material_stiffness->total(f, p, t, dm);
    (*_material_expansion)(p, t, at); _material_expansion->total(f, p, t, dat);
    
    m = m1;
    
    m.right_multiply(at);
    m.scale(dh);
    
    m1.right_multiply(dat);
    dm.right_multiply(at);
    m1.add(1., dm);
    
    m.add(h, m1);
}




MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionBMatrix::
SectionIntegratedThermalExpansionBMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat_stiff,
                                         MAST::FieldFunction<DenseMatrix<Real> > *mat_expansion,
                                         MAST::FieldFunction<Real> * h,
                                         MAST::FieldFunction<Real> *off):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedThermalExpansionBMatrix2D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_h(h),
_off(off) {
    _functions.insert(mat_stiff);
    _functions.insert(mat_expansion);
    _functions.insert(off);
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedThermalExpansionBMatrix::operator() (const Point& p,
                                                      const Real t,
                                                      DenseMatrix<Real>& m) const {
    DenseMatrix<Real> at;
    Real h, off;
    (*_h)(p, t, h);
    (*_off)(p, t, off);
    (*_material_stiffness)(p, t, m);
    (*_material_expansion)(p, t, at);
    
    m.right_multiply(at);
    m.scale(h*off);
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedThermalExpansionBMatrix::partial (const MAST::FieldFunctionBase& f,
                                                   const Point& p,
                                                   const Real t,
                                                   DenseMatrix<Real>& m) const {
    DenseMatrix<Real> m1, at, dm, dat;
    Real h, dh, off, doff;
    (*_h)(p, t, h); _h->partial(f, p, t, dh);
    (*_off)(p, t, off); _off->partial(f, p, t, doff);
    (*_material_stiffness)(p, t, m1); _material_stiffness->partial(f, p, t, dm);
    (*_material_expansion)(p, t, at); _material_expansion->partial(f, p, t, dat);
    
    m = m1;
    
    m.right_multiply(at);
    m.scale(dh*off+h*doff);
    
    m1.right_multiply(dat);
    dm.right_multiply(at);
    m1.add(1., dm);
    
    m.add(h*off, m1);
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedThermalExpansionBMatrix::total (const MAST::FieldFunctionBase& f,
                                                 const Point& p,
                                                 const Real t,
                                                 DenseMatrix<Real>& m) const {
    DenseMatrix<Real> m1, at, dm, dat;
    Real h, dh, off, doff;
    (*_h)(p, t, h); _h->total(f, p, t, dh);
    (*_off)(p, t, off); _off->total(f, p, t, doff);
    (*_material_stiffness)(p, t, m1); _material_stiffness->total(f, p, t, dm);
    (*_material_expansion)(p, t, at); _material_expansion->total(f, p, t, dat);
    
    m = m1;
    
    m.right_multiply(at);
    m.scale(dh*off+h*doff);
    
    m1.right_multiply(dat);
    dm.right_multiply(at);
    m1.add(1., dm);
    
    m.add(h*off, m1);
}




MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressAMatrix::
SectionIntegratedPrestressAMatrix(MAST::FieldFunction<DenseMatrix<Real> > *prestress,
                                  MAST::FieldFunction<DenseMatrix<Real> > *T,
                                  MAST::FieldFunction<Real> * h):
MAST::SectionIntegratedPrestressMatrixBase("SectionIntegratedPrestressAMatrix2D"),
_prestress(prestress),
_T(T),
_h(h) {
    _functions.insert(prestress);
    _functions.insert(h);
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedPrestressAMatrix::operator() (const Point& p,
                                               const Real t,
                                               DenseMatrix<Real>& m) const {
    DenseMatrix<Real> s, T;
    m.resize(2, 2);
    Real h;
    (*_h)(p, t, h);
    (*_prestress)(p, t, s);
    (*_T)(p, t, T);

    // convert the stress to the local coordinate
    s.right_multiply(T);
    s.left_multiply_transpose(T);
    
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(i,j) = s(i,j)*h;
}





void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedPrestressAMatrix::partial (const MAST::FieldFunctionBase& f,
                                            const Point& p,
                                            const Real t,
                                            DenseMatrix<Real>& m) const {
    DenseMatrix<Real> s, ds, T, dT;
    m.resize(2, 2);
    Real h, dh;
    (*_h)(p, t, h); _h->partial(f, p, t, dh);
    (*_prestress)(p, t, s); _prestress->partial(f, p, t, ds);
    (*_T)(p, t, T); _T->partial(f, p, t, dT);

    // convert the stress to the local coordinate
    s.right_multiply(T);
    s.left_multiply_transpose(T);
    
    // ds =  dT^T s T + T^T s dT + T^T ds T
    DenseMatrix<Real> tmp;
    ds.right_multiply(T);
    ds.left_multiply_transpose(T);
    
    tmp = s;
    tmp.right_multiply(dT);
    tmp.left_multiply_transpose(T);
    ds += tmp;

    tmp = s;
    tmp.right_multiply(T);
    tmp.left_multiply_transpose(dT);
    ds += tmp;

    

    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(i,j) = ds(i,j)*h + s(i,j)*dh;
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedPrestressAMatrix::total (const MAST::FieldFunctionBase& f,
                                         const Point& p,
                                         const Real t,
                                         DenseMatrix<Real>& m) const {
    DenseMatrix<Real> s, ds, T, dT;
    m.resize(2, 2);
    Real h, dh;
    (*_h)(p, t, h); _h->total(f, p, t, dh);
    (*_prestress)(p, t, s); _prestress->total(f, p, t, ds);
    (*_T)(p, t, T); _T->total(f, p, t, dT);
    
    // convert the stress to the local coordinate
    s.right_multiply(T);
    s.left_multiply_transpose(T);
    
    // ds =  dT^T s T + T^T s dT + T^T ds T
    DenseMatrix<Real> tmp;
    ds.right_multiply(T);
    ds.left_multiply_transpose(T);
    
    tmp = s;
    tmp.right_multiply(dT);
    tmp.left_multiply_transpose(T);
    ds += tmp;
    
    tmp = s;
    tmp.right_multiply(T);
    tmp.left_multiply_transpose(dT);
    ds += tmp;
    
    
    
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(i,j) = ds(i,j)*h + s(i,j)*dh;
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedPrestressAMatrix::convert_to_vector(const DenseMatrix<Real> &m,
                                                     DenseVector<Real> &v) const {
    libmesh_assert_equal_to(m.m(), 2);
    libmesh_assert_equal_to(m.n(), 2);
    v.resize(3);
    v(0) = m(0,0);  // sigma x
    v(1) = m(1,1);  // sigma y
    v(2) = m(0,1);  // tau xy
}




MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressBMatrix::
SectionIntegratedPrestressBMatrix(MAST::FieldFunction<DenseMatrix<Real> > *prestress,
                                  MAST::FieldFunction<DenseMatrix<Real> > *T,
                                  MAST::FieldFunction<Real> * h,
                                  MAST::FieldFunction<Real> *off):
MAST::SectionIntegratedPrestressMatrixBase("SectionIntegratedPrestressBMatrix2D"),
_prestress(prestress),
_T(T),
_h(h),
_off(off) {
    _functions.insert(prestress);
    _functions.insert(h);
    _functions.insert(off);
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedPrestressBMatrix::operator() (const Point& p,
                                               const Real t,
                                               DenseMatrix<Real>& m) const {
    DenseMatrix<Real> s, T;
    m.resize(2, 2);
    Real h, off;
    (*_h)(p, t, h);
    (*_off)(p, t, off);
    (*_prestress)(p, t, s);
    (*_T)(p, t, T);
    
    // convert the stress to the local coordinate
    s.right_multiply(T);
    s.left_multiply_transpose(T);
    
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(i,j) = s(i,j)*(h*off);
}





void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedPrestressBMatrix::partial (const MAST::FieldFunctionBase& f,
                                            const Point& p,
                                            const Real t,
                                            DenseMatrix<Real>& m) const {
    DenseMatrix<Real> s, ds, T, dT;
    m.resize(2, 2);
    Real h, dh, off, doff;
    (*_h)(p, t, h); _h->partial(f, p, t, dh);
    (*_off)(p, t, off); _off->partial(f, p, t, doff);
    (*_prestress)(p, t, s); _prestress->partial(f, p, t, ds);
    (*_T)(p, t, T); _T->partial(f, p, t, dT);
    
    // convert the stress to the local coordinate
    s.right_multiply(T);
    s.left_multiply_transpose(T);
    
    // ds =  dT^T s T + T^T s dT + T^T ds T
    DenseMatrix<Real> tmp;
    ds.right_multiply(T);
    ds.left_multiply_transpose(T);
    
    tmp = s;
    tmp.right_multiply(dT);
    tmp.left_multiply_transpose(T);
    ds += tmp;
    
    tmp = s;
    tmp.right_multiply(T);
    tmp.left_multiply_transpose(dT);
    ds += tmp;
    
    
    
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(i,j) = ds(i,j)*(h*off) + s(i,j)*(dh*off+h*doff);
}



void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedPrestressBMatrix::total (const MAST::FieldFunctionBase& f,
                                          const Point& p,
                                          const Real t,
                                          DenseMatrix<Real>& m) const {
    DenseMatrix<Real> s, ds, T, dT;
    m.resize(2, 2);
    Real h, dh, off, doff;
    (*_h)(p, t, h); _h->total(f, p, t, dh);
    (*_off)(p, t, off); _off->total(f, p, t, doff);
    (*_prestress)(p, t, s); _prestress->total(f, p, t, ds);
    (*_T)(p, t, T); _T->total(f, p, t, dT);
    
    // convert the stress to the local coordinate
    s.right_multiply(T);
    s.left_multiply_transpose(T);
    
    // ds =  dT^T s T + T^T s dT + T^T ds T
    DenseMatrix<Real> tmp;
    ds.right_multiply(T);
    ds.left_multiply_transpose(T);
    
    tmp = s;
    tmp.right_multiply(dT);
    tmp.left_multiply_transpose(T);
    ds += tmp;
    
    tmp = s;
    tmp.right_multiply(T);
    tmp.left_multiply_transpose(dT);
    ds += tmp;
    
    
    
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(i,j) = ds(i,j)*(h*off) + s(i,j)*(dh*off+h*doff);
}




void
MAST::Solid2DSectionElementPropertyCard::
SectionIntegratedPrestressBMatrix::convert_to_vector(const DenseMatrix<Real> &m,
                                                     DenseVector<Real> &v) const {
    // nothing to be done for a symmetric section
    v.resize(3);
}


std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > >
MAST::Solid2DSectionElementPropertyCard::get_property(MAST::ElemenetPropertyMatrixType t,
                                                      const MAST::StructuralElementBase& e) const {
    
    std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > > rval;
    
    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            rval.reset(new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionStiffnessMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 2).release(),
                        this->get<MAST::FieldFunction<Real> >("h").clone().release()));
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
            rval.reset(new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedExtensionBendingStiffnessMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 2).release(),
                        this->get<MAST::FieldFunction<Real> >("h").clone().release(),
                        this->get<MAST::FieldFunction<Real> >("off").clone().release()));
            break;

        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX:
            rval.reset(new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedBendingStiffnessMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 2).release(),
                        this->get<MAST::FieldFunction<Real> >("h").clone().release(),
                        this->get<MAST::FieldFunction<Real> >("off").clone().release()));
            break;

        case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
            rval.reset(new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedInertiaMatrix
                       (_material->get<MAST::FieldFunction<Real> >("rho").clone().release(),
                        this->get<MAST::FieldFunction<Real> >("h").clone().release(),
                        this->get<MAST::FieldFunction<Real> >("off").clone().release()));
            break;

        case MAST::SECTION_INTEGRATED_PRESTRESS_A_MATRIX:
            rval.reset(new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressAMatrix
                       (this->get<MAST::FieldFunction<DenseMatrix<Real> > >
                        ("prestress").clone().release(),
                        e.local_elem().T_matrix_function().release(),
                        this->get<MAST::FieldFunction<Real> >("h").clone().release()));
            break;
            
        case MAST::SECTION_INTEGRATED_PRESTRESS_B_MATRIX:
            rval.reset(new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedPrestressBMatrix
                       (this->get<MAST::FieldFunction<DenseMatrix<Real> > >
                        ("prestress").clone().release(),
                        e.local_elem().T_matrix_function().release(),
                        this->get<MAST::FieldFunction<Real> >("h").clone().release(),
                        this->get<MAST::FieldFunction<Real> >("off").clone().release()));
            break;

        case MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX:
            rval.reset(new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedTransverseStiffnessMatrix
                       (_material->get_property(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX, *this, 2).release(),
                        this->get<MAST::FieldFunction<Real> >("h").clone().release()));
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX:
            rval.reset(new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionAMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 2).release(),
                        _material->get_property(MAST::MATERIAL_THERMAL_EXPANSION_MATRIX, *this, 2).release(),
                        this->get<MAST::FieldFunction<Real> >("h").clone().release()));
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_B_MATRIX:
            rval.reset(new MAST::Solid2DSectionElementPropertyCard::SectionIntegratedThermalExpansionBMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 2).release(),
                        _material->get_property(MAST::MATERIAL_THERMAL_EXPANSION_MATRIX, *this, 2).release(),
                        this->get<MAST::FieldFunction<Real> >("h").clone().release(),
                        this->get<MAST::FieldFunction<Real> >("off").clone().release()));
            break;

        case MAST::SECTION_INTEGRATED_MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
    
    return rval;
}



