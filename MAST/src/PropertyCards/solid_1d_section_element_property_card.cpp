//
//  solid_1d_section_element_property_card.cpp
//  MAST
//
//  Created by Manav Bhatia on 1/30/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

// MAST includes
#include "PropertyCards/solid_1d_section_element_property_card.h"
#include "StructuralElems/bending_structural_elem.h"


MAST::Solid1DSectionElementPropertyCard::SectionIntegratedExtensionStiffnessMatrix::
SectionIntegratedExtensionStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                          MAST::FieldFunction<Real> * A,
                                          MAST::FieldFunction<Real> * J):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedExtensionStiffnessMatrix1D"),
_material_stiffness(mat),
_A(A),
_J(J) {
    _functions.insert(mat);
    _functions.insert(A);
    _functions.insert(J);
}





void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedExtensionStiffnessMatrix::operator() (const Point& p,
                                                       const Real t,
                                                       DenseMatrix<Real>& m) const {
    // [C]*h
    Real A, J;
    (*_A)(p, t, A);
    (*_J)(p, t, J);
    (*_material_stiffness)(p, t, m);
    m.scale_row(0, A);
    m.scale_row(1, J);
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedExtensionStiffnessMatrix::partial (const MAST::FieldFunctionBase& f,
                                                    const Point& p,
                                                    const Real t,
                                                    DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(2,2);
    Real A, J, dA, dJ;
    (*_A)(p, t, A); _A->partial(f, p, t, dA);
    (*_J)(p, t, J); _J->partial(f, p, t, dJ);
    (*_material_stiffness)(p, t, m); _material_stiffness->partial(f, p, t, dm);
    
    // [C]*dh
    m.scale_row(0, dA);
    m.scale_row(1, dJ);
    
    // += [dC]*h
    dm.scale_row(0, A);
    dm.scale_row(1, J);
    m.add(1., dm);
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedExtensionStiffnessMatrix::total (const MAST::FieldFunctionBase& f,
                                                  const Point& p,
                                                  const Real t,
                                                  DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(2,2);
    Real A, J, dA, dJ;
    (*_A)(p, t, A); _A->total(f, p, t, dA);
    (*_J)(p, t, J); _J->total(f, p, t, dJ);
    (*_material_stiffness)(p, t, m); _material_stiffness->total(f, p, t, dm);
    
    // [C]*dh
    m.scale_row(0, dA);
    m.scale_row(1, dJ);
    
    // += [dC]*h
    dm.scale_row(0, A);
    dm.scale_row(1, J);
    m.add(1., dm);
}






MAST::Solid1DSectionElementPropertyCard::SectionIntegratedExtensionBendingStiffnessMatrix::
SectionIntegratedExtensionBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                                 MAST::FieldFunction<Real> * A_y_moment,
                                                 MAST::FieldFunction<Real> * A_z_moment):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedExtensionBendingStiffnessMatrix1D"),
_material_stiffness(mat),
_A_y_moment(A_y_moment),
_A_z_moment(A_z_moment) {
    _functions.insert(mat);
    _functions.insert(A_y_moment);
    _functions.insert(A_z_moment);
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedExtensionBendingStiffnessMatrix::operator() (const Point& p,
                                                              const Real t,
                                                              DenseMatrix<Real>& m) const {
    Real Ay, Az;
    (*_A_y_moment)(p, t, Ay);
    (*_A_z_moment)(p, t, Az);
    (*_material_stiffness)(p, t, m);
    
    m(0,1)  = m(0,0)*Ay;  // coupling of u and w bending (== theta_y)
    m(0,0) *= Az;        // coupling of u and v bending (== theta_z)
    
    m.scale_row(1, 0); // no coupling for torsion for symmetic sections
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedExtensionBendingStiffnessMatrix::partial (const MAST::FieldFunctionBase& f,
                                                           const Point& p,
                                                           const Real t,
                                                           DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    Real Ay, Az, dAy, dAz;
    (*_A_y_moment)(p, t, Ay); _A_y_moment->partial(f, p, t, dAy);
    (*_A_z_moment)(p, t, Az); _A_z_moment->partial(f, p, t, dAz);
    (*_material_stiffness)(p, t, m); _material_stiffness->partial(f, p, t, dm);
    
    m(0,1)  = m(0,0)*dAy;  // coupling of u and w bending (== theta_y)
    m(0,0) *= dAz;        // coupling of u and v bending (== theta_z)
    m.scale_row(1, 0);     // no coupling for torsion for symmetic sections
    
    dm(0,1)  = dm(0,0)*Ay;
    dm(0,0) *= Az;
    dm.scale_row(1, 0.);
    m.add(1., dm);
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedExtensionBendingStiffnessMatrix::total (const MAST::FieldFunctionBase& f,
                                                         const Point& p,
                                                         const Real t,
                                                         DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    Real Ay, Az, dAy, dAz;
    (*_A_y_moment)(p, t, Ay); _A_y_moment->total(f, p, t, dAy);
    (*_A_z_moment)(p, t, Az); _A_z_moment->total(f, p, t, dAz);
    (*_material_stiffness)(p, t, m); _material_stiffness->total(f, p, t, dm);
    
    m(0,1)  = m(0,0)*dAy;  // coupling of u and w bending (== theta_y)
    m(0,0) *= dAz;        // coupling of u and v bending (== theta_z)
    m.scale_row(1, 0);     // no coupling for torsion for symmetic sections
    
    dm(0,1)  = dm(0,0)*Ay;
    dm(0,0) *= Az;
    dm.scale_row(1, 0.);
    m.add(1., dm);
}




MAST::Solid1DSectionElementPropertyCard::SectionIntegratedBendingStiffnessMatrix::
SectionIntegratedBendingStiffnessMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat,
                                        MAST::FieldFunction<DenseMatrix<Real> > *I):
MAST::FieldFunction<DenseMatrix<Real> > ("SectionIntegratedBendingStiffnessMatrix1D"),
_material_stiffness(mat),
_I(I) {
    _functions.insert(mat);
    _functions.insert(I);
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedBendingStiffnessMatrix::operator() (const Point& p,
                                                     const Real t,
                                                     DenseMatrix<Real>& m) const {
    DenseMatrix<Real> mat;
    (*_I)(p, t, m);
    (*_material_stiffness)(p, t, mat);
    
    // E*I
    m.scale(mat(0,0)); // scale the inertia matrix with modulus of elasticity
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedBendingStiffnessMatrix::partial (const MAST::FieldFunctionBase& f,
                                                  const Point& p,
                                                  const Real t,
                                                  DenseMatrix<Real>& m) const {
    DenseMatrix<Real> mat, dmat, dm;
    (*_I)(p, t, m); _I->partial(f, p, t, dm);
    (*_material_stiffness)(p, t, mat); _material_stiffness->partial(f, p, t, dmat);
    
    // dE*I
    m.scale(dmat(0,0)); // scale the inertia matrix with modulus of elasticity
    
    // E*dI
    m.add(mat(0,0), dm); // scale the inertia matrix with modulus of elasticity
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedBendingStiffnessMatrix::total (const MAST::FieldFunctionBase& f,
                                                const Point& p,
                                                const Real t,
                                                DenseMatrix<Real>& m) const {
    DenseMatrix<Real> mat, dmat, dm;
    (*_I)(p, t, m); _I->total(f, p, t, dm);
    (*_material_stiffness)(p, t, mat); _material_stiffness->total(f, p, t, dmat);
    
    // dE*I
    m.scale(dmat(0,0)); // scale the inertia matrix with modulus of elasticity
    
    // E*dI
    m.add(mat(0,0), dm); // scale the inertia matrix with modulus of elasticity
}





MAST::Solid1DSectionElementPropertyCard::SectionIntegratedInertiaMatrix::
SectionIntegratedInertiaMatrix(MAST::FieldFunction<Real> * rho,
                               MAST::FieldFunction<Real> * A,
                               MAST::FieldFunction<Real> * A_y_moment,
                               MAST::FieldFunction<Real> * A_z_moment,
                               MAST::FieldFunction<Real> * Ip,
                               MAST::FieldFunction<DenseMatrix<Real> >* I):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedInertiaMatrix1D"),
_rho(rho),
_A(A),
_A_y_moment(A_y_moment),
_A_z_moment(A_z_moment),
_Ip(Ip),
_I(I) {
    _functions.insert(rho);
    _functions.insert(A);
    _functions.insert(A_y_moment);
    _functions.insert(A_z_moment);
    _functions.insert(Ip);
    _functions.insert(I);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedInertiaMatrix::operator() (const Point& p,
                                            const Real t,
                                            DenseMatrix<Real>& m) const {
    m.resize(6, 6);
    DenseMatrix<Real> I;
    Real rho, A, Ay, Az, Ip;
    (*_rho)(p, t, rho);
    (*_A)(p, t, A);
    (*_A_y_moment)(p, t, Ay);
    (*_A_z_moment)(p, t, Az);
    (*_Ip)(p, t, Ip);
    (*_I)(p, t, I);
    
    // translation velocities
    m(0,0) = A; m(1,1) = A; m(2,2) = A;
    
    // torsion
    m(3,3) = Ip;
    
    // rotational velocities
    m(0,4) = Ay;  m(4,0) = Ay;   // w-displacement
    m(0,5) = -Az; m(5,0) = -Az;  // v-displacement
    
    // bending rotation inertia
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++)
            m(4+i,4+j) = I(i,j);

    // reduce the rotation inertia component
    for (unsigned int i=0; i<3; i++)
        m(i+3,i+3) *= 1.0e-16;
    
    m.scale(rho);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedInertiaMatrix::partial (const MAST::FieldFunctionBase& f,
                                         const Point& p,
                                         const Real t,
                                         DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(6, 6); dm.resize(6, 6);
    DenseMatrix<Real> I, dI;
    Real rho, A, Ay, Az, Ip, drho, dA, dAy, dAz, dIp;
    (*_rho)(p, t, rho); _rho->partial(f, p, t, drho);
    (*_A)(p, t, A); _A->partial(f, p, t, dA);
    (*_A_y_moment)(p, t, Ay); _A_y_moment->partial(f, p, t, dAy);
    (*_A_z_moment)(p, t, Az); _A_z_moment->partial(f, p, t, dAz);
    (*_Ip)(p, t, Ip); _Ip->partial(f, p, t, dIp);
    (*_I)(p, t, I); _I->partial(f, p, t, dI);
    
    // translation velocities
    m(0,0) = A;  m(1,1) = A;  m(2,2) = A;
    dm(0,0) = A; dm(1,1) = A; dm(2,2) = A;
    
    // torsion
    m(3,3) = Ip;
    dm(3,3) = dIp;
    
    // rotational velocities
    m(0,4) = Ay;  m(4,0) = Ay;   // w-displacement
    dm(0,4) = dAy;  dm(4,0) = dAy;   // w-displacement
    m(0,5) = -Az; m(5,0) = -Az;  // v-displacement
    dm(0,5) = -dAz; m(5,0) = -dAz;  // v-displacement
    
    // bending rotation inertia
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++) {
            m(4+i,4+j) = I(i,j);
            dm(4+i,4+j) = dI(i,j);
        }
    
    // reduce the rotation inertia component
    for (unsigned int i=0; i<3; i++)
        m(i+3,i+3) *= 1.0e-16;

    m.scale(drho);
    m.add(rho, dm);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedInertiaMatrix::total (const MAST::FieldFunctionBase& f,
                                       const Point& p,
                                       const Real t,
                                       DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(6, 6); dm.resize(6, 6);
    DenseMatrix<Real> I, dI;
    Real rho, A, Ay, Az, Ip, drho, dA, dAy, dAz, dIp;
    (*_rho)(p, t, rho); _rho->total(f, p, t, drho);
    (*_A)(p, t, A); _A->total(f, p, t, dA);
    (*_A_y_moment)(p, t, Ay); _A_y_moment->total(f, p, t, dAy);
    (*_A_z_moment)(p, t, Az); _A_z_moment->total(f, p, t, dAz);
    (*_Ip)(p, t, Ip); _Ip->total(f, p, t, dIp);
    (*_I)(p, t, I); _I->total(f, p, t, dI);
    
    // translation velocities
    m(0,0) = A;  m(1,1) = A;  m(2,2) = A;
    dm(0,0) = A; dm(1,1) = A; dm(2,2) = A;
    
    // torsion
    m(3,3) = Ip;
    dm(3,3) = dIp;
    
    // rotational velocities
    m(0,4) = Ay;  m(4,0) = Ay;   // w-displacement
    dm(0,4) = dAy;  dm(4,0) = dAy;   // w-displacement
    m(0,5) = -Az; m(5,0) = -Az;  // v-displacement
    dm(0,5) = -dAz; m(5,0) = -dAz;  // v-displacement
    
    // bending rotation inertia
    for (unsigned int i=0; i<2; i++)
        for (unsigned int j=0; j<2; j++) {
            m(4+i,4+j) = I(i,j);
            dm(4+i,4+j) = dI(i,j);
        }
    
    // reduce the rotation inertia component
    for (unsigned int i=0; i<3; i++)
        m(i+3,i+3) *= 1.0e-16;

    m.scale(drho);
    m.add(rho, dm);
}




MAST::Solid1DSectionElementPropertyCard::SectionIntegratedThermalExpansionAMatrix::
SectionIntegratedThermalExpansionAMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat_stiff,
                                        MAST::FieldFunction<DenseMatrix<Real> > *mat_expansion,
                                        MAST::FieldFunction<Real> * A):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedThermalExpansionAMatrix1D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_A(A) {
    _functions.insert(mat_stiff);
    _functions.insert(mat_expansion);
    _functions.insert(A);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedThermalExpansionAMatrix::operator() (const Point& p,
                                                      const Real t,
                                                      DenseMatrix<Real>& m) const {
    Real A;
    DenseMatrix<Real> at;
    (*_A)(p, t, A);
    (*_material_stiffness)(p, t, m);
    (*_material_expansion)(p, t, at);
    
    m.right_multiply(at);
    m.scale(A);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedThermalExpansionAMatrix::partial (const MAST::FieldFunctionBase& f,
                                                  const Point& p,
                                                  const Real t,
                                                  DenseMatrix<Real>& m) const {
    Real A, dA;
    DenseMatrix<Real> m1, at, dat, dm;
    (*_A)(p, t, A); _A->partial(f, p, t, dA);
    (*_material_stiffness)(p, t, m1); _material_stiffness->partial(f, p, t, dm);
    (*_material_expansion)(p, t, at); _material_expansion->partial(f, p, t, dat);
    
    m=m1;
    
    m.right_multiply(at);
    m.scale(dA);
    
    m1.right_multiply(dat);
    dm.right_multiply(at);
    m1.add(1., dm);
    
    m.add(A, m1);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedThermalExpansionAMatrix::total (const MAST::FieldFunctionBase& f,
                                                const Point& p,
                                                const Real t,
                                                DenseMatrix<Real>& m) const {
    Real A, dA;
    DenseMatrix<Real> m1, at, dat, dm;
    (*_A)(p, t, A); _A->total(f, p, t, dA);
    (*_material_stiffness)(p, t, m1); _material_stiffness->total(f, p, t, dm);
    (*_material_expansion)(p, t, at); _material_expansion->total(f, p, t, dat);
    
    m=m1;
    
    m.right_multiply(at);
    m.scale(dA);
    
    m1.right_multiply(dat);
    dm.right_multiply(at);
    m1.add(1., dm);
    
    m.add(A, m1);
}




MAST::Solid1DSectionElementPropertyCard::SectionIntegratedThermalExpansionBMatrix::
SectionIntegratedThermalExpansionBMatrix(MAST::FieldFunction<DenseMatrix<Real> > *mat_stiff,
                                         MAST::FieldFunction<DenseMatrix<Real> > *mat_expansion,
                                         MAST::FieldFunction<Real> * A_y_moment,
                                         MAST::FieldFunction<Real> * A_z_moment):
MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedThermalExpansionBMatrix1D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion),
_A_y_moment(A_y_moment),
_A_z_moment(A_z_moment) {
    _functions.insert(mat_stiff);
    _functions.insert(mat_expansion);
    _functions.insert(A_y_moment);
    _functions.insert(A_z_moment);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedThermalExpansionBMatrix::operator() (const Point& p,
                                                      const Real t,
                                                      DenseMatrix<Real>& m) const {
    Real Ay, Az;
    DenseMatrix<Real> at;
    (*_A_y_moment)(p, t, Ay);
    (*_A_z_moment)(p, t, Az);
    (*_material_stiffness)(p, t, m);
    (*_material_expansion)(p, t, at);
    
    m.right_multiply(at);
    m(1,0)  = Ay * m(0,0);
    m(0,0) *= Az;
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedThermalExpansionBMatrix::partial (const MAST::FieldFunctionBase& f,
                                                   const Point& p,
                                                   const Real t,
                                                   DenseMatrix<Real>& m) const {
    Real Ay, Az, dAy, dAz;
    DenseMatrix<Real> at, dat, m1, dm;
    (*_A_y_moment)(p, t, Ay); _A_y_moment->partial(f, p, t, dAy);
    (*_A_z_moment)(p, t, Az); _A_z_moment->partial(f, p, t, dAz);
    (*_material_stiffness)(p, t, m1); _material_stiffness->partial(f, p, t, dm);
    (*_material_expansion)(p, t, at); _material_expansion->partial(f, p, t, dat);
    
    m = m1;
    m.right_multiply(at);
    m(1,0)  = dAy * m(0,0);
    m(0,0) *= dAz;
    
    m1.right_multiply(dat);
    dm.right_multiply(at);
    m1.add(1., dm);
    m1(1,0)  = Ay * m1(0,0);
    m1(0,0) *= Az;
    
    m.add(1., m1);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedThermalExpansionBMatrix::total (const MAST::FieldFunctionBase& f,
                                                 const Point& p,
                                                 const Real t,
                                                 DenseMatrix<Real>& m) const {
    Real Ay, Az, dAy, dAz;
    DenseMatrix<Real> at, dat, m1, dm;
    (*_A_y_moment)(p, t, Ay); _A_y_moment->total(f, p, t, dAy);
    (*_A_z_moment)(p, t, Az); _A_z_moment->total(f, p, t, dAz);
    (*_material_stiffness)(p, t, m1); _material_stiffness->total(f, p, t, dm);
    (*_material_expansion)(p, t, at); _material_expansion->total(f, p, t, dat);
    
    m = m1;
    m.right_multiply(at);
    m(1,0)  = dAy * m(0,0);
    m(0,0) *= dAz;
    
    m1.right_multiply(dat);
    dm.right_multiply(at);
    m1.add(1., dm);
    m1(1,0)  = Ay * m1(0,0);
    m1(0,0) *= Az;
    
    m.add(1., m1);
}




MAST::Solid1DSectionElementPropertyCard::SectionIntegratedPrestressAMatrix::
SectionIntegratedPrestressAMatrix(MAST::FieldFunction<DenseMatrix<Real> > *prestress,
                                  MAST::FieldFunction<DenseMatrix<Real> > *T,
                                  MAST::FieldFunction<Real> * A):
MAST::SectionIntegratedPrestressMatrixBase("SectionIntegratedPrestressAMatrix1D"),
_prestress(prestress),
_T(T),
_A(A) {
    _functions.insert(prestress);
    _functions.insert(A);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedPrestressAMatrix::operator() (const Point& p,
                                               const Real t,
                                               DenseMatrix<Real>& m) const {
    DenseMatrix<Real> s, T;
    m.resize(2, 2);
    Real A;
    (*_A)(p, t, A);
    (*_prestress)(p, t, s);
    (*_T)(p, t, T);
    
    // convert the stress to the local coordinate
    s.right_multiply(T);
    s.left_multiply_transpose(T);
    
    m(0,0) = s(0,0)*A; // only sigma_xx is applied, and torsion is neglected
}





void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedPrestressAMatrix::partial (const MAST::FieldFunctionBase& f,
                                            const Point& p,
                                            const Real t,
                                            DenseMatrix<Real>& m) const {
    DenseMatrix<Real> s, ds, T, dT;
    m.resize(2, 2);
    Real A, dA;
    (*_A)(p, t, A); _A->partial(f, p, t, dA);
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
    
    m(0,0) = s(0,0)*dA + ds(0,0)*A; // only sigma_xx is applied, and torsion is neglected
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedPrestressAMatrix::total (const MAST::FieldFunctionBase& f,
                                          const Point& p,
                                          const Real t,
                                          DenseMatrix<Real>& m) const {
    DenseMatrix<Real> s, ds, T, dT;
    m.resize(2, 2);
    Real A, dA;
    (*_A)(p, t, A); _A->total(f, p, t, dA);
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
    
    m(0,0) = s(0,0)*dA + ds(0,0)*A; // only sigma_xx is applied, and torsion is neglected
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedPrestressAMatrix::convert_to_vector(const DenseMatrix<Real> &m,
                                                     DenseVector<Real> &v)  const {
    libmesh_assert_equal_to(m.m(), 2);
    libmesh_assert_equal_to(m.n(), 2);
    v.resize(2);
    v(0) = m(0,0);
}




MAST::Solid1DSectionElementPropertyCard::SectionIntegratedPrestressBMatrix::
SectionIntegratedPrestressBMatrix(MAST::FieldFunction<DenseMatrix<Real> > *prestress,
                                  MAST::FieldFunction<DenseMatrix<Real> > *T,
                                  MAST::FieldFunction<Real> * A_y_moment,
                                  MAST::FieldFunction<Real> * A_z_moment):
MAST::SectionIntegratedPrestressMatrixBase("SectionIntegratedPrestressBMatrix1D"),
_prestress(prestress),
_T(T),
_A_y_moment(A_y_moment),
_A_z_moment(A_z_moment) {
    _functions.insert(prestress);
    _functions.insert(A_y_moment);
    _functions.insert(A_z_moment);
}




void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedPrestressBMatrix::operator() (const Point& p,
                                               const Real t,
                                               DenseMatrix<Real>& m) const {
    DenseMatrix<Real> s, T;
    m.resize(2, 2);
    Real Ay, Az;
    (*_A_y_moment)(p, t, Ay);
    (*_A_z_moment)(p, t, Az);
    (*_prestress)(p, t, s);
    (*_T)(p, t, T);
    
    // convert the stress to the local coordinate
    s.right_multiply(T);
    s.left_multiply_transpose(T);
    
    // only sigma_xx is applied, and torsion is neglected
    m(0,0) =  s(0,0)*Az;
    m(0,1) =  s(0,0)*Ay;
}





void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedPrestressBMatrix::partial (const MAST::FieldFunctionBase& f,
                                            const Point& p,
                                            const Real t,
                                            DenseMatrix<Real>& m) const {
    DenseMatrix<Real> s, ds, T, dT;
    m.resize(2, 2);
    Real Ay, Az, dAy, dAz;
    (*_A_y_moment)(p, t, Ay); _A_y_moment->partial(f, p, t, dAy);
    (*_A_z_moment)(p, t, Az); _A_z_moment->partial(f, p, t, dAz);
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
    
    // only sigma_xx is applied, and torsion is neglected
    m(0,0) =  (s(0,0)*dAz + ds(0,0)*Az);
    m(0,1) =  s(0,0)*dAy + ds(0,0)*Ay;;
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedPrestressBMatrix::total (const MAST::FieldFunctionBase& f,
                                          const Point& p,
                                          const Real t,
                                          DenseMatrix<Real>& m) const {
    DenseMatrix<Real> s, ds, T, dT;
    m.resize(2, 2);
    Real Ay, Az, dAy, dAz;
    (*_A_y_moment)(p, t, Ay); _A_y_moment->total(f, p, t, dAy);
    (*_A_z_moment)(p, t, Az); _A_z_moment->total(f, p, t, dAz);
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
    
    // only sigma_xx is applied, and torsion is neglected
    m(0,0) =  (s(0,0)*dAz + ds(0,0)*Az);
    m(0,1) =  s(0,0)*dAy + ds(0,0)*Ay;
}



void
MAST::Solid1DSectionElementPropertyCard::
SectionIntegratedPrestressBMatrix::convert_to_vector(const DenseMatrix<Real> &m,
                                                     DenseVector<Real> &v)  const {
    libmesh_assert_equal_to(m.m(), 2);
    libmesh_assert_equal_to(m.n(), 2);
    v.resize(2);
    v(0) = m(0,0);
    v(1) = m(0,1);
}




Real
MAST::Solid1DSectionElementPropertyCard::value(const std::string& val) const {
    Point p; // dummy value
    Real v=0.;
    
    if (val == "A") {
        MAST::Solid1DSectionElementPropertyCard::Area
        A(this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
          this->get<MAST::FieldFunction<Real> >("hz").clone().release());
        A(p, 0., v);
    }
    else if (val == "J") {
        MAST::Solid1DSectionElementPropertyCard::TorsionalConstant
        J(this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
          this->get<MAST::FieldFunction<Real> >("hz").clone().release());
        J(p, 0., v);
    }
    else if (val == "IYY" ||
             val == "IZZ" ||
             val == "IYZ") {
        DenseMatrix<Real> Imat;
        MAST::Solid1DSectionElementPropertyCard::AreaInertiaMatrix
        I(this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
          this->get<MAST::FieldFunction<Real> >("hz").clone().release(),
          this->get<MAST::FieldFunction<Real> >("hy_offset").clone().release(),
          this->get<MAST::FieldFunction<Real> >("hz_offset").clone().release());
        I(p, 0., Imat);
        if (val == "IYY")
            v = Imat(0,0);
        else if (val == "IZZ")
            v = Imat(1,1);
        else if (val == "IYZ")
            v = Imat(1,0);
        else
            libmesh_error(); // should not get here
    }
    else
        libmesh_error(); // should not get here
    return v;
}




template <>
std::auto_ptr<MAST::FieldFunction<Real> >
MAST::Solid1DSectionElementPropertyCard::section_property<MAST::FieldFunction<Real> >
(const std::string& val) const {
    std::auto_ptr<MAST::FieldFunction<Real> > rval;
    
    if (val == "A")
        rval.reset( new MAST::Solid1DSectionElementPropertyCard::Area
                   (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                    this->get<MAST::FieldFunction<Real> >("hz").clone().release()));
    else if (val == "J")
        rval.reset(new MAST::Solid1DSectionElementPropertyCard::TorsionalConstant
                   (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                    this->get<MAST::FieldFunction<Real> >("hz").clone().release()));
    else
        libmesh_error(); // should not get here
    
    return rval;
}



std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > >
MAST::Solid1DSectionElementPropertyCard::get_property(MAST::ElemenetPropertyMatrixType t,
                                                      const MAST::StructuralElementBase& e) const {
    
    std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > > rval;
    
    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            rval.reset(new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedExtensionStiffnessMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 1).release(),
                        new MAST::Solid1DSectionElementPropertyCard::Area
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release()),
                        new MAST::Solid1DSectionElementPropertyCard::TorsionalConstant
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release())));
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
            rval.reset(new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedExtensionBendingStiffnessMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 1).release(),
                        new MAST::Solid1DSectionElementPropertyCard::AreaYMoment
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz_offset").clone().release()),
                        new MAST::Solid1DSectionElementPropertyCard::AreaZMoment
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hy_offset").clone().release())));
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX:
            rval.reset(new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedBendingStiffnessMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 1).release(),
                        new MAST::Solid1DSectionElementPropertyCard::AreaInertiaMatrix
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hy_offset").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz_offset").clone().release())));
            break;

        case MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX:
            rval.reset(new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedTransverseStiffnessMatrix
                       (_material->get_property(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX, *this, 1).release(),
                        new MAST::Solid1DSectionElementPropertyCard::Area
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release())));
            break;

        case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
            rval.reset(new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedInertiaMatrix
                       (_material->get<MAST::FieldFunction<Real> >("rho").clone().release(),
                        new MAST::Solid1DSectionElementPropertyCard::Area
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release()),
                        new MAST::Solid1DSectionElementPropertyCard::AreaYMoment
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz_offset").clone().release()),
                        new MAST::Solid1DSectionElementPropertyCard::AreaZMoment
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hy_offset").clone().release()),
                        new MAST::Solid1DSectionElementPropertyCard::PolarInertia
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hy_offset").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz_offset").clone().release()),
                        new MAST::Solid1DSectionElementPropertyCard::AreaInertiaMatrix
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hy_offset").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz_offset").clone().release())));
            break;

        case MAST::SECTION_INTEGRATED_PRESTRESS_A_MATRIX:
            rval.reset(new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedPrestressAMatrix
                       (this->get<MAST::FieldFunction<DenseMatrix<Real> > >
                        ("prestress").clone().release(),
                        e.local_elem().T_matrix_function().release(),
                        new MAST::Solid1DSectionElementPropertyCard::Area
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release())));
            break;
            
        case MAST::SECTION_INTEGRATED_PRESTRESS_B_MATRIX:
            rval.reset(new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedPrestressBMatrix
                       (this->get<MAST::FieldFunction<DenseMatrix<Real> > >
                        ("prestress").clone().release(),
                        e.local_elem().T_matrix_function().release(),
                        new MAST::Solid1DSectionElementPropertyCard::AreaYMoment
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz_offset").clone().release()),
                         new MAST::Solid1DSectionElementPropertyCard::AreaZMoment
                         (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                          this->get<MAST::FieldFunction<Real> >("hz").clone().release(),
                          this->get<MAST::FieldFunction<Real> >("hz_offset").clone().release())));
            break;

        case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX:
            rval.reset(new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedThermalExpansionAMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 1).release(),
                        _material->get_property(MAST::MATERIAL_THERMAL_EXPANSION_MATRIX, *this, 1).release(),
                        new MAST::Solid1DSectionElementPropertyCard::Area
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release())));
            break;

        case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_B_MATRIX:
            rval.reset(new MAST::Solid1DSectionElementPropertyCard::SectionIntegratedThermalExpansionBMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 1).release(),
                        _material->get_property(MAST::MATERIAL_THERMAL_EXPANSION_MATRIX, *this, 1).release(),
                        new MAST::Solid1DSectionElementPropertyCard::AreaYMoment
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz_offset").clone().release()),
                        new MAST::Solid1DSectionElementPropertyCard::AreaZMoment
                        (this->get<MAST::FieldFunction<Real> >("hy").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz").clone().release(),
                         this->get<MAST::FieldFunction<Real> >("hz_offset").clone().release())));
            break;

        case MAST::SECTION_INTEGRATED_MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
    
    return rval;
}

