//
//  isotropic_material_property_card.cpp
//  MAST
//
//  Created by Manav Bhatia on 1/30/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//


// MAST includes
#include "PropertyCards/isotropic_material_property_card.h"
#include "PropertyCards/element_property_card_2D.h"



MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix1D::StiffnessMatrix1D(MAST::FieldFunction<libMesh::Real>* E,
                                     MAST::FieldFunction<libMesh::Real>* nu ):
MAST::FieldFunction<libMesh::DenseMatrix<libMesh::Real> >("StiffnessMatrix1D"),
_E(E),
_nu(nu)
{
    _functions.insert(E);
    _functions.insert(nu);
}



MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix1D::~StiffnessMatrix1D() {
    delete _E;
    delete _nu;
}


void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix1D::operator() (const libMesh::Point& p,
                               const libMesh::Real t,
                               libMesh::DenseMatrix<libMesh::Real>& m) const {
    m.resize(2,2);
    libMesh::Real E, nu, G;
    (*_E)(p, t, E); (*_nu)(p, t, nu);
    G = E/2./(1.+nu);
    m(0,0) = E;
    m(1,1) = G;
}


void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix1D::partial (const MAST::FieldFunctionBase& f,
                            const libMesh::Point& p,
                            const libMesh::Real t,
                            libMesh::DenseMatrix<libMesh::Real>& m) const {
    libMesh::DenseMatrix<libMesh::Real> dm;
    m.resize(2,2); dm.resize(2,2);
    libMesh::Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E); _E->partial(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->partial(f, p, t, dnudf);
    
    // parM/parE * parE/parf
    dm(0,0) = 1.;
    dm(1,1) = 1./2./(1.+nu);
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    dm(0,0) = 0.;
    dm(1,1) = -E/2./pow(1.+nu,2);
    m.add(dnudf, dm);
}


void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix1D::total (const MAST::FieldFunctionBase& f,
                          const libMesh::Point& p, const libMesh::Real t, libMesh::DenseMatrix<libMesh::Real>& m) const {
    
    
    libMesh::DenseMatrix<libMesh::Real> dm;
    m.resize(2,2); dm.resize(2,2);
    libMesh::Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E); _E->total(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->total(f, p, t, dnudf);
    
    // parM/parE * parE/parf
    dm(0,0) = 1.;
    dm(1,1) = 1./2./(1.+nu);
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    dm(0,0) = 0.;
    dm(1,1) = -E/2./pow(1.+nu,2);
    m.add(dnudf, dm);
}


MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::TransverseShearStiffnessMatrix( MAST::FieldFunction<libMesh::Real> * E,
                                                               MAST::FieldFunction<libMesh::Real> * nu,
                                                               MAST::FieldFunction<libMesh::Real> * kappa):
MAST::FieldFunction<libMesh::DenseMatrix<libMesh::Real> >("TransverseShearStiffnessMatrix"),
_E(E),
_nu(nu),
_kappa(kappa)
{
    _functions.insert(E);
    _functions.insert(nu);
    _functions.insert(kappa);
}



MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::~TransverseShearStiffnessMatrix() {
    delete _E;
    delete _nu;
    delete _kappa;
}


void
MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::operator() (const libMesh::Point& p,
                                            const libMesh::Real t,
                                            libMesh::DenseMatrix<libMesh::Real>& m) const {
    m.resize(2,2);
    libMesh::Real E, nu, kappa, G;
    (*_E)(p, t, E); (*_nu)(p, t, nu); (*_kappa)(p, t, kappa);
    G = E/2./(1.+nu);
    m(0,0) = G*kappa;
    m(1,1) = m(0,0);
}


void
MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::partial (const MAST::FieldFunctionBase& f,
                                         const libMesh::Point& p,
                                         const libMesh::Real t,
                                         libMesh::DenseMatrix<libMesh::Real>& m) const {
    
    libMesh::DenseMatrix<libMesh::Real> dm;
    m.resize(2,2); dm.resize(2, 2);
    libMesh::Real E, nu, kappa, dEdf, dnudf, dkappadf, G;
    (*_E)(p, t, E); _E->partial(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->partial(f, p, t, dnudf);
    (*_kappa)(p, t, kappa); _kappa->partial(f, p, t, dkappadf);
    G = E/2./(1.+nu);
    
    
    // parM/parE * parE/parf
    dm(0,0) = 1./2./(1.+nu)*kappa;
    dm(1,1) = dm(0,0);
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    dm(0,0) = -E/2./pow(1.+nu,2)*kappa;
    dm(1,1) = dm(0,0);
    m.add(dnudf, dm);
    
    // parM/parnu * parkappa/parf
    
    dm(0,0) = G; dm(1,1) = G;
    dm.add(dkappadf, dm);
}




void
MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::total (const MAST::FieldFunctionBase& f,
                                       const libMesh::Point& p,
                                       const libMesh::Real t,
                                       libMesh::DenseMatrix<libMesh::Real>& m) const {
    libMesh::DenseMatrix<libMesh::Real> dm;
    m.resize(2,2); dm.resize(2, 2);
    libMesh::Real E, nu, kappa, dEdf, dnudf, dkappadf, G;
    (*_E)(p, t, E); _E->total(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->total(f, p, t, dnudf);
    (*_kappa)(p, t, kappa); _kappa->total(f, p, t, dkappadf);
    G = E/2./(1.+nu);
    
    
    // parM/parE * parE/parf
    dm(0,0) = 1./2./(1.+nu)*kappa;
    dm(1,1) = dm(0,0);
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    dm(0,0) = -E/2./pow(1.+nu,2)*kappa;
    dm(1,1) = dm(0,0);
    m.add(dnudf, dm);
    
    // parM/parnu * parkappa/parf
    
    dm(0,0) = G; dm(1,1) = G;
    dm.add(dkappadf, dm);
}




MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix2D::StiffnessMatrix2D(MAST::FieldFunction<libMesh::Real> * E,
                                     MAST::FieldFunction<libMesh::Real> * nu ,
                                     bool plane_stress ):
MAST::FieldFunction<libMesh::DenseMatrix<libMesh::Real> >("StiffnessMatrix2D"),
_E(E),
_nu(nu),
_plane_stress(plane_stress)
{
    _functions.insert(E);
    _functions.insert(nu);
}



MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix2D::~StiffnessMatrix2D() {
    delete _E;
    delete _nu;
}



void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix2D::operator() (const libMesh::Point& p,
                               const libMesh::Real t,
                               libMesh::DenseMatrix<libMesh::Real>& m) const {
    libmesh_assert(_plane_stress); // currently only implemented for plane stress
    m.resize(3,3);
    libMesh::Real E, nu;
    (*_E)(p, t, E); (*_nu)(p, t, nu);
    for (unsigned int i=0; i<2; i++) {
        for (unsigned int j=0; j<2; j++)
            if (i == j) // diagonal: direct stress
                m(i,i) = E/(1.-nu*nu);
            else // offdiagonal: direct stress
                m(i,j) = E*nu/(1.-nu*nu);
    }
    m(2,2) = E/2./(1.+nu); // diagonal: shear stress
}



void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix2D::partial (const MAST::FieldFunctionBase& f,
                            const libMesh::Point& p,
                            const libMesh::Real t,
                            libMesh::DenseMatrix<libMesh::Real>& m) const {
    libmesh_assert(_plane_stress); // currently only implemented for plane stress
    libMesh::DenseMatrix<libMesh::Real> dm;
    m.resize(3,3); dm.resize(3, 3);
    libMesh::Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E); _E->partial(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->partial(f, p, t, dnudf);
    
    // parM/parE * parE/parf
    for (unsigned int i=0; i<2; i++) {
        for (unsigned int j=0; j<2; j++)
            if (i == j) // diagonal: direct stress
                dm(i,i) = 1./(1.-nu*nu);
            else // offdiagonal: direct stress
                dm(i,j) = 1.*nu/(1.-nu*nu);
    }
    dm(2,2) = 1./2./(1.+nu); // diagonal: shear stress
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    for (unsigned int i=0; i<2; i++) {
        for (unsigned int j=0; j<2; j++)
            if (i == j) // diagonal: direct stress
                dm(i,i) = E/pow(1.-nu*nu, 2)*2.*nu;
            else // offdiagonal: direct stress
                dm(i,j) = E/(1.-nu*nu) + E*nu/pow(1.-nu*nu,2)*2.*nu;
    }
    dm(2,2) = -E/2./pow(1.+nu,2); // diagonal: shear stress
    m.add(dnudf, dm);
}




void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix2D::total (const MAST::FieldFunctionBase& f,
                          const libMesh::Point& p, const libMesh::Real t, libMesh::DenseMatrix<libMesh::Real>& m) const {
    libmesh_assert(_plane_stress); // currently only implemented for plane stress
    libMesh::DenseMatrix<libMesh::Real> dm;
    m.resize(3,3); dm.resize(3, 3);
    libMesh::Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E); _E->total(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->total(f, p, t, dnudf);
    
    // parM/parE * parE/parf
    for (unsigned int i=0; i<2; i++) {
        for (unsigned int j=0; j<2; j++)
            if (i == j) // diagonal: direct stress
                dm(i,i) = 1./(1.-nu*nu);
            else // offdiagonal: direct stress
                dm(i,j) = 1.*nu/(1.-nu*nu);
    }
    dm(2,2) = 1./2./(1.+nu); // diagonal: shear stress
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    for (unsigned int i=0; i<2; i++) {
        for (unsigned int j=0; j<2; j++)
            if (i == j) // diagonal: direct stress
                dm(i,i) = E/pow(1.-nu*nu, 2)*2.*nu;
            else // offdiagonal: direct stress
                dm(i,j) = E/(1.-nu*nu) + E*nu/pow(1.-nu*nu,2)*2.*nu;
    }
    dm(2,2) = -E/2./pow(1.+nu,2); // diagonal: shear stress
    m.add(dnudf, dm);
}



MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix3D::StiffnessMatrix3D(MAST::FieldFunction<libMesh::Real> * E,
                                     MAST::FieldFunction<libMesh::Real> * nu):
MAST::FieldFunction<libMesh::DenseMatrix<libMesh::Real> >("StiffnessMatrix2D"),
_E(E),
_nu(nu)
{
    _functions.insert(E);
    _functions.insert(nu);
}



MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix3D::~StiffnessMatrix3D() {
    delete _E;
    delete _nu;
}


void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix3D::operator() (const libMesh::Point& p,
                               const libMesh::Real t,
                               libMesh::DenseMatrix<libMesh::Real>& m) const {
    m.resize(3,3);
    libMesh::Real E, nu;
    (*_E)(p, t, E); (*_nu)(p, t, nu);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++)
            if (i == j) // diagonal: direct stress
                m(i,i) = E*(1.-nu)/(1.-nu-2.*nu*nu);
            else // offdiagonal: direct stress
                m(i,j) = E*nu/(1.-nu-2.*nu*nu);
        m(i+3,i+3) = E/2./(1.+nu); // diagonal: shear stress
    }
}



void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix3D::partial (const MAST::FieldFunctionBase& f,
                            const libMesh::Point& p,
                            const libMesh::Real t,
                            libMesh::DenseMatrix<libMesh::Real>& m) const {
    libMesh::DenseMatrix<libMesh::Real> dm;
    m.resize(3,3); dm.resize(3, 3);
    libMesh::Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E); _E->partial(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->partial(f, p, t, dnudf);
    
    // parM/parE * parE/parf
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++)
            if (i == j) // diagonal: direct stress
                dm(i,i) = (1.-nu)/(1.-nu-2.*nu*nu);
            else // offdiagonal: direct stress
                dm(i,j) = nu/(1.-nu-2.*nu*nu);
        dm(i+3,i+3) = 1./2./(1.+nu); // diagonal: shear stress
    }
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++)
            if (i == j) // diagonal: direct stress
                dm(i,i) = -E/(1.-nu-2.*nu*nu) + E*(1.-nu)/pow(1.-nu-2.*nu*nu,2)*(1.+4.*nu);
            else // offdiagonal: direct stress
                dm(i,j) =  E/(1.-nu-2.*nu*nu) + E*nu/pow(1.-nu-2.*nu*nu,2)*(1.+4.*nu);
        dm(i+3,i+3) = -E/2./pow(1.+nu,2); // diagonal: shear stress
    }
    m.add(dnudf, dm);
}




void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix3D::total (const MAST::FieldFunctionBase& f,
                          const libMesh::Point& p,
                          const libMesh::Real t,
                          libMesh::DenseMatrix<libMesh::Real>& m) const {
    libMesh::DenseMatrix<libMesh::Real> dm;
    m.resize(3,3); dm.resize(3, 3);
    libMesh::Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E); _E->total(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->total(f, p, t, dnudf);
    
    // parM/parE * parE/parf
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++)
            if (i == j) // diagonal: direct stress
                dm(i,i) = (1.-nu)/(1.-nu-2.*nu*nu);
            else // offdiagonal: direct stress
                dm(i,j) = nu/(1.-nu-2.*nu*nu);
        dm(i+3,i+3) = 1./2./(1.+nu); // diagonal: shear stress
    }
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++)
            if (i == j) // diagonal: direct stress
                dm(i,i) = -E/(1.-nu-2.*nu*nu) + E*(1.-nu)/pow(1.-nu-2.*nu*nu,2)*(1.+4.*nu);
            else // offdiagonal: direct stress
                dm(i,j) =  E/(1.-nu-2.*nu*nu) + E*nu/pow(1.-nu-2.*nu*nu,2)*(1.+4.*nu);
        dm(i+3,i+3) = -E/2./pow(1.+nu,2); // diagonal: shear stress
    }
    m.add(dnudf, dm);
}



std::auto_ptr<MAST::FieldFunction<libMesh::DenseMatrix<libMesh::Real> > >
MAST::IsotropicMaterialPropertyCard::get_property(MAST::MaterialPropertyMatrixType t,
                                                  const MAST::ElementPropertyCardBase& p,
                                                  const unsigned int dim) const  {
    MAST::FieldFunction<libMesh::DenseMatrix<libMesh::Real> > *rval = NULL;
    
    switch (t) {
        case MAST::MATERIAL_STIFFNESS_MATRIX:
            switch (dim) {
                case 1:
                    rval = new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix1D
                    (this->get<MAST::FieldFunction<libMesh::Real> >("E").clone().release(),
                     this->get<MAST::FieldFunction<libMesh::Real> >("nu").clone().release());
                    break;
                    
                case 2:
                    rval = new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix2D
                    (this->get<MAST::FieldFunction<libMesh::Real> >("E").clone().release(),
                     this->get<MAST::FieldFunction<libMesh::Real> >("nu").clone().release(),
                     dynamic_cast<const MAST::ElementPropertyCard2D&>(p).plane_stress());
                    break;
                    
                case 3:
                    rval = new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix3D
                    (this->get<MAST::FieldFunction<libMesh::Real> >("E").clone().release(),
                     this->get<MAST::FieldFunction<libMesh::Real> >("nu").clone().release());
                    break;
            }
            break;
            
        case MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX:
            rval = new MAST::IsotropicMaterialPropertyCard::TransverseShearStiffnessMatrix
            (this->get<MAST::FieldFunction<libMesh::Real> >("E").clone().release(),
             this->get<MAST::FieldFunction<libMesh::Real> >("nu").clone().release(),
             this->get<MAST::FieldFunction<libMesh::Real> >("kappa").clone().release());
            break;
            
        case MAST::MATERIAL_THERMAL_EXPANSION_MATRIX:
            rval = new MAST::IsotropicMaterialPropertyCard::ThermalExpansionMatrix
            (dim,
             this->get<MAST::FieldFunction<libMesh::Real> >("alpha").clone().release());
            break;
            
        case MAST::MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
    
    // make sure that this is not null
    libmesh_assert(rval);
    
    return std::auto_ptr<MAST::FieldFunction<libMesh::DenseMatrix<libMesh::Real> > >(rval);
}
