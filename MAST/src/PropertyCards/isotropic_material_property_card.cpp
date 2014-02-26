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
StiffnessMatrix1D::StiffnessMatrix1D(MAST::FieldFunction<Real>* E,
                                     MAST::FieldFunction<Real>* nu ):
MAST::FieldFunction<DenseMatrix<Real> >("StiffnessMatrix1D"),
_E(E),
_nu(nu)
{
    _functions.insert(E);
    _functions.insert(nu);
    
    // function to calculate the matrix
    _m = new std::function<void(const Real, const Real, DenseMatrix<Real>&)>
    ([&](const Real E, const Real nu, DenseMatrix<Real>& m)
     {
         m.zero();
         Real G = E/2./(1.+nu);
         m(0,0) = E;
         m(1,1) = G;
     });
    
    // function to calculate the matrix sensitivity wrt E
    _parm_parE = new std::function<void(const Real, const Real, DenseMatrix<Real>&)>
    ([&](const Real E, const Real nu, DenseMatrix<Real>& m)
     {
         m.zero();
         m(0,0) = 1.;
         m(1,1) = 1./2./(1.+nu);
     });
    
    // function to calculate the matrix sensitivity wrt nu
    _parm_parnu = new std::function<void(const Real, const Real, DenseMatrix<Real>&)>
    ([&](const Real E, const Real nu, DenseMatrix<Real>& m)
     {
         m.zero();
         m(0,0) = 0.;
         m(1,1) = -E/2./pow(1.+nu,2);
     });
}



MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix1D::~StiffnessMatrix1D() {
    delete _E;
    delete _nu;
    delete _m;
    delete _parm_parE;
    delete _parm_parnu;
}


void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix1D::operator() (const Point& p,
                               const Real t,
                               DenseMatrix<Real>& m) const {
    m.resize(2,2);
    Real E, nu;
    (*_E)(p, t, E); (*_nu)(p, t, nu);
    (*_m)(E, nu, m);
}


void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix1D::partial (const MAST::FieldFunctionBase& f,
                            const Point& p,
                            const Real t,
                            DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(2,2); dm.resize(2,2);
    Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E); _E->partial(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->partial(f, p, t, dnudf);
    
    // parM/parE * parE/parf
    (*_parm_parE)(E, nu, dm);
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    (*_parm_parnu)(E, nu, dm);
    m.add(dnudf, dm);
}


void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix1D::total (const MAST::FieldFunctionBase& f,
                          const Point& p, const Real t, DenseMatrix<Real>& m) const {
    
    
    DenseMatrix<Real> dm;
    m.resize(2,2); dm.resize(2,2);
    Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E); _E->total(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->total(f, p, t, dnudf);
    
    // parM/parE * parE/parf
    (*_parm_parE)(E, nu, dm);
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    (*_parm_parnu)(E, nu, dm);
    m.add(dnudf, dm);
}


MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::TransverseShearStiffnessMatrix( MAST::FieldFunction<Real> * E,
                                                               MAST::FieldFunction<Real> * nu,
                                                               MAST::FieldFunction<Real> * kappa):
MAST::FieldFunction<DenseMatrix<Real> >("TransverseShearStiffnessMatrix"),
_E(E),
_nu(nu),
_kappa(kappa)
{
    _functions.insert(E);
    _functions.insert(nu);
    _functions.insert(kappa);
    
    // function to calculate the matrix
    _m = new std::function<void(const Real, const Real, const Real, DenseMatrix<Real>&)>
    ([&](const Real E, const Real nu, const Real kappa, DenseMatrix<Real>& m)
     {
         m.zero();
         Real G = E/2./(1.+nu);
         m(0,0) = G*kappa;
         m(1,1) = m(0,0);
     });
    
    // function to calculate the matrix sensitivity wrt E
    _parm_parE = new std::function<void(const Real, const Real, const Real, DenseMatrix<Real>&)>
    ([&](const Real E, const Real nu, const Real kappa, DenseMatrix<Real>& m)
     {
         m.zero();
         m(0,0) = 1./2./(1.+nu)*kappa;
         m(1,1) = m(0,0);
     });
    
    // function to calculate the matrix sensitivity wrt nu
    _parm_parnu = new std::function<void(const Real, const Real, const Real, DenseMatrix<Real>&)>
    ([&](const Real E, const Real nu, const Real kappa, DenseMatrix<Real>& m)
     {
         m.zero();
         m(0,0) = -E/2./pow(1.+nu,2)*kappa;
         m(1,1) = m(0,0);
     });
    
    // function to calculate the matrix sensitivity wrt kappa
    _parm_parkappa = new std::function<void(const Real, const Real, const Real, DenseMatrix<Real>&)>
    ([&](const Real E, const Real nu, const Real kappa, DenseMatrix<Real>& m)
     {
         m.zero();
         Real G = E/2./(1.+nu);
         m(0,0) = G; m(1,1) = G;
     });
}



MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::~TransverseShearStiffnessMatrix() {
    delete _E;
    delete _nu;
    delete _kappa;
    delete _m;
    delete _parm_parE;
    delete _parm_parnu;
    delete _parm_parkappa;
}


void
MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::operator() (const Point& p,
                                            const Real t,
                                            DenseMatrix<Real>& m) const {
    m.resize(2,2);
    Real E, nu, kappa;
    (*_E)(p, t, E); (*_nu)(p, t, nu); (*_kappa)(p, t, kappa);
    (*_m)(E, nu, kappa, m);
}


void
MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::partial (const MAST::FieldFunctionBase& f,
                                         const Point& p,
                                         const Real t,
                                         DenseMatrix<Real>& m) const {
    
    DenseMatrix<Real> dm;
    m.resize(2,2); dm.resize(2, 2);
    Real E, nu, kappa, dEdf, dnudf, dkappadf;
    (*_E)(p, t, E); _E->partial(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->partial(f, p, t, dnudf);
    (*_kappa)(p, t, kappa); _kappa->partial(f, p, t, dkappadf);
    
    
    
    // parM/parE * parE/parf
    (*_parm_parE)(E, nu, kappa, dm);
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    (*_parm_parnu)(E, nu, kappa, dm);
    m.add(dnudf, dm);
    
    // parM/parnu * parkappa/parf
    (*_parm_parkappa)(E, nu, kappa, dm);
    m.add(dkappadf, dm);
}




void
MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::total (const MAST::FieldFunctionBase& f,
                                       const Point& p,
                                       const Real t,
                                       DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(2,2); dm.resize(2, 2);
    Real E, nu, kappa, dEdf, dnudf, dkappadf;
    (*_E)(p, t, E); _E->total(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->total(f, p, t, dnudf);
    (*_kappa)(p, t, kappa); _kappa->total(f, p, t, dkappadf);
    
    
    
    // parM/parE * parE/parf
    (*_parm_parE)(E, nu, kappa, dm);
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    (*_parm_parnu)(E, nu, kappa, dm);
    m.add(dnudf, dm);
    
    // parM/parnu * parkappa/parf
    (*_parm_parkappa)(E, nu, kappa, dm);
    m.add(dkappadf, dm);
}




MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix2D::StiffnessMatrix2D(MAST::FieldFunction<Real> * E,
                                     MAST::FieldFunction<Real> * nu ,
                                     bool plane_stress ):
MAST::FieldFunction<DenseMatrix<Real> >("StiffnessMatrix2D"),
_E(E),
_nu(nu),
_plane_stress(plane_stress)
{
    _functions.insert(E);
    _functions.insert(nu);
    
    // function to calculate the matrix
    if (plane_stress) {
        _m = new std::function<void(const Real, const Real, DenseMatrix<Real>&)>
        ([&](const Real E, const Real nu, DenseMatrix<Real>& m)
         {
             m.zero();
             for (unsigned int i=0; i<2; i++) {
                 for (unsigned int j=0; j<2; j++)
                     if (i == j) // diagonal: direct stress
                         m(i,i) = E/(1.-nu*nu);
                     else // offdiagonal: direct stress
                         m(i,j) = E*nu/(1.-nu*nu);
             }
             m(2,2) = E/2./(1.+nu); // diagonal: shear stress
         });
        
        // function to calculate the matrix sensitivity wrt E
        _parm_parE = new std::function<void(const Real, const Real, DenseMatrix<Real>&)>
        ([&](const Real E, const Real nu, DenseMatrix<Real>& m)
         {
             m.zero();
             for (unsigned int i=0; i<2; i++) {
                 for (unsigned int j=0; j<2; j++)
                     if (i == j) // diagonal: direct stress
                         m(i,i) = 1./(1.-nu*nu);
                     else // offdiagonal: direct stress
                         m(i,j) = 1.*nu/(1.-nu*nu);
             }
             m(2,2) = 1./2./(1.+nu); // diagonal: shear stress
         });
        
        // function to calculate the matrix sensitivity wrt nu
        _parm_parnu = new std::function<void(const Real, const Real, DenseMatrix<Real>&)>
        ([&](const Real E, const Real nu, DenseMatrix<Real>& m)
         {
             m.zero();
             for (unsigned int i=0; i<2; i++) {
                 for (unsigned int j=0; j<2; j++)
                     if (i == j) // diagonal: direct stress
                         m(i,i) = E/pow(1.-nu*nu, 2)*2.*nu;
                     else // offdiagonal: direct stress
                         m(i,j) = E/(1.-nu*nu) + E*nu/pow(1.-nu*nu,2)*2.*nu;
             }
             m(2,2) = -E/2./pow(1.+nu,2); // diagonal: shear stress
         });
    }
    else
        libmesh_error(); // plane strain not yet implemented
}



MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix2D::~StiffnessMatrix2D() {
    delete _E;
    delete _nu;
    delete _m;
    delete _parm_parE;
    delete _parm_parnu;
}



void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix2D::operator() (const Point& p,
                               const Real t,
                               DenseMatrix<Real>& m) const {
    m.resize(3,3);
    Real E, nu;
    (*_E)(p, t, E); (*_nu)(p, t, nu);
    (*_m)(E, nu, m);
}



void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix2D::partial (const MAST::FieldFunctionBase& f,
                            const Point& p,
                            const Real t,
                            DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(3,3); dm.resize(3, 3);
    Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E); _E->partial(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->partial(f, p, t, dnudf);
    
    // parM/parE * parE/parf
    (*_parm_parE)(E, nu, dm);
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    (*_parm_parnu)(E, nu, dm);
    m.add(dnudf, dm);
}




void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix2D::total (const MAST::FieldFunctionBase& f,
                          const Point& p, const Real t, DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(3,3); dm.resize(3, 3);
    Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E); _E->total(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->total(f, p, t, dnudf);
    
    // parM/parE * parE/parf
    (*_parm_parE)(E, nu, dm);
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    (*_parm_parnu)(E, nu, dm);
    m.add(dnudf, dm);
}



MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix3D::StiffnessMatrix3D(MAST::FieldFunction<Real> * E,
                                     MAST::FieldFunction<Real> * nu):
MAST::FieldFunction<DenseMatrix<Real> >("StiffnessMatrix2D"),
_E(E),
_nu(nu)
{
    _functions.insert(E);
    _functions.insert(nu);
    
    _functions.insert(E);
    _functions.insert(nu);
    
    // function to calculate the matrix
    _m = new std::function<void(const Real, const Real, DenseMatrix<Real>&)>
    ([&](const Real E, const Real nu, DenseMatrix<Real>& m)
     {
         m.zero();
         for (unsigned int i=0; i<3; i++) {
             for (unsigned int j=0; j<3; j++)
                 if (i == j) // diagonal: direct stress
                     m(i,i) = E*(1.-nu)/(1.-nu-2.*nu*nu);
                 else // offdiagonal: direct stress
                     m(i,j) = E*nu/(1.-nu-2.*nu*nu);
             m(i+3,i+3) = E/2./(1.+nu); // diagonal: shear stress
         }
     });
    
    // function to calculate the matrix sensitivity wrt E
    _parm_parE = new std::function<void(const Real, const Real, DenseMatrix<Real>&)>
    ([&](const Real E, const Real nu, DenseMatrix<Real>& m)
     {
         m.zero();
         for (unsigned int i=0; i<3; i++) {
             for (unsigned int j=0; j<3; j++)
                 if (i == j) // diagonal: direct stress
                     m(i,i) = (1.-nu)/(1.-nu-2.*nu*nu);
                 else // offdiagonal: direct stress
                     m(i,j) = nu/(1.-nu-2.*nu*nu);
             m(i+3,i+3) = 1./2./(1.+nu); // diagonal: shear stress
         }
     });
    
    // function to calculate the matrix sensitivity wrt nu
    _parm_parnu = new std::function<void(const Real, const Real, DenseMatrix<Real>&)>
    ([&](const Real E, const Real nu, DenseMatrix<Real>& m)
     {
         m.zero();
         for (unsigned int i=0; i<3; i++) {
             for (unsigned int j=0; j<3; j++)
                 if (i == j) // diagonal: direct stress
                     m(i,i) = -E/(1.-nu-2.*nu*nu) + E*(1.-nu)/pow(1.-nu-2.*nu*nu,2)*(1.+4.*nu);
                 else // offdiagonal: direct stress
                     m(i,j) =  E/(1.-nu-2.*nu*nu) + E*nu/pow(1.-nu-2.*nu*nu,2)*(1.+4.*nu);
             m(i+3,i+3) = -E/2./pow(1.+nu,2); // diagonal: shear stress
         }
     });
}



MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix3D::~StiffnessMatrix3D() {
    delete _E;
    delete _nu;
    delete _m;
    delete _parm_parE;
    delete _parm_parnu;
}


void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix3D::operator() (const Point& p,
                               const Real t,
                               DenseMatrix<Real>& m) const {
    m.resize(3,3);
    Real E, nu;
    (*_E)(p, t, E); (*_nu)(p, t, nu);
    (*_m)(E, nu, m);
}



void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix3D::partial (const MAST::FieldFunctionBase& f,
                            const Point& p,
                            const Real t,
                            DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(3,3); dm.resize(3, 3);
    Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E); _E->partial(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->partial(f, p, t, dnudf);
    
    // parM/parE * parE/parf
    (*_parm_parE)(E, nu, dm);
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    (*_parm_parnu)(E, nu, dm);
    m.add(dnudf, dm);
}




void
MAST::IsotropicMaterialPropertyCard::
StiffnessMatrix3D::total (const MAST::FieldFunctionBase& f,
                          const Point& p,
                          const Real t,
                          DenseMatrix<Real>& m) const {
    DenseMatrix<Real> dm;
    m.resize(3,3); dm.resize(3, 3);
    Real E, nu, dEdf, dnudf;
    (*_E)(p, t, E); _E->total(f, p, t, dEdf);
    (*_nu)(p, t, nu); _nu->total(f, p, t, dnudf);
    
    // parM/parE * parE/parf
    (*_parm_parE)(E, nu, dm);
    m.add(dEdf, dm);
    
    
    // parM/parnu * parnu/parf
    (*_parm_parnu)(E, nu, dm);
    m.add(dnudf, dm);
}



std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > >
MAST::IsotropicMaterialPropertyCard::get_property(MAST::MaterialPropertyMatrixType t,
                                                  const MAST::ElementPropertyCardBase& p,
                                                  const unsigned int dim) const  {
    MAST::FieldFunction<DenseMatrix<Real> > *rval = NULL;
    
    switch (t) {
        case MAST::MATERIAL_STIFFNESS_MATRIX:
            switch (dim) {
                case 1:
                    rval = new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix1D
                    (this->get<MAST::FieldFunction<Real> >("E").clone().release(),
                     this->get<MAST::FieldFunction<Real> >("nu").clone().release());
                    break;
                    
                case 2:
                    rval = new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix2D
                    (this->get<MAST::FieldFunction<Real> >("E").clone().release(),
                     this->get<MAST::FieldFunction<Real> >("nu").clone().release(),
                     dynamic_cast<const MAST::ElementPropertyCard2D&>(p).plane_stress());
                    break;
                    
                case 3:
                    rval = new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix3D
                    (this->get<MAST::FieldFunction<Real> >("E").clone().release(),
                     this->get<MAST::FieldFunction<Real> >("nu").clone().release());
                    break;
            }
            break;
            
        case MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX:
            rval = new MAST::IsotropicMaterialPropertyCard::TransverseShearStiffnessMatrix
            (this->get<MAST::FieldFunction<Real> >("E").clone().release(),
             this->get<MAST::FieldFunction<Real> >("nu").clone().release(),
             this->get<MAST::FieldFunction<Real> >("kappa").clone().release());
            break;
            
        case MAST::MATERIAL_THERMAL_EXPANSION_MATRIX:
            rval = new MAST::IsotropicMaterialPropertyCard::ThermalExpansionMatrix
            (dim,
             this->get<MAST::FieldFunction<Real> >("alpha").clone().release());
            break;
            
        case MAST::MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
    
    // make sure that this is not null
    libmesh_assert(rval);
    
    return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real> > >(rval);
}
