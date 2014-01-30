//
//  isotropic_material_property_card.cpp
//  MAST
//
//  Created by Manav Bhatia on 1/30/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//


// MAST includes
#include "PropertyCards/isotropic_material_property_card.h"




MAST::IsotropicMaterialPropertyCard::StiffnessMatrix1D::StiffnessMatrix1D(MAST::FieldFunction<Real>& E,
                                                                          MAST::FieldFunction<Real>& nu ):
MAST::FieldFunction<DenseMatrix<Real> >("StiffnessMatrix1D"),
_E(E),
_nu(nu)
{
    _functions.insert(E.master());
    _functions.insert(nu.master());
}



void
MAST::IsotropicMaterialPropertyCard::StiffnessMatrix1D::operator() (const Point& p,
                                                                    const Real t,
                                                                    DenseMatrix<Real>& m) const {
    m.resize(2,2);
    Real E, nu, G;
    _E(p, t, E); _nu(p, t, nu);
    G = E/2./(1.+nu);
    m(0,0) = E;
    m(1,1) = G;
}


void
MAST::IsotropicMaterialPropertyCard::StiffnessMatrix1D::partial_derivative (const MAST::SensitivityParameters& par,
                                                                            const Point& p,
                                                                            const Real t,
                                                                            DenseMatrix<Real>& m) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    const MAST::FunctionBase& f = par.get_first_order_derivative_parameter();
    
    m.resize(2,2);
    Real E, nu;
    _E(p, t, E); _nu(p, t, nu);
    if (_E.is_equal(f)) {
        m(0,0) = 1.;
        m(1,1) = 1./2./(1.+nu);
    }
    else if (_nu.is_equal(f)) {
        m(0,0) = 0.;
        m(1,1) = -E/2./pow(1.+nu,2);
    }
}


void
MAST::IsotropicMaterialPropertyCard::StiffnessMatrix1D::total_derivative (const MAST::SensitivityParameters& par,
                                                                          const Point& p, const Real t, DenseMatrix<Real>& m) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    
    DenseMatrix<Real> dm;
    m.resize(2,2); dm.resize(2,2);
    Real df;
    
    // total derivative of E wrt parameter
    _E.total_derivative(par, p, t, df);
    if (fabs(df) > 0.) {
        MAST::SensitivityParameters dE; dE.add_parameter(&_E, 1);
        this->partial_derivative(dE, p, t, dm);
        m.add(df, dm);
    }
    
    // total derivative of nu wrt parameter
    _nu.total_derivative(par, p, t, df);
    if (fabs(df) > 0.) {
        MAST::SensitivityParameters dnu; dnu.add_parameter(&_nu, 1);
        this->partial_derivative(dnu, p, t, dm);
        m.add(df, dm);
    }
}


MAST::IsotropicMaterialPropertyCard::TransverseShearStiffnessMatrix::TransverseShearStiffnessMatrix( MAST::FieldFunction<Real>& E,
                                                                                                    MAST::FieldFunction<Real>& nu,
                                                                                                    MAST::FieldFunction<Real>& kappa):
MAST::FieldFunction<DenseMatrix<Real> >("TransverseShearStiffnessMatrix"),
_E(E),
_nu(nu),
_kappa(kappa)
{
    _functions.insert(E.master());
    _functions.insert(nu.master());
    _functions.insert(kappa.master());
}




void
MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::operator() (const Point& p,
                                            const Real t,
                                            DenseMatrix<Real>& m) const {
    m.resize(2,2);
    Real E, nu, kappa, G;
    _E(p, t, E); _nu(p, t, nu); _kappa(p, t, kappa);
    G = E/2./(1.+nu);
    for (unsigned int i=0; i<2; i++)
        m(i,i) = G*kappa;
}


void
MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                    const Point& p,
                                                    const Real t,
                                                    DenseMatrix<Real>& m) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    const MAST::FunctionBase& f = par.get_first_order_derivative_parameter();
    
    m.resize(2,2);
    Real E, nu, G, kappa;
    _E(p, t, E); _nu(p, t, nu); _kappa(p, t, kappa);
    G = E/2./(1.+nu);
    if (_E.is_equal(f)) {
        for (unsigned int i=0; i<2; i++)
            m(i,i) = 1./2./(1.+nu)*kappa;
    }
    else if (_nu.is_equal(f)) {
        for (unsigned int i=0; i<2; i++)
            m(i,i) = -E/2./pow(1.+nu,2)*kappa;
    }
    else if (_kappa.is_equal(f)) {
        for (unsigned int i=0; i<2; i++)
            m(i,i) = G;
    }
}




void
MAST::IsotropicMaterialPropertyCard::
TransverseShearStiffnessMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                  const Point& p,
                                                  const Real t,
                                                  DenseMatrix<Real>& m) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    
    DenseMatrix<Real> dm;
    m.resize(2,2); dm.resize(2,2);
    Real df;
    
    // total derivative of E wrt parameter
    _E.total_derivative(par, p, t, df);
    if (fabs(df) > 0.) {
        MAST::SensitivityParameters dE; dE.add_parameter(&_E, 1);
        this->partial_derivative(dE, p, t, dm);
        m.add(df, dm);
    }
    
    // total derivative of nu wrt parameter
    _nu.total_derivative(par, p, t, df);
    if (fabs(df) > 0.) {
        MAST::SensitivityParameters dnu; dnu.add_parameter(&_nu, 1);
        this->partial_derivative(dnu, p, t, dm);
        m.add(df, dm);
    }
    
    // total derivative of kappa wrt parameter
    _kappa.total_derivative(par, p, t, df);
    if (fabs(df) > 0.) {
        MAST::SensitivityParameters dkappa; dkappa.add_parameter(&_kappa, 1);
        this->partial_derivative(dkappa, p, t, dm);
        m.add(df, dm);
    }
}


MAST::IsotropicMaterialPropertyCard::InertiaMatrix::InertiaMatrix( MAST::FieldFunction<Real>& rho):
MAST::FieldFunction<DenseMatrix<Real> >("InertiaMatrix"),
_rho(rho)
{
    _functions.insert(rho.master());
}



void
MAST::IsotropicMaterialPropertyCard::InertiaMatrix::operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const {
    m.resize(3,3);
    Real rho;
    _rho(p, t, rho);
    for (unsigned int i=0; i<3; i++) // for u, v, w
        m(i,i) = rho;
}



void
MAST::IsotropicMaterialPropertyCard::InertiaMatrix::partial_derivative (const MAST::SensitivityParameters& par,
                                                                        const Point& p, const Real t, DenseMatrix<Real>& m) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    const MAST::FunctionBase& f = par.get_first_order_derivative_parameter();
    
    m.resize(3,3);
    if (_rho.is_equal(f)) {
        for (unsigned int i=0; i<3; i++)
            m(i,i) = 1.;
    }
}



void
MAST::IsotropicMaterialPropertyCard::InertiaMatrix::total_derivative (const MAST::SensitivityParameters& par,
                                                                      const Point& p, const Real t, DenseMatrix<Real>& m) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    
    DenseMatrix<Real> dm;
    m.resize(6,6); dm.resize(6,6);
    Real df;
    
    // total derivative of rho wrt parameter
    _rho.total_derivative(par, p, t, df);
    if (fabs(df) > 0.) {
        MAST::SensitivityParameters drho; drho.add_parameter(&_rho, 1);
        this->partial_derivative(drho, p, t, dm);
        m.add(df, dm);
    }
}



MAST::IsotropicMaterialPropertyCard::StiffnessMatrix2D::StiffnessMatrix2D(MAST::FieldFunction<Real>& E,
                                                                          MAST::FieldFunction<Real>& nu ,
                                                                          bool plane_stress ):
MAST::FieldFunction<DenseMatrix<Real> >("StiffnessMatrix2D"),
_E(E),
_nu(nu),
_plane_stress(plane_stress)
{
    _functions.insert(E.master());
    _functions.insert(nu.master());
}




void
MAST::IsotropicMaterialPropertyCard::StiffnessMatrix2D::set_plane_stress(bool val) {
    _plane_stress = val;
}



void
MAST::IsotropicMaterialPropertyCard::StiffnessMatrix2D::operator() (const Point& p,
                                                                    const Real t,
                                                                    DenseMatrix<Real>& m) const {
    m.resize(3,3);
    Real E, nu, G;
    _E(p, t, E); _nu(p, t, nu);
    G = E/2./(1.+nu);
    
    if (_plane_stress) {
        for (unsigned int i=0; i<2; i++) {
            for (unsigned int j=0; j<2; j++)
                if (i == j) // diagonal: direct stress
                    m(i,i) = E/(1.-nu*nu);
                else // offdiagonal: direct stress
                    m(i,j) = E*nu/(1.-nu*nu);
        }
        m(2,2) = E/2./(1.+nu); // diagonal: shear stress
    }
    else
        libmesh_error(); // not implemented for plane strain yet
}



void
MAST::IsotropicMaterialPropertyCard::StiffnessMatrix2D::partial_derivative (const MAST::SensitivityParameters& par,
                                                                            const Point& p,
                                                                            const Real t,
                                                                            DenseMatrix<Real>& m) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    const MAST::FunctionBase& f = par.get_first_order_derivative_parameter();
    
    m.resize(3,3);
    Real E, nu;
    _E(p, t, E); _nu(p, t, nu);
    
    if (_plane_stress) {
        if (_E.is_equal(f)) {
            for (unsigned int i=0; i<2; i++) {
                for (unsigned int j=0; j<2; j++)
                    if (i == j) // diagonal: direct stress
                        m(i,i) = 1./(1.-nu*nu);
                    else // offdiagonal: direct stress
                        m(i,j) = 1.*nu/(1.-nu*nu);
            }
            m(2,2) = 1./2./(1.+nu); // diagonal: shear stress
        }
        else if (_nu.is_equal(f)) {
            for (unsigned int i=0; i<2; i++) {
                for (unsigned int j=0; j<2; j++)
                    if (i == j) // diagonal: direct stress
                        m(i,i) = E/pow(1.-nu*nu, 2)*2.*nu;
                    else // offdiagonal: direct stress
                        m(i,j) = E/(1.-nu*nu) + E*nu/pow(1.-nu*nu,2)*2.*nu;
            }
            m(2,2) = -E/2./pow(1.+nu,2); // diagonal: shear stress
        }
        else
            m.zero();
    }
    else
        libmesh_error(); // to be implemented
}




void
MAST::IsotropicMaterialPropertyCard::StiffnessMatrix2D::total_derivative (const MAST::SensitivityParameters& par,
                                                                          const Point& p, const Real t, DenseMatrix<Real>& m) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    
    DenseMatrix<Real> dm;
    m.resize(2,2); dm.resize(2,2);
    Real df;
    
    // total derivative of E wrt parameter
    _E.total_derivative(par, p, t, df);
    if (fabs(df) > 0.) {
        MAST::SensitivityParameters dE; dE.add_parameter(&_E, 1);
        this->partial_derivative(dE, p, t, dm);
        m.add(df, dm);
    }
    
    // total derivative of nu wrt parameter
    _nu.total_derivative(par, p, t, df);
    if (fabs(df) > 0.) {
        MAST::SensitivityParameters dnu; dnu.add_parameter(&_nu, 1);
        this->partial_derivative(dnu, p, t, dm);
        m.add(df, dm);
    }
}



MAST::IsotropicMaterialPropertyCard::StiffnessMatrix3D::StiffnessMatrix3D(MAST::FieldFunction<Real>& E,
                                                                          MAST::FieldFunction<Real>& nu):
MAST::FieldFunction<DenseMatrix<Real> >("StiffnessMatrix2D"),
_E(E),
_nu(nu)
{
    _functions.insert(E.master());
    _functions.insert(nu.master());
}


void
MAST::IsotropicMaterialPropertyCard::StiffnessMatrix3D::operator() (const Point& p,
                                                                    const Real t,
                                                                    DenseMatrix<Real>& m) const {
    m.resize(3,3);
    Real E, nu;
    _E(p, t, E); _nu(p, t, nu);
    
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
MAST::IsotropicMaterialPropertyCard::StiffnessMatrix3D::partial_derivative (const MAST::SensitivityParameters& par,
                                                                            const Point& p,
                                                                            const Real t,
                                                                            DenseMatrix<Real>& m) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    const MAST::FunctionBase& f = par.get_first_order_derivative_parameter();
    
    m.resize(3,3);
    Real E, nu;
    _E(p, t, E); _nu(p, t, nu);
    
    if (_E.is_equal(f)) {
        for (unsigned int i=0; i<3; i++) {
            for (unsigned int j=0; j<3; j++)
                if (i == j) // diagonal: direct stress
                    m(i,i) = (1.-nu)/(1.-nu-2.*nu*nu);
                else // offdiagonal: direct stress
                    m(i,j) = nu/(1.-nu-2.*nu*nu);
            m(i+3,i+3) = 1./2./(1.+nu); // diagonal: shear stress
        }
    }
    else if (_nu.is_equal(f)) {
        for (unsigned int i=0; i<3; i++) {
            for (unsigned int j=0; j<3; j++)
                if (i == j) // diagonal: direct stress
                    m(i,i) = -E/(1.-nu-2.*nu*nu) + E*(1.-nu)/pow(1.-nu-2.*nu*nu,2)*(1.+4.*nu);
                else // offdiagonal: direct stress
                    m(i,j) =  E/(1.-nu-2.*nu*nu) + E*nu/pow(1.-nu-2.*nu*nu,2)*(1.+4.*nu);
            m(i+3,i+3) = -E/2./pow(1.+nu,2); // diagonal: shear stress
        }
    }
}




void
MAST::IsotropicMaterialPropertyCard::StiffnessMatrix3D::total_derivative (const MAST::SensitivityParameters& par,
                                                                          const Point& p,
                                                                          const Real t,
                                                                          DenseMatrix<Real>& m) const {
    libmesh_assert_equal_to(par.total_order(), 1);
    
    DenseMatrix<Real> dm;
    m.resize(2,2); dm.resize(2,2);
    Real df;
    
    // total derivative of E wrt parameter
    _E.total_derivative(par, p, t, df);
    if (fabs(df) > 0.) {
        MAST::SensitivityParameters dE; dE.add_parameter(&_E, 1);
        this->partial_derivative(dE, p, t, dm);
        m.add(df, dm);
    }
    
    // total derivative of nu wrt parameter
    _nu.total_derivative(par, p, t, df);
    if (fabs(df) > 0.) {
        MAST::SensitivityParameters dnu; dnu.add_parameter(&_nu, 1);
        this->partial_derivative(dnu, p, t, dm);
        m.add(df, dm);
    }
}



std::auto_ptr<MAST::FunctionBase>
MAST::IsotropicMaterialPropertyCard::get_property(MAST::MaterialPropertyMatrixType t,
                                                  const unsigned int dim)  {
    std::auto_ptr<MAST::FunctionBase> rval;
    
    switch (t) {
        case MAST::MATERIAL_STIFFNESS_MATRIX:
            switch (dim) {
                case 1:
                    rval.reset
                    (new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix1D
                     ( this->get<MAST::FieldFunction<Real> >("E"),
                      this->get<MAST::FieldFunction<Real> >("nu") ));
                    break;
                    
                case 2:
                    rval.reset
                    (new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix2D
                     ( this->get<MAST::FieldFunction<Real> >("E"),
                      this->get<MAST::FieldFunction<Real> >("nu") ));
                    break;
                    
                case 3:
                    rval.reset
                    (new MAST::IsotropicMaterialPropertyCard::StiffnessMatrix3D
                     ( this->get<MAST::FieldFunction<Real> >("E"),
                      this->get<MAST::FieldFunction<Real> >("nu") ));
                    break;
                    
            }
            break;
            
        case MAST::MATERIAL_INERTIA_MATRIX:
            rval.reset
            (new MAST::IsotropicMaterialPropertyCard::InertiaMatrix
             ( this->get<MAST::FieldFunction<Real> >("rho")));
            break;
            
            
        case MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX:
            rval.reset
            (new MAST::IsotropicMaterialPropertyCard::TransverseShearStiffnessMatrix
             ( this->get<MAST::FieldFunction<Real> >("E"),
              this->get<MAST::FieldFunction<Real> >("nu"),
              this->get<MAST::FieldFunction<Real> >("kappa")));
            break;
            
        case MAST::MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
}
