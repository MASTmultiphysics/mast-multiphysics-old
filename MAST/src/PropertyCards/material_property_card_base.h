//
//  material_property_card_base.h
//  MAST
//
//  Created by Manav Bhatia on 10/15/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_material_property_card_base_h__
#define __MAST_material_property_card_base_h__

// MAST includes
#include "PropertyCards/property_card_base.h"

namespace MAST
{
    enum MaterialPropertyMatrixType {
        MATERIAL_STIFFNESS_MATRIX,
        MATERIAL_DAMPING_MATRIX,
        MATERIAL_INERTIA_MATRIX,
        MATERIAL_THERMAL_EXPANSION_MATRIX,
        MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX
    };

    
    class  MaterialPropertyCardBase: public MAST::PropertyCardBase {
        
    public:
        
        MaterialPropertyCardBase(unsigned int pid):
        MAST::PropertyCardBase (),
        _pid(pid)
        { }
        
        /*!
         *    returns the id for this card
         */
        unsigned int id() const {
            return _pid;
        }

        /*!
         *   calculates the material matrix in \par m of type \par t for a
         *   1D element.
         */
        virtual void calculate_1d_matrix(MAST::MaterialPropertyMatrixType t,
                                         DenseMatrix<Real>& m) const = 0;

        /*!
         *   calculates the material matrix in \par m of type \par t for a
         *   2D element.
         */
        virtual void calculate_2d_matrix(MAST::MaterialPropertyMatrixType t,
                                         DenseMatrix<Real>& m,
                                         bool if_plane_stress) const = 0;

        
        /*!
         *   calculates the material matrix in \par m of type \par t for a 
         *   3D element.
         */
        virtual void calculate_3d_matrix(MAST::MaterialPropertyMatrixType t,
                                         DenseMatrix<Real>& m) const = 0;

        
        /*!
         *   calculates the sensitivity of material matrix in \par m of
         *   type \par t for a 1D element.
         */
        virtual void calculate_1d_matrix_sensitivity(MAST::MaterialPropertyMatrixType t,
                                                     DenseMatrix<Real>& m,
                                                     const MAST::SensitivityParameters& p) const = 0;
        
        /*!
         *   calculates the sensitivity of material matrix in \par m of
         *   type \par t for a 2D element.
         */
        virtual void calculate_2d_matrix_sensitivity(MAST::MaterialPropertyMatrixType t,
                                                     DenseMatrix<Real>& m,
                                                     bool if_plane_stress,
                                                     const MAST::SensitivityParameters& p) const = 0;
        
        
        /*!
         *   calculates the sensitivity of material matrix in \par m of
         *   type \par t for a 3D element.
         */
        virtual void calculate_3d_matrix_sensitivity(MAST::MaterialPropertyMatrixType t,
                                                     DenseMatrix<Real>& m,
                                                     const MAST::SensitivityParameters& p) const = 0;

        
    protected:
        
        /*!
         *   property card id
         */
        unsigned int _pid;
        
    };
    
    
    

    
    class IsotropicMaterialPropertyCard: public MAST::MaterialPropertyCardBase {
        
    public:
        
        IsotropicMaterialPropertyCard(unsigned int pid):
        MAST::MaterialPropertyCardBase (pid)
        { }
        
        /*!
         *   calculates the material matrix in \par m of type \par t for a
         *   1D element.
         */
        virtual void calculate_1d_matrix(MAST::MaterialPropertyMatrixType t,
                                         DenseMatrix<Real>& m) const;

        /*!
         *   calculates the material matrix in \par m of type \par t.
         */
        virtual void calculate_2d_matrix(MAST::MaterialPropertyMatrixType t,
                                         DenseMatrix<Real>& m,
                                         bool if_plane_stress) const;

        /*!
         *   calculates the material matrix in \par m of type \par t.
         */
        virtual void calculate_3d_matrix(MAST::MaterialPropertyMatrixType t,
                                         DenseMatrix<Real>& m) const;
        
        /*!
         *   calculates the sensitivity of material matrix in \par m of
         *   type \par t for a 1D element.
         */
        virtual void calculate_1d_matrix_sensitivity(MAST::MaterialPropertyMatrixType t,
                                                     DenseMatrix<Real>& m,
                                                     const MAST::SensitivityParameters& p) const;
        
        /*!
         *   calculates the sensitivity of material matrix in \par m of
         *   type \par t for a 2D element.
         */
        virtual void calculate_2d_matrix_sensitivity(MAST::MaterialPropertyMatrixType t,
                                                     DenseMatrix<Real>& m,
                                                     bool if_plane_stress,
                                                     const MAST::SensitivityParameters& p) const;
        
        
        /*!
         *   calculates the sensitivity of material matrix in \par m of
         *   type \par t for a 3D element.
         */
        virtual void calculate_3d_matrix_sensitivity(MAST::MaterialPropertyMatrixType t,
                                                     DenseMatrix<Real>& m,
                                                     const MAST::SensitivityParameters& p) const;
        
    protected:
        
    };

}




inline
void
MAST::IsotropicMaterialPropertyCard::calculate_1d_matrix(MAST::MaterialPropertyMatrixType t,
                                                         DenseMatrix<Real> &m) const {
    
    switch (t) {
        case MAST::MATERIAL_STIFFNESS_MATRIX: {
            m.resize(2,2);
            Real E = this->get<Real>("E")(),
            nu = this->get<Real>("nu")(),
            G = E/2./(1.+nu);
            m(0,0) = E;
            m(1,1) = G;
        }
            break;

        case MAST::MATERIAL_INERTIA_MATRIX: {
            m.resize(6,6);
            Real rho = this->get<Real>("rho")();
            for (unsigned int i=0; i<6; i++)
                m(i,i) = rho;
        }
            break;
            
        case MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX: {
            m.resize(2,2);
            Real E = this->get<Real>("E")(),
            nu = this->get<Real>("nu")(),
            kappa = this->get<Real>("kappa")(),
            G = E/2./(1.+nu);
            for (unsigned int i=0; i<2; i++)
                m(i,i) = G*kappa;
        }
            break;

        case MAST::MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
}




inline
void
MAST::IsotropicMaterialPropertyCard::calculate_1d_matrix_sensitivity
(MAST::MaterialPropertyMatrixType t,
 DenseMatrix<Real> &m,
 const MAST::SensitivityParameters& p) const {
    
    // only first order sensitivities are calculated at this point
    libmesh_assert_equal_to(p.total_order(), 1);
    const MAST::FunctionBase& f = p.get_first_order_derivative_parameter();
    
    const bool depends = this->depends_on(f);

    switch (t) {

        case MAST::MATERIAL_STIFFNESS_MATRIX: {
            m.resize(2,2);
            Real E = this->get<Real>("E")(),
            nu = this->get<Real>("nu")();
            if (depends && f.name() == "E") {
                m(0,0) = 1.;
                m(1,1) = 1./2./(1.+nu);
            }
            else if (depends && f.name() == "nu") {
                m(0,0) = 0.;
                m(1,1) = -E/2./pow(1.+nu,2);
            }
            else
                m.zero();
        }
            break;

        case MAST::MATERIAL_INERTIA_MATRIX: {
            m.resize(6,6);
            if (depends && f.name() == "rho")
                for (unsigned int i=0; i<6; i++)
                    m(i,i) = 1.;
            else
                for (unsigned int i=0; i<6; i++)
                    m(i,i) = 0.;
        }
            break;
            
        case MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX: {
            m.resize(2,2);
            Real E = this->get<Real>("E")(),
            nu = this->get<Real>("nu")(),
            kappa = this->get<Real>("kappa")(),
            G = E/2./(1.+nu);
            if (depends && f.name() == "E") {
                for (unsigned int i=0; i<2; i++)
                    m(i,i) = 1./2./(1.+nu)*kappa;
            }
            else if (depends && f.name() == "nu") {
                for (unsigned int i=0; i<2; i++)
                    m(i,i) = -E/2./pow(1.+nu,2)*kappa;
            }
            else if (depends && f.name() == "kappa") {
                for (unsigned int i=0; i<2; i++)
                    m(i,i) = G;
            }
            else
                m.zero();
        }
            break;
            
        case MAST::MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
}




inline
void
MAST::IsotropicMaterialPropertyCard::calculate_2d_matrix(MAST::MaterialPropertyMatrixType t,
                                                         DenseMatrix<Real> &m,
                                                         bool if_plane_stress) const {
    
    switch (t) {
        case MAST::MATERIAL_STIFFNESS_MATRIX: {
            m.resize(3,3);
            Real E = this->get<Real>("E")(),
            nu = this->get<Real>("nu")();
            if (if_plane_stress) {
                for (unsigned int i=0; i<2; i++) {
                    for (unsigned int j=0; j<2; j++)
                        if (i == j) // diagonal: direct stress
                            m(i,i) = E/(1.-nu*nu);
                        else // offdiagonal: direct stress
                            m(i,j) = E*nu/(1.-nu*nu);
                }
                m(2,2) = E/2./(1.+nu); // diagonal: shear stress
            }
            else {
                libmesh_error(); // to be implemented
            }
        }
            break;
            
        case MAST::MATERIAL_INERTIA_MATRIX: {
            m.resize(6,6);
            Real rho = this->get<Real>("rho")();
            for (unsigned int i=0; i<6; i++)
                m(i,i) = rho;
        }
            break;
            
        case MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX: {
            m.resize(2,2);
            Real E = this->get<Real>("E")(),
            nu = this->get<Real>("nu")(),
            kappa = this->get<Real>("kappa")(),
            G = E/2./(1.+nu);
            for (unsigned int i=0; i<2; i++)
                m(i,i) = G*kappa;
        }
            break;
            
            
        case MAST::MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
}





inline
void
MAST::IsotropicMaterialPropertyCard::calculate_2d_matrix_sensitivity
(MAST::MaterialPropertyMatrixType t,
 DenseMatrix<Real> &m,
 bool if_plane_stress,
 const MAST::SensitivityParameters& p) const {
    
    // only first order sensitivities are calculated at this point
    libmesh_assert_equal_to(p.total_order(), 1);
    const MAST::FunctionBase& f = p.get_first_order_derivative_parameter();

    const bool depends = this->depends_on(f);

    switch (t) {
        case MAST::MATERIAL_STIFFNESS_MATRIX: {
            m.resize(3,3);
            Real E = this->get<Real>("E")(),
            nu = this->get<Real>("nu")();
            if (if_plane_stress) {
                if (depends && f.name() == "E") {
                    for (unsigned int i=0; i<2; i++) {
                        for (unsigned int j=0; j<2; j++)
                            if (i == j) // diagonal: direct stress
                                m(i,i) = 1./(1.-nu*nu);
                            else // offdiagonal: direct stress
                                m(i,j) = 1.*nu/(1.-nu*nu);
                    }
                    m(2,2) = 1./2./(1.+nu); // diagonal: shear stress
                }
                else if (depends && f.name() == "nu") {
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
            else {
                libmesh_error(); // to be implemented
            }
        }
            break;
            
        case MAST::MATERIAL_INERTIA_MATRIX: {
            m.resize(6,6);
            if (depends && f.name() == "rho")
                for (unsigned int i=0; i<6; i++)
                    m(i,i) = 1.;
            else
                m.zero();
        }
            break;
            
        case MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX: {
            m.resize(2,2);
            Real E = this->get<Real>("E")(),
            nu = this->get<Real>("nu")(),
            kappa = this->get<Real>("kappa")(),
            G = E/2./(1.+nu);
            if (depends && f.name() == "E")
                for (unsigned int i=0; i<2; i++)
                    m(i,i) = 1./2./(1.+nu)*kappa;
            else if (depends && f.name() == "nu")
                for (unsigned int i=0; i<2; i++)
                    m(i,i) = -E/2./pow(1.+nu,2)*kappa;
            else if (depends && f.name() == "kappa")
                for (unsigned int i=0; i<2; i++)
                    m(i,i) = G;
            else
                m.zero();
        }
            break;
            
            
        case MAST::MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
}


inline
void
MAST::IsotropicMaterialPropertyCard::calculate_3d_matrix(MAST::MaterialPropertyMatrixType t,
                                                         DenseMatrix<Real> &m) const {
    
    switch (t) {
        case MAST::MATERIAL_STIFFNESS_MATRIX: {
            m.resize(6,6);
            Real E = this->get<Real>("E")(),
            nu = this->get<Real>("nu")();
            for (unsigned int i=0; i<3; i++) {
                for (unsigned int j=0; j<3; j++)
                    if (i == j) // diagonal: direct stress
                        m(i,i) = E*(1.-nu)/(1.-nu-2.*nu*nu);
                    else // offdiagonal: direct stress
                        m(i,j) = E*nu/(1.-nu-2.*nu*nu);
                m(i+3,i+3) = E/2./(1.+nu); // diagonal: shear stress
            }
        }
            break;
            
        case MAST::MATERIAL_INERTIA_MATRIX: {
            m.resize(6,6);
            Real rho = this->get<Real>("rho")();
            for (unsigned int i=0; i<6; i++)
                m(i,i) = rho;
        }
            break;

        case MAST::MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
}




inline
void
MAST::IsotropicMaterialPropertyCard::calculate_3d_matrix_sensitivity
(MAST::MaterialPropertyMatrixType t,
 DenseMatrix<Real> &m,
 const MAST::SensitivityParameters& p) const {
    
    // only first order sensitivities are calculated at this point
    libmesh_assert_equal_to(p.total_order(), 1);
    const MAST::FunctionBase& f = p.get_first_order_derivative_parameter();

    const bool depends = this->depends_on(f);

    switch (t) {
        case MAST::MATERIAL_STIFFNESS_MATRIX: {
            m.resize(6,6);
            Real E = this->get<Real>("E")(),
            nu = this->get<Real>("nu")();
            if (depends && f.name() == "E") {
                for (unsigned int i=0; i<3; i++) {
                    for (unsigned int j=0; j<3; j++)
                        if (i == j) // diagonal: direct stress
                            m(i,i) = (1.-nu)/(1.-nu-2.*nu*nu);
                        else // offdiagonal: direct stress
                            m(i,j) = nu/(1.-nu-2.*nu*nu);
                    m(i+3,i+3) = 1./2./(1.+nu); // diagonal: shear stress
                }
            }
            else if (depends && f.name() == "nu") {
                for (unsigned int i=0; i<3; i++) {
                    for (unsigned int j=0; j<3; j++)
                        if (i == j) // diagonal: direct stress
                            m(i,i) = -E/(1.-nu-2.*nu*nu) + E*(1.-nu)/pow(1.-nu-2.*nu*nu,2)*(1.+4.*nu);
                        else // offdiagonal: direct stress
                            m(i,j) =  E/(1.-nu-2.*nu*nu) + E*nu/pow(1.-nu-2.*nu*nu,2)*(1.+4.*nu);
                    m(i+3,i+3) = -E/2./pow(1.+nu,2); // diagonal: shear stress
                }
            }
            else
                m.zero();
        }
            break;
            
        case MAST::MATERIAL_INERTIA_MATRIX: {
            m.resize(6,6);
            if (depends && f.name() == "rho")
                for (unsigned int i=0; i<6; i++)
                    m(i,i) = 1.;
            else
                m.zero();
        }
            break;
            
        case MAST::MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
}


#endif // __MAST_material_property_card_base_h__
