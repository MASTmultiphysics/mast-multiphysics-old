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
        MATERIAL_THERMAL_EXPANSION_MATRIX
    };

    
    class  MaterialPropertyCardBase: public MAST::PropertyCardBase {
        
    public:

        MaterialPropertyCardBase():
        MAST::PropertyCardBase ()
        { }
        
        /*!
         *   calculates the material matrix in \par m of type \par t for a
         *   3D element.
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

        
    protected:
        
    };

    
    class IsotropicMaterialPropertyCard: public MAST::MaterialPropertyCardBase {
        
    public:
        
        IsotropicMaterialPropertyCard():
        MAST::MaterialPropertyCardBase ()
        { }
        
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
        
    protected:
        
    };

}


inline
void
MAST::IsotropicMaterialPropertyCard::calculate_2d_matrix(MAST::MaterialPropertyMatrixType t,
                                                         DenseMatrix<Real> &m,
                                                         bool if_plane_stress) const {
    switch (t) {
        case MAST::MATERIAL_STIFFNESS_MATRIX: {
            m.resize(3,3);
            double E = this->get<Real>("E")(),
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
            double rho = this->get<Real>("rho")();
            for (unsigned int i=0; i<3; i++)
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
MAST::IsotropicMaterialPropertyCard::calculate_3d_matrix(MAST::MaterialPropertyMatrixType t,
                                                         DenseMatrix<Real> &m) const {
    switch (t) {
        case MAST::MATERIAL_STIFFNESS_MATRIX: {
            m.resize(6,6);
            double E = this->get<Real>("E")(),
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
            double rho = this->get<Real>("rho")();
            for (unsigned int i=0; i<3; i++)
                m(i,i) = rho;
        }
            break;
            
        case MAST::MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
}


#endif // __MAST_material_property_card_base_h__
