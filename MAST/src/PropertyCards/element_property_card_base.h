//
//  element_property_card_base.h
//  MAST
//
//  Created by Manav Bhatia on 10/15/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_element_property_card_base_h__
#define __MAST_element_property_card_base_h__

// MAST includes
#include "PropertyCards/property_card_base.h"

// libMesh includes
#include "libmesh/elem.h"


namespace MAST
{
    enum ElemenetPropertyMatrixType {
        MATERIAL_STIFFNESS_MATRIX,
        MATERIAL_DAMPING_MATRIX,
        MATERIAL_INERTIA_MATRIX,
        SECTION_INTEGRATED_MATERIAL_STIFFNESS_MATRIX,
        SECTION_INTEGRATED_MATERIAL_DAMPING_MATRIX,
        SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX
    };
    
    
    class ElementPropertyCardBase: public MAST::PropertyCardBase {
        
    public:
        ElementPropertyCardBase():
        MAST::PropertyCardBase()
        { }
        
        
        /*!
         *   calculates the material matrix in \par m of type \par t
         */
        virtual void calculate_matrix(const Elem& elem,
                                      MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const;
        
        
    protected:
        
        /*!
         *    pointer to the material property card
         */
        MAST::MaterialPropertyCardBase* _material;
        
    };
    
    
    inline void
    MAST::ElementPropertyCardBase::calculate_matrix(const libMesh::Elem &elem,
                                                    MAST::ElemenetPropertyMatrixType t,
                                                    DenseMatrix<Real>& m) const
    {
        switch (elem.dim()) {
            case 1:
                switch (t) {
                    case MAST::MATERIAL_STIFFNESS_MATRIX:
                    case MAST::MATERIAL_DAMPING_MATRIX:
                    case MAST::MATERIAL_INERTIA_MATRIX:
                    case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_MATRIX:
                    case MAST::SECTION_INTEGRATED_MATERIAL_DAMPING_MATRIX:
                    case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
                    default:
                        libmesh_error();
                        break;
                }
                break;
                
            case 2:
                switch (t) {
                    case MAST::MATERIAL_STIFFNESS_MATRIX:
                    case MAST::MATERIAL_DAMPING_MATRIX:
                    case MAST::MATERIAL_INERTIA_MATRIX:
                    case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_MATRIX:
                    case MAST::SECTION_INTEGRATED_MATERIAL_DAMPING_MATRIX:
                    case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
                    default:
                        libmesh_error();
                        break;
                }
                break;
                
            case 3:
                switch (t) {
                    case MAST::MATERIAL_STIFFNESS_MATRIX: {
                        m.resize(6,6);
                        double E = this->get<Real>("E")(),
                        nu = this->get<Real>("nu")();
                        for (unsigned int i=0; i<3; i++) {
                            for (unsigned int j=0; i<3; j++)
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
                    case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_MATRIX:
                    case MAST::SECTION_INTEGRATED_MATERIAL_DAMPING_MATRIX:
                    case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
                    default:
                        libmesh_error();
                        break;
                }
                break;
                
            default:
                libmesh_error();
                break;
        }
    }
    
}



#endif // __MAST_element_property_card_base_h__
