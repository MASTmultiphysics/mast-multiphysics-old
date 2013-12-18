//
//  element_property_card_2D.h
//  MAST
//
//  Created by Manav Bhatia on 10/23/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_element_property_card_2D_h__
#define __MAST_element_property_card_2D_h__

// MAST includes
#include "PropertyCards/element_property_card_base.h"
#include "PropertyCards/material_property_card_base.h"


namespace MAST
{
    
    class ElementPropertyCard2D: public MAST::ElementPropertyCardBase {
        
    public:
        ElementPropertyCard2D():
        MAST::ElementPropertyCardBase(),
        _bending_model(MAST::DEFAULT_BENDING),
        _if_plane_stress(true)
        { }
        
        /*!
         *   virtual destructor
         */
        ~ElementPropertyCard2D() { }
        
        /*!
         *   returns the bending model to be used for the 2D element
         */
        void set_bending_model(MAST::BendingOperatorType b)  {
            _bending_model = b;
        }
        
        
        /*!
         *   returns the bending model to be used for the 2D element.
         */
        MAST::BendingOperatorType bending_model(const Elem& elem,
                                         const FEType& fe) const;
        
        
        /*!
         *    returns the extra quadrature order (on top of the system) that
         *    this element should use. This is elevated by two orders for a DKT
         *    element
         */
        virtual unsigned int extra_quadrature_order(const Elem& elem,
                                                    const FEType& fe) const {
            if (this->bending_model(elem, fe) == MAST::DKT)
                return 2;
            else
                return 0;
        }
        
        /*!
         *   sets the flag for plane stress.
         */
        void set_plane_stress(bool val) {
            _if_plane_stress = val;
        }
        
        /*!
         *   returns the flag for plane stress.
         */
        bool plane_stress() const {
            return _if_plane_stress;
        }
        
        
    protected:
        
        /*!
         *    initializes the vector to the prestress in the element
         */
        virtual void _prestress_vector(const DenseMatrix<Real>& T,
                                       DenseVector<Real>& v) const {
            v.resize(3); // zero, if the stress has not been defined
            if (_prestress.m() != 0) {

                DenseMatrix<Real> mat;
                _prestress_matrix(T, mat);
                v(0) = mat(0,0); // sigma_xx
                v(1) = mat(1,1); // sigma_yy
                v(2) = mat(0,1); // sigma_xy
            }
        }
        
        
        /*!
         *    initializes the matrix to the prestress in the element
         */
        virtual void _prestress_matrix(const DenseMatrix<Real>& T,
                                       DenseMatrix<Real>& m) const {
            m.resize(2, 2);
            if (_prestress.m() != 0) {
                DenseMatrix<Real> mat; mat = _prestress;
                mat.right_multiply_transpose(T);
                mat.left_multiply(T);
                
                for (unsigned int i=0; i<2; i++)
                    for (unsigned int j=0; j<2; j++)
                        m(i,j) = mat(i,j);
            }
        }
        

        
        /*!
         *   material property card. By default this chooses DKT for 3 noded
         *   triangles and Mindling for all other elements
         */
        MAST::BendingOperatorType _bending_model;
        
        /*!
         *   if the analysis is plne stress, otherwise it is plane strain.
         *   Note that this is true by default
         */
        bool _if_plane_stress;
        
    };
    
    
    
    
    class Solid2DSectionElementPropertyCard : public MAST::ElementPropertyCard2D {
    public:
        Solid2DSectionElementPropertyCard():
        MAST::ElementPropertyCard2D(),
        _material(NULL)
        { }
        
        /*!
         *   virtual destructor
         */
        ~Solid2DSectionElementPropertyCard() { }
        
        /*!
         *    sets the material card
         */
        void set_material(MAST::MaterialPropertyCardBase& mat) {
            _material = &mat;
        }
        
        
        /*!
         *    returns a reference to the material
         */
        const MAST::MaterialPropertyCardBase& get_material() const {
            libmesh_assert(_material); // make sure it has already been set
            return *_material;
        }
        

        /*!
         *   calculates the matrix in \par m of type \par t.
         */
        virtual void calculate_matrix(const Elem& elem,
                                      MAST::ElemenetPropertyMatrixType t,
                                      DenseMatrix<Real>& m) const;

        /*!
         *   calculates the sensitivity of matrix in \par m of type \par t.
         */
        virtual void calculate_matrix_sensitivity(const Elem& elem,
                                                  MAST::ElemenetPropertyMatrixType t,
                                                  DenseMatrix<Real>& m,
                                                  const MAST::SensitivityParameters& p) const;

        /*!
         *    initializes the vector to the prestress in the element
         */
        virtual void prestress_vector(MAST::ElemenetPropertyMatrixType t,
                                      const DenseMatrix<Real>& T,
                                      DenseVector<Real>& v) const;
        
        
        /*!
         *    initializes the matrix to the prestress in the element
         */
        virtual void prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                      const DenseMatrix<Real>& T,
                                      DenseMatrix<Real>& m) const;
        
        /*!
         *  returns true if the property card depends on the function \p f
         */
        virtual bool depends_on(const FunctionBase& f) const {
            return _material->depends_on(f) ||            // check if the material property depends on the function
            MAST::ElementPropertyCardBase::depends_on(f); // check with this property card
        }
        
        
        /*!
         *  returns true if the property card depends on the functions in \p p
         */
        virtual bool depends_on(const MAST::SensitivityParameters& p) const {
            return _material->depends_on(p) ||            // check if the material property depends on the function
            MAST::ElementPropertyCardBase::depends_on(p); // check with this property card
        }

    protected:
        
        /*!
         *   material property card
         */
        MAST::MaterialPropertyCardBase *_material;
    };
}



inline
MAST::BendingOperatorType
MAST::ElementPropertyCard2D::bending_model(const Elem& elem,
                                           const FEType& fe) const {
    // for a TRI3 element, default bending is DKT. For all other elements
    // the default is Mindlin. Otherwise it returns the model set for
    // this card.
    switch (elem.type()) {
        case TRI3:
            if ((fe.family == LAGRANGE) &&
                (fe.order  == FIRST) &&
                (_bending_model == MAST::DEFAULT_BENDING))
                return MAST::DKT;
            else
                return MAST::TIMOSHENKO;
            break;
            
        default:
            if (_bending_model == MAST::DEFAULT_BENDING)
                return MAST::MINDLIN;
            else
                return _bending_model;
            break;
    }
}



inline void
MAST::Solid2DSectionElementPropertyCard::calculate_matrix(const libMesh::Elem &elem,
                                                          MAST::ElemenetPropertyMatrixType t,
                                                          DenseMatrix<Real>& m) const
{
    libmesh_assert(_material); // should have been set
    
    switch (elem.dim()) {

        case 2: {
            double h = this->get<Real>("h")();
            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
                    _material->calculate_2d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m, _if_plane_stress);
                    m.scale(h);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
                    // for solid sections with isotropic material this is zero
                    m.resize(3,3);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX:
                    _material->calculate_2d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m, _if_plane_stress);
                    m.scale(pow(h,3)/12.);
                    break;

                case MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX:
                    _material->calculate_2d_matrix(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                   m, _if_plane_stress);
                    m.scale(h);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
                    _material->calculate_2d_matrix(MAST::MATERIAL_INERTIA_MATRIX,
                                                   m, _if_plane_stress);
                    m.scale(h);
                    // now scale the rotation dofs with small factors
                    for (unsigned int i=0; i<3; i++) {
                        m(i+3, i+3) *= 1.0e-6;
                    }
                    break;

                default:
                    libmesh_error();
                    break;
            }
        }
            break;
        
        case 1:
        case 3:
        default:
            libmesh_error();
            break;
    }
}




inline void
MAST::Solid2DSectionElementPropertyCard::calculate_matrix_sensitivity(const libMesh::Elem &elem,
                                                                      MAST::ElemenetPropertyMatrixType t,
                                                                      DenseMatrix<Real>& m,
                                                                      const MAST::SensitivityParameters& p) const
{
    libmesh_assert(_material); // should have been set
    
    // only first order sensitivities are calculated at this point
    libmesh_assert_equal_to(p.total_order(), 1);
    
    const MAST::SensitivityParameters::ParameterMap& p_map = p.get_map();
    MAST::SensitivityParameters::ParameterMap::const_iterator it, end;
    it = p_map.begin(); end = p_map.end();

    DenseMatrix<Real> dm;
    const MAST::FunctionBase& f = *(it->first);

    switch (elem.dim()) {
            
        case 2: {
            double h = this->get<Real>("h")();
            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX: {
                    _material->calculate_2d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m, _if_plane_stress);
                    _material->calculate_2d_matrix_sensitivity(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                               dm, _if_plane_stress, p);
                    if (f.name() != "h")
                        m.scale(0.);
                    m.add(h, dm);
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
                    // for solid sections with isotropic material this is zero
                    m.resize(3,3);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX: {
                    _material->calculate_2d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m, _if_plane_stress);
                    _material->calculate_2d_matrix_sensitivity(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                               dm, _if_plane_stress, p);
                    if (f.name() == "h")
                        m.scale(pow(h,2)/4.);
                    else
                        m.scale(0.);
                    m.add(pow(h,3)/12., dm);
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX: {
                    _material->calculate_2d_matrix(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                   m, _if_plane_stress);
                    _material->calculate_2d_matrix_sensitivity(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                               dm, _if_plane_stress, p);
                    if (f.name() != "h")
                        m.scale(0.);
                    m.add(h, dm);
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX: {
                    _material->calculate_2d_matrix(MAST::MATERIAL_INERTIA_MATRIX,
                                                   m, _if_plane_stress);
                    _material->calculate_2d_matrix_sensitivity(MAST::MATERIAL_INERTIA_MATRIX,
                                                               dm, _if_plane_stress, p);
                    if (f.name() != "h")
                        m.scale(0.);
                    m.add(h, dm);
                    // now scale the rotation dofs with small factors
                    for (unsigned int i=0; i<3; i++) {
                        m(i+3, i+3) *= 1.0e-6;
                    }
                }
                    break;
                    
                default:
                    libmesh_error();
                    break;
            }
        }
            break;
            
        case 1:
        case 3:
        default:
            libmesh_error();
            break;
    }
}






inline void
MAST::Solid2DSectionElementPropertyCard::prestress_vector(MAST::ElemenetPropertyMatrixType t,
                                                          const DenseMatrix<Real>& T,
                                                          DenseVector<Real>& v) const {
    double h = this->get<Real>("h")();
    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            _prestress_vector(T, v);
            v.scale(h);
            break;
                    
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
            // for solid sections with isotropic material this is zero
            v.resize(3);
            break;
                    
        default:
            libmesh_error();
            break;
    }
}





inline void
MAST::Solid2DSectionElementPropertyCard::prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                                          const DenseMatrix<Real>& T,
                                                          DenseMatrix<Real>& m) const {
    double h = this->get<Real>("h")();
    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            _prestress_matrix(T, m);
            m.scale(h);
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
            // for solid sections with isotropic material this is zero
            m.resize(2,2);
            break;
            
        default:
            libmesh_error();
            break;
    }
}



#endif // __MAST_element_property_card_2D_h__
