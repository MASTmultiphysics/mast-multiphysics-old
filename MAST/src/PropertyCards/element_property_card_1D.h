//
//  element_property_card_1D.h
//  MAST
//
//  Created by Manav Bhatia on 10/25/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_element_property_card_1D_h__
#define __MAST_element_property_card_1D_h__


// MAST includes
#include "PropertyCards/element_property_card_base.h"
#include "PropertyCards/material_property_card_base.h"


namespace MAST
{
    class ElementPropertyCard1D: public MAST::ElementPropertyCardBase {
        
    public:
        ElementPropertyCard1D():
        MAST::ElementPropertyCardBase(),
        _bending_model(MAST::DEFAULT_BENDING)
        { }
        
        /*!
         *   virtual destructor
         */
        ~ElementPropertyCard1D() { }
        
        /*!
         *   returns the bending model to be used for the 2D element
         */
        void set_bending_model(MAST::BendingOperatorType b)  {
            _bending_model = b;
        }
        
        
        /*!
         *   returns the bending model to be used for the 2D element.
         */
        virtual MAST::BendingOperatorType bending_model(const Elem& elem,
                                                 const FEType& fe) const;
        
        
        /*!
         *    returns the extra quadrature order (on top of the system) that
         *    this element should use. This is elevated by two orders for a DKT
         *    element
         */
        virtual unsigned int extra_quadrature_order(const Elem& elem,
                                                    const FEType& fe) const {
            if (this->bending_model(elem, fe) == MAST::BERNOULLI)
                return 2;
            else
                return 0;
        }
        
    protected:
        
        /*!
         *    initializes the vector to the prestress in the element
         */
        virtual void _prestress_vector(DenseVector<Real>& v) const {
            v.resize(2);
            if (_prestress.size() == 6)
                v(0) = _prestress(0);
        }
        
        
        /*!
         *    initializes the matrix to the prestress in the element
         */
        virtual void _prestress_matrix(DenseMatrix<Real>& m) const {
            m.resize(2, 2);
            if (_prestress.size() == 6)
                m(0,0) = _prestress(0);
        }

        /*!
         *   material property card. By default this chooses DKT for 3 noded
         *   triangles and Mindling for all other elements
         */
        MAST::BendingOperatorType _bending_model;
        
    };
    
    
    class Solid1DSectionElementPropertyCard : public MAST::ElementPropertyCard1D {

    public:
        
        Solid1DSectionElementPropertyCard():
        MAST::ElementPropertyCard1D(),
        _material(NULL)
        { }
        
        /*!
         *   virtual destructor
         */
        ~Solid1DSectionElementPropertyCard() { }
        
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
MAST::ElementPropertyCard1D::bending_model(const Elem& elem,
                                           const FEType& fe) const {
    // for an EDGE2 element, default bending is Bernoulli. For all other elements
    // the default is Timoshenko. Otherwise it returns the model set for
    // this card.
    switch (elem.type()) {
        case EDGE2:
            if ((fe.family == LAGRANGE) &&
                (fe.order  == FIRST) &&
                (_bending_model == MAST::DEFAULT_BENDING))
                return MAST::BERNOULLI;
            else
                return MAST::TIMOSHENKO;
            break;
            
        default:
            if (_bending_model == MAST::DEFAULT_BENDING)
                return MAST::TIMOSHENKO;
            else
                return _bending_model;
            break;
    }
}



inline void
MAST::Solid1DSectionElementPropertyCard::calculate_matrix(const libMesh::Elem &elem,
                                                          MAST::ElemenetPropertyMatrixType t,
                                                          DenseMatrix<Real>& m) const
{
    libmesh_assert(_material); // should have been set
    
    switch (elem.dim()) {
            
        case 1: {
            Real h = this->get<Real>("h")(), // section height
            b = this->get<Real>("b")(),        // section width
            Area = b*h, Iyy = b*pow(h,3)/12., Izz = h*pow(b,3)/12.,
            Iyz = 0., J=1.;
            
            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX: {
                    _material->calculate_1d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m);
                    m.scale_row(0, Area);
                    m.scale_row(1, J);
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
                    // for solid sections with isotropic material this is zero
                    // m(0, 0) and m(0, 1) identify coupling of extension and bending
                    // m(1, 0) and m(1, 1) identify coupling of torsion and bending
                    // m(0,0) = E Lz (y2^2- y1^2)
                    // m(0,1) = E Ly (z2^2- z1^2)
                    m.resize(2,2);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX: {
                    m.resize(2,2);
                    Real E = _material->get<Real>("E")();
                    m(0,0) = E*Iyy;
                    m(0,1) = E*Iyz;
                    m(1,0) = E*Iyz;
                    m(1,1) = E*Izz;
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX:
                    _material->calculate_1d_matrix(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                   m);
                    m.scale(Area);
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
                    _material->calculate_1d_matrix(MAST::MATERIAL_INERTIA_MATRIX,
                                                   m);
                    m.scale(Area);
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
            
        case 2:
        case 3:
        default:
            libmesh_error();
            break;
    }
}



inline void
MAST::Solid1DSectionElementPropertyCard::calculate_matrix_sensitivity(const libMesh::Elem &elem,
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
    
    const MAST::FunctionBase& f = *(it->first);

    DenseMatrix<Real> dm;
    
    switch (elem.dim()) {
            
        case 1: {
            double h = this->get<Real>("h")(), // section height
            b = this->get<Real>("b")(),        // section width
            Area = b*h, Iyy = b*pow(h,3)/12., Izz = h*pow(b,3)/12.,
            dAreadb = h, dAreadh = b,
            dIyydb = pow(h,3)/12., dIyydh = b*pow(h,2)/4.,
            dIzzdb = h*pow(b,2)/4., dIzzdh =  pow(b,3)/12.,
            Iyz = 0., dIyzdb = 0., dIyzdh = 0.,
            J = 0., dJdh = 0., dJdb = 0.;

            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX: {
                    _material->calculate_1d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m);
                    _material->calculate_1d_matrix_sensitivity(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                               dm, p);
                    dm.scale_row(0, Area);
                    dm.scale_row(1, J);
                    
                    if (f.name() == "h") {
                        m.scale_row(0, dAreadh);
                        m.scale_row(1, dJdh);
                    }
                    else if (f.name() == "b") {
                        m.scale_row(0, dAreadb);
                        m.scale_row(1, dJdb);
                    }
                    else
                        m.zero();
                    
                    m.add(1., dm);
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
                    // for solid sections with isotropic material this is zero
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX: {
                    m.resize(2,2); dm.resize(2,2);
                    Real E = _material->get<Real>("E")(),
                    dEdp = _material->get<Real>("E").sensitivity(p);
                    dm(0,0) = dEdp*Iyy;
                    dm(0,1) = dEdp*Iyz;
                    dm(1,0) = dEdp*Iyz;
                    dm(1,1) = dEdp*Izz;

                    if (f.name() == "h") {
                        m(0,0) = E*dIyydh;
                        m(0,1) = E*dIyzdh;
                        m(1,0) = E*dIyzdh;
                        m(1,1) = E*dIzzdh;
                    }
                    else if (f.name() == "b") {
                        m(0,0) = E*dIyydb;
                        m(0,1) = E*dIyzdb;
                        m(1,0) = E*dIyzdb;
                        m(1,1) = E*dIzzdb;
                    }
                    else {
                        m.zero();
                    }
                    
                    m.add(1., dm);
                }
                    break;
                    
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX: {
                    _material->calculate_1d_matrix(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                   m);
                    _material->calculate_1d_matrix_sensitivity(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX,
                                                               dm, p);
                    if (f.name() == "h")
                        m.scale(dAreadh);
                    else if (f.name() == "b")
                        m.scale(dAreadb);
                    else
                        m.zero();
                    
                    m.add(Area, dm);
                }
                    break;
                    
                default:
                    libmesh_error();
                    break;
            }
        }
            break;
            
        case 2:
        case 3:
        default:
            libmesh_error();
            break;
    }
}



inline void
MAST::Solid1DSectionElementPropertyCard::prestress_vector(MAST::ElemenetPropertyMatrixType t,
                                                          const DenseMatrix<Real>& T,
                                                          DenseVector<Real>& v) const {
    double h = this->get<Real>("h")(), // section height
    b = this->get<Real>("b")(),        // section width
    Area = b*h;
    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            _prestress_vector(v);
            v.scale(Area);
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
            // for solid sections with isotropic material this is zero
            v.resize(2);
            break;
            
        default:
            libmesh_error();
            break;
    }
}





inline void
MAST::Solid1DSectionElementPropertyCard::prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                                          const DenseMatrix<Real>& T,
                                                          DenseMatrix<Real>& m) const {
    double h = this->get<Real>("h")(), // section height
    b = this->get<Real>("b")(),        // section width
    Area = b*h;
    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            _prestress_matrix(m);
            m.scale(Area);
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


#endif  // __MAST_element_property_card_1D_h__
