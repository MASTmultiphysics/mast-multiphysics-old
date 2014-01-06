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
        ElementPropertyCard1D(unsigned int pid):
        MAST::ElementPropertyCardBase(pid),
        _bending_model(MAST::DEFAULT_BENDING)
        { }
        
        /*!
         *   virtual destructor
         */
        ~ElementPropertyCard1D() { }
        
        
        /*!
         *   dimension of the element for which this property is defined
         */
        virtual unsigned int dim() const {
            return 1;
        }

        
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
        
        
        /*!
         *   returns value of the property \par val. The string values for 
         *   \par val are IYY, IZZ, IYZ
         */
        virtual Real value(const std::string& val) const = 0;
        
        /*!
         *   vector in the x-y plane of the element. This should not be the same
         *   as the element x-axis.
         */
        Point& y_vector() {
            return _local_y;
        }
        

        /*!
         *   constant reference to vector in the x-y plane of the element. 
         *   This should not be the same as the element x-axis.
         */
        const Point& y_vector() const {
            return _local_y;
        }
        
        
    protected:
        
        /*!
         *    initializes the vector to the prestress in the element
         */
        virtual void _prestress_vector(const DenseMatrix<Real>& T,
                                       DenseVector<Real>& v) const {
            v.resize(2); // zero, if the stress has not been defined

            if (_prestress.m() != 0) {
                
                DenseMatrix<Real> mat;
                _prestress_matrix(T, mat);
                v(0) = mat(0,0); // sigma_xx
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
                mat.right_multiply(T);
                mat.left_multiply_transpose(T);
                
                m(0,0) = mat(0,0);
            }
        }

        /*!
         *   material property card. By default this chooses DKT for 3 noded
         *   triangles and Mindling for all other elements
         */
        MAST::BendingOperatorType _bending_model;
        
        /*!
         *   vector in the x-y plane.
         */
        Point _local_y;
        
    };
    
    
    class Solid1DSectionElementPropertyCard : public MAST::ElementPropertyCard1D {

    public:
        
        Solid1DSectionElementPropertyCard(unsigned int pid):
        MAST::ElementPropertyCard1D(pid),
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
         *   return true if the property is isotropic
         */
        virtual bool if_isotropic() const {
            return true;
        }
        
        
        /*!
         *    returns a reference to the material
         */
        virtual const MAST::MaterialPropertyCardBase& get_material() const {
            libmesh_assert(_material); // make sure it has already been set
            return *_material;
        }
        
        
        /*!
         *   returns value of the property \par val. The string values for
         *   \par val are IYY, IZZ, IYZ
         */
        virtual Real value(const std::string& val) const;
        
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
         *    initializes the vector to the sensitivity of prestress in the element
         */
        virtual void prestress_vector_sensitivity(MAST::ElemenetPropertyMatrixType t,
                                                  const DenseMatrix<Real>& T,
                                                  DenseVector<Real>& v,
                                                  const MAST::SensitivityParameters& p) const;
        
        /*!
         *    initializes the matrix to the prestress in the element
         */
        virtual void prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                      const DenseMatrix<Real>& T,
                                      DenseMatrix<Real>& m) const;
        
        /*!
         *    initializes the matrix to the sensitivity of prestress in the element
         */
        virtual void prestress_matrix_sensitivity(MAST::ElemenetPropertyMatrixType t,
                                                  const DenseMatrix<Real>& T,
                                                  DenseMatrix<Real>& m,
                                                  const MAST::SensitivityParameters& p) const;
        
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



inline Real
MAST::Solid1DSectionElementPropertyCard::value(const std::string& val) const {
    
    Real h_y = this->get<Real>("h_y")(), // section height
    h_z = this->get<Real>("h_z")(),        // section width
    off_h_y = 0., off_h_z = 0.;
    
    if (_properties.count("off_h_y"))
        off_h_y = this->get<Real>("off_h_y")();
    if (_properties.count("off_h_z"))
        off_h_z = this->get<Real>("off_h_z")();

    Real Area = h_y*h_z,
    Iyy = h_y*(pow(h_z,3)/12. + pow(off_h_z,2)*h_z),
    Izz = h_z*(pow(h_y,3)/12. + pow(off_h_y,2)*h_y),
    Iyz = (off_h_y*off_h_z*h_y*h_z),
    J=1.;

    if (val == "A")
        return Area;
    else if (val == "IYY")
        return Iyy;
    else if (val == "IZZ")
        return Izz;
    else if (val == "IYZ")
        return Iyz;
    else if (val == "J")
        return J;
    
    // should not get here
    libmesh_error();
    return 0.;
}




inline void
MAST::Solid1DSectionElementPropertyCard::calculate_matrix(const libMesh::Elem &elem,
                                                          MAST::ElemenetPropertyMatrixType t,
                                                          DenseMatrix<Real>& m) const
{
    libmesh_assert(_material); // should have been set
    
    switch (elem.dim()) {
            
        case 1: {
            Real h_y = this->get<Real>("h_y")(), // section height
            h_z = this->get<Real>("h_z")(),      // section width
            off_h_y = 0., off_h_z = 0.;          // offset values of mid-plane from longitudinal axis
            
            if (_properties.count("off_h_y"))
                off_h_y = this->get<Real>("off_h_y")();
            if (_properties.count("off_h_z"))
                off_h_z = this->get<Real>("off_h_z")();
            
            Real Area = h_y*h_z,
            Iyy = h_y*(pow(h_z,3)/12. + pow(off_h_z,2)*h_z),   // used for w-bending
            Izz = h_z*(pow(h_y,3)/12. + pow(off_h_y,2)*h_y),   // used for v-bending
            Iyz = (off_h_y*off_h_z*h_y*h_z),
            J=1.;
            
            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX: {
                    _material->calculate_1d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m);
                    m.scale_row(0, Area);
                    m.scale_row(1, J);
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX: {
                    // for solid sections with isotropic material this is zero
                    // m(0, 0) and m(0, 1) identify coupling of extension and bending
                    // m(1, 0) and m(1, 1) identify coupling of torsion and bending
                    Real E = _material->get<Real>("E")();
                    m.resize(2,2);
                    m(0,0) = E*h_y*h_z*off_h_y;  // u <-> v
                    m(0,1) = E*h_y*h_z*off_h_z;  // u <-> w
                    m(1,0) = 0.;           // tx <-> v
                    m(1,1) = 0.;           // tx <-> w
                }
                    break;
                    
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX: {
                    m.resize(2,2);
                    Real E = _material->get<Real>("E")();
                    m(0,0) = E*Izz;
                    m(0,1) = E*Iyz;
                    m(1,0) = E*Iyz;
                    m(1,1) = E*Iyy;
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
    
    const bool depends = this->depends_on(f);
    
    switch (elem.dim()) {
            
        case 1: {
            Real h_y = this->get<Real>("h_y")(), // section height
            h_z = this->get<Real>("h_z")(),        // section width
            off_h_y = 0., off_h_z = 0.;          // offset values of mid-plane from longitudinal axis
            
            if (_properties.count("off_h_y"))
                off_h_y = this->get<Real>("off_h_y")();
            if (_properties.count("off_h_z"))
                off_h_z = this->get<Real>("off_h_z")();
            
            Real Area = h_y*h_z,
            Iyy = h_y*(pow(h_z,3)/12. + pow(off_h_z,2)*h_z),   // used for w-bending
            Izz = h_z*(pow(h_y,3)/12. + pow(off_h_y,2)*h_y),   // used for v-bending
            Iyz = (off_h_y*off_h_z*h_y*h_z),
            
            dAreadhy = h_z, dAreadhz = h_y,
            
            dIyydhy = (pow(h_z,3)/12. + pow(off_h_z,2)*h_z),
            dIyydhz = h_y*(pow(h_z,2)/4. + pow(off_h_z,2)),
            dIzzdhy = h_z*(pow(h_y,2)/4. + pow(off_h_y,2)),
            dIzzdhz = (pow(h_y,3)/12. + pow(off_h_y,2)*h_y),
            
            dIyzdhy = off_h_y*off_h_z*h_z,
            dIyzdhz = off_h_y*off_h_z*h_y,
            
            J = 1., dJdhy = 0., dJdhz = 0.;

            switch (t) {
                case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX: {
                    _material->calculate_1d_matrix(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                   m);
                    _material->calculate_1d_matrix_sensitivity(MAST::MATERIAL_STIFFNESS_MATRIX,
                                                               dm, p);
                    dm.scale_row(0, Area);
                    dm.scale_row(1, J);
                    
                    if (depends && f.name() == "h_y") {
                        m.scale_row(0, dAreadhy);
                        m.scale_row(1, dJdhy);
                    }
                    else if (depends && f.name() == "h_z") {
                        m.scale_row(0, dAreadhz);
                        m.scale_row(1, dJdhz);
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

                    if (depends && f.name() == "h_y") {
                        m(0,0) = E*dIyydhy;
                        m(0,1) = E*dIyzdhy;
                        m(1,0) = E*dIyzdhy;
                        m(1,1) = E*dIzzdhy;
                    }
                    else if (depends && f.name() == "h_z") {
                        m(0,0) = E*dIyydhz;
                        m(0,1) = E*dIyzdhz;
                        m(1,0) = E*dIyzdhz;
                        m(1,1) = E*dIzzdhz;
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
                    if (depends && f.name() == "h_y")
                        m.scale(dAreadhy);
                    else if (depends && f.name() == "h_z")
                        m.scale(dAreadhz);
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
    // prestress is independent of cross-section.

    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            _prestress_vector(T, v);
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
MAST::Solid1DSectionElementPropertyCard::prestress_vector_sensitivity(MAST::ElemenetPropertyMatrixType t,
                                                                      const DenseMatrix<Real>& T,
                                                                      DenseVector<Real>& v,
                                                                      const MAST::SensitivityParameters& p) const {
    // only first order sensitivities are calculated at this point
    libmesh_assert_equal_to(p.total_order(), 1);
    // prestress is independent of cross-section.
    v.resize(2);
}





inline void
MAST::Solid1DSectionElementPropertyCard::prestress_matrix(MAST::ElemenetPropertyMatrixType t,
                                                          const DenseMatrix<Real>& T,
                                                          DenseMatrix<Real>& m) const {
    // prestress is independent of cross-section.
    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            _prestress_matrix(T, m);
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



inline void
MAST::Solid1DSectionElementPropertyCard::prestress_matrix_sensitivity(MAST::ElemenetPropertyMatrixType t,
                                                                      const DenseMatrix<Real>& T,
                                                                      DenseMatrix<Real>& m,
                                                                      const MAST::SensitivityParameters& p) const {
    // only first order sensitivities are calculated at this point
    libmesh_assert_equal_to(p.total_order(), 1);
    
    // prestress is independent of cross-section.
    m.resize(2,2);
}


#endif  // __MAST_element_property_card_1D_h__
