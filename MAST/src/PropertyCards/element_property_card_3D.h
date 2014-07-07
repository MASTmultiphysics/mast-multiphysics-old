/*
 * MAST: Multidisciplinary-design Adaptation and Sensitivity Toolkit
 * Copyright (C) 2013-2014  Manav Bhatia
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef __MAST_element_property_card_3D_h__
#define __MAST_element_property_card_3D_h__

// MAST includes
#include "PropertyCards/element_property_card_base.h"
#include "PropertyCards/material_property_card_base.h"


namespace MAST
{
    
    class ElementPropertyCard3D: public MAST::ElementPropertyCardBase {
        
    public:
    public:
        ElementPropertyCard3D(unsigned int pid):
        MAST::ElementPropertyCardBase(pid),
        _material(NULL)
        { }
        
        /*!
         *   virtual destructor
         */
        virtual ~ElementPropertyCard3D() { }
        
        
        
        class SectionIntegratedStiffnessMatrix: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            SectionIntegratedStiffnessMatrix(MAST::FieldFunction<DenseRealMatrix > *mat);
            
            SectionIntegratedStiffnessMatrix(const MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _material_stiffness(f._material_stiffness->clone().release()) {
                _functions.insert(_material_stiffness);
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix(*this));
            }

            virtual ~SectionIntegratedStiffnessMatrix() { delete _material_stiffness;}
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                             const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                           const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseRealMatrix > *_material_stiffness;
        };
        
        
        
        class SectionIntegratedInertiaMatrix: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            SectionIntegratedInertiaMatrix(MAST::FieldFunction<DenseRealMatrix > *mat);
            
            SectionIntegratedInertiaMatrix(const MAST::ElementPropertyCard3D::SectionIntegratedInertiaMatrix& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _material_inertia(f._material_inertia->clone().release()) {
                _functions.insert(_material_inertia);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::ElementPropertyCard3D::SectionIntegratedInertiaMatrix(*this));
            }

            virtual ~SectionIntegratedInertiaMatrix() { delete _material_inertia;}
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                             const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                           const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseRealMatrix > *_material_inertia;
        };
        
        
        
        class SectionIntegratedThermalExpansionMatrix: public MAST::FieldFunction<DenseRealMatrix > {
        public:
            SectionIntegratedThermalExpansionMatrix(MAST::FieldFunction<DenseRealMatrix > *mat_stiff,
                                                    MAST::FieldFunction<DenseRealMatrix > *mat_expansion);
            
            SectionIntegratedThermalExpansionMatrix(const MAST::ElementPropertyCard3D::
                                                    SectionIntegratedThermalExpansionMatrix& f):
            MAST::FieldFunction<DenseRealMatrix >(f),
            _material_stiffness(f._material_stiffness->clone().release()),
            _material_expansion(f._material_expansion->clone().release()) {
                _functions.insert(_material_stiffness);
                _functions.insert(_material_expansion);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::ElementPropertyCard3D::SectionIntegratedThermalExpansionMatrix(*this));
            }

            virtual ~SectionIntegratedThermalExpansionMatrix() {
                delete _material_stiffness;
                delete _material_expansion;
            }
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                             const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                           const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
        protected:
            
            MAST::FieldFunction<DenseRealMatrix > *_material_stiffness;
            MAST::FieldFunction<DenseRealMatrix > *_material_expansion;
        };
        
        
        
        
        class SectionIntegratedPrestressAMatrix: public MAST::SectionIntegratedPrestressMatrixBase {
        public:
            SectionIntegratedPrestressAMatrix(MAST::FieldFunction<DenseRealMatrix > *prestress);
            
            SectionIntegratedPrestressAMatrix(const MAST::ElementPropertyCard3D::SectionIntegratedPrestressAMatrix& f):
            MAST::SectionIntegratedPrestressMatrixBase(f),
            _prestress(f._prestress->clone().release()) {
                _functions.insert(_prestress);
            }

            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
                (new MAST::ElementPropertyCard3D::SectionIntegratedPrestressAMatrix(*this));
            }

            virtual ~SectionIntegratedPrestressAMatrix() { delete _prestress;}
            
            virtual void operator() (const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                             const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                           const libMesh::Point& p, const Real t, DenseRealMatrix& m) const;
            
            virtual void convert_to_vector(const DenseRealMatrix& m, DenseRealVector& v) const;
            
        protected:
            
            MAST::FieldFunction<DenseRealMatrix > *_prestress;
        };

        
        /*!
         *   dimension of the element for which this property is defined
         */
        virtual unsigned int dim() const {
            return 3;
        }

        /*!
         *   return true if the property is isotropic
         */
        virtual bool if_isotropic() const {
            return true;
        }

        
        /*!
         *    sets the material card
         */
        virtual void set_material(MAST::MaterialPropertyCardBase& mat) {
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
         *   returns a function to evaluate the specified quantitys
         *   type \par t.
         */
        virtual std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
        get_property(MAST::ElemenetPropertyMatrixType t,
                     const MAST::StructuralElementBase& e) const;

        /*!
         *  returns true if the property card depends on the function \p f
         */
        virtual bool depends_on(const MAST::FieldFunctionBase& f) const {
            return _material->depends_on(f) ||            // check if the material property depends on the function
            MAST::ElementPropertyCardBase::depends_on(f); // check with this property card
        }
        
        
        protected:
        
        /*!
         *    pointer to the material property card
         */
        MAST::MaterialPropertyCardBase* _material;
    };
}


#endif // __MAST_element_property_card_3D_h__
