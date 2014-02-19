//
//  multilayer_1d_section_element_property_card.h
//  MAST
//
//  Created by Manav Bhatia on 2/18/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

#ifndef __MAST_multilayer_1d_section_element_property_card_h__
#define __MAST_multilayer_1d_section_element_property_card_h__


// MAST includes
#include "PropertyCards/element_property_card_1D.h"
#include "PropertyCards/solid_1d_section_element_property_card.h"


namespace MAST {
    class Multilayer1DSectionElementPropertyCard : public MAST::ElementPropertyCard1D {
        
    public:
        
        Multilayer1DSectionElementPropertyCard(unsigned int pid):
        MAST::ElementPropertyCard1D(pid)
        { }
        
        
        /*!
         *   virtual destructor
         */
        virtual ~Multilayer1DSectionElementPropertyCard() {
            // delete the layer offset functions
            for (unsigned int i=0; i<_layer_offsets.size(); i++)
                delete _layer_offsets[i];
        }
        
        
        class LayerOffset: public MAST::FieldFunction<Real> {
        public:
            LayerOffset(unsigned int layer_num,
                        std::vector<MAST::FieldFunction<Real>*>& layer_hz):
            MAST::FieldFunction<Real>("hz_offset"),
            _layer_num(layer_num),
            _layer_hz(layer_hz) {
                for (unsigned int i=0; i < _layer_hz.size(); i++)
                    _functions.insert(_layer_hz[i]);
            }
            
            LayerOffset(const MAST::Multilayer1DSectionElementPropertyCard::LayerOffset &f):
            MAST::FieldFunction<Real>(f),
            _layer_num(f._layer_num)
            {
                // initialize the vector
                _layer_hz.resize(f._layer_hz.size());
                for (unsigned int i=0; i < _layer_hz.size(); i++) {
                    _layer_hz[i] = f._layer_hz[i]->clone().release();
                    _functions.insert(_layer_hz[i]);
                }
            }
                
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<Real> > clone() const {
                return std::auto_ptr<MAST::FieldFunction<Real> >
                (new MAST::Multilayer1DSectionElementPropertyCard::LayerOffset(*this));
            }
            
            virtual ~LayerOffset() {
                // delete all the layer functions
                for (unsigned int i=0; i<_layer_hz.size(); i++)
                    delete _layer_hz[i];
            }

            virtual void operator() (const Point& p, const Real t, Real& m) const {
                Real val = 0., hi=0.;
                for (unsigned int i=0; i<_layer_num; i++) {
                    (*_layer_hz[i])(p, t, val);
                    hi += val; // currently the offset is chosen from h=0;
                }
                // finally, add half of the current layer thickness
                (*_layer_hz[_layer_num])(p, t, val);
                hi += 0.5*val;
            }
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, Real& m) const {
                libmesh_error(); // to be implemented
            }
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, Real& m) const {
                libmesh_error(); // to be implemented
            }
            
        protected:
            
            const unsigned int _layer_num;
            std::vector<MAST::FieldFunction<Real>*> _layer_hz;
        };
        
        
        class SectionIntegratedMatrix: public MAST::FieldFunction<DenseMatrix<Real> > {
        public:
            SectionIntegratedMatrix(std::vector<MAST::FieldFunction<DenseMatrix<Real> >*>& layer_mats):
            MAST::FieldFunction<DenseMatrix<Real> >("SectionIntegratedMatrix1D"),
            _layer_mats(layer_mats) {
                for (unsigned int i=0; i < _layer_mats.size(); i++) {
                    _functions.insert(_layer_mats[i]);
                }
            }
            
            SectionIntegratedMatrix(const MAST::Multilayer1DSectionElementPropertyCard::SectionIntegratedMatrix &f):
            MAST::FieldFunction<DenseMatrix<Real> >(f) {
                // initialize the vector
                _layer_mats.resize(f._layer_mats.size());
                for (unsigned int i=0; i < _layer_mats.size(); i++) {
                    _layer_mats[i] = f._layer_mats[i]->clone().release();
                    _functions.insert(_layer_mats[i]);
                }
            }
            
            /*!
             *   @returns a clone of the function
             */
            virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>> clone() const {
                return std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
                (new MAST::Multilayer1DSectionElementPropertyCard::SectionIntegratedMatrix(*this));
            }
            
            virtual ~SectionIntegratedMatrix() {
                // delete all the layer functions
                for (unsigned int i=0; i<_layer_mats.size(); i++)
                    delete _layer_mats[i];
            }
            
            virtual void operator() (const Point& p, const Real t, DenseMatrix<Real>& m) const {
                // add the values of each matrix to get the integrated value
                DenseMatrix<Real> mi;
                m.resize(3,3);
                for (unsigned int i=0; i<_layer_mats.size(); i++) {
                    (*_layer_mats[i])(p, t, mi);
                    m.add(1., mi);
                }
            }
            
            virtual void partial (const MAST::FieldFunctionBase& f,
                                  const Point& p, const Real t, DenseMatrix<Real>& m) const {
                // add the values of each matrix to get the integrated value
                DenseMatrix<Real> mi;
                m.resize(3,3);
                for (unsigned int i=0; i<_layer_mats.size(); i++) {
                    _layer_mats[i]->partial(f, p, t, mi);
                    m.add(1., mi);
                }
            }
            
            virtual void total (const MAST::FieldFunctionBase& f,
                                const Point& p, const Real t, DenseMatrix<Real>& m) const {
                // add the values of each matrix to get the integrated value
                DenseMatrix<Real> mi;
                m.resize(3,3);
                for (unsigned int i=0; i<_layer_mats.size(); i++) {
                    _layer_mats[i]->total(f, p, t, mi);
                    m.add(1., mi);
                }
            }

            
        protected:
            
            std::vector<MAST::FieldFunction<DenseMatrix<Real>>*> _layer_mats;
        };
        
        
        /*!
         *    sets the layers of this section
         */
        void set_layers(std::vector<MAST::Solid1DSectionElementPropertyCard*>& layers) {
            
            // make sure that this has not been already set
            libmesh_assert(_layers.size() == 0);
            
            // now create the vector of offsets for each later
            std::vector<MAST::FieldFunction<Real>*> layer_hz(_layers.size());
            _layer_offsets.resize(_layers.size());
            
            for (unsigned int i=0; i<_layers.size(); i++)
                layer_hz[i] = _layers[i]->get<MAST::FieldFunction<Real> >("hz").clone().release();
            
            for (unsigned int i=0; i<_layers.size(); i++) {
                // create the offset function
                _layer_offsets[i] = new MAST::Multilayer1DSectionElementPropertyCard::LayerOffset
                (i, layer_hz);
                // tell the layer about the offset
                _layers[i]->add(*_layer_offsets[i]);
            }
        }
        
        /*!
         *   return true if the property is isotropic
         */
        virtual bool if_isotropic() const {
            return false;
        }
        
        
        /*!
         *   returns value of the property \par val. The string values for
         *   \par val are A, J, IYY, IZZ, IYZ
         */
        virtual Real value(const std::string& val) const {
            libmesh_assert(false); // not implemented for multilayered sections yet
        }
        
        /*!
         *   returns a function to evaluate the specified quantitys
         *   type \par t.
         */
        virtual std::auto_ptr<MAST::FieldFunction<DenseMatrix<Real>>>
        get_property(MAST::ElemenetPropertyMatrixType t,
                     const MAST::StructuralElementBase& e) const;
        
        /*!
         *  returns true if the property card depends on the function \p f
         */
        virtual bool depends_on(const MAST::FieldFunctionBase& f) const {
            // ask each layer for the dependence
            for (unsigned int i=0; i<_layers.size(); i++)
                if (_layers[i]->depends_on(f))
                    return true;

            // ask each offset for the dependence
            for (unsigned int i=0; i<_layer_offsets.size(); i++)
                if (_layer_offsets[i]->depends_on(f))
                    return true;

            // if it gets here, then there is no dependence
            return false;
        }
        
    protected:

        std::vector<MAST::FieldFunction<Real>*> _layer_offsets;
        
        /*!
         *   vector of thickness function for each layer
         */
        std::vector<MAST::Solid1DSectionElementPropertyCard*> _layers;
    };
    
}




#endif // __MAST_multilayer_1d_section_element_property_card_h__
