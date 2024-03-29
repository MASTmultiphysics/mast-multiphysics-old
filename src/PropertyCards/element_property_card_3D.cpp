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

// MAST includes
#include "PropertyCards/element_property_card_3D.h"


MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix::
SectionIntegratedStiffnessMatrix(MAST::FieldFunction<DenseRealMatrix > *mat):
MAST::FieldFunction<DenseRealMatrix > ("SectionIntegratedStiffnessMatrix3D"),
_material_stiffness(mat) {
    _functions.insert(mat);
}


void
MAST::ElementPropertyCard3D::
SectionIntegratedStiffnessMatrix::operator() (const libMesh::Point& p,
                                              const Real t,
                                              DenseRealMatrix& m) const {
    // this only returns the material stiffness
    (*_material_stiffness)
    (p, t, m);
}


void
MAST::ElementPropertyCard3D::
SectionIntegratedStiffnessMatrix::partial (const MAST::FieldFunctionBase& f,
                                           const libMesh::Point& p,
                                           const Real t,
                                           DenseRealMatrix& m) const {
    // this only returns the material stiffness
    _material_stiffness->partial(f, p, t, m);
}


void
MAST::ElementPropertyCard3D::
SectionIntegratedStiffnessMatrix::total (const MAST::FieldFunctionBase& f,
                                         const libMesh::Point& p,
                                         const Real t,
                                         DenseRealMatrix& m) const {
    // this only returns the material stiffness
    _material_stiffness->total(f, p, t, m);
}




MAST::ElementPropertyCard3D::
SectionIntegratedInertiaMatrix::SectionIntegratedInertiaMatrix(MAST::FieldFunction<DenseRealMatrix > *mat):
MAST::FieldFunction<DenseRealMatrix >("SectionIntegratedInertiaMatrix3D"),
_material_inertia(mat) {
    _functions.insert(mat);
}



void
MAST::ElementPropertyCard3D::
SectionIntegratedInertiaMatrix::operator() (const libMesh::Point& p,
                                            const Real t,
                                            DenseRealMatrix& m) const {
    DenseRealMatrix mat;
    m.resize(6, 6);
    // this only returns the material inertia
    (*_material_inertia)(p, t, mat);
    for (unsigned int i=0; i<3; i++) {
        for (unsigned int j=0; j<3; j++) {
            m(i,j) = mat(i,j);
        }
        m(i+3,i+3) = mat(i,i) * 1.0e-12;
    }
}



void
MAST::ElementPropertyCard3D::
SectionIntegratedInertiaMatrix::partial (const MAST::FieldFunctionBase& f,
                                         const libMesh::Point& p,
                                         const Real t,
                                         DenseRealMatrix& m) const {
    DenseRealMatrix mat;
    m.resize(6, 6);
    // this only returns the material inertia
    // sensitivity of rotary inertia is assumed to be zero
    _material_inertia->partial(f, p, t, mat);
    for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<3; j++)
            m(i,j) = mat(i,j);
}


void
MAST::ElementPropertyCard3D::
SectionIntegratedInertiaMatrix::total (const MAST::FieldFunctionBase& f,
                                       const libMesh::Point& p,
                                       const Real t,
                                       DenseRealMatrix& m) const {
    DenseRealMatrix mat;
    m.resize(6, 6);
    // this only returns the material inertia
    // sensitivity of rotary inertia is assumed to be zero
    _material_inertia->total(f, p, t, mat);
    for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<3; j++)
            m(i,j) = mat(i,j);
}




MAST::ElementPropertyCard3D::SectionIntegratedThermalExpansionMatrix::
SectionIntegratedThermalExpansionMatrix(MAST::FieldFunction<DenseRealMatrix > *mat_stiff,
                                        MAST::FieldFunction<DenseRealMatrix > *mat_expansion):
MAST::FieldFunction<DenseRealMatrix >("SectionIntegratedThermalExpansionMatrix3D"),
_material_stiffness(mat_stiff),
_material_expansion(mat_expansion) {
    _functions.insert(mat_stiff);
    _functions.insert(mat_expansion);
}




void
MAST::ElementPropertyCard3D::SectionIntegratedThermalExpansionMatrix::operator() (const libMesh::Point& p,
                                                                                  const Real t,
                                                                                  DenseRealMatrix& m) const {
    DenseRealMatrix mat;
    (*_material_stiffness)(p, t, m);
    (*_material_expansion)(p, t, mat);
    m.right_multiply(mat);
}




void
MAST::ElementPropertyCard3D::SectionIntegratedThermalExpansionMatrix::partial (const MAST::FieldFunctionBase& f,
                                                                               const libMesh::Point& p,
                                                                               const Real t,
                                                                               DenseRealMatrix& m) const {
    libmesh_error();
}




void
MAST::ElementPropertyCard3D::SectionIntegratedThermalExpansionMatrix::total (const MAST::FieldFunctionBase& f,
                                                                             const libMesh::Point& p,
                                                                             const Real t,
                                                                             DenseRealMatrix& m) const {
    libmesh_error();
}




MAST::ElementPropertyCard3D::SectionIntegratedPrestressAMatrix::
SectionIntegratedPrestressAMatrix(MAST::FieldFunction<DenseRealMatrix > *prestress):
MAST::SectionIntegratedPrestressMatrixBase("SectionIntegratedPrestressAMatrix3D"),
_prestress(prestress){
    _functions.insert(prestress);
}




void
MAST::ElementPropertyCard3D::
SectionIntegratedPrestressAMatrix::operator() (const libMesh::Point& p,
                                               const Real t,
                                               DenseRealMatrix& m) const {
    (*_prestress)(p, t, m);
}





void
MAST::ElementPropertyCard3D::
SectionIntegratedPrestressAMatrix::partial (const MAST::FieldFunctionBase& f,
                                            const libMesh::Point& p,
                                            const Real t,
                                            DenseRealMatrix& m) const {
    _prestress->partial(f, p, t, m);
}



void
MAST::ElementPropertyCard3D::
SectionIntegratedPrestressAMatrix::total (const MAST::FieldFunctionBase& f,
                                          const libMesh::Point& p,
                                          const Real t,
                                          DenseRealMatrix& m) const {
    _prestress->total(f, p, t, m);
}



void
MAST::ElementPropertyCard3D::
SectionIntegratedPrestressAMatrix::convert_to_vector(const DenseRealMatrix &m,
                                                     DenseRealVector &v) const {
    libmesh_error();
}


std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > >
MAST::ElementPropertyCard3D::get_property(MAST::ElemenetPropertyMatrixType t,
                                          const MAST::StructuralElementBase& e) const {
    
    std::auto_ptr<MAST::FieldFunction<DenseRealMatrix > > rval;

    switch (t) {
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_A_MATRIX:
            rval.reset(new MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 3).release()));
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_B_MATRIX:
            rval.reset(new MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 3).release()));
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_STIFFNESS_D_MATRIX:
            rval.reset(new MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 3).release()));
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX:
            rval.reset(new MAST::ElementPropertyCard3D::SectionIntegratedStiffnessMatrix
                       (_material->get_property(MAST::MATERIAL_TRANSVERSE_SHEAR_STIFFNESS_MATRIX, *this, 3).release()));
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_INERTIA_MATRIX:
        case MAST::SECTION_INTEGRATED_PRESTRESS_A_MATRIX:
        case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_A_MATRIX:
        case MAST::SECTION_INTEGRATED_MATERIAL_THERMAL_EXPANSION_B_MATRIX:
            rval.reset(new MAST::ElementPropertyCard3D::SectionIntegratedThermalExpansionMatrix
                       (_material->get_property(MAST::MATERIAL_STIFFNESS_MATRIX, *this, 3).release(),
                        _material->get_property(MAST::MATERIAL_THERMAL_EXPANSION_MATRIX, *this, 3).release()));
            break;
            
        case MAST::SECTION_INTEGRATED_MATERIAL_DAMPING_MATRIX:
        default:
            libmesh_error();
            break;
    }
    
    return rval;
}



