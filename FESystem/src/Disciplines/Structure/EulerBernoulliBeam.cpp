//
//  EulerBernoulliBeam.cpp
//  FESystem
//
//  Created by Manav Bhatia on 8/13/12.
//
//

// FESystem includes
#include "Disciplines/Structure/EulerBernoulliBeam.h"
#include "Mesh/ElemBase.h"
#include "Numerics/MatrixBase.h"


FESystem::Structures::EulerBernoulliBeam::EulerBernoulliBeam():
FESystem::Structures::LinearBeamElementBase()
{
}


FESystem::Structures::EulerBernoulliBeam::~EulerBernoulliBeam()
{
    
}



void
FESystem::Structures::EulerBernoulliBeam::calculateConsistentMassMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    
}



void
FESystem::Structures::EulerBernoulliBeam::calculateStiffnessMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    
    FESystemDouble length=this->geometric_elem->getElementSize(*(this->finite_element), *(this->quadrature)), EI_trans = 0.0, EI_chord = 0.0,
    EIl_term_trans = EI_trans / pow(length,3), EIl_term_chord = EI_chord / pow(length,3);
    
    //     v1 : chordwise bending
    mat.setVal(2,2, EIl_term_chord * 12.0);
    mat.setVal(2,6, EIl_term_chord * 6.0 * length);
    mat.setVal(2,8, EIl_term_chord * (-12.0));
    mat.setVal(2,12, EIl_term_chord * 6.0 * length);
    
    //     w1 : spanwise bending
    mat.setVal(3,3, EIl_term_trans * 12.0);
    mat.setVal(3,5, EIl_term_trans * (-6.0) * length);
    mat.setVal(3,9, EIl_term_trans * (-12.0));
    mat.setVal(3,11, EIl_term_trans * (-6.0) * length);
    
    
    //     thetay1 : spanwise bending
    mat.setVal(5,3, EIl_term_trans * (-6.0) * length);
    mat.setVal(5,5, EIl_term_trans * 4.0 * pow(length,2));
    mat.setVal(5,9, EIl_term_trans * 6.0 * length);
    mat.setVal(5,11, EIl_term_trans * 2.0 * pow(length,2));
    
    //     thetaz1 : chordwise bending
    mat.setVal(6,2, EIl_term_chord * 6.0 * length);
    mat.setVal(6,6, EIl_term_chord * 4.0 * pow(length,2));
    mat.setVal(6,8, EIl_term_chord * (-6.0) * length);
    mat.setVal(6,12, EIl_term_chord * 2.0 * pow(length,2));
    
    
    //     v2 :  chordwise bending
    mat.setVal(8,2, EIl_term_chord * (-12.0));
    mat.setVal(8,6, EIl_term_chord * (-6.0) * length);
    mat.setVal(8,8, EIl_term_chord * 12.0);
    mat.setVal(8,12, EIl_term_chord * (-6.0) * length);
    
    //     w2 :  spanwise bending
    mat.setVal(9,3, EIl_term_trans * (-12.0));
    mat.setVal(9,5, EIl_term_trans * 6.0 * length);
    mat.setVal(9,9, EIl_term_trans * 12.0);
    mat.setVal(9,11, EIl_term_trans * 6.0 * length);
    
    
    //     thetay2 : spanwise bending
    mat.setVal(11,3, EIl_term_trans * (-6.0) * length);
    mat.setVal(11,5, EIl_term_trans * 2.0 * pow(length,2));
    mat.setVal(11,9, EIl_term_trans * 6.0 * length);
    mat.setVal(11,11, EIl_term_trans * 4.0 * pow(length,2));
    
    //     thetaz2 : chordwise bending
    mat.setVal(12,2, EIl_term_chord * 6.0 * length);
    mat.setVal(12,6, EIl_term_chord * 2.0 * pow(length,2));
    mat.setVal(12,8, EIl_term_chord * (-6.0) * length);
    mat.setVal(12,12, EIl_term_chord * 4.0 * pow(length,2));
}




