//
//  flutter_solution.cpp
//  RealSolver
//
//  Created by Manav Bhatia on 4/15/14.
//  Copyright (c) 2014 Manav Bhatia. All rights reserved.
//

// C++ includes
#include <iomanip>

// MAST includes
#include "Aeroelasticity/flutter_solution.h"
#include "Aeroelasticity/flutter_solver_base.h"



void
MAST::TimeDomainFlutterRoot::init(const libMesh::Real ref_val, const libMesh::Real b_ref,
                                  const libMesh::Complex num,
                                  const libMesh::Complex den,
                                  const RealMatrixX& Bmat,
                                  const ComplexVectorX& eig_vec)
{
    V = ref_val;
    
    if (std::abs(den) > 0.)
    {
        root = num/den;
        if (std::real(root) > 0.)
        {
            V     = sqrt(1./std::real(root));
            g     = std::imag(root)/std::real(root);
            omega = k_ref*V/b_ref;
            if_nonphysical_root = false;
        }
        else
        {
            V     = 0.;
            g     = 0.;
            omega = 0.;
            if_nonphysical_root = true;
        }
    }
    
    // calculate the modal participation vector
    const unsigned int nvals = (int)Bmat.rows();
    mode = eig_vec;
    ComplexVectorX k_q = Bmat * mode;
    modal_participation.resize(nvals, 1);
    for (unsigned int i=0; i<nvals; i++)
        modal_participation(i) =  std::abs(std::conj(mode(i)) * k_q(i));
    modal_participation *= (1./modal_participation.sum());
}


void
MAST::FlutterSolutionBase::sort(const MAST::FlutterSolutionBase& sol)
{
    const unsigned int nvals = this->n_roots();
    libmesh_assert_equal_to(nvals, sol.n_roots());
    
    // two roots with highest modal_participation dot product are sorted
    // in the same serial order
    for (unsigned int i=0; i<nvals-1; i++)
    {
        const MAST::FlutterRootBase& r = sol.get_root(i);
        libMesh::Real max_val = 0., val = 0.;
        unsigned int max_val_root = nvals+1;
        for (unsigned int j=i; j<nvals; j++) {
            val = _roots[j]->modal_participation.dot(r.modal_participation);
            if (val > max_val) {
                max_val = val;
                max_val_root = j;
            }
        }
        
        // now we should have the one with highest dot product
        std::swap(_roots[i], _roots[max_val_root]);
    }
}





void
MAST::FrequencyDomainFlutterSolution::init (const libMesh::Real ref_val, const libMesh::Real bref,
                                            const LAPACK_ZGGEV& eig_sol)
{
    // make sure that it hasn't already been initialized
    libmesh_assert(!_roots.size());
    
    _ref_val = ref_val;
    // iterate over the roots and initialize the vector
    const ComplexMatrixX &B = eig_sol.B(), &VR = eig_sol.right_eigenvectors();
    const ComplexVectorX &num = eig_sol.alphas(), &den = eig_sol.betas();
    unsigned int nvals = (int)B.rows();

    _roots.resize(nvals);
    for (unsigned int i=0; i<nvals; i++) {
        _roots[i] =  _solver.build_flutter_root();
        dynamic_cast<MAST::FrequencyDomainFlutterRoot*>
        (_roots[i])->init(ref_val, bref, num(i), den(i), B, VR.col(i));
    }
}




void
MAST::TimeDomainFlutterSolution::init (const libMesh::Real ref_val, const libMesh::Real bref,
                                       const LAPACK_DGGEV& eig_sol)
{
    // make sure that it hasn't already been initialized
    libmesh_assert(!_roots.size());
    
    _ref_val = ref_val;
    // iterate over the roots and initialize the vector
    const RealMatrixX &B = eig_sol.B();
    const ComplexMatrixX &VR = eig_sol.right_eigenvectors();
    const ComplexVectorX &num = eig_sol.alphas(), &den = eig_sol.betas();
    unsigned int nvals = (int)B.rows();
    
    _roots.resize(nvals);
    for (unsigned int i=0; i<nvals; i++) {
        _roots[i] =  _solver.build_flutter_root();
        dynamic_cast<MAST::TimeDomainFlutterRoot*>
        (_roots[i])->init(ref_val, bref, num(i), den(i), B, VR.col(i));
    }
}




void
MAST::FlutterSolutionBase::print(std::ostream &output,
                                 std::ostream &mode_output)
{
    const unsigned int nvals = this->n_roots();
    libmesh_assert(nvals);
    
    output << "k = " << _ref_val << std::endl
    << std::setw(5) << "#"
    << std::setw(15) << "Re"
    << std::setw(15) << "Im"
    << std::setw(15) << "g"
    << std::setw(15) << "V"
    << std::setw(15) << "omega";
    
    // output the headers for flutter mode participation
    for (unsigned int i=0; i<nvals; i++)
        output
        << std::setw(2) << " "
        << std::setw(5) << "Mode "
        << std::setw(5) << i
        << std::setw(3) << " ";
    output << std::endl;
    
    // output the headers for flutter mode
    mode_output << "k = " << _ref_val << "  (Right Eigenvectors)  " << std::endl;
    mode_output << std::setw(5) << "#";
    for (unsigned int i=0; i<nvals; i++)
        mode_output
        << std::setw(10) << " "
        << std::setw(5) << "Mode "
        << std::setw(5) << i
        << std::setw(10) << " ";
    mode_output << std::endl;
    
    for (unsigned int i=0; i<nvals; i++)
    {
        const MAST::FlutterRootBase& root = this->get_root(i);
        
        if (root.if_nonphysical_root)
        {
            std::stringstream oss;
            oss << "**" << i;
            output << std::setw(5) << oss.str();
            mode_output << std::setw(5) << oss.str();
        }
        else
        {
            output << std::setw(5) << i;
            mode_output << std::setw(5) << i;
        }
        
        // flutter root details
        output
        << std::setw(15) << std::real(root.root)
        << std::setw(15) << std::imag(root.root)
        << std::setw(15) << root.g
        << std::setw(15) << root.V
        << std::setw(15) << root.omega;
        
        // now write the modal participation
        for (unsigned int j=0; j<nvals; j++)
            output << std::setw(15) << root.modal_participation(j);
        output << std::endl;
        
        // now write the flutter mode
        for (unsigned int j=0; j<nvals; j++)
        {
            mode_output
            << std::setw(13) << std::real(root.mode(j))
            << std::setw(4) << " "
            << std::setw(13) << std::imag(root.mode(j));
        }
        mode_output << std::endl;
    }
}



void
MAST::FlutterRootCrossoverBase::print(std::ostream &output) const
{
    const MAST::FlutterRootBase
    &lower = crossover_solutions.first->get_root(root_num),
    &upper = crossover_solutions.second->get_root(root_num);
    
    output
    << " Lower Root: " << std::endl
    << "    k : " << std::setw(15) << lower.k_ref << std::endl
    << "    g : " << std::setw(15) << lower.g << std::endl
    << "    V : " << std::setw(15) << lower.V << std::endl
    << "omega : " << std::setw(15) << lower.omega << std::endl
    << " Upper Root: " << std::endl
    << "    k : " << std::setw(15) << upper.k_ref << std::endl
    << "    g : " << std::setw(15) << upper.g << std::endl
    << "    V : " << std::setw(15) << upper.V << std::endl
    << "omega : " << std::setw(15) << upper.omega << std::endl;
    
    if (root)
        output
        << "Critical Root: " << std::endl
        << "    k : " << std::setw(15) << root->k_ref << std::endl
        << "    g : " << std::setw(15) << root->g << std::endl
        << "    V : " << std::setw(15) << root->V << std::endl
        << "omega : " << std::setw(15) << root->omega << std::endl;
    else
        output << "Critical root not yet calculated." << std::endl;
}
