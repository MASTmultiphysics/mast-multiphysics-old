//
//  ug_flutter_solver.cpp
//  MAST
//
//  Created by Manav Bhatia on 9/16/13.
//  Copyright (c) 2013 Manav Bhatia. All rights reserved.
//

// C++ includes
#include <iomanip>

// MAST includes
#include "Aeroelasticity//ug_flutter_solver.h"


void
UGFlutterSolver::UGFlutterRoot::init(const Real k, const Real b_ref,
                                     const Complex num,
                                     const Complex den,
                                     const ComplexMatrixX& Bmat,
                                     const ComplexVectorX& eig_vec)
{
    k_ref = k;

    if (std::abs(den) > 0.)
    {
        root = num/den;
        if (std::real(root) > 0.)
        {
            V     = sqrt(1./std::real(root));
            g     = std::imag(root)/std::real(root);
            omega = k_ref*V/b_ref;
            if_imag_V = false;
        }
        else
        {
            V     = 0.;
            g     = 0.;
            omega = 0.;
            if_imag_V = true;
        }
    }

    // calculate the modal participation vector
    const unsigned int nvals = (int)Bmat.rows();
    mode = eig_vec;
    ComplexVectorX k_q = Bmat * mode;
    RealVectorX modal_participation(nvals);
    for (unsigned int i=0; i<nvals; i++)
        modal_participation(i) =  std::abs(std::conj(mode(i)) * k_q(i));
    modal_participation *= (1./modal_participation.sum());
}




void
UGFlutterSolver::UGFlutterSolution::init (const Real kref, const Real bref,
                                          const LAPACK_ZGGEV& eig_sol)
{
    _k_ref = kref;
    // iterate over the roots and initialize the vector
    const ComplexMatrixX &B = eig_sol.B(), &VR = eig_sol.right_eigenvectors();
    const ComplexVectorX &num = eig_sol.alphas(), &den = eig_sol.betas();
    unsigned int nvals = (int)B.rows();
    
    for (unsigned int i=0; i<nvals; i++)
        _roots[i]->init(kref, bref, num(i), den(i), B, VR.col(i));
}




void
UGFlutterSolver::UGFlutterSolution::sort(const UGFlutterSolver::UGFlutterSolution& sol)
{
    const unsigned int nvals = this->n_roots();
    libmesh_assert_equal_to(nvals, sol.n_roots());
    
    for (unsigned int i=0; i<nvals-1; i++)
    {
        const UGFlutterSolver::UGFlutterRoot& r = sol.get_root(i);
        Real max_val = 0., val = 0.;
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



UGFlutterSolver::~UGFlutterSolver()
{
    std::map<Real, UGFlutterSolver::UGFlutterSolution*>::iterator it =
    _flutter_solutions.begin();
    
    for ( ; it != _flutter_solutions.end(); it++)
        delete it->second;
}



bool UGFlutterSolver::find_next_root()
{
    // active range in which to search
    std::pair<Real, Real> active_k_ref_range(0.,0.);
    
    // if this is the first solve, then the active range is same as the
    // k_ref range
    if (n_roots_found() == 0)
        active_k_ref_range = k_ref_range;
    else
        active_k_ref_range = std::pair<Real, Real>
        (k_ref_range.first, _previous_k_ref);
    
    // now march from the upper limit to the lower to find the roots
    Real current_k_ref = active_k_ref_range.second,
    delta_k_ref = (active_k_ref_range.second-active_k_ref_range.first)/n_k_divs;
    bool continue_finding = true, found_root = false;
    
    std::vector<Real> k_vals(n_k_divs+1);
    for (unsigned int i=0; i<n_k_divs+1; i++)
    {
        k_vals[i] = current_k_ref;
        current_k_ref -= delta_k_ref;
    }
    k_vals[n_k_divs] = active_k_ref_range.first; // to get around finite-precision arithmetic
    
    ComplexMatrixX m, k;
    
    UGFlutterSolver::UGFlutterSolution* prev_sol = NULL;
    
    for (unsigned int i=0; i<n_k_divs+1; i++)
    {
        current_k_ref = k_vals[i];
        initialize_matrices(current_k_ref, m, k);
        LAPACK_ZGGEV ges;
        ges.compute(m, k);
        ges.scale_eigenvectors_to_identity_innerproduct();
        
        // now insert the root
        std::pair<Real, UGFlutterSolver::UGFlutterSolution*>
        val(current_k_ref, new UGFlutterSolver::UGFlutterSolution((int)m.rows()));
        val.second->init(current_k_ref, flight_condition->ref_chord, ges);
        bool success = _flutter_solutions.insert(val).second;
        libmesh_assert (success);
        
        if (i > 0)
            val.second->sort(*prev_sol);
        prev_sol = val.second;
    }

    
    //print_roots(current_k_ref, ges);

    return found_root;
}



bool UGFlutterSolver::bisection_search(const std::pair<Real, Real>& active_k_ref_range)
{
    // to be programmed
    libmesh_assert(false);
}


void UGFlutterSolver::initialize_matrices(Real k_ref,
                                          ComplexMatrixX& m, // mass & aero
                                          ComplexMatrixX& k) // stiffness
{
    bool has_matrix = false;
    RealMatrixX mat_r;
    
    // stiffness matrix forms the rhs of the eigenvalue problem
    has_matrix = aero_structural_model->get_structural_stiffness_matrix(mat_r);
    k = mat_r.cast<Complex>();
    libmesh_assert(has_matrix);
    
    // combination of mass and aero matrix forms lhs of the eigenvalue problem
    has_matrix = aero_structural_model->get_structural_mass_matrix(mat_r);
    libmesh_assert(has_matrix);
    
    
    has_matrix = aero_structural_model->get_aero_operator_matrix(k_ref, m);
    libmesh_assert(has_matrix);
    
    m *= 0.5 * flight_condition->gas_property.rho;
    mat_r *= pow(k_ref/flight_condition->ref_chord, 2);
    m += mat_r.cast<Complex>();
}


void
UGFlutterSolver::evaluate_roots(const Real k_ref,
                                const LAPACK_ZGGEV& ges)
{
    const ComplexVectorX &num = ges.alphas(), &den = ges.betas();
    const ComplexMatrixX &vr = ges.right_eigenvectors(), &B = ges.B();
    
    unsigned int nvals = std::max(num.rows(), num.cols());
    
    
    Real g, V, omega;
    Complex eig;
    bool if_imag_v = false;
    
    _output << "k = " << k_ref << std::endl
    << std::setw(5) << "#"
    << std::setw(15) << "Re"
    << std::setw(15) << "Im"
    << std::setw(15) << "g"
    << std::setw(15) << "V"
    << std::setw(15) << "omega";
    
    // output the headers for flutter mode participation
    for (unsigned int i=0; i<nvals; i++)
        _output
        << std::setw(2) << " "
        << std::setw(5) << "Mode "
        << std::setw(5) << i
        << std::setw(3) << " ";
    _output << std::endl;
    
    // output the headers for flutter mode
    _mode_output << "k = " << k_ref << "  (Right Eigenvectors)  " << std::endl;
    _mode_output << std::setw(5) << "#";
    for (unsigned int i=0; i<nvals; i++)
        _mode_output
        << std::setw(10) << " "
        << std::setw(5) << "Mode "
        << std::setw(5) << i
        << std::setw(10) << " ";
    _mode_output << std::endl;
    
    for (unsigned int i=0; i<nvals; i++)
    {
        g = 0.; V= 0., omega = 0.;
        eig = 0.;
        if (std::abs(den(i)) > 0.)
        {
            eig = num(i)/den(i);
            if (std::real(eig) > 0.)
            {
                V     = sqrt(1./std::real(eig));
                g     = std::imag(eig)/std::real(eig);
                omega = k_ref*V/flight_condition->ref_chord;
            }
            else
                if_imag_v = true;
        }
        
        
        if (if_imag_v)
        {
            std::stringstream oss;
            oss << "**" << i;
            _output << std::setw(5) << oss.str();
            _mode_output << std::setw(5) << oss.str();
        }
        else
        {
            _output << std::setw(5) << i;
            _mode_output << std::setw(5) << i;
        }
        
        // flutter root details
        _output
        << std::setw(15) << std::real(eig)
        << std::setw(15) << std::imag(eig)
        << std::setw(15) << g
        << std::setw(15) << V
        << std::setw(15) << omega;
        
        // now write the modal participation
//        for (unsigned int j=0; j<nvals; j++)
//            _output << std::setw(15) << modal_participation(j);
//        _output << std::endl;
        
        // now write the flutter mode
        for (unsigned int j=0; j<nvals; j++)
        {
            _mode_output
            << std::setw(13) << std::real(vr(j,i))
            << std::setw(4) << " "
            << std::setw(13) << std::imag(vr(j,i));
        }
        _mode_output << std::endl;
        
        if_imag_v = false;
    }
    
}

