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

// C++ includes
#include <iomanip>

// MAST includes
#include "Aeroelasticity/ug_flutter_solver.h"



void
MAST::UGFlutterRoot::init(const Real k, const Real b_ref,
                          const Complex num,
                          const Complex den,
                          const ComplexMatrixX& Bmat,
                          const ComplexVectorX& evec_right,
                          const ComplexVectorX& evec_left)
{
    k_red = k;
    
    if (std::abs(den) > 0.)
    {
        root = num/den;
        if (std::real(root) > 0.)
        {
            V     = sqrt(1./std::real(root));
            g     = std::imag(root)/std::real(root);
            omega = k_red*V/b_ref;
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
    eig_vec_right = evec_right;
    eig_vec_left  = evec_left;
    ComplexVectorX k_q = Bmat * evec_right;
    modal_participation.resize(nvals, 1);
    for (unsigned int i=0; i<nvals; i++)
        modal_participation(i) =  std::abs(std::conj(evec_right(i)) * k_q(i));
    modal_participation *= (1./modal_participation.sum());
}





MAST::UGFlutterSolver::~UGFlutterSolver()
{ }



std::auto_ptr<MAST::FlutterRootBase>
MAST::UGFlutterSolver::build_flutter_root() {
    return std::auto_ptr<MAST::FlutterRootBase>(new MAST::UGFlutterRoot);
}



void MAST::UGFlutterSolver::_identify_crossover_points()
{
    // if the initial scanning has not been done, then do it now
    const Real tol = 1.0e-5, max_allowable_g = 0.75;
    
    const unsigned int nvals = _flutter_solutions.begin()->second->n_roots();
    // make sure that the solution has been generated
    libmesh_assert(nvals);
    
    //
    // for some cases some roots trail along the g=0 axis
    // and should not be considered as flutter. These are simply
    // modes where aerodynamics do not provide any damping.
    // These modes will not have a damping more than tolerance
    //
    
    std::vector<bool> modes_to_neglect(nvals);
    std::fill(modes_to_neglect.begin(),
              modes_to_neglect.end(), false);
    
    // look for the max g val for a mode, which will indicate if the
    // mode is undamped or not
    for (unsigned int i=0; i<nvals; i++) {
        std::map<Real, MAST::FlutterSolutionBase*>::const_iterator
        sol_it    = _flutter_solutions.begin(),
        sol_end   = _flutter_solutions.end();
        Real max_g_val = 0., val = 0.;
        for ( ; sol_it!=sol_end; sol_it++) {
            val = fabs(sol_it->second->get_root(i).g);
            if (val > max_g_val)
                max_g_val = val;
        }
        // check the maximum damping seen for this mode
        if (max_g_val < tol)
            modes_to_neglect[i] = true;
    }
    
    // identify the flutter cross-overs. For this, move from
    // higher k to lower k, and handle the k=0 cases separately
    {
        std::map<Real, MAST::FlutterSolutionBase*>::const_iterator
        sol_it    = _flutter_solutions.begin(), // first of the pair
        sol_end   = _flutter_solutions.end();
        if (sol_it == sol_end)
            return;
        
        // check if k=0 exists, and identify
        if (fabs(sol_it->second->ref_val()) < tol) { // k = 0
            
            // k=0 makes sense only for divergence roots. Do not use them
            // for crossover points if a finite damping was seen. Hence,
            // do nothing here, and move to the next iterator
            for (unsigned int i=0; i<nvals; i++) {
                if (!sol_it->second->get_root(i).if_nonphysical_root &&
                    fabs(sol_it->second->get_root(i).g) < tol) {
                    MAST::FlutterRootCrossoverBase* cross =
                    new MAST::FlutterRootCrossoverBase;
                    cross->crossover_solutions.first = sol_it->second;
                    cross->crossover_solutions.second = sol_it->second;
                    cross->root = &sol_it->second->get_root(i);
                    cross->root_num = i;
                    std::pair<Real, MAST::FlutterRootCrossoverBase*>
                    val( cross->root->V, cross);
                    _flutter_crossovers.insert(val);
                }
            }
        }
    }
    
    // now look for oscillatory roots crossover points in decreasing
    // order of k_red
    for (unsigned int i=0; i<nvals; i++) {
        std::map<Real, MAST::FlutterSolutionBase*>::const_reverse_iterator
        sol_rit    = _flutter_solutions.rbegin(), // first of the pair
        sol_ritp1    = _flutter_solutions.rbegin(), // first of the pair
        sol_rend   = _flutter_solutions.rend();
        if (sol_rit == sol_rend)
            return;
        
        sol_ritp1++; // increment for the next pair of results
        while (sol_ritp1 != sol_rend) {
            // do not use k_red = 0, or if the root is invalid
            if (sol_rit->second->get_root(i).if_nonphysical_root ||
                sol_ritp1->second->get_root(i).if_nonphysical_root ||
                fabs(sol_rit->second->ref_val()) < tol ||
                fabs(sol_ritp1->second->ref_val()) < tol ||
                fabs(sol_rit->second->get_root(i).g) > max_allowable_g ||
                fabs(sol_ritp1->second->get_root(i).g) > max_allowable_g) {
                // do nothing
            }
            else if (!modes_to_neglect[i]) { // look for the flutter roots
                MAST::FlutterSolutionBase *lower = sol_rit->second,
                *upper = sol_ritp1->second;
                
                if ((lower->get_root(i).g <= 0.) &&
                    (upper->get_root(i).g > 0.)) {
                    MAST::FlutterRootCrossoverBase* cross =
                    new MAST::FlutterRootCrossoverBase;
                    cross->crossover_solutions.first  = lower; // -ve g
                    cross->crossover_solutions.second = upper; // +ve g
                    cross->root_num = i;
                    std::pair<Real, MAST::FlutterRootCrossoverBase*>
                    val( lower->get_root(i).V, cross);
                    _flutter_crossovers.insert(val);
                }
                else if ((lower->get_root(i).g > 0.) &&
                         (upper->get_root(i).g <= 0.)) {
                    MAST::FlutterRootCrossoverBase* cross =
                    new MAST::FlutterRootCrossoverBase;
                    cross->crossover_solutions.first  = upper; // -ve g
                    cross->crossover_solutions.second = lower; // +ve g
                    cross->root_num = i;
                    std::pair<Real, MAST::FlutterRootCrossoverBase*>
                    val( upper->get_root(i).V, cross);
                    _flutter_crossovers.insert(val);
                }
            }
            
            // increment the pointers for next pair of roots
            sol_rit++;
            sol_ritp1++;
        }
    }
}



MAST::FlutterSolutionBase*
MAST::UGFlutterSolver::analyze(const Real k_red,
                               const Real v_ref,
                               const MAST::FlutterSolutionBase* prev_sol) {
    ComplexMatrixX m, k;
    
    libMesh::out
    << " ====================================================" << std::endl
    << "UG Solution for k_red = "
    << std::setw(10) << k_red << std::endl;
    
    initialize_matrices(k_red, v_ref, m, k);
    LAPACK_ZGGEV ges;
    ges.compute(m, k);
    ges.scale_eigenvectors_to_identity_innerproduct();

    // now insert the root
    std::pair<Real, MAST::FlutterSolutionBase*>
    val(k_red, new MAST::FrequencyDomainFlutterSolution(*this));
    bool success = _flutter_solutions.insert(val).second;
    libmesh_assert (success); // make sure that it was successfully added
    dynamic_cast<MAST::FrequencyDomainFlutterSolution*>(val.second)->init
    (k_red, flight_condition->ref_chord, ges);
    val.second->print(_output, _mode_output);
    
    if (prev_sol)
        val.second->sort(*prev_sol);
    
    libMesh::out
    << "Finished UG Solution" << std::endl
    << " ====================================================" << std::endl;
    
    
    return val.second;
}




void
MAST::UGFlutterSolver::calculate_sensitivity(MAST::FlutterRootBase& root,
                                             const libMesh::ParameterVector& params,
                                             const unsigned int i) {
    // make sure that the aero_structural_model is a valid pointer
    libmesh_assert(aero_structural_model);
    
    Complex eig = root.root, sens = 0., k_sens = 0., den = 0.;
    Real par_g_par_alpha = 0., par_g_par_kref = 0., par_k_par_alpha = 0.,
    V_sens=0.;
    
    // get the sensitivity of the matrices
    ComplexMatrixX mat_A, mat_B, mat_A_sens, mat_B_sens;
    ComplexVectorX v;

    // initialize the baseline matrices
    initialize_matrices(root.k_red, root.V, mat_A, mat_B);

    // calculate the eigenproblem sensitivity
    initialize_matrix_sensitivity_for_param(params, i,
                                            root.k_red,
                                            root.V,
                                            mat_A_sens,
                                            mat_B_sens);

    // the eigenproblem is     A x - lambda B x = 0
    // therefore, the denominator is obtained from the inner product of
    // x^T B x
    // sensitivity is
    //   -dlambda/dp x^T B x = - x^T (dA/dp - lambda dB/dp)
    // or
    //   dlambda/dp = [x^T (dA/dp - lambda dB/dp)]/(x^T B x)
    
    // now calculate the quotient for sensitivity
    // numerator =  ( dA/dp - lambda dB/dp)
    mat_B_sens *= -eig;
    mat_B_sens += mat_A_sens;
    v = mat_B_sens*root.eig_vec_right;
    den = root.eig_vec_left.dot(mat_B*root.eig_vec_right);
    sens = root.eig_vec_left.dot(v)/den;
    
    // now add the correction from sensitivity of g(k) = 0
    par_g_par_alpha =
    sens.imag()/eig.real() - eig.imag()/pow(eig.real(),2) * sens.real();

    
    // next we need the sensitivity of k_red before we can calculate
    // the sensitivity of flutter eigenvalue
    initialize_matrix_sensitivity_for_reduced_freq(root.k_red,
                                                   root.V,
                                                   mat_A_sens,
                                                   mat_B_sens);
    
    // now calculate the quotient for sensitivity wrt k_red
    // calculate numerator
    mat_B_sens *= -eig;
    mat_B_sens += mat_A_sens;
    v = mat_B_sens*root.eig_vec_right;
    k_sens = root.eig_vec_left.dot(v) / den;
    
    // use this to calculate the partial derivative of g wrt k_red
    par_g_par_kref =
    k_sens.imag()/eig.real() - eig.imag()/pow(eig.real(),2) * k_sens.real();
    
    // use this to calculate the sensitivity of k_red wrt alpha
    par_k_par_alpha = -par_g_par_alpha / par_g_par_kref;

    // finally add the correction to the flutter sensitivity
    sens += k_sens * par_k_par_alpha;

    // finally, the flutter speed sensitivity
    V_sens = -.5*sens.real()/pow(eig.real(), 1.5);
    
    // set value in the return root
    root.has_sensitivity_data = true;
    root.root_sens  = sens;
    root.k_red_sens = k_sens;
    root.V_sens     = V_sens;
}




void
MAST::UGFlutterSolver::initialize_matrices(const Real k_red,
                                           const Real v_ref,
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
    
    
    has_matrix =
    aero_structural_model->get_aero_operator_matrix(k_red, v_ref, m);
    libmesh_assert(has_matrix);
    
    m *= 0.5 * flight_condition->gas_property.rho;
    mat_r *= pow(k_red/flight_condition->ref_chord, 2);
    m += mat_r.cast<Complex>();
}


void
MAST::UGFlutterSolver::
initialize_matrix_sensitivity_for_param(const libMesh::ParameterVector& params,
                                        unsigned int p,
                                        const Real k_red,
                                        const Real v_ref,
                                        ComplexMatrixX& m, // mass & aero
                                        ComplexMatrixX& k) { // stiffness
    bool has_matrix = false;
    RealMatrixX mat_r;
    
    // stiffness matrix forms the rhs of the eigenvalue problem
    has_matrix =
    aero_structural_model->get_structural_stiffness_matrix_sensitivity(params,
                                                                       p,
                                                                       mat_r);
    k = mat_r.cast<Complex>();
    libmesh_assert(has_matrix);
    
    // combination of mass and aero matrix forms lhs of the eigenvalue problem
    has_matrix =
    aero_structural_model->get_structural_mass_matrix_sensitivity(params,
                                                                  p,
                                                                  mat_r);
    libmesh_assert(has_matrix);
    
    
    has_matrix =
    aero_structural_model->get_aero_operator_matrix_sensitivity(params,
                                                                p,
                                                                k_red,
                                                                v_ref,
                                                                m);
    libmesh_assert(has_matrix);
    
    m *= 0.5 * flight_condition->gas_property.rho;
    mat_r *= pow(k_red/flight_condition->ref_chord, 2);
    m += mat_r.cast<Complex>();
}



void
MAST::UGFlutterSolver::
initialize_matrix_sensitivity_for_reduced_freq(const Real k_red,
                                               const Real v_ref,
                                               ComplexMatrixX& m, // mass & aero
                                               ComplexMatrixX& k) { // stiffness
    bool has_matrix = false;
    RealMatrixX mat_r;

    // combination of mass and aero matrix forms lhs of the eigenvalue problem
    has_matrix =
    aero_structural_model->get_structural_mass_matrix(mat_r);
    libmesh_assert(has_matrix);
    
    
    has_matrix =
    aero_structural_model->get_aero_operator_matrix_sensitivity_for_reduced_freq(k_red,
                                                                                 v_ref,
                                                                                 m);
    libmesh_assert(has_matrix);
    
    m *= 0.5 * flight_condition->gas_property.rho;
    mat_r *= 2.*k_red*pow(1./flight_condition->ref_chord, 2);
    m += mat_r.cast<Complex>();
    
    k.setZero(m.rows(), m.cols());
}

