//
//  euler_elem_base.cpp
//  FESystem
//
//  Created by Manav Bhatia on 3/6/13.
//
//


// FESystem includes
#include "euler/euler_elem_base.h"


// C++ includes
#include <iomanip>

// Basic include files
#include "libmesh/getpot.h"
#include "libmesh/mesh.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/parameters.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;




PrimitiveSolution::PrimitiveSolution()
{
    this->zero();
}


void
PrimitiveSolution::zero()
{
    this->primitive_sol.zero();
    dimension = 0;
    cp = 0.; cv = 0.;
    rho = 0.; u1 = 0.; u2 = 0.; u3 = 0.; T = 0.; p = 0.; a = 0.; e_tot = 0.; k = 0.;
    entropy = 0.; mach = 0.;
}


void PrimitiveSolution::init(const unsigned int dim, const DenseVector<Real> &conservative_sol, const Real cp_val, const Real cv_val)
{
    dimension = dim;
    const unsigned int n1 = dim+2;
    cp = cp_val; cv = cv_val;
    const Real R = cp-cv, gamma = cp/cv;
    primitive_sol.resize(n1);
    
    rho = conservative_sol(0);
    primitive_sol(0) = rho;
    
    u1 = conservative_sol(1)/rho;
    primitive_sol(1) = u1;
    k = u1*u1;
    
    if (dim > 1)
    {
        u2 = conservative_sol(2)/rho;
        primitive_sol(2) = u2;
        k += u2*u2;
    }
    
    if (dim > 2)
    {
        u3 = conservative_sol(3)/rho;
        primitive_sol(3) = u3;
        k += u3*u3;
    }
    k *= 0.5;
    
    e_tot = conservative_sol(n1-1)/rho; // cv*T+k;

    T = (e_tot - k)/cv;
    primitive_sol(n1-1) = T;
    
    p = R*T*rho;
    a = sqrt(gamma*R*T);
    mach = sqrt(2.0*k)/a;
    entropy = log(p/pow(rho,gamma));
}


Real
PrimitiveSolution::c_pressure(const Real p0, const Real q0)
{
    return (p-p0)/q0;
}


void
PrimitiveSolution::print(std::ostream& out) const
{
    out << "Primitive Solution:" << std::endl;
    primitive_sol.print(out);
    std::cout
    << std::setw(15) <<  " rho: " << rho << std::endl
    << std::setw(15) <<  " u1: " << u1 << std::endl
    << std::setw(15) <<  " u2: " << u2 << std::endl
    << std::setw(15) <<  " u3: " << u3 << std::endl
    << std::setw(15) <<  " mach: " << mach << std::endl
    << std::setw(15) <<  " a: " << a << std::endl
    << std::setw(15) <<  " T: " << T << std::endl
    << std::setw(15) <<  " p: " << p << std::endl
    << std::setw(15) <<  " e_tot: " << e_tot << std::endl
    << std::setw(15) <<  " k: " << k << std::endl
    << std::setw(15) <<  " entropy: " << entropy << std::endl << std::endl;
}




template <typename ValType>
SmallPerturbationPrimitiveSolution<ValType>::SmallPerturbationPrimitiveSolution()
{
    this->zero();
}


template <typename ValType>
void SmallPerturbationPrimitiveSolution<ValType>::zero()
{
    perturb_primitive_sol.zero();
    drho = 0.; du1 = 0.; du2 = 0.; du3 = 0.; dT = 0.; dp = 0.; da = 0.; de_tot = 0.; dk = 0.;
    dentropy = 0.; dmach = 0.;
}


template <typename ValType>
void
SmallPerturbationPrimitiveSolution<ValType>::init(const PrimitiveSolution& sol, const DenseVector<ValType>& delta_sol)
{
    const unsigned int n1 = sol.dimension+2;
    const Real R = sol.cp-sol.cv, gamma = sol.cp/sol.cv;
    perturb_primitive_sol.resize(n1);
    
    drho = delta_sol(0);
    perturb_primitive_sol(0) = drho;
    
    du1 = (delta_sol(1) - drho * sol.u1)/sol.rho;
    perturb_primitive_sol(1) = du1;
    dk = sol.u1*du1;
    
    if (sol.dimension > 1)
    {
        du2 = (delta_sol(2) - drho * sol.u2)/sol.rho;
        perturb_primitive_sol(2) = du2;
        dk += sol.u2*du2;
    }
    
    if (sol.dimension > 2)
    {
        du3 = (delta_sol(3) - drho * sol.u3)/sol.rho;
        perturb_primitive_sol(3) = du3;
        dk += sol.u3*du3;
    }
    
    de_tot = (delta_sol(n1-1) - drho * sol.e_tot)/sol.rho;

    dT = (de_tot - dk)/sol.cv;
    perturb_primitive_sol(n1-1) = dT;
    
    dp = R*(dT*sol.rho + sol.T*drho);
    da = 0.5*sqrt(gamma*R/sol.T)*dT;
    dmach =  dk/sqrt(2.0*sol.k)/sol.a - sqrt(2.*sol.k)/pow(sol.a,2) * da;
    dentropy = (dp/pow(sol.rho,gamma) - gamma*sol.p/pow(sol.rho,gamma+1.)*drho) / (sol.p/pow(sol.rho,gamma)) ;
}


template <typename ValType>
void
SmallPerturbationPrimitiveSolution<ValType>::print(std::ostream& out) const
{
    out << "Small Perturbation Primitive Solution:" << std::endl;
    perturb_primitive_sol.print(out);
    std::cout
    << std::setw(15) <<  " drho: " << drho << std::endl
    << std::setw(15) <<  " du1: " << du1 << std::endl
    << std::setw(15) <<  " du2: " << du2 << std::endl
    << std::setw(15) <<  " du3: " << du3 << std::endl
    << std::setw(15) <<  " dmach: " << dmach << std::endl
    << std::setw(15) <<  " da: " << da << std::endl
    << std::setw(15) <<  " dT: " << dT << std::endl
    << std::setw(15) <<  " dp: " << dp << std::endl
    << std::setw(15) <<  " de_tot: " << de_tot << std::endl
    << std::setw(15) <<  " dk: " << dk << std::endl
    << std::setw(15) <<  " dentropy: " << dentropy << std::endl << std::endl;
}
    

template <typename ValType>
ValType
SmallPerturbationPrimitiveSolution<ValType>::c_pressure(const Real p0, const Real q0)
{
    
}
    




void EulerElemBase::init_data ()
{
    // Check the input file for Reynolds number, application type,
    // approximation type
    const Real pi = acos(-1.);

    GetPot infile("euler.in");
    aoa = infile("aoa",0.0);
    rho_inf = infile("rho",1.05);
    mach_inf = infile("mach",0.5);
    temp_inf = infile("temp",300.0);
    cp = infile("cp",1.003e3);
    cv = infile("cv",0.716e3);
    
    gamma = cp/cv;
    R = cp-cv;
    a_inf = sqrt(gamma*R*temp_inf);
    
    u1_inf = mach_inf*a_inf*cos(aoa*pi/180.0);
    u2_inf = mach_inf*a_inf*sin(aoa*pi/180.0);
    u3_inf = 0.0;
    q0_inf = 0.5*rho_inf*(u1_inf*u1_inf+u2_inf*u2_inf+u3_inf*u3_inf);
    p_inf = R*rho_inf*temp_inf;
}





void EulerElemBase::get_infinity_vars( DenseVector<Real>& vars_inf )
{
    Real k = 0.0;
    vars_inf(0) = rho_inf;
    vars_inf(1) = rho_inf*u1_inf; k += u1_inf*u1_inf;
    if (dim > 1)
    {
        vars_inf(2) = rho_inf*u2_inf;
        k += u2_inf*u2_inf;
    }
    if (dim > 2)
    {
        vars_inf(3) = rho_inf*u3_inf;
        k += u3_inf*u3_inf;
    }
    vars_inf(dim+2-1) = rho_inf*(cv*temp_inf+0.5*k);
}




void EulerElemBase::update_solution_at_quadrature_point( const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c, 
                                                        const bool if_elem_domain, const DenseVector<Real>& elem_solution,
                                                        DenseVector<Real>& conservative_sol, PrimitiveSolution& primitive_sol,
                                                        DenseMatrix<Real>& B_mat)
{
    conservative_sol.zero();
    B_mat.zero();
    
    FEBase* fe;

    for ( unsigned int i_var=0; i_var<c.n_vars(); i_var++ )
    {
        if (if_elem_domain)
            c.get_element_fe(vars[i_var], fe);
        else
            c.get_side_fe(vars[i_var], fe);
        
        const std::vector<std::vector<Real> >& phi = fe->get_phi();
        const unsigned int n_phi = phi.size();

        for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
            B_mat( i_var, i_var*n_phi+i_nd ) = phi[i_nd][qp];
    }
    
    B_mat.vector_mult( conservative_sol, elem_solution );
    
    primitive_sol.zero();
    primitive_sol.init(dim, conservative_sol, cp, cv);
}




void EulerElemBase::update_jacobian_at_quadrature_point( const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c, PrimitiveSolution& primitive_sol,
                                                        std::vector<DenseMatrix<Real> >& dB_mat,
                                                        std::vector<DenseMatrix<Real> >& Ai_advection,
                                                        DenseMatrix<Real>& A_entropy, DenseMatrix<Real>& A_inv_entropy )
{
    for ( unsigned int i_dim=0; i_dim<dim; i_dim++ )
        dB_mat[ i_dim ].zero();
    
    FEBase* fe;

    for ( unsigned int i_var=0; i_var<c.n_vars(); i_var++ )
    {
        c.get_element_fe(vars[0], fe); // assuming that all variables have same interpolation
        
        const std::vector<std::vector<RealVectorValue> >& dphi = fe->get_dphi();
        const unsigned int n_phi = dphi.size();

        for ( unsigned int i_dim=0; i_dim<dim; i_dim++ )
            for ( unsigned int i_nd=0; i_nd<n_phi; i_nd++ )
                dB_mat[i_dim]( i_var, i_var*n_phi+i_nd ) = dphi[i_nd][qp](i_dim);
    }
    
    for ( unsigned int i_dim=0; i_dim < dim; i_dim++)
        this->calculate_advection_flux_jacobian ( i_dim, primitive_sol, Ai_advection[i_dim] );
    
    this->calculate_entropy_variable_jacobian ( primitive_sol, A_entropy, A_inv_entropy );
}



void
EulerElemBase::calculate_advection_flux(const unsigned int calculate_dim,
                                        const PrimitiveSolution& sol,
                                        DenseVector<Real>& flux)
{
    const unsigned int n1 = 2 + dim;
    
    const Real rho = sol.rho,
    u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    p = sol.p,
    e_tot = sol.e_tot;
    
    flux.zero();
    
    // calculate the flux using given flow parameters at this point
    switch (calculate_dim)
    {
        case 0:
        {
            flux(0) =  rho * u1;
            flux(n1-1) =  u1 * (rho * e_tot + p);
            switch (dim)
            {
                case 3:
                    flux(3) =  rho * u1 * u3;
                case 2:
                    flux(2) =  rho * u1 * u2;
                case 1:
                    flux(1) =  rho * u1 * u1 + p;
            }
        }
            break;
            
        case 1:
        {
            flux(0) =  rho * u2;
            flux(n1-1) =  u2 * (rho * e_tot + p);
            switch (dim)
            {
                case 3:
                    flux(3) =  rho * u2 * u3;
                case 2:
                    flux(2) =  rho * u2 * u2 + p;
                case 1:
                    flux(1) =  rho * u2 * u1;
            }
        }
            break;
            
        case 2:
        {
            flux(0) =  rho * u3;
            flux(n1-1) =  u3 * (rho * e_tot + p);
            switch (dim)
            {
                case 3:
                    flux(3) =  rho * u3 * u3 + p;
                case 2:
                    flux(2) =  rho * u3 * u2;
                case 1:
                    flux(1) =  rho * u3 * u1;
            }
        }
            break;
    }
}




void
EulerElemBase::calculate_advection_flux_jacobian(const unsigned int calculate_dim,
                                                 const PrimitiveSolution& sol,
                                                 DenseMatrix<Real>& mat)
{
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector
    
    const unsigned int n1 = 2 + dim;
    //    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    //
    //    FESystemAssert4((s.first == n1) && (s.second == n1), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n1);
    //    FESystemAssert0(div_coord < dim, FESystem::Exception::InvalidValue);
    
    mat.zero();
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    e_tot = sol.e_tot,
    T = sol.T;
    
    switch (calculate_dim)
    {
        case 0:
        {
            switch (dim)
            {
                case 3:
                {
                    mat(1, 3) = -u3*R/cv;
                    
                    mat(3, 0) = -u1*u3;
                    mat(3, 1) = u3;
                    mat(3, 3) = u1;
                    
                    mat(n1-1, 3) = -u1*u3*R/cv;
                }
                    
                case 2:
                {
                    mat(1, 2) = -u2*R/cv;
                    
                    mat(2, 0) = -u1*u2;
                    mat(2, 1) = u2;
                    mat(2, 2) = u1;
                    
                    mat(n1-1, 2) = -u1*u2*R/cv;
                }
                    
                case 1:
                {
                    mat(0, 1) = 1.0; // d U / d (rho u1)
                    
                    mat(1, 0) = -u1*u1+R*k/cv;
                    mat(1, 1) = u1*(2.0-R/cv);
                    mat(1, n1-1) = R/cv;
                    
                    mat(n1-1, 0) = u1*(R*(-e_tot+2.0*k)-e_tot*cv)/cv;
                    mat(n1-1, 1) = e_tot+R*(T-u1*u1/cv);
                    mat(n1-1, n1-1) = u1*gamma;
                }
                    break;
            }
        }
            break;
            
        case 1:
        {
            switch (dim)
            {
                case 3:
                {
                    mat(2, 3) = -u3*R/cv;
                    
                    mat(3, 0) = -u2*u3;
                    mat(3, 2) = u3;
                    mat(3, 3) = u2;
                    
                    mat(n1-1, 3) = -u2*u3*R/cv;
                }
                    
                case 2:
                {
                    mat(0, 2) = 1.0; // d U / d (rho u2)
                    
                    mat(1, 0) = -u1*u2;
                    mat(1, 1) = u2;
                    mat(1, 2) = u1;
                    
                    mat(2, 0) = -u2*u2+R*k/cv;
                    mat(2, 1) = -u1*R/cv;
                    mat(2, 2) = u2*(2.0-R/cv);
                    mat(2, n1-1) = R/cv;
                    
                    mat(n1-1, 0) = u2*(R*(-e_tot+2.0*k)-e_tot*cv)/cv;
                    mat(n1-1, 1) = -u1*u2*R/cv;
                    mat(n1-1, 2) = e_tot+R*(T-u2*u2/cv);
                    mat(n1-1, n1-1) = u2*gamma;
                }
                    break;
                    
                case 1:
                    // if second coordinate divergence is being asked for, then the element is atleast 2D
                    libmesh_assert_msg(false, "invalid dim");
                    break;
            }
        }
            break;
            
        case 2:
        {
            mat(0, 3) = 1.0; // d U / d (rho u3)
            
            mat(1, 0) = -u1*u3;
            mat(1, 1) = u3;
            mat(1, 3) = u1;
            
            mat(2, 0) = -u2*u3;
            mat(2, 2) = u3;
            mat(2, 3) = u2;
            
            mat(3, 0) = -u3*u3+R*k/cv;
            mat(3, 1) = -u1*R/cv;
            mat(3, 2) = -u2*R/cv;
            mat(3, 3) = u3*(2.0-R/cv);
            mat(3, n1-1) = R/cv;
            
            mat(n1-1, 0) = u3*(R*(-e_tot+2.0*k)-e_tot*cv)/cv;
            mat(n1-1, 1) = -u1*u3*R/cv;
            mat(n1-1, 2) = -u2*u3*R/cv;
            mat(n1-1, 3) = e_tot+R*(T-u3*u3/cv);
            mat(n1-1, n1-1) = u3*gamma;
        }
            break;
            
        default:
            libmesh_assert_msg(false, "invalid dim");
            break;
    }
}


void
EulerElemBase::calculate_advection_left_eigenvector_and_inverse_for_normal(const PrimitiveSolution& sol,
                                                                           const Point& normal, DenseMatrix<Real>& eig_vals,
                                                                           DenseMatrix<Real>& l_eig_mat, DenseMatrix<Real>& l_eig_mat_inv_tr)
{
    const unsigned int n1 = 2 + dim;
    
    //    const std::pair<FESystemUInt, FESystemUInt> s_eig_val = eig_vals.getSize(), s_l_eig_mat = l_eig_mat.getSize(), s_l_eig_mat_inv_tr = l_eig_mat_inv_tr.getSize();
    //
    //    FESystemAssert4((s_eig_val.first == n1) && (s_eig_val.first == n1), FESystem::Numerics::MatrixSizeMismatch, s_eig_val.first, s_eig_val.second, n1, n1);
    //    FESystemAssert4((s_l_eig_mat.first == n1) && (s_l_eig_mat.first == n1), FESystem::Numerics::MatrixSizeMismatch, s_l_eig_mat.first, s_l_eig_mat.second, n1, n1);
    //    FESystemAssert4((s_l_eig_mat_inv_tr.first == n1) && (s_l_eig_mat_inv_tr.first == n1), FESystem::Numerics::MatrixSizeMismatch, s_l_eig_mat_inv_tr.first, s_l_eig_mat_inv_tr.second, n1, n1);
    
    eig_vals.zero(); l_eig_mat.zero(); l_eig_mat_inv_tr.zero();
    
    Real nx=0., ny=0., nz=0., u=0.;
    unsigned int dim_for_eig_vec=100; // initializing with arbitrarily high value
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    T = sol.T,
    a = sol.a;
    
    // initialize the values
    switch (dim)
    {
        case 3:
        {
            nz = normal(2);
            u += u3*nz;
        }
            
        case 2:
        {
            ny = normal(1);
            u += u2*ny;
        }
            
        case 1:
        {
            nx = normal(0);
            u += u1*nx;
        }
    }
    
    // select the largest value of surface normal component, so that the appropriate section can be chosen
    if ((fabs(nx)>=fabs(ny)) && (fabs(nx)>=fabs(nz)))
        dim_for_eig_vec = 1;
    else if ((fabs(ny)>=fabs(nx)) && (fabs(ny)>=fabs(nz)))
        dim_for_eig_vec = 2;
    else if ((fabs(nz)>=fabs(nx)) && (fabs(nz)>=fabs(ny)))
        dim_for_eig_vec = 3;
    
    // set eigenvalues
    switch (dim)
    {
        case 3:
            eig_vals(2, 2) = u;
        case 2:
            eig_vals(1, 1) = u;
        case 1:
        {
            eig_vals(0, 0) = u;
            eig_vals(n1-2, n1-2) = u-a;
            eig_vals(n1-1, n1-1) = u+a;
        }
    }
    
    // set last two columns of the eigenvector matrices, and all columns of the eigenvector inverse.
    // Note that the dim_for_eig_vec column of the inverse matrix will be overwritten by the appropriate matrix
    switch (dim)
    {
        case 3:
        {
            // for u-a
            l_eig_mat(3, n1-2) = -u3-nz*a/(gamma-1.0);
            
            // for u+a
            l_eig_mat(3, n1-1) = -u3+nz*a/(gamma-1.0);
            
            // for u
            l_eig_mat_inv_tr(3, 0) = (gamma-1.0)*u1*u3/(a*a)-nx*nz;
            
            // for u
            l_eig_mat_inv_tr(3, 1) = (gamma-1.0)*u2*u3/(a*a)-ny*nz;
            
            // for u
            l_eig_mat_inv_tr(0, 2) = (gamma-1.0)*u3/(a*a);
            l_eig_mat_inv_tr(1, 2) = (gamma-1.0)*u3*u1/(a*a)-nz*nx;
            l_eig_mat_inv_tr(2, 2) = (gamma-1.0)*u3*u2/(a*a)-nz*ny;
            l_eig_mat_inv_tr(3, 2) = (gamma-1.0)*u3*u3/(a*a)+1.0-nz*nz;
            l_eig_mat_inv_tr(n1-1, 2) = (gamma-1.0)*k*u3/(a*a)+u3-nz*u;
            
            // overwrite for u
            l_eig_mat_inv_tr(3, dim_for_eig_vec-1) = u3;
            
            // for u-a
            l_eig_mat_inv_tr(3, n1-2) = 0.5*(gamma-1.0)*(-nz+u3/a)/a;
            
            // for u+a
            l_eig_mat_inv_tr(3, n1-1) = 0.5*(gamma-1.0)*(nz+u3/a)/a;
            
        }
            
        case 2:
        {
            // for u-a
            l_eig_mat(2, n1-2) = -u2-ny*a/(gamma-1.0);
            
            // for u+a
            l_eig_mat(2, n1-1) = -u2+ny*a/(gamma-1.0);
            
            // for u
            l_eig_mat_inv_tr(2, 0) = (gamma-1.0)*u1*u2/(a*a)-nx*ny;
            
            // for u
            l_eig_mat_inv_tr(0, 1) = (gamma-1.0)*u2/(a*a);
            l_eig_mat_inv_tr(1, 1) = (gamma-1.0)*u2*u1/(a*a)-ny*nx;
            l_eig_mat_inv_tr(2, 1) = (gamma-1.0)*u2*u2/(a*a)+1.0-ny*ny;
            l_eig_mat_inv_tr(n1-1, 1) = (gamma-1.0)*k*u2/(a*a)+u2-ny*u;
            
            // overwrite for u
            l_eig_mat_inv_tr(2, dim_for_eig_vec-1) = u2;
            
            // for u-a
            l_eig_mat_inv_tr(2, n1-2) = 0.5*(gamma-1.0)*(-ny+u2/a)/a;
            
            // for u+a
            l_eig_mat_inv_tr(2, n1-1) = 0.5*(gamma-1.0)*(ny+u2/a)/a;
        }
            
        case 1:
        {
            // for u-a
            l_eig_mat(0, n1-2) = k+u*a/(gamma-1.0);
            l_eig_mat(1, n1-2) = -u1-nx*a/(gamma-1.0);
            l_eig_mat(n1-1, n1-2) = 1.0;
            
            // for u+a
            l_eig_mat(0, n1-1) = k-u*a/(gamma-1.0);
            l_eig_mat(1, n1-1) = -u1+nx*a/(gamma-1.0);
            l_eig_mat(n1-1, n1-1) = 1.0;
            
            // for u
            l_eig_mat_inv_tr(0, 0) = (gamma-1.0)*u1/(a*a);
            l_eig_mat_inv_tr(1, 0) = (gamma-1.0)*u1*u1/(a*a)+1.0-nx*nx;
            l_eig_mat_inv_tr(n1-1, 0) = (gamma-1.0)*k*u1/(a*a)+u1-nx*u;
            
            // overwrite for u
            l_eig_mat_inv_tr(0, dim_for_eig_vec-1) = 1.0;
            l_eig_mat_inv_tr(1, dim_for_eig_vec-1) = u1;
            l_eig_mat_inv_tr(n1-1, dim_for_eig_vec-1) = k;
            l_eig_mat_inv_tr.scale_column(dim_for_eig_vec-1, -1.0/(cv*T*gamma));
            
            // for u-a
            l_eig_mat_inv_tr(0, n1-2) = 1.0/(2.0*cv*T*gamma);
            l_eig_mat_inv_tr(1, n1-2) = 0.5*(gamma-1.0)*(-nx+u1/a)/a;
            l_eig_mat_inv_tr(n1-1, n1-2) = 0.5*(1.0+(gamma-1.0)*(-u+k/a)/a);
            
            // for u+a
            l_eig_mat_inv_tr(0, n1-1) = 1.0/(2.0*cv*T*gamma);
            l_eig_mat_inv_tr(1, n1-1) = 0.5*(gamma-1.0)*(nx+u1/a)/a;
            l_eig_mat_inv_tr(n1-1, n1-1) = 0.5*(1.0+(gamma-1.0)*(u+k/a)/a);
            
        }
    }
    
    
    switch (dim_for_eig_vec)
    {
        case 1:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            switch (dim)
            {
                case 3:
                {
                    // for u
                    l_eig_mat(0, 2) = nz*u1/nx-u3;
                    l_eig_mat(1, 2) = -nz/nx;
                    l_eig_mat(3, 2) = 1.0;
                }
                    
                case 2:
                {
                    // for u
                    l_eig_mat(0, 1) = ny*u1/nx-u2;
                    l_eig_mat(1, 1) = -ny/nx;
                    l_eig_mat(2, 1) = 1.0;
                }
                    
                case 1:
                {
                    // for u
                    l_eig_mat(0, 0) = u*u1/nx-(cv*T*gamma+k);
                    l_eig_mat(1, 0) = -u/nx;
                    l_eig_mat(n1-1, 0) = 1.0;
                }
            }
        }
            break;
            
        case 2:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            switch (dim)
            {
                case 3:
                {
                    // for u
                    l_eig_mat(0, 2) = nz*u2/ny-u3;
                    l_eig_mat(2, 2) = -nz/ny;
                    l_eig_mat(3, 2) = 1.0;
                }
                    
                case 2:
                case 1:
                {
                    // for u
                    l_eig_mat(0, 0) = nx*u2/ny-u1;
                    l_eig_mat(1, 0) = 1.0;
                    l_eig_mat(2, 0) = -nx/ny;
                    
                    // for u
                    l_eig_mat(0, 1) = u*u2/ny-(cv*T*gamma+k);
                    l_eig_mat(2, 1) = -u/ny;
                    l_eig_mat(n1-1, 1) = 1.0;
                }
            }
        }
            break;
            
        case 3:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            
            // for u
            l_eig_mat(0, 0) = nx*u3/nz-u1;
            l_eig_mat(1, 0) = 1.0;
            l_eig_mat(3, 0) = -nx/nz;
            
            // for u
            l_eig_mat(0, 1) = ny*u3/nz-u2;
            l_eig_mat(2, 1) = 1.0;
            l_eig_mat(3, 1) = -ny/nz;
            
            // for u
            l_eig_mat(0, 2) = u*u3/nz-(cv*T*gamma+k);
            l_eig_mat(3, 2) = -u/nz;
            l_eig_mat(n1-1, 2) = 1.0;
        }
            break;
    }
}



void
EulerElemBase::calculate_advection_flux_jacobian_for_moving_solid_wall_boundary(const PrimitiveSolution& sol,
                                                                                const Real xi_ni, const Point& nvec, DenseMatrix<Real>& mat)
{
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector
    
    const unsigned int n1 = 2 + dim;
    //    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    //
    //    FESystemAssert4((s.first == n1) && (s.second == n1), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n1);
    
    mat.zero();
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k;
    
    switch (dim)
    {
        case 3:
        {
            mat(1, 3) = -R/cv*nvec(0)*u3;
            
            mat(2, 3) = -R/cv*nvec(1)*u3;
            
            mat(3, 0) = nvec(2)*R/cv*k;
            mat(3, 1) = -R/cv*nvec(2)*u1;
            mat(3, 2) = -R/cv*nvec(2)*u2;
            mat(3, 3) = xi_ni-R/cv*nvec(2)*u3;
            mat(3, n1-1) = R/cv*nvec(2);
            
            mat(n1-1, 3) = xi_ni*R/cv*u3;
        }
            
        case 2:
        {
            mat(1, 2) = -R/cv*nvec(0)*u2;
            
            mat(2, 0) = nvec(1)*R/cv*k;
            mat(2, 1) = -R/cv*nvec(1)*u1;
            mat(2, 2) = xi_ni-R/cv*nvec(1)*u2;
            mat(2, n1-1) = R/cv*nvec(1);
            
            mat(n1-1, 2) = xi_ni*R/cv*u2;
        }
            
        case 1:
        {
            mat(0, 1) = xi_ni; // d U / d (rho u1)
            
            mat(1, 0) = nvec(0)*R/cv*k;
            mat(1, 1) = xi_ni-R/cv*nvec(0)*u1;
            mat(1, n1-1) = R/cv*nvec(0);
            
            mat(n1-1, 0) = xi_ni*R/cv*k;
            mat(n1-1, 1) = xi_ni*R/cv*u1;
            mat(n1-1, n1-1) = xi_ni*(R+cv)/cv;
        }
            break;
    }
}




void
EulerElemBase::calculate_entropy_variable_jacobian(const PrimitiveSolution& sol,
                                                   DenseMatrix<Real>& dUdV, DenseMatrix<Real>& dVdU)
{
    // calculates dU/dV where V is the Entropy variable vector
    
    // calculate A0 = d U / d Y , where U = conservation variable vector, and Y = unknown variable vector
    // note that for conservation variables as the unknown, this is an identity matrix
    
    const unsigned int n1 = 2 + dim;
    const Real rho = sol.rho,
    u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    T = sol.T,
    e_tot = sol.e_tot;
    
    dUdV.zero(); dVdU.zero();
    
    // du/dv
    switch (dim)
    {
        case 3:
        {
            dUdV(0, 3) = u3;
            
            dUdV(1, 3) = u1*u3;
            
            dUdV(2, 3) = u2*u3;
            
            dUdV(3, 0) = dUdV(0, 3);
            dUdV(3, 1) = dUdV(1, 3);
            dUdV(3, 2) = dUdV(2, 3);
            dUdV(3, 3) = u3*u3+cv*T*(gamma-1.0);
            dUdV(3, n1-1) = u3*(cv*T*gamma+k);
            
            dUdV(n1-1, 3) = dUdV(3, n1-1);
        }
            
        case 2:
        {
            dUdV(0, 2) = u2;
            
            dUdV(1, 2) = u1*u2;
            
            dUdV(2, 0) = dUdV(0, 2);
            dUdV(2, 1) = dUdV(1, 2);
            dUdV(2, 2) = u2*u2+cv*T*(gamma-1.0);
            dUdV(2, n1-1) = u2*(cv*T*gamma+k);
            
            dUdV(n1-1, 2) = dUdV(2, n1-1);
        }
            
        case 1:
        {
            dUdV(0, 0) = 1.0;
            dUdV(0, 1) = u1;
            dUdV(0, n1-1) = e_tot;
            
            dUdV(1, 0) = dUdV(0, 1);
            dUdV(1, 1) = u1*u1+cv*T*(gamma-1.0);
            dUdV(1, n1-1) = u1*(cv*T*gamma+k);
            
            dUdV(n1-1, 0) = dUdV(0, n1-1);
            dUdV(n1-1, 1) = dUdV(1, n1-1);
            dUdV(n1-1, n1-1) = k*k+gamma*cv*T*(cv*T+2*k);
            
        }
            break;
            
        default:
            break;
    }
    
    dUdV.scale(rho/(gamma-1.0));
    
    
    // dv/du
    switch (dim)
    {
        case 3:
        {
            dVdU(0, 3) = -u3*k;
            
            dVdU(1, 3) = u1*u3;
            
            dVdU(2, 3) = u2*u3;
            
            dVdU(3, 0) = dVdU(0, 3);
            dVdU(3, 1) = dVdU(1, 3);
            dVdU(3, 2) = dVdU(2, 3);
            dVdU(3, 3) = u3*u3+cv*T;
            dVdU(3, n1-1) = -u3;
            
            dVdU(n1-1, 3) = dVdU(3, n1-1);
        }
            
        case 2:
        {
            dVdU(0, 2) = -u2*k;
            
            dVdU(1, 2) = u1*u2;
            
            dVdU(2, 0) = dVdU(0, 2);
            dVdU(2, 1) = dVdU(1, 2);
            dVdU(2, 2) = u2*u2+cv*T;
            dVdU(2, n1-1) = -u2;
            
            dVdU(n1-1, 2) = dVdU(2, n1-1);
        }
            
        case 1:
        {
            dVdU(0, 0) = k*k+cv*cv*T*T*gamma;
            dVdU(0, 1) = -u1*k;
            dVdU(0, n1-1) = -e_tot+2.0*k;
            
            dVdU(1, 0) = dVdU(0, 1);
            dVdU(1, 1) = u1*u1+cv*T;
            dVdU(1, n1-1) = -u1;
            
            dVdU(n1-1, 0) = dVdU(0, n1-1);
            dVdU(n1-1, 1) = dVdU(1, n1-1);
            dVdU(n1-1, n1-1) = 1.0;
        }
            break;
            
        default:
            break;
    }
    
    dVdU.scale(1.0/(rho*cv*cv*T*T));
}




void
EulerElemBase::calculate_artificial_diffusion_operator(const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c,
                                                       const PrimitiveSolution& sol,
                                                       DenseMatrix<Real>& streamline_operator)
{
    const unsigned int n1 = 2 + dim;
    
    FEBase* fe;
    c.get_element_fe(vars[0], fe);
    std::vector<std::vector<RealVectorValue> > dphi = fe->get_dphi(); // assuming that all variables have the same interpolation
    
    DenseVector<Real> u, dN;
    u.resize(dim); dN.resize(dim);
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    a = sol.a,
    dt = c.get_deltat_value();
    
    streamline_operator.zero();
    
    // calculate the gradients
    switch (dim)
    {
        case 3:
            u(2) = u3;
            
        case 2:
            u(1) = u2;
            
        case 1:
            u(0) = u1;
            break;
            
        default:
            break;
    }
    
    // calculate the dot product of velocity times gradient of shape function
    Real h = 0, u_val = u.l2_norm(), tau_rho, tau_m, tau_e;
    u.scale(1.0/u_val);
    
    for (unsigned int i_nodes=0; i_nodes<dphi.size(); i_nodes++)
    {
        // set value of shape function gradient
        for (unsigned int i_dim=0; i_dim<dim; i_dim++)
            dN(i_dim) = dphi[i_nodes][qp](i_dim);
        
        h += fabs(dN.dot(u));
    }
    
    h = 2.0/h;
    
    // now set the value of streamwise dissipation
    tau_rho = 1.0/sqrt(pow(2.0/dt, 2)+ pow(2.0/h*(u_val+a), 2));
    tau_m = 1.0/sqrt(pow(2.0/dt, 2)+ pow(2.0/h*(u_val+a), 2));
    tau_e = 1.0/sqrt(pow(2.0/dt, 2)+ pow(2.0/h*(u_val+a), 2));
    
    streamline_operator(0, 0) = tau_rho;
    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        streamline_operator(1+i_dim, 1+i_dim) = tau_m;
    streamline_operator(n1-1, n1-1) = tau_e;
}


void
EulerElemBase::calculate_dxidX (const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c,
                                DenseMatrix<Real>& dxi_dX)
{
    // initialize dxi_dX
    dxi_dX.zero();
    Real val=0.;
    FEBase* fe;
    c.get_element_fe(vars[0], fe); // assuming that all elements have the same interpolation fe
    
    for (unsigned int i_dim=0; i_dim<dim; i_dim++)
        for (unsigned int j_dim=0; j_dim<dim; j_dim++)
        {
            switch (i_dim)
            {
                case 0:
                {
                    switch (j_dim)
                    {
                        case 0:
                            val = fe->get_dxidx()[qp];
                            break;
                        case 1:
                            val = fe->get_dxidy()[qp];
                            break;
                        case 2:
                            val = fe->get_dxidz()[qp];
                            break;
                    }
                }
                    break;
                case 1:
                {
                    switch (j_dim)
                    {
                        case 0:
                            val = fe->get_detadx()[qp];
                            break;
                        case 1:
                            val = fe->get_detady()[qp];
                            break;
                        case 2:
                            val = fe->get_detadz()[qp];
                            break;
                    }
                }
                    break;
                case 2:
                {
                    switch (j_dim)
                    {
                        case 0:
                            val = fe->get_dzetadx()[qp];
                            break;
                        case 1:
                            val = fe->get_dzetady()[qp];
                            break;
                        case 2:
                            val = fe->get_dzetadz()[qp];
                            break;
                    }
                }
                    break;
            }
            dxi_dX(i_dim, j_dim) = val;
        }
}


void EulerElemBase::calculate_differential_operator_matrix(const std::vector<unsigned int>& vars, const unsigned int qp, FEMContext& c, const DenseVector<Real>& elem_solution, const PrimitiveSolution& sol,
                                                           const DenseMatrix<Real>& B_mat, const std::vector<DenseMatrix<Real> >& dB_mat,
                                                           const std::vector<DenseMatrix<Real> >& Ai_advection, const DenseMatrix<Real>& Ai_Bi_advection,
                                                           const DenseMatrix<Real>& A_inv_entropy,
                                                           DenseMatrix<Real>& LS_operator, Real& discontinuity_val )
{
    const unsigned int n1 = 2 + dim, n2 = B_mat.n();
    
    std::vector<DenseVector<Real> > diff_vec(3);
    DenseMatrix<Real> tmp_mat, tmp_mat_n1n1, dxi_dX;
    DenseVector<Real> vec1, vec2;
    tmp_mat.resize(n1, n2); tmp_mat_n1n1.resize(n1, n1); dxi_dX.resize(dim, dim);
    vec1.resize(n1); vec2.resize(n1);
    for (unsigned int i=0; i<dim; i++) diff_vec[i].resize(n1);
    
    this->calculate_dxidX (vars, qp, c, dxi_dX);
    
    // contribution of unsteady term
    LS_operator.zero();
    
    vec2.zero();
    
    // contribution of advection flux term
    for (unsigned int i=0; i<dim; i++)
    {
        Ai_advection[i].get_transpose(tmp_mat);
        tmp_mat.right_multiply(dB_mat[i]);
        LS_operator.add(1.0, tmp_mat);  // A_i^T dB/dx_i
        
        dB_mat[i].vector_mult(diff_vec[i], elem_solution);
        
        Ai_advection[i].vector_mult(vec1, diff_vec[i]);
        
        vec2.add(1.0, vec1); // sum A_i dU/dx_i
    }
    
    // add the velocity and calculate the numerator of the discontinuity capturing term coefficient
    //vec2 += c.elem_solution; // add velocity TODO: how to get the velocity for all calls to this method
    A_inv_entropy.vector_mult(vec1, vec2);
    discontinuity_val = vec1.dot(vec2); // this is the numerator term
    
    // now evaluate the dissipation factor for the discontinuity capturing term
    // this is the denominator term
    
    Real val1 = 0.0;
    for (unsigned int i=0; i<dim; i++)
    {
        vec1.zero();
        tmp_mat.zero();
        
        for (unsigned int j=0; j<dim; j++)
            vec1.add(dxi_dX(i, j), diff_vec[j]);
        
        // calculate the value of denominator
        A_inv_entropy.vector_mult(vec2, vec1);
        val1 += vec1.dot(vec2);
    }
    
    //    // now calculate the discontinuity capturing operator
    if ((fabs(val1) > 0.0) &&  (fabs(discontinuity_val) > 0.0))
        discontinuity_val = sqrt(discontinuity_val/val1);
    else
        discontinuity_val = 0.0;
    
    // scale the LS matrix with the correct factor
    this->calculate_artificial_diffusion_operator(vars, qp, c, sol, tmp_mat_n1n1);
    LS_operator.left_multiply(tmp_mat_n1n1);
}



template class SmallPerturbationPrimitiveSolution<Number>;
