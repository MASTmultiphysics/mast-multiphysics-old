

// MAST includes
#include "FluidElems/fluid_elem_base.h"


void
FluidElemBase::
calculate_advection_flux_jacobian_rho_derivative(const unsigned int calculate_dim,
                                                 const PrimitiveSolution& sol,
                                                 DenseRealMatrix& mat)
{
    // calculate dAi/drho = d^2 F_adv / d x_i drho, where F_adv is the Euler advection flux vector
    
    mat.zero();
}


void
FluidElemBase::
calculate_advection_flux_jacobian_u1_derivative(const unsigned int calculate_dim,
                                                const PrimitiveSolution& sol,
                                                DenseRealMatrix& mat)
{
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector
    
    const unsigned int n1 = 2 + dim;
    
    mat.zero();
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    T = sol.T,
    gamma = flight_condition->gas_property.gamma,
    R = flight_condition->gas_property.R,
    cv = flight_condition->gas_property.cv;
    
    switch (calculate_dim)
    {
        case 0:
        {
            switch (dim)
            {
                case 3:
                {
                    mat(3, 0) = -u3;
                    mat(3, 3) = 1.;
                    
                    mat(n1-1, 3) = -u3*R/cv;
                }
                    
                case 2:
                {
                    mat(2, 0) = -u2;
                    mat(2, 2) = 1.;
                    
                    mat(n1-1, 2) = -u2*R/cv;
                }
                    
                case 1:
                {
                    mat(1, 0) = (-2+R/cv)*u1;
                    mat(1, 1) = (2.0-R/cv);
                    
                    mat(n1-1, 0) = (-2*pow(cv,2)*T + R*(3*pow(u1,2) + pow(u2,2) + pow(u3,2)) -
                                    cv*(2*R*T + 3*pow(u1,2) + pow(u2,2) + pow(u3,2)))/(2.*cv);
                    mat(n1-1, 1) = u1 - (2*R*u1)/cv;
                    mat(n1-1, n1-1) = gamma;
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
                }
                    
                case 2:
                {
                    mat(1, 0) = -u2;
                    mat(1, 2) = 1.;
                    
                    mat(2, 0) = R*u1/cv;
                    mat(2, 1) = -R/cv;
                    
                    mat(n1-1, 0) = ((-cv + R)*u1*u2)/cv;
                    mat(n1-1, 1) = -u2*R/cv;
                    mat(n1-1, 2) = u1;
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
            mat(1, 0) = -u3;
            mat(1, 3) = 1.;
            
            mat(3, 0) = R*u1/cv;
            mat(3, 1) = -R/cv;
            
            mat(n1-1, 0) = ((-cv + R)*u1*u3)/cv;
            mat(n1-1, 1) = -u3*R/cv;
            mat(n1-1, 3) = u1;
        }
            break;
            
        default:
            libmesh_assert_msg(false, "invalid dim");
            break;
    }
}




void
FluidElemBase::
calculate_advection_flux_jacobian_u2_derivative(const unsigned int calculate_dim,
                                                const PrimitiveSolution& sol,
                                                DenseRealMatrix& mat)
{
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector
    
    const unsigned int n1 = 2 + dim;
    
    mat.zero();
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    T = sol.T,
    gamma = flight_condition->gas_property.gamma,
    R = flight_condition->gas_property.R,
    cv = flight_condition->gas_property.cv;
    
    switch (calculate_dim)
    {
        case 0:
        {
            switch (dim)
            {
                case 3:
                {
                }
                    
                case 2:
                {
                    mat(1, 2) = -R/cv;
                    
                    mat(2, 0) = -u1;
                    mat(2, 1) = 1.;
                    
                    mat(n1-1, 2) = -u1*R/cv;
                }
                    
                case 1:
                {
                    mat(1, 0) = R*u2/cv;
                    
                    mat(n1-1, 0) = ((-cv + R)*u1*u2)/cv;
                    mat(n1-1, 1) = u2;
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
                    mat(3, 0) = -u3;
                    mat(3, 3) = 1.;
                    
                    mat(n1-1, 3) = -u3*R/cv;
                }
                    
                case 2:
                {
                    mat(1, 0) = -u1;
                    mat(1, 1) = 1.;
                    
                    mat(2, 0) = (-2+R/cv)*u2;
                    mat(2, 2) = (2.0-R/cv);
                    
                    mat(n1-1, 0) = (-2*pow(cv,2)*T + R*(pow(u1,2) + 3*pow(u2,2) + pow(u3,2)) -
                                    cv*(2*R*T + pow(u1,2) + 3*pow(u2,2) + pow(u3,2)))/(2.*cv);
                    mat(n1-1, 1) = -u1*R/cv;
                    mat(n1-1, 2) = u2-2.*R*u2/cv;
                    mat(n1-1, n1-1) = gamma;
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
            mat(2, 0) = -u3;
            mat(2, 3) = 1.;
            
            mat(3, 0) = R*u2/cv;
            mat(3, 2) = -R/cv;
            
            mat(n1-1, 0) = ((-cv + R)*u2*u3)/cv;
            mat(n1-1, 2) = -u3*R/cv;
            mat(n1-1, 3) = u2;
        }
            break;
            
        default:
            libmesh_assert_msg(false, "invalid dim");
            break;
    }
}



void
FluidElemBase::
calculate_advection_flux_jacobian_u3_derivative(const unsigned int calculate_dim,
                                                const PrimitiveSolution& sol,
                                                DenseRealMatrix& mat)
{
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector
    
    const unsigned int n1 = 2 + dim;
    
    mat.zero();
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    T = sol.T,
    gamma = flight_condition->gas_property.gamma,
    R = flight_condition->gas_property.R,
    cv = flight_condition->gas_property.cv;
    
    switch (calculate_dim)
    {
        case 0:
        {
            switch (dim)
            {
                case 3:
                {
                    mat(3, 0) = -u1;

                    mat(3, 1) = 1.;

                    mat(1, 3) = -R/cv;
                    mat(n1-1, 3) = -u1*R/cv;
                }
                    
                case 2:
                {
                }
                    
                case 1:
                {
                    mat(1, 0) = R*u3/cv;
                    mat(n1-1, 0) = ((-cv + R)*u1*u3)/cv;

                    mat(n1-1, 1) = u3;
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
                    mat(2, 3) = -R/cv;
                    
                    mat(3, 0) = -u2;
                    mat(3, 2) = 1.;
                    
                    mat(n1-1, 3) = -u2*R/cv;
                }
                    
                case 2:
                {
                    mat(2, 0) = R*u3/cv;
                    
                    mat(n1-1, 0) = ((-cv + R)*u2*u3)/cv;
                    mat(n1-1, 2) = u3;
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
            mat(1, 0) = -u1;
            mat(1, 1) = 1.;
            
            mat(2, 0) = -u2;
            mat(2, 2) = 1.;
            
            mat(3, 0) = (-2+R/cv)*u3;
            mat(3, 3) = (2.0-R/cv);
            
            mat(n1-1, 0) = (-2*pow(cv,2)*T + R*(pow(u1,2) + pow(u2,2) + 3*pow(u3,2)) -
                            cv*(2*R*T + pow(u1,2) + pow(u2,2) + 3*pow(u3,2)))/(2.*cv);
            mat(n1-1, 1) = -u1*R/cv;
            mat(n1-1, 2) = -u2*R/cv;
            mat(n1-1, 3) = u3-2.*R*u3/cv;
            mat(n1-1, n1-1) = gamma;
        }
            break;
            
        default:
            libmesh_assert_msg(false, "invalid dim");
            break;
    }
}



void
FluidElemBase::
calculate_advection_flux_jacobian_T_derivative(const unsigned int calculate_dim,
                                               const PrimitiveSolution& sol,
                                               DenseRealMatrix& mat)
{
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector
    
    const unsigned int n1 = 2 + dim;
    
    mat.zero();
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    R = flight_condition->gas_property.R,
    cv = flight_condition->gas_property.cv;
    
    switch (calculate_dim)
    {
        case 0:
        {
            switch (dim)
            {
                case 3:
                {
                }
                    
                case 2:
                {
                }
                    
                case 1:
                {
                    mat(n1-1, 0) = -(cv+R)*u1;
                    mat(n1-1, 1) = cv+R;
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
                }
                    
                case 2:
                {
                    mat(n1-1, 0) = -(cv+R)*u2;
                    mat(n1-1, 2) = cv+R;
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
            mat(n1-1, 0) = -(cv+R)*u3;
            mat(n1-1, 3) = cv+R;
        }
            break;
            
        default:
            libmesh_assert_msg(false, "invalid dim");
            break;
    }
}




void
FluidElemBase::
calculate_advection_left_eigenvector_and_inverse_rho_derivative_for_normal(const PrimitiveSolution& sol,
                                                                           const libMesh::Point& normal,
                                                                           DenseRealMatrix& eig_vals,
                                                                           DenseRealMatrix& l_eig_mat,
                                                                           DenseRealMatrix& l_eig_mat_inv_tr)
{
    eig_vals.zero(); l_eig_mat.zero(); l_eig_mat_inv_tr.zero();
}





void
FluidElemBase::
calculate_advection_left_eigenvector_and_inverse_u1_derivative_for_normal(const PrimitiveSolution& sol,
                                                                          const libMesh::Point& normal,
                                                                          DenseRealMatrix& eig_vals,
                                                                          DenseRealMatrix& l_eig_mat,
                                                                          DenseRealMatrix& l_eig_mat_inv_tr)
{
    const unsigned int n1 = 2 + dim;
    
    eig_vals.zero(); l_eig_mat.zero(); l_eig_mat_inv_tr.zero();
    
    Real nx=0., ny=0., nz=0., u=0.;
    unsigned int dim_for_eig_vec=100; // initializing with arbitrarily high value
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    T = sol.T,
    a = sol.a,
    gamma = flight_condition->gas_property.gamma,
    cv = flight_condition->gas_property.cv,
    factor1 = cv*T*gamma*(gamma-1.),
    factor2 = cv*T*gamma,
    factor1_sqrt = sqrt(factor1);
    
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
            eig_vals(2, 2) = nx;
        case 2:
            eig_vals(1, 1) = nx;
        case 1:
        {
            eig_vals(0, 0) = nx;
            eig_vals(n1-2, n1-2) = nx;
            eig_vals(n1-1, n1-1) = nx;
        }
    }
    
    // set last two columns of the eigenvector matrices, and all columns of the eigenvector inverse.
    // Note that the dim_for_eig_vec column of the inverse matrix will be overwritten by the appropriate matrix
    switch (dim)
    {
        case 3:
        {
            // for u
            l_eig_mat_inv_tr(3, 0) = u3/factor2;
            
            // for u
            l_eig_mat_inv_tr(1, 2) = u3/factor2;
            l_eig_mat_inv_tr(n1-1, 2) = (u1*u3)/factor2-nx*nz;
            
            // overwrite for u
            l_eig_mat_inv_tr(3, dim_for_eig_vec-1) = 0;
        }
            
        case 2:
        {
            // for u
            l_eig_mat_inv_tr(2, 0) = u2/factor2;
            
            // for u
            l_eig_mat_inv_tr(1, 1) = u2/factor2;
            l_eig_mat_inv_tr(n1-1, 1) = (u1*u2)/factor2-nx*ny;
            
            // overwrite for u
            l_eig_mat_inv_tr(2, dim_for_eig_vec-1) = 0;
        }
            
        case 1:
        {
            // for u-a
            l_eig_mat(0, n1-2) = (nx*factor1_sqrt + (-1 + gamma)*u1)/(-1 + gamma);
            l_eig_mat(1, n1-2) = -1.;
            
            // for u+a
            l_eig_mat(0, n1-1) = (-(nx*factor1_sqrt) + (-1 + gamma)*u1)/(-1 + gamma);
            l_eig_mat(1, n1-1) = -1.;
            
            // for u
            l_eig_mat_inv_tr(0, 0) = 1./factor2;
            l_eig_mat_inv_tr(1, 0) = 2.*u1/factor2;
            l_eig_mat_inv_tr(n1-1, 0) = 1.-nx*nx+(k+u1*u1)/factor2;
            
            // overwrite for u
            l_eig_mat_inv_tr(0, dim_for_eig_vec-1) = 0;
            l_eig_mat_inv_tr(1, dim_for_eig_vec-1) = 1;
            l_eig_mat_inv_tr(n1-1, dim_for_eig_vec-1) = u1;
            l_eig_mat_inv_tr.scale_column(dim_for_eig_vec-1, -1.0/factor2);
            
            // for u-a
            l_eig_mat_inv_tr(1, n1-2) = 1./(2.*factor2);
            l_eig_mat_inv_tr(n1-1, n1-2) = (-1.+gamma)*(-factor1*nx + factor1_sqrt*u1)/(2.*factor1*factor1_sqrt);
            
            // for u+a
            l_eig_mat_inv_tr(1, n1-1) = 1./(2.*factor2);
            l_eig_mat_inv_tr(n1-1, n1-1) = (-1.+gamma)*(factor1*nx + factor1_sqrt*u1)/(2.*factor1*factor1_sqrt);
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
                    l_eig_mat(0, 2) = nz/nx;
                }
                    
                case 2:
                {
                    // for u
                    l_eig_mat(0, 1) = ny/nx;
                }
                    
                case 1:
                {
                    // for u
                    l_eig_mat(0, 0) = u/nx;
                    l_eig_mat(1, 0) = -1.;
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
                }
                    
                case 2:
                case 1:
                {
                    // for u
                    l_eig_mat(0, 0) = -1.;
                    
                    // for u
                    l_eig_mat(0, 1) = -u1+nx*u2/ny;
                    l_eig_mat(2, 1) = -nx/ny;
                }
            }
        }
            break;
            
        case 3:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            
            // for u
            l_eig_mat(0, 0) = -1;
            
            // for u
            l_eig_mat(0, 2) = -u1+nx*u3/nz;
            l_eig_mat(3, 2) = -nx/nz;
        }
            break;
    }
}





void
FluidElemBase::
calculate_advection_left_eigenvector_and_inverse_u2_derivative_for_normal(const PrimitiveSolution& sol,
                                                                          const libMesh::Point& normal,
                                                                          DenseRealMatrix& eig_vals,
                                                                          DenseRealMatrix& l_eig_mat,
                                                                          DenseRealMatrix& l_eig_mat_inv_tr)
{
    const unsigned int n1 = 2 + dim;
    
    eig_vals.zero(); l_eig_mat.zero(); l_eig_mat_inv_tr.zero();
    
    Real nx=0., ny=0., nz=0., u=0.;
    unsigned int dim_for_eig_vec=100; // initializing with arbitrarily high value
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    T = sol.T,
    a = sol.a,
    gamma = flight_condition->gas_property.gamma,
    cv = flight_condition->gas_property.cv,
    factor1 = cv*T*gamma*(gamma-1.),
    factor2 = cv*T*gamma,
    factor1_sqrt = sqrt(factor1);
    
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
            eig_vals(2, 2) = ny;
        case 2:
            eig_vals(1, 1) = ny;
        case 1:
        {
            eig_vals(0, 0) = ny;
            eig_vals(n1-2, n1-2) = ny;
            eig_vals(n1-1, n1-1) = ny;
        }
    }
    
    // set last two columns of the eigenvector matrices, and all columns of the eigenvector inverse.
    // Note that the dim_for_eig_vec column of the inverse matrix will be overwritten by the appropriate matrix
    switch (dim)
    {
        case 3:
        {
            // for u
            l_eig_mat_inv_tr(3, 1) = u3/factor2;
            
            // for u
            l_eig_mat_inv_tr(2, 2) = u3/factor2;
            l_eig_mat_inv_tr(n1-1, 2) = (u2*u3)/factor2-ny*nz;
            
            // overwrite for dim_for_eig_vec-1
            l_eig_mat_inv_tr(3, dim_for_eig_vec-1) = 0.;
        }
            
        case 2:
        {
            // for u-a
            l_eig_mat(2, n1-2) = -1.;
            
            // for u+a
            l_eig_mat(2, n1-1) = -1.;
            
            // for u
            l_eig_mat_inv_tr(2, 0) = u1/factor2;

            // for u
            l_eig_mat_inv_tr(0, 1) = 1./factor2;
            l_eig_mat_inv_tr(1, 1) = u1/factor2;
            l_eig_mat_inv_tr(2, 1) = 2.*u2/factor2;
            l_eig_mat_inv_tr(n1-1, 1) = 1.-ny*ny+(k+u2*u2)/factor2;
            
            // overwrite for u
            l_eig_mat_inv_tr(0, dim_for_eig_vec-1) = 0.;
            l_eig_mat_inv_tr(1, dim_for_eig_vec-1) = 0.;
            l_eig_mat_inv_tr(2, dim_for_eig_vec-1) = 1.;
            
            // for u-a
            l_eig_mat_inv_tr(2, n1-2) = 1./(2.*factor2);
            
            // for u+a
            l_eig_mat_inv_tr(2, n1-1) = 1./(2.*factor2);
        }
            
        case 1:
        {
            // for u-a
            l_eig_mat(0, n1-2) = (ny*factor1_sqrt + (-1 + gamma)*u2)/(-1 + gamma);
            
            // for u+a
            l_eig_mat(0, n1-1) = (-ny*factor1_sqrt + (-1 + gamma)*u2)/(-1 + gamma);
            
            // for u
            l_eig_mat_inv_tr(n1-1, 0) = (u1*u2)/factor2-nx*ny;
            
            // overwrite for u
            l_eig_mat_inv_tr(n1-1, dim_for_eig_vec-1) = u2;
            l_eig_mat_inv_tr.scale_column(dim_for_eig_vec-1, -1.0/factor2);
            
            // for u-a
            l_eig_mat_inv_tr(n1-1, n1-2) = (-1.+gamma)*(-factor1*ny + factor1_sqrt*u2)/(2.*factor1*factor1_sqrt);
            
            // for u+a
            l_eig_mat_inv_tr(n1-1, n1-1) = (-1.+gamma)*(factor1*ny + factor1_sqrt*u2)/(2.*factor1*factor1_sqrt);
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
                }
                    
                case 2:
                {
                    // for u
                    l_eig_mat(0, 1) = -1.;
                }
                    
                case 1:
                {
                    // for u
                    l_eig_mat(0, 0) = ny*u1/nx-u2;
                    l_eig_mat(1, 0) = -ny/nx;
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
                    l_eig_mat(0, 2) = nz/ny;
                }
                    
                case 2:
                case 1:
                {
                    // for u
                    l_eig_mat(0, 0) = nx/ny;
                    
                    // for u
                    l_eig_mat(0, 1) = u/ny;
                    l_eig_mat(2, 1) = -1.;
                }
            }
        }
            break;
            
        case 3:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            
            // for u
            l_eig_mat(0, 1) = -1;
            
            // for u
            l_eig_mat(0, 2) = -u2+ny*u3/nz;
            l_eig_mat(3, 2) = -ny/nz;
        }
            break;
    }
}





void
FluidElemBase::
calculate_advection_left_eigenvector_and_inverse_u3_derivative_for_normal(const PrimitiveSolution& sol,
                                                                          const libMesh::Point& normal,
                                                                          DenseRealMatrix& eig_vals,
                                                                          DenseRealMatrix& l_eig_mat,
                                                                          DenseRealMatrix& l_eig_mat_inv_tr)
{
    const unsigned int n1 = 2 + dim;
    
    eig_vals.zero(); l_eig_mat.zero(); l_eig_mat_inv_tr.zero();
    
    Real nx=0., ny=0., nz=0., u=0.;
    unsigned int dim_for_eig_vec=100; // initializing with arbitrarily high value
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    T = sol.T,
    a = sol.a,
    gamma = flight_condition->gas_property.gamma,
    cv = flight_condition->gas_property.cv,
    factor1 = cv*T*gamma*(gamma-1.),
    factor2 = cv*T*gamma,
    factor1_sqrt = sqrt(factor1);
    
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
            eig_vals(2, 2) = nz;
        case 2:
            eig_vals(1, 1) = nz;
        case 1:
        {
            eig_vals(0, 0) = nz;
            eig_vals(n1-2, n1-2) = nz;
            eig_vals(n1-1, n1-1) = nz;
        }
    }
    
    // set last two columns of the eigenvector matrices, and all columns of the eigenvector inverse.
    // Note that the dim_for_eig_vec column of the inverse matrix will be overwritten by the appropriate matrix
    switch (dim)
    {
        case 3:
        {
            // for u-a
            l_eig_mat(3, n1-2) = -1.;
            
            // for u+a
            l_eig_mat(3, n1-1) = -1.;
            
            // for u
            l_eig_mat_inv_tr(3, 0) = u1/factor2;
            
            // for u
            l_eig_mat_inv_tr(3, 1) = u2/factor2;
            
            // for u
            l_eig_mat_inv_tr(0, 2) = 1./factor2;
            l_eig_mat_inv_tr(1, 2) = u1/factor2;
            l_eig_mat_inv_tr(2, 2) = u2/factor2;
            l_eig_mat_inv_tr(3, 2) = 2.*u3/factor2;
            l_eig_mat_inv_tr(n1-1, 2) = 1.-nz*nz+(k+u3*u3)/factor2;
            
            // overwrite for u
            l_eig_mat_inv_tr(0, dim_for_eig_vec-1) = 0.;
            l_eig_mat_inv_tr(1, dim_for_eig_vec-1) = 0.;
            l_eig_mat_inv_tr(2, dim_for_eig_vec-1) = 0.;
            l_eig_mat_inv_tr(3, dim_for_eig_vec-1) = 1.;
            
            // for u-a
            l_eig_mat_inv_tr(3, n1-2) = 1.0/(2.0*factor2);
            
            // for u+a
            l_eig_mat_inv_tr(3, n1-1) = 1.0/(2.0*factor2);
        }
            
        case 2:
        {
            // for u
            l_eig_mat_inv_tr(n1-1, 1) = (u2*u3)/factor2-ny*nz;
        }
            
        case 1:
        {
            // for u-a
            l_eig_mat(0, n1-2) = (nz*factor1_sqrt + (-1 + gamma)*u3)/(-1 + gamma);
            
            // for u+a
            l_eig_mat(0, n1-1) = (-(nz*factor1_sqrt) + (-1 + gamma)*u3)/(-1 + gamma);
            
            // for u
            l_eig_mat_inv_tr(n1-1, 0) = (u1*u3)/factor2-nx*nz;
            
            // overwrite for u
            l_eig_mat_inv_tr(n1-1, dim_for_eig_vec-1) = u3;
            l_eig_mat_inv_tr.scale_column(dim_for_eig_vec-1, -1.0/factor2);
            
            // for u-a
            l_eig_mat_inv_tr(n1-1, n1-2) = (-1 + gamma)*(-factor1*nz + factor1_sqrt*u3)/(2.*factor1*factor1_sqrt);
            
            // for u+a
            l_eig_mat_inv_tr(n1-1, n1-1) = (-1 + gamma)*(factor1*nz + factor1_sqrt*u3)/(2.*factor1*factor1_sqrt);
            
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
                    l_eig_mat(0, 2) = -1.;
                }
                    
                case 2:
                {
                }
                    
                case 1:
                {
                    // for u
                    l_eig_mat(0, 0) = nz*u1/nx-u3;
                    l_eig_mat(1, 0) = -nz/nx;
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
                    l_eig_mat(0, 2) = -1;
                }
                    
                case 2:
                case 1:
                {
                    // for u
                    l_eig_mat(0, 1) = nz*u2/ny-u3;
                    l_eig_mat(2, 1) = -nz/ny;
                }
            }
        }
            break;
            
        case 3:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            
            // for u
            l_eig_mat(0, 0) = nx/nz;
            
            // for u
            l_eig_mat(0, 1) = ny/nz;
            
            // for u
            l_eig_mat(0, 2) = u/nz;
            l_eig_mat(3, 2) = -1.;
        }
            break;
    }
}





void
FluidElemBase::
calculate_advection_left_eigenvector_and_inverse_T_derivative_for_normal(const PrimitiveSolution& sol,
                                                                         const libMesh::Point& normal,
                                                                         DenseRealMatrix& eig_vals,
                                                                         DenseRealMatrix& l_eig_mat,
                                                                         DenseRealMatrix& l_eig_mat_inv_tr)
{
    const unsigned int n1 = 2 + dim;
    
    eig_vals.zero(); l_eig_mat.zero(); l_eig_mat_inv_tr.zero();
    
    Real nx=0., ny=0., nz=0., u=0.;
    unsigned int dim_for_eig_vec=100; // initializing with arbitrarily high value
    
    const Real u1 = sol.u1,
    u2 = sol.u2,
    u3 = sol.u3,
    k = sol.k,
    T = sol.T,
    a = sol.a,
    gamma = flight_condition->gas_property.gamma,
    R     = flight_condition->gas_property.R,
    cv = flight_condition->gas_property.cv,
    factor1 = cv*T*gamma*(gamma-1.),
    factor2 = cv*T*T*gamma,
    factor1_sqrt = sqrt(factor1);
    
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
            eig_vals(2, 2) = 0.;
        case 2:
            eig_vals(1, 1) = 0.;
        case 1:
        {
            eig_vals(0, 0) = 0.;
            eig_vals(n1-2, n1-2) = -sqrt(gamma*R/T)/2.;
            eig_vals(n1-1, n1-1) =  sqrt(gamma*R/T)/2.;
        }
    }
    
    // set last two columns of the eigenvector matrices, and all columns of the eigenvector inverse.
    // Note that the dim_for_eig_vec column of the inverse matrix will be overwritten by the appropriate matrix
    switch (dim)
    {
        case 3:
        {
            // for u-a
            l_eig_mat(3, n1-2) = -(cv*gamma*nz)/(2.*factor1_sqrt);
            
            // for u+a
            l_eig_mat(3, n1-1) =  (cv*gamma*nz)/(2.*factor1_sqrt);
            
            // for u
            l_eig_mat_inv_tr(3, 0) = -u1*u3/factor2;
            
            // for u
            l_eig_mat_inv_tr(3, 1) = -u2*u3/factor2;
            
            // for u
            l_eig_mat_inv_tr(0, 2) = -u3/factor2;
            l_eig_mat_inv_tr(1, 2) = -u1*u3/factor2;
            l_eig_mat_inv_tr(2, 2) = -u2*u3/factor2;
            l_eig_mat_inv_tr(3, 2) = -u3*u3/factor2;
            l_eig_mat_inv_tr(n1-1, 2) = -u3*k/factor2;
            
            // overwrite for u
            l_eig_mat_inv_tr(3, dim_for_eig_vec-1) = u3;
            
            // for u-a
            l_eig_mat_inv_tr(3, n1-2) = (factor1*nz - 2*factor1_sqrt*u3)/(4.*factor2*factor1_sqrt);
            
            // for u+a
            l_eig_mat_inv_tr(3, n1-1) = (-factor1*nz - 2*factor1_sqrt*u3)/(4.*factor2*factor1_sqrt);
            
        }
            
        case 2:
        {
            // for u-a
            l_eig_mat(2, n1-2) = -(cv*gamma*ny)/(2.*factor1_sqrt);
            
            // for u+a
            l_eig_mat(2, n1-1) =  (cv*gamma*ny)/(2.*factor1_sqrt);
            
            // for u
            l_eig_mat_inv_tr(2, 0) = -u1*u2/factor2;
            
            // for u
            l_eig_mat_inv_tr(0, 1) = -u2/factor2;
            l_eig_mat_inv_tr(1, 1) = -u1*u2/factor2;
            l_eig_mat_inv_tr(2, 1) = -u2*u2/factor2;
            l_eig_mat_inv_tr(n1-1, 1) = -u2*k/factor2;
            
            // overwrite for u
            l_eig_mat_inv_tr(2, dim_for_eig_vec-1) = u2;
            
            // for u-a
            l_eig_mat_inv_tr(2, n1-2) = (factor1*ny - 2*factor1_sqrt*u2)/(4.*factor2*factor1_sqrt);
            
            // for u+a
            l_eig_mat_inv_tr(2, n1-1) = (-factor1*ny - 2*factor1_sqrt*u2)/(4.*factor2*factor1_sqrt);
        }
            
        case 1:
        {
            // for u-a
            l_eig_mat(0, n1-2) =(cv*gamma*u)/(2.*factor1_sqrt);
            l_eig_mat(1, n1-2) = -(cv*gamma*nx)/(2.*factor1_sqrt);
            
            // for u+a
            l_eig_mat(0, n1-1) = -(cv*gamma*u)/(2.*factor1_sqrt);
            l_eig_mat(1, n1-1) = (cv*gamma*nx)/(2.*factor1_sqrt);
            
            // for u
            l_eig_mat_inv_tr(0, 0) = -u1/factor2;
            l_eig_mat_inv_tr(1, 0) = -u1*u1/factor2;
            l_eig_mat_inv_tr(n1-1, 0) = -u1*k/factor2;
            
            // overwrite for u
            l_eig_mat_inv_tr(0, dim_for_eig_vec-1) = 1.0;
            l_eig_mat_inv_tr(1, dim_for_eig_vec-1) = u1;
            l_eig_mat_inv_tr(n1-1, dim_for_eig_vec-1) = k;
            l_eig_mat_inv_tr.scale_column(dim_for_eig_vec-1, 1.0/factor2);
            
            // for u-a
            l_eig_mat_inv_tr(0, n1-2) = -1.0/(2.0*factor2);
            l_eig_mat_inv_tr(1, n1-2) = (factor1*nx - 2*factor1_sqrt*u1)/(4.*factor2*factor1_sqrt);
            l_eig_mat_inv_tr(n1-1, n1-2) = (factor1*u - factor1_sqrt*2.*k)/(4.*factor2*factor1_sqrt);
            
            // for u+a
            l_eig_mat_inv_tr(0, n1-1) = -1.0/(2.0*factor2);
            l_eig_mat_inv_tr(1, n1-1) = (-factor1*nx - 2*factor1_sqrt*u1)/(4.*factor2*factor1_sqrt);
            l_eig_mat_inv_tr(n1-1, n1-1) = (-factor1*u - factor1_sqrt*2.*k)/(4.*factor2*factor1_sqrt);
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
                }
                    
                case 2:
                {
                }
                    
                case 1:
                {
                    // for u
                    l_eig_mat(0, 0) = -cv*gamma;
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
                }
                    
                case 2:
                {
                    // for u
                    l_eig_mat(0, 1) = -cv*gamma;
                }
            }
        }
            break;
            
        case 3:
        {
            // set values in the left eigenvector matrix and eigenvalue matrix
            
            // for u
            l_eig_mat(0, 2) = -cv*gamma;
        }
            break;
    }
}







