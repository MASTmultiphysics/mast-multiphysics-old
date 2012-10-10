//
//  FluidElementBase.cpp
//  FESystem
//
//  Created by Manav Bhatia on 9/19/12.
//
//


// FESystem includes
#include "Disciplines/Fluid/FluidElementBase.h"
#include "Base/FESystemExceptions.h"
#include "Numerics/DenseMatrix.h"
#include "Numerics/LocalVector.h"
#include "Mesh/ElemBase.h"
#include "FiniteElems/FiniteElementBase.h"
#include "Quadrature/QuadratureBase.h"
#include "Geom/Point.h"



FESystem::Fluid::FluidElementBase::FluidElementBase():
if_initialized(false),
geometric_elem(NULL),
quadrature(NULL),
finite_element(NULL),
if_include_diffusion_flux(false),
interpolated_sol(NULL),
solution(NULL),
velocity(NULL),
dt(0.0),
cp(0.0),
cv(0.0),
gamma(0.0),
R(0.0),
s0(0.0),
p0(0.0),
T0(0.0),
mu(0.0),
rho(0.0),
p(0.0),
T(0.0),
u1(0.0),
u2(0.0),
u3(0.0),
e_tot(0.0),
s(0.0),
v(0.0),
a(0.0),
d(0.0),
h(0.0),
e(0.0),
alpha_p(0.0),
beta_T(0.0),
gamma_bar(0.0),
k(0.0),
e_1(0.0),
e_1_bar(0.0),
dp_drho_T(0.0),
dp_dT_rho(0.0),
de_drho_T(0.0),
de_dT_rho(0.0),
drho_dp_T(0.0),
drho_dT_p(0.0),
de_dp_T(0.0),
de_dT_p(0.0),
e_rho_1(0.0),
e_rho_2(0.0),
e_rho_3(0.0),
e_rho_4(0.0),
e_p_1(0.0),
e_p_2(0.0),
e_p_3(0.0),
e_p_4(0.0),
s1(0.0),
s2(0.0),
e_c_1(0.0),
e_c_2(0.0),
e_c_3(0.0),
e_c_4(0.0)
{
    this->interpolated_sol = new FESystem::Numerics::LocalVector<FESystemDouble>;
}




FESystem::Fluid::FluidElementBase::~FluidElementBase()
{
    this->clear();
    if (this->interpolated_sol != NULL) delete this->interpolated_sol;
}



void
FESystem::Fluid::FluidElementBase::clear()
{
    this->if_initialized = false;
    this->geometric_elem = NULL;
    this->quadrature = NULL;
    this->finite_element = NULL;
    this->if_include_diffusion_flux = false;
    this->solution = NULL;
    this->velocity = NULL;
    this->dt = 0.0;
    this->cp = 0.0;
    this->cv = 0.0;
    this->gamma = 0.0;
    this->R = 0.0;
    this->s0 = 0.0;
    this->p0 = 0.0;
    this->T0 = 0.0;
    this->mu = 0.0;
    this->rho = 0.0;
    this->p = 0.0;
    this->T = 0.0;
    this->u1 = 0.0;
    this->u2 = 0.0;
    this->u3 = 0.0;
    this->e_tot = 0.0;
    this->s = 0.0;
    this->v = 0.0;
    this->a = 0.0;
    this->d = 0.0;
    this->h = 0.0;
    this->e = 0.0;
    this->alpha_p = 0.0;
    this->beta_T = 0.0;
    this->gamma_bar = 0.0;
    this->k = 0.0;
    this->e_1 = 0.0;
    this->e_1_bar = 0.0;
    this->dp_drho_T = 0.0;
    this->dp_dT_rho = 0.0;
    this->de_drho_T = 0.0;
    this->de_dT_rho = 0.0;
    this->drho_dp_T = 0.0;
    this->drho_dT_p = 0.0;
    this->de_dp_T = 0.0;
    this->de_dT_p = 0.0;
    this->e_rho_1 = 0.0;
    this->e_rho_2 = 0.0;
    this->e_rho_3 = 0.0;
    this->e_rho_4 = 0.0;
    this->e_p_1 = 0.0;
    this->e_p_2 = 0.0;
    this->e_p_3 = 0.0;
    this->e_p_4 = 0.0;
    this->s1 = 0.0;
    this->s2 = 0.0;
    this->e_c_1 = 0.0;
    this->e_c_2 = 0.0;
    this->e_c_3 = 0.0;
    this->e_c_4 = 0.0;
}



void
FESystem::Fluid::FluidElementBase::initialize(const FESystem::Mesh::ElemBase& elem, const FESystem::FiniteElement::FiniteElementBase& fe, const FESystem::Quadrature::QuadratureBase& q_rule,
                                              FESystemDouble dt_val, FESystemDouble cp_val, FESystemDouble cv_val)
{
    FESystemAssert0(!this->if_initialized, FESystem::Exception::InvalidState);
    this->dt = dt_val;
    this->cp = cp_val;
    this->cv = cv_val;
    this->R = cp - cv;
    this->gamma = cp/cv;
    
    this->geometric_elem = &elem;
    this->finite_element = &fe;
    this->quadrature = &q_rule;
    
    // arbitrary reference values
    this->s0 = 0.0;
    this->p0 = 1.01335e3; // STP
    this->T0 = 273.0;  // STP
    
    FESystemUInt dim = this->geometric_elem->getDimension(), n1 = 2 + dim;
    this->interpolated_sol->resize(n1);
    
    this->if_initialized = true;
}




void
FESystem::Fluid::FluidElementBase::calculateFluxBoundaryCondition(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, const FESystem::Numerics::VectorBase<FESystemDouble>& mass_flux,
                                                                  const FESystem::Numerics::MatrixBase<FESystemDouble>& momentum_flux_tensor, const FESystem::Numerics::VectorBase<FESystemDouble>& energy_flux,
                                                                  FESystem::Numerics::VectorBase<FESystemDouble>& bc_vec)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    std::pair<FESystemUInt, FESystemUInt> s1 = momentum_flux_tensor.getSize();
    
    FESystemAssert2(mass_flux.getSize() == dim, FESystem::Exception::DimensionsDoNotMatch, mass_flux.getSize(), dim);
    FESystemAssert2(energy_flux.getSize() == dim, FESystem::Exception::DimensionsDoNotMatch, energy_flux.getSize(), dim);
    FESystemAssert4((s1.first == dim) && (s1.second == dim), FESystem::Numerics::MatrixSizeMismatch, s1.first, s1.second, dim, dim);
    
    FESystemAssert2(bc_vec.getSize() == n2, FESystem::Exception::DimensionsDoNotMatch, bc_vec.getSize(), n2);
    
    
    static FESystem::Numerics::LocalVector<FESystemDouble>  Nvec, normal, normal_local, tmp_vec1;
    Nvec.resize(n); normal.resize(3); normal_local.resize(dim); tmp_vec1.resize(dim);
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = q_boundary.getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = q_boundary.getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    bc_vec.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValueForBoundary(b_id, *(q_pts[i]));
        
        this->geometric_elem->calculateBoundaryNormal(b_id, normal);
        normal_local.setSubVectorVals(0, dim-1, 0, dim-1, normal);
        this->finite_element->getShapeFunctionForBoundary(b_id, *(q_pts[i]), Nvec);
        
        // add the mass flux
        bc_vec.addSubVectorVals(0, n-1, 0, n-1, -q_weight[i]*jac*mass_flux.dotProduct(normal_local), Nvec);
        
        // next, add the momentum flux
        momentum_flux_tensor.rightVectorMultiply(normal_local, tmp_vec1);
        for (FESystemUInt ii = 0; ii<dim; ii++)
            bc_vec.addSubVectorVals(n*(ii+1), n*(ii+2)-1, 0, n-1, -q_weight[i]*jac*tmp_vec1.getVal(ii), Nvec);
        
        // next, add the energy flux
        bc_vec.addSubVectorVals(n*(dim+1), n2-1, 0, n-1, -q_weight[i]*jac*energy_flux.dotProduct(normal_local), Nvec);
    }
}



void
FESystem::Fluid::FluidElementBase::calculateResidualVector(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, const FESystem::Numerics::VectorBase<FESystemDouble>& vel,
                                                           FESystem::Numerics::VectorBase<FESystemDouble>& res)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    
    FESystemAssert2(res.getSize() == n2, FESystem::Exception::DimensionsDoNotMatch, res.getSize(), n2);
    
    this->solution = &sol;
    this->velocity = &vel;
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> A, B_mat, B_matdx, LS_mat, diff1, diff2;
    static FESystem::Numerics::LocalVector<FESystemDouble> flux, tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2;
    B_mat.resize(n1, n2); B_matdx.resize(n1, n2); LS_mat.resize(n1, n2); diff1.resize(n1, n1); diff2.resize(n1, n1);
    flux.resize(n1); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1); tmp_vec3_n2.resize(n2); A.resize(n1, n1);
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    res.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        
        // shape functions
        this->calculateOperatorMatrix(*(q_pts[i]), B_mat, false, 0);

        // first update the variables at the current quadrature point
        this->updateVariablesAtQuadraturePoint(B_mat);

        this->calculateDifferentialOperatorMatrix(*(q_pts[i]), LS_mat);
        this->calculateArtificialDiffusionOperator(*(q_pts[i]), diff1, diff2);
        diff1.matrixTransposeRightMultiply(1.0, LS_mat, B_matdx); // LS^T tau
        LS_mat.copyMatrix(B_matdx);

        // contribution from unsteady term
        // interpolate velocity for current point
//        B_mat.rightVectorMultiply(vel, tmp_vec1_n1); // interpolated velocity
//        this->calculateConservationVariableJacobian(A);
//        A.rightVectorMultiply(tmp_vec1_n1, tmp_vec2_n1);
//        
//        // Galerkin contribution of velocity
//        B_mat.leftVectorMultiply(tmp_vec2_n1, tmp_vec3_n2);
//        res.add(q_weight[i]*jac, tmp_vec3_n2); // Bw^T u_dot
//        
//        // LS contribution of velocity
//        LS_mat.leftVectorMultiply(tmp_vec2_n1, tmp_vec3_n2); // LS^T tau A U_dot
//        res.add(q_weight[i]*jac, tmp_vec3_n2);
        
        for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
        {
            this->calculateOperatorMatrix(*(q_pts[i]), B_matdx, true, i_dim); // dBw/dx_i

            // Galerkin contribution from the advection flux terms
            this->calculateAdvectionFlux(i_dim, flux); // F^adv_i
            B_matdx.leftVectorMultiply(flux, tmp_vec3_n2); // dBw/dx_i F^adv_i
            res.add(q_weight[i]*jac, tmp_vec3_n2);
            
            // Least square contribution from flux
            this->calculateAdvectionFluxSpatialDerivative(i_dim, B_matdx, &flux, NULL); // d F^adv_i / dxi
            LS_mat.leftVectorMultiply(flux, tmp_vec3_n2); // LS^T tau F^adv_i
            res.add(-q_weight[i]*jac, tmp_vec3_n2);
        }
    }
}




void
FESystem::Fluid::FluidElementBase::calculateTangentMatrix(const FESystem::Numerics::VectorBase<FESystemDouble>& sol, const FESystem::Numerics::VectorBase<FESystemDouble>& vel,
                                                          FESystem::Numerics::MatrixBase<FESystemDouble>& dres_dx, FESystem::Numerics::MatrixBase<FESystemDouble>& dres_dxdot)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    
    const std::pair<FESystemUInt, FESystemUInt> s_mat1 = dres_dx.getSize(), s_mat2 = dres_dxdot.getSize();
    
    FESystemAssert4((s_mat1.first == n2) && (s_mat1.first == n2), FESystem::Numerics::MatrixSizeMismatch, s_mat1.first, s_mat1.first, n2, n2);
    FESystemAssert4((s_mat2.first == n2) && (s_mat2.first == n2), FESystem::Numerics::MatrixSizeMismatch, s_mat2.first, s_mat2.first, n2, n2);
    
    this->solution = &sol;
    this->velocity = &vel;

    static FESystem::Numerics::DenseMatrix<FESystemDouble> A, B_mat, B_matdx, LS_mat, diff1, diff2, tmp_mat1_n2n2, tmp_mat2_n1n2, tmp_mat3_n1n1;
    static FESystem::Numerics::LocalVector<FESystemDouble> flux, tmp_vec1_n1;
    B_mat.resize(n1, n2); B_matdx.resize(n1, n2); LS_mat.resize(n1, n2); diff1.resize(n1, n1); diff2.resize(n1, n1);
    A.resize(n1, n1); flux.resize(n1); tmp_vec1_n1.resize(n1); tmp_mat1_n2n2.resize(n2,n2); tmp_mat2_n1n2.resize(n1, n2); tmp_mat3_n1n1.resize(n1, n1); A.resize(n1, n1);
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble jac=0.0;
    dres_dx.zero(); dres_dxdot.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        // first update the variables at the current quadrature point
        jac = this->finite_element->getJacobianValue(*(q_pts[i]));
        
        // shape functions
        this->calculateOperatorMatrix(*(q_pts[i]), B_mat, false, 0);
        
        // update all variables at the quadrature point for given shape functions
        this->updateVariablesAtQuadraturePoint(B_mat);

        // get the other matrices of interest
        this->calculateDifferentialOperatorMatrix(*(q_pts[i]), LS_mat);
        this->calculateArtificialDiffusionOperator(*(q_pts[i]), diff1, diff2);
        diff1.matrixTransposeRightMultiply(1.0, LS_mat, tmp_mat2_n1n2); // LS^T tau
        LS_mat.copyMatrix(tmp_mat2_n1n2);
        
        // contribution from unsteady term
        // Galerkin contribution of velocity
        this->calculateConservationVariableJacobian(A);
        A.matrixRightMultiply(1.0, B_mat, tmp_mat2_n1n2); // A Bmat
        B_mat.matrixTransposeRightMultiply(1.0, tmp_mat2_n1n2, tmp_mat1_n2n2);
        dres_dxdot.add(q_weight[i]*jac, tmp_mat1_n2n2);
        
        // LS contribution of velocity
        LS_mat.matrixTransposeRightMultiply(1.0, tmp_mat2_n1n2, tmp_mat1_n2n2); // LS^T tau A Bmat
        dres_dxdot.add(q_weight[i]*jac, tmp_mat1_n2n2);
        
        for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
        {
            this->calculateOperatorMatrix(*(q_pts[i]), B_matdx, true, i_dim); // dBw/dx_i

            // Galerkin contribution from the advection flux terms
            this->calculateAdvectionFluxJacobian(i_dim, A); // TODO: Generalize this to group FEM
            A.matrixRightMultiply(1.0, B_mat, tmp_mat2_n1n2);
            B_matdx.matrixTransposeRightMultiply(1.0, tmp_mat2_n1n2, tmp_mat1_n2n2); // dBw/dx_i^T  dF^adv_i/ dU
            dres_dx.add(q_weight[i]*jac, tmp_mat1_n2n2);
            
            // Least square contribution from flux
            this->calculateAdvectionFluxSpatialDerivative(i_dim, B_matdx, NULL, &tmp_mat2_n1n2); // d^2 F^adv_i / dxi dU
            LS_mat.matrixTransposeRightMultiply(1.0, tmp_mat2_n1n2, tmp_mat1_n2n2); // LS^T tau d^2F^adv_i / dx dU
            dres_dx.add(-q_weight[i]*jac, tmp_mat1_n2n2);
        }
    }
}





void
FESystem::Fluid::FluidElementBase::calculateConservationVariableJacobian(FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    // calculate A0 = d U / d Y , where U = conservation variable vector, and Y = unknown variable vector
    // note that for conservation variables as the unknown, this is an identity matrix

    FESystemUInt dim = this->geometric_elem->getDimension(), n1 = 2 + dim;
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4((s.first == n1) && (s.second == n1), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n1);
    
    mat.setToIdentity(); // for conservative formulation
}




void
FESystem::Fluid::FluidElementBase::calculateAdvectionFluxJacobian(FESystemUInt div_coord, FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector

    FESystemUInt dim = this->geometric_elem->getDimension(), n1 = 2 + dim;
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4((s.first == n1) && (s.second == n1), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n1);
    FESystemAssert0(div_coord < dim, FESystem::Exception::InvalidValue);

    mat.zero();
    FESystemUInt energy_i = n1-1;
    
    switch (div_coord)
    {
        case 0:
        {
            switch (dim)
            {
                case 3:
                {
                    mat.setVal(1, 3, -u3*gamma_bar);
                    
                    mat.setVal(3, 0, -u1*u3);
                    mat.setVal(3, 1, u3);
                    mat.setVal(3, 3, u1);
                    
                    mat.setVal(energy_i, 3, -u1*u3*gamma_bar);
                }

                case 2:
                {
                    mat.setVal(1, 2, -u2*gamma_bar);
                    
                    mat.setVal(2, 0, -u1*u2);
                    mat.setVal(2, 1, u2);
                    mat.setVal(2, 2, u1);
                    
                    mat.setVal(energy_i, 2, -u1*u2*gamma_bar);
                }
                    
                case 1:
                {
                    mat.setVal(0, 1, 1.0); // d U / d (rho u1)
                    
                    mat.setVal(1, 0, a*a-u1*u1-e_1_bar*gamma_bar);
                    mat.setVal(1, 1, u1*(2.0-gamma_bar));
                    mat.setVal(1, energy_i, gamma_bar);
                    
                    mat.setVal(energy_i, 0, u1*e_c_2);
                    mat.setVal(energy_i, 1, e_c_3-u1*u1*gamma_bar);
                    mat.setVal(energy_i, energy_i, u1*e_c_4);
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
                    mat.setVal(2, 3, -u3*gamma_bar);
                    
                    mat.setVal(3, 0, -u2*u3);
                    mat.setVal(3, 2, u3);
                    mat.setVal(3, 3, u2);
                    
                    mat.setVal(energy_i, 3, -u2*u3*gamma_bar);
                }
                    
                case 2:
                {
                    mat.setVal(0, 2, 1.0); // d U / d (rho u2)
                    
                    mat.setVal(1, 0, -u1*u2);
                    mat.setVal(1, 1, u2);
                    mat.setVal(1, 2, u1);

                    mat.setVal(2, 0, a*a-u2*u2-e_1_bar*gamma_bar);
                    mat.setVal(2, 1, -u1*gamma_bar);
                    mat.setVal(2, 2, u2*(2.0-gamma_bar));
                    mat.setVal(2, energy_i, gamma_bar);
                                        
                    mat.setVal(energy_i, 0, u2*e_c_2);
                    mat.setVal(energy_i, 1, -u1*u2*gamma_bar);
                    mat.setVal(energy_i, 2, e_c_3-u2*u2*gamma_bar);
                    mat.setVal(energy_i, energy_i, u2*e_c_4);
                }
                    break;
                    
                case 1:
                    // if second coordinate divergence is being asked for, then the element is atleast 2D
                    FESystemAssert0(false, FESystem::Exception::InvalidValue);
                    break;
            }
        }
            break;

        case 2:
        {
            mat.setVal(0, 3, 1.0); // d U / d (rho u3)
            
            mat.setVal(1, 0, -u1*u3);
            mat.setVal(1, 1, u3);
            mat.setVal(1, 3, u1);
            
            mat.setVal(2, 0, -u2*u3);
            mat.setVal(2, 2, u3);
            mat.setVal(2, 3, u2);
            
            mat.setVal(3, 0, a*a-u3*u3-e_1_bar*gamma_bar);
            mat.setVal(3, 1, -u1*gamma_bar);
            mat.setVal(3, 2, -u2*gamma_bar);
            mat.setVal(3, 3, u3*(2.0-gamma_bar));
            mat.setVal(3, energy_i, gamma_bar);

            mat.setVal(energy_i, 0, u3*e_c_2);
            mat.setVal(energy_i, 1, -u1*u3*gamma_bar);
            mat.setVal(energy_i, 2, -u2*u3*gamma_bar);
            mat.setVal(energy_i, 3, e_c_3-u3*u3*gamma_bar);
            mat.setVal(energy_i, energy_i, u3*e_c_4);
        }
            break;

        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
}




void
FESystem::Fluid::FluidElementBase::calculateAdvectionFlux(const FESystemUInt i, FESystem::Numerics::VectorBase<FESystemDouble>& flux)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemUInt dim = this->geometric_elem->getDimension(), n1 = 2 + dim;
    
    FESystemAssert2(flux.getSize() == n1, FESystem::Exception::DimensionsDoNotMatch, flux.getSize(), n1);
    FESystemAssert0(i<dim, FESystem::Exception::InvalidValue);
    
    FESystemUInt energy_i = n1-1;
    
    // calculate the flux using given flow parameters at this point
    switch (i)
    {
        case 0:
        {
            flux.setVal(0, rho * u1);
            flux.setVal(energy_i, u1 * (rho * e_tot + p));
            switch (dim)
            {
                case 3:
                    flux.setVal(3, rho * u1 * u3);
                case 2:
                    flux.setVal(2, rho * u1 * u2);
                case 1:
                    flux.setVal(1, rho * u1 * u1 + p);
            }
        }
            break;
            
        case 1:
        {
            flux.setVal(0, rho * u2);
            flux.setVal(energy_i, u2 * (rho * e_tot + p));
            switch (dim)
            {
                case 3:
                    flux.setVal(3, rho * u2 * u3);
                case 2:
                    flux.setVal(2, rho * u2 * u2 + p);
                case 1:
                    flux.setVal(1, rho * u2 * u1);
            }
        }
            break;
            
        case 2:
        {
            flux.setVal(0, rho * u3);
            flux.setVal(energy_i, u3 * (rho * e_tot + p));
            switch (dim)
            {
                case 3:
                    flux.setVal(3, rho * u3 * u3 + p);
                case 2:
                    flux.setVal(2, rho * u3 * u2);
                case 1:
                    flux.setVal(1, rho * u3 * u1);
            }
        }
            break;
            
        default:
            FESystemAssert0(false, FESystem::Exception::InvalidValue);
            break;
    }
}





void
FESystem::Fluid::FluidElementBase::calculateAdvectionFluxSpatialDerivative(const FESystemUInt i, const FESystem::Numerics::MatrixBase<FESystemDouble>& dBmat_dx,
                                                                           FESystem::Numerics::VectorBase<FESystemDouble>* flux, FESystem::Numerics::MatrixBase<FESystemDouble>* dflux_dU)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemUInt dim = this->geometric_elem->getDimension(), n1 = 2 + dim, n=this->geometric_elem->getNNodes(), n2 = n*n1;
    const std::pair<FESystemUInt, FESystemUInt> s = dBmat_dx.getSize();
    
    FESystemAssert0(( flux!=NULL ) || (dflux_dU != NULL), FESystem::Exception::InvalidFunctionCall); // both can't be null
    FESystemAssert4((s.first == n1) && (s.second == n2), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n2);
    FESystemAssert0(i < dim, FESystem::Exception::InvalidValue);
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> A;
    static FESystem::Numerics::LocalVector<FESystemDouble> tmp_vec;
    A.resize(n1, n1); tmp_vec.resize(n1);
    
    // calculate the Jacobian
    this->calculateAdvectionFluxJacobian(i, A);
    
    if (flux != NULL)
    {
        FESystemAssert2(flux->getSize() == n1, FESystem::Exception::DimensionsDoNotMatch, flux->getSize(), n1);
        dBmat_dx.rightVectorMultiply(*(this->solution), tmp_vec);
        A.rightVectorMultiply(tmp_vec, *flux);
    }

    if (dflux_dU != NULL)
    {
        FESystemAssert4((s.first == n1) && (s.second == n2), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n2);
        A.matrixRightMultiply(1.0, dBmat_dx, *dflux_dU); // conventional interpolation
    }
}






void
FESystem::Fluid::FluidElementBase::calculateDiffusiveFluxJacobian(FESystemUInt div_coord, FESystemUInt flux_coord, FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    // calculate Aij = d F_diff_j / d x_i, where F_diff_j = Kij d U / dxj is the Euler advection flux vector

    FESystemUInt dim = this->geometric_elem->getDimension(), n1 = 2 + dim;
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4((s.first == n1) && (s.second == n1), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n1);
    FESystemAssert0(div_coord < dim, FESystem::Exception::InvalidValue);
    FESystemAssert0(flux_coord < dim, FESystem::Exception::InvalidValue);

    
}





void
FESystem::Fluid::FluidElementBase::calculateArtificialDiffusionOperator(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& streamline_operator, FESystem::Numerics::MatrixBase<FESystemDouble>& discontinuity_operator)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n = this->geometric_elem->getNNodes(), n1 = 2 + dim;
    const std::pair<FESystemUInt, FESystemUInt> s_mat1 = streamline_operator.getSize(), s_mat2 = discontinuity_operator.getSize();
    
    FESystemAssert4((s_mat1.first == n1) && (s_mat1.second == n1), FESystem::Numerics::MatrixSizeMismatch, s_mat1.first, s_mat1.second, n1, n1);
    FESystemAssert4((s_mat2.first == n1) && (s_mat2.second == n1), FESystem::Numerics::MatrixSizeMismatch, s_mat2.first, s_mat2.second, n1, n1);
    
    static FESystem::Numerics::LocalVector<FESystemDouble> N, dNdx, dNdy, dNdz, u, dN;
    static std::vector<FESystemUInt> deriv;
    N.resize(n); dNdx.resize(n); dNdy.resize(n); dNdz.resize(n); u.resize(dim); dN.resize(dim);
    deriv.resize(dim);
    
    streamline_operator.zero();
    discontinuity_operator.zero();
    
    // calculate the gradients
    switch (dim)
    {
        case 3:
        {
            std::fill(deriv.begin(), deriv.end(), 0);
            deriv[2] = 1;
            this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(deriv, pt, dNdz);
            u.setVal(2, this->u3);
        }

        case 2:
        {
            std::fill(deriv.begin(), deriv.end(), 0);
            deriv[1] = 1;
            this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(deriv, pt, dNdy);
            u.setVal(1, this->u2);
        }

        case 1:
        {
            std::fill(deriv.begin(), deriv.end(), 0);
            deriv[0] = 1;
            this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(deriv, pt, dNdx);
            u.setVal(0, this->u1);
        }
            break;

        default:
            break;
    }

    // calculate the dot product of velocity times gradient of shape function
    FESystemDouble h = 0, u_val = u.getL2Norm(), tau_rho, tau_m, tau_e;
    u.scaleToUnitLength();

    for (FESystemUInt i_nodes=0; i_nodes<n; i_nodes++)
    {
        // set value of shape function gradient
        switch (dim)
        {
            case 3:
                dN.setVal(2, dNdz.getVal(i_nodes));
            case 2:
                dN.setVal(1, dNdy.getVal(i_nodes));
            case 1:
                dN.setVal(0, dNdx.getVal(i_nodes));
                break;
                
            default:
                break;
        }
        
        h += fabs(dN.dotProduct(u));
    }
    
    h = 2.0/h;
    
    // now set the value of streamwise dissipation
    tau_rho = 1.0/sqrt(pow(2.0/this->dt, 2)+ pow(2.0/h*(u_val+this->a), 2));
    tau_m = 1.0/sqrt(pow(2.0/this->dt, 2)+ pow(2.0/h*(u_val+this->a), 2));
    tau_e = 1.0/sqrt(pow(2.0/this->dt, 2)+ pow(2.0/h*(u_val+this->a), 2));
    
//    streamline_operator.setVal(0, 0, tau_rho);
//    switch (dim)
//    {
//        case 3:
//            streamline_operator.setVal(3, 3, tau_m);
//        case 2:
//            streamline_operator.setVal(2, 2, tau_m);
//        default:
//            streamline_operator.setVal(1, 1, tau_m);
//            break;
//    }
//    streamline_operator.setVal(n1-1, n1-1, tau_e);
    
}





void
FESystem::Fluid::FluidElementBase::calculateDifferentialOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4((s.first == n1) && (s.second == n2), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n2);
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> A, B_mat, tmp_mat;
    A.resize(n1, n1); B_mat.resize(n1, n2); tmp_mat.resize(n1, n2);
    
    mat.zero();
        
    // contribution of unsteady term
    this->calculateConservationVariableJacobian(A);
    this->calculateOperatorMatrix(pt, B_mat, false, 0);
    A.matrixRightMultiply(1.0, B_mat, mat);

    // contribution of advection flux term
    for (FESystemUInt i=0; i<dim; i++)
    {
        this->calculateAdvectionFluxJacobian(i, A);
        this->calculateOperatorMatrix(pt, B_mat, true, i);
        A.matrixRightMultiply(1.0, B_mat, tmp_mat);
        
        mat.add(1.0, tmp_mat);
    }
}





void
FESystem::Fluid::FluidElementBase::calculateOperatorMatrix(const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat, FESystemBoolean if_strain, FESystemUInt deriv_dim)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    FESystemAssert4((s.first == n1) && (s.second == n2), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n2);
    
    static std::vector<FESystemUInt> derivatives;

    static FESystem::Numerics::LocalVector<FESystemDouble> Nvec;
    Nvec.resize(n); Nvec.zero();
    B_mat.zero();
    
    // prepare the shape function derivative vector
    if (if_strain)
    {
        FESystemAssert0(deriv_dim < dim, FESystem::Exception::InvalidValue);
        
        derivatives.resize(dim);
        std::fill(derivatives.begin(), derivatives.end(), 0);
        derivatives[deriv_dim] = 1;
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinates(derivatives, pt, Nvec);
    }
    else
        this->finite_element->getShapeFunction(pt, Nvec);
        
    // now put these values in the operator matrix
    for (FESystemUInt i=0; i<n1; i++)
        B_mat.setRowVals(i, i*n, (i+1)*n-1, Nvec);
}




void
FESystem::Fluid::FluidElementBase::updateVariablesAtQuadraturePoint(const FESystem::Numerics::MatrixBase<FESystemDouble>& Bmat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;

    const std::pair<FESystemUInt, FESystemUInt> s = Bmat.getSize();
    
    FESystemAssert4((s.first == n1) && (s.second == n2), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n2);
    
    // set the value of interpolated sol through interpolation
    Bmat.rightVectorMultiply(*(this->solution), *(this->interpolated_sol));
    
    this->rho = this->interpolated_sol->getVal(0);

    switch (dim)
    {
        case 3:
            this->u3 = this->interpolated_sol->getVal(3)/rho;
        case 2:
            this->u2 = this->interpolated_sol->getVal(2)/rho;
        case 1:
            this->u1 = this->interpolated_sol->getVal(1)/rho;
    }

    this->e_tot = this->interpolated_sol->getVal(n1-1)/rho; // total energy: defined as (e + k)
    this->k = 0.5 * (u1*u1 + u2*u2 + u3*u3); // kinetic energy
    this->e =  e_tot - k;// internal energy: defined as (cv * T)
    this->T =  e / cv; // Temperature (K)
    this->p =  R * T * rho; // pressure (using ideal gas relationship)
    this->mu =  cp * T * (1.0 - log(T / T0)) + R * T * log(p / p0) - T * s0;// equation of state
    this->s =  cp * log(T / T0) - R * log(p / p0) + s0;// entropy
    this->v =  1.0 / rho;// specific volume
    this->h =  cp * T;// enthalpy
    this->alpha_p = 1.0 / T;// expansivity (ideal gas)
    this->beta_T =  1.0 / p; // isothermal compressibility (ideal gas)
    this->a =  sqrt(v * cp / cv / beta_T);// speed of sound (isentropic speed of sound)
    this->d =  v * alpha_p * T/ beta_T;//
    this->gamma_bar = v * alpha_p / beta_T / cv; //
    this->e_1 = h + k;
    this->e_1_bar = h - k;
    this->dp_drho_T = 1.0 / rho / beta_T;
    this->dp_dT_rho = alpha_p / beta_T;
    this->de_drho_T = -(T * alpha_p / beta_T - p) / pow(rho,2);
    this->de_dT_rho = cv;
    this->drho_dp_T = rho * beta_T;
    this->drho_dT_p = - rho * alpha_p;
    this->de_dp_T = (beta_T * p - alpha_p * T) / rho;
    this->de_dT_p = cp - p * alpha_p / rho;
    this->e_rho_1 = rho * de_drho_T + e_tot;
    this->e_rho_2 = e_rho_1 + dp_drho_T;
    this->e_rho_3 = rho * e_tot + p;
    this->e_rho_4 = rho * de_dT_rho + dp_dT_rho;
    this->e_p_1 = drho_dp_T * e_tot + rho * de_dp_T;
    this->e_p_2 = e_p_1 + 1.0;
    this->e_p_3 = rho * e_tot + p;
    this->e_p_4 = drho_dT_p * e_tot + rho * de_dT_p;
    this->s1 = dp_drho_T / rho / T;
    this->s2 = (T / rho * dp_dT_rho - h) / pow(T, 2) +  k / pow(T,2);
    this->e_c_1 = (2.0 * k - e_rho_1) / rho / cv;
    this->e_c_3 = h + k;
    this->e_c_2 = dp_drho_T - e_c_3 + e_c_1 * dp_dT_rho;
    this->e_c_4 = 1.0 + dp_dT_rho / rho / cv;
}


