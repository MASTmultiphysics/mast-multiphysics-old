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
interpolated_vel(NULL),
solution(NULL),
velocity(NULL),
dX_dxi(NULL),
dxi_dX(NULL),
B_mat(NULL),
A_entropy(NULL),
A_inv_entropy(NULL),
Ai_Bi_advection(NULL),
N_vec(NULL),
jac(0.0),
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
    this->interpolated_vel = new FESystem::Numerics::LocalVector<FESystemDouble>;
    this->dX_dxi = new FESystem::Numerics::DenseMatrix<FESystemDouble>;
    this->dxi_dX = new FESystem::Numerics::DenseMatrix<FESystemDouble>;
    this->B_mat = new FESystem::Numerics::DenseMatrix<FESystemDouble>;
    this->A_entropy = new FESystem::Numerics::DenseMatrix<FESystemDouble>;
    this->A_inv_entropy = new FESystem::Numerics::DenseMatrix<FESystemDouble>;
    this->Ai_Bi_advection = new FESystem::Numerics::DenseMatrix<FESystemDouble>;
    this->N_vec = new FESystem::Numerics::LocalVector<FESystemDouble>;
    
    this->N_vec_dx.resize(3);
    this->B_mat_dxi.resize(3);
    this->Ai_advection.resize(3);
    for (FESystemUInt i=0; i<3; i++)
    {
        this->N_vec_dx[i] = new FESystem::Numerics::LocalVector<FESystemDouble>;
        this->B_mat_dxi[i] = new FESystem::Numerics::DenseMatrix<FESystemDouble>;
        this->Ai_advection[i] = new FESystem::Numerics::DenseMatrix<FESystemDouble>;
    }
}




FESystem::Fluid::FluidElementBase::~FluidElementBase()
{
    this->clear();
    if (this->interpolated_sol != NULL) delete this->interpolated_sol;
    if (this->interpolated_vel != NULL) delete this->interpolated_vel;
    if (this->dX_dxi != NULL) delete this->dX_dxi;
    if (this->dxi_dX != NULL) delete this->dxi_dX;
    if (this->B_mat != NULL) delete this->B_mat;
    if (this->A_entropy != NULL) delete this->A_entropy;
    if (this->A_inv_entropy != NULL) delete this->A_inv_entropy;
    if (this->Ai_Bi_advection != NULL) delete this->Ai_Bi_advection;
    if (this->N_vec != NULL) delete this->N_vec;
    for (FESystemUInt i=0; i<3; i++)
    {
        if (this->N_vec_dx[i] != NULL) delete this->N_vec_dx[i];
        if (this->B_mat_dxi[i] != NULL) delete this->B_mat_dxi[i];
        if (this->Ai_advection[i] != NULL) delete this->Ai_advection[i];
    }
}



void
FESystem::Fluid::FluidElementBase::clear()
{
    this->if_initialized = false;
    this->interpolated_sol->zero();
    this->interpolated_vel->zero();
    this->dX_dxi->zero();
    this->dxi_dX->zero();
    this->B_mat->zero();
    this->A_entropy->zero();
    this->A_inv_entropy->zero();
    this->Ai_Bi_advection->zero();
    this->N_vec->zero();
    for (FESystemUInt i=0; i<3; i++)
    {
        this->N_vec_dx[i]->zero();
        this->B_mat_dxi[i]->zero();
        this->Ai_advection[i]->zero();
    }
    this->geometric_elem = NULL;
    this->quadrature = NULL;
    this->finite_element = NULL;
    this->if_include_diffusion_flux = false;
    this->solution = NULL;
    this->velocity = NULL;
    this->jac = 0.0;
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
                                              FESystemDouble dt_val, FESystemDouble cp_val, FESystemDouble cv_val, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
                                              const FESystem::Numerics::VectorBase<FESystemDouble>& vel)
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
    
    this->solution = &sol;
    this->velocity = &vel;

    // arbitrary reference values
    this->s0 = 0.0;
    this->p0 = 1.01335e3; // STP
    this->T0 = 273.0;  // STP
    
    FESystemUInt dim = this->geometric_elem->getDimension(), n1 = 2 + dim, n2=n1*this->geometric_elem->getNNodes();
    this->interpolated_sol->resize(n1);
    this->interpolated_vel->resize(n1);

    this->dX_dxi->resize(dim, dim);
    this->dxi_dX->resize(dim, dim);
    this->B_mat->resize(n1, n2);
    this->A_entropy->resize(n1,n1);
    this->A_inv_entropy->resize(n1, n1);
    this->Ai_Bi_advection->resize(n1, n2);
    this->N_vec->resize(n2);
    for (FESystemUInt i=0; i<dim; i++)
    {
        this->N_vec_dx[i]->resize(n2);
        this->B_mat_dxi[i]->resize(n1, n2);
        this->Ai_advection[i]->resize(n1, n1);
    }

    
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
    
    bc_vec.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        this->jac = this->finite_element->getJacobianValueForBoundary(b_id, *(q_pts[i]));
        
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
FESystem::Fluid::FluidElementBase::calculateSolidWallFluxBoundaryCondition(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, FESystem::Numerics::VectorBase<FESystemDouble>& bc_vec)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    
    FESystemAssert2(bc_vec.getSize() == n2, FESystem::Exception::DimensionsDoNotMatch, bc_vec.getSize(), n2);
    
    static FESystem::Numerics::LocalVector<FESystemDouble>  normal, normal_local, tmp_vec1, flux;
    normal.resize(3); normal_local.resize(dim); tmp_vec1.resize(n2); flux.resize(n1);
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = q_boundary.getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = q_boundary.getQuadraturePointWeights();
    
    bc_vec.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        this->geometric_elem->calculateBoundaryNormal(b_id, normal);
        normal_local.setSubVectorVals(0, dim-1, 0, dim-1, normal);
        
        // first update the variables at the current quadrature point
        this->updateVariablesAtQuadraturePointForBoundary(b_id, *(q_pts[i]));
        
        // now calculate the flux
        for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
        {
            flux.zero();
            flux.setVal(i_dim+1, this->p); // only pressure is set
            this->B_mat->leftVectorMultiply(flux, tmp_vec1);
            bc_vec.add(-q_weight[i]*jac*normal_local.getVal(i_dim), tmp_vec1);
        }
    }
}





void
FESystem::Fluid::FluidElementBase::calculateTangentMatrixForSolidWallFluxBoundaryCondition(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, FESystem::Numerics::MatrixBase<FESystemDouble>& dfdx)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    
    const std::pair<FESystemUInt, FESystemUInt> s_mat1 = dfdx.getSize();
    
    FESystemAssert4((s_mat1.first == n2) && (s_mat1.first == n2), FESystem::Numerics::MatrixSizeMismatch, s_mat1.first, s_mat1.first, n2, n2);
    
    static FESystem::Numerics::LocalVector<FESystemDouble>  normal, normal_local;
    static FESystem::Numerics::DenseMatrix<FESystemDouble>  A, tmp_mat_n1n2, tmp_mat2_n2n2;
    normal.resize(3); normal_local.resize(dim);
    A.resize(n1, n1); tmp_mat_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = q_boundary.getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = q_boundary.getQuadraturePointWeights();
    
    dfdx.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        this->geometric_elem->calculateBoundaryNormal(b_id, normal);
        normal_local.setSubVectorVals(0, dim-1, 0, dim-1, normal);

        this->updateVariablesAtQuadraturePointForBoundary(b_id, *(q_pts[i]));

        // now calculate the flux
        for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
        {
            this->calculatePressureFluxJacobianOnSolidWall(i_dim, A);
            A.matrixRightMultiply(1.0, *(this->B_mat), tmp_mat_n1n2);
            this->B_mat->matrixTransposeRightMultiply(1.0, tmp_mat_n1n2, tmp_mat2_n2n2);
            dfdx.add(-q_weight[i]*jac*normal_local.getVal(i_dim), tmp_mat2_n2n2);
        }
    }
}





void
FESystem::Fluid::FluidElementBase::calculateFluxBoundaryConditionUsingLocalSolution(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, FESystem::Numerics::VectorBase<FESystemDouble>& bc_vec)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    
    FESystemAssert2(bc_vec.getSize() == n2, FESystem::Exception::DimensionsDoNotMatch, bc_vec.getSize(), n2);
    
    static FESystem::Numerics::LocalVector<FESystemDouble>  normal, normal_local, tmp_vec1, flux;
    normal.resize(3); normal_local.resize(dim); tmp_vec1.resize(n2); flux.resize(n1);
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = q_boundary.getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = q_boundary.getQuadraturePointWeights();
    
    bc_vec.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {        
        this->geometric_elem->calculateBoundaryNormal(b_id, normal);
        normal_local.setSubVectorVals(0, dim-1, 0, dim-1, normal);

        this->updateVariablesAtQuadraturePointForBoundary(b_id, *(q_pts[i]));

        // now calculate the flux
        for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
        {
            this->calculateAdvectionFlux(i_dim, flux); // F^adv_i
            this->B_mat->leftVectorMultiply(flux, tmp_vec1);
            bc_vec.add(-q_weight[i]*jac*normal_local.getVal(i_dim), tmp_vec1);
        }
    }
}




void
FESystem::Fluid::FluidElementBase::calculateTangentMatrixForFluxBoundaryConditionUsingLocalSolution(const FESystemUInt b_id, const FESystem::Quadrature::QuadratureBase& q_boundary, FESystem::Numerics::MatrixBase<FESystemDouble>& dfdx)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    
    const std::pair<FESystemUInt, FESystemUInt> s_mat1 = dfdx.getSize();
    
    FESystemAssert4((s_mat1.first == n2) && (s_mat1.first == n2), FESystem::Numerics::MatrixSizeMismatch, s_mat1.first, s_mat1.first, n2, n2);
    
    static FESystem::Numerics::LocalVector<FESystemDouble>  normal, normal_local;
    static FESystem::Numerics::DenseMatrix<FESystemDouble>  tmp_mat_n1n2, tmp_mat2_n2n2;
    normal.resize(3); normal_local.resize(dim);
    tmp_mat_n1n2.resize(n1, n2); tmp_mat2_n2n2.resize(n2, n2);
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = q_boundary.getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = q_boundary.getQuadraturePointWeights();
    
    dfdx.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        this->geometric_elem->calculateBoundaryNormal(b_id, normal);
        normal_local.setSubVectorVals(0, dim-1, 0, dim-1, normal);

        this->updateVariablesAtQuadraturePointForBoundary(b_id, *(q_pts[i]));

        // now calculate the flux
        for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
        {
            this->Ai_advection[i_dim]->matrixRightMultiply(1.0, *(this->B_mat), tmp_mat_n1n2);
            this->B_mat->matrixTransposeRightMultiply(1.0, tmp_mat_n1n2, tmp_mat2_n2n2);
            dfdx.add(-q_weight[i]*jac*normal_local.getVal(i_dim), tmp_mat2_n2n2);
        }
    }
}



void
FESystem::Fluid::FluidElementBase::calculateResidualVector(FESystem::Numerics::VectorBase<FESystemDouble>& res)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    
    FESystemAssert2(res.getSize() == n2, FESystem::Exception::DimensionsDoNotMatch, res.getSize(), n2);
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> LS_mat;
    static FESystem::Numerics::LocalVector<FESystemDouble> flux, tmp_vec1_n1, tmp_vec2_n1, tmp_vec3_n2, diff_sens;
    LS_mat.resize(n1, n2);
    flux.resize(n1); tmp_vec1_n1.resize(n1); tmp_vec2_n1.resize(n1); tmp_vec3_n2.resize(n2); diff_sens.resize(n2);
    
    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble diff_val=0.0;
    res.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        // first update the variables at the current quadrature point
        this->updateVariablesAtQuadraturePoint(*(q_pts[i]));

        this->calculateDifferentialOperatorMatrix(LS_mat, diff_val, diff_sens);

        for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
        {
            // Galerkin contribution from the advection flux terms
            this->calculateAdvectionFlux(i_dim, flux); // F^adv_i
            this->B_mat_dxi[i_dim]->leftVectorMultiply(flux, tmp_vec3_n2); // dBw/dx_i F^adv_i
            res.add(q_weight[i]*jac, tmp_vec3_n2);
            
            // discontinuity capturing operator
            this->B_mat_dxi[i_dim]->rightVectorMultiply(*(this->solution), flux);
            this->B_mat_dxi[i_dim]->leftVectorMultiply(flux, tmp_vec3_n2);
            res.add(-q_weight[i]*jac*diff_val, tmp_vec3_n2);
        }
        
        // Least square contribution from flux
        this->Ai_Bi_advection->rightVectorMultiply(*(this->solution), flux); // d F^adv_i / dxi
        LS_mat.leftVectorMultiply(flux, tmp_vec3_n2); // LS^T tau F^adv_i
        res.add(-q_weight[i]*jac, tmp_vec3_n2);
    }
}




void
FESystem::Fluid::FluidElementBase::calculateTangentMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& dres_dx, FESystem::Numerics::MatrixBase<FESystemDouble>& dres_dxdot)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    
    const std::pair<FESystemUInt, FESystemUInt> s_mat1 = dres_dx.getSize(), s_mat2 = dres_dxdot.getSize();
    
    FESystemAssert4((s_mat1.first == n2) && (s_mat1.first == n2), FESystem::Numerics::MatrixSizeMismatch, s_mat1.first, s_mat1.first, n2, n2);
    FESystemAssert4((s_mat2.first == n2) && (s_mat2.first == n2), FESystem::Numerics::MatrixSizeMismatch, s_mat2.first, s_mat2.first, n2, n2);
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> LS_mat, diff2, tmp_mat1_n2n2, tmp_mat2_n1n2;
    static FESystem::Numerics::LocalVector<FESystemDouble> flux, tmp_vec1_n1, tmp_vec2_n2, diff_sens;
    LS_mat.resize(n1, n2); diff2.resize(n1, n1);
    tmp_vec1_n1.resize(n1); tmp_mat1_n2n2.resize(n2,n2); tmp_mat2_n1n2.resize(n1, n2);
    flux.resize(n1); diff_sens.resize(n2); tmp_vec2_n2.resize(n2);

    const std::vector<FESystem::Geometry::Point*>& q_pts = this->quadrature->getQuadraturePoints();
    const std::vector<FESystemDouble>& q_weight = this->quadrature->getQuadraturePointWeights();
    
    FESystemDouble diff_val=0.0;
    dres_dx.zero(); dres_dxdot.zero();
    
    for (FESystemUInt i=0; i<q_pts.size(); i++)
    {
        // update all variables at the quadrature point for given shape functions
        this->updateVariablesAtQuadraturePoint(*(q_pts[i]));

        // get the other matrices of interest
        this->calculateDifferentialOperatorMatrix(LS_mat, diff_val, diff_sens);
        
        // contribution from unsteady term
        // Galerkin contribution of velocity
        this->B_mat->matrixTransposeRightMultiply(1.0, *(this->B_mat), tmp_mat1_n2n2);
        dres_dxdot.add(q_weight[i]*jac, tmp_mat1_n2n2);
        
        // LS contribution of velocity
        LS_mat.matrixTransposeRightMultiply(1.0, *(this->B_mat), tmp_mat1_n2n2); // LS^T tau A Bmat
        dres_dxdot.add(q_weight[i]*jac, tmp_mat1_n2n2);
        
        for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
        {
            // Galerkin contribution from the advection flux terms
            this->Ai_advection[i_dim]->matrixRightMultiply(1.0, *(this->B_mat), tmp_mat2_n1n2);
            this->B_mat_dxi[i_dim]->matrixTransposeRightMultiply(1.0, tmp_mat2_n1n2, tmp_mat1_n2n2); // dBw/dx_i^T  dF^adv_i/ dU
            dres_dx.add(q_weight[i]*jac, tmp_mat1_n2n2);
            
            // discontinuity capturing term
            this->B_mat_dxi[i_dim]->matrixTransposeRightMultiply(1.0, *(this->B_mat_dxi[i_dim]), tmp_mat1_n2n2);
            dres_dx.add(-q_weight[i]*jac*diff_val, tmp_mat1_n2n2);
            
//            // discontinuity capturing term Jac due to coefficient dependence on solution
//            B_matdx.rightVectorMultiply(*(this->solution), flux);
//            B_matdx.leftVectorMultiply(flux, tmp_vec2_n2);
//            for (FESystemUInt ii=0; ii<n2; ii++)
//                for (FESystemUInt jj=0; jj<n2; jj++)
//                    tmp_mat1_n2n2.setVal(ii, jj, tmp_vec2_n2.getVal(ii)*diff_sens.getVal(jj));
//            dres_dx.add(-q_weight[i]*jac, tmp_mat1_n2n2);
        }

        // Lease square contribution of flux gradient
        LS_mat.matrixTransposeRightMultiply(1.0, *(this->Ai_Bi_advection), tmp_mat1_n2n2); // LS^T tau d^2F^adv_i / dx dU
        dres_dx.add(-q_weight[i]*jac, tmp_mat1_n2n2);
    }
}





void
FESystem::Fluid::FluidElementBase::calculateEntropyVariableJacobian(FESystem::Numerics::MatrixBase<FESystemDouble>& dUdV, FESystem::Numerics::MatrixBase<FESystemDouble>& dVdU)
{
    // calculates dU/dV where V is the Entropy variable vector
    
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    // calculate A0 = d U / d Y , where U = conservation variable vector, and Y = unknown variable vector
    // note that for conservation variables as the unknown, this is an identity matrix
    
    FESystemUInt dim = this->geometric_elem->getDimension(), n1 = 2 + dim;
    const std::pair<FESystemUInt, FESystemUInt> s = dUdV.getSize();
    
    FESystemAssert4((s.first == n1) && (s.second == n1), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n1);

    
    // du/dv
    switch (dim)
    {
        case 3:
        {
            dUdV.setVal(0, 3, u3);

            dUdV.setVal(1, 3, u1*u3);

            dUdV.setVal(2, 3, u2*u3);

            dUdV.setVal(3, 0, dUdV.getVal(0, 3));
            dUdV.setVal(3, 1, dUdV.getVal(1, 3));
            dUdV.setVal(3, 2, dUdV.getVal(2, 3));
            dUdV.setVal(3, 3, u3*u3+v/beta_T);
            dUdV.setVal(3, n1-1, u3*(h+k-v*(alpha_p*T+1)/beta_T));
            
            dUdV.setVal(n1-1, 3, dUdV.getVal(3, n1-1));
        }
            
        case 2:
        {
            dUdV.setVal(0, 2, u2);
            
            dUdV.setVal(1, 2, u1*u2);

            dUdV.setVal(2, 0, dUdV.getVal(0, 2));
            dUdV.setVal(2, 1, dUdV.getVal(1, 2));
            dUdV.setVal(2, 2, u2*u2+v/beta_T);
            dUdV.setVal(2, n1-1, u2*(h+k-v*(alpha_p*T+1)/beta_T));
            
            dUdV.setVal(n1-1, 2, dUdV.getVal(2, n1-1));
        }

        case 1:
        {
            dUdV.setVal(0, 0, 1.0);
            dUdV.setVal(0, 1, u1);
            dUdV.setVal(0, n1-1, h+k-v*alpha_p*T/beta_T);

            dUdV.setVal(1, 0, dUdV.getVal(0, 1));
            dUdV.setVal(1, 1, u1*u1+v/beta_T);
            dUdV.setVal(1, n1-1, u1*(h+k-v*(alpha_p*T+1)/beta_T));

            dUdV.setVal(n1-1, 0, dUdV.getVal(0, n1-1));
            dUdV.setVal(n1-1, 1, dUdV.getVal(1, n1-1));
            dUdV.setVal(n1-1, n1-1, pow(h+k,2)+v/beta_T*(cp*T-2*h*alpha_p*T-2*k*(alpha_p*T-1)));

        }
            break;

        default:
            break;
    }
    
    dUdV.scale(beta_T*T/v/v);

    
    // dv/du
    switch (dim)
    {
        case 3:
        {
            dVdU.setVal(0, 3, u3);
            
            dVdU.setVal(1, 3, u1*u3);
            
            dVdU.setVal(2, 3, u2*u3);
            
            dVdU.setVal(3, 0, dVdU.getVal(0, 3));
            dVdU.setVal(3, 1, dVdU.getVal(1, 3));
            dVdU.setVal(3, 2, dVdU.getVal(2, 3));
            dVdU.setVal(3, 3, u3*u3+v/beta_T);
            dVdU.setVal(3, n1-1, u3*(h+k-v*(alpha_p*T+1)/beta_T));
            
            dVdU.setVal(n1-1, 3, dVdU.getVal(3, n1-1));
        }
            
        case 2:
        {
            dVdU.setVal(0, 2, u2);
            
            dVdU.setVal(1, 2, u1*u2);
            
            dVdU.setVal(2, 0, dVdU.getVal(0, 2));
            dVdU.setVal(2, 1, dVdU.getVal(1, 2));
            dVdU.setVal(2, 2, u2*u2+v/beta_T);
            dVdU.setVal(2, n1-1, u2*(h+k-v*(alpha_p*T+1)/beta_T));
            
            dVdU.setVal(n1-1, 2, dVdU.getVal(2, n1-1));
        }
            
        case 1:
        {
            dVdU.setVal(0, 0, 1.0);
            dVdU.setVal(0, 1, u1);
            dVdU.setVal(0, n1-1, h+k-v*alpha_p*T/beta_T);
            
            dVdU.setVal(1, 0, dVdU.getVal(0, 1));
            dVdU.setVal(1, 1, u1*u1+v/beta_T);
            dVdU.setVal(1, n1-1, u1*(h+k-v*(alpha_p*T+1)/beta_T));
            
            dVdU.setVal(n1-1, 0, dVdU.getVal(0, n1-1));
            dVdU.setVal(n1-1, 1, dVdU.getVal(1, n1-1));
            dVdU.setVal(n1-1, n1-1, pow(h+k,2)+v/beta_T*(cp*T-2*h*alpha_p*T-2*k*(alpha_p*T-1)));
        }
            break;
            
        default:
            break;
    }
    
    dVdU.scale(beta_T*T/v/v);
    dVdU.zero();
}




void
FESystem::Fluid::FluidElementBase::calculatePressureFluxJacobianOnSolidWall(FESystemUInt div_coord, FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    // calculate Ai = d F_adv / d x_i, where F_adv is the Euler advection flux vector
    
    FESystemUInt dim = this->geometric_elem->getDimension(), n1 = 2 + dim;
    const std::pair<FESystemUInt, FESystemUInt> s = mat.getSize();
    
    FESystemAssert4((s.first == n1) && (s.second == n1), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n1);
    FESystemAssert0(div_coord < dim, FESystem::Exception::InvalidValue);
    
    mat.zero();
    FESystemUInt energy_i = n1-1;
    
    switch (dim)
    {
        case 3:
            mat.setVal(div_coord+1, 3, -u3);
            
        case 2:
            mat.setVal(div_coord+1, 2, -u2);
            
        case 1:
        {
            mat.setVal(div_coord+1, 0, k);
            mat.setVal(div_coord+1, 1, -u1);
            mat.setVal(div_coord+1, energy_i, 1.0);
        }
            break;
    }

    mat.scale(this->gamma_bar);
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
FESystem::Fluid::FluidElementBase::calculateAdvectionFluxSpatialDerivative(const FESystemUInt i, FESystem::Numerics::VectorBase<FESystemDouble>* flux, FESystem::Numerics::MatrixBase<FESystemDouble>* dflux_dU)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    
    FESystemUInt dim = this->geometric_elem->getDimension(), n1 = 2 + dim, n=this->geometric_elem->getNNodes(), n2 = n*n1;
    
    FESystemAssert0(( flux!=NULL ) || (dflux_dU != NULL), FESystem::Exception::InvalidFunctionCall); // both can't be null
    FESystemAssert0(i < dim, FESystem::Exception::InvalidValue);
    
    static FESystem::Numerics::LocalVector<FESystemDouble> tmp_vec;
    tmp_vec.resize(n1);
    
    // calculate the Jacobian
    if (flux != NULL)
    {
        FESystemAssert2(flux->getSize() == n1, FESystem::Exception::DimensionsDoNotMatch, flux->getSize(), n1);
        this->B_mat_dxi[i]->rightVectorMultiply(*(this->solution), tmp_vec);
        this->Ai_advection[i]->rightVectorMultiply(tmp_vec, *flux);
    }

    if (dflux_dU != NULL)
    {
        const std::pair<FESystemUInt, FESystemUInt> s1 = dflux_dU->getSize();
        FESystemAssert4((s1.first == n1) && (s1.second == n2), FESystem::Numerics::MatrixSizeMismatch, s1.first, s1.second, n1, n2);
        this->Ai_advection[i]->matrixRightMultiply(1.0, *(this->B_mat_dxi[i]), *dflux_dU); // conventional interpolation
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
FESystem::Fluid::FluidElementBase::calculateArtificialDiffusionOperator(FESystem::Numerics::MatrixBase<FESystemDouble>& streamline_operator)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n = this->geometric_elem->getNNodes(), n1 = 2 + dim;
    const std::pair<FESystemUInt, FESystemUInt> s_mat1 = streamline_operator.getSize();
    
    FESystemAssert4((s_mat1.first == n1) && (s_mat1.second == n1), FESystem::Numerics::MatrixSizeMismatch, s_mat1.first, s_mat1.second, n1, n1);
    
    static FESystem::Numerics::LocalVector<FESystemDouble> u, dN;
    u.resize(dim); dN.resize(dim);
    
    streamline_operator.zero();
    
    // calculate the gradients
    switch (dim)
    {
        case 3:
            u.setVal(2, this->u3);

        case 2:
            u.setVal(1, this->u2);

        case 1:
            u.setVal(0, this->u1);
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
        for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
            dN.setVal(i_dim, this->N_vec_dx[i_dim]->getVal(i_nodes));
        
        h += fabs(dN.dotProduct(u));
    }
    
    h = 2.0/h;
    
    // now set the value of streamwise dissipation
    tau_rho = 1.0/sqrt(pow(2.0/this->dt, 2)+ pow(2.0/h*(u_val+this->a), 2));
    tau_m = 1.0/sqrt(pow(2.0/this->dt, 2)+ pow(2.0/h*(u_val+this->a), 2));
    tau_e = 1.0/sqrt(pow(2.0/this->dt, 2)+ pow(2.0/h*(u_val+this->a), 2));
    
    streamline_operator.setVal(0, 0, tau_rho);
    for (FESystemUInt i_dim=0; i_dim<dim; i_dim++)
        streamline_operator.setVal(1+i_dim, 1+i_dim, tau_m);
    streamline_operator.setVal(n1-1, n1-1, tau_e);
}





void
FESystem::Fluid::FluidElementBase::calculateDifferentialOperatorMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& LS_operator, FESystemDouble& discontinuity_val, FESystem::Numerics::VectorBase<FESystemDouble>& discont_operator_sens)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    const std::pair<FESystemUInt, FESystemUInt> s = LS_operator.getSize();
    
    FESystemAssert4((s.first == n1) && (s.second == n2), FESystem::Numerics::MatrixSizeMismatch, s.first, s.second, n1, n2);
    
    static std::vector<FESystem::Numerics::LocalVector<FESystemDouble> > diff_vec(3);
    static FESystem::Numerics::DenseMatrix<FESystemDouble> tmp_mat, tmp_mat_n1n1, diff_operator;
    static FESystem::Numerics::LocalVector<FESystemDouble> vec1, vec2, vec3, vec4, vec5;
    tmp_mat.resize(n1, n2); tmp_mat_n1n1.resize(n1, n1); diff_operator.resize(n1, n2);
    vec1.resize(n1); vec2.resize(n1); vec3.resize(n1); vec4.resize(n2); vec5.resize(n2);
    for (FESystemUInt i=0; i<dim; i++) diff_vec[i].resize(n1);
    
    // contribution of unsteady term
    LS_operator.copyMatrix(*(this->B_mat));

    diff_operator.copyMatrix(*(this->B_mat));
    diff_operator.add(1.0, *(this->Ai_Bi_advection));

    FESystemDouble val1 = 0.0;

    vec2.zero();

    // contribution of advection flux term
    for (FESystemUInt i=0; i<dim; i++)
    {
        this->Ai_advection[i]->matrixTransposeRightMultiply(1.0, *(this->B_mat_dxi[i]), tmp_mat);
        LS_operator.add(1.0, tmp_mat);  // (B + A_i^T B_i)

        this->B_mat_dxi[i]->rightVectorMultiply(*(this->solution), diff_vec[i]);
        this->Ai_advection[i]->rightVectorMultiply(diff_vec[i], vec1);

        vec2.add(1.0, vec1); // sum A_i dU/dx_i
    }

    // add the velocity and calculate the numerator of the discontinuity capturing term coefficient
    vec2.add(1.0, *(this->interpolated_vel)); // add velocity
    this->A_inv_entropy->rightVectorMultiply(vec2, vec1);
    discontinuity_val = vec1.dotProduct(vec2); // this is the numerator term
    
    // calculate the sensitivity of the numerator
    vec3.copyVector(vec1); // A0inv v
    this->A_inv_entropy->leftVectorMultiply(vec2, vec1); // A0inv^T v
    vec3.add(1.0, vec1);
    diff_operator.leftVectorMultiply(vec3, discont_operator_sens); // v^T (A0inv^T + A0inv) (B + A_i B_i)

    // now evaluate the dissipation factor for the discontinuity capturing term
    // this is the denominator term
    
    val1 = 0.0;
    vec5.zero();
    for (FESystemUInt i=0; i<dim; i++)
    {
        vec1.zero();
        tmp_mat.zero();

        for (FESystemUInt j=0; j<dim; j++)
        {
            vec1.add(this->dxi_dX->getVal(i, j), diff_vec[j]);
            tmp_mat.add(this->dxi_dX->getVal(i, j), *(this->B_mat_dxi[j]));
        }
        
        // calculate the value of denominator
        this->A_inv_entropy->rightVectorMultiply(vec1, vec2);
        val1 += vec1.dotProduct(vec2);
        
        // calculate the sensitivity of the denominator
        vec3.copyVector(vec2); // A0inv v
        this->A_inv_entropy->leftVectorMultiply(vec1, vec2); // A0inv^T v
        vec3.add(1.0, vec2); // v^T (A0inv + A0inv^T)
        tmp_mat.leftVectorMultiply(vec3, vec4);
        vec5.add(1.0, vec4);
    }
    
    //    // now calculate the discontinuity capturing operator
//    if ((fabs(val1) > FESystem::Base::getMachineEpsilon<FESystemDouble>()) &&  (fabs(discontinuity_val) > FESystem::Base::getMachineEpsilon<FESystemDouble>()))
    if ((fabs(val1) > 0.0) &&  (fabs(discontinuity_val) > 0.0))
    {
        discont_operator_sens.scale(1.0/discontinuity_val);
        discont_operator_sens.add(-1.0/val1, vec5);
        discontinuity_val /= val1;
        discontinuity_val = sqrt(discontinuity_val);
        discont_operator_sens.scale(0.5*discontinuity_val);
        discont_operator_sens.zero();
    }
    else
    {
        discontinuity_val = 0.0;
        discont_operator_sens.zero();
    }
    
    // scale the LS matrix with the correct factor
    this->calculateArtificialDiffusionOperator(tmp_mat_n1n1);
    tmp_mat_n1n1.matrixRightMultiply(1.0, LS_operator, tmp_mat);
    LS_operator.copyMatrix(tmp_mat);
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
FESystem::Fluid::FluidElementBase::calculateOperatorMatrixForBoundary(const FESystemUInt b_id, const FESystem::Geometry::Point& pt, FESystem::Numerics::MatrixBase<FESystemDouble>& B_mat, FESystemBoolean if_strain, FESystemUInt deriv_dim)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    const std::pair<FESystemUInt, FESystemUInt> s = B_mat.getSize();
    
    FESystemAssert2(pt.getSize() == dim-1, FESystem::Exception::DimensionsDoNotMatch, dim-1, pt.getSize());
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
        this->finite_element->getShapeFunctionDerivativeForPhysicalCoordinatesForBoundary(b_id, derivatives, pt, Nvec);
    }
    else
        this->finite_element->getShapeFunctionForBoundary(b_id, pt, Nvec);
    
    // now put these values in the operator matrix
    for (FESystemUInt i=0; i<n1; i++)
        B_mat.setRowVals(i, i*n, (i+1)*n-1, Nvec);
}



void
FESystem::Fluid::FluidElementBase::updateVariablesAtQuadraturePoint(const FESystem::Geometry::Point& pt)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> tmp_mat;
    tmp_mat.resize(n1, n2);
    
    // update the shape function and Jacobian matrices
    this->calculateOperatorMatrix(pt, *(this->B_mat), false, 0);
    this->B_mat->getRowVals(0, 0, n, *(this->N_vec));
    for (FESystemUInt i=0; i<dim; i++)
    {
        this->calculateOperatorMatrix(pt, *(this->B_mat_dxi[i]), true, i);
        this->B_mat_dxi[i]->getRowVals(0, 0, n, *(this->N_vec_dx[i]));
    }

    // shape function Jacobian matrix
    this->finite_element->getJacobianMatrix(pt, *(this->dX_dxi));
    this->dX_dxi->getInverse(*(this->dxi_dX));

    // determinant of the Jacobian for numerical integration
    this->jac = this->dX_dxi->getDeterminant();
    
    this->updateVariablesForInterpolationOperator(*(this->B_mat));
}



void
FESystem::Fluid::FluidElementBase::updateVariablesAtQuadraturePointForBoundary(const FESystemUInt b_id, const FESystem::Geometry::Point& pt)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> tmp_mat;
    tmp_mat.resize(n1, n2);
    
    // update the shape function and Jacobian matrices
    this->calculateOperatorMatrixForBoundary(b_id, pt, *(this->B_mat), false, 0);
    this->B_mat->getRowVals(0, 0, n, *(this->N_vec));
    for (FESystemUInt i=0; i<dim; i++)
    {
        this->calculateOperatorMatrixForBoundary(b_id, pt, *(this->B_mat_dxi[i]), true, i);
        this->B_mat_dxi[i]->getRowVals(0, 0, n, *(this->N_vec_dx[i]));
    }

    // shape function Jacobian matrix
    this->dX_dxi->zero();  // not used for boundary integrals
    this->dxi_dX->zero();  // not used for boundary integrals
    this->jac = this->finite_element->getJacobianValueForBoundary(b_id, pt);

    this->updateVariablesForInterpolationOperator(*(this->B_mat));
}

    
void
FESystem::Fluid::FluidElementBase::updateVariablesForInterpolationOperator(const FESystem::Numerics::MatrixBase<FESystemDouble>& Bmat)
{
    FESystemAssert0(this->if_initialized, FESystem::Exception::InvalidState);
    FESystemUInt dim = this->geometric_elem->getDimension(), n=this->geometric_elem->getNNodes(), n1 = 2 + dim, n2 = n1*n;
    
    static FESystem::Numerics::DenseMatrix<FESystemDouble> tmp_mat;
    tmp_mat.resize(n1, n2);
    
    // set the value of interpolated sol through interpolation
    this->B_mat->rightVectorMultiply(*(this->solution), *(this->interpolated_sol));
    this->B_mat->rightVectorMultiply(*(this->velocity), *(this->interpolated_vel));
    
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
    
    // update the Jacobian matrices
    this->Ai_Bi_advection->zero();
    for (FESystemUInt i=0; i<dim; i++)
    {
        this->calculateAdvectionFluxJacobian(i, *(this->Ai_advection[i]));
        this->Ai_advection[i]->matrixRightMultiply(1.0, *(this->B_mat_dxi[i]), tmp_mat);
        this->Ai_Bi_advection->add(1.0, tmp_mat);
    }
    this->calculateEntropyVariableJacobian(*(this->A_entropy), *(this->A_inv_entropy));
}



FESystemDouble
FESystem::Fluid::FluidElementBase::estimateJacobianSpectralRadius(const FESystem::Numerics::MatrixBase<FESystemDouble>& mat)
{
    const std::pair<FESystemUInt, FESystemUInt> size = mat.getSize();
    
    FESystemAssert0(size.first == size.second, FESystem::Numerics::MatrixMustBeSquareForOperation);
    
    static FESystem::Numerics::LocalVector<FESystemDouble> vec1, vec2;
    vec1.resize(size.first); vec2.resize(size.first);
    
    vec1.setAllVals(1.0);
    FESystemDouble val1=100.0, val2=0.0;
    FESystemUInt n_iters = 10;
    
    // use power iteration to calculate the largest eigenvalue
    while ((fabs((val1-val2)/val1) > FESystem::Base::getMachineEpsilon<FESystemDouble>()) && (n_iters > 0))
    {
        val2 = val1;
        mat.rightVectorMultiply(vec1, vec2);
        val1 = vec2.dotProduct(vec1)/vec1.getL2Norm();
        vec1.copyVector(vec2);
        vec1.scaleToUnitLength();
        n_iters--;
    }
    return val1;
}



