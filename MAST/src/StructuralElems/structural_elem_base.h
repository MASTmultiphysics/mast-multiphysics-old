

#ifndef __mast_structural_element_base_h__
#define __mast_structural_element_base_h__

// C++ includes
#include <string>
#include <vector>


// MAST includes
#include "StructuralElems/structural_system_base.h"

// libMesh includes
#include "libmesh/quadrature.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"


class StructuralElementBase
{
public:
    /*!
     *   Constructor
     */
    StructuralElementBase();
    
    virtual ~StructuralElementBase();
    
//    virtual void getStressTensor(const FESystem::Geometry::Point& pt, const FESystem::Numerics::VectorBase<FESystemDouble>& sol,
//                                 FESystem::Numerics::MatrixBase<FESystemDouble>& mat) = 0;
//    
//    virtual void transformMatrixToGlobalSystem(const FESystem::Numerics::MatrixBase<FESystemDouble>& elem_mat, FESystem::Numerics::MatrixBase<FESystemDouble>& global_mat);
//    
//    virtual void transformVectorToGlobalSystem(const FESystem::Numerics::VectorBase<FESystemDouble>& elem_vec, FESystem::Numerics::VectorBase<FESystemDouble>& global_vec);
    
    virtual void initialize(const Elem& elem, const StructuralSystemBase& system);

    
    virtual bool element_time_derivative (bool request_jacobian,
                                          DenseVector<Real>& res,
                                          DenseMatrix<Real>& jac) = 0;
    
    virtual bool mass_residual (bool request_jacobian,
                                DenseVector<Real>& res,
                                DenseMatrix<Real>& jac) = 0;
    
    
    
protected:
    
    virtual void clear();
    
    //void calculateDeformationTransformationMatrix(FESystem::Numerics::MatrixBase<FESystemDouble>& mat);
    
    bool _if_initialized;
    
    const Elem* _elem;
    
    const StructuralSystemBase* _system;
};



#endif // __mast_structural_element_base_h__
