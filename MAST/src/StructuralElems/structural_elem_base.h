

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

#ifndef LIBMESH_USE_COMPLEX_NUMBERS


enum StructuralElemType {
    // linear elements
    AXIAL_BAR,
    TORSION_BAR,
    EULER_BERNOULLI_BEAM,
    TIMOSHENKO_BEAM,
    MEMBRANE,
    MINDLIN_PLATE,
    DKT_PLATE,
    SOLID,
    // nonlinear elements
    VON_KARMAN_STRAIN
};


class StructuralElementBase
{
public:
    /*!
     *   Constructor
     */
    StructuralElementBase();
    
    virtual ~StructuralElementBase();
    
    
    void initialize();
    
    
    virtual void getStressTensor(const Point& pt,
                                 const DenseVector<Real>& sol,
                                 DenseMatrix<Real>& mat) = 0;
    
    virtual void transform_matrix_to_global_system(const DenseMatrix<Real>& local_mat,
                                                   DenseMatrix<Real>& global_mat);
    
    virtual void transform_vector_to_local_system(const DenseVector<Real>& global_vec,
                                                   DenseVector<Real>& local_vec);

    virtual void transform_vector_to_global_system(const DenseVector<Real>& local_vec,
                                                   DenseVector<Real>& global_vec);
    
    virtual bool element_time_derivative (bool request_jacobian,
                                          DenseVector<Real>& res,
                                          DenseMatrix<Real>& jac) = 0;
    
    virtual bool mass_residual (bool request_jacobian,
                                DenseVector<Real>& res,
                                DenseMatrix<Real>& jac) = 0;
    
    
    
protected:
    
    bool _if_initialized;
    
    const Elem* _elem;
};

#endif // LIBMESH_USE_COMPLEX_NUMBERS


#endif // __mast_structural_element_base_h__
