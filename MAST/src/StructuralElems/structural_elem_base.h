//
//  structural_elememt_base.h
//  MAST
//
//  Created by Manav Bhatia on 10/16/13.
//
//

#ifndef __mast_structural_element_base_h__
#define __mast_structural_element_base_h__

// C++ includes
#include <memory>

// libMesh includes
#include "libmesh/system.h"
#include "libmesh/elem.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"

// MAST includes


namespace MAST {

    // forward declerations
    class ElementPropertyCardBase;
    class Temperature;
    
    enum StructuralElemType {
        ELASTIC_ELEMENT_1D,
        ELASTIC_ELEMENT_2D,
        ELASTIC_ELEMENT_3D
    };
    
    
    class StructuralElementBase
    {
    public:
        /*!
         *   Constructor
         */
        StructuralElementBase(System& sys,
                              const Elem& elem,
                              const MAST::ElementPropertyCardBase& p);
        
        virtual ~StructuralElementBase();
        
        
        virtual bool internal_force (bool request_jacobian,
                                     DenseVector<Real>& f,
                                     DenseMatrix<Real>& jac) = 0;
        
        virtual bool damping_force (bool request_jacobian,
                                    DenseVector<Real>& f,
                                    DenseMatrix<Real>& jac);
        
        virtual bool inertial_force (bool request_jacobian,
                                     DenseVector<Real>& f,
                                     DenseMatrix<Real>& jac);
        
        virtual bool side_external_force (bool request_jacobian,
                                          DenseVector<Real>& f,
                                          DenseMatrix<Real>& jac);
        
        virtual bool prestress_force (bool request_jacobian,
                                      DenseVector<Real>& f,
                                      DenseMatrix<Real>& jac) = 0;
        
        virtual bool volume_external_force (bool request_jacobian,
                                            DenseVector<Real>& f,
                                            DenseMatrix<Real>& jac);
        
        /*!
         *   element solution in the local coordinate system
         */
        DenseVector<Real> local_solution;
        
        /*!
         *   element velocity in the local coordinate system
         */
        DenseVector<Real> local_velocity;
        
        /*!
         *   element acceleration in the local coordinate system
         */
        DenseVector<Real> local_acceleration;
        
    protected:
        
        virtual void _transform_to_global_system(const DenseMatrix<Real>& local_mat,
                                                 DenseMatrix<Real>& global_mat) const;
        
        virtual void _transform_to_local_system(const DenseVector<Real>& global_vec,
                                                DenseVector<Real>& local_vec) const;
        
        virtual void _transform_to_global_system(const DenseVector<Real>& local_vec,
                                                 DenseVector<Real>& global_vec) const;

        /*!
         *   matrix that transforms the global dofs to the local element coordinate
         *   system
         */
        virtual const DenseMatrix<Real>& _transformation_matrix() const = 0;
        
        /*!
         *   returns the quadrature and finite element for element integration. 
         *   These are raw pointers created using new. The pointers must be 
         *   deleted at the end of scope.
         */
        virtual void _init_fe_and_qrule(const Elem& e);

        /*!
         *   returns the quadrature and finite element for element side 
         *   integration. These are raw pointers created using new. The 
         *   pointers must be deleted at the end of scope.
         */
        virtual void _get_side_fe_and_qrule(unsigned int s,
                                            std::auto_ptr<FEBase>& fe,
                                            std::auto_ptr<QBase>& qrule);


        /*!
         *    Calculates the force vector and Jacobian due to thermal stresses. 
         *    this should be implemented for each element type
         */
        virtual bool thermal_force(bool request_jacobian,
                                   DenseVector<Real>& f,
                                   DenseMatrix<Real>& jac) = 0;

        /*!
         *    System to which this system belongs
         */
        System& _system;
        
        /*!
         *   geometric element for which the computations are performed
         */
        const Elem& _elem;
        
        /*!
         *   element property
         */
        const MAST::ElementPropertyCardBase& _property;
        
        /*!
         *    function that defines the temperature in the element domain
         */
        MAST::Temperature* _temperature;
        
        /*!
         *   element finite element for computations
         */
        std::auto_ptr<FEBase> _fe;

        /*!
         *   element quadrature rule for computations
         */
        std::auto_ptr<QBase> _qrule;
    };

    /*!
     *    builds the structural element for the specified element type
     */
    std::auto_ptr<MAST::StructuralElementBase>
    build_structural_element(System& sys,
                             const Elem& elem,
                             const MAST::ElementPropertyCardBase& p);
}


#endif // __mast_structural_element_base_h__


