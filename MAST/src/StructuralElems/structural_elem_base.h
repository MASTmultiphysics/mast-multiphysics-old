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


// MAST includes
#include "StructuralElems/structural_system_base.h"

// libMesh includes
#include "libmesh/elem.h"



namespace MAST {

    // forward declerations
    class StrainOperator;
    class StressTensor;
    
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
        
        
        /*!
         *   clears the data structure for reinitialization
         */
        void clear();
        
        
        /*!
         *   initializes the data structure for the specified system and element
         */
        virtual void initialize(System& sys, Elem* elem);
        
        
        virtual std::auto_ptr<MAST::StressTensor> get_stress_tensor() = 0;
        
        virtual void transform_to_global_system(const DenseMatrix<Real>& local_mat,
                                                DenseMatrix<Real>& global_mat);
        
        virtual void transform_to_local_system(const DenseVector<Real>& global_vec,
                                               DenseVector<Real>& local_vec);
        
        virtual void transform_to_global_system(const DenseVector<Real>& local_vec,
                                                DenseVector<Real>& global_vec);
        
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
        
        virtual bool volume_external_force (bool request_jacobian,
                                            DenseVector<Real>& f,
                                            DenseMatrix<Real>& jac);
        
    protected:
        
        /*!
         *    initialize the finite element in the local coordinate system
         */
        virtual void initialize_local_element() = 0;
        
        
        /*!
         *   matrix that transforms the global dofs to the local element coordinate
         *   system
         */
        void transformation_matrix(DenseMatrix<Real>& mat);
        
        /*!
         *   returns the quadrature and finite element for element integration. 
         *   These are raw pointers created using new. The pointers must be 
         *   deleted at the end of scope.
         */
        virtual void get_fe_and_qrule(FEBase* fe, QBase* qrule);

        /*!
         *   returns the quadrature and finite element for element side 
         *   integration. These are raw pointers created using new. The 
         *   pointers must be deleted at the end of scope.
         */
        virtual void get_side_fe_and_qrule(unsigned int s,
                                           FEBase* fe, QBase* qrule);

        
        /*!
         *   boolean flag to identify if the object is initialized
         */
        bool _if_initialized;
        
        /*!
         *    System to which this system belongs
         */
        System* _system;
        
        /*!
         *   geometric element for which the computations are performed
         */
        const Elem* _elem;
        
        /*!
         *   element property
         */
        const MAST::ElementPropertyCardBase* _property;
        
        /*!
         *   element solution in the local coordinate system
         */
        DenseVector<Real> _local_solution;
        
        /*!
         *   element velocity in the local coordinate system
         */
        DenseVector<Real> _local_velocity;

        /*!
         *   element acceleration in the local coordinate system
         */
        DenseVector<Real> _local_acceleration;
    };
    
    /*!
     *   returns the matrix operator to calculate strain for the specified point
     */
    std::auto_ptr<MAST::StrainOperator> build_strain_operator();
}


#endif // __mast_structural_element_base_h__


