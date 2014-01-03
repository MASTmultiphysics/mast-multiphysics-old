//
//  structural_elememt_base.h
//  MAST
//
//  Created by Manav Bhatia on 10/16/13.
//
//

#ifndef __MAST_structural_element_base_h__
#define __MAST_structural_element_base_h__

// C++ includes
#include <memory>
#include <map>

// libMesh includes
#include "libmesh/system.h"
#include "libmesh/elem.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"

// MAST includes


namespace MAST {

    // forward declerations
    class ElementPropertyCardBase;
    class FunctionBase;
    class Temperature;
    class SensitivityParameters;
    class BoundaryCondition;
    
    enum StructuralElemType {
        ELASTIC_ELEMENT_1D,
        ELASTIC_ELEMENT_2D,
        ELASTIC_ELEMENT_3D
    };
    
    
    /*!
     *   stress tensor
     */
    class Stress: public DenseMatrix<Real> {
    public:
        Stress():
        DenseMatrix<Real>()
        {
            this->resize(3,3);
        }
        
        virtual ~Stress() { }

        /*!
         *    calculates and returns the von Mises stress for this stress
         *    tensor
         */
        Real von_mises_stress() const {
            
            const DenseMatrix<Real>& s = *this;
            
            Real val =
            pow(s(0,0) - s(1,1), 2) +
            pow(s(1,1) - s(2,2), 2) +
            pow(s(2,2) - s(0,0), 2) +
            6.* (pow(s(0,1), 2) +
                 pow(s(1,2), 2) +
                 pow(s(2,0), 2));
            
            val = sqrt(.5 * val);
            
            return val;
        }
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
        
        
        /*!
         *   returns a constant reference to the finite element object
         */
        System& system()  {
            return _system;
        }

        /*!
         *   returns a constant reference to the finite element object
         */
        const MAST::ElementPropertyCardBase& elem_property()  {
            return _property;
        }

        
        /*!
         *   returns a constant reference to the element
         */
        const Elem& elem() const {
            return _elem;
        }

        
        /*!
         *   returns a constant reference to the element in local coordinate system
         */
        virtual const Elem& local_elem() const = 0;

        
        /*!
         *   returns a constant reference to the quadrature rule
         */
        QBase& quadrature_rule()  {
            return *_qrule;
        }

        
        
        /*!
         *   returns a constant reference to the finite element object
         */
        FEBase& fe()  {
            return *_fe;
        }

        
        /*!
         *   internal force contribution to system residual
         */
        virtual bool internal_force (bool request_jacobian,
                                     DenseVector<Real>& f,
                                     DenseMatrix<Real>& jac) = 0;
        
        /*!
         *   damping force contribution to system residual
         */
        virtual bool damping_force (bool request_jacobian,
                                    DenseVector<Real>& f,
                                    DenseMatrix<Real>& jac);
        
        /*!
         *   inertial force contribution to system residual
         */
        virtual bool inertial_force (bool request_jacobian,
                                     DenseVector<Real>& f,
                                     DenseMatrix<Real>& jac);
        
        /*!
         *   side external force contribution to system residual
         */
        virtual bool side_external_force (bool request_jacobian,
                                          DenseVector<Number>& f,
                                          DenseMatrix<Number>& jac,
                                          std::multimap<boundary_id_type, MAST::BoundaryCondition*>& bc);
        
        /*!
         *   prestress force contribution to system residual
         */
        virtual bool prestress_force (bool request_jacobian,
                                      DenseVector<Real>& f,
                                      DenseMatrix<Real>& jac) = 0;
        
        /*!
         *   volume external force contribution to system residual
         */
        virtual bool volume_external_force (bool request_jacobian,
                                            DenseVector<Number>& f,
                                            DenseMatrix<Number>& jac,
                                            std::multimap<subdomain_id_type, MAST::BoundaryCondition*>& bc);

        /*!
         *   sensitivity of the internal force contribution to system residual
         */
        virtual bool internal_force_sensitivity (bool request_jacobian,
                                                 DenseVector<Real>& f,
                                                 DenseMatrix<Real>& jac) = 0;
        /*!
         *   sensitivity of the damping force contribution to system residual
         */
        virtual bool damping_force_sensitivity (bool request_jacobian,
                                                DenseVector<Real>& f,
                                                DenseMatrix<Real>& jac);
        
        /*!
         *   sensitivity of the inertial force contribution to system residual
         */
        virtual bool inertial_force_sensitivity (bool request_jacobian,
                                                 DenseVector<Real>& f,
                                                 DenseMatrix<Real>& jac);
        
        /*!
         *   sensitivity of the side external force contribution to system residual
         */
        virtual bool side_external_force_sensitivity (bool request_jacobian,
                                                      DenseVector<Number>& f,
                                                      DenseMatrix<Number>& jac,
                                                      std::multimap<boundary_id_type, MAST::BoundaryCondition*>& bc);
        
        /*!
         *   sensitivity of the prestress force contribution to system residual
         */
        virtual bool prestress_force_sensitivity (bool request_jacobian,
                                                  DenseVector<Real>& f,
                                                  DenseMatrix<Real>& jac) = 0;
        
        /*!
         *   sensitivity of the volume external force contribution to system residual
         */
        virtual bool volume_external_force_sensitivity (bool request_jacobian,
                                                        DenseVector<Number>& f,
                                                        DenseMatrix<Number>& jac,
                                                        std::multimap<subdomain_id_type, MAST::BoundaryCondition*>& bc);

        /*!
         *   returns the value of maximum von Mises stress over the element
         */
        virtual Real max_von_mises_stress() = 0;
        
        
        /*!
         *   returns the sensitivity of maximum von Mises stress over the element
         */
        virtual Real max_von_mises_stress_sensitivity() = 0;
        
        
        /*!
         *    flag for follower forces
         */
        bool follower_forces;
        
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
        
        /*!
         *   parameters for which sensitivity has to be calculated.
         */
        const MAST::SensitivityParameters* sensitivity_params;
        
        
        template <typename ValType>
        void transform_to_global_system(const DenseMatrix<ValType>& local_mat,
                                                DenseMatrix<ValType>& global_mat) const;
        
        template <typename ValType>
        void transform_to_local_system(const DenseVector<ValType>& global_vec,
                                       DenseVector<ValType>& local_vec) const;
        
        template <typename ValType>
        void transform_to_global_system(const DenseVector<ValType>& local_vec,
                                        DenseVector<ValType>& global_vec) const;
        
    protected:

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
        virtual void _get_side_fe_and_qrule(const Elem& e,
                                            unsigned int s,
                                            std::auto_ptr<FEBase>& fe,
                                            std::auto_ptr<QBase>& qrule);


        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_pressure_force(bool request_jacobian,
                                            DenseVector<Real>& f,
                                            DenseMatrix<Real>& jac,
                                            const unsigned int side,
                                            MAST::BoundaryCondition& p);

        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure 
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual bool surface_pressure_force(bool request_jacobian,
                                            DenseVector<Real>& f,
                                            DenseMatrix<Real>& jac,
                                            MAST::BoundaryCondition& p);

        
        /*!
         *    Calculates the force vector and Jacobian due to small 
         *    perturbation surface pressure.
         */
        virtual bool small_disturbance_surface_pressure_force(bool request_jacobian,
                                                              DenseVector<Number>& f,
                                                              DenseMatrix<Number>& jac,
                                                              const unsigned int side,
                                                              MAST::BoundaryCondition& p);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual bool small_disturbance_surface_pressure_force(bool request_jacobian,
                                                              DenseVector<Number>& f,
                                                              DenseMatrix<Number>& jac,
                                                              MAST::BoundaryCondition& p);

        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool surface_pressure_force_sensitivity(bool request_jacobian,
                                                        DenseVector<Real>& f,
                                                        DenseMatrix<Real>& jac,
                                                        const unsigned int side,
                                                        MAST::BoundaryCondition& p);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_pressure_force_sensitivity(bool request_jacobian,
                                                        DenseVector<Real>& f,
                                                        DenseMatrix<Real>& jac,
                                                        MAST::BoundaryCondition& p);

        
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


