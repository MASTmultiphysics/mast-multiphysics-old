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
    class FieldFunctionBase;
    class Temperature;
    class BoundaryCondition;
    class LocalElemBase;
    
    enum StructuralElemType {
        ELASTIC_ELEMENT_1D,
        ELASTIC_ELEMENT_2D,
        ELASTIC_ELEMENT_3D
    };
    
    
    /*!
     *   stress tensor
     */
    class Stress: public libMesh::DenseMatrix<libMesh::Real> {
    public:
        Stress():
        libMesh::DenseMatrix<libMesh::Real>()
        {
            this->resize(3,3);
        }
        
        virtual ~Stress() { }

        /*!
         *    calculates and returns the von Mises stress for this stress
         *    tensor
         */
        libMesh::Real von_mises_stress() const {
            
            const libMesh::DenseMatrix<libMesh::Real>& s = *this;
            
            libMesh::Real val =
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
        StructuralElementBase(libMesh::System& sys,
                              const libMesh::Elem& elem,
                              const MAST::ElementPropertyCardBase& p);
        
        virtual ~StructuralElementBase();
        
        
        /*!
         *   returns a constant reference to the finite element object
         */
        libMesh::System& system()  {
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
        const libMesh::Elem& elem() const {
            return _elem;
        }

        
        /*!
         *   returns a constant reference to the element in local coordinate system
         */
        virtual const MAST::LocalElemBase& local_elem() const = 0;

        
        /*!
         *    Local elements are defined for 1D and 2D elements that exist in 
         *    3D space. These elements have a local coordinate system associated
         *    with the local coordinate. This method accepts the point defined 
         *    the local coordinate system as the input and maps it to the 
         *    global coordinate system.
         */
        virtual void global_coordinates(const libMesh::Point& local,
                                        libMesh::Point& global) const = 0;
        
        
        /*!
         *   returns a constant reference to the element in local coordinate system
         */
        virtual const libMesh::Elem& get_elem_for_quadrature() const = 0;

        
        /*!
         *   returns a constant reference to the quadrature rule
         */
        libMesh::QBase& quadrature_rule()  {
            return *_qrule;
        }

        
        
        /*!
         *   returns a constant reference to the finite element object
         */
        libMesh::FEBase& fe()  {
            return *_fe;
        }

        
        /*!
         *   internal force contribution to system residual
         */
        virtual bool internal_force (bool request_jacobian,
                                     libMesh::DenseVector<libMesh::Real>& f,
                                     libMesh::DenseMatrix<libMesh::Real>& jac,
                                     bool if_ignore_ho_jac) = 0;
        
        /*!
         *   damping force contribution to system residual
         */
        virtual bool damping_force (bool request_jacobian,
                                    libMesh::DenseVector<libMesh::Real>& f,
                                    libMesh::DenseMatrix<libMesh::Real>& jac);
        
        /*!
         *   inertial force contribution to system residual
         */
        virtual bool inertial_force (bool request_jacobian,
                                     libMesh::DenseVector<libMesh::Real>& f,
                                     libMesh::DenseMatrix<libMesh::Real>& jac);
        
        /*!
         *   side external force contribution to system residual
         */
        virtual bool side_external_force (bool request_jacobian,
                                          libMesh::DenseVector<libMesh::Number>& f,
                                          libMesh::DenseMatrix<libMesh::Number>& jac,
                                          std::multimap<libMesh::boundary_id_type, MAST::BoundaryCondition*>& bc);
        
        /*!
         *   prestress force contribution to system residual
         */
        virtual bool prestress_force (bool request_jacobian,
                                      libMesh::DenseVector<libMesh::Real>& f,
                                      libMesh::DenseMatrix<libMesh::Real>& jac) = 0;
        
        /*!
         *   volume external force contribution to system residual
         */
        virtual bool volume_external_force (bool request_jacobian,
                                            libMesh::DenseVector<libMesh::Number>& f,
                                            libMesh::DenseMatrix<libMesh::Number>& jac,
                                            std::multimap<libMesh::subdomain_id_type, MAST::BoundaryCondition*>& bc);

        /*!
         *   sensitivity of the internal force contribution to system residual
         */
        virtual bool internal_force_sensitivity (bool request_jacobian,
                                                 libMesh::DenseVector<libMesh::Real>& f,
                                                 libMesh::DenseMatrix<libMesh::Real>& jac,
                                                 bool if_ignore_ho_jac) = 0;
        /*!
         *   sensitivity of the damping force contribution to system residual
         */
        virtual bool damping_force_sensitivity (bool request_jacobian,
                                                libMesh::DenseVector<libMesh::Real>& f,
                                                libMesh::DenseMatrix<libMesh::Real>& jac);
        
        /*!
         *   sensitivity of the inertial force contribution to system residual
         */
        virtual bool inertial_force_sensitivity (bool request_jacobian,
                                                 libMesh::DenseVector<libMesh::Real>& f,
                                                 libMesh::DenseMatrix<libMesh::Real>& jac);
        
        /*!
         *   sensitivity of the side external force contribution to system residual
         */
        virtual bool side_external_force_sensitivity (bool request_jacobian,
                                                      libMesh::DenseVector<libMesh::Number>& f,
                                                      libMesh::DenseMatrix<libMesh::Number>& jac,
                                                      std::multimap<libMesh::boundary_id_type, MAST::BoundaryCondition*>& bc);
        
        /*!
         *   sensitivity of the prestress force contribution to system residual
         */
        virtual bool prestress_force_sensitivity (bool request_jacobian,
                                                  libMesh::DenseVector<libMesh::Real>& f,
                                                  libMesh::DenseMatrix<libMesh::Real>& jac) = 0;
        
        /*!
         *   sensitivity of the volume external force contribution to system residual
         */
        virtual bool volume_external_force_sensitivity (bool request_jacobian,
                                                        libMesh::DenseVector<libMesh::Number>& f,
                                                        libMesh::DenseMatrix<libMesh::Number>& jac,
                                                        std::multimap<libMesh::subdomain_id_type, MAST::BoundaryCondition*>& bc);

        /*!
         *   returns the value of maximum von Mises stress over the element
         */
        virtual libMesh::Real max_von_mises_stress() = 0;
        
        
        /*!
         *   returns the sensitivity of maximum von Mises stress over the element
         */
        virtual libMesh::Real max_von_mises_stress_sensitivity() = 0;
        
        
        /*!
         *    flag for follower forces
         */
        bool follower_forces;
        
        /*!
         *   element solution, and sensitivity in the local coordinate system
         */
        libMesh::DenseVector<libMesh::Real> local_solution, local_solution_sens;
        
        /*!
         *   element velocity, and sensitivity in the local coordinate system
         */
        libMesh::DenseVector<libMesh::Real> local_velocity, local_velocity_sens;
        
        /*!
         *   element acceleration, and sensitivity in the local coordinate system
         */
        libMesh::DenseVector<libMesh::Real> local_acceleration, local_acceleration_sens;
        
        /*!
         *   parameter for which sensitivity has to be calculated.
         */
        const MAST::FieldFunctionBase* sensitivity_param;
        
        
        template <typename ValType>
        void transform_to_global_system(const libMesh::DenseMatrix<ValType>& local_mat,
                                                libMesh::DenseMatrix<ValType>& global_mat) const;
        
        template <typename ValType>
        void transform_to_local_system(const libMesh::DenseVector<ValType>& global_vec,
                                       libMesh::DenseVector<ValType>& local_vec) const;
        
        template <typename ValType>
        void transform_to_global_system(const libMesh::DenseVector<ValType>& local_vec,
                                        libMesh::DenseVector<ValType>& global_vec) const;
        
    protected:

        /*!
         *   matrix that transforms the global dofs to the local element coordinate
         *   system
         */
        virtual const libMesh::DenseMatrix<libMesh::Real>& _transformation_matrix() const = 0;
        
        /*!
         *   returns the quadrature and finite element for element integration. 
         *   These are raw pointers created using new. The pointers must be 
         *   deleted at the end of scope.
         */
        virtual void _init_fe_and_qrule(const libMesh::Elem& e);

        /*!
         *   returns the quadrature and finite element for element side 
         *   integration. These are raw pointers created using new. The 
         *   pointers must be deleted at the end of scope.
         */
        virtual void _get_side_fe_and_qrule(const libMesh::Elem& e,
                                            unsigned int s,
                                            std::auto_ptr<libMesh::FEBase>& fe,
                                            std::auto_ptr<libMesh::QBase>& qrule);


        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_pressure_force(bool request_jacobian,
                                            libMesh::DenseVector<libMesh::Real>& f,
                                            libMesh::DenseMatrix<libMesh::Real>& jac,
                                            const unsigned int side,
                                            MAST::BoundaryCondition& p);

        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure 
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual bool surface_pressure_force(bool request_jacobian,
                                            libMesh::DenseVector<libMesh::Real>& f,
                                            libMesh::DenseMatrix<libMesh::Real>& jac,
                                            MAST::BoundaryCondition& p);

        
        /*!
         *    Calculates the force vector and Jacobian due to small 
         *    perturbation surface pressure.
         */
        virtual bool small_disturbance_surface_pressure_force(bool request_jacobian,
                                                              libMesh::DenseVector<libMesh::Number>& f,
                                                              libMesh::DenseMatrix<libMesh::Number>& jac,
                                                              const unsigned int side,
                                                              MAST::BoundaryCondition& p);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure
         *    applied on the entire element domain. This is applicable for
         *    only 1D and 2D elements.
         */
        virtual bool small_disturbance_surface_pressure_force(bool request_jacobian,
                                                              libMesh::DenseVector<libMesh::Number>& f,
                                                              libMesh::DenseMatrix<libMesh::Number>& jac,
                                                              MAST::BoundaryCondition& p);

        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         */
        virtual bool surface_pressure_force_sensitivity(bool request_jacobian,
                                                        libMesh::DenseVector<libMesh::Real>& f,
                                                        libMesh::DenseMatrix<libMesh::Real>& jac,
                                                        const unsigned int side,
                                                        MAST::BoundaryCondition& p);
        
        
        /*!
         *    Calculates the force vector and Jacobian due to surface pressure.
         *    this should be implemented for each element type
         */
        virtual bool surface_pressure_force_sensitivity(bool request_jacobian,
                                                        libMesh::DenseVector<libMesh::Real>& f,
                                                        libMesh::DenseMatrix<libMesh::Real>& jac,
                                                        MAST::BoundaryCondition& p);

        
        /*!
         *    Calculates the force vector and Jacobian due to thermal stresses. 
         *    this should be implemented for each element type
         */
        virtual bool thermal_force(bool request_jacobian,
                                   libMesh::DenseVector<libMesh::Real>& f,
                                   libMesh::DenseMatrix<libMesh::Real>& jac,
                                   MAST::BoundaryCondition& p) = 0;

        /*!
         *    Calculates the sensitivity of force vector and Jacobian due to 
         *    thermal stresses. this should be implemented for each element type
         */
        virtual bool thermal_force_sensitivity(bool request_jacobian,
                                               libMesh::DenseVector<libMesh::Real>& f,
                                               libMesh::DenseMatrix<libMesh::Real>& jac,
                                               MAST::BoundaryCondition& p) = 0;

        /*!
         *    System to which this system belongs
         */
        libMesh::System& _system;
        
        /*!
         *   geometric element for which the computations are performed
         */
        const libMesh::Elem& _elem;
        
        /*!
         *   element property
         */
        const MAST::ElementPropertyCardBase& _property;
        
        /*!
         *   element finite element for computations
         */
        std::auto_ptr<libMesh::FEBase> _fe;

        /*!
         *   element quadrature rule for computations
         */
        std::auto_ptr<libMesh::QBase> _qrule;
    };

    /*!
     *    builds the structural element for the specified element type
     */
    std::auto_ptr<MAST::StructuralElementBase>
    build_structural_element(libMesh::System& sys,
                             const libMesh::Elem& elem,
                             const MAST::ElementPropertyCardBase& p);
}


#endif // __mast_structural_element_base_h__


