#include "fvm_ghost_boundary.h"

#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void GhostCreationFromMesh::addGhostParticleAndSetInConfiguration()
{
    ghost_bound_.second = ghost_bound_.first;

    for (size_t index_i = 0; index_i != particles_->TotalRealParticles(); ++index_i)
    {
        for (size_t neighbor_index = 0; neighbor_index != mesh_topology_[index_i].size(); ++neighbor_index)
        {
            size_t boundary_type = mesh_topology_[index_i][neighbor_index][1];
            if (mesh_topology_[index_i][neighbor_index][1] != 2)
            {
                mutex_create_ghost_particle_.lock();
                size_t ghost_particle_index = ghost_bound_.second;
                ghost_bound_.second++;
                ghost_boundary_.checkWithinGhostSize(ghost_bound_);

                particles_->updateGhostParticle(ghost_particle_index, index_i);
                size_t node1_index = mesh_topology_[index_i][neighbor_index][2];
                size_t node2_index = mesh_topology_[index_i][neighbor_index][3];
                Vecd node1_position = node_coordinates_[node1_index];
                Vecd node2_position = node_coordinates_[node2_index];
                Vecd ghost_particle_position = 0.5 * (node1_position + node2_position);

                mesh_topology_[index_i][neighbor_index][0] = ghost_particle_index + 1;
                pos_[ghost_particle_index] = ghost_particle_position;
                mutex_create_ghost_particle_.unlock();

                std::vector<std::vector<size_t>> new_element;
                // Add (corresponding_index_i,boundary_type,node1_index,node2_index) to the new element
                std::vector<size_t> sub_element1 = {index_i + 1, boundary_type, node1_index, node2_index};
                new_element.push_back(sub_element1);
                // Add (corresponding_index_i,boundary_type,node1_index,node2_index) to the new element
                std::vector<size_t> sub_element2 = {index_i + 1, boundary_type, node1_index, node2_index};
                new_element.push_back(sub_element2);
                // Add (corresponding_index_i,boundary_type,node1_index,node2_index) to the new element
                std::vector<size_t> sub_element3 = {index_i + 1, boundary_type, node1_index, node2_index};
                new_element.push_back(sub_element3);
                // Add the new element to mesh_topology_
                mesh_topology_.push_back(new_element);
                // creating the boundary files with ghost particle index
                each_boundary_type_with_all_ghosts_index_[boundary_type].push_back(ghost_particle_index);
                // creating the boundary files with contact real particle index
                each_boundary_type_contact_real_index_[boundary_type].push_back(index_i);

                // creating the boundary files with ghost eij
                Vecd interface_area_vector = node1_position - node2_position;
                Real interface_area_size = interface_area_vector.norm();
                Vecd unit_vector = interface_area_vector / interface_area_size;
                // normal unit vector
                Vecd normal_vector = Vecd(unit_vector[1], -unit_vector[0]);
                // judge the direction
                Vecd particle_position = pos_[index_i];
                Vecd node1_to_center_direction = particle_position - node1_position;
                if (node1_to_center_direction.dot(normal_vector) < 0)
                {
                    normal_vector = -normal_vector;
                };
                each_boundary_type_with_all_ghosts_eij_[boundary_type].push_back(normal_vector);
            }
        }
    }
}
//=================================================================================================//
void BoundaryConditionSetupInFVM::resetBoundaryConditions()
{
    for (size_t boundary_type = 0; boundary_type < each_boundary_type_with_all_ghosts_index_.size(); ++boundary_type)
    {
        if (!each_boundary_type_with_all_ghosts_index_[boundary_type].empty())
        {
            for (size_t ghost_number = 0; ghost_number != each_boundary_type_with_all_ghosts_index_[boundary_type].size(); ++ghost_number)
            {
                size_t ghost_index = each_boundary_type_with_all_ghosts_index_[boundary_type][ghost_number];
                size_t index_i = each_boundary_type_contact_real_index_[boundary_type][ghost_number];
                Vecd e_ij = each_boundary_type_with_all_ghosts_eij_[boundary_type][ghost_number];

                // Dispatch the appropriate boundary condition
                switch (boundary_type)
                {
                case 3: // this refer to the different types of wall boundary conditions
                    applyNonSlipWallBoundary(ghost_index, index_i);
                    applyReflectiveWallBoundary(ghost_index, index_i, e_ij);
                    break;
                case 4:
                    applyTopBoundary(ghost_index, index_i);
                    break;
                case 5:
                    applyPressureOutletBC(ghost_index, index_i);
                    break;
                case 7:
                    applySymmetryBoundary(ghost_index, index_i, e_ij);
                    break;
                case 9:
                    applyFarFieldBoundary(ghost_index);
                    break;
                case 10:
                    applyGivenValueInletFlow(ghost_index);
                    applyVelocityInletFlow(ghost_index, index_i);
                    break;
                case 36:
                    applyOutletBoundary(ghost_index, index_i);
                    break;
                }
            }
        }
    }
}
//=================================================================================================//
} // namespace SPH
