/*
Copyright (c) 2005-2019, University of Oxford.
Copyright (c) 2025, Oliver M. Drozdowski and Ulrich S. Schwarz (Heidelberg University)

All rights reserved.
University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.
This file is part of Chaste.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "SurfaceTensionForce.hpp"
#include "LinearSystem.hpp"
#include "MonolayerVertexElement.hpp"
#include "PetscVecTools.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
SurfaceTensionForce<DIM>::SurfaceTensionForce()
        : AbstractForce<DIM>(),
          mSurfaceTensionParameters(), // Initialize as empty vector, then use 1.0 as standard value
          mProjectOnVolumeConservingForces(true),
          mPerformActiveT1Swaps(false),
          mRelativeNoiseStrength(0.0),
          mAnnealingDecayTime(1.0),
          mInitialTemperature(1e-10),
          mApicalStandardTension(0.0),
          mBasalStandardTension(0.0),
          mLateralStandardTension(0.0),
          mT1AreTemperatureDependent(true),
          mT1TransitionRate(100.0),
          mModeOfSurfaceTensionCalculation(0),
          mIndividualNodeDampingActivated(false),
          mpSimulation(nullptr)
{
}

template <unsigned DIM>
SurfaceTensionForce<DIM>::~SurfaceTensionForce()
{
}

template <unsigned DIM>
void SurfaceTensionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    ///\todo: check whether this line influences profiling tests - if so, we should remove it.
    if (dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("SurfaceTensionForce is to be used with a 3D MonolayerVertexBasedCellPopulation only");
    }

    // Define some helper variables
    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableMonolayerVertexMesh<DIM, DIM>& r_mesh = p_cell_population->rGetMesh();
    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();
    unsigned num_faces = p_cell_population->rGetMesh().GetNumFaces();

    // First we update the surface tensions, which is necessary after cell division or remeshing
    // with different mutation states
    UpdateSurfaceTensions(&(rCellPopulation));

    // If we do not have any Surface Tensions use a standard value
    if (mSurfaceTensionParameters.empty())
    {
        mSurfaceTensionParameters = std::vector<double>(num_faces, 1.0);
    }

    // Begin by computing the areas and volume of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_volumes(num_elements);
    std::map<MonolayerVertexElement<DIM - 1, DIM>*, double> faces_areas_old;
    std::vector<double> target_volumes(num_elements);
    std::vector<double> cell_types(num_elements);

    // For volume conservation we determine a volume conservation matrix
    // This matrix contains the sum of scalar products of valume gradients
    // as entries and defines the linear system to determine the force projections
    LinearSystem* p_system_volume_correction = new LinearSystem(num_elements, num_elements);
    LinearSystem* p_system_volume_conservation = new LinearSystem(num_elements, num_elements);
    LinearSystem* p_system_random_force = new LinearSystem(num_elements, num_elements);

    // Matrix consists of symmetric scalar products
    // p_system_volume_correction->SetMatrixIsSymmetric();
    // p_system_volume_conservation->SetMatrixIsSymmetric();

    // p_system_volume_correction->AssembleIntermediateLinearSystem();
    // p_system_volume_conservation->AssembleIntermediateLinearSystem();

    // p_system_volume_correction->ZeroLinearSystem();
    // p_system_volume_conservation->ZeroLinearSystem();

    for (typename MonolayerVertexMesh<DIM, DIM>::MonolayerVertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        element_volumes[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
        try
        {
            // If we haven't specified a growth modifier, there won't be any target areas in the CellData array and CellData
            // will throw an exception that it doesn't have "target area" entries.  We add this piece of code to give a more
            // understandable message. There is a slight chance that the exception is thrown although the error is not about the
            // target areas.
            target_volumes[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target volume");
            // cell_types[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("cell type");
        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetVolumeModifier to the simulation in order to use a SurfaceTensionForce");
        }

        // Add the excess volume to the RHS of system for volume correction
        double volume_excess = element_volumes[elem_index] - target_volumes[elem_index];
        p_system_volume_correction->AddToRhsVectorElement(elem_index, -volume_excess);
    }

    // Uncorrected forces need to be saved for the projection
    std::vector<c_vector<double, DIM> > forces_on_nodes;

    // Save all volume gradients, outer index nodes, inner index containing cells
    std::vector<std::vector<c_vector<double, DIM> > > vector_volume_gradients_all_nodes;

    // Store the initial node positions (these may be needed when applying random forces)
    std::map<Node<DIM>*, c_vector<double, DIM> > old_node_locations;

    // Random forces are saved
    std::map<Node<DIM>*, c_vector<double, DIM> > random_forces;

    // Set up everything for the simulated annealing
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();
    double temperature;
    if (current_time < mAnnealingDecayTime)
    {
        temperature = mInitialTemperature * (1.0 - current_time / mAnnealingDecayTime);
    }
    else
    {
        temperature = 1.0e-10;
    }
    double relative_strength = mRelativeNoiseStrength;
    RandomNumberGenerator* p_random = RandomNumberGenerator::Instance();

    // Set the probability for Active T1s
    if (DoPerformActiveT1Swaps())
    {
        double time_step = p_simulation_time->GetTimeStep();
        double probability = 0.0;
        if (IsT1TransitonRateTemperatureDependent())
            probability = temperature * time_step * mT1TransitionRate / r_mesh.GetNumFaces();
        else
            probability = time_step * mT1TransitionRate / r_mesh.GetNumFaces();
        r_mesh.SetActiveT1SwapProbability(probability);
    }
    else
    {
        r_mesh.SetActiveT1SwapProbability(0.0);
    }

    // Iterate over vertices in the cell population
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
        old_node_locations[p_this_node] = p_this_node->rGetLocation();

        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the CellPopulation, evaluated at the position of the vertex.
         *
         *
         *
         *
         *
         * Note that since the movement of this Node only affects the free energy
         * of the CellPtrs containing it, we can just consider the contributions
         * to the free energy gradient from each of these CellPtrs.
         */
        c_vector<double, DIM> surface_tension_contribution = zero_vector<double>(DIM);
        std::vector<c_vector<double, DIM> > vector_volume_gradients;

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            MonolayerVertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned elem_index = p_element->GetIndex();
            // unsigned num_nodes_elem = p_element->GetNumNodes();
            unsigned num_faces_elem = p_element->GetNumFaces();
            unsigned local_node_index = p_element->GetNodeLocalIndex(node_index);

            // For each face in the elements which contain the node/vertex
            // we calculate the term from the area gradients
            for (unsigned index_face = 0; index_face < num_faces_elem; ++index_face)
            {
                MonolayerVertexElement<DIM - 1, DIM>* p_face = p_element->GetFace(index_face);
                unsigned face_index = p_face->GetIndex();
                unsigned local_node_index_in_face = p_face->GetNodeLocalIndex(node_index);
                // If the node is not even in the face, we just continue to the next face
                if (local_node_index_in_face == UINT_MAX)
                {
                    continue;
                }

                double surface_tension_parameter = mSurfaceTensionParameters[face_index];

                // Note that the 1/2 from counting faces twice is part of the surface tension parameter
                surface_tension_contribution -= surface_tension_parameter * r_mesh.GetAreaGradientOfFaceAtNode(p_face, local_node_index_in_face);
            }

            // We then determine the volume gradient, which we need to
            // calculate the volume conservation matrix
            c_vector<double, DIM> volume_gradient = r_mesh.GetVolumeGradientAtNode(elem_index, local_node_index);
            vector_volume_gradients.push_back(volume_gradient);
        }

        c_vector<double, DIM> force_on_node = surface_tension_contribution;
        forces_on_nodes.push_back(force_on_node);
        vector_volume_gradients_all_nodes.push_back(vector_volume_gradients);

        // We now create a random force
        c_vector<double, DIM> random_force_on_node = zero_vector<double>(DIM);

        for (unsigned i_dim = 0; i_dim < DIM; i_dim++)
        {
            // Add Gaussian noise to force with sigma = rel_strength*force
            double std_rand = p_random->StandardNormalRandomDeviate();
            std_rand *= relative_strength;
            random_force_on_node[i_dim] += std_rand;
        }

        random_forces[p_this_node] = random_force_on_node;

        // We check whether the node is in small area and should not
        // be moved fully in the volume correction step
        double small_area_correction_matrix_coefficient = 1.0;
        if (mIndividualNodeDampingActivated)
        {
            if (mMapNodeToDampingThresholds.find(p_cell_population->GetNode(node_index))
                != mMapNodeToDampingThresholds.end())
            {
                small_area_correction_matrix_coefficient = mMapNodeToDampingThresholds[p_this_node].volumeCorrectionDamping;
            }
        }

        // Set up the linear systems for volume conservation
        // Iterate over containing elements
        unsigned index_volume_gradient_larger = 0;
        for (std::set<unsigned>::iterator iter_elem_larger = containing_elem_indices.begin();
             iter_elem_larger != containing_elem_indices.end();
             ++iter_elem_larger)
        {
            MonolayerVertexElement<DIM, DIM>* p_element_larger = p_cell_population->GetElement(*iter_elem_larger);
            unsigned global_elem_index_larger = p_element_larger->GetIndex();

            c_vector<double, DIM> volume_gradient_1 = vector_volume_gradients[index_volume_gradient_larger];

            unsigned index_volume_gradient_smaller = 0;
            std::set<unsigned>::iterator iter_elem_larger_next = iter_elem_larger;
            iter_elem_larger_next++;
            for (std::set<unsigned>::iterator iter_elem_smaller = containing_elem_indices.begin();
                 iter_elem_smaller != iter_elem_larger_next;
                 ++iter_elem_smaller)
            {
                MonolayerVertexElement<DIM, DIM>* p_element_smaller = p_cell_population->GetElement(*iter_elem_smaller);
                unsigned global_elem_index_smaller = p_element_smaller->GetIndex();

                c_vector<double, DIM> volume_gradient_2 = vector_volume_gradients[index_volume_gradient_smaller];
                double product_volume_gradients = inner_prod(volume_gradient_1, volume_gradient_2);
                // std::cout << "\n" << global_elem_index_larger << " , " << global_elem_index_smaller << std::flush;

                // p_system_volume_correction->SwitchWriteModeLhsMatrix();
                // p_system_volume_conservation->SwitchWriteModeLhsMatrix();

                p_system_volume_correction->AddToMatrixElement(global_elem_index_larger, global_elem_index_smaller,
                                                               small_area_correction_matrix_coefficient * product_volume_gradients);
                if (global_elem_index_smaller != global_elem_index_larger)
                {
                    p_system_volume_correction->AddToMatrixElement(global_elem_index_smaller, global_elem_index_larger,
                                                                   small_area_correction_matrix_coefficient * product_volume_gradients);
                }

                p_system_volume_conservation->AddToMatrixElement(global_elem_index_larger, global_elem_index_smaller, product_volume_gradients);
                p_system_random_force->AddToMatrixElement(global_elem_index_larger, global_elem_index_smaller, product_volume_gradients);
                if (global_elem_index_smaller != global_elem_index_larger)
                {
                    p_system_volume_conservation->AddToMatrixElement(global_elem_index_smaller, global_elem_index_larger, product_volume_gradients);
                    p_system_random_force->AddToMatrixElement(global_elem_index_smaller, global_elem_index_larger, product_volume_gradients);
                }

                index_volume_gradient_smaller++;
            }
            // Perform diagonal element where iter_larger == iter_smaller
            /*if(index_volume_gradient_smaller == index_volume_gradient_larger)
            {
              MonolayerVertexElement<DIM, DIM>* p_element_smaller = p_cell_population->GetElement(*iter_elem_smaller);
              unsigned global_elem_index_smaller = p_element_smaller->GetIndex();

              c_vector<double, DIM> volume_gradient_2 = vector_volume_gradient[index_volume_gradient_smaller];
              double product_volume_gradients = inner_prod(volume_gradient_1, volume_gradient_2);

              p_system_volume_correction->AddToMatrixElement(global_elem_index_larger, global_elem_index_smaller, product_volume_gradients);
              p_system_volume_correction->AddToMatrixElement(global_elem_index_smaller, global_elem_index_larger, product_volume_gradients);

              p_system_volume_conservation->AddToMatrixElement(global_elem_index_larger, global_elem_index_smaller, product_volume_gradients);
              p_system_volume_conservation->AddToMatrixElement(global_elem_index_smaller, global_elem_index_larger, product_volume_gradients);
            }
            else
            {
              EXCEPTION("Something went wrong!");
            }*/

            // add force to rhs
            double product_force_vol_gradient = inner_prod(force_on_node, volume_gradient_1);
            p_system_volume_conservation->AddToRhsVectorElement(global_elem_index_larger, product_force_vol_gradient);
            // add random force to rhs
            double product_random_force_vol_gradient = inner_prod(random_force_on_node, volume_gradient_1);
            p_system_random_force->AddToRhsVectorElement(global_elem_index_larger, product_random_force_vol_gradient);

            // Flush entries
            p_system_volume_correction->AssembleIntermediateLinearSystem();
            p_system_volume_conservation->AssembleIntermediateLinearSystem();
            p_system_random_force->AssembleIntermediateLinearSystem();

            index_volume_gradient_larger++;
        }
        /*
        for(unsigned index_volume_gradient=0; index_volume_gradient < vector_volume_gradients.size(); ++index_volume_gradient)
        {
          c_vector<double, DIM> volume_gradient_1 = vector_volume_gradient[index_volume_gradient];

          for(unsigned index_smaller_volume_gradient=0;
              index_smaller_volume_gradient <= index_valume_gradient;
              ++index_smaller_volume_gradient)
          {
              c_vector<double, DIM> volume_gradient_2 = vector_volume_gradient[index_smaller_volume_gradient];
              double product_volume_gradients = inner_prod(volume_gradient_1, volume_gradient_2);

              system_volume_correction->AddToMatrixElement(index_volume_gradient, index_smaller_volume_gradient, product_volume_gradients);
              system_volume_correction->AddToMatrixElement(index_smaller_volume_gradient, index_volume_gradient, product_volume_gradients);

              system_volume_conservation->AddToMatrixElement(index_volume_gradient, index_smaller_volume_gradient, product_volume_gradients);
              system_volume_conservation->AddToMatrixElement(index_smaller_volume_gradient, index_volume_gradient, product_volume_gradients);
          }
          // add force to rhs
          double product_force_vol_gradient = inner_prod(force_on_node, volume_gradient_1);
          system_volume_conservation->AddToRhsVectorElement(index_volume_gradient, product_force_vol_gradient);
        }
        */
    }

    // Now that all forces and gradients have been calculated and inserted into the linear systems
    // we can solve for the parameters, which determine the volume correction (in the direction of
    // volume gradients) and the force projected on the conservation space (orthogonal to vol. grads.)
    p_system_volume_correction->AssembleFinalLinearSystem();
    p_system_volume_conservation->AssembleFinalLinearSystem();
    p_system_random_force->AssembleFinalLinearSystem();

    Vec parameters_correction = p_system_volume_correction->Solve();
    Vec parameters_conservation = p_system_volume_conservation->Solve();
    Vec parameters_random_force = nullptr;
    if (mRelativeNoiseStrength > 0.0)
    {
        parameters_random_force = p_system_random_force->Solve();
    }

    // Iterate over vertices in the cell population to calculate and save the force contributions
    // and perform the single-step volume correction
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        c_vector<double, DIM> correction_vector = zero_vector<double>(DIM);
        c_vector<double, DIM> force_on_node = forces_on_nodes[node_index];
        c_vector<double, DIM> random_force_on_node = random_forces[p_cell_population->GetNode(node_index)];

        // We perform the volume correction and the force projection in one step
        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        unsigned local_cell_index_for_node = 0;
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            MonolayerVertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned global_elem_index = p_element->GetIndex();

            // Determine cell-indexed vector for correction
            double coefficient_cell = PetscVecTools::GetElement(parameters_correction, global_elem_index);
            c_vector<double, DIM> gradient_vector_cell_node = vector_volume_gradients_all_nodes[node_index][local_cell_index_for_node];

            /* If we have node-specific damping, we need to threshold the correction with
             * the corresponding node-specific threshold, which comes from the
             * map of nodes to damping values. First we check whether we have node-
             * specific damping and then we reduce the magnitude of the correction
             * Note that this can break volume conservation! Since this only occurs
             * for small areas this should still converge to the correct volume
             * for time steps small enough
             */

            double threshold_correction_value = 1.0;
            if (mIndividualNodeDampingActivated)
            {
                if (mMapNodeToDampingThresholds.find(p_cell_population->GetNode(node_index))
                    != mMapNodeToDampingThresholds.end())
                {
                    threshold_correction_value = mMapNodeToDampingThresholds[p_cell_population->GetNode(node_index)].volumeCorrectionDamping;
                }
            }

            // Correct the node position
            c_vector<double, DIM>& node_position_vector = p_cell_population->GetNode(node_index)->rGetModifiableLocation();
            node_position_vector += threshold_correction_value * coefficient_cell * gradient_vector_cell_node;
            old_node_locations[p_cell_population->GetNode(node_index)] = node_position_vector;

            // Determine cell-indexed vectors for projection of forces
            double coefficient_forces = PetscVecTools::GetElement(parameters_conservation, global_elem_index);
            double coefficient_random_forces = 0.0;
            if (mRelativeNoiseStrength > 0.0)
            {
                coefficient_random_forces = PetscVecTools::GetElement(parameters_random_force, global_elem_index);
            }

            // Project the force
            force_on_node -= coefficient_forces * gradient_vector_cell_node;
            random_force_on_node -= coefficient_random_forces * gradient_vector_cell_node;

            local_cell_index_for_node++;
        }
        // Save the force into the node
        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);
        random_forces[p_cell_population->GetNode(node_index)] = random_force_on_node;
    }

    /* For simulated annealing we now have to calculate the energy difference we get from the
     * random force. If this energy difference is negative, we accept the step, if not we
     * accept with the Boltzmann probability with the current annealing temperature.
     *
     * For this we first calculate the energy, move the nodes according to our random
     * forces, calculate the new energy, move the nodes back and finally save the
     * force contribution if this step is accepted
     */

    if ((mRelativeNoiseStrength > 0.0) & (mpSimulation != nullptr))
    {
        double energy_before_random = CalculateTotalEnergy(rCellPopulation);
        double timestep = mpSimulation->GetDt();
        // Move the nodes (assuming Forward Euler Step)
        for (unsigned node_index = 0; node_index < num_nodes; node_index++)
        {
            Node<DIM>* p_node = p_cell_population->GetNode(node_index);
            double damping = p_cell_population->GetDampingConstant(node_index);
            c_vector<double, DIM> random_force_on_node = random_forces[p_node];
            c_vector<double, DIM>& location_node = p_node->rGetModifiableLocation();
            location_node += random_force_on_node * timestep / damping;
        }
        // Compute the new energy
        double energy_after_random = CalculateTotalEnergy(rCellPopulation);

        // Reset the nodes
        for (unsigned node_index = 0; node_index < num_nodes; node_index++)
        {
            Node<DIM>* p_node = p_cell_population->GetNode(node_index);
            c_vector<double, DIM>& location_node = p_node->rGetModifiableLocation();
            location_node = old_node_locations[p_node];
        }
        // Add the force contribution, depending on energy
        double uniform_rand = p_random->ranf();
        double boltzmann_factor = exp(-(energy_after_random - energy_before_random) / temperature);
        if (uniform_rand < boltzmann_factor)
        {
            for (unsigned node_index = 0; node_index < num_nodes; node_index++)
            {
                Node<DIM>* p_node = p_cell_population->GetNode(node_index);
                p_node->AddAppliedForceContribution(random_forces[p_node]);
            }
        }
    }

    /* If we have node-specific damping, we need to threshold the force with
     * the corresponding node-specific threshold, which comes from the
     * map of nodes to damping values. First we check whether we have node-
     * specific damping and then we reduce the magnitude of the force
     * Note that this can break volume conservation! Since this only occurs
     * for small areas this is however not such a large effect and we
     * correct for it in the correction step of the next iteration.
     */

    if (mIndividualNodeDampingActivated)
    {
        for (auto iter = mMapNodeToDampingThresholds.begin(); iter != mMapNodeToDampingThresholds.end(); ++iter)
        {
            Node<DIM>* p_node = iter->first;
            double threshold = iter->second.forceDamping;
            c_vector<double, DIM> force = p_node->rGetAppliedForce();
            double norm_force = norm_2(force);
            if (norm_force > threshold)
            {
                p_node->ClearAppliedForce();
                p_node->AddAppliedForceContribution(threshold / norm_force * force);
            }
        }
    }

    // Cleanup
    delete p_system_volume_conservation;
    delete p_system_volume_correction;
    delete p_system_random_force;

    // Must release the memory from the Petsc Vector to not have memory leaks
    if (parameters_correction)
    {
        PetscTools::Destroy(parameters_correction);
    }
    if (parameters_conservation)
    {
        PetscTools::Destroy(parameters_conservation);
    }
    if (parameters_random_force)
    {
        PetscTools::Destroy(parameters_random_force);
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SurfaceTensionForce<1>::AddForceContribution(AbstractCellPopulation<1>& rCellPopulation)
{
    EXCEPTION("SurfaceTensionForce not possible in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SurfaceTensionForce<1>::CreateSurfaceTensionParametersForCells(double apical_surface_tension,
                                                                    double basal_surface_tension,
                                                                    double lateral_surface_tension,
                                                                    MonolayerVertexMesh<1, 1>* p_mesh)
{
    EXCEPTION("CreateSurfaceTensionParametersForCells not implemented in one dimension.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
void SurfaceTensionForce<DIM>::CreateSurfaceTensionParametersForCells(double apical_surface_tension,
                                                                      double basal_surface_tension,
                                                                      double lateral_surface_tension,
                                                                      MonolayerVertexMesh<DIM, DIM>* p_mesh)
{
    // Set the mode to uniform
    mModeOfSurfaceTensionCalculation = 1;

    unsigned num_faces = p_mesh->GetNumFaces();
    if (mSurfaceTensionParameters.size() != num_faces)
    {
        mSurfaceTensionParameters = std::vector<double>(num_faces, 1.0);
    }

    for (unsigned faceIndex = 0; faceIndex < num_faces; ++faceIndex)
    {
        MonolayerVertexElement<DIM - 1, DIM>* p_face = p_mesh->GetFace(faceIndex);
        bool boundary_face = p_face->IsBoundaryFace();
        double boundary_factor = boundary_face ? 2.0 : 1.0;
        if (p_face->GetFaceType() == MonolayerVertexElementType::Apical)
        {
            mSurfaceTensionParameters[faceIndex] = apical_surface_tension;
        }
        else if (p_face->GetFaceType() == MonolayerVertexElementType::Basal)
        {
            mSurfaceTensionParameters[faceIndex] = basal_surface_tension;
        }
        else if (p_face->GetFaceType() == MonolayerVertexElementType::Lateral)
        {
            mSurfaceTensionParameters[faceIndex] = boundary_factor * lateral_surface_tension / 2.0;
        }
        else
        {
            mSurfaceTensionParameters[faceIndex] = boundary_factor * 0.0;
        }
    }
    // Save into internal variables
    mApicalStandardTension = apical_surface_tension;
    mBasalStandardTension = basal_surface_tension;
    mLateralStandardTension = lateral_surface_tension;
}

template <unsigned DIM>
double SurfaceTensionForce<DIM>::GetApicalStandardTension()
{
    return mApicalStandardTension;
}

template <unsigned DIM>
double SurfaceTensionForce<DIM>::GetBasalStandardTension()
{
    return mBasalStandardTension;
}

template <unsigned DIM>
double SurfaceTensionForce<DIM>::GetLateralStandardTension()
{
    return mLateralStandardTension;
}

template <unsigned DIM>
double SurfaceTensionForce<DIM>::GetSurfaceTensionParameter(unsigned face)
{
    if (face >= mSurfaceTensionParameters.size())
    {
        EXCEPTION("Error! Index to large in GetSurfaceTensionParameter");
    }
    return mSurfaceTensionParameters[face];
}

template <unsigned DIM>
std::vector<double> SurfaceTensionForce<DIM>::GetSurfaceTensionParameters()
{
    return mSurfaceTensionParameters;
}

template <unsigned DIM>
bool SurfaceTensionForce<DIM>::DoProjectOnVolumeConservingForces()
{
    return true;
}

template <unsigned DIM>
bool SurfaceTensionForce<DIM>::DoPerformActiveT1Swaps()
{
    return mPerformActiveT1Swaps;
}

template <unsigned DIM>
double SurfaceTensionForce<DIM>::GetVolumeElasticityParameter()
{
    return mVolumeElasticityParameter;
}

template <unsigned DIM>
void SurfaceTensionForce<DIM>::SetSurfaceTensionParameter(unsigned face, double surfaceTensionParameter)
{
    // Set mode to matrix
    mModeOfSurfaceTensionCalculation = 2;

    if (face >= mSurfaceTensionParameters.size())
    {
        EXCEPTION("Error! Index to large in GetSurfaceTensionParameter");
    }
    mSurfaceTensionParameters[face] = surfaceTensionParameter;
}

template <unsigned DIM>
void SurfaceTensionForce<DIM>::SetSurfaceTensionParametersByMutation(std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3> > mappingMutationTension)
{
    // Set mode to dictionary
    mModeOfSurfaceTensionCalculation = 3;

    mMappingMutationTension = mappingMutationTension;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SurfaceTensionForce<1>::UpdateSurfaceTensions(AbstractCellPopulation<1>* pCellPopulation)
{
    EXCEPTION("SurfaceTensionForce not possible in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
void SurfaceTensionForce<DIM>::UpdateSurfaceTensions(AbstractCellPopulation<DIM>* pCellPopulation)
{
    if (mModeOfSurfaceTensionCalculation == 1)
    {
        MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(pCellPopulation);

        // Update with the old values.
        CreateSurfaceTensionParametersForCells(mApicalStandardTension, mBasalStandardTension,
                                               mLateralStandardTension, &(p_cell_population->rGetMesh()));
    }
    else if (mModeOfSurfaceTensionCalculation == 3)
    {
        MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(pCellPopulation);

        // Update with the old mapping to tensions
        UpdateSurfaceTensionsByMutation(p_cell_population);
    }
    // In the matrix case we cannot consistently update.
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SurfaceTensionForce<1>::UpdateSurfaceTensionsByMutation(AbstractCellPopulation<1>* pCellPopulation)
{
    EXCEPTION("SurfaceTensionForce not possible in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
void SurfaceTensionForce<DIM>::UpdateSurfaceTensionsByMutation(AbstractCellPopulation<DIM>* pCellPopulation)
{
    assert(mModeOfSurfaceTensionCalculation == 3);

    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(pCellPopulation);

    // We initialize an empty surface tension parameter vector
    unsigned num_faces = p_cell_population->rGetMesh().GetNumFaces();
    mSurfaceTensionParameters = std::vector<double>(num_faces, 0.0);

    for (typename MonolayerVertexMesh<DIM, DIM>::MonolayerVertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        // Get the tensions for the cell with the correpsonding mutation
        CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(elem_iter->GetIndex());
        boost::shared_ptr<AbstractCellMutationState> p_mutation = p_cell->GetMutationState();
        std::array<double, 3> surface_tensions_cell = mMappingMutationTension[p_mutation];
        double apical_surface_tension = surface_tensions_cell[0];
        double basal_surface_tension = surface_tensions_cell[1];
        double lateral_surface_tension = surface_tensions_cell[2];

        unsigned num_faces_in_element = elem_iter->GetNumFaces();
        for (unsigned face_index = 0; face_index < num_faces_in_element; ++face_index)
        {
            MonolayerVertexElement<DIM - 1, DIM>* p_face = elem_iter->GetFace(face_index);
            bool boundary_face = p_face->IsBoundaryFace();
            unsigned global_face_index = p_face->GetIndex();
            assert(global_face_index < num_faces);

            // If it is a shared face, we take the average of the cell-specific face values
            // In the sum of forces we also go twice over this face. Therefore multiply by 1/4.
            // For boundaries we remultiply by 4, because this face is only shared by ONE element!
            // Therefore we need to save it as 1.0*x and don't average
            double boundary_factor = boundary_face ? 4.0 : 1.0;

            if (p_face->GetFaceType() == MonolayerVertexElementType::Apical)
            {
                mSurfaceTensionParameters[global_face_index] += apical_surface_tension;
            }
            else if (p_face->GetFaceType() == MonolayerVertexElementType::Basal)
            {
                mSurfaceTensionParameters[global_face_index] += basal_surface_tension;
            }
            else if (p_face->GetFaceType() == MonolayerVertexElementType::Lateral)
            {
                // We sum over each face twice and divide by 2.0 to get the average.
                mSurfaceTensionParameters[global_face_index] += boundary_factor * lateral_surface_tension / 4.0;
            }
            else
            {
                mSurfaceTensionParameters[global_face_index] += boundary_factor * 0.0;
            }
        }
    }
}

template <unsigned DIM>
void SurfaceTensionForce<DIM>::SetSurfaceTensionParameters(std::vector<double> surfaceTensionParameters)
{
    mSurfaceTensionParameters = surfaceTensionParameters;
}

template <unsigned DIM>
void SurfaceTensionForce<DIM>::SetVolumeElasticityParameter(double volumeElasticityParameter)
{
    mVolumeElasticityParameter = volumeElasticityParameter;
}

template <unsigned DIM>
void SurfaceTensionForce<DIM>::SetT1TransitionParameters(double T1TransitionRate, bool T1AreTemperatureDependent)
{
    mT1TransitionRate = T1TransitionRate;
    mT1AreTemperatureDependent = mT1AreTemperatureDependent;
}

template <unsigned DIM>
double SurfaceTensionForce<DIM>::GetT1TransitionRate()
{
    return mT1TransitionRate;
}

template <unsigned DIM>
bool SurfaceTensionForce<DIM>::IsT1TransitonRateTemperatureDependent()
{
    return mT1AreTemperatureDependent;
}

template <unsigned DIM>
void SurfaceTensionForce<DIM>::SetProjectOnVolumeConservingForces(bool doProjection)
{
    if (!doProjection)
    {
        EXCEPTION("We have not yet implement SurfaceTensionForce without Volume conservation projection!");
    }
}

template <unsigned DIM>
void SurfaceTensionForce<DIM>::SetPerformActiveT1Swaps(bool doActiveT1Swaps)
{
    mPerformActiveT1Swaps = doActiveT1Swaps;
}

template <>
double SurfaceTensionForce<1>::CalculateTotalEnergy(AbstractCellPopulation<1>& rCellPopulation)
{
    EXCEPTION("SurfaceTensionForce not possible in 1D.");
}

template <unsigned DIM>
double SurfaceTensionForce<DIM>::CalculateTotalEnergy(AbstractCellPopulation<DIM>& rCellPopulation)
{
    double total_energy = 0.0;

    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_faces = p_cell_population->rGetMesh().GetNumFaces();

    // If we do not have any Surface Tensions use a standard value
    if (mSurfaceTensionParameters.empty())
    {
        mSurfaceTensionParameters = std::vector<double>(num_faces, 1.0);
    }

    for (typename MonolayerVertexMesh<DIM, DIM>::MonolayerVertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned num_faces_in_element = elem_iter->GetNumFaces();
        for (unsigned face_index = 0; face_index < num_faces_in_element; ++face_index)
        {
            MonolayerVertexElement<DIM - 1, DIM>* p_face = elem_iter->GetFace(face_index);
            unsigned global_face_index = p_face->GetIndex();
            double face_area = p_cell_population->rGetMesh().CalculateAreaOfFace(p_face);
            total_energy += mSurfaceTensionParameters[global_face_index] * face_area;
        }
    }

    return total_energy;
}

template <unsigned DIM>
void SurfaceTensionForce<DIM>::SetSimulationInstance(AbstractCellBasedSimulation<DIM, DIM>* pSimulation)
{
    mpSimulation = pSimulation;
}

template <unsigned DIM>
void SurfaceTensionForce<DIM>::SetSimulatedAnnealingParameters(double relativeNoiseStrength,
                                                               double annealingDecayTime,
                                                               double initialTemperature)
{
    mRelativeNoiseStrength = relativeNoiseStrength;
    mAnnealingDecayTime = annealingDecayTime;
    mInitialTemperature = initialTemperature;
}

template <unsigned DIM>
void SurfaceTensionForce<DIM>::SetMapNodeToDampingThresholds(
    std::map<Node<DIM>*, DampingCoefficientPair> mapNodeToDampingThreshold)
{
    mMapNodeToDampingThresholds = mapNodeToDampingThreshold;
    mIndividualNodeDampingActivated = true;
}

template <unsigned DIM>
void SurfaceTensionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    if (mSurfaceTensionParameters.empty())
    {
        *rParamsFile << "\t\t\t<SurfaceTensionParameterZero>" << 1.0 << "</SurfaceTensionParametersZero>\n";
    }
    else
    {
        *rParamsFile << "\t\t\t<SurfaceTensionParameterZero>" << mSurfaceTensionParameters[0] << "</SurfaceTensionParametersZero>\n";
    }
    *rParamsFile << "\t\t\t<ModeOfSurfaceTensionCalculation>" << mModeOfSurfaceTensionCalculation << "</ModeOfSurfaceTensionCalculation>\n";

    *rParamsFile << "\t\t\t<StandardApicalSurfaceTension>" << mApicalStandardTension << "</StandardApicalSurfaceTension>\n";
    *rParamsFile << "\t\t\t<StandardBasalSurfaceTension>" << mBasalStandardTension << "</StandardBasalSurfaceTension>\n";
    *rParamsFile << "\t\t\t<StandardLateralSurfaceTension>" << mLateralStandardTension << "</StandardLateralSurfaceTension>\n";

    *rParamsFile << "\t\t\t<MappingMutationTension>\n";
    for (auto it = mMappingMutationTension.begin(); it != mMappingMutationTension.end(); ++it)
    {
        *rParamsFile << "\t\t\t\t<MutationCellType>" << it->first->GetColour() << "</MutationCellType>\n";
        *rParamsFile << "\t\t\t\t<ApicalSurfaceTension>" << it->second[0] << "</ApicalSurfaceTension>\n";
        *rParamsFile << "\t\t\t\t<BasalSurfaceTension>" << it->second[1] << "</BasalSurfaceTension>\n";
        *rParamsFile << "\t\t\t\t<LateralSurfaceTension>" << it->second[2] << "</LateralSurfaceTension>\n";
    }
    *rParamsFile << "\t\t\t<\nMappingMutationTension>\n";

    *rParamsFile << "\t\t\t<T1AreTemperatureDependent>" << mT1AreTemperatureDependent << "</T1AreTemperatureDependent>\n";
    *rParamsFile << "\t\t\t<T1TransitionRate>" << mT1TransitionRate << "</T1TransitionRate>\n";

    *rParamsFile << "\t\t\t<SimulatedAnnealingNoiseStrength>" << mRelativeNoiseStrength << "</SimulatedAnnealingNoiseStrength>\n";
    *rParamsFile << "\t\t\t<SimulatedAnnealingDecayTime>" << mAnnealingDecayTime << "</SimulatedAnnealingDecayTime>\n";
    *rParamsFile << "\t\t\t<SimulatedAnnealingInitialTemperature>" << mInitialTemperature << "</SimulatedAnnealingInitialTemperature>\n";

    *rParamsFile << "\t\t\t<VolumeElasticityParameter>" << mVolumeElasticityParameter << "</VolumeElasticityParameter>\n";
    *rParamsFile << "\t\t\t<ProjectOnVolumeConservingForces>" << mProjectOnVolumeConservingForces << "</ProjectOnVolumeConservingForces>\n";
    *rParamsFile << "\t\t\t<ActiveT1Swaps>" << mPerformActiveT1Swaps << "</ActiveT1Swaps>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class SurfaceTensionForce<1>;
template class SurfaceTensionForce<2>;
template class SurfaceTensionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"