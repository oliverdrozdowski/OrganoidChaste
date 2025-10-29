#include "GeneralizedVolumeConservingForce.hpp"
#include <algorithm>
#include "LinearSystem.hpp"
#include "MonolayerVertexElement.hpp"
#include "PetscTools.hpp"
#include "PetscVecTools.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
GeneralizedVolumeConservingForce<DIM>::GeneralizedVolumeConservingForce()
        : AbstractForce<DIM>(),
          mSubForceCollection(),
          mProjectOnVolumeConservingForces(true),
          mIndividualNodeDampingActivated(false),
          mMapNodeToDampingThresholds(),
          mpSimulation(nullptr),
          mConstantLumenVolume(false)
{
}

template <unsigned DIM>
GeneralizedVolumeConservingForce<DIM>::~GeneralizedVolumeConservingForce()
{
}

template <unsigned DIM>
void GeneralizedVolumeConservingForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("GeneralizedVolumeConservingForce is to be used with a 3D MonolayerVertexBasedCellPopulation only");
    }
    if (mpSimulation == nullptr)
    {
        EXCEPTION("The pointer to the simulation instance must be set manually for GeneralizedVolumeConservingForce.");
    }

    // Define some helper variables
    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableMonolayerVertexMesh<DIM, DIM>& r_mesh = p_cell_population->rGetMesh();
    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();

    // The forces are determined from the subforces
    // Map from node pointers to force vectors
    NodeSubForceMap<DIM> forces_on_nodes;
    for (auto it_force = mSubForceCollection.begin(); it_force != mSubForceCollection.end(); ++it_force)
    {
        NodeSubForceMap<DIM> node_force_map = (*it_force)->GetForceContributions(*p_cell_population, mpSimulation);
        for (auto it_nd = node_force_map.begin(); it_nd != node_force_map.end(); ++it_nd)
        {
            // If node does not have a force, add the new force
            if (forces_on_nodes.find(it_nd->first) == forces_on_nodes.end())
            {
                forces_on_nodes[it_nd->first] = it_nd->second;
            }
            // otherwise just add to the vector
            else
            {
                forces_on_nodes[it_nd->first] += it_nd->second;
            }
        }
    }

    // We compute all the element volumes.
    // We also compute the volume gradients at each node
    std::vector<double> element_volumes(num_elements);
    // std::map< MonolayerVertexElement<DIM, DIM>*, double > element_volume_map;
    std::vector<double> target_volumes(num_elements);
    std::vector<double> cell_types(num_elements);

    // For volume conservation we determine a volume conservation matrix
    // This matrix contains the sum of scalar products of volume gradients
    // as entries and defines the linear system to determine the force projections
    unsigned linear_system_size = num_elements;

    if (mConstantLumenVolume)
    {
        linear_system_size += 1;
    }

    LinearSystem* p_system_volume_correction = new LinearSystem(linear_system_size, linear_system_size);
    LinearSystem* p_system_volume_conservation = new LinearSystem(linear_system_size, linear_system_size);
    // LinearSystem* p_system_random_force = new LinearSystem(num_elements, num_elements);

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
        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetVolumeModifier to the simulation in order to use a GeneralizedVolumeConservingForce");
        }

        // Add the excess volume to the RHS of system for volume correction
        double volume_excess = element_volumes[elem_index] - target_volumes[elem_index];
        p_system_volume_correction->AddToRhsVectorElement(elem_index, -volume_excess);
    }

    // Add excess Lumen Volume if wanted
    if (mConstantLumenVolume)
    {
        if (p_cell_population->GetTargetLumenVolume() == 0.0)
        {
            EXCEPTION("You need to set a target lumen volume.");
        }
        double volume_excess = r_mesh.GetLumenVolume() - p_cell_population->GetTargetLumenVolume();
        // std::cout << "Volume Excess: " << volume_excess << std::endl;
        p_system_volume_correction->AddToRhsVectorElement(num_elements, -volume_excess);
    }

    // Iterate over vertices in the cell population to determine volumes and volume gradients
    // Save all volume gradients, outer index nodes, inner index containing cells
    std::vector<std::vector<c_vector<double, DIM> > > vector_volume_gradients_all_nodes;

    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
        std::vector<c_vector<double, DIM> > vector_volume_gradients;

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        if (mConstantLumenVolume)
        {
            containing_elem_indices.insert(num_elements);
        }

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            if (*iter == num_elements) // this also makes sure that mConstantLumenVolume == true
            {
                c_vector<double, DIM> lumen_gradient = r_mesh.CalculateLumenVolGradient(node_index);
                vector_volume_gradients.push_back(lumen_gradient);
            }
            else
            {
                // Get this element, its index and its number of nodes
                MonolayerVertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
                unsigned elem_index = p_element->GetIndex(); // same as *iter

                unsigned local_node_index = p_element->GetNodeLocalIndex(node_index);

                // We then determine the volume gradient, which we need to
                // calculate the volume conservation matrix
                c_vector<double, DIM> volume_gradient = r_mesh.GetVolumeGradientAtNode(elem_index, local_node_index);
                vector_volume_gradients.push_back(volume_gradient);
            }
        }

        vector_volume_gradients_all_nodes.push_back(vector_volume_gradients);

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
            unsigned global_elem_index_larger;
            if (*iter_elem_larger == num_elements) // this also makes sure that mConstantLumenVolume == true
            {
                global_elem_index_larger = num_elements; // index for lumen
            }
            else
            {
                MonolayerVertexElement<DIM, DIM>* p_element_larger = p_cell_population->GetElement(*iter_elem_larger);
                global_elem_index_larger = p_element_larger->GetIndex();
            }

            c_vector<double, DIM> volume_gradient_1 = vector_volume_gradients[index_volume_gradient_larger];

            unsigned index_volume_gradient_smaller = 0;
            std::set<unsigned>::iterator iter_elem_larger_next = iter_elem_larger;
            iter_elem_larger_next++;
            for (std::set<unsigned>::iterator iter_elem_smaller = containing_elem_indices.begin();
                 iter_elem_smaller != iter_elem_larger_next;
                 ++iter_elem_smaller)
            {
                unsigned global_elem_index_smaller;
                if (*iter_elem_smaller == num_elements) // this also makes sure that mConstantLumenVolume == true
                {
                    global_elem_index_smaller = num_elements; // index for lumen
                }
                else
                {
                    MonolayerVertexElement<DIM, DIM>* p_element_smaller = p_cell_population->GetElement(*iter_elem_smaller);
                    global_elem_index_smaller = p_element_smaller->GetIndex();
                }

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
                // p_system_random_force->AddToMatrixElement(global_elem_index_larger, global_elem_index_smaller, product_volume_gradients);
                if (global_elem_index_smaller != global_elem_index_larger)
                {
                    p_system_volume_conservation->AddToMatrixElement(global_elem_index_smaller, global_elem_index_larger, product_volume_gradients);
                    // p_system_random_force->AddToMatrixElement(global_elem_index_smaller, global_elem_index_larger, product_volume_gradients);
                }

                index_volume_gradient_smaller++;
            }

            // add force to rhs
            c_vector<double, DIM> force_on_node = forces_on_nodes[p_this_node];
            double product_force_vol_gradient = inner_prod(force_on_node, volume_gradient_1);
            p_system_volume_conservation->AddToRhsVectorElement(global_elem_index_larger, product_force_vol_gradient);
            // add random force to rhs
            // double product_random_force_vol_gradient = inner_prod(random_force_on_node, volume_gradient_1);
            // p_system_random_force->AddToRhsVectorElement(global_elem_index_larger, product_random_force_vol_gradient);

            // Flush entries
            p_system_volume_correction->AssembleIntermediateLinearSystem();
            p_system_volume_conservation->AssembleIntermediateLinearSystem();
            // p_system_random_force->AssembleIntermediateLinearSystem();

            index_volume_gradient_larger++;
        }
    }

    // Now that all forces and gradients have been calculated and inserted into the linear systems
    // we can solve for the parameters, which determine the volume correction (in the direction of
    // volume gradients) and the force projected on the conservation space (orthogonal to vol. grads.)
    p_system_volume_correction->AssembleFinalLinearSystem();
    p_system_volume_conservation->AssembleFinalLinearSystem();
    // p_system_random_force->AssembleFinalLinearSystem();

    Vec parameters_correction = nullptr;
    Vec parameters_conservation = nullptr;
    // Vec parameters_random_force = nullptr;

    parameters_correction = p_system_volume_correction->Solve();
    parameters_conservation = p_system_volume_conservation->Solve();
    // parameters_random_force;
    /*
    if(mRelativeNoiseStrength > 0.0)
    {
      parameters_random_force = p_system_random_force->Solve();
    }
    */

    // Iterate over vertices in the cell population to calculate and save the force contributions
    // and perform the single-step volume correction
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        c_vector<double, DIM> correction_vector = zero_vector<double>(DIM);
        c_vector<double, DIM> force_on_node = forces_on_nodes[p_cell_population->GetNode(node_index)];
        // c_vector<double, DIM> random_force_on_node = random_forces[p_cell_population->GetNode(node_index)];

        // We perform the volume correction and the force projection in one step
        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        if (mConstantLumenVolume)
        {
            containing_elem_indices.insert(num_elements);
        }

        // Iterate over these elements
        unsigned local_cell_index_for_node = 0;
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            unsigned global_elem_index;
            if (*iter == num_elements)
            {
                global_elem_index = num_elements;
            }
            else
            {
                MonolayerVertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
                global_elem_index = p_element->GetIndex();
            }

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
            // old_node_locations[p_cell_population->GetNode(node_index)] = node_position_vector;

            // Determine cell-indexed vectors for projection of forces
            double coefficient_forces = PetscVecTools::GetElement(parameters_conservation, global_elem_index);
            /*
            double coefficient_random_forces=0.0;
            if(mRelativeNoiseStrength > 0.0)
            {
              coefficient_random_forces = PetscVecTools::GetElement(parameters_random_force, global_elem_index);
            }
            */
            // Project the force
            force_on_node -= coefficient_forces * gradient_vector_cell_node;
            // random_force_on_node -= coefficient_random_forces * gradient_vector_cell_node;

            local_cell_index_for_node++;
        }
        // Save the force into the node
        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);
        // random_forces[p_cell_population->GetNode(node_index)] = random_force_on_node;
    }

    /* For simulated annealing we now have to calculate the energy difference we get from the
     * random force. If this energy difference is negative, we accept the step, if not we
     * accept with the Boltzmann probability with the current annealing temperature.
     *
     * For this we first calculate the energy, move the nodes according to our random
     * forces, calculate the new energy, move the nodes back and finally save the
     * force contribution if this step is accepted
     */
    /*
    if((mRelativeNoiseStrength > 0.0) & (mpSimulation != nullptr))
    {
      double energy_before_random = CalculateTotalEnergy(rCellPopulation);
      double timestep = mpSimulation->GetDt();
      // Move the nodes (assuming Forward Euler Step)
      for (unsigned node_index=0; node_index<num_nodes; node_index++)
      {
        Node<DIM>* p_node = p_cell_population->GetNode(node_index);
        double damping = p_cell_population -> GetDampingConstant(node_index);
        c_vector<double,DIM> random_force_on_node = random_forces[p_node];
        c_vector<double, DIM>& location_node = p_node->rGetModifiableLocation();
        location_node += random_force_on_node * timestep / damping;
      }
      // Compute the new energy
      double energy_after_random = CalculateTotalEnergy(rCellPopulation);

      // Reset the nodes
      for (unsigned node_index=0; node_index<num_nodes; node_index++)
      {
        Node<DIM>* p_node = p_cell_population->GetNode(node_index);
        c_vector<double, DIM>& location_node = p_node->rGetModifiableLocation();
        location_node = old_node_locations[p_node];
      }
      // Add the force contribution, depending on energy
      double uniform_rand = p_random->ranf();
      double boltzmann_factor = exp(-(energy_after_random - energy_before_random)/temperature);
      if(uniform_rand < boltzmann_factor)
      {
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
          Node<DIM>* p_node = p_cell_population->GetNode(node_index);
          p_node -> AddAppliedForceContribution(random_forces[p_node]);
        }
      }
    }
    */
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
    // delete p_system_random_force;

    // Must release the memory from the Petsc Vector to not have memory leaks
    if (parameters_correction)
    {
        PetscTools::Destroy(parameters_correction);
    }
    if (parameters_conservation)
    {
        PetscTools::Destroy(parameters_conservation);
    }
    // if (parameters_random_force)
    //{
    //     PetscTools::Destroy(parameters_random_force);
    // }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void GeneralizedVolumeConservingForce<1>::AddForceContribution(AbstractCellPopulation<1>& rCellPopulation)
{
    EXCEPTION("GeneralizedVolumeConservingForce not possible in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
bool GeneralizedVolumeConservingForce<DIM>::DoProjectOnVolumeConservingForces()
{
    return true;
}

template <unsigned DIM>
void GeneralizedVolumeConservingForce<DIM>::SetProjectOnVolumeConservingForces(bool doProjection)
{
    if (!doProjection)
    {
        EXCEPTION("We have not yet implement GeneralizedVolumeConservingForce without Volume conservation projection!");
    }
}

template <unsigned DIM>
void GeneralizedVolumeConservingForce<DIM>::SetSimulationInstance(AbstractCellBasedSimulation<DIM, DIM>* pSimulation)
{
    mpSimulation = pSimulation;
}

template <unsigned DIM>
AbstractCellBasedSimulation<DIM, DIM>* GeneralizedVolumeConservingForce<DIM>::GetSimulationInstance()
{
    return mpSimulation;
}

template <unsigned DIM>
void GeneralizedVolumeConservingForce<DIM>::SetMapNodeToDampingThresholds(
    std::map<Node<DIM>*, DampingCoefficientPair> mapNodeToDampingThreshold)
{
    mMapNodeToDampingThresholds = mapNodeToDampingThreshold;
    mIndividualNodeDampingActivated = true;
}

template <unsigned DIM>
void GeneralizedVolumeConservingForce<DIM>::AddSubForce(boost::shared_ptr<AbstractSubForce<DIM> > pForce)
{
    mSubForceCollection.push_back(pForce);
}

template <unsigned DIM>
void GeneralizedVolumeConservingForce<DIM>::RemoveSubForce(boost::shared_ptr<AbstractSubForce<DIM> > pForce)
{
    auto it_del = std::find(mSubForceCollection.begin(), mSubForceCollection.end(), pForce);
    if (it_del != mSubForceCollection.end())
    {
        mSubForceCollection.erase(it_del);
    }
}

template <unsigned DIM>
boost::shared_ptr<SurfaceTensionSubForce<DIM> > GeneralizedVolumeConservingForce<DIM>::GetSurfaceTensionSubForce()
{
    for (auto it_force = mSubForceCollection.begin(); it_force != mSubForceCollection.end(); ++it_force)
    {
        if (dynamic_cast<SurfaceTensionSubForce<DIM>*>(&(*(*it_force))) != nullptr)
            return boost::dynamic_pointer_cast<SurfaceTensionSubForce<DIM> >(*it_force);
    }
    return nullptr;
}

template <unsigned DIM>
void GeneralizedVolumeConservingForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ProjectOnVolumeConservingForces>" << mProjectOnVolumeConservingForces << "</ProjectOnVolumeConservingForces>\n";
    *rParamsFile << "\t\t\t<IndividualNodeDampingActivated>" << mIndividualNodeDampingActivated << "</IndividualNodeDampingActivated>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

template <unsigned DIM>
void GeneralizedVolumeConservingForce<DIM>::SetConstantLumenVolumeFlag(bool keepLumenVolumeConstant)
{
    mConstantLumenVolume = keepLumenVolumeConstant;
}

// Explicit instantiation
template class GeneralizedVolumeConservingForce<1>;
template class GeneralizedVolumeConservingForce<2>;
template class GeneralizedVolumeConservingForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"