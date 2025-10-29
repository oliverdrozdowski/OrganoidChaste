#include "MovingBoundarySubForce.hpp"
#include "LinearSystem.hpp"
#include "MonolayerVertexElement.hpp"
#include "PetscVecTools.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
MovingBoundarySubForce<DIM>::MovingBoundarySubForce()
        : AbstractSubForce<DIM>()
{
}

template <unsigned DIM>
MovingBoundarySubForce<DIM>::~MovingBoundarySubForce()
{
}

template <unsigned DIM>
NodeSubForceMap<DIM> MovingBoundarySubForce<DIM>::GetForceContributions(
    MonolayerVertexBasedCellPopulation<DIM>& rCellPopulation,
    AbstractCellBasedSimulation<DIM, DIM>* pSimulation)
{
    // Define some helper variables
    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_nodes = p_cell_population->GetNumNodes();

    // clear for read out
    for (auto boundary : p_cell_population->GetMovingBoundaries())
    {
        boundary->ResetAppliedForceOnPopulation();
    }

    // Map of nodes to forces which will be returned
    NodeSubForceMap<DIM> map_forces_on_nodes;

    // Iterate over vertices in the cell population
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);

        c_vector<double, DIM> boundary_contribution = zero_vector<double>(DIM);

        for (auto boundary : p_cell_population->GetMovingBoundaries())
        {
            c_vector<double, DIM> boundary_force = boundary->CalculateForceOnNode(p_this_node->rGetLocation(), p_cell_population->GetDampingConstantNormal(), pSimulation->GetDt());

            // needed for read out
            boundary->AddAppliedForceOnPopulation(boundary_force);

            boundary_contribution += boundary_force;
        }

        map_forces_on_nodes[p_this_node] = boundary_contribution;
    }
    return map_forces_on_nodes;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
NodeSubForceMap<1> MovingBoundarySubForce<1>::GetForceContributions(
    MonolayerVertexBasedCellPopulation<1>& rCellPopulation,
    AbstractCellBasedSimulation<1, 1>* pSimulation)
{
    EXCEPTION("MovingBoundarySubForce not possible in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
void MovingBoundarySubForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
}

// Explicit instantiation
template class MovingBoundarySubForce<1>;
template class MovingBoundarySubForce<2>;
template class MovingBoundarySubForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"