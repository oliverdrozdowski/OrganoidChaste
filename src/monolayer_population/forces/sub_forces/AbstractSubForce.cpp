#include "AbstractSubForce.hpp"

template <unsigned DIM>
AbstractSubForce<DIM>::AbstractSubForce()
{
}

template <unsigned DIM>
AbstractSubForce<DIM>::~AbstractSubForce()
{
}

template <unsigned DIM>
void AbstractSubForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("Standard implementation of AbstractSubForce only allows for MonolayerVertexBasedCellPopulations.");
    }

    // Get forces from GetForceContributions and write force to nodes.
    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    NodeSubForceMap<DIM> node_force_map = GetForceContributions((*p_cell_population), mpSimulation);
    for (auto it = node_force_map.begin(); it != node_force_map.end(); ++it)
    {
        it->first->AddAppliedForceContribution(it->second);
    }
}

template <unsigned DIM>
void AbstractSubForce<DIM>::SetSimulationInstance(AbstractCellBasedSimulation<DIM, DIM>* pSimulation)
{
    mpSimulation = pSimulation;
}

template <unsigned DIM>
AbstractCellBasedSimulation<DIM, DIM>* AbstractSubForce<DIM>::GetSimulationInstance()
{
    return mpSimulation;
}

template <unsigned DIM>
void AbstractSubForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to output
}

// Explicit instantiation
template class AbstractSubForce<1>;
template class AbstractSubForce<2>;
template class AbstractSubForce<3>;
