#include "LumenPressureSubForce.hpp"
#include "LinearSystem.hpp"
#include "MonolayerVertexElement.hpp"
#include "PetscVecTools.hpp"

template <unsigned DIM>
LumenPressureSubForce<DIM>::LumenPressureSubForce(double pressure_constant)
        : AbstractSubForce<DIM>(), mPressureConstant(pressure_constant)
{
}

template <unsigned DIM>
LumenPressureSubForce<DIM>::~LumenPressureSubForce()
{
}

template <unsigned DIM>
NodeSubForceMap<DIM> LumenPressureSubForce<DIM>::GetForceContributions(
    MonolayerVertexBasedCellPopulation<DIM>& rCellPopulation,
    AbstractCellBasedSimulation<DIM, DIM>* pSimulation)
{
    // Define some helper variables
    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_nodes = p_cell_population->GetNumNodes();

    MutableMonolayerVertexMesh<DIM, DIM>& r_mesh = p_cell_population->rGetMesh();

    // Calculated current pressure
    // double pressure = mPressureConstant / r_mesh.GetLumenVolume();
    double pressure = mPressureConstant;

    // Map of nodes to forces which will be returned
    NodeSubForceMap<DIM> map_forces_on_nodes;

    // Iterate over vertices in the cell population
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
        map_forces_on_nodes[p_this_node] = pressure * r_mesh.CalculateLumenVolGradient(node_index);
    }
    return map_forces_on_nodes;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
NodeSubForceMap<1> LumenPressureSubForce<1>::GetForceContributions(
    MonolayerVertexBasedCellPopulation<1>& rCellPopulation,
    AbstractCellBasedSimulation<1, 1>* pSimulation)
{
    EXCEPTION("LumenPressureSubForce not possible in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
void LumenPressureSubForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PressureParameter>"
                 << mPressureConstant
                 << "</PressureParameter>\n";
}

template <unsigned DIM>
void LumenPressureSubForce<DIM>::SetPressureConstant(double pressure_constant)
{
    mPressureConstant = pressure_constant;
}

template <unsigned DIM>
double LumenPressureSubForce<DIM>::GetPressureConstant()
{
    return mPressureConstant;
}

// Explicit instantiation
template class LumenPressureSubForce<1>;
template class LumenPressureSubForce<2>;
template class LumenPressureSubForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"