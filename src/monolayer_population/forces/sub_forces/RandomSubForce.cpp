#include "RandomSubForce.hpp"
#include "LinearSystem.hpp"
#include "MonolayerVertexElement.hpp"
#include "PetscVecTools.hpp"
#include "RandomNumberGenerator.hpp"
#include "Temperature.hpp"

template <unsigned DIM>
RandomSubForce<DIM>::RandomSubForce(double strength, boost::shared_ptr<Temperature> temp) : AbstractSubForce<DIM>(), mStrength(strength), mTemperature(temp)
{
}

template <unsigned DIM>
RandomSubForce<DIM>::~RandomSubForce()
{
}

template <unsigned DIM>
NodeSubForceMap<DIM> RandomSubForce<DIM>::GetForceContributions(
    MonolayerVertexBasedCellPopulation<DIM>& rCellPopulation,
    AbstractCellBasedSimulation<DIM, DIM>* pSimulation)
{
    // Define some helper variables
    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_nodes = p_cell_population->GetNumNodes();

    // Create pointer to random number generator
    RandomNumberGenerator* p_random = RandomNumberGenerator::Instance();

    // Calculate current strength
    double current_strength = mStrength * mTemperature->GetCurrentTemperature();

    // Map of nodes to forces which will be returned
    NodeSubForceMap<DIM> map_forces_on_nodes;

    // Iterate over vertices in the cell population
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);

        c_vector<double, DIM> random_contribution = zero_vector<double>(DIM);

        for (unsigned i_dim = 0; i_dim < DIM; i_dim++)
        {
            // Add Gaussian noise
            double std_rand = p_random->StandardNormalRandomDeviate();

            random_contribution[i_dim] += current_strength * std_rand;
        }

        for (auto boundary : p_cell_population->GetMovingBoundaries())
        {
            double boundary_z_position = boundary->GetZPosition();

            double z_distance = boundary_z_position - p_cell_population->GetNode(node_index)->rGetLocation()[DIM - 1];

            if (abs(z_distance < 0.5) & ((boundary_z_position * random_contribution[DIM - 1]) > 0))
            {
                random_contribution[DIM - 1] = 0.0;
            }
        }

        map_forces_on_nodes[p_this_node] = random_contribution;
    }
    return map_forces_on_nodes;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
NodeSubForceMap<1> RandomSubForce<1>::GetForceContributions(
    MonolayerVertexBasedCellPopulation<1>& rCellPopulation,
    AbstractCellBasedSimulation<1, 1>* pSimulation)
{
    EXCEPTION("RandomSubForce not possible in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
void RandomSubForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<RandomStrengthParameter>"
                 << mStrength
                 << "</RandomStrengthParameter>\n";
}

template <unsigned DIM>
void RandomSubForce<DIM>::ResetTemperature()
{
    mTemperature->SetStartTimeToNow();
}

template <unsigned DIM>
void RandomSubForce<DIM>::SetStrength(double strength)
{
    mStrength = strength;
}

template <unsigned DIM>
void RandomSubForce<DIM>::SetTemperatureDecay(double decay_time)
{
    mTemperature->SetRelativeDecayTime(decay_time);
}

template <unsigned DIM>
void RandomSubForce<DIM>::SetTemperatureToConstant()
{
    mTemperature->SetTemperatureToConstant();
}

// Explicit instantiation
template class RandomSubForce<1>;
template class RandomSubForce<2>;
template class RandomSubForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"