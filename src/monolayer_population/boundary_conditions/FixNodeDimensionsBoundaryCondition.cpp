#include "FixNodeDimensionsBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MonolayerVertexElement.hpp"
#include "MutableElement.hpp"
#include "VertexBasedCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FixNodeDimensionsBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::FixNodeDimensionsBoundaryCondition(
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
    std::vector<Node<SPACE_DIM>*> nodesToConsider,
    bool keepConstantX,
    bool keepConstantY,
    bool keepConstantZ)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(pCellPopulation),
          mNodesToConsider(nodesToConsider),
          mKeepConstantX(keepConstantX),
          mKeepConstantY(keepConstantY),
          mKeepConstantZ(keepConstantZ)
{
    // First we save the node positions for later
    for (auto it = nodesToConsider.begin(); it != nodesToConsider.end(); ++it)
    {
        mOriginalNodePositions[(*it)] = (*it)->rGetLocation();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FixNodeDimensionsBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::ImposeBoundaryCondition(
    const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations)
{
    if (SPACE_DIM == 3)
    {
        for (auto it = mNodesToConsider.begin(); it != mNodesToConsider.end(); ++it)
        {
            if (mKeepConstantX)
            {
                (*it)->rGetModifiableLocation()[0] = mOriginalNodePositions[(*it)][0];
            }
            if (mKeepConstantY)
            {
                (*it)->rGetModifiableLocation()[1] = mOriginalNodePositions[(*it)][1];
            }
            if (mKeepConstantZ)
            {
                (*it)->rGetModifiableLocation()[2] = mOriginalNodePositions[(*it)][2];
            }
        }
    }
    else
    {
        EXCEPTION("FixNodeDimensionsBoundaryCondition only implemented in 3D!");
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool FixNodeDimensionsBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    if (SPACE_DIM != 3)
    {
        EXCEPTION("FixNodeDimensionsBoundaryCondition is only implemented in 3D");
    }
    else
    {
        return condition_satisfied;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FixNodeDimensionsBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<KeepConstantX>";
    *rParamsFile << mKeepConstantX;
    *rParamsFile << "</KeepConstantX>\n";

    *rParamsFile << "\t\t\t<KeepConstantY>";
    *rParamsFile << mKeepConstantY;
    *rParamsFile << "</KeepConstantY>\n";

    *rParamsFile << "\t\t\t<KeepConstantZ>";
    *rParamsFile << mKeepConstantZ;
    *rParamsFile << "</KeepConstantZ>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class FixNodeDimensionsBoundaryCondition<1, 1>;
template class FixNodeDimensionsBoundaryCondition<1, 2>;
template class FixNodeDimensionsBoundaryCondition<2, 2>;
template class FixNodeDimensionsBoundaryCondition<1, 3>;
template class FixNodeDimensionsBoundaryCondition<2, 3>;
template class FixNodeDimensionsBoundaryCondition<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FixNodeDimensionsBoundaryCondition)
