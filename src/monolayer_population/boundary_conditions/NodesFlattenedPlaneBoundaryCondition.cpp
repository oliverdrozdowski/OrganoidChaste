#include "NodesFlattenedPlaneBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MonolayerVertexElement.hpp"
#include "MutableElement.hpp"
#include "VertexBasedCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NodesFlattenedPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::NodesFlattenedPlaneBoundaryCondition(
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
    c_vector<double, SPACE_DIM> midNormal)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(pCellPopulation),
          mUseJiggledNodesOnSubstrate(false)
{
    assert(norm_2(midNormal) > 0.0);
    mEdgeNormal = midNormal / norm_2(midNormal);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& NodesFlattenedPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetEdgeNormal() const
{
    return mEdgeNormal;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodesFlattenedPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::SetUseJiggledNodesOnSubstrate(bool useJiggledNodesOnSubstrate)
{
    mUseJiggledNodesOnSubstrate = useJiggledNodesOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool NodesFlattenedPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::GetUseJiggledNodesOnSubstrate()
{
    return mUseJiggledNodesOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodesFlattenedPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::SetEdgeNodes(std::vector<Node<SPACE_DIM>*> edgeNodes)
{
    mEdgeNodes = edgeNodes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodesFlattenedPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::ImposeBoundaryCondition(
    const std::map<Node<SPACE_DIM>*,
                   c_vector<double, SPACE_DIM> >& rOldLocations)
{
    ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("NodesFlattenedPlaneBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }

    assert((dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation))
           || (SPACE_DIM == ELEMENT_DIM && (dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation))));

    // This is a magic number
    double max_jiggle = 1e-4;

    if (SPACE_DIM != 1)
    {
        assert(SPACE_DIM == ELEMENT_DIM); // LCOV_EXCL_LINE
        assert(dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation));

        // First calculate center of mass of edge nodes
        c_vector<double, SPACE_DIM> position_com = zero_vector<double>(SPACE_DIM);
        unsigned num_nodes = 0;
        for (auto node_iter = mEdgeNodes.begin(); node_iter != mEdgeNodes.end(); ++node_iter)
        {
            Node<SPACE_DIM>* p_node = (*node_iter);
            position_com += p_node->rGetLocation();
            num_nodes++;
        }
        position_com /= num_nodes;

        // Then move nodes accordingly
        for (auto node_iter = mEdgeNodes.begin(); node_iter != mEdgeNodes.end(); ++node_iter)
        {
            Node<SPACE_DIM>* p_node = (*node_iter);
            c_vector<double, SPACE_DIM>& r_position = p_node->rGetModifiableLocation();

            double signed_distance = inner_prod(r_position - position_com, mEdgeNormal);
            c_vector<double, SPACE_DIM> distance_vector = signed_distance * mEdgeNormal;
            if (mUseJiggledNodesOnSubstrate)
            {
                distance_vector += max_jiggle * RandomNumberGenerator::Instance()->ranf() * mEdgeNormal;
            }
            r_position -= distance_vector;
        }
    }
    else
    {
        // DIM == 1
        NEVER_REACHED;
        // NodesFlattenedPlaneBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool NodesFlattenedPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::VerifyBoundaryCondition()
{
    return true;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void NodesFlattenedPlaneBoundaryCondition<1, 1>::ImposeBoundaryCondition(const std::map<Node<1>*, c_vector<double, 1> >& rOldLocations)
{
    EXCEPTION("NodesFlattenedPlaneBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void NodesFlattenedPlaneBoundaryCondition<1, 2>::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    EXCEPTION("NodesFlattenedPlaneBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void NodesFlattenedPlaneBoundaryCondition<1, 3>::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{
    EXCEPTION("NodesFlattenedPlaneBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void NodesFlattenedPlaneBoundaryCondition<2, 2>::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    EXCEPTION("NodesFlattenedPlaneBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void NodesFlattenedPlaneBoundaryCondition<2, 3>::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{
    EXCEPTION("NodesFlattenedPlaneBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool NodesFlattenedPlaneBoundaryCondition<1, 1>::VerifyBoundaryCondition()
{
    EXCEPTION("NodesFlattenedPlaneBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool NodesFlattenedPlaneBoundaryCondition<1, 2>::VerifyBoundaryCondition()
{
    EXCEPTION("NodesFlattenedPlaneBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool NodesFlattenedPlaneBoundaryCondition<1, 3>::VerifyBoundaryCondition()
{
    EXCEPTION("NodesFlattenedPlaneBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool NodesFlattenedPlaneBoundaryCondition<2, 2>::VerifyBoundaryCondition()
{
    EXCEPTION("NodesFlattenedPlaneBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool NodesFlattenedPlaneBoundaryCondition<2, 3>::VerifyBoundaryCondition()
{
    EXCEPTION("NodesFlattenedPlaneBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NodesFlattenedPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<EdgeNormal>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mEdgeNormal[index] << ",";
    }
    *rParamsFile << mEdgeNormal[SPACE_DIM - 1] << "</EdgeNormal>\n";

    *rParamsFile << "\t\t\t<UseJiggledNodesOnSubstrate>" << mUseJiggledNodesOnSubstrate << "</UseJiggledNodesOnSubstrate>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class NodesFlattenedPlaneBoundaryCondition<1, 1>;
template class NodesFlattenedPlaneBoundaryCondition<1, 2>;
template class NodesFlattenedPlaneBoundaryCondition<2, 2>;
template class NodesFlattenedPlaneBoundaryCondition<1, 3>;
template class NodesFlattenedPlaneBoundaryCondition<2, 3>;
template class NodesFlattenedPlaneBoundaryCondition<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NodesFlattenedPlaneBoundaryCondition)
