#include "StraightEdgeBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MonolayerVertexElement.hpp"
#include "MutableElement.hpp"
#include "VertexBasedCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
StraightEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::StraightEdgeBoundaryCondition(
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
    c_vector<double, SPACE_DIM> midNormal)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(pCellPopulation),
          mUseJiggledNodesOnSubstrate(false)
{
    assert(norm_2(midNormal) > 0.0);
    mEdgeNormal = midNormal / norm_2(midNormal);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& StraightEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetEdgeNormal() const
{
    return mEdgeNormal;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StraightEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::SetUseJiggledNodesOnSubstrate(bool useJiggledNodesOnSubstrate)
{
    mUseJiggledNodesOnSubstrate = useJiggledNodesOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool StraightEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::GetUseJiggledNodesOnSubstrate()
{
    return mUseJiggledNodesOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StraightEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::SetEdgeNodes(std::vector<Node<SPACE_DIM>*> edgeNodes)
{
    mEdgeNodes = edgeNodes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StraightEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::ImposeBoundaryCondition(
    const std::map<Node<SPACE_DIM>*,
                   c_vector<double, SPACE_DIM> >& rOldLocations)
{
    ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("StraightEdgeBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }

    assert((dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation))
           || (SPACE_DIM == ELEMENT_DIM && (dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation))));

    // This is a magic number
    double max_jiggle = 1e-4;

    if (SPACE_DIM != 1)
    {
        assert(SPACE_DIM == ELEMENT_DIM); // LCOV_EXCL_LINE
        assert(dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation));
        MonolayerVertexBasedCellPopulation<SPACE_DIM>* p_cell_population = dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation);

        MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& r_mesh = p_cell_population->rGetMesh();

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

        // Loop over lateral faces
        for (auto face_iter = r_mesh.GetFaceIteratorBegin(); face_iter != r_mesh.GetFaceIteratorEnd(); ++face_iter)
        {
            // Only check lateral faces, which are boundary faces
            // Note that the mesh MUST have properly set boundary faces!
            if (face_iter->GetFaceType() != MonolayerVertexElementType::Lateral || !face_iter->IsBoundaryFace())
            {
                continue;
            }

            // We now go along the face to find an apico-basal edge
            unsigned num_nodes = face_iter->GetNumNodes();
            for (unsigned local_index = 0; local_index < num_nodes; local_index++)
            {
                Node<SPACE_DIM>* p_current_node = face_iter->GetNode(local_index);
                Node<SPACE_DIM>* p_next_node = face_iter->GetNode((local_index + 1) % num_nodes);

                MonolayerVertexElementType current_node_type = face_iter->GetNodeType(local_index);
                MonolayerVertexElementType next_node_type = face_iter->GetNodeType((local_index + 1) % num_nodes);

                Node<SPACE_DIM>* p_basal_node = nullptr;
                Node<SPACE_DIM>* p_apical_node = nullptr;
                // Name nodes correctly
                if (current_node_type == MonolayerVertexElementType::Basal && next_node_type == MonolayerVertexElementType::Apical)
                {
                    p_basal_node = p_current_node;
                    p_apical_node = p_next_node;
                }
                else if (current_node_type == MonolayerVertexElementType::Apical && next_node_type == MonolayerVertexElementType::Basal)
                {
                    p_apical_node = p_current_node;
                    p_basal_node = p_next_node;
                }
                else
                {
                    continue;
                }

                // Check if both nodes should actually be constrained
                // If no list of nodes was supplied, we just assume alle edges are constrained
                // If first element is nullptr, we assume no nodes are constrained
                if (!mEdgeNodes.empty())
                {
                    if (mEdgeNodes[0] == nullptr)
                    {
                        continue;
                    }
                    if (std::find(mEdgeNodes.begin(), mEdgeNodes.end(), p_apical_node) == mEdgeNodes.end())
                    {
                        continue;
                    }
                    if (std::find(mEdgeNodes.begin(), mEdgeNodes.end(), p_basal_node) == mEdgeNodes.end())
                    {
                        continue;
                    }
                }

                // Determine distance from edge center plane
                c_vector<double, SPACE_DIM>& basal_point = p_basal_node->rGetModifiableLocation();
                c_vector<double, SPACE_DIM>& apical_point = p_apical_node->rGetModifiableLocation();
                double signed_distance = inner_prod((apical_point + basal_point) / 2.0 - position_com, mEdgeNormal);
                c_vector<double, SPACE_DIM> distance_vector = signed_distance * mEdgeNormal;
                if (mUseJiggledNodesOnSubstrate)
                {
                    distance_vector += max_jiggle * RandomNumberGenerator::Instance()->ranf() * mEdgeNormal;
                }

                apical_point -= distance_vector;
                basal_point -= distance_vector;
            }
        }
    }
    else
    {
        // DIM == 1
        NEVER_REACHED;
        // StraightEdgeBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool StraightEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::VerifyBoundaryCondition()
{
    return true;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void StraightEdgeBoundaryCondition<1, 1>::ImposeBoundaryCondition(const std::map<Node<1>*, c_vector<double, 1> >& rOldLocations)
{
    EXCEPTION("StraightEdgeBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void StraightEdgeBoundaryCondition<1, 2>::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    EXCEPTION("StraightEdgeBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void StraightEdgeBoundaryCondition<1, 3>::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{
    EXCEPTION("StraightEdgeBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void StraightEdgeBoundaryCondition<2, 2>::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    EXCEPTION("StraightEdgeBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void StraightEdgeBoundaryCondition<2, 3>::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{
    EXCEPTION("StraightEdgeBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool StraightEdgeBoundaryCondition<1, 1>::VerifyBoundaryCondition()
{
    EXCEPTION("StraightEdgeBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool StraightEdgeBoundaryCondition<1, 2>::VerifyBoundaryCondition()
{
    EXCEPTION("StraightEdgeBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool StraightEdgeBoundaryCondition<1, 3>::VerifyBoundaryCondition()
{
    EXCEPTION("StraightEdgeBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool StraightEdgeBoundaryCondition<2, 2>::VerifyBoundaryCondition()
{
    EXCEPTION("StraightEdgeBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool StraightEdgeBoundaryCondition<2, 3>::VerifyBoundaryCondition()
{
    EXCEPTION("StraightEdgeBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void StraightEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
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
template class StraightEdgeBoundaryCondition<1, 1>;
template class StraightEdgeBoundaryCondition<1, 2>;
template class StraightEdgeBoundaryCondition<2, 2>;
template class StraightEdgeBoundaryCondition<1, 3>;
template class StraightEdgeBoundaryCondition<2, 3>;
template class StraightEdgeBoundaryCondition<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(StraightEdgeBoundaryCondition)
