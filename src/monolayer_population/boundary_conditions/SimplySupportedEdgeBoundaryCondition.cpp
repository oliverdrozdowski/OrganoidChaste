#include "SimplySupportedEdgeBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MonolayerVertexElement.hpp"
#include "MutableElement.hpp"
#include "VertexBasedCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SimplySupportedEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::SimplySupportedEdgeBoundaryCondition(
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
    c_vector<double, SPACE_DIM> midPoint,
    c_vector<double, SPACE_DIM> midNormal)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(pCellPopulation),
          mMidPlanePoint(midPoint),
          mUseJiggledNodesOnSubstrate(false)
{
    assert(norm_2(midNormal) > 0.0);
    mMidPlaneNormal = midNormal / norm_2(midNormal);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& SimplySupportedEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetMidPlanePoint() const
{
    return mMidPlanePoint;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& SimplySupportedEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetMidPlaneNormal() const
{
    return mMidPlaneNormal;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SimplySupportedEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::SetUseJiggledNodesOnSubstrate(bool useJiggledNodesOnSubstrate)
{
    mUseJiggledNodesOnSubstrate = useJiggledNodesOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool SimplySupportedEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::GetUseJiggledNodesOnSubstrate()
{
    return mUseJiggledNodesOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SimplySupportedEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::SetMidPlaneNodes(std::vector<Node<SPACE_DIM>*> midPlaneNodes)
{
    mMidPlaneNodes = midPlaneNodes;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SimplySupportedEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::ImposeBoundaryCondition(
    const std::map<Node<SPACE_DIM>*,
                   c_vector<double, SPACE_DIM> >& rOldLocations)
{
    ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("SimplySupportedEdgeBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
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

        // Loop over lateral faces
        // First we orthogonalize the edges. In a second step we correct the midplane.
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

                // Determine midplane distance
                c_vector<double, SPACE_DIM>& basal_point = p_basal_node->rGetModifiableLocation();
                c_vector<double, SPACE_DIM>& apical_point = p_apical_node->rGetModifiableLocation();

                // Determine average basal & apical normals for elements that contain the edge
                // If edges are shared by elements, we do this procedure multiple times; no harm done
                std::set<unsigned> elements_of_basal = p_basal_node->rGetContainingElementIndices();
                c_vector<double, SPACE_DIM> face_normal = zero_vector<double>(SPACE_DIM);
                for (auto it_elt = elements_of_basal.begin(); it_elt != elements_of_basal.end(); ++it_elt)
                {
                    c_vector<double, SPACE_DIM> nrml = zero_vector<double>(SPACE_DIM);

                    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = r_mesh.GetFaceOfType(*it_elt, MonolayerVertexElementType::Apical);
                    r_mesh.CalculateUnitNormalToFaceWithArea(p_face, nrml);
                    nrml *= inner_prod(face_normal, nrml) < 0.0 ? -1.0 : 1.0;
                    face_normal += nrml;

                    nrml = zero_vector<double>(SPACE_DIM);
                    p_face = r_mesh.GetFaceOfType(*it_elt, MonolayerVertexElementType::Basal);
                    r_mesh.CalculateUnitNormalToFaceWithArea(p_face, nrml);
                    nrml *= inner_prod(face_normal, nrml) < 0.0 ? -1.0 : 1.0;
                    face_normal += nrml;
                }
                face_normal /= norm_2(face_normal);

                // Find orthogonal component to averaged normal, align edge and shift to midplane
                double parallel_component = inner_prod((apical_point - basal_point), face_normal);
                c_vector<double, 3> orthogonal_vector = (apical_point - basal_point) - parallel_component * face_normal;
                apical_point -= orthogonal_vector / 2.0;
                basal_point += orthogonal_vector / 2.0;
            }
        }
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
                if (!mMidPlaneNodes.empty())
                {
                    if (mMidPlaneNodes[0] == nullptr)
                    {
                        continue;
                    }
                    if (std::find(mMidPlaneNodes.begin(), mMidPlaneNodes.end(), p_apical_node) == mMidPlaneNodes.end())
                    {
                        continue;
                    }
                    if (std::find(mMidPlaneNodes.begin(), mMidPlaneNodes.end(), p_basal_node) == mMidPlaneNodes.end())
                    {
                        continue;
                    }
                }

                // Determine midplane distance
                c_vector<double, SPACE_DIM>& basal_point = p_basal_node->rGetModifiableLocation();
                c_vector<double, SPACE_DIM>& apical_point = p_apical_node->rGetModifiableLocation();
                double signed_distance = inner_prod((apical_point + basal_point) / 2.0 - mMidPlanePoint, mMidPlaneNormal);
                c_vector<double, SPACE_DIM> distance_vector = signed_distance * mMidPlaneNormal;
                if (mUseJiggledNodesOnSubstrate)
                {
                    distance_vector += max_jiggle * RandomNumberGenerator::Instance()->ranf() * mMidPlaneNormal;
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
        // SimplySupportedEdgeBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool SimplySupportedEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::VerifyBoundaryCondition()
{
    // bool condition_satisfied = true;
    return true;
    /*****
    if (SPACE_DIM == 1)
    {
        EXCEPTION("SimplySupportedEdgeBoundaryCondition is not implemented in 1D");
    }
    else
    {
        assert(dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation));
        MonolayerVertexBasedCellPopulation<SPACE_DIM>* p_cell_population = dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation);

        // This is a magic number
        double max_jiggle = 1e-4;

        MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& r_mesh = p_cell_population->rGetMesh();

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

                // Determine midplane distance
                c_vector<double, SPACE_DIM>& basal_point = p_basal_node->rGetModifiableLocation();
                c_vector<double, SPACE_DIM>& apical_point = p_apical_node->rGetModifiableLocation();
                double signed_distance = inner_prod((apical_point + basal_point) / 2.0 - mMidPlanePoint, mMidPlaneNormal);

                if (abs(signed_distance) > max_jiggle)
                {
                    condition_satisfied = false;
                }

                // Determine average basal & apical normals for elements that contain the edge
                // If edges are shared by elements, we do this procedure multiple times; no harm done
                std::set<unsigned> elements_of_basal = p_basal_node->rGetContainingElementIndices();
                c_vector<double, SPACE_DIM> face_normal = zero_vector<double>(3);
                for (auto it_elt = elements_of_basal.begin(); it_elt != elements_of_basal.end(); ++it_elt)
                {
                    c_vector<double, SPACE_DIM> nrml = zero_vector<double>(SPACE_DIM);

                    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = r_mesh.GetFaceOfType(*it_elt, MonolayerVertexElementType::Apical);
                    r_mesh.CalculateUnitNormalToFaceWithArea(p_face, nrml);
                    nrml *= inner_prod(face_normal, nrml) < 0.0 ? -1.0 : 1.0;
                    face_normal += nrml;

                    nrml = zero_vector<double>(3);
                    p_face = r_mesh.GetFaceOfType(*it_elt, MonolayerVertexElementType::Basal);
                    r_mesh.CalculateUnitNormalToFaceWithArea(p_face, nrml);
                    nrml *= inner_prod(face_normal, nrml) < 0.0 ? -1.0 : 1.0;
                    face_normal += nrml;
                }
                face_normal /= norm_2(face_normal);

                // Find orthogonal component to averaged normal, align edge and shift to midplane
                double parallel_component = inner_prod((apical_point - basal_point), face_normal);
                c_vector<double, 3> orthogonal_vect = (apical_point - basal_point) - parallel_component * face_normal;

                if (norm_2(orthogonal_vect) > max_jiggle * 10.0)
                {
                    condition_satisfied = false;
                }
            }
        }
    }
    *****/
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SimplySupportedEdgeBoundaryCondition<1, 1>::ImposeBoundaryCondition(const std::map<Node<1>*, c_vector<double, 1> >& rOldLocations)
{
    EXCEPTION("SimplySupportedEdgeBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SimplySupportedEdgeBoundaryCondition<1, 2>::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    EXCEPTION("SimplySupportedEdgeBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SimplySupportedEdgeBoundaryCondition<1, 3>::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{
    EXCEPTION("SimplySupportedEdgeBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SimplySupportedEdgeBoundaryCondition<2, 2>::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    EXCEPTION("SimplySupportedEdgeBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SimplySupportedEdgeBoundaryCondition<2, 3>::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{
    EXCEPTION("SimplySupportedEdgeBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool SimplySupportedEdgeBoundaryCondition<1, 1>::VerifyBoundaryCondition()
{
    EXCEPTION("SimplySupportedEdgeBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool SimplySupportedEdgeBoundaryCondition<1, 2>::VerifyBoundaryCondition()
{
    EXCEPTION("SimplySupportedEdgeBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool SimplySupportedEdgeBoundaryCondition<1, 3>::VerifyBoundaryCondition()
{
    EXCEPTION("SimplySupportedEdgeBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool SimplySupportedEdgeBoundaryCondition<2, 2>::VerifyBoundaryCondition()
{
    EXCEPTION("SimplySupportedEdgeBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool SimplySupportedEdgeBoundaryCondition<2, 3>::VerifyBoundaryCondition()
{
    EXCEPTION("SimplySupportedEdgeBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SimplySupportedEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MidPlanePoint>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mMidPlanePoint[index] << ",";
    }
    *rParamsFile << mMidPlanePoint[SPACE_DIM - 1] << "</MidPlanePoint>\n";

    *rParamsFile << "\t\t\t<MidPlaneNormal>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mMidPlaneNormal[index] << ",";
    }
    *rParamsFile << mMidPlaneNormal[SPACE_DIM - 1] << "</MidPlaneNormal>\n";

    *rParamsFile << "\t\t\t<UseJiggledNodesOnSubstrate>" << mUseJiggledNodesOnSubstrate << "</UseJiggledNodesOnSubstrate>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class SimplySupportedEdgeBoundaryCondition<1, 1>;
template class SimplySupportedEdgeBoundaryCondition<1, 2>;
template class SimplySupportedEdgeBoundaryCondition<2, 2>;
template class SimplySupportedEdgeBoundaryCondition<1, 3>;
template class SimplySupportedEdgeBoundaryCondition<2, 3>;
template class SimplySupportedEdgeBoundaryCondition<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SimplySupportedEdgeBoundaryCondition)
