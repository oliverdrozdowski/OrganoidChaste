#include "RectangularEdgeSubForce.hpp"
#include "LinearSystem.hpp"
#include "MonolayerVertexElement.hpp"
#include "PetscVecTools.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
RectangularEdgeSubForce<DIM>::RectangularEdgeSubForce()
        : AbstractSubForce<DIM>(),
          mpSimulation(nullptr)
{
    mLeftEdgeForce = zero_vector<double>(DIM);
    mRightEdgeForce = zero_vector<double>(DIM);
    mTopEdgeForce = zero_vector<double>(DIM);
    mBottomEdgeForce = zero_vector<double>(DIM);
}

template <unsigned DIM>
RectangularEdgeSubForce<DIM>::~RectangularEdgeSubForce()
{
}

template <unsigned DIM>
NodeSubForceMap<DIM> RectangularEdgeSubForce<DIM>::GetForceContributions(
    MonolayerVertexBasedCellPopulation<DIM>& rCellPopulation, AbstractCellBasedSimulation<DIM, DIM>* pSimulation)
{
    MutableMonolayerVertexMesh<DIM, DIM>& r_mesh = rCellPopulation.rGetMesh();

    // First we transform into internal coordinates if need be
    c_vector<double, DIM> left_force = zero_vector<double>(DIM);
    c_vector<double, DIM> top_force = zero_vector<double>(DIM);
    c_vector<double, DIM> right_force = zero_vector<double>(DIM);
    c_vector<double, DIM> bottom_force = zero_vector<double>(DIM);
    if (mUseInternalCoordinates)
    {
        // Left edge
        c_vector<double, DIM> population_center = rCellPopulation.GetCentroidOfCellPopulation();
        c_vector<double, DIM> average_edge_normal_l = zero_vector<double>(DIM);
        c_vector<double, DIM> average_apical_normal_l = zero_vector<double>(DIM);
        unsigned num_edge_nodes_l = 0;
        for (auto it = mLeftEdgeNodes.begin(); it != mLeftEdgeNodes.end(); ++it)
        {
            // Edge normal
            c_vector<double, DIM> edge_center = (*it)->rGetLocation();
            c_vector<double, DIM> edge_nrml = population_center - edge_center;
            edge_nrml /= norm_2(edge_nrml);
            average_edge_normal_l += edge_nrml;
            num_edge_nodes_l++;

            // Apical vector
            unsigned first_element_index = *((*it)->ContainingElementsBegin());
            MonolayerVertexElement<DIM - 1, DIM>* p_first_elements_face = r_mesh.GetFaceOfType(first_element_index, MonolayerVertexElementType::Apical);
            c_vector<double, DIM> temp_nrml_vect = zero_vector<double>(DIM);
            r_mesh.CalculateUnitNormalToFaceWithArea(p_first_elements_face, temp_nrml_vect);
            temp_nrml_vect /= norm_2(temp_nrml_vect);
            average_apical_normal_l += temp_nrml_vect;
        }
        average_edge_normal_l /= norm_2(average_edge_normal_l);
        average_apical_normal_l /= norm_2(average_apical_normal_l);

        // RIGHT EDGE
        // Left edge
        c_vector<double, DIM> average_edge_normal_r = zero_vector<double>(DIM);
        c_vector<double, DIM> average_apical_normal_r = zero_vector<double>(DIM);
        unsigned num_edge_nodes_r = 0;
        for (auto it = mRightEdgeNodes.begin(); it != mRightEdgeNodes.end(); ++it)
        {
            // Edge normal
            c_vector<double, DIM> edge_center = (*it)->rGetLocation();
            c_vector<double, DIM> edge_nrml = population_center - edge_center;
            edge_nrml /= norm_2(edge_nrml);
            average_edge_normal_r += edge_nrml;
            num_edge_nodes_r++;

            // Apical vector
            unsigned first_element_index = *((*it)->ContainingElementsBegin());
            MonolayerVertexElement<DIM - 1, DIM>* p_first_elements_face = r_mesh.GetFaceOfType(first_element_index, MonolayerVertexElementType::Apical);
            c_vector<double, DIM> temp_nrml_vect = zero_vector<double>(DIM);
            r_mesh.CalculateUnitNormalToFaceWithArea(p_first_elements_face, temp_nrml_vect);
            temp_nrml_vect /= norm_2(temp_nrml_vect);
            average_apical_normal_r += temp_nrml_vect;
        }
        average_edge_normal_r /= norm_2(average_edge_normal_r);
        average_apical_normal_r /= norm_2(average_apical_normal_r);

        // If we set the normals to parallel externally, we have to correct for this here
        if (mParallelizeOppositeEdgeNormals)
        {
            c_vector<double, DIM> average_edge_normal = zero_vector<double>(DIM);
            c_vector<double, DIM> average_apical_normal = zero_vector<double>(DIM);

            // weighted averages, negative because of opposite directions
            average_edge_normal = num_edge_nodes_l * average_edge_normal_l - num_edge_nodes_r * average_edge_normal_r;
            average_apical_normal = num_edge_nodes_l * average_apical_normal_l + num_edge_nodes_r * average_apical_normal_r;

            average_edge_normal /= norm_2(average_edge_normal);
            average_apical_normal /= norm_2(average_apical_normal);

            // We orthogonalize the edge normal w.r.t. the apical normal
            average_edge_normal -= inner_prod(average_edge_normal, average_apical_normal) * average_apical_normal;
            average_edge_normal /= norm_2(average_edge_normal);

            c_vector<double, DIM> tangential_vector = VectorProduct(average_edge_normal, average_apical_normal);
            tangential_vector /= norm_2(tangential_vector);

            left_force = mLeftEdgeForce[0] * average_edge_normal + mLeftEdgeForce[1] * average_apical_normal + mLeftEdgeForce[2] * tangential_vector;

            right_force = -mRightEdgeForce[0] * average_edge_normal + mRightEdgeForce[1] * average_apical_normal + -mRightEdgeForce[2] * tangential_vector;
        }
        else
        {
            // Left edge
            // We orthogonalize the edge normal w.r.t. the apical normal
            average_edge_normal_l -= inner_prod(average_edge_normal_l, average_apical_normal_l) * average_apical_normal_l;
            average_edge_normal_l /= norm_2(average_edge_normal_l);

            c_vector<double, DIM> tangential_vector_l = VectorProduct(average_edge_normal_l, average_apical_normal_l);
            tangential_vector_l /= norm_2(tangential_vector_l);

            left_force = mLeftEdgeForce[0] * average_edge_normal_l + mLeftEdgeForce[1] * average_apical_normal_l + mLeftEdgeForce[2] * tangential_vector_l;

            // Right edge
            // We orthogonalize the edge normal w.r.t. the apical normal
            average_edge_normal_r -= inner_prod(average_edge_normal_r, average_apical_normal_r) * average_apical_normal_r;
            average_edge_normal_r /= norm_2(average_edge_normal_r);

            c_vector<double, DIM> tangential_vector_r = VectorProduct(average_edge_normal_r, average_apical_normal_r);
            tangential_vector_r /= norm_2(tangential_vector_r);

            right_force = mRightEdgeForce[0] * average_edge_normal_r + mRightEdgeForce[1] * average_apical_normal_r + mRightEdgeForce[2] * tangential_vector_r;
        }

        // Top edge
        c_vector<double, DIM> average_edge_normal_tp = zero_vector<double>(DIM);
        c_vector<double, DIM> average_apical_normal_tp = zero_vector<double>(DIM);
        unsigned num_edge_nodes_tp = 0;
        for (auto it = mTopEdgeNodes.begin(); it != mTopEdgeNodes.end(); ++it)
        {
            // Edge normal
            c_vector<double, DIM> edge_center = (*it)->rGetLocation();
            c_vector<double, DIM> edge_nrml = population_center - edge_center;
            edge_nrml /= norm_2(edge_nrml);
            average_edge_normal_tp += edge_nrml;
            num_edge_nodes_tp++;

            // Apical vector
            unsigned first_element_index = *((*it)->ContainingElementsBegin());
            MonolayerVertexElement<DIM - 1, DIM>* p_first_elements_face = r_mesh.GetFaceOfType(first_element_index, MonolayerVertexElementType::Apical);
            c_vector<double, DIM> temp_nrml_vect = zero_vector<double>(DIM);
            r_mesh.CalculateUnitNormalToFaceWithArea(p_first_elements_face, temp_nrml_vect);
            temp_nrml_vect /= norm_2(temp_nrml_vect);
            average_apical_normal_tp += temp_nrml_vect;
        }
        average_edge_normal_tp /= norm_2(average_edge_normal_tp);
        average_apical_normal_tp /= norm_2(average_apical_normal_tp);

        // Bottom edge
        c_vector<double, DIM> average_edge_normal_b = zero_vector<double>(DIM);
        c_vector<double, DIM> average_apical_normal_b = zero_vector<double>(DIM);
        unsigned num_edge_nodes_b = 0;
        for (auto it = mBottomEdgeNodes.begin(); it != mBottomEdgeNodes.end(); ++it)
        {
            // Edge normal
            c_vector<double, DIM> edge_center = (*it)->rGetLocation();
            c_vector<double, DIM> edge_nrml = population_center - edge_center;
            edge_nrml /= norm_2(edge_nrml);
            average_edge_normal_b += edge_nrml;
            num_edge_nodes_b++;

            // Apical vector
            unsigned first_element_index = *((*it)->ContainingElementsBegin());
            MonolayerVertexElement<DIM - 1, DIM>* p_first_elements_face = r_mesh.GetFaceOfType(first_element_index, MonolayerVertexElementType::Apical);
            c_vector<double, DIM> temp_nrml_vect = zero_vector<double>(DIM);
            r_mesh.CalculateUnitNormalToFaceWithArea(p_first_elements_face, temp_nrml_vect);
            temp_nrml_vect /= norm_2(temp_nrml_vect);
            average_apical_normal_b += temp_nrml_vect;
        }
        average_edge_normal_b /= norm_2(average_edge_normal_b);
        average_apical_normal_b /= norm_2(average_apical_normal_b);

        // If we set the normals to parallel externally, we have to correct for this here
        if (mParallelizeOppositeEdgeNormals)
        {
            c_vector<double, DIM> average_edge_normal = zero_vector<double>(DIM);
            c_vector<double, DIM> average_apical_normal = zero_vector<double>(DIM);

            // weighted averages, negative because of opposite directions
            average_edge_normal = num_edge_nodes_tp * average_edge_normal_tp - num_edge_nodes_b * average_edge_normal_b;
            average_apical_normal = num_edge_nodes_tp * average_apical_normal_tp + num_edge_nodes_b * average_apical_normal_b;

            average_edge_normal /= norm_2(average_edge_normal);
            average_apical_normal /= norm_2(average_apical_normal);

            // We orthogonalize the edge normal w.r.t. the apical normal
            average_edge_normal -= inner_prod(average_edge_normal, average_apical_normal) * average_apical_normal;
            average_edge_normal /= norm_2(average_edge_normal);

            c_vector<double, DIM> tangential_vector = VectorProduct(average_edge_normal, average_apical_normal);
            tangential_vector /= norm_2(tangential_vector);

            top_force = mTopEdgeForce[0] * average_edge_normal + mTopEdgeForce[1] * average_apical_normal + mTopEdgeForce[2] * tangential_vector;

            bottom_force = -mBottomEdgeForce[0] * average_edge_normal + mBottomEdgeForce[1] * average_apical_normal + -mBottomEdgeForce[2] * tangential_vector;
        }
        else
        {
            // Top edge
            // We orthogonalize the edge normal w.r.t. the apical normal
            average_edge_normal_tp -= inner_prod(average_edge_normal_tp, average_apical_normal_tp) * average_apical_normal_tp;
            average_edge_normal_tp /= norm_2(average_edge_normal_tp);

            c_vector<double, DIM> tangential_vector_tp = VectorProduct(average_edge_normal_tp, average_apical_normal_tp);
            tangential_vector_tp /= norm_2(tangential_vector_tp);

            top_force = mTopEdgeForce[0] * average_edge_normal_tp + mTopEdgeForce[1] * average_apical_normal_tp + mTopEdgeForce[2] * tangential_vector_tp;

            // Bottom edge
            // We orthogonalize the edge normal w.r.t. the apical normal
            average_edge_normal_b -= inner_prod(average_edge_normal_b, average_apical_normal_b) * average_apical_normal_b;
            average_edge_normal_b /= norm_2(average_edge_normal_b);

            c_vector<double, DIM> tangential_vector_b = VectorProduct(average_edge_normal_b, average_apical_normal_b);
            tangential_vector_b /= norm_2(tangential_vector_b);

            bottom_force = mBottomEdgeForce[0] * average_edge_normal_b + mBottomEdgeForce[1] * average_apical_normal_b + mBottomEdgeForce[2] * tangential_vector_b;
        }
    }
    else
    {
        left_force = mLeftEdgeForce;
        top_force = mTopEdgeForce;
        right_force = mRightEdgeForce;
        bottom_force = mBottomEdgeForce;
    }

    NodeSubForceMap<DIM> map_forces;
    for (auto it = mLeftEdgeNodes.begin(); it != mLeftEdgeNodes.end(); ++it)
    {
        std::set<Node<DIM>*> nodes_already_added;
        if (map_forces.find(*it) == map_forces.end())
        {
            map_forces[(*it)] = left_force;
            nodes_already_added.insert(*it);
        }
        else
        {
            if (nodes_already_added.find(*it) == nodes_already_added.end())
            {
                map_forces[(*it)] += left_force;
            }
        }
    }
    for (auto it = mRightEdgeNodes.begin(); it != mRightEdgeNodes.end(); ++it)
    {
        std::set<Node<DIM>*> nodes_already_added;
        if (map_forces.find(*it) == map_forces.end())
        {
            map_forces[(*it)] = right_force;
            nodes_already_added.insert(*it);
        }
        else
        {
            if (nodes_already_added.find(*it) == nodes_already_added.end())
            {
                map_forces[(*it)] += right_force;
            }
        }
    }
    for (auto it = mTopEdgeNodes.begin(); it != mTopEdgeNodes.end(); ++it)
    {
        std::set<Node<DIM>*> nodes_already_added;
        if (map_forces.find(*it) == map_forces.end())
        {
            map_forces[(*it)] = top_force;
            nodes_already_added.insert(*it);
        }
        else
        {
            if (nodes_already_added.find(*it) == nodes_already_added.end())
            {
                map_forces[(*it)] += top_force;
            }
        }
    }
    for (auto it = mBottomEdgeNodes.begin(); it != mBottomEdgeNodes.end(); ++it)
    {
        std::set<Node<DIM>*> nodes_already_added;
        if (map_forces.find(*it) == map_forces.end())
        {
            map_forces[(*it)] = bottom_force;
            nodes_already_added.insert(*it);
        }
        else
        {
            if (nodes_already_added.find(*it) == nodes_already_added.end())
            {
                map_forces[(*it)] += bottom_force;
            }
        }
    }
    return map_forces;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
NodeSubForceMap<1> RectangularEdgeSubForce<1>::GetForceContributions(
    MonolayerVertexBasedCellPopulation<1>& rCellPopulation, AbstractCellBasedSimulation<1, 1>* pSimulation)
{
    EXCEPTION("RectangularEdgeSubForce not possible in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
std::vector<Node<DIM>*> RectangularEdgeSubForce<DIM>::GetLeftEdgeNodes()
{
    return mLeftEdgeNodes;
}

template <unsigned DIM>
std::vector<Node<DIM>*> RectangularEdgeSubForce<DIM>::GetRightEdgeNodes()
{
    return mRightEdgeNodes;
}

template <unsigned DIM>
std::vector<Node<DIM>*> RectangularEdgeSubForce<DIM>::GetTopEdgeNodes()
{
    return mTopEdgeNodes;
}

template <unsigned DIM>
std::vector<Node<DIM>*> RectangularEdgeSubForce<DIM>::GetBottomEdgeNodes()
{
    return mBottomEdgeNodes;
}

template <unsigned DIM>
void RectangularEdgeSubForce<DIM>::SetLeftEdgeNodes(std::vector<Node<DIM>*> edgeNodes)
{
    mLeftEdgeNodes = edgeNodes;
}

template <unsigned DIM>
void RectangularEdgeSubForce<DIM>::SetRightEdgeNodes(std::vector<Node<DIM>*> edgeNodes)
{
    mRightEdgeNodes = edgeNodes;
}

template <unsigned DIM>
void RectangularEdgeSubForce<DIM>::SetTopEdgeNodes(std::vector<Node<DIM>*> edgeNodes)
{
    mTopEdgeNodes = edgeNodes;
}

template <unsigned DIM>
void RectangularEdgeSubForce<DIM>::SetBottomEdgeNodes(std::vector<Node<DIM>*> edgeNodes)
{
    mBottomEdgeNodes = edgeNodes;
}

template <unsigned DIM>
c_vector<double, DIM> RectangularEdgeSubForce<DIM>::GetLeftEdgeForce()
{
    return mLeftEdgeForce;
}

template <unsigned DIM>
c_vector<double, DIM> RectangularEdgeSubForce<DIM>::GetRightEdgeForce()
{
    return mRightEdgeForce;
}

template <unsigned DIM>
c_vector<double, DIM> RectangularEdgeSubForce<DIM>::GetTopEdgeForce()
{
    return mTopEdgeForce;
}

template <unsigned DIM>
c_vector<double, DIM> RectangularEdgeSubForce<DIM>::GetBottomEdgeForce()
{
    return mBottomEdgeForce;
}

template <unsigned DIM>
void RectangularEdgeSubForce<DIM>::SetLeftEdgeForce(c_vector<double, DIM> edgeForce)
{
    mLeftEdgeForce = edgeForce;
}

template <unsigned DIM>
void RectangularEdgeSubForce<DIM>::SetRightEdgeForce(c_vector<double, DIM> edgeForce)
{
    mRightEdgeForce = edgeForce;
}

template <unsigned DIM>
void RectangularEdgeSubForce<DIM>::SetTopEdgeForce(c_vector<double, DIM> edgeForce)
{
    mTopEdgeForce = edgeForce;
}

template <unsigned DIM>
void RectangularEdgeSubForce<DIM>::SetBottomEdgeForce(c_vector<double, DIM> edgeForce)
{
    mBottomEdgeForce = edgeForce;
}

template <unsigned DIM>
void RectangularEdgeSubForce<DIM>::UseInternalCoordinates(bool useInternalCoordinates)
{
    mUseInternalCoordinates = useInternalCoordinates;
}

template <unsigned DIM>
void RectangularEdgeSubForce<DIM>::ParallelizeOppositeEdgeNormals(bool parallelizeOppositeEdgeNormals)
{
    mParallelizeOppositeEdgeNormals = parallelizeOppositeEdgeNormals;
}

template <unsigned DIM>
void RectangularEdgeSubForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<LeftEdgeForce>";
    for (unsigned index = 0; index != DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mLeftEdgeForce[index] << ",";
    }
    *rParamsFile << "</LeftEdgeForce>\n";

    *rParamsFile << "\t\t\t<RightEdgeForce>";
    for (unsigned index = 0; index != DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mRightEdgeForce[index] << ",";
    }
    *rParamsFile << "</RightEdgeForce>\n";

    *rParamsFile << "\t\t\t<TopEdgeForce>";
    for (unsigned index = 0; index != DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mTopEdgeForce[index] << ",";
    }
    *rParamsFile << "</TopEdgeForce>\n";

    *rParamsFile << "\t\t\t<BottomEdgeForce>";
    for (unsigned index = 0; index != DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mBottomEdgeForce[index] << ",";
    }
    *rParamsFile << "</BottomEdgeForce>\n";

    *rParamsFile << "\t\t\t<UseInternalCoordinates>" << mUseInternalCoordinates << "<UseInternalCoordinates>\n";

    *rParamsFile << "\t\t\t<ParallelizeOppositeEdgeNormals>" << mParallelizeOppositeEdgeNormals << "<ParallelizeOppositeEdgeNormals>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class RectangularEdgeSubForce<1>;
template class RectangularEdgeSubForce<2>;
template class RectangularEdgeSubForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"