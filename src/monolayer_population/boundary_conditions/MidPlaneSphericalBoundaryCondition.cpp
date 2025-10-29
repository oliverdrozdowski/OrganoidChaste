#include "MidPlaneSphericalBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MonolayerVertexElement.hpp"
#include "MutableElement.hpp"
#include "VertexBasedCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MidPlaneSphericalBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::MidPlaneSphericalBoundaryCondition(
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
    c_vector<double, SPACE_DIM> midPoint,
    double radius)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(pCellPopulation),
          mMidPoint(midPoint)
{
    assert(radius > 0.0);
    mRadius = radius;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& MidPlaneSphericalBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetMidPoint() const
{
    return mMidPoint;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MidPlaneSphericalBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::GetRadius() const
{
    return mRadius;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MidPlaneSphericalBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::ImposeBoundaryCondition(
    const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations)
{
    ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("MidPlaneSphericalBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }

    assert((dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation))
           || (SPACE_DIM == ELEMENT_DIM && (dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation))));

    if (SPACE_DIM != 1)
    {
        assert(SPACE_DIM == ELEMENT_DIM); // LCOV_EXCL_LINE
        assert(dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation));
        MonolayerVertexBasedCellPopulation<SPACE_DIM>* p_cell_population = dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation);

        MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& r_mesh = p_cell_population->rGetMesh();

        // Loop over lateral faces
        // We correct the midplane.
        for (auto face_iter = r_mesh.GetFaceIteratorBegin(); face_iter != r_mesh.GetFaceIteratorEnd(); ++face_iter)
        {
            // Only check lateral faces
            if (face_iter->GetFaceType() != MonolayerVertexElementType::Lateral)
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

                c_vector<double, SPACE_DIM> mid_point = (basal_point + apical_point) / 2.0;
                c_vector<double, SPACE_DIM> distance_vector = (mid_point - mMidPoint);
                double correction = mRadius - norm_2(distance_vector);

                // Then correct the edge
                basal_point += correction * distance_vector / norm_2(distance_vector);
                apical_point += correction * distance_vector / norm_2(distance_vector);
            }
        }
    }
    else
    {
        // DIM == 1
        NEVER_REACHED;
        // MidPlaneSphericalBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MidPlaneSphericalBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::VerifyBoundaryCondition()
{
    return true;
    if (SPACE_DIM == 1)
    {
        EXCEPTION("MidPlaneSphericalBoundaryCondition is not implemented in 1D");
    }
    else
    {
        // NO CHECK!
        return true;
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MidPlaneSphericalBoundaryCondition<1, 1>::ImposeBoundaryCondition(const std::map<Node<1>*, c_vector<double, 1> >& rOldLocations)
{
    EXCEPTION("MidPlaneSphericalBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MidPlaneSphericalBoundaryCondition<1, 2>::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    EXCEPTION("MidPlaneSphericalBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MidPlaneSphericalBoundaryCondition<1, 3>::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{
    EXCEPTION("MidPlaneSphericalBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MidPlaneSphericalBoundaryCondition<2, 2>::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    EXCEPTION("MidPlaneSphericalBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MidPlaneSphericalBoundaryCondition<2, 3>::ImposeBoundaryCondition(const std::map<Node<3>*, c_vector<double, 3> >& rOldLocations)
{
    EXCEPTION("MidPlaneSphericalBoundaryCondition only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool MidPlaneSphericalBoundaryCondition<1, 1>::VerifyBoundaryCondition()
{
    EXCEPTION("MidPlaneSphericalBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool MidPlaneSphericalBoundaryCondition<1, 2>::VerifyBoundaryCondition()
{
    EXCEPTION("MidPlaneSphericalBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool MidPlaneSphericalBoundaryCondition<1, 3>::VerifyBoundaryCondition()
{
    EXCEPTION("MidPlaneSphericalBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool MidPlaneSphericalBoundaryCondition<2, 2>::VerifyBoundaryCondition()
{
    EXCEPTION("MidPlaneSphericalBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool MidPlaneSphericalBoundaryCondition<2, 3>::VerifyBoundaryCondition()
{
    EXCEPTION("MidPlaneSphericalBoundaryCondition only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MidPlaneSphericalBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MidPoint>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mMidPoint[index] << ",";
    }
    *rParamsFile << mMidPoint[SPACE_DIM - 1] << "</MidPoint>\n";

    *rParamsFile << "\t\t\t<Radius>" << mRadius << "</Radius>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class MidPlaneSphericalBoundaryCondition<1, 1>;
template class MidPlaneSphericalBoundaryCondition<1, 2>;
template class MidPlaneSphericalBoundaryCondition<2, 2>;
template class MidPlaneSphericalBoundaryCondition<1, 3>;
template class MidPlaneSphericalBoundaryCondition<2, 3>;
template class MidPlaneSphericalBoundaryCondition<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MidPlaneSphericalBoundaryCondition)
