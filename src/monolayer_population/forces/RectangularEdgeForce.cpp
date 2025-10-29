/*
Copyright (c) 2005-2019, University of Oxford.
Copyright (c) 2025, Oliver M. Drozdowski and Ulrich S. Schwarz (Heidelberg University)

All rights reserved.
University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.
This file is part of Chaste.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "RectangularEdgeForce.hpp"
#include "LinearSystem.hpp"
#include "MonolayerVertexElement.hpp"
#include "PetscVecTools.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
RectangularEdgeForce<DIM>::RectangularEdgeForce()
        : AbstractForce<DIM>(),
          mpSimulation(nullptr)
{
    mLeftEdgeForce = zero_vector<double>(DIM);
    mRightEdgeForce = zero_vector<double>(DIM);
    mTopEdgeForce = zero_vector<double>(DIM);
    mBottomEdgeForce = zero_vector<double>(DIM);
}

template <unsigned DIM>
RectangularEdgeForce<DIM>::~RectangularEdgeForce()
{
}

template <unsigned DIM>
void RectangularEdgeForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    ///\todo: check whether this line influences profiling tests - if so, we should remove it.
    if (dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("RectangularEdgeForce is to be used with a 3D MonolayerVertexBasedCellPopulation only");
    }

    // Define some helper variables
    /*MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableMonolayerVertexMesh<DIM,DIM>& r_mesh = p_cell_population->rGetMesh();
    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();
    unsigned num_faces = p_cell_population->rGetMesh().GetNumFaces();
    */

    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableMonolayerVertexMesh<DIM, DIM>& r_mesh = p_cell_population->rGetMesh();

    // First we transform into internal coordinates if need be
    c_vector<double, DIM> left_force = zero_vector<double>(DIM);
    c_vector<double, DIM> top_force = zero_vector<double>(DIM);
    c_vector<double, DIM> right_force = zero_vector<double>(DIM);
    c_vector<double, DIM> bottom_force = zero_vector<double>(DIM);
    if (mUseInternalCoordinates)
    {
        // Left edge
        c_vector<double, DIM> population_center = rCellPopulation.GetCentroidOfCellPopulation();
        c_vector<double, DIM> average_edge_normal = zero_vector<double>(DIM);
        c_vector<double, DIM> average_apical_normal = zero_vector<double>(DIM);
        unsigned num_edge_nodes = 0;
        for (auto it = mLeftEdgeNodes.begin(); it != mLeftEdgeNodes.end(); ++it)
        {
            // Edge normal
            c_vector<double, DIM> edge_center = (*it)->rGetLocation();
            c_vector<double, DIM> edge_nrml = population_center - edge_center;
            edge_nrml /= norm_2(edge_nrml);
            average_edge_normal += edge_nrml;
            num_edge_nodes++;

            // Apical vector
            unsigned first_element_index = *((*it)->ContainingElementsBegin());
            MonolayerVertexElement<DIM - 1, DIM>* p_first_elements_face = r_mesh.GetFaceOfType(first_element_index, MonolayerVertexElementType::Apical);
            c_vector<double, DIM> temp_nrml_vect = zero_vector<double>(DIM);
            r_mesh.CalculateUnitNormalToFaceWithArea(p_first_elements_face, temp_nrml_vect);
            temp_nrml_vect /= norm_2(temp_nrml_vect);
            average_apical_normal += temp_nrml_vect;
        }
        average_edge_normal /= norm_2(average_edge_normal);
        average_apical_normal /= norm_2(average_apical_normal);

        // We orthogonalize the edge normal w.r.t. the apical normal
        average_edge_normal -= inner_prod(average_edge_normal, average_apical_normal) * average_apical_normal;
        average_edge_normal /= norm_2(average_edge_normal);

        c_vector<double, DIM> tangential_vector = VectorProduct(average_edge_normal, average_apical_normal);
        tangential_vector /= norm_2(tangential_vector);

        left_force = mLeftEdgeForce[0] * average_edge_normal + mLeftEdgeForce[1] * average_apical_normal + mLeftEdgeForce[2] * tangential_vector;

        // Top edge
        average_edge_normal = zero_vector<double>(DIM);
        average_apical_normal = zero_vector<double>(DIM);
        num_edge_nodes = 0;
        for (auto it = mTopEdgeNodes.begin(); it != mTopEdgeNodes.end(); ++it)
        {
            // Edge normal
            c_vector<double, DIM> edge_center = (*it)->rGetLocation();
            c_vector<double, DIM> edge_nrml = population_center - edge_center;
            edge_nrml /= norm_2(edge_nrml);
            average_edge_normal += edge_nrml;
            num_edge_nodes++;

            // Apical vector
            unsigned first_element_index = *((*it)->ContainingElementsBegin());
            MonolayerVertexElement<DIM - 1, DIM>* p_first_elements_face = r_mesh.GetFaceOfType(first_element_index, MonolayerVertexElementType::Apical);
            c_vector<double, DIM> temp_nrml_vect = zero_vector<double>(DIM);
            r_mesh.CalculateUnitNormalToFaceWithArea(p_first_elements_face, temp_nrml_vect);
            temp_nrml_vect /= norm_2(temp_nrml_vect);
            average_apical_normal += temp_nrml_vect;
        }
        average_edge_normal /= norm_2(average_edge_normal);
        average_apical_normal /= norm_2(average_apical_normal);

        // We orthogonalize the edge normal w.r.t. the apical normal
        average_edge_normal -= inner_prod(average_edge_normal, average_apical_normal) * average_apical_normal;
        average_edge_normal /= norm_2(average_edge_normal);

        tangential_vector = VectorProduct(average_edge_normal, average_apical_normal);
        tangential_vector /= norm_2(tangential_vector);

        top_force = mTopEdgeForce[0] * average_edge_normal + mTopEdgeForce[1] * average_apical_normal + mTopEdgeForce[2] * tangential_vector;

        // Right edge
        average_edge_normal = zero_vector<double>(DIM);
        average_apical_normal = zero_vector<double>(DIM);
        num_edge_nodes = 0;
        for (auto it = mRightEdgeNodes.begin(); it != mRightEdgeNodes.end(); ++it)
        {
            // Edge normal
            c_vector<double, DIM> edge_center = (*it)->rGetLocation();
            c_vector<double, DIM> edge_nrml = population_center - edge_center;
            edge_nrml /= norm_2(edge_nrml);
            average_edge_normal += edge_nrml;
            num_edge_nodes++;

            // Apical vector
            unsigned first_element_index = *((*it)->ContainingElementsBegin());
            MonolayerVertexElement<DIM - 1, DIM>* p_first_elements_face = r_mesh.GetFaceOfType(first_element_index, MonolayerVertexElementType::Apical);
            c_vector<double, DIM> temp_nrml_vect = zero_vector<double>(DIM);
            r_mesh.CalculateUnitNormalToFaceWithArea(p_first_elements_face, temp_nrml_vect);
            temp_nrml_vect /= norm_2(temp_nrml_vect);
            average_apical_normal += temp_nrml_vect;
        }
        average_edge_normal /= norm_2(average_edge_normal);
        average_apical_normal /= norm_2(average_apical_normal);

        // We orthogonalize the edge normal w.r.t. the apical normal
        average_edge_normal -= inner_prod(average_edge_normal, average_apical_normal) * average_apical_normal;
        average_edge_normal /= norm_2(average_edge_normal);

        tangential_vector = VectorProduct(average_edge_normal, average_apical_normal);
        tangential_vector /= norm_2(tangential_vector);

        right_force = mRightEdgeForce[0] * average_edge_normal + mRightEdgeForce[1] * average_apical_normal + mRightEdgeForce[2] * tangential_vector;

        // Bottom edge
        average_edge_normal = zero_vector<double>(DIM);
        average_apical_normal = zero_vector<double>(DIM);
        num_edge_nodes = 0;
        for (auto it = mBottomEdgeNodes.begin(); it != mBottomEdgeNodes.end(); ++it)
        {
            // Edge normal
            c_vector<double, DIM> edge_center = (*it)->rGetLocation();
            c_vector<double, DIM> edge_nrml = population_center - edge_center;
            edge_nrml /= norm_2(edge_nrml);
            average_edge_normal += edge_nrml;
            num_edge_nodes++;

            // Apical vector
            unsigned first_element_index = *((*it)->ContainingElementsBegin());
            MonolayerVertexElement<DIM - 1, DIM>* p_first_elements_face = r_mesh.GetFaceOfType(first_element_index, MonolayerVertexElementType::Apical);
            c_vector<double, DIM> temp_nrml_vect = zero_vector<double>(DIM);
            r_mesh.CalculateUnitNormalToFaceWithArea(p_first_elements_face, temp_nrml_vect);
            temp_nrml_vect /= norm_2(temp_nrml_vect);
            average_apical_normal += temp_nrml_vect;
        }
        average_edge_normal /= norm_2(average_edge_normal);
        average_apical_normal /= norm_2(average_apical_normal);

        // We orthogonalize the edge normal w.r.t. the apical normal
        average_edge_normal -= inner_prod(average_edge_normal, average_apical_normal) * average_apical_normal;
        average_edge_normal /= norm_2(average_edge_normal);

        tangential_vector = VectorProduct(average_edge_normal, average_apical_normal);
        tangential_vector /= norm_2(tangential_vector);

        bottom_force = mBottomEdgeForce[0] * average_edge_normal + mBottomEdgeForce[1] * average_apical_normal + mBottomEdgeForce[2] * tangential_vector;
    }
    else
    {
        left_force = mLeftEdgeForce;
        top_force = mTopEdgeForce;
        right_force = mRightEdgeForce;
        bottom_force = mBottomEdgeForce;
    }

    for (auto it = mLeftEdgeNodes.begin(); it != mLeftEdgeNodes.end(); ++it)
    {
        (*it)->AddAppliedForceContribution(left_force);
    }
    for (auto it = mRightEdgeNodes.begin(); it != mRightEdgeNodes.end(); ++it)
    {
        (*it)->AddAppliedForceContribution(right_force);
    }
    for (auto it = mTopEdgeNodes.begin(); it != mTopEdgeNodes.end(); ++it)
    {
        (*it)->AddAppliedForceContribution(top_force);
    }
    for (auto it = mBottomEdgeNodes.begin(); it != mBottomEdgeNodes.end(); ++it)
    {
        (*it)->AddAppliedForceContribution(bottom_force);
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void RectangularEdgeForce<1>::AddForceContribution(AbstractCellPopulation<1>& rCellPopulation)
{
    EXCEPTION("RectangularEdgeForce not possible in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
std::vector<Node<DIM>*> RectangularEdgeForce<DIM>::GetLeftEdgeNodes()
{
    return mLeftEdgeNodes;
}

template <unsigned DIM>
std::vector<Node<DIM>*> RectangularEdgeForce<DIM>::GetRightEdgeNodes()
{
    return mRightEdgeNodes;
}

template <unsigned DIM>
std::vector<Node<DIM>*> RectangularEdgeForce<DIM>::GetTopEdgeNodes()
{
    return mTopEdgeNodes;
}

template <unsigned DIM>
std::vector<Node<DIM>*> RectangularEdgeForce<DIM>::GetBottomEdgeNodes()
{
    return mBottomEdgeNodes;
}

template <unsigned DIM>
void RectangularEdgeForce<DIM>::SetLeftEdgeNodes(std::vector<Node<DIM>*> edgeNodes)
{
    mLeftEdgeNodes = edgeNodes;
}

template <unsigned DIM>
void RectangularEdgeForce<DIM>::SetRightEdgeNodes(std::vector<Node<DIM>*> edgeNodes)
{
    mRightEdgeNodes = edgeNodes;
}

template <unsigned DIM>
void RectangularEdgeForce<DIM>::SetTopEdgeNodes(std::vector<Node<DIM>*> edgeNodes)
{
    mTopEdgeNodes = edgeNodes;
}

template <unsigned DIM>
void RectangularEdgeForce<DIM>::SetBottomEdgeNodes(std::vector<Node<DIM>*> edgeNodes)
{
    mBottomEdgeNodes = edgeNodes;
}

template <unsigned DIM>
c_vector<double, DIM> RectangularEdgeForce<DIM>::GetLeftEdgeForce()
{
    return mLeftEdgeForce;
}

template <unsigned DIM>
c_vector<double, DIM> RectangularEdgeForce<DIM>::GetRightEdgeForce()
{
    return mRightEdgeForce;
}

template <unsigned DIM>
c_vector<double, DIM> RectangularEdgeForce<DIM>::GetTopEdgeForce()
{
    return mTopEdgeForce;
}

template <unsigned DIM>
c_vector<double, DIM> RectangularEdgeForce<DIM>::GetBottomEdgeForce()
{
    return mBottomEdgeForce;
}

template <unsigned DIM>
void RectangularEdgeForce<DIM>::SetLeftEdgeForce(c_vector<double, DIM> edgeForce)
{
    mLeftEdgeForce = edgeForce;
}

template <unsigned DIM>
void RectangularEdgeForce<DIM>::SetRightEdgeForce(c_vector<double, DIM> edgeForce)
{
    mRightEdgeForce = edgeForce;
}

template <unsigned DIM>
void RectangularEdgeForce<DIM>::SetTopEdgeForce(c_vector<double, DIM> edgeForce)
{
    mTopEdgeForce = edgeForce;
}

template <unsigned DIM>
void RectangularEdgeForce<DIM>::SetBottomEdgeForce(c_vector<double, DIM> edgeForce)
{
    mBottomEdgeForce = edgeForce;
}

template <unsigned DIM>
void RectangularEdgeForce<DIM>::UseInternalCoordinates(bool useInternalCoordinates)
{
    mUseInternalCoordinates = useInternalCoordinates;
}

template <unsigned DIM>
void RectangularEdgeForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
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

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class RectangularEdgeForce<1>;
template class RectangularEdgeForce<2>;
template class RectangularEdgeForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"