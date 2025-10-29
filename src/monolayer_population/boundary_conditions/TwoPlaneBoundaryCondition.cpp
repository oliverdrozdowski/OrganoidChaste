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

#include "TwoPlaneBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MonolayerVertexElement.hpp"
#include "MutableElement.hpp"
#include "VertexBasedCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::TwoPlaneBoundaryCondition(
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
    c_vector<double, SPACE_DIM> basalPoint,
    c_vector<double, SPACE_DIM> basalNormal,
    c_vector<double, SPACE_DIM> apicalPoint,
    c_vector<double, SPACE_DIM> apicalNormal)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(pCellPopulation),
          mBasalPointOnSubstrate(basalPoint),
          mApicalPointOnSubstrate(apicalPoint),
          mUseJiggledNodesOnSubstrate(false)
{
    assert(norm_2(basalNormal) > 0.0);
    mBasalNormalToSubstrate = basalNormal / norm_2(basalNormal);
    assert(norm_2(apicalNormal) > 0.0);
    mApicalNormalToSubstrate = apicalNormal / norm_2(apicalNormal);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetBasalPointOnSubstrate() const
{
    return mBasalPointOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetBasalNormalToSubstrate() const
{
    return mBasalNormalToSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetApicalPointOnSubstrate() const
{
    return mApicalPointOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetApicalNormalToSubstrate() const
{
    return mApicalNormalToSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::SetUseJiggledNodesOnSubstrate(bool useJiggledNodesOnSubstrate)
{
    mUseJiggledNodesOnSubstrate = useJiggledNodesOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::GetUseJiggledNodesOnSubstrate()
{
    return mUseJiggledNodesOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::ImposeBoundaryCondition(
    const std::map<Node<SPACE_DIM>*,
                   c_vector<double, SPACE_DIM> >& rOldLocations)
{
    ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("TwoPlaneBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
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

        // Iterate over all nodes and update their positions according to the boundary conditions
        unsigned num_nodes = this->mpCellPopulation->GetNumNodes();
        for (unsigned node_index = 0; node_index < num_nodes; node_index++)
        {
            Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

            // Only do something for boundary nodes
            if (!p_node->IsBoundaryNode())
            {
                continue;
            }

            c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();
            std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
            unsigned first_element = *(containing_elements.begin());
            unsigned local_node_index = p_cell_population->GetElement(first_element)->GetNodeLocalIndex(node_index);
            MonolayerVertexElementType node_type = p_cell_population->GetElement(first_element)->GetNodeType(local_node_index);

            if (node_type == MonolayerVertexElementType::Basal)
            {
                double signed_basal_distance = inner_prod(node_location - mBasalPointOnSubstrate, mBasalNormalToSubstrate);
                if (abs(signed_basal_distance) > max_jiggle)
                {
                    // For the closest point on the plane we travel from node_location the signed_distance in the direction
                    // of -mNormalToSubstrate
                    c_vector<double, SPACE_DIM> nearest_point;
                    if (mUseJiggledNodesOnSubstrate)
                    {
                        nearest_point = node_location - (signed_basal_distance + max_jiggle * RandomNumberGenerator::Instance()->ranf()) * mBasalNormalToSubstrate;
                    }
                    else
                    {
                        nearest_point = node_location - signed_basal_distance * mBasalNormalToSubstrate;
                    }
                    p_node->rGetModifiableLocation() = nearest_point;
                }
            }
            else if (node_type == MonolayerVertexElementType::Apical)
            {
                double signed_apical_distance = inner_prod(node_location - mApicalPointOnSubstrate, mApicalNormalToSubstrate);
                if (abs(signed_apical_distance) > max_jiggle)
                {
                    // For the closest point on the plane we travel from node_location the signed_distance in the direction
                    // of -mNormalToSubstrate
                    c_vector<double, SPACE_DIM> nearest_point;
                    if (mUseJiggledNodesOnSubstrate)
                    {
                        nearest_point = node_location - (signed_apical_distance + max_jiggle * RandomNumberGenerator::Instance()->ranf()) * mApicalNormalToSubstrate;
                    }
                    else
                    {
                        nearest_point = node_location - signed_apical_distance * mApicalNormalToSubstrate;
                    }
                    p_node->rGetModifiableLocation() = nearest_point;
                }
            }
        }
    }
    else
    {
        // DIM == 1
        NEVER_REACHED;
        // TwoPlaneBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    if (SPACE_DIM == 1)
    {
        EXCEPTION("TwoPlaneBoundaryCondition is not implemented in 1D");
    }
    else
    {
        assert(dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation));
        MonolayerVertexBasedCellPopulation<SPACE_DIM>* p_cell_population = dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation);

        // This is a magic number
        double max_jiggle = 1e-4;

        unsigned num_nodes = p_cell_population->GetNumNodes();
        for (unsigned node_index = 0; node_index < num_nodes; node_index++)
        {
            Node<SPACE_DIM>* p_node = p_cell_population->GetNode(node_index);

            // Boundary condition only applies to boundary nodes
            if (!p_node->IsBoundaryNode())
            {
                continue;
            }

            c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();
            std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
            unsigned first_element = *(containing_elements.begin());
            unsigned local_node_index = p_cell_population->GetElement(first_element)->GetNodeLocalIndex(node_index);
            MonolayerVertexElementType node_type = p_cell_population->GetElement(first_element)->GetNodeType(local_node_index);

            if (node_type == MonolayerVertexElementType::Basal)
            {
                if (abs(inner_prod(node_location - mBasalPointOnSubstrate, mBasalNormalToSubstrate)) > max_jiggle)
                {
                    condition_satisfied = false;
                    break;
                }
            }
            else if (node_type == MonolayerVertexElementType::Apical)
            {
                if (abs(inner_prod(node_location - mApicalPointOnSubstrate, mApicalNormalToSubstrate)) > max_jiggle)
                {
                    condition_satisfied = false;
                    break;
                }
            }
        }
    }
    return condition_satisfied;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TwoPlaneBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<BasalPointOnSubstrate>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mBasalPointOnSubstrate[index] << ",";
    }
    *rParamsFile << mBasalPointOnSubstrate[SPACE_DIM - 1] << "</BasalPointOnSubstrate>\n";

    *rParamsFile << "\t\t\t<BasalNormalToSubstrate>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mBasalNormalToSubstrate[index] << ",";
    }
    *rParamsFile << mBasalNormalToSubstrate[SPACE_DIM - 1] << "</BasalNormalToSubstrate>\n";

    *rParamsFile << "\t\t\t<ApicalPointOnSubstrate>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mApicalPointOnSubstrate[index] << ",";
    }
    *rParamsFile << mApicalPointOnSubstrate[SPACE_DIM - 1] << "</ApicalPointOnSubstrate>\n";

    *rParamsFile << "\t\t\t<ApicalNormalToSubstrate>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mApicalNormalToSubstrate[index] << ",";
    }
    *rParamsFile << mApicalNormalToSubstrate[SPACE_DIM - 1] << "</ApicalNormalToSubstrate>\n";

    *rParamsFile << "\t\t\t<UseJiggledNodesOnSubstrate>" << mUseJiggledNodesOnSubstrate << "</UseJiggledNodesOnSubstrate>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class TwoPlaneBoundaryCondition<1, 1>;
template class TwoPlaneBoundaryCondition<1, 2>;
template class TwoPlaneBoundaryCondition<2, 2>;
template class TwoPlaneBoundaryCondition<1, 3>;
template class TwoPlaneBoundaryCondition<2, 3>;
template class TwoPlaneBoundaryCondition<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(TwoPlaneBoundaryCondition)
