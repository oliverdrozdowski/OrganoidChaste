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

#include "SubstrateBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MonolayerVertexElement.hpp"
#include "MutableElement.hpp"
#include "VertexBasedCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SubstrateBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::SubstrateBoundaryCondition(
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
    c_vector<double, SPACE_DIM> point,
    c_vector<double, SPACE_DIM> normal)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(pCellPopulation),
          mPointOnSubstrate(point),
          mUseJiggledNodesOnSubstrate(false)
{
    assert(norm_2(normal) > 0.0);
    mNormalToSubstrate = normal / norm_2(normal);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& SubstrateBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetPointOnSubstrate() const
{
    return mPointOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& SubstrateBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::rGetNormalToSubstrate() const
{
    return mNormalToSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SubstrateBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::SetUseJiggledNodesOnSubstrate(bool useJiggledNodesOnSubstrate)
{
    mUseJiggledNodesOnSubstrate = useJiggledNodesOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool SubstrateBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::GetUseJiggledNodesOnSubstrate()
{
    return mUseJiggledNodesOnSubstrate;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SubstrateBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::ImposeBoundaryCondition(
    const std::map<Node<SPACE_DIM>*,
                   c_vector<double, SPACE_DIM> >& rOldLocations)
{
    ///\todo Move this to constructor. If this is in the constructor then Exception always throws.
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>*>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("SubstrateBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
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
            c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();
            std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
            unsigned first_element = *(containing_elements.begin());
            unsigned local_node_index = p_cell_population->GetElement(first_element)->GetNodeLocalIndex(node_index);
            MonolayerVertexElementType node_type = p_cell_population->GetElement(first_element)->GetNodeType(local_node_index);

            if (node_type != MonolayerVertexElementType::Basal)
            {
                continue;
            }

            double signed_distance = inner_prod(node_location - mPointOnSubstrate, mNormalToSubstrate);
            if (signed_distance != 0.0)
            {
                // For the closest point on the plane we travel from node_location the signed_distance in the direction of -mNormalToSubstrate
                c_vector<double, SPACE_DIM> nearest_point;
                if (mUseJiggledNodesOnSubstrate)
                {
                    nearest_point = node_location - (signed_distance + max_jiggle * RandomNumberGenerator::Instance()->ranf()) * mNormalToSubstrate;
                }
                else
                {
                    nearest_point = node_location - signed_distance * mNormalToSubstrate;
                }
                p_node->rGetModifiableLocation() = nearest_point;
            }
        }
    }
    else
    {
        // DIM == 1
        NEVER_REACHED;
        // SubstrateBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool SubstrateBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    if (SPACE_DIM == 1)
    {
        EXCEPTION("SubstrateBoundaryCondition is not implemented in 1D");
    }
    else
    {
        assert(dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation));
        MonolayerVertexBasedCellPopulation<SPACE_DIM>* p_cell_population = dynamic_cast<MonolayerVertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation);

        unsigned num_nodes = p_cell_population->GetNumNodes();
        for (unsigned node_index = 0; node_index < num_nodes; node_index++)
        {
            Node<SPACE_DIM>* p_node = p_cell_population->GetNode(node_index);
            c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();
            std::set<unsigned> containing_elements = p_node->rGetContainingElementIndices();
            unsigned first_element = *(containing_elements.begin());
            unsigned local_node_index = p_cell_population->GetElement(first_element)->GetNodeLocalIndex(node_index);
            MonolayerVertexElementType node_type = p_cell_population->GetElement(first_element)->GetNodeType(local_node_index);

            if (node_type != MonolayerVertexElementType::Basal)
            {
                continue;
            }
            if (inner_prod(node_location - mPointOnSubstrate, mNormalToSubstrate) > 0.0)
            {
                condition_satisfied = false;
                break;
            }
        }
    }

    return condition_satisfied;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SubstrateBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PointOnSubstrate>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mPointOnSubstrate[index] << ",";
    }
    *rParamsFile << mPointOnSubstrate[SPACE_DIM - 1] << "</PointOnSubstrate>\n";

    *rParamsFile << "\t\t\t<NormalToSubstrate>";
    for (unsigned index = 0; index != SPACE_DIM - 1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mNormalToSubstrate[index] << ",";
    }
    *rParamsFile << mNormalToSubstrate[SPACE_DIM - 1] << "</NormalToSubstrate>\n";
    *rParamsFile << "\t\t\t<UseJiggledNodesOnSubstrate>" << mUseJiggledNodesOnSubstrate << "</UseJiggledNodesOnSubstrate>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class SubstrateBoundaryCondition<1, 1>;
template class SubstrateBoundaryCondition<1, 2>;
template class SubstrateBoundaryCondition<2, 2>;
template class SubstrateBoundaryCondition<1, 3>;
template class SubstrateBoundaryCondition<2, 3>;
template class SubstrateBoundaryCondition<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SubstrateBoundaryCondition)
