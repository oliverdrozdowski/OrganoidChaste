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

#include "DecayingRandomForce.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
DecayingRandomForce<DIM>::DecayingRandomForce()
        : AbstractForce<DIM>(),
          mRelativeStrengthParameter(0.01),
          mDecayTime(1.0)
{
}

template <unsigned DIM>
DecayingRandomForce<DIM>::~DecayingRandomForce()
{
}

template <unsigned DIM>
void DecayingRandomForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    ///\todo: check whether this line influences profiling tests - if so, we should remove it.
    if (dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("DecayingRandomForce is to be used with a 3D MonolayerVertexBasedCellPopulation only");
    }

    // Define some helper variables
    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    // MonolayerVertexMesh<DIM,DIM>& r_mesh = p_cell_population->rGetMesh();
    unsigned num_nodes = p_cell_population->GetNumNodes();

    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();

    double relative_strength = mRelativeStrengthParameter * exp(-current_time / mDecayTime);

    RandomNumberGenerator* p_random = RandomNumberGenerator::Instance();

    // Iterate over vertices in the cell population
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
        // c_vector<double, DIM>& applied_force_vector = p_this_node->rGetAppliedForce();
        c_vector<double, DIM> force_on_node = zero_vector<double>(DIM);

        for (unsigned i_dim = 0; i_dim < DIM; i_dim++)
        {
            // Add Gaussian noise to force with sigma = rel_strength*force
            double std_rand = p_random->StandardNormalRandomDeviate();
            // std_rand *= relative_strength * applied_force_vector[i_dim];
            std_rand *= relative_strength;
            force_on_node[i_dim] += std_rand;
        }
        p_this_node->AddAppliedForceContribution(force_on_node);
    }
}

template <unsigned DIM>
double DecayingRandomForce<DIM>::GetRelativeStrengthParameter()
{
    return mRelativeStrengthParameter;
}

template <unsigned DIM>
double DecayingRandomForce<DIM>::GetDecayTimeParameter()
{
    return mDecayTime;
}

template <unsigned DIM>
void DecayingRandomForce<DIM>::SetRelativeStrengthParameter(double relativeStrengthParameter)
{
    mRelativeStrengthParameter = relativeStrengthParameter;
}

template <unsigned DIM>
void DecayingRandomForce<DIM>::SetDecayTimeParameter(double decayTimeParameter)
{
    mDecayTime = decayTimeParameter;
}

template <unsigned DIM>
void DecayingRandomForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DecayTimeParameter>" << mDecayTime << "</DecayTimeParameter>\n";
    *rParamsFile << "\t\t\t<RelativeStrengthParameter>" << mRelativeStrengthParameter << "</RelativeStrengthParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class DecayingRandomForce<1>;
template class DecayingRandomForce<2>;
template class DecayingRandomForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"