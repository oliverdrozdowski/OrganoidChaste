/*

Copyright (c) 2005-2023, University of Oxford.
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

#include "PopulationEdgeLengthWriter.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PopulationEdgeLengthWriter<ELEMENT_DIM, SPACE_DIM>::PopulationEdgeLengthWriter()
        : AbstractMonolayerVertexPopulationWriter<ELEMENT_DIM, SPACE_DIM>("cellpopulationmidplanelength.dat") {}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationEdgeLengthWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MonolayerVertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("PopulationEdgeLengthWriter is to be used with a 3D MonolayerVertexBasedCellPopulation only");
}

template <>
void PopulationEdgeLengthWriter<3, 3>::Visit(MonolayerVertexBasedCellPopulation<3>* pCellPopulation)
{

    MutableMonolayerVertexMesh<3, 3>& r_mesh = pCellPopulation->rGetMesh();

    // Loop over lateral faces
    // Calculated mid edge length for every lateral face
    for (auto face_iter = r_mesh.GetFaceIteratorBegin(); face_iter != r_mesh.GetFaceIteratorEnd(); ++face_iter)
    {
        // Only check lateral faces
        if (face_iter->GetFaceType() != MonolayerVertexElementType::Lateral)
        {
            continue;
        }

        // We approximated the midplane edge length by the mean of the basal and apical edge lengths
        double apical_length{};
        double basal_length{};

        unsigned num_nodes = face_iter->GetNumNodes();
        assert(num_nodes == 4);

        for (unsigned local_index = 0; local_index < num_nodes; local_index++)
        {
            Node<3>* p_current_node = face_iter->GetNode(local_index);
            Node<3>* p_next_node = face_iter->GetNode((local_index + 1) % num_nodes);

            MonolayerVertexElementType current_node_type = face_iter->GetNodeType(local_index);
            MonolayerVertexElementType next_node_type = face_iter->GetNodeType((local_index + 1) % num_nodes);

            // Name nodes correctly
            if (current_node_type == MonolayerVertexElementType::Apical && next_node_type == MonolayerVertexElementType::Apical)
            {
                // we found the apical edge
                apical_length = r_mesh.GetDistanceBetweenNodes(p_current_node->GetIndex(), p_next_node->GetIndex());
            }
            else if (current_node_type == MonolayerVertexElementType::Basal && next_node_type == MonolayerVertexElementType::Basal)
            {
                // we found the basal edge
                basal_length = r_mesh.GetDistanceBetweenNodes(p_current_node->GetIndex(), p_next_node->GetIndex());
            }
        }
        double midplane_length = (apical_length + basal_length) / 2.0;

        *this->mpOutStream << midplane_length << " ";
    }
}

// Explicit instantiation
template class PopulationEdgeLengthWriter<1, 1>;
template class PopulationEdgeLengthWriter<1, 2>;
template class PopulationEdgeLengthWriter<2, 2>;
template class PopulationEdgeLengthWriter<1, 3>;
template class PopulationEdgeLengthWriter<2, 3>;
template class PopulationEdgeLengthWriter<3, 3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PopulationEdgeLengthWriter)
