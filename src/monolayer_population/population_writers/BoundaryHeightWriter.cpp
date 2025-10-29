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

#include "BoundaryHeightWriter.hpp"
#include "AbstractMovingBoundary.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundaryHeightWriter<ELEMENT_DIM, SPACE_DIM>::BoundaryHeightWriter()
        : AbstractMonolayerVertexPopulationWriter<ELEMENT_DIM, SPACE_DIM>(
              "movingboundaryheight.dat") {}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BoundaryHeightWriter<ELEMENT_DIM, SPACE_DIM>::Visit(
    MonolayerVertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE

    // double max_z_position = pCellPopulation->GetMaxZPosition();

    std::vector<boost::shared_ptr<AbstractMovingBoundary<SPACE_DIM> > > boundaries = pCellPopulation->GetMovingBoundaries();

    for (auto boundary : boundaries)
    {
        *this->mpOutStream << boundary->GetZPosition() << " ";
    }
}

// Explicit instantiation
template class BoundaryHeightWriter<1, 1>;
template class BoundaryHeightWriter<1, 2>;
template class BoundaryHeightWriter<2, 2>;
template class BoundaryHeightWriter<1, 3>;
template class BoundaryHeightWriter<2, 3>;
template class BoundaryHeightWriter<3, 3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(BoundaryHeightWriter)
