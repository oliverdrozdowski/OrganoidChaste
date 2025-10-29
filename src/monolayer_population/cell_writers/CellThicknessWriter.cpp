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

#include "CellThicknessWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellThicknessWriter<ELEMENT_DIM, SPACE_DIM>::CellThicknessWriter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellthickness.dat")
{
    this->mVtkCellDataName = "Cell thickness";
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellThicknessWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CellThicknessWriter is to be used with a 3D MonolayerVertexBasedCellPopulation only");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellThicknessWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population)
{
    EXCEPTION("CellThicknessWriter is to be used with a 3D MonolayerVertexBasedCellPopulation only");
}

template <>
double CellThicknessWriter<3, 3>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<3, 3>* p_Cell_Population)
{
    // AbstractOffLatticeCellPopulation<3>* p_CellPopulation = dynamic_cast<AbstractOffLatticeCellPopulation<3>*>(p_Cell_Population);
    if (dynamic_cast<MonolayerVertexBasedCellPopulation<3>*>(p_Cell_Population) == nullptr)
    {
        EXCEPTION("CellThicknessWriter is to be used with a 3D MonolayerVertexBasedCellPopulation only");
    }
    MonolayerVertexBasedCellPopulation<3>* pCellPopulation = static_cast<MonolayerVertexBasedCellPopulation<3>*>(p_Cell_Population);

    double thickness = pCellPopulation->GetThicknessOfCell(pCell);
    return thickness;
}

template <>
void CellThicknessWriter<3, 3>::VisitCell(CellPtr pCell, AbstractCellPopulation<3, 3>* p_Cell_Population)
{
    // AbstractCellPopulation<3>* p_CellPopulation = static_cast<AbstractCellPopulation<3>*>(p_Cell_Population);
    if (dynamic_cast<MonolayerVertexBasedCellPopulation<3>*>(p_Cell_Population) == nullptr)
    {
        EXCEPTION("CellThicknessWriter is to be used with a 3D MonolayerVertexBasedCellPopulation only");
    }
    MonolayerVertexBasedCellPopulation<3>* pCellPopulation = static_cast<MonolayerVertexBasedCellPopulation<3>*>(p_Cell_Population);

    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    unsigned cell_id = pCell->GetCellId();
    c_vector<double, 3> centre_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    double thickness = pCellPopulation->GetThicknessOfCell(pCell);

    if (thickness < DBL_MAX) // Only write cells with finite thickness
    {
        *this->mpOutStream << location_index << " " << cell_id << " ";
        for (unsigned i = 0; i < 3; i++)
        {
            *this->mpOutStream << centre_location[i] << " ";
        }

        *this->mpOutStream << thickness << " ";
    }
}

// Explicit instantiation
template class CellThicknessWriter<1, 1>;
template class CellThicknessWriter<1, 2>;
template class CellThicknessWriter<2, 2>;
template class CellThicknessWriter<1, 3>;
template class CellThicknessWriter<2, 3>;
template class CellThicknessWriter<3, 3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellThicknessWriter)
