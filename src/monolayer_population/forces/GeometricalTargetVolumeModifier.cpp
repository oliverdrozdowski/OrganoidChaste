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

#include "GeometricalTargetVolumeModifier.hpp"
#include <cmath>
#include "AbstractPhaseBasedCellCycleModel.hpp"
#include "ApoptoticCellProperty.hpp"

template <unsigned DIM>
GeometricalTargetVolumeModifier<DIM>::GeometricalTargetVolumeModifier(MonolayerVertexBasedCellPopulation<DIM>* pPopulation)
        : AbstractTargetVolumeModifier<DIM>(),
          mGrowthDuration(DOUBLE_UNSET),
          mT1AdaptationDuration(0.0),
          mpPopulation(pPopulation)
{
    mpPopulation->SetSaveT1TimeAndVolumeToCells(true);
}

template <unsigned DIM>
GeometricalTargetVolumeModifier<DIM>::GeometricalTargetVolumeModifier()
        : AbstractTargetVolumeModifier<DIM>(),
          mGrowthDuration(DOUBLE_UNSET),
          mT1AdaptationDuration(0.0),
          mpPopulation(nullptr)
{
}

template <unsigned DIM>
GeometricalTargetVolumeModifier<DIM>::~GeometricalTargetVolumeModifier()
{
}

template <unsigned DIM>
void GeometricalTargetVolumeModifier<DIM>::UpdateTargetVolumeOfCell(CellPtr pCell)
{
    // Get target area A of a healthy cell in S, G2 or M phase
    double cell_target_area = this->mReferenceTargetVolume;

    // We take the growth duration that is slowest
    double growth_duration = mGrowthDuration;
    if (growth_duration == DOUBLE_UNSET)
    {
        if (dynamic_cast<AbstractPhaseBasedCellCycleModel*>(pCell->GetCellCycleModel()) == nullptr)
        {
            EXCEPTION("If SetGrowthDuration() has not been called, a subclass of AbstractPhaseBasedCellCycleModel must be used");
        }
        AbstractPhaseBasedCellCycleModel* p_model = static_cast<AbstractPhaseBasedCellCycleModel*>(pCell->GetCellCycleModel());

        growth_duration = p_model->GetG1Duration();

        // If the cell is differentiated then its G1 duration is infinite
        if (growth_duration == DBL_MAX)
        {
            // This is just for fixed cell-cycle models, need to work out how to find the g1 duration
            growth_duration = p_model->GetTransitCellG1Duration();
        }
    }

    if (pCell->HasCellProperty<ApoptoticCellProperty>())
    {
        // Age of cell when apoptosis begins
        if (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime() < growth_duration)
        {
            cell_target_area *= 0.5 * (1 + (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime()) / growth_duration);
        }

        // The target area of an apoptotic cell decreases linearly to zero
        double time_spent_apoptotic = SimulationTime::Instance()->GetTime() - pCell->GetStartOfApoptosisTime();

        cell_target_area *= 1.0 - 0.5 / (pCell->GetApoptosisTime()) * time_spent_apoptotic;
        if (cell_target_area < 0)
        {
            cell_target_area = 0;
        }
    }
    else
    {
        double cell_age = pCell->GetAge();

        // The target area of a proliferating cell increases linearly from A/2 to A over the course of the prescribed duration
        if (cell_age < growth_duration)
        {
            cell_target_area *= 0.5 * (1 + cell_age / growth_duration);
        }
        else
        {
            /**
             * At division, daughter cells inherit the cell data array from the mother cell.
             * Here, we assign the target area that we want daughter cells to have to cells
             * that we know to divide in this time step.
             *
             * \todo This is a little hack that we might want to clean up in the future.
             */
            if (pCell->ReadyToDivide())
            {
                cell_target_area = 0.5 * this->mReferenceTargetVolume;
            }
        }
    }

    // If we have had T1 transitions/protorosette formation and no instantaneous T1 relaxation,
    // we also need to calculate the cell target volume for this case
    double cell_target_volume_t1 = this->mReferenceTargetVolume;
    std::vector<std::string> vector_of_cell_data_keys = pCell->GetCellData()->GetKeys();
    if ((std::find(vector_of_cell_data_keys.begin(), vector_of_cell_data_keys.end(), "time last T1")
         != vector_of_cell_data_keys.end())
        && mT1AdaptationDuration > 0.0)
    {
        double time_last_t1 = pCell->GetCellData()->GetItem("time last T1");
        double volume_at_t1 = pCell->GetCellData()->GetItem("volume last T1");

        double time_since_t1 = SimulationTime::Instance()->GetTime() - time_last_t1;
        if (time_since_t1 < mT1AdaptationDuration)
        {
            // We linearly interpolate between the volume after T1 and the reference volume
            cell_target_volume_t1 = volume_at_t1 + (this->mReferenceTargetVolume - volume_at_t1) * (time_since_t1) / mT1AdaptationDuration;
        }
    }

    // Since we want both effects, we always take the slower volume adaptation
    double diff_cell_target_area = std::abs(cell_target_area - this->mReferenceTargetVolume);
    double diff_cell_target_volume_t1 = std::abs(cell_target_volume_t1 - this->mReferenceTargetVolume);
    double final_cell_target_volume = diff_cell_target_area > diff_cell_target_volume_t1 ? cell_target_area : cell_target_volume_t1;

    // Set cell data
    pCell->GetCellData()->SetItem("target volume", final_cell_target_volume);
}

template <unsigned DIM>
double GeometricalTargetVolumeModifier<DIM>::GetGrowthDuration()
{
    return mGrowthDuration;
}

template <unsigned DIM>
void GeometricalTargetVolumeModifier<DIM>::SetGrowthDuration(double growthDuration)
{
    assert(growthDuration >= 0.0);
    mGrowthDuration = growthDuration;
}

template <unsigned DIM>
double GeometricalTargetVolumeModifier<DIM>::GetT1AdaptationDuration()
{
    return mT1AdaptationDuration;
}

template <unsigned DIM>
void GeometricalTargetVolumeModifier<DIM>::SetT1AdaptationDuration(double t1AdaptationDuration)
{
    assert(t1AdaptationDuration >= 0.0);
    mT1AdaptationDuration = t1AdaptationDuration;
}

template <unsigned DIM>
void GeometricalTargetVolumeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<GrowthDuration>" << mGrowthDuration << "</GrowthDuration>\n";
    *rParamsFile << "\t\t\t<T1AdaptationDuration>" << mT1AdaptationDuration << "</T1AdaptationDuration>\n";

    // Next, call method on direct parent class
    AbstractTargetVolumeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class GeometricalTargetVolumeModifier<1>;
template class GeometricalTargetVolumeModifier<2>;
template class GeometricalTargetVolumeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GeometricalTargetVolumeModifier)
