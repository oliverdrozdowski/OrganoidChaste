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

#include "GeometricalGrowthTargetVolumeModifier.hpp"
#include <cmath>
#include "AbstractPhaseBasedCellCycleModel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

template <unsigned DIM>
GeometricalGrowthTargetVolumeModifier<DIM>::GeometricalGrowthTargetVolumeModifier(MonolayerVertexBasedCellPopulation<DIM>* pPopulation)
        : AbstractTargetVolumeModifier<DIM>(),
          mAgeToStartGrowing(DOUBLE_UNSET),
          mGrowthRate(DOUBLE_UNSET),
          mT1AdaptationDuration(0.0),
          mpPopulation(pPopulation)
{
    mpPopulation->SetSaveT1TimeAndVolumeToCells(true);
}

template <unsigned DIM>
GeometricalGrowthTargetVolumeModifier<DIM>::GeometricalGrowthTargetVolumeModifier()
        : AbstractTargetVolumeModifier<DIM>(),
          mAgeToStartGrowing(DOUBLE_UNSET),
          mGrowthRate(DOUBLE_UNSET),
          mT1AdaptationDuration(0.0),
          mpPopulation(nullptr)
{
}

template <unsigned DIM>
GeometricalGrowthTargetVolumeModifier<DIM>::~GeometricalGrowthTargetVolumeModifier()
{
}

template <unsigned DIM>
void GeometricalGrowthTargetVolumeModifier<DIM>::UpdateTargetVolumeOfCell(CellPtr pCell)
{
    // Get target area A of a healthy cell in S, G2 or M phase
    double cell_target_volume = this->mReferenceTargetVolume;
    bool use_double_reference_volume_t1 = false;

    // The target area of an apoptotic cell decreases linearly to zero
    if (pCell->HasCellProperty<ApoptoticCellProperty>())
    {
        ///\todo: which cells are apoptotic? if they get apoptotic during G2-phase then this line has to be changed
        double time_spent_apoptotic = SimulationTime::Instance()->GetTime() - pCell->GetStartOfApoptosisTime();
        cell_target_volume = cell_target_volume - 0.5 * cell_target_volume / (pCell->GetApoptosisTime()) * time_spent_apoptotic;

        // Don't allow a negative target area
        if (cell_target_volume < 0)
        {
            cell_target_volume = 0;
        }
    }

    else
    {
        bool cell_is_differentiated = pCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>();
        if (!cell_is_differentiated)
        {
            /*
             * If not using a phase-based cell-cycle model, the target area of a proliferating cell increases
             * linearly from mReferenceTargetVolume as soon as the cell's age exceeds mAgeToStartGrowing, with
             * growth rate mGrowthRate.
             */
            double growth_start_time = mAgeToStartGrowing;
            double growth_rate = mGrowthRate;
            if (growth_start_time == DOUBLE_UNSET)
            {
                if (dynamic_cast<AbstractPhaseBasedCellCycleModel*>(pCell->GetCellCycleModel()) == nullptr)
                {
                    EXCEPTION("If SetAgeToStartGrowing() has not been called, a subclass of AbstractPhaseBasedCellCycleModel must be used");
                }

                /*
                 * If using a phase-based cell-cycle model, the target area of a proliferating cell increases
                 * linearly from mReferenceTargetVolume to 2*mReferenceTargetVolume during the cell's G2 phase.
                 */
                AbstractPhaseBasedCellCycleModel* p_model = static_cast<AbstractPhaseBasedCellCycleModel*>(pCell->GetCellCycleModel());
                growth_start_time = p_model->GetMDuration() + p_model->GetG1Duration() + p_model->GetSDuration();

                double g2_duration = p_model->GetG2Duration();
                growth_rate = cell_target_volume / g2_duration;
            }
            else if (mGrowthRate == DOUBLE_UNSET)
            {
                EXCEPTION("If SetAgeToStartGrowing() has been called, then SetGrowthRate() must also be called");
            }

            double time_spent_growing = pCell->GetAge() - growth_start_time;

            // The target area of a proliferating cell increases linearly
            if (time_spent_growing > 0)
            {
                cell_target_volume += time_spent_growing * growth_rate;
                use_double_reference_volume_t1 = true;
            }

            /**
             * At division, daughter cells inherit the cell data array from the mother cell.
             * Here, we assign the target area that we want daughter cells to have to cells
             * that we know to divide in this time step.
             *
             * \todo This is a little hack that we might want to clean up in the future.
             */
            if (pCell->ReadyToDivide())
            {
                cell_target_volume = this->mReferenceTargetVolume;
            }
        }
    }

    // If we have had T1 transitions/protorosette formation and no instantaneous T1 relaxation,
    // we also need to calculate the cell target volume for this case
    double t1_target_factor = use_double_reference_volume_t1 ? 2.0 : 1.0;
    double cell_target_volume_t1 = t1_target_factor * this->mReferenceTargetVolume;
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
            cell_target_volume_t1 = volume_at_t1 + (t1_target_factor * this->mReferenceTargetVolume - volume_at_t1) * (time_since_t1) / mT1AdaptationDuration;
        }
    }

    // Since we want both effects, we always take the slower volume adaptation
    double diff_cell_target_volume = std::abs(cell_target_volume - t1_target_factor * this->mReferenceTargetVolume);
    double diff_cell_target_volume_t1 = std::abs(cell_target_volume_t1 - t1_target_factor * this->mReferenceTargetVolume);
    double final_cell_target_volume = diff_cell_target_volume > diff_cell_target_volume_t1 ? cell_target_volume : cell_target_volume_t1;

    // Set cell data
    pCell->GetCellData()->SetItem("target volume", final_cell_target_volume);
}

template <unsigned DIM>
double GeometricalGrowthTargetVolumeModifier<DIM>::GetAgeToStartGrowing()
{
    return mAgeToStartGrowing;
}

template <unsigned DIM>
void GeometricalGrowthTargetVolumeModifier<DIM>::SetAgeToStartGrowing(double ageToStartGrowing)
{
    assert(ageToStartGrowing >= 0.0);
    mAgeToStartGrowing = ageToStartGrowing;
}

template <unsigned DIM>
double GeometricalGrowthTargetVolumeModifier<DIM>::GetGrowthRate()
{
    return mGrowthRate;
}

template <unsigned DIM>
void GeometricalGrowthTargetVolumeModifier<DIM>::SetGrowthRate(double growthRate)
{
    assert(growthRate >= 0.0);
    mGrowthRate = growthRate;
}

template <unsigned DIM>
double GeometricalGrowthTargetVolumeModifier<DIM>::GetT1AdaptationDuration()
{
    return mT1AdaptationDuration;
}

template <unsigned DIM>
void GeometricalGrowthTargetVolumeModifier<DIM>::SetT1AdaptationDuration(double t1AdaptationDuration)
{
    assert(t1AdaptationDuration >= 0.0);
    mT1AdaptationDuration = t1AdaptationDuration;
}

template <unsigned DIM>
void GeometricalGrowthTargetVolumeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AgeToStartGrowing>" << mAgeToStartGrowing << "</AgeToStartGrowing>\n";
    *rParamsFile << "\t\t\t<GrowthRate>" << mGrowthRate << "</GrowthRate>\n";
    *rParamsFile << "\t\t\t<T1AdaptationDuration>" << mT1AdaptationDuration << "</T1AdaptationDuration>\n";

    // Next, call method on direct parent class
    AbstractTargetVolumeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class GeometricalGrowthTargetVolumeModifier<1>;
template class GeometricalGrowthTargetVolumeModifier<2>;
template class GeometricalGrowthTargetVolumeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GeometricalGrowthTargetVolumeModifier)
