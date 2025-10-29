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

#include "FiniteThicknessSimulation3d.hpp"
#include "GeometricalGrowthTargetVolumeModifier.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "SmallAreaDampingModifier.hpp"
#include "SurfaceTensionForce.hpp"

FiniteThicknessSimulation3d::FiniteThicknessSimulation3d(AbstractCellPopulation<3>& rCellPopulation,
                                                         bool deleteCellPopulationInDestructor,
                                                         bool initialiseCells)
        : OffLatticeSimulation<3>(rCellPopulation,
                                  deleteCellPopulationInDestructor,
                                  initialiseCells)
{
    /*
     * Throw an exception message if not using a  MonolayerVertexBasedCellPopulation.
     */
    if (dynamic_cast<MonolayerVertexBasedCellPopulation<3>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("FiniteThicknessSimulation3d is to be used with MonolayerVertexBasedCellPopulation (or subclasses) only");
    }
}

FiniteThicknessSimulation3d::~FiniteThicknessSimulation3d()
{
}

void FiniteThicknessSimulation3d::UpdateCellLocationsAndTopology()
{
    // First call UpdateTargetVolumeAfterReMesh()
    this->UpdateTargetVolumeAfterReMesh();

    // Then call method on base class
    OffLatticeSimulation<3>::UpdateCellLocationsAndTopology();
}

void FiniteThicknessSimulation3d::UpdateCellPopulation()
{
    unsigned deaths_before = this->mNumDeaths;
    unsigned births_before = this->mNumBirths;

    // First call method on base class
    OffLatticeSimulation<3>::UpdateCellPopulation();

    // Then call UpdateModifiersAfterBirthsAndDeaths()
    bool birthsOrDeathsOccured = ((this->mNumDeaths > deaths_before) || (this->mNumBirths > births_before));
    this->UpdateModifiersAfterBirthsAndDeaths(birthsOrDeathsOccured);
}

void FiniteThicknessSimulation3d::UpdateModifiersAfterBirthsAndDeaths(bool birthsOrDeathsOccured)
{
    // Go through each force
    for (typename std::vector<boost::shared_ptr<AbstractForce<3, 3> > >::iterator iter = this->mForceCollection.begin();
         iter != this->mForceCollection.end();
         ++iter)
    {
        // If this force is a surface tension force, we re-set the surface tensions
        if (dynamic_cast<SurfaceTensionForce<3>*>(&(*(*iter))) != nullptr)
        {
            SurfaceTensionForce<3>* p_force = static_cast<SurfaceTensionForce<3>*>(&(*(*iter)));
            // double apical_surface_tension = p_force->GetApicalStandardTension();
            // double basal_surface_tension = p_force->GetBasalStandardTension();
            // double lateral_surface_tension = p_force->GetLateralStandardTension();

            if (true)
            {
                MonolayerVertexBasedCellPopulation<3>* p_population = static_cast<MonolayerVertexBasedCellPopulation<3>*>(&(mrCellPopulation));

                // Set the surface tensions
                // p_force->CreateSurfaceTensionParametersForCells(apical_surface_tension, basal_surface_tension,
                //																								lateral_surface_tension, &(p_population->rGetMesh()) );
                p_force->UpdateSurfaceTensions(p_population);
            }
        }
        // On the other hand, if we have a GeneralizedVolumeConservingForce, we must check for
        // a surface tension subforce and update
        else if (dynamic_cast<GeneralizedVolumeConservingForce<3>*>(&(*(*iter))) != nullptr)
        {
            GeneralizedVolumeConservingForce<3>* p_force = static_cast<GeneralizedVolumeConservingForce<3>*>(&(*(*iter)));
            boost::shared_ptr<SurfaceTensionSubForce<3> > p_subforce = p_force->GetSurfaceTensionSubForce();
            if (p_subforce)
            {
                MonolayerVertexBasedCellPopulation<3>* p_population = static_cast<MonolayerVertexBasedCellPopulation<3>*>(&(mrCellPopulation));
                p_subforce->UpdateSurfaceTensions(p_population);
            }
        }
    }
    // this->UpdateTargetVolumeAfterReMesh();
}

void FiniteThicknessSimulation3d::UpdateTargetVolumeAfterReMesh()
{
    // Go through each modifier
    for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<3, 3> > >::iterator iter = this->mSimulationModifiers.begin();
         iter != this->mSimulationModifiers.end();
         ++iter)
    {
        // If this modifier is indeed a GeometricalTargetVolumeModifier we update the target volumes
        // This should be generalized for other AbstractTargetVolumeModifiers at some point
        // And maybe with a morge general procedure than just update target volumes...
        if (dynamic_cast<GeometricalTargetVolumeModifier<3>*>(&(*(*iter))) != nullptr)
        {
            GeometricalTargetVolumeModifier<3>* p_modifier = static_cast<GeometricalTargetVolumeModifier<3>*>(&(*(*iter)));
            p_modifier->UpdateTargetVolumes(this->mrCellPopulation);
        }
        // Also for GeometricalGrowthTargetVolumeModifier
        if (dynamic_cast<GeometricalGrowthTargetVolumeModifier<3>*>(&(*(*iter))) != nullptr)
        {
            GeometricalGrowthTargetVolumeModifier<3>* p_modifier = static_cast<GeometricalGrowthTargetVolumeModifier<3>*>(&(*(*iter)));
            p_modifier->UpdateTargetVolumes(this->mrCellPopulation);
        }
        // Same thing with the SmallAreaDampingModifier
        if (dynamic_cast<SmallAreaDampingModifier<3>*>(&(*(*iter))) != nullptr)
        {
            SmallAreaDampingModifier<3>* p_modifier = static_cast<SmallAreaDampingModifier<3>*>(&(*(*iter)));
            p_modifier->UpdateAtEndOfTimeStep(this->mrCellPopulation);
        }
    }
}

void FiniteThicknessSimulation3d::SetupSolve()
{
    // Go through each modifier
    for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<3, 3> > >::iterator iter = this->mSimulationModifiers.begin();
         iter != this->mSimulationModifiers.end();
         ++iter)
    {
        // If this modifier is indeed a GeometricalTargetVolumeModifier we initialize the cell data
        if (dynamic_cast<GeometricalTargetVolumeModifier<3>*>(&(*(*iter))) != nullptr)
        {
            GeometricalTargetVolumeModifier<3>* p_modifier = static_cast<GeometricalTargetVolumeModifier<3>*>(&(*(*iter)));
            MonolayerVertexBasedCellPopulation<3>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<3>*>(&(this->mrCellPopulation));
            // initialize to target volume and negative time, such that we stay at target volume.
            p_cell_population->InitializeT1TimeAndVolumeToCells(-10.0, p_modifier->GetReferenceTargetVolume());
        }
        // Same for GeometricalGrowthTargetVolumeModifier
        if (dynamic_cast<GeometricalGrowthTargetVolumeModifier<3>*>(&(*(*iter))) != nullptr)
        {
            GeometricalGrowthTargetVolumeModifier<3>* p_modifier = static_cast<GeometricalGrowthTargetVolumeModifier<3>*>(&(*(*iter)));
            MonolayerVertexBasedCellPopulation<3>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<3>*>(&(this->mrCellPopulation));
            // initialize to target volume and negative time, such that we stay at target volume.
            p_cell_population->InitializeT1TimeAndVolumeToCells(-10.0, p_modifier->GetReferenceTargetVolume());
        }
    }

    // Then call method on base class
    OffLatticeSimulation<3>::SetupSolve();
}

void FiniteThicknessSimulation3d::OutputSimulationParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    OffLatticeSimulation<3>::OutputSimulationParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FiniteThicknessSimulation3d)