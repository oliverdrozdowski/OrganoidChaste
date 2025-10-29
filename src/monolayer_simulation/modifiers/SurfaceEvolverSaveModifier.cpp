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

#include "SurfaceEvolverSaveModifier.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"

template <unsigned DIM>
SurfaceEvolverSaveModifier<DIM>::SurfaceEvolverSaveModifier()
        : AbstractCellBasedSimulationModifier<DIM>(),
          mOutputDirectory("")
{
}

template <unsigned DIM>
SurfaceEvolverSaveModifier<DIM>::~SurfaceEvolverSaveModifier()
{
}

template <unsigned DIM>
MonolayerVertexMeshSurfaceEvolverWriter<DIM, DIM>* SurfaceEvolverSaveModifier<DIM>::CreateSurfaceEvolverWriter(std::string filename)
{
    MonolayerVertexMeshSurfaceEvolverWriter<DIM, DIM>* p_writer = new MonolayerVertexMeshSurfaceEvolverWriter<DIM, DIM>(mOutputDirectory, filename, false);
    if (mpSurfaceTensionSubForce)
    {
        p_writer->SetSurfaceTensionSubForce(mpSurfaceTensionSubForce);
    }
    if (mUseRandomizedVolumes)
    {
        p_writer->SetUseRandomizedVolumes(mUseRandomizedVolumes, mUniformDistributionBoundary);
    }
    if (!mMapTensionToColor.empty())
    {
        p_writer->SetMapTensionToColor(mMapTensionToColor);
    }
    if (mWriteFaceTypeIntoFile)
    {
        p_writer->SetWriteFaceTypeIntoFile(true);
    }
    if (mFixBoundaryNodes)
    {
        p_writer->SetFixBoundaryNodes(true);
    }
    if (mConsiderLumenAsCell)
    {
        p_writer->SetConsiderLumenAsCell(true);
    }
    return p_writer;
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::SetSaveAtEnd(bool saveAtEnd)
{
    mSaveAtEnd = saveAtEnd;
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::SetSaveEveryTime(bool saveEveryTime)
{
    mSaveEveryTime = saveEveryTime;
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::SetSurfaceTensionSubForce(boost::shared_ptr<SurfaceTensionSubForce<DIM> > pForce)
{
    mpSurfaceTensionSubForce = pForce;
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::SetMapTensionToColor(std::map<double, int> map)
{
    mMapTensionToColor = map;
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::SetFixBoundaryNodes(bool fixBoundaries)
{
    mFixBoundaryNodes = fixBoundaries;
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::SetWriteFaceTypeIntoFile(bool writeFaceType)
{
    mWriteFaceTypeIntoFile = writeFaceType;
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::SetConsiderLumenAsCell(bool considerLumenAsCell)
{
    mConsiderLumenAsCell = considerLumenAsCell;
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::SetUseRandomizedVolumes(bool useRandomVol, double distributionBoundary)
{
    mUseRandomizedVolumes = useRandomVol;
    mUniformDistributionBoundary = distributionBoundary;
}

template <unsigned DIM>
bool SurfaceEvolverSaveModifier<DIM>::GetUseRandomizedVolumes()
{
    return mUseRandomizedVolumes;
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    if (mSaveEveryTime)
    {
        if (bool(dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation)))
        {
            MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&(rCellPopulation));
            // MonolayerVertexMesh<DIM,DIM> rMesh = static_cast<MonolayerVertexMesh<DIM,DIM>>(rCellPopulation.rGetMesh());

            unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
            std::stringstream filename;
            filename << "results_se_" << num_timesteps;

            MonolayerVertexMeshSurfaceEvolverWriter<DIM, DIM>* surfaceEvolverWriter = CreateSurfaceEvolverWriter(filename.str());
            surfaceEvolverWriter->WriteFilesUsingMesh(p_cell_population->rGetMesh());
            delete surfaceEvolverWriter;
        }
        else
        {
            EXCEPTION("SurfaceEvolverSaveModifier only works with MonolayerVertexBasedCellPopulation.");
        }
    }
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    if (bool(dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation)))
    {
        MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&(rCellPopulation));
        // MonolayerVertexMesh<DIM,DIM> rMesh = static_cast<MonolayerVertexMesh<DIM,DIM>>(rCellPopulation.rGetMesh());
        MonolayerVertexMeshSurfaceEvolverWriter<DIM, DIM>* surfaceEvolverWriter = CreateSurfaceEvolverWriter("results_se_fin");
        surfaceEvolverWriter->WriteFilesUsingMesh(p_cell_population->rGetMesh());
        delete surfaceEvolverWriter;
    }
    else
    {
        EXCEPTION("SurfaceEvolverSaveModifier only works with MonolayerVertexBasedCellPopulation.");
    }
}

template <unsigned DIM>
void SurfaceEvolverSaveModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class SurfaceEvolverSaveModifier<1>;
template class SurfaceEvolverSaveModifier<2>;
template class SurfaceEvolverSaveModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SurfaceEvolverSaveModifier)
