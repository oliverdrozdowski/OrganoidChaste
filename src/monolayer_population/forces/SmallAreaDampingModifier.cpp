
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

#include "SmallAreaDampingModifier.hpp"

template <unsigned DIM>
SmallAreaDampingModifier<DIM>::SmallAreaDampingModifier()
        : AbstractCellBasedSimulationModifier<DIM>(),
          mSmallAreaThreshold(0.0),
          mForceDampingRatio(1.0),
          mVolumeDampingRatio(1.0),
          mpPopulation(nullptr),
          mpForce(nullptr)
{
}

template <unsigned DIM>
SmallAreaDampingModifier<DIM>::SmallAreaDampingModifier(MonolayerVertexBasedCellPopulation<DIM>* pPopulation,
                                                        AbstractForce<DIM>* pForce)
        : AbstractCellBasedSimulationModifier<DIM>(),
          mSmallAreaThreshold(0.0),
          mForceDampingRatio(1.0),
          mVolumeDampingRatio(1.0),
          mpPopulation(pPopulation),
          mpForce(pForce)
{
    if (dynamic_cast<SurfaceTensionForce<DIM>*>(pForce) == nullptr && dynamic_cast<GeneralizedVolumeConservingForce<DIM>*>(pForce) == nullptr)
    {
        EXCEPTION("SmallAreaDampingModifier is only implemented with SurfaceTensionForce and GeneralizedVolumeConservingForce for now.");
    }
    if (dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(pPopulation) == nullptr)
    {
        EXCEPTION("SmallAreaDampingModifier is to be used with a 3D MonolayerVertexBasedCellPopulation only");
    }
}

template <unsigned DIM>
SmallAreaDampingModifier<DIM>::~SmallAreaDampingModifier()
{
}

template <unsigned DIM>
void SmallAreaDampingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    if (dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("SmallAreaDampingModifier is to be used with a 3D MonolayerVertexBasedCellPopulation only");
    }
    else
    {
        MonolayerVertexBasedCellPopulation<DIM>* p_cellPopulation = dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        GiveMapOfDampedNodesToForce(p_cellPopulation);
    }
}

template <unsigned DIM>
void SmallAreaDampingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We already determine the damping factor in the beginning if the initial data already contains
     * relevant small faces
     */
    if (dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("SmallAreaDampingModifier is to be used with a 3D MonolayerVertexBasedCellPopulation only");
    }
    else
    {
        MonolayerVertexBasedCellPopulation<DIM>* p_cellPopulation = dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        GiveMapOfDampedNodesToForce(p_cellPopulation);
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void SmallAreaDampingModifier<1>::GiveMapOfDampedNodesToForce(MonolayerVertexBasedCellPopulation<1>* pCellPopulation)
{
    EXCEPTION("SmallAreaDampingModifier does not work in 1D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned DIM>
void SmallAreaDampingModifier<DIM>::GiveMapOfDampedNodesToForce(MonolayerVertexBasedCellPopulation<DIM>* pCellPopulation)
{
    std::map<Node<DIM>*, DampingCoefficientPair> mapNodeToDampingRate;
    MonolayerVertexMesh<DIM, DIM>& r_mesh = pCellPopulation->rGetMesh();
    // Loop over faces, because if we look at cell's faces, we go over the same faces multiple times
    unsigned num_faces = r_mesh.GetNumFaces();
    for (unsigned index_face = 0; index_face < num_faces; ++index_face)
    {
        MonolayerVertexElement<DIM - 1, DIM>* p_face = r_mesh.GetFace(index_face);
        double face_area = r_mesh.CalculateAreaOfFace(p_face);
        if (face_area < mSmallAreaThreshold)
        {
            unsigned num_nodes = p_face->GetNumNodes();
            double damping_ratio_force = mForceDampingRatio * sqrt(face_area);
            double damping_ratio_volume = mVolumeDampingRatio * sqrt(face_area);
            DampingCoefficientPair damping_ratios;
            damping_ratios.volumeCorrectionDamping = damping_ratio_volume;
            damping_ratios.forceDamping = damping_ratio_force;
            for (unsigned node_index = 0; node_index < num_nodes; ++node_index)
            {
                Node<DIM>* p_node = p_face->GetNode(node_index);
                // We assume that we only have one such small face and thus
                // just write the rate
                mapNodeToDampingRate[p_node] = damping_ratios;
            }
        }
    }
    // We implement the case of SurfaceTensionForce, which has the corresponding mapping method
    if (dynamic_cast<SurfaceTensionForce<DIM>*>(mpForce) != nullptr)
    {
        SurfaceTensionForce<DIM>* p_force = dynamic_cast<SurfaceTensionForce<DIM>*>(mpForce);
        p_force->SetMapNodeToDampingThresholds(mapNodeToDampingRate);
    }
    else if (dynamic_cast<GeneralizedVolumeConservingForce<DIM>*>(mpForce) != nullptr)
    {
        GeneralizedVolumeConservingForce<DIM>* p_force = dynamic_cast<GeneralizedVolumeConservingForce<DIM>*>(mpForce);
        p_force->SetMapNodeToDampingThresholds(mapNodeToDampingRate);
    }
}

template <unsigned DIM>
void SmallAreaDampingModifier<DIM>::SetSmallAreaThreshold(double smallAreaThreshold)
{
    assert(smallAreaThreshold >= 0.0);
    mSmallAreaThreshold = smallAreaThreshold;
}

template <unsigned DIM>
void SmallAreaDampingModifier<DIM>::SetForceDampingRatio(double forceDampingRatio)
{
    assert(forceDampingRatio > 0.0);
    mForceDampingRatio = forceDampingRatio;
}

template <unsigned DIM>
void SmallAreaDampingModifier<DIM>::SetVolumeDampingRatio(double volumeDampingRatio)
{
    assert(volumeDampingRatio > 0.0);
    mVolumeDampingRatio = volumeDampingRatio;
}

template <unsigned DIM>
double SmallAreaDampingModifier<DIM>::GetSmallAreaThreshold()
{
    return mSmallAreaThreshold;
}

template <unsigned DIM>
double SmallAreaDampingModifier<DIM>::GetForceDampingRatio()
{
    return mForceDampingRatio;
}

template <unsigned DIM>
double SmallAreaDampingModifier<DIM>::GetVolumeDampingRatio()
{
    return mVolumeDampingRatio;
}

template <unsigned DIM>
void SmallAreaDampingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<SmallAreaThreshold>" << mSmallAreaThreshold << "</SmallAreaThreshold>\n";
    *rParamsFile << "\t\t\t<VolumeDampingRatio>" << mVolumeDampingRatio << "</VolumeDampingRatio>\n";
    *rParamsFile << "\t\t\t<ForceDampingRatio>" << mForceDampingRatio << "</ForceDampingRatio>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class SmallAreaDampingModifier<1>;
template class SmallAreaDampingModifier<2>;
template class SmallAreaDampingModifier<3>;