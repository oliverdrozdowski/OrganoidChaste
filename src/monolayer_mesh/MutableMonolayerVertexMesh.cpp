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

#include "MutableMonolayerVertexMesh.hpp"
#include "LogFile.hpp"
#include "MonolayerVertexElement.hpp"
#include "UblasCustomFunctions.hpp"
#include "Warnings.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MutableMonolayerVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                                                               std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> faces,
                                                                               std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements,
                                                                               double cellRearrangementThreshold,
                                                                               double t2Threshold,
                                                                               double cellRearrangementRatio,
                                                                               double protorosetteFormationProbability,
                                                                               double protorosetteResolutionProbabilityPerTimestep,
                                                                               double rosetteResolutionProbabilityPerTimestep)
        : mCellRearrangementThreshold(cellRearrangementThreshold),
          mCellRearrangementRatio(cellRearrangementRatio),
          mT2Threshold(t2Threshold),
          mActiveT1SwapProbability(0.0),
          mProtorosetteFormationProbability(protorosetteFormationProbability),
          mProtorosetteResolutionProbabilityPerTimestep(protorosetteResolutionProbabilityPerTimestep),
          mRosetteResolutionProbabilityPerTimestep(rosetteResolutionProbabilityPerTimestep),
          mCheckForInternalIntersections(false),
          mDistanceForT3SwapChecking(5.0),
          mLengthDependentActiveT1BoltzmannParameter(0.0),
          mIsActiveT1ProbPerEdge(false),
          mPassiveT1TransitionsCounter(0),
          mActiveT1TransitionsCounter(0)
{
    // Threshold parameters must be strictly positive
    // assert(cellRearrangementThreshold >= 0.0); // we do not assert this. Setting this to neg. allows zero-faces
    assert(t2Threshold > 0.0);
    assert(protorosetteFormationProbability >= 0.0);
    assert(protorosetteFormationProbability <= 1.0);
    assert(protorosetteResolutionProbabilityPerTimestep >= 0.0);
    assert(protorosetteResolutionProbabilityPerTimestep <= 1.0);
    assert(rosetteResolutionProbabilityPerTimestep >= 0.0);
    assert(rosetteResolutionProbabilityPerTimestep <= 1.0);

    // Reset member variables and clear mNodes and mElements
    Clear();

    // Populate mNodes and mElements an mFaces
    for (unsigned node_index = 0; node_index < nodes.size(); node_index++)
    {
        Node<SPACE_DIM>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }
    for (unsigned elem_index = 0; elem_index < vertexElements.size(); elem_index++)
    {
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = vertexElements[elem_index];
        this->mElements.push_back(p_temp_vertex_element);
    }
    for (unsigned face_index = 0; face_index < faces.size(); face_index++)
    {
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_temp_face = faces[face_index];
        this->mFaces.push_back(p_temp_face);
    }

    // Register elements with nodes
    for (unsigned index = 0; index < this->mElements.size(); index++)
    {
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = this->mElements[index];
        for (unsigned node_index = 0; node_index < p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<SPACE_DIM>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }

    this->mMeshChangesDuringSimulation = true;
    MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::UpdateElementsFacesMap();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MutableMonolayerVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                                                               std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements,
                                                                               double cellRearrangementThreshold,
                                                                               double t2Threshold,
                                                                               double cellRearrangementRatio,
                                                                               double protorosetteFormationProbability,
                                                                               double protorosetteResolutionProbabilityPerTimestep,
                                                                               double rosetteResolutionProbabilityPerTimestep)
        : mCellRearrangementThreshold(cellRearrangementThreshold),
          mCellRearrangementRatio(cellRearrangementRatio),
          mT2Threshold(t2Threshold),
          mActiveT1SwapProbability(0.0),
          mProtorosetteFormationProbability(protorosetteFormationProbability),
          mProtorosetteResolutionProbabilityPerTimestep(protorosetteResolutionProbabilityPerTimestep),
          mRosetteResolutionProbabilityPerTimestep(rosetteResolutionProbabilityPerTimestep),
          mCheckForInternalIntersections(false),
          mDistanceForT3SwapChecking(5.0),
          mLengthDependentActiveT1BoltzmannParameter(0.0),
          mIsActiveT1ProbPerEdge(false),
          mPassiveT1TransitionsCounter(0),
          mActiveT1TransitionsCounter(0)
{
    // Threshold parameters must be strictly positive
    // assert(cellRearrangementThreshold >= 0.0);
    assert(t2Threshold > 0.0);
    assert(protorosetteFormationProbability >= 0.0);
    assert(protorosetteFormationProbability <= 1.0);
    assert(protorosetteResolutionProbabilityPerTimestep >= 0.0);
    assert(protorosetteResolutionProbabilityPerTimestep <= 1.0);
    assert(rosetteResolutionProbabilityPerTimestep >= 0.0);
    assert(rosetteResolutionProbabilityPerTimestep <= 1.0);

    // Reset member variables and clear mNodes and mElements
    Clear();

    // Populate mNodes and mElements
    for (unsigned node_index = 0; node_index < nodes.size(); node_index++)
    {
        Node<SPACE_DIM>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }
    for (unsigned elem_index = 0; elem_index < vertexElements.size(); elem_index++)
    {
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = vertexElements[elem_index];
        this->mElements.push_back(p_temp_vertex_element);
    }

    // If in 3D, then also populate mFaces
    if (SPACE_DIM == 3)
    {
        // Use a std::set to keep track of which faces have been added to mFaces
        std::set<unsigned> faces_counted;

        // Loop over mElements
        for (unsigned elem_index = 0; elem_index < this->mElements.size(); elem_index++)
        {
            // Loop over faces of this element
            for (unsigned face_index = 0; face_index < this->mElements[elem_index]->GetNumFaces(); face_index++)
            {
                MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = this->mElements[elem_index]->GetFace(face_index);

                // If this face is not already contained in mFaces, then add it and update faces_counted
                if (faces_counted.find(p_face->GetIndex()) == faces_counted.end())
                {
                    this->mFaces.push_back(p_face);
                    faces_counted.insert(p_face->GetIndex());
                }
            }
        }
    }

    // Register elements with nodes
    for (unsigned index = 0; index < this->mElements.size(); index++)
    {
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = this->mElements[index];
        for (unsigned node_index = 0; node_index < p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<SPACE_DIM>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }

    this->mMeshChangesDuringSimulation = true;
    MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::UpdateElementsFacesMap();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MutableMonolayerVertexMesh()
        : mCellRearrangementThreshold(0.01),
          mCellRearrangementRatio(1.5),
          mT2Threshold(0.001),
          mActiveT1SwapProbability(0.0),
          mProtorosetteFormationProbability(0.0),
          mProtorosetteResolutionProbabilityPerTimestep(0.0),
          mRosetteResolutionProbabilityPerTimestep(0.0),
          mCheckForInternalIntersections(false),
          mDistanceForT3SwapChecking(5.0)
{
    // Note that the member variables initialised above will be overwritten as soon as archiving is complete
    this->mMeshChangesDuringSimulation = true;
    Clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::~MutableMonolayerVertexMesh()
{
    Clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCellRearrangementThreshold() const
{
    return mCellRearrangementThreshold;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetT2Threshold() const
{
    return mT2Threshold;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCellRearrangementRatio() const
{
    return mCellRearrangementRatio;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetProtorosetteFormationProbability() const
{
    return this->mProtorosetteFormationProbability;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetActiveT1SwapProbability() const
{
    return this->mActiveT1SwapProbability;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetProtorosetteResolutionProbabilityPerTimestep() const
{
    return this->mProtorosetteResolutionProbabilityPerTimestep;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetRosetteResolutionProbabilityPerTimestep() const
{
    return this->mRosetteResolutionProbabilityPerTimestep;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetDistanceForT3SwapChecking(double distanceForT3SwapChecking)
{
    mDistanceForT3SwapChecking = distanceForT3SwapChecking;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetDistanceForT3SwapChecking() const
{
    return mDistanceForT3SwapChecking;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCheckForInternalIntersections() const
{
    return mCheckForInternalIntersections;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetCellRearrangementThreshold(double cellRearrangementThreshold)
{
    mCellRearrangementThreshold = cellRearrangementThreshold;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetT2Threshold(double t2Threshold)
{
    mT2Threshold = t2Threshold;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetCellRearrangementRatio(double cellRearrangementRatio)
{
    mCellRearrangementRatio = cellRearrangementRatio;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetProtorosetteFormationProbability(double protorosetteFormationProbability)
{
    // Check that the new value is in [0,1]
    if (protorosetteFormationProbability < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }
    if (protorosetteFormationProbability > 1.0)
    {
        EXCEPTION("Attempting to assign a probability greater than one.");
    }

    // Assign the new value
    mProtorosetteFormationProbability = protorosetteFormationProbability;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetActiveT1SwapProbability(double activeT1SwapProbability, bool isActiveT1ProbabilityPerEdge)
{
    // Check that the new value is in [0,1]
    if (activeT1SwapProbability < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }

    // Assign the new value
    mActiveT1SwapProbability = activeT1SwapProbability;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetProtorosetteResolutionProbabilityPerTimestep(double protorosetteResolutionProbabilityPerTimestep)
{
    // Check that the new value is in [0,1]
    if (protorosetteResolutionProbabilityPerTimestep < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }
    if (protorosetteResolutionProbabilityPerTimestep > 1.0)
    {
        EXCEPTION("Attempting to assign a probability greater than one.");
    }

    // Assign the new value
    mProtorosetteResolutionProbabilityPerTimestep = protorosetteResolutionProbabilityPerTimestep;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetRosetteResolutionProbabilityPerTimestep(double rosetteResolutionProbabilityPerTimestep)
{
    // Check that the new value is in [0,1]
    if (rosetteResolutionProbabilityPerTimestep < 0.0)
    {
        EXCEPTION("Attempting to assign a negative probability.");
    }
    if (rosetteResolutionProbabilityPerTimestep > 1.0)
    {
        EXCEPTION("Attempting to assign a probability greater than one.");
    }

    // Assign the new value
    mRosetteResolutionProbabilityPerTimestep = rosetteResolutionProbabilityPerTimestep;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetCheckForInternalIntersections(bool checkForInternalIntersections)
{
    mCheckForInternalIntersections = checkForInternalIntersections;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    mDeletedNodeIndices.clear();
    mDeletedElementIndices.clear();

    MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size() - mDeletedNodeIndices.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return this->mElements.size() - mDeletedElementIndices.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetLocationsOfT1Swaps()
{
    return mLocationsOfT1Swaps;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetLastT2SwapLocation()
{
    return mLastT2SwapLocation;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetLocationsOfT3Swaps()
{
    return mLocationsOfT3Swaps;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElementsThatUnderwentT1Transitions()
{
    return mElementsThatUnderwentT1Transitions;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::ClearLocationsOfT1Swaps()
{
    mLocationsOfT1Swaps.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::ClearLocationsOfT3Swaps()
{
    mLocationsOfT3Swaps.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::ClearMapOfProtorosettes()
{
    mMapOfProtorosettes.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::ClearSwappedFaceIndices()
{
    mSwappedFaceIndices.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNewNode)
{
    if (mDeletedNodeIndices.empty())
    {
        pNewNode->SetIndex(this->mNodes.size());
        this->mNodes.push_back(pNewNode);
    }
    else
    {
        unsigned index = mDeletedNodeIndices.back();
        pNewNode->SetIndex(index);
        mDeletedNodeIndices.pop_back();
        delete this->mNodes[index];
        this->mNodes[index] = pNewNode;
    }
    return pNewNode->GetIndex();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::AddFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pNewFace)
{
    if (mDeletedFaceIndices.empty())
    {
        pNewFace->SetIndex(this->mFaces.size());
        this->mFaces.push_back(pNewFace);
    }
    else
    {
        unsigned index = mDeletedFaceIndices.back();
        pNewFace->SetIndex(index);
        mDeletedFaceIndices.pop_back();
        delete this->mFaces[index];
        this->mFaces[index] = pNewFace;
    }
    return pNewFace->GetIndex();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::AddElement(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pNewElement)
{
    unsigned new_element_index = pNewElement->GetIndex();

    if (new_element_index == this->mElements.size())
    {
        this->mElements.push_back(pNewElement);
    }
    else
    {
        this->mElements[new_element_index] = pNewElement;
    }
    pNewElement->RegisterWithNodes();
    return pNewElement->GetIndex();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point)
{
    this->mNodes[nodeIndex]->SetPoint(point);
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MutableMonolayerVertexMesh<1, 1>::DivideElementAlongGivenAxis(MonolayerVertexElement<1, 1>* pElement,
                                                                       c_vector<double, 1> axisOfDivision,
                                                                       bool placeOriginalElementInDirection, unsigned iterative_call_counter)
{
    EXCEPTION("DivideElementAlongGivenAxis only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MutableMonolayerVertexMesh<1, 2>::DivideElementAlongGivenAxis(MonolayerVertexElement<1, 2>* pElement,
                                                                       c_vector<double, 2> axisOfDivision,
                                                                       bool placeOriginalElementInDirection, unsigned iterative_call_counter)
{
    EXCEPTION("DivideElementAlongGivenAxis only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MutableMonolayerVertexMesh<1, 3>::DivideElementAlongGivenAxis(MonolayerVertexElement<1, 3>* pElement,
                                                                       c_vector<double, 3> axisOfDivision,
                                                                       bool placeOriginalElementInDirection, unsigned iterative_call_counter)
{
    EXCEPTION("DivideElementAlongGivenAxis only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MutableMonolayerVertexMesh<2, 2>::DivideElementAlongGivenAxis(MonolayerVertexElement<2, 2>* pElement,
                                                                       c_vector<double, 2> axisOfDivision,
                                                                       bool placeOriginalElementInDirection, unsigned iterative_call_counter)
{
    EXCEPTION("DivideElementAlongGivenAxis only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MutableMonolayerVertexMesh<2, 3>::DivideElementAlongGivenAxis(MonolayerVertexElement<2, 3>* pElement,
                                                                       c_vector<double, 3> axisOfDivision,
                                                                       bool placeOriginalElementInDirection, unsigned iterative_call_counter)
{
    EXCEPTION("DivideElementAlongGivenAxis only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElementAlongGivenAxis(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                                         c_vector<double, SPACE_DIM> axisOfDivision,
                                                                                         bool placeOriginalElementInDirection, unsigned iterative_call_counter)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE

    /* We now determine basal edges through which the plane orthogonal to axisOfDivision passes
     * end then we separate the lateral faces such that the distance of the new apical nodes to
     * the plane is minimized, but the new edge length is at least 1% the old long edge length.
     *             /o-------o\        /o--o----o\         /o--o----o\
     *            / :_ _ _ _: \      / :_ _\_ _: \       / :_ |\ _ : \
     *           o /         \ o => o /         \ o =>  o /   | |   \ o
     *           |\           /|    |\           /|     |\    | |    /|
     *            \\o-------o/ /     \\o--o----o/ /      \\o--o----o/ /
     *             \|_______| /       \|___\___| /        \|___\|__| /
     */
    unsigned element_index = pElement->GetIndex();
    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_basal_face = this->GetFaceOfType(element_index, MonolayerVertexElementType::Basal);

    // left means smaller index w.r.t. ordering in face in basal
    // in apical they are opposite nodes
    Node<SPACE_DIM>* p_basal_node_1_left = nullptr;
    Node<SPACE_DIM>* p_basal_node_1_right = nullptr;
    Node<SPACE_DIM>* p_basal_node_2_left = nullptr;
    Node<SPACE_DIM>* p_basal_node_2_right = nullptr;

    // Find the nodes at the basal edges which intersect the plane
    c_vector<double, SPACE_DIM> centroid = this->GetCentroidOfElement(element_index);
    Node<SPACE_DIM>* current_node = p_basal_face->GetNode(0);
    unsigned num_nodes = p_basal_face->GetNumNodes();
    for (unsigned ind_node = 0; ind_node < num_nodes; ind_node++)
    {
        Node<SPACE_DIM>* next_node = p_basal_face->GetNode((ind_node + 1) % num_nodes);

        c_vector<double, SPACE_DIM> centr_to_curr = this->GetVectorFromAtoB(centroid, current_node->rGetLocation());
        c_vector<double, SPACE_DIM> centr_to_nxt = this->GetVectorFromAtoB(centroid, next_node->rGetLocation());

        // If the scalar product with the normal changes sign, we switched the side of the plane
        if (inner_prod(centr_to_curr, axisOfDivision) * inner_prod(centr_to_nxt, axisOfDivision) < 0.0)
        {
            if (p_basal_node_1_left == nullptr)
            {
                p_basal_node_1_left = current_node;
                p_basal_node_1_right = next_node;
            }
            else
            {
                p_basal_node_2_left = current_node;
                p_basal_node_2_right = next_node;
            }
        }
        current_node = next_node;
    }
    if (p_basal_node_2_left == nullptr || p_basal_node_1_left == nullptr)
    {
        // If the axis defines division plane which does not intersect any basal edge, we try again with
        // an axis, which is only given by the long axis of the basal side. If we have tried this and it
        // still does not work, we throw the exception
        if (iterative_call_counter > 1)
            EXCEPTION("Tried to divide element with an axis, where the division plane does not intersect any basal edges. Even the basal-correction failed!");

        c_vector<double, SPACE_DIM> newAxisOfDivision = this->GetBasalLongAxisOfElement(element_index);
        return DivideElementAlongGivenAxis(pElement, newAxisOfDivision, placeOriginalElementInDirection, iterative_call_counter + 1);
    }

    // Add the two new apical nodes
    bool boundary_1 = p_basal_node_1_left->IsBoundaryNode() && p_basal_node_1_right->IsBoundaryNode();
    c_vector<double, SPACE_DIM> dir_ba_1 = this->GetVectorFromAtoB(p_basal_node_1_left->rGetLocation(), p_basal_node_1_right->rGetLocation());
    c_vector<double, SPACE_DIM> centr_to_ba_1_l = this->GetVectorFromAtoB(centroid, p_basal_node_1_left->rGetLocation());

    /* If c centroid, p point 1, l vector from point 1 to point 2 and n normal, we get the intersection distance
     *  d= <(p-c),n> / <l,n>
     * the intersection is then at p+d*l
     */
    double fraction_1 = abs((inner_prod(centr_to_ba_1_l, axisOfDivision)) / (inner_prod(dir_ba_1, axisOfDivision)));
    c_vector<double, SPACE_DIM> pos_basal_new_1 = p_basal_node_1_left->rGetLocation() + fraction_1 * dir_ba_1;
    Node<SPACE_DIM>* p_basal_new_1 = new Node<SPACE_DIM>(0, boundary_1, pos_basal_new_1[0], pos_basal_new_1[1], pos_basal_new_1[2]);
    this->AddNode(p_basal_new_1);

    bool boundary_2 = p_basal_node_2_left->IsBoundaryNode() && p_basal_node_2_right->IsBoundaryNode();
    c_vector<double, SPACE_DIM> dir_ba_2 = this->GetVectorFromAtoB(p_basal_node_2_left->rGetLocation(), p_basal_node_2_right->rGetLocation());
    c_vector<double, SPACE_DIM> centr_to_ba_2_l = this->GetVectorFromAtoB(centroid, p_basal_node_2_left->rGetLocation());

    double fraction_2 = abs((inner_prod(centr_to_ba_2_l, axisOfDivision)) / (inner_prod(dir_ba_2, axisOfDivision)));
    c_vector<double, SPACE_DIM> pos_basal_new_2 = p_basal_node_2_left->rGetLocation() + fraction_2 * dir_ba_2;
    Node<SPACE_DIM>* p_basal_new_2 = new Node<SPACE_DIM>(0, boundary_2, pos_basal_new_2[0], pos_basal_new_2[1], pos_basal_new_2[2]);
    this->AddNode(p_basal_new_2);

    /*
     * Determine the apical nodes
     */
    Node<SPACE_DIM>* p_apical_node_1_left = nullptr;
    Node<SPACE_DIM>* p_apical_node_1_right = nullptr;
    Node<SPACE_DIM>* p_apical_node_2_left = nullptr;
    Node<SPACE_DIM>* p_apical_node_2_right = nullptr;

    unsigned num_faces = pElement->GetNumFaces();
    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*p_lateral_face_1{}, *p_lateral_face_2{};
    std::set<unsigned> apical_edge_1, apical_edge_2;

    // Go through the lateral faces to find the corresponding basal edges for the two faces
    // we call the basal node which is connected to the 'left' apical node 'left'.
    std::set<unsigned> basal_edge_1 = { p_basal_node_1_left->GetIndex(), p_basal_node_1_right->GetIndex() };
    std::set<unsigned> basal_edge_2 = { p_basal_node_2_left->GetIndex(), p_basal_node_2_right->GetIndex() };
    for (unsigned ind_face = 0; ind_face < num_faces; ind_face++)
    {
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = pElement->GetFace(ind_face);
        if (p_face->GetFaceType() != MonolayerVertexElementType::Lateral)
            continue;
        unsigned num_nodes_face = p_face->GetNumNodes();
        Node<SPACE_DIM>* current_node = p_face->GetNode(0);
        unsigned global_index_current = current_node->GetIndex();
        Node<SPACE_DIM>*p_temp_apical_1{}, *p_temp_apical_2{};
        bool relevant_face_1 = false;
        bool relevant_face_2 = false;
        bool ordering_basal_left_right = false;
        // Iterate through nodes and search for the faces with the two edges
        for (unsigned ind_node = 0; ind_node < num_nodes_face; ind_node++)
        {
            Node<SPACE_DIM>* next_node = p_face->GetNode((ind_node + 1) % num_nodes_face);
            unsigned global_index_next = next_node->GetIndex();
            // Save the basal nodes if we need them
            if (p_face->GetNodeType(ind_node) == MonolayerVertexElementType::Apical && p_face->GetNodeType((ind_node + 1) % num_nodes_face) == MonolayerVertexElementType::Apical)
            {
                p_temp_apical_1 = current_node;
                p_temp_apical_2 = next_node;
            }
            else
            {
                current_node = next_node;
                global_index_current = global_index_next;
                continue;
            }
        }
        for (unsigned ind_node = 0; ind_node < num_nodes_face; ind_node++)
        {
            Node<SPACE_DIM>* next_node = p_face->GetNode((ind_node + 1) % num_nodes_face);
            unsigned global_index_next = next_node->GetIndex();

            // Check if this is one of the relevant faces
            std::set<unsigned> edge = { global_index_current, global_index_next };
            if (edge == basal_edge_1)
            {
                if (p_basal_node_1_left == current_node) // is left the first in ordering?
                    ordering_basal_left_right = true;
                p_lateral_face_1 = p_face;
                relevant_face_1 = true;
            }
            else if (edge == basal_edge_2)
            {
                if (p_basal_node_2_left == current_node) // is left the first in ordering?
                    ordering_basal_left_right = true;
                p_lateral_face_2 = p_face;
                relevant_face_2 = true;
            }
            current_node = next_node;
            global_index_current = global_index_next;
        }

        if (relevant_face_1)
        {
            p_apical_node_1_left = ordering_basal_left_right ? p_temp_apical_2 : p_temp_apical_1;
            p_apical_node_1_right = ordering_basal_left_right ? p_temp_apical_1 : p_temp_apical_2;
        }
        else if (relevant_face_2)
        {
            p_apical_node_2_left = ordering_basal_left_right ? p_temp_apical_2 : p_temp_apical_1;
            p_apical_node_2_right = ordering_basal_left_right ? p_temp_apical_1 : p_temp_apical_2;
        }
    }
    assert(p_apical_node_2_left != nullptr && p_apical_node_1_left != nullptr);

    // Go through the apical nodes and check whether the edge is divided by the division plane
    // MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_apical_face = this->GetFaceOfType(element_index, MonolayerVertexElementType::Apical);

    // Find the nodes at the apical edges which intersect the plane or minimize the distance to division plane
    // EDGE 1
    c_vector<double, SPACE_DIM> centr_to_apical_1_l = this->GetVectorFromAtoB(centroid, p_apical_node_1_left->rGetLocation());
    c_vector<double, SPACE_DIM> centr_to_apical_1_r = this->GetVectorFromAtoB(centroid, p_apical_node_1_right->rGetLocation());
    c_vector<double, SPACE_DIM> pos_apical_new_1 = zero_vector<double>(3);
    // If this edge is divided...
    if (inner_prod(centr_to_apical_1_l, axisOfDivision) * inner_prod(centr_to_apical_1_r, axisOfDivision) < 0.0)
    {
        // ... we create the new basal node at the division point
        c_vector<double, SPACE_DIM> dir_ap_1 = this->GetVectorFromAtoB(p_apical_node_1_right->rGetLocation(), p_apical_node_1_left->rGetLocation());
        c_vector<double, SPACE_DIM> centr_to_ap_1_r = this->GetVectorFromAtoB(centroid, p_apical_node_1_right->rGetLocation());

        /* If c centroid, p point 1, l vector from point 1 to point 2 and n normal, we get the intersection distance
         *  d= <(p-c),n> / <l,n>
         * the intersection is then at p+d*l
         */
        double fraction_1 = abs((inner_prod(centr_to_ap_1_r, axisOfDivision)) / (inner_prod(dir_ap_1, axisOfDivision)));
        pos_apical_new_1 = p_apical_node_1_right->rGetLocation() + fraction_1 * dir_ap_1;
    }
    else
    {
        // ... otherwise we minimize the distance to the plane, taking a minimum distance of the old node
        //     of 1 % of the old length of the edge.
        bool right_is_closer = abs(inner_prod(centr_to_apical_1_l, axisOfDivision)) > abs(inner_prod(centr_to_apical_1_r, axisOfDivision));
        c_vector<double, SPACE_DIM> dir_ap_1 = this->GetVectorFromAtoB(p_apical_node_1_right->rGetLocation(), p_apical_node_1_left->rGetLocation());
        if (right_is_closer)
            pos_apical_new_1 = p_apical_node_1_right->rGetLocation() + 0.01 * dir_ap_1;
        else
            pos_apical_new_1 = p_apical_node_1_left->rGetLocation() - 0.01 * dir_ap_1;
    }
    boundary_1 = p_apical_node_1_left->IsBoundaryNode() && p_apical_node_1_right->IsBoundaryNode();
    Node<SPACE_DIM>* p_apical_new_1 = new Node<SPACE_DIM>(0, boundary_1, pos_apical_new_1[0], pos_apical_new_1[1], pos_apical_new_1[2]);
    this->AddNode(p_apical_new_1);

    // EDGE 2
    c_vector<double, SPACE_DIM> centr_to_apical_2_l = this->GetVectorFromAtoB(centroid, p_apical_node_2_left->rGetLocation());
    c_vector<double, SPACE_DIM> centr_to_apical_2_r = this->GetVectorFromAtoB(centroid, p_apical_node_2_right->rGetLocation());
    c_vector<double, SPACE_DIM> pos_apical_new_2 = zero_vector<double>(3);
    // If this edge is divided...
    if (inner_prod(centr_to_apical_2_l, axisOfDivision) * inner_prod(centr_to_apical_2_r, axisOfDivision) < 0.0)
    {
        // ... we create the new basal node at the division point
        c_vector<double, SPACE_DIM> dir_ap_2 = this->GetVectorFromAtoB(p_apical_node_2_right->rGetLocation(), p_apical_node_2_left->rGetLocation());
        c_vector<double, SPACE_DIM> centr_to_ap_2_r = this->GetVectorFromAtoB(centroid, p_apical_node_2_right->rGetLocation());

        /* If c centroid, p point 1, l vector from point 1 to point 2 and n normal, we get the intersection distance
         *  d= <(p-c),n> / <l,n>
         * the intersection is then at p+d*l
         */
        double fraction_2 = abs((inner_prod(centr_to_ap_2_r, axisOfDivision)) / (inner_prod(dir_ap_2, axisOfDivision)));
        pos_apical_new_2 = p_apical_node_2_right->rGetLocation() + fraction_2 * dir_ap_2;
    }
    else
    {
        // ... otherwise we minimize the distance to the plane, taking a minimum distance of the old node
        //     of 1 % of the old length of the edge.
        bool right_is_closer = abs(inner_prod(centr_to_apical_2_l, axisOfDivision)) > abs(inner_prod(centr_to_apical_2_r, axisOfDivision));
        c_vector<double, SPACE_DIM> dir_ap_2 = this->GetVectorFromAtoB(p_apical_node_2_right->rGetLocation(), p_apical_node_2_left->rGetLocation());
        if (right_is_closer)
            pos_apical_new_2 = p_apical_node_2_right->rGetLocation() + 0.01 * dir_ap_2;
        else
            pos_apical_new_2 = p_apical_node_2_left->rGetLocation() - 0.01 * dir_ap_2;
    }
    boundary_2 = p_apical_node_2_left->IsBoundaryNode() && p_apical_node_2_right->IsBoundaryNode();
    Node<SPACE_DIM>* p_apical_new_2 = new Node<SPACE_DIM>(0, boundary_2, pos_apical_new_2[0], pos_apical_new_2[1], pos_apical_new_2[2]);
    this->AddNode(p_apical_new_2);

    // Now divide the lateral faces, where the function also takes care of apical and basal sides in shared elements
    this->DivideLateralFaceWithNodes(p_lateral_face_1, p_apical_new_1, p_basal_new_1);
    this->DivideLateralFaceWithNodes(p_lateral_face_2, p_apical_new_2, p_basal_new_2);

    // Create the array of (ordered) node indices where to divide and use the division function
    std::array<unsigned, 4> division_node_indices;
    division_node_indices[0] = pElement->GetNodeLocalIndex(p_basal_new_1->GetIndex());
    division_node_indices[1] = pElement->GetNodeLocalIndex(p_apical_new_1->GetIndex());
    division_node_indices[2] = pElement->GetNodeLocalIndex(p_apical_new_2->GetIndex());
    division_node_indices[3] = pElement->GetNodeLocalIndex(p_basal_new_2->GetIndex());

    unsigned new_element_index = this->DivideElement(pElement, division_node_indices, placeOriginalElementInDirection);

    return new_element_index;
}

/*
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::AddDivisionNodeToNeighbourApicalBasal(MonolayerVertexElement<ELEMENT_DIM,SPACE_DIM>* pElementToDivide,
                                                                                                                                                Node<SPACE_DIM>* pNewNode, Node<SPACE_DIM>* pOldNode1, Node<SPACE_DIM>* pOldNode2)
{
        // Get elements that contain both nodes
  std::set<unsigned> elements_of_node_1 = pOldNode1->rGetContainingElementIndices();
  std::set<unsigned> elements_of_node_2 = pOldNode1->rGetContainingElementIndices();

  std::set<unsigned> elements_sharing_face;
  std::set_intersection(elements_of_node_1.begin(), elements_of_node_1.end(),
                                                elements_of_node_2.begin(), elements_of_node_2.end(),
                                                std::inserter(elements_sharing_face, elements_sharing_face.begin()));

        // Go through these elements
        for(auto it=elements_sharing_face.begin(); it!=elements_sharing_face.end(); ++it)
        {
                MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_element = this->GetElement(*it);
                if(p_temp_element==pElementToDivide) // We ignore the element, which we divide and only look at neighbours
                                continue;
                unsigned num_faces = p_temp_element->GetNumFaces();
                for(unsigned ind_face = 0; ind_faces<num_faces; ind_faces++)
                {
                        // We only add it to apical/basal faces
                        if(p_temp_element->GetFaceType(ind_face)!= MonolayerVertexElementType::Apical &&
                                 p_temp_element->GetFaceType(ind_face)!= MonolayerVertexElementType::Basal)
                        {
                                continue;
                        }
                        unsigned num_nodes = p_temp_element->GetFace(ind_face)->GetNumNodes();
                        Node<SPACE_DIM>* p_current_node = p_temp_element->GetFace(ind_face)->GetNode(0);
                        for(unsigned ind_node = 0; ind_node < num_nodes; ind_node++)
                        {
                                Node<SPAE_DIM>* p_next_node = p_temp_element->GetFace(ind_face)->GetNode((ind_node+1)%num_nodes);
                                if((p_current_node==pOldNode1 && p_next_node==pOldNode2) || (p_current_node==pOldNode2 && p_next_node==pOldNode1))
                                {
                                        // If the new node is inbetween the current and next node, we add it inbetween to the face
                                        MonolayerVertexElementType node_type = p_temp_element->GetFace(ind_face)->GetNodeType(ind_node);
                                        p_temp_element->GetFace(ind_face)->AddNode(pNewNode, ind_node, node_type);
                                        break;
                                }
                        }
                }
        }

}
*/

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MutableMonolayerVertexMesh<1, 1>::DivideLateralFaceWithNodes(MonolayerVertexElement<0, 1>* pFace, Node<1>* pNewApicalNode, Node<1>* pNewBasalNode)
{
    EXCEPTION("DivideLateralFaceWithNodes only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MutableMonolayerVertexMesh<1, 2>::DivideLateralFaceWithNodes(MonolayerVertexElement<0, 2>* pFace, Node<2>* pNewApicalNode, Node<2>* pNewBasalNode)
{
    EXCEPTION("DivideLateralFaceWithNodes only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MutableMonolayerVertexMesh<1, 3>::DivideLateralFaceWithNodes(MonolayerVertexElement<0, 3>* pFace, Node<3>* pNewApicalNode, Node<3>* pNewBasalNode)
{
    EXCEPTION("DivideLateralFaceWithNodes only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MutableMonolayerVertexMesh<2, 2>::DivideLateralFaceWithNodes(MonolayerVertexElement<1, 2>* pFace, Node<2>* pNewApicalNode, Node<2>* pNewBasalNode)
{
    EXCEPTION("DivideLateralFaceWithNodes only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MutableMonolayerVertexMesh<2, 3>::DivideLateralFaceWithNodes(MonolayerVertexElement<1, 3>* pFace, Node<3>* pNewApicalNode, Node<3>* pNewBasalNode)
{
    EXCEPTION("DivideLateralFaceWithNodes only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideLateralFaceWithNodes(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace,
                                                                                    Node<SPACE_DIM>* pNewApicalNode, Node<SPACE_DIM>* pNewBasalNode)
{
    // Get Pointers to the apical and basal nodes
    // B comes after A in internal ordering
    Node<SPACE_DIM>* pApicalNodeA = nullptr;
    Node<SPACE_DIM>* pApicalNodeB = nullptr;
    Node<SPACE_DIM>* pBasalNodeA = nullptr;
    Node<SPACE_DIM>* pBasalNodeB = nullptr;
    unsigned num_nodes_face = pFace->GetNumNodes();
    assert(num_nodes_face == 4);

    // unsigned number_found_basal_nodes_in_pFace = 0;
    for (unsigned index_node = 0; index_node < num_nodes_face; index_node++)
    {
        if (pFace->GetNodeType(index_node) == MonolayerVertexElementType::Apical && pFace->GetNodeType((index_node + 1) % num_nodes_face) == MonolayerVertexElementType::Apical)
        {
            pApicalNodeA = pFace->GetNode(index_node);
            pApicalNodeB = pFace->GetNode((index_node + 1) % num_nodes_face);
        }
        else if (pFace->GetNodeType(index_node) == MonolayerVertexElementType::Basal && pFace->GetNodeType((index_node + 1) % num_nodes_face) == MonolayerVertexElementType::Basal)
        {
            pBasalNodeA = pFace->GetNode(index_node);
            pBasalNodeB = pFace->GetNode((index_node + 1) % num_nodes_face);
        }
    }
    if (pApicalNodeA == nullptr || pApicalNodeB == nullptr || pBasalNodeA == nullptr || pBasalNodeB == nullptr)
    {
        EXCEPTION("Tried to divide lateral face without two adjacent apical/basal nodes. This is not implemented.");
    }

    // We change the nodes in lateral face to contain new and apical A, basal B
    unsigned node_basal_A_local_face_index = pFace->GetNodeLocalIndex(pBasalNodeA->GetIndex());
    unsigned node_apical_B_local_face_index = pFace->GetNodeLocalIndex(pApicalNodeB->GetIndex());
    pFace->UpdateNode(node_basal_A_local_face_index, pNewBasalNode, MonolayerVertexElementType::Basal);
    pFace->UpdateNode(node_apical_B_local_face_index, pNewApicalNode, MonolayerVertexElementType::Apical);

    // Create new lateral face with same orientation
    std::vector<Node<SPACE_DIM>*> face_nodes{ pNewApicalNode, pApicalNodeB, pBasalNodeA, pNewBasalNode };
    std::vector<MonolayerVertexElementType> face_node_types{ MonolayerVertexElementType::Apical, MonolayerVertexElementType::Apical,
                                                             MonolayerVertexElementType::Basal, MonolayerVertexElementType::Basal };
    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_new_face_lateral = new MonolayerVertexElement<2, 3>(0, MonolayerVertexElementType::Lateral, face_nodes, face_node_types);
    this->AddFace(p_new_face_lateral);

    // Determine the elements which share this to-be-divided face
    std::set<unsigned> basal_nodeA_elem_indices = pBasalNodeA->rGetContainingElementIndices();
    std::set<unsigned> basal_nodeB_elem_indices = pBasalNodeB->rGetContainingElementIndices();
    std::set<unsigned> elements_sharing_face;
    std::set_intersection(basal_nodeA_elem_indices.begin(), basal_nodeA_elem_indices.end(),
                          basal_nodeB_elem_indices.begin(), basal_nodeB_elem_indices.end(),
                          std::inserter(elements_sharing_face, elements_sharing_face.begin()));

    std::set<unsigned> apical_edge{ pApicalNodeA->GetIndex(), pApicalNodeB->GetIndex() };
    std::set<unsigned> basal_edge{ pBasalNodeA->GetIndex(), pBasalNodeB->GetIndex() };

    for (auto it = elements_sharing_face.begin(); it != elements_sharing_face.end(); ++it)
    {
        // Add new lateral face
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement = this->mElements[*it];
        unsigned local_index_old_face = pElement->GetFaceLocalIndex(pFace);
        bool orientation_old_face = pElement->FaceIsOrientatedClockwise(local_index_old_face);
        pElement->AddNode(pNewApicalNode, pElement->GetNumNodes() - 1, MonolayerVertexElementType::Apical);
        pElement->AddNode(pNewBasalNode, pElement->GetNumNodes() - 1, MonolayerVertexElementType::Basal);
        pElement->AddFace(p_new_face_lateral, MonolayerVertexElementType::Lateral, orientation_old_face);

        // Add new apical node to apical face
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_apical = this->GetFaceOfType(*it, MonolayerVertexElementType::Apical);
        unsigned num_apical_nodes = p_face_apical->GetNumNodes();
        for (unsigned ind_node = 0; ind_node < num_apical_nodes; ind_node++)
        {
            std::set<unsigned> edge{ p_face_apical->GetNode(ind_node)->GetIndex(),
                                     p_face_apical->GetNode((ind_node + 1) % num_apical_nodes)->GetIndex() };
            if (edge == apical_edge)
            {
                p_face_apical->AddNode(pNewApicalNode, ind_node, MonolayerVertexElementType::Apical);
                break;
            }
        }

        // Add new basal node to basal face
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_basal = this->GetFaceOfType(*it, MonolayerVertexElementType::Basal);
        unsigned num_basal_nodes = p_face_basal->GetNumNodes();
        for (unsigned ind_node = 0; ind_node < num_basal_nodes; ind_node++)
        {
            std::set<unsigned> edge{ p_face_basal->GetNode(ind_node)->GetIndex(),
                                     p_face_basal->GetNode((ind_node + 1) % num_basal_nodes)->GetIndex() };
            if (edge == basal_edge)
            {
                p_face_basal->AddNode(pNewBasalNode, ind_node, MonolayerVertexElementType::Basal);
                break;
            }
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElementAlongShortAxis(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                                         bool placeOriginalElementBelow)
{
    EXCEPTION("DivideElementAlongShortAxis is not implemented");
    return 0;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MutableMonolayerVertexMesh<1, 1>::DivideElement(MonolayerVertexElement<1, 1>* pElement,
                                                         std::array<unsigned, 4> nodeIndices,
                                                         bool placeOriginalElementOutside)
{
    EXCEPTION("DivideElement only in 3D.");
    return 0;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MutableMonolayerVertexMesh<1, 2>::DivideElement(MonolayerVertexElement<1, 2>* pElement,
                                                         std::array<unsigned, 4> nodeIndices,
                                                         bool placeOriginalElementOutside)
{
    EXCEPTION("DivideElement only in 3D.");
    return 0;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MutableMonolayerVertexMesh<1, 3>::DivideElement(MonolayerVertexElement<1, 3>* pElement,
                                                         std::array<unsigned, 4> nodeIndices,
                                                         bool placeOriginalElementOutside)
{
    EXCEPTION("DivideElement only in 3D.");
    return 0;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MutableMonolayerVertexMesh<2, 2>::DivideElement(MonolayerVertexElement<2, 2>* pElement,
                                                         std::array<unsigned, 4> nodeIndices,
                                                         bool placeOriginalElementOutside)
{
    EXCEPTION("DivideElement only in 3D.");
    return 0;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MutableMonolayerVertexMesh<2, 3>::DivideElement(MonolayerVertexElement<2, 3>* pElement,
                                                         std::array<unsigned, 4> nodeIndices,
                                                         bool placeOriginalElementOutside)
{
    EXCEPTION("DivideElement only in 3D.");
    return 0;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideElement(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                                                           std::array<unsigned, 4> nodeIndices,
                                                                           bool placeOriginalElementOutside)

{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE

    // First create the new lateral face.
    std::vector<Node<SPACE_DIM>*> lateral_nodes;
    std::vector<MonolayerVertexElementType> lateral_types;
    for (auto it = nodeIndices.begin(); it != nodeIndices.end(); ++it)
    {
        lateral_nodes.push_back(pElement->GetNode(*it));
        lateral_types.push_back(pElement->GetNodeType(*it));
    }
    MonolayerVertexElement<SPACE_DIM - 1, SPACE_DIM>* p_face_shared_lateral = new MonolayerVertexElement<SPACE_DIM - 1, SPACE_DIM>(0, MonolayerVertexElementType::Lateral, lateral_nodes, lateral_types);
    this->AddFace(p_face_shared_lateral);

    /* We take care of the apical side.
     *      4___3____x
     *     /    :    \
     *    /     :     \
     *   0      :      x
     *    \     :     /
     *     \1___2___x/
     * We iterate through the apical side and jump over all to-be-deleted elements (x), which are
     * added to the new apical face, starting from the first shared node we encounter (2)
     * This gives two faces with the same anti-clockwise orientation
     * We also remember the edges to identify the lateral faces we need to keep/delete
     */
    unsigned old_element_index = pElement->GetIndex();
    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_old_apical_face = this->GetFaceOfType(old_element_index, MonolayerVertexElementType::Apical);

    // Iterate through nodes, remember edges and remove nodes from new element
    unsigned num_nodes = p_old_apical_face->GetNumNodes();
    Node<SPACE_DIM>* current_node = p_old_apical_face->GetNode(0);

    /* We need to check, whether the node we start with on the apical side in our analysis is part of the
     * new or old cell based on our choice, whether the old element is to be put outside or inside the new lateral face
     */
    c_vector<double, SPACE_DIM> nrml_shared_lateral = zero_vector<double>(SPACE_DIM);
    this->CalculateUnitNormalToFaceWithArea(p_face_shared_lateral, nrml_shared_lateral);

    // If the node is the edge at which we share...
    Node<SPACE_DIM>* node_to_check_side = current_node;
    if (std::find(lateral_nodes.begin(), lateral_nodes.end(), current_node) != lateral_nodes.end())
    {
        // We take the next one
        node_to_check_side = p_old_apical_face->GetNode(1);
    }
    c_vector<double, SPACE_DIM> passive_center_shared = this->GetPassiveCenterOfFace(p_face_shared_lateral);
    c_vector<double, SPACE_DIM> passive_center_to_node = this->GetVectorFromAtoB(passive_center_shared, node_to_check_side->rGetLocation());

    bool first_element_is_in_new_element = true;
    // If the orientation to the first node does not match the choice if the old element is outside
    // we can conclude that the first node has to be new.
    if ((inner_prod(passive_center_to_node, nrml_shared_lateral) > 0.0) == placeOriginalElementOutside)
    {
        first_element_is_in_new_element = false;
    }

    std::set<std::set<unsigned> > set_edges_apical_old;
    std::set<unsigned> set_nodes_removed_old;

    std::vector<Node<SPACE_DIM>*> nodes_new_apical;
    std::vector<Node<SPACE_DIM>*> nodes_to_delete_apical;
    bool node_for_new_element = first_element_is_in_new_element;
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        Node<SPACE_DIM>* next_node = p_old_apical_face->GetNode((node_index + 1) % num_nodes);
        std::set<unsigned> edge = { current_node->GetIndex(), next_node->GetIndex() };

        // Do we hit a division node?
        if (std::find(lateral_nodes.begin(), lateral_nodes.end(), current_node) != lateral_nodes.end())
        {
            // The current node is part of both faces.
            nodes_new_apical.push_back(current_node);

            // If we come from the new element, the next edge is part of the old element
            // Which is incorrect if we started at the boundary
            if (node_for_new_element == (node_index != 0))
            {
                set_edges_apical_old.insert(edge);
            }

            // Change the status of new/old node if we did not start at the separation element
            if (node_index != 0)
            {
                node_for_new_element = !node_for_new_element;
            }
        }
        // If this node is part of the old element...
        else if (!node_for_new_element)
        {
            // remmeber edge for old
            set_edges_apical_old.insert(edge);
        }
        else // if not...
        {
            // Remove node from old and add it to new apical face
            nodes_new_apical.push_back(current_node);
            // We delete later to not change the indices.
            nodes_to_delete_apical.push_back(current_node);
            set_nodes_removed_old.insert(current_node->GetIndex());
        }
        current_node = next_node;
    }

    // Delete the nodes
    for (auto it = nodes_to_delete_apical.begin(); it != nodes_to_delete_apical.end(); ++it)
    {
        unsigned local_index = p_old_apical_face->GetNodeLocalIndex((*it)->GetIndex());
        p_old_apical_face->DeleteNode(local_index);
    }

    // Create the new apical face
    std::vector<MonolayerVertexElementType> apical_node_types(nodes_new_apical.size(), MonolayerVertexElementType::Apical);
    MonolayerVertexElement<SPACE_DIM - 1, SPACE_DIM>* p_face_apical_new = new MonolayerVertexElement<SPACE_DIM - 1, SPACE_DIM>(0, MonolayerVertexElementType::Apical, nodes_new_apical, apical_node_types);
    this->AddFace(p_face_apical_new);

    /*
     * We now look at the lateral faces to determine which faces are part
     * of the new/old cell and to save the lateral edges connecting new/old
     * apical and basal nodes. This we need, because the starting node could
     * be on the wrong side of the normal if the shared face is non-coplanar
     */

    std::set<std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> > apicobasal_edges_lateral_new;

    // Vector of new faces
    std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> faces_new_element;
    std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> faces_delete_old_element;
    std::vector<bool> orientations_new_element;

    // Iterate through the outer lateral faces, which we need to distribute between the elements
    unsigned num_faces_old = pElement->GetNumFaces();
    for (unsigned face_index = 0; face_index < num_faces_old; face_index++)
    {
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = pElement->GetFace(face_index);
        if (p_face->GetFaceType() != MonolayerVertexElementType::Lateral)
            continue;
        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> apicobasal_edge_1 = std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>(nullptr, nullptr);
        std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> apicobasal_edge_2 = std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>(nullptr, nullptr);

        unsigned num_nodes_in_face = p_face->GetNumNodes();
        // Check whether face has edge belonging to new apical side of old element
        bool face_for_new = true;
        for (unsigned index_node_in_face = 0; index_node_in_face < num_nodes_in_face; index_node_in_face++)
        {
            unsigned ind1 = p_face->GetNode(index_node_in_face)->GetIndex();
            unsigned ind2 = p_face->GetNode((index_node_in_face + 1) % num_nodes_in_face)->GetIndex();
            std::set<unsigned> edge = { ind1, ind2 };
            if (std::find(set_edges_apical_old.begin(), set_edges_apical_old.end(), edge) != set_edges_apical_old.end())
            {
                // This face is kept and nothing happens if it is in old
                face_for_new = false;
            }
            // Remember apico-basal edge if it is one
            if (p_face->GetNodeType(index_node_in_face) != p_face->GetNodeType((index_node_in_face + 1) % num_nodes_in_face))
            {
                Node<SPACE_DIM>* apical_node = (p_face->GetNodeType(index_node_in_face) == MonolayerVertexElementType::Apical) ? p_face->GetNode(index_node_in_face) : p_face->GetNode((index_node_in_face + 1) % num_nodes_in_face);
                Node<SPACE_DIM>* basal_node = (p_face->GetNodeType(index_node_in_face) == MonolayerVertexElementType::Apical) ? p_face->GetNode((index_node_in_face + 1) % num_nodes_in_face) : p_face->GetNode(index_node_in_face);
                if (apicobasal_edge_1.first == nullptr)
                    apicobasal_edge_1 = std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>(apical_node, basal_node);
                else
                    apicobasal_edge_2 = std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>(apical_node, basal_node);
            }
        }
        if (face_for_new)
        {
            // Save face as new if it is
            faces_new_element.push_back(p_face);
            orientations_new_element.push_back(pElement->FaceIsOrientatedClockwise(face_index));
            // We delete later to not change the indices.
            faces_delete_old_element.push_back(p_face);

            // Save apicobasal_edges as new, if it belongs to new face
            // n.b.: They can be shared by old and new cell!
            apicobasal_edges_lateral_new.insert(apicobasal_edge_1);
            apicobasal_edges_lateral_new.insert(apicobasal_edge_2);
        }
    }

    /*
     * Take care of the basal side
     */
    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_old_basal_face = this->GetFaceOfType(old_element_index, MonolayerVertexElementType::Basal);
    std::vector<Node<SPACE_DIM>*> nodes_new_basal;
    std::vector<Node<SPACE_DIM>*> nodes_to_delete_basal;
    current_node = p_old_basal_face->GetNode(0);

    /* We need to check, whether the node we start with on the apical side in our analysis is part of the
     * new or old cell based on our choice, whether the old element is to be put outside or inside the new lateral face
     * For this we have to check whether the corresponding edge is in a new face.
     */
    // If the node is the edge at which we share...
    node_to_check_side = current_node;
    if (std::find(lateral_nodes.begin(), lateral_nodes.end(), current_node) != lateral_nodes.end())
    {
        // We take the next one
        node_to_check_side = p_old_apical_face->GetNode(1);
    }

    // Go through all the new edges and check whether this node has a corresponding new apical node
    first_element_is_in_new_element = false;
    for (auto it = apicobasal_edges_lateral_new.begin(); it != apicobasal_edges_lateral_new.end(); ++it)
    {
        if ((*it).second == node_to_check_side)
        {
            first_element_is_in_new_element = true;
            break;
        }
    }

    node_for_new_element = first_element_is_in_new_element;
    // Same number as apical nodes needed!
    assert(num_nodes == p_old_basal_face->GetNumNodes());
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        Node<SPACE_DIM>* next_node = p_old_basal_face->GetNode((node_index + 1) % num_nodes);

        // Do we hit a division node?
        if (std::find(lateral_nodes.begin(), lateral_nodes.end(), current_node) != lateral_nodes.end())
        {
            // The current node is part of both faces.
            nodes_new_basal.push_back(current_node);

            // Change the status of new/old node if we did not start at the separation element
            if (node_index != 0)
            {
                node_for_new_element = !node_for_new_element;
            }
        }
        // If this node is part of the old element...
        else if (!node_for_new_element)
        {
            // Do nothing
        }
        else // if not...
        {
            // Remove node from old and add it to new basal face
            nodes_new_basal.push_back(current_node);
            // We delete later to not change the indices.
            nodes_to_delete_basal.push_back(current_node);
            set_nodes_removed_old.insert(current_node->GetIndex());
        }
        current_node = next_node;
    }

    // Delete the nodes
    for (auto it = nodes_to_delete_basal.begin(); it != nodes_to_delete_basal.end(); ++it)
    {
        unsigned local_index = p_old_basal_face->GetNodeLocalIndex((*it)->GetIndex());
        p_old_basal_face->DeleteNode(local_index);
    }

    // Create the new basal face
    assert(nodes_new_basal.size() == nodes_new_apical.size());
    std::vector<MonolayerVertexElementType> basal_node_types(nodes_new_basal.size(), MonolayerVertexElementType::Basal);
    MonolayerVertexElement<SPACE_DIM - 1, SPACE_DIM>* p_face_basal_new = new MonolayerVertexElement<SPACE_DIM - 1, SPACE_DIM>(0, MonolayerVertexElementType::Basal, nodes_new_basal, basal_node_types);
    this->AddFace(p_face_basal_new);

    /*
     * Create the element and delete old faces and nodes
     */

    // Add apical and basal faces to new element
    faces_new_element.push_back(p_face_apical_new);
    faces_new_element.push_back(p_face_basal_new);
    orientations_new_element.push_back(false); // always wrong orientation == false
    orientations_new_element.push_back(false);

    // Delete the faces from element where they do not belong
    for (auto it = faces_delete_old_element.begin(); it != faces_delete_old_element.end(); ++it)
    {
        unsigned local_index = pElement->GetFaceLocalIndex(*it);
        pElement->DeleteFace(local_index);
    }

    // Delete deleted nodes from old element
    std::vector<Node<SPACE_DIM>*> unique_nodes_to_delete;
    for (auto it = set_nodes_removed_old.begin(); it != set_nodes_removed_old.end(); ++it)
    {
        // We need to work with the node pointers as the indices will change
        unique_nodes_to_delete.push_back(this->GetNode(*it));
    }
    for (auto it = unique_nodes_to_delete.begin(); it != unique_nodes_to_delete.end(); ++it)
    {
        unsigned local_index = pElement->GetNodeLocalIndex((*it)->GetIndex());
        pElement->DeleteNode(local_index);
    }

    /*
     * Now we want to take care of the lateral face which is shared.
     * Calculate the normal of the new shared lateral face and check
     * whether it points towards the new or old element.
     * This orientation will be used to add the face to the two elements.
     */
    // nrml_shared_lateral = zero_vector<double>(SPACE_DIM);
    // this->CalculateUnitNormalToFaceWithArea(p_face_shared_lateral, nrml_shared_lateral);

    c_vector<double, SPACE_DIM> old_center = this->GetCentroidOfElement(pElement->GetIndex());
    // passive_center_shared = this->GetPassiveCenterOfFace(p_face_shared_lateral);
    c_vector<double, SPACE_DIM> face_to_old_center = this->GetVectorFromAtoB(passive_center_shared, old_center);

    bool points_outside_old = inner_prod(nrml_shared_lateral, face_to_old_center) < 0.0;

    pElement->AddFace(p_face_shared_lateral, MonolayerVertexElementType::Lateral, !points_outside_old);

    faces_new_element.push_back(p_face_shared_lateral);
    orientations_new_element.push_back(points_outside_old);

    // Get the index of the new element
    unsigned new_element_index;
    if (mDeletedElementIndices.empty())
    {
        new_element_index = this->mElements.size();
    }
    else
    {
        new_element_index = mDeletedElementIndices.back();
        mDeletedElementIndices.pop_back();
        delete this->mElements[new_element_index];
    }

    // Create and Add the new element to the mesh
    this->AddElement(new MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>(new_element_index, MonolayerVertexElementType::Undetermined, faces_new_element, orientations_new_element));

    // Now swap elements if necessary
    /*if(placeOriginalElementOutside != !points_outside_old)
    {
            unsigned old_element_index = pElement->GetIndex();
            MonolayerVertexElement<ELEMENT_DIM,SPACE_DIM>* p_temp = this->mElements[new_element_index];
            this->mElements[new_element_index] = this->mElements[old_element_index];
            this->mElements[old_element_index] = p_temp;
            this->mElements[new_element_index]->SetIndex(new_element_index);
            this->mElements[old_element_index]->SetIndex(old_element_index);
    }*/

    return new_element_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::DeleteElementPriorToReMesh(unsigned index)
{
    if (SPACE_DIM == 2)
    {
        // Mark any nodes that are contained only in this element as deleted
        for (unsigned i = 0; i < this->mElements[index]->GetNumNodes(); i++)
        {
            Node<SPACE_DIM>* p_node = this->mElements[index]->GetNode(i);

            if (p_node->rGetContainingElementIndices().size() == 1)
            {
                DeleteNodePriorToReMesh(p_node->GetIndex());
            }

            // Mark all the nodes contained in the removed element as boundary nodes
            p_node->SetAsBoundaryNode(true);
        }

        // Mark this element as deleted
        this->mElements[index]->MarkAsDeleted();
        mDeletedElementIndices.push_back(index);
    }
    else
    {
        EXCEPTION("DeleteElementPriorToReMesh only implemented in two dimensions.");
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::DeleteFacePriorToReMesh(unsigned index)
{
    // Mark this element as deleted
    this->mFaces[index]->MarkAsDeleted();
    mDeletedFaceIndices.push_back(index);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::DeleteNodePriorToReMesh(unsigned index)
{
    this->mNodes[index]->MarkAsDeleted();
    mDeletedNodeIndices.push_back(index);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::DivideEdge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB)
{
    EXCEPTION("DivideEdge is not implemented");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::RemoveDeletedNodesAndElements(VertexElementMap& rElementMap)
{
    // Make sure the map is big enough.  Each entry will be set in the loop below.
    rElementMap.Resize(this->GetNumAllElements());

    // Remove any elements that have been marked for deletion and store all other elements in a temporary structure
    std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> live_elements;
    for (unsigned i = 0; i < this->mElements.size(); i++)
    {
        if (this->mElements[i]->IsDeleted())
        {
            delete this->mElements[i];
            rElementMap.SetDeleted(i);
        }
        else
        {
            live_elements.push_back(this->mElements[i]);
            rElementMap.SetNewIndex(i, (unsigned)(live_elements.size() - 1));
        }
    }

    // Sanity check
    assert(mDeletedElementIndices.size() == this->mElements.size() - live_elements.size());

    // Repopulate the elements vector and reset the list of deleted element indices
    mDeletedElementIndices.clear();
    this->mElements = live_elements;

    // Finally, reset the element indices to run from zero
    for (unsigned i = 0; i < this->mElements.size(); i++)
    {
        this->mElements[i]->ResetIndex(i);
    }

    // Remove deleted nodes
    RemoveDeletedNodes();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::RemoveDeletedNodes()
{
    // Remove any nodes that have been marked for deletion and store all other nodes in a temporary structure
    std::vector<Node<SPACE_DIM>*> live_nodes;
    for (unsigned i = 0; i < this->mNodes.size(); i++)
    {
        if (this->mNodes[i]->IsDeleted())
        {
            delete this->mNodes[i];
        }
        else
        {
            live_nodes.push_back(this->mNodes[i]);
        }
    }

    // Sanity check
    assert(mDeletedNodeIndices.size() == this->mNodes.size() - live_nodes.size());

    // Repopulate the nodes vector and reset the list of deleted node indices
    this->mNodes = live_nodes;
    mDeletedNodeIndices.clear();

    // Finally, reset the node indices to run from zero
    for (unsigned i = 0; i < this->mNodes.size(); i++)
    {
        this->mNodes[i]->SetIndex(i);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::RemoveDeletedFaces()
{
    EXCEPTION("We cannot remove deleted faces in these dimensions.");
}

template <>
void MutableMonolayerVertexMesh<3, 3>::RemoveDeletedFaces()
{
    // Remove any faces that have been marked for deletion and store all other nodes in a temporary structure
    std::vector<MonolayerVertexElement<2, 3>*> live_faces;
    for (unsigned i = 0; i < this->mFaces.size(); i++)
    {
        if (this->mFaces[i]->IsDeleted())
        {
            delete this->mFaces[i];
        }
        else
        {
            live_faces.push_back(this->mFaces[i]);
        }
    }

    // Sanity check
    assert(mDeletedFaceIndices.size() == this->mFaces.size() - live_faces.size());

    // Repopulate the faces vector and reset the list of deleted node indices
    this->mFaces = live_faces;
    mDeletedFaceIndices.clear();

    // Finally, reset the face indices to run from zero
    for (unsigned i = 0; i < this->mFaces.size(); i++)
    {
        this->mFaces[i]->SetIndex(i);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::IsFaceConvexAndNonSelfIntersecting(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace)
{
    unsigned num_nodes = pFace->GetNumNodes();
    c_vector<double, SPACE_DIM> passive_center = this->GetPassiveCenterOfFace(pFace);
    c_vector<double, SPACE_DIM> this_node = pFace->GetNode(0)->rGetLocation();
    c_vector<double, SPACE_DIM> next_node = pFace->GetNode((0 + 1) % num_nodes)->rGetLocation();

    c_vector<double, SPACE_DIM> vector_pc_this = this->GetVectorFromAtoB(passive_center, this_node);
    c_vector<double, SPACE_DIM> vector_pc_next = this->GetVectorFromAtoB(passive_center, next_node);

    c_vector<double, SPACE_DIM> vector_prev_edge = this->GetVectorFromAtoB(this_node, next_node);

    c_vector<double, SPACE_DIM> cross_triangle_previous = VectorProduct(vector_pc_this, vector_pc_next);
    c_vector<double, SPACE_DIM> cross_edges_previous = zero_vector<double>(3);

    /**
     * We need to go to num_nodes+1, because
     *                                              ind ind+1
     * in step num_nodes-1 we check the edges (n-2, n-1, 0) and (n-3, n-2, n-1)
     * in step num_nodes   we check the edges (n-1, 0,   1) and (n-2, n-1, 0)
     * in step num_nodes+1 we check the edges (0  , 1,   2) and (n-1, 0,   1)
     */
    for (unsigned index = 1; index < num_nodes + 2; index++)
    {
        this_node = pFace->GetNode((index) % num_nodes)->rGetLocation();
        next_node = pFace->GetNode((index + 1) % num_nodes)->rGetLocation();

        vector_pc_this = this->GetVectorFromAtoB(passive_center, this_node);
        vector_pc_next = this->GetVectorFromAtoB(passive_center, next_node);

        c_vector<double, SPACE_DIM> vector_this_edge = this->GetVectorFromAtoB(this_node, next_node);

        c_vector<double, SPACE_DIM> cross_triangle_now = VectorProduct(vector_pc_this, vector_pc_next);
        c_vector<double, SPACE_DIM> cross_edges_now = VectorProduct(vector_prev_edge, vector_this_edge);

        // Now we check whether the triangle normals next to each other have the same sign
        double triangle_scalar = inner_prod(cross_triangle_previous, cross_triangle_now);
        if (triangle_scalar < 0)
            return false;

        double edge_scalar = inner_prod(cross_edges_previous, cross_edges_now);
        // If we have already two cross products from edges we can check the sign
        if (index > 1 && edge_scalar < 0)
            return false;

        // Now save to previous quantities
        vector_prev_edge = vector_this_edge;
        cross_triangle_previous = cross_triangle_now;
        cross_edges_previous = cross_edges_now;
    }
    // After checking all normals this is (more or less) regularly oriented
    return true;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh(VertexElementMap& rElementMap)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE
    assert(ELEMENT_DIM == SPACE_DIM); // LCOV_EXCL_LINE

    if (SPACE_DIM == 3)
    {
        // Make sure the map is big enough
        rElementMap.Resize(this->GetNumAllElements());

        /*
         * To begin the remeshing process, we do not need to call Clear() and remove all current data,
         * since cell birth, rearrangement and death result only in local remeshing of a vertex-based
         * mesh. Instead, we just remove any deleted elements and nodes.
         */
        RemoveDeletedNodesAndElements(rElementMap);
        RemoveDeletedFaces();
        ClearSwappedFaceIndices();
        ClearElementsThatUnderwentT1Transitions();

        // At first we resolve (proto)rosettes from prior steps, if necessary
        if (mProtorosetteFormationProbability > 0.0)
        {
            this->CheckForRosettes();
        }

        bool recheck_mesh = true;
        while (recheck_mesh == true)
        {
            // We check for any short edges and perform swaps if necessary and possible.
            recheck_mesh = CheckForSwapsFromShortEdges();
        }

        RemoveDeletedNodes();
        /*
         * This is handled in a separate method to allow child classes to implement additional ReMeshing functionality
         * (see #2664).
         */
        // this->CheckForRosettes();

        if ((mActiveT1SwapProbability > 0.0) || (mLengthDependentActiveT1BoltzmannParameter > 0.0))
        {
            // \TODO FIX THIS WORKAROUND!!!
            // double temp = mCellRearrangementThreshold;
            // mCellRearrangementThreshold = 0.3;
            PerformActiveT1Swaps();
            // mCellRearrangementThreshold = temp;
        }

        RemoveDeletedNodesAndElements(rElementMap);
        RemoveDeletedFaces();
        ClearSwappedFaceIndices();

        // Update the map from elements to faces
        // As we do not use it during the ReMesh, we only update at the end
        this->UpdateElementsFacesMap();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::ReMesh()
{
    VertexElementMap map(GetNumElements());
    ReMesh(map);
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool MutableMonolayerVertexMesh<1, 1>::CheckForSwapsFromShortEdges()
{
    EXCEPTION("Check for Swaps only for 3D elements.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool MutableMonolayerVertexMesh<1, 2>::CheckForSwapsFromShortEdges()
{
    EXCEPTION("Check for Swaps only for 3D elements.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool MutableMonolayerVertexMesh<1, 3>::CheckForSwapsFromShortEdges()
{
    EXCEPTION("Check for Swaps only for 3D elements.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::CheckForSwapsFromShortEdges()
{
    // Loop over lateral faces, where we allow T1 transitions from short edges
    for (typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator face_iter = this->GetFaceIteratorBegin();
         face_iter != this->GetFaceIteratorEnd();
         ++face_iter)
    {
        // Only check lateral faces
        if (face_iter->GetFaceType() != MonolayerVertexElementType::Lateral)
        {
            continue;
        }

        unsigned num_nodes = face_iter->GetNumNodes();
        assert(num_nodes > 0);

        Node<SPACE_DIM>* p_basal_node_a = nullptr;
        Node<SPACE_DIM>* p_basal_node_b = nullptr;

        // Caveat: We assume that we only have one apical and one basal edge
        // If both edges are too short, we allows swaps
        unsigned num_found_short_edges = 0;
        for (unsigned local_index = 0; local_index < num_nodes; local_index++)
        {
            Node<SPACE_DIM>* p_current_node = face_iter->GetNode(local_index);
            Node<SPACE_DIM>* p_next_node = face_iter->GetNode((local_index + 1) % num_nodes);

            MonolayerVertexElementType current_node_type = face_iter->GetNodeType(local_index);
            MonolayerVertexElementType next_node_type = face_iter->GetNodeType((local_index + 1) % num_nodes);

            // Check if neighbouring nodes are either both apical or basal and too short
            if (current_node_type == MonolayerVertexElementType::Apical && next_node_type == MonolayerVertexElementType::Apical)
            {
                double distance_apical_nodes = this->GetDistanceBetweenNodes(p_current_node->GetIndex(), p_next_node->GetIndex());
                if (distance_apical_nodes >= mCellRearrangementThreshold)
                {
                    break;
                }
                else // If distance too small, we have found a short edge
                {
                    num_found_short_edges++;
                }
            }
            else if (current_node_type == MonolayerVertexElementType::Basal && next_node_type == MonolayerVertexElementType::Basal)
            {
                double distance_basal_nodes = this->GetDistanceBetweenNodes(p_current_node->GetIndex(), p_next_node->GetIndex());
                if (distance_basal_nodes >= mCellRearrangementThreshold)
                {
                    break;
                }
                else // If distance too small, we have found a short edge
                {
                    num_found_short_edges++;
                    p_basal_node_a = p_current_node;
                    p_basal_node_b = p_next_node;
                }
            }
        }
        // If indeed both edges are too short...
        // Here we should also implement the alternative of too small area...
        if (num_found_short_edges == 2)
        {
            // ...then check if any triangular elements are shared by this face...
            std::set<unsigned> elements_of_node_a = p_basal_node_a->rGetContainingElementIndices();
            std::set<unsigned> elements_of_node_b = p_basal_node_b->rGetContainingElementIndices();

            std::set<unsigned> shared_elements;
            std::set_intersection(elements_of_node_a.begin(), elements_of_node_a.end(),
                                  elements_of_node_b.begin(), elements_of_node_b.end(),
                                  std::inserter(shared_elements, shared_elements.begin()));

            bool both_nodes_share_triangular_element = false;
            for (std::set<unsigned>::const_iterator it = shared_elements.begin();
                 it != shared_elements.end();
                 ++it)
            {
                // We check whether one of the containg elements has a basal side with <=3 nodes
                MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(*it);
                unsigned num_faces_in_element = p_element->GetNumFaces();
                for (unsigned face_index_in_element = 0; face_index_in_element < num_faces_in_element; face_index_in_element++)
                {
                    MonolayerVertexElementType face_in_element_type = p_element->GetFace(face_index_in_element)->GetFaceType();
                    if (face_in_element_type == MonolayerVertexElementType::Basal)
                    {
                        if (p_element->GetFace(face_index_in_element)->GetNumNodes() <= 3)
                        {
                            both_nodes_share_triangular_element = true;
                            break;
                        }
                    }
                }
                // ...and if none are, then perform the required type of swap and halt the search, returning true
                if (!both_nodes_share_triangular_element)
                {
                    mPassiveT1TransitionsCounter++;
                    IdentifySwapType(p_basal_node_a, p_basal_node_b, &(*face_iter));
                    // Since this face was passively moved (if possible), we need not actively move it
                    mSwappedFaceIndices.insert(face_iter->GetIndex());
                    return true;
                }
            }
        }
    }
    return false;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MutableMonolayerVertexMesh<1, 1>::PerformActiveT1Swaps()
{
    EXCEPTION("Active T1 swaps only for 3D elements.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MutableMonolayerVertexMesh<1, 2>::PerformActiveT1Swaps()
{
    EXCEPTION("Active T1 swaps only for 3D elements.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MutableMonolayerVertexMesh<1, 3>::PerformActiveT1Swaps()
{
    EXCEPTION("Active T1 swaps only for 3D elements.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformActiveT1Swaps()
{
    // Loop over lateral faces, where we allow T1 transitions with the applicable probability
    for (typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator face_iter = this->GetFaceIteratorBegin();
         face_iter != this->GetFaceIteratorEnd();
         ++face_iter)
    {
        // Only check lateral faces, which were not actively considered and are note deleted
        if (face_iter->GetFaceType() != MonolayerVertexElementType::Lateral || mSwappedFaceIndices.find(face_iter->GetIndex()) != mSwappedFaceIndices.end() || face_iter->IsDeleted())
        {
            continue;
        }

        unsigned num_nodes = face_iter->GetNumNodes();
        assert(num_nodes > 0);

        Node<SPACE_DIM>* p_basal_node_a = nullptr;
        Node<SPACE_DIM>* p_basal_node_b = nullptr;

        double distance_apical_nodes = DOUBLE_UNSET;
        double distance_basal_nodes = DOUBLE_UNSET;

        // Caveat: We assume that we only have one apical and one basal edge
        for (unsigned local_index = 0; local_index < num_nodes; local_index++)
        {
            Node<SPACE_DIM>* p_current_node = face_iter->GetNode(local_index);
            Node<SPACE_DIM>* p_next_node = face_iter->GetNode((local_index + 1) % num_nodes);

            MonolayerVertexElementType current_node_type = face_iter->GetNodeType(local_index);
            MonolayerVertexElementType next_node_type = face_iter->GetNodeType((local_index + 1) % num_nodes);

            // Get length between apical nodes
            if (current_node_type == MonolayerVertexElementType::Apical && next_node_type == MonolayerVertexElementType::Apical)
            {
                distance_apical_nodes = this->GetDistanceBetweenNodes(p_current_node->GetIndex(), p_next_node->GetIndex());
            }

            // Check if neighbouring nodes are both basal and name them
            if (current_node_type == MonolayerVertexElementType::Basal && next_node_type == MonolayerVertexElementType::Basal)
            {
                p_basal_node_a = p_current_node;
                p_basal_node_b = p_next_node;
                // Get length between basal nodes
                distance_basal_nodes = this->GetDistanceBetweenNodes(p_current_node->GetIndex(), p_next_node->GetIndex());
                break;
            }
        }
        assert(p_basal_node_b != nullptr && p_basal_node_a != nullptr);

        // If one of the nodes is marked for deletion we should not perform any T1 swaps
        if (p_basal_node_b->IsDeleted() || p_basal_node_a->IsDeleted())
        {
            continue;
        }

        // We only allow active T1s for elements without protorosettes.
        // Determine the element indices
        std::set<unsigned> elements_of_node_a = p_basal_node_a->rGetContainingElementIndices();
        std::set<unsigned> elements_of_node_b = p_basal_node_b->rGetContainingElementIndices();

        std::set<unsigned> all_indices, temp_union_set;
        std::set_union(elements_of_node_a.begin(), elements_of_node_a.end(),
                       elements_of_node_b.begin(), elements_of_node_b.end(),
                       std::inserter(temp_union_set, temp_union_set.begin()));
        all_indices.swap(temp_union_set); // temp_set will be deleted, all_indices now contains all
                                          // the indices of elements
        // that touch the potentially swapping nodes

        // bool element_already_has_protorosette = false;
        //  Check all elements which contain the face
        /*
                          for (std::set<unsigned>::const_iterator it = all_indices.begin();
             it != all_indices.end();
             ++it)
        {
                                  // We check whether one of the containg elements has a basal side with <=3 nodes
                                  MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(*it);
                                  if(DoesElementHaveProtorosette(p_element))
                                  {
                                          //element_already_has_protorosette = true;
                                          break;
                                  }
                          }

                          // If one element has a protorosette, skip this face
                          if(element_already_has_protorosette)
                          {
                                  continue;
                          }
                          */
        // Check if any triangular elements are shared by this face...
        std::set<unsigned> shared_elements;
        std::set_intersection(elements_of_node_a.begin(), elements_of_node_a.end(),
                              elements_of_node_b.begin(), elements_of_node_b.end(),
                              std::inserter(shared_elements, shared_elements.begin()));

        bool both_nodes_share_triangular_element = false;
        for (std::set<unsigned>::const_iterator it = shared_elements.begin();
             it != shared_elements.end();
             ++it)
        {
            // We check whether one of the containg elements has a basal side with <=3 nodes
            MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(*it);
            unsigned num_faces_in_element = p_element->GetNumFaces();
            for (unsigned face_index_in_element = 0; face_index_in_element < num_faces_in_element; face_index_in_element++)
            {
                MonolayerVertexElementType face_in_element_type = p_element->GetFace(face_index_in_element)->GetFaceType();
                if (face_in_element_type == MonolayerVertexElementType::Basal)
                {
                    if (p_element->GetFace(face_index_in_element)->GetNumNodes() <= 3)
                    {
                        both_nodes_share_triangular_element = true;
                        break;
                    }
                }
            }
        }
        // ...and if none are, then perform the required type of swap with the appropriate probability

        // calculate the probability of performing the active T1 swap. Return mActiveT1SwapProbability if the probability should not be length dependent
        double active_t1_swap_probability = CalculateActiveT1SwapProbability(distance_apical_nodes, distance_basal_nodes);

        if ((!both_nodes_share_triangular_element) && (active_t1_swap_probability > RandomNumberGenerator::Instance()->ranf()))
        {
            mActiveT1TransitionsCounter++;
            IdentifySwapType(p_basal_node_a, p_basal_node_b, &(*face_iter), true);
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::CheckForT2Swaps(VertexElementMap& rElementMap)
{
    EXCEPTION("CheckForT2Swaps is not implemented");
    return false;
}

template <>
void MutableMonolayerVertexMesh<1, 1>::IdentifySwapType(Node<1>* pNodeA,
                                                        Node<1>* pNodeB,
                                                        MonolayerVertexElement<0, 1>* pFace,
                                                        bool performActiveT1Swaps)
{
    EXCEPTION("IdentifySwapType is not implemented in these dimensions");
}

template <>
void MutableMonolayerVertexMesh<1, 2>::IdentifySwapType(Node<2>* pNodeA,
                                                        Node<2>* pNodeB,
                                                        MonolayerVertexElement<0, 2>* pFace,
                                                        bool performActiveT1Swaps)
{
    EXCEPTION("IdentifySwapType is not implemented in these dimensions");
}

template <>
void MutableMonolayerVertexMesh<1, 3>::IdentifySwapType(Node<3>* pNodeA,
                                                        Node<3>* pNodeB,
                                                        MonolayerVertexElement<0, 3>* pFace,
                                                        bool performActiveT1Swaps)
{
    EXCEPTION("IdentifySwapType is not implemented in these dimensions");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::IdentifySwapType(Node<SPACE_DIM>* pNodeA,
                                                                          Node<SPACE_DIM>* pNodeB,
                                                                          MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace,
                                                                          bool activeT1Swap)
{
    // Find the sets of elements containing nodes A and B
    std::set<unsigned> nodeA_elem_indices = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> nodeB_elem_indices = pNodeB->rGetContainingElementIndices();

    // Form the set union
    std::set<unsigned> all_indices, temp_union_set;
    std::set_union(nodeA_elem_indices.begin(), nodeA_elem_indices.end(),
                   nodeB_elem_indices.begin(), nodeB_elem_indices.end(),
                   std::inserter(temp_union_set, temp_union_set.begin()));
    all_indices.swap(temp_union_set); // temp_set will be deleted, all_indices now contains all the indices of elements
                                      // that touch the potentially swapping nodes

    if ((nodeA_elem_indices.size() > 3) || (nodeB_elem_indices.size() > 3))
    {
        /*
         * Looks like
         *
         *  \
         *   \ A   B
         * ---o---o---
         *   /
         *  /
         *
         */

        /*
         * This case is handled in a separate method to allow child classes to implement different
         * functionality for high-order-junction remodelling events (see #2664).
         * For now not implemented in the finite thickness case
         */
        if (activeT1Swap)
            return;
        // EXCEPTION("rosette formation");
        // std::cout << "Rosette formation!" << std::flush;
        // PerformNodeMerge(pNodeA, pNodeB, pFace);
        // RemoveDeletedNodes();
        // RemoveDeletedFaces();
        this->HandleHighOrderJunctions(pNodeA, pNodeB, pFace);
    }
    else // each node is contained in at most three elements
    {
        switch (all_indices.size())
        {
            case 1:
            {
                /*
                 * Each node is contained in a single element, so the nodes must lie on the boundary
                 * of the mesh, as shown below. In this case, we merge the nodes and tidy up node
                 * indices through calls to PerformNodeMerge() and RemoveDeletedNodes().
                 *
                 *    A   B
                 * ---o---o---
                 */
                if (activeT1Swap)
                    return;
                assert(pNodeA->IsBoundaryNode());
                assert(pNodeB->IsBoundaryNode());
                assert(pFace->IsBoundaryFace());

                PerformNodeMerge(pNodeA, pNodeB, pFace);
                RemoveDeletedNodes();
                // RemoveDeletedFaces();
                break;
            }
            case 2:
            {
                if (nodeA_elem_indices.size() == 2 && nodeB_elem_indices.size() == 2)
                {
                    if (pNodeA->IsBoundaryNode() && pNodeB->IsBoundaryNode())
                    {
                        /*
                         * The node configuration is as shown below, with voids on either side. In this case
                         * we perform a T1 swap, which separates the elements.
                         *
                         *   \   /
                         *    \ / Node A
                         * (1) |   (2)      (element number in brackets)
                         *    / \ Node B
                         *   /   \
                         */
                        PerformT1Swap(pNodeA, pNodeB, pFace);
                    }
                    else if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
                    {
                        /*
                         * The node configuration is as shown below, with a void on one side. We should not
                         * be able to reach this case at present, since we allow only for three-way junctions
                         * or boundaries, so we throw an exception.
                         *
                         *   \   /
                         *    \ / Node A
                         * (1) |   (2)      (element number in brackets)
                         *     x Node B
                         *     |
                         */
                        EXCEPTION("There is a non-boundary node contained only in two elements; something has gone wrong.");
                    }
                    else
                    {
                        /*
                         * Each node is contained in two elements, so the nodes lie on an internal edge, as shown below.
                         * We should not be able to reach this case at present, since we allow only for three-way junctions
                         * or boundaries, so we throw an exception.
                         *
                         *    A   B
                         * ---o---o---
                         */
                        EXCEPTION("There are non-boundary nodes contained only in two elements; something has gone wrong.");
                    }
                } // from [if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)]
                else
                {
                    /*
                     * The node configuration looks like that shown below. In this case, we merge the nodes
                     * and tidy up node indices through calls to PerformNodeMerge() and  RemoveDeletedNodes().
                     *
                     * Outside
                     *         /
                     *   --o--o (2)
                     *     (1) \
                     *
                     * ///\todo this should be a T1 swap (see #1263 and #2401)
                     * Referring to the todo: this should probably stay a node-merge. If this is a T1 swap then
                     * the single boundary node will travel from element 1 to element 2, but still remain a single node.
                     * I.e. we would not reduce the total number of nodes in this situation.
                     */
                    if (activeT1Swap)
                        return;
                    PerformNodeMerge(pNodeA, pNodeB, pFace);
                    RemoveDeletedNodes();
                    // RemoveDeletedFaces();
                }
                break;
            }
            case 3:
            {
                if (nodeA_elem_indices.size() == 1 || nodeB_elem_indices.size() == 1)
                {
                    /*
                     * One node is contained in one element and the other node is contained in three elements.
                     * We should not be able to reach this case at present, since we allow each boundary node
                     * to be contained in at most two elements, so we throw an exception.
                     *
                     *    A   B
                     *
                     *  empty   /
                     *         / (3)
                     * ---o---o-----   (element number in brackets)
                     *  (1)    \ (2)
                     *          \
                     */
                    assert(pNodeA->IsBoundaryNode());
                    assert(pNodeB->IsBoundaryNode());

                    EXCEPTION("There is a boundary node contained in three elements something has gone wrong.");
                }
                else if (nodeA_elem_indices.size() == 2 && nodeB_elem_indices.size() == 2)
                {
                    // The short edge must be at the boundary. We need to check whether this edge is
                    // adjacent to a triangular void before we swap. If it is a triangular void, we perform a T2-type swap.
                    // If not, then we perform a normal T1 swap. I.e. in detail we need to check whether the
                    // element in nodeA_elem_indices which is not in nodeB_elem_indices contains a shared node
                    // with the element in nodeB_elem_indices which is not in nodeA_elem_indices.

                    std::set<unsigned> element_A_not_B, temp_set;
                    std::set_difference(all_indices.begin(), all_indices.end(), nodeB_elem_indices.begin(),
                                        nodeB_elem_indices.end(), std::inserter(temp_set, temp_set.begin()));
                    element_A_not_B.swap(temp_set);

                    // There must be only one such element
                    assert(element_A_not_B.size() == 1);

                    std::set<unsigned> element_B_not_A;
                    std::set_difference(all_indices.begin(), all_indices.end(), nodeA_elem_indices.begin(),
                                        nodeA_elem_indices.end(), std::inserter(temp_set, temp_set.begin()));
                    element_B_not_A.swap(temp_set);

                    // There must be only one such element
                    assert(element_B_not_A.size() == 1);

                    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_A_not_B = this->mElements[*element_A_not_B.begin()];
                    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_B_not_A = this->mElements[*element_B_not_A.begin()];

                    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_basal_A_not_B = this->GetFaceOfType(*element_A_not_B.begin(),
                                                                                                              MonolayerVertexElementType::Basal);
                    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_basal_B_not_A = this->GetFaceOfType(*element_B_not_A.begin(),
                                                                                                              MonolayerVertexElementType::Basal);

                    unsigned local_index_1 = p_basal_A_not_B->GetNodeLocalIndex(pNodeA->GetIndex());
                    unsigned next_node_1 = p_basal_A_not_B->GetNodeGlobalIndex((local_index_1 + 1) % (p_basal_A_not_B->GetNumNodes()));
                    unsigned previous_node_1 = p_basal_A_not_B->GetNodeGlobalIndex(
                        (local_index_1 + p_basal_A_not_B->GetNumNodes() - 1) % (p_basal_A_not_B->GetNumNodes()));
                    unsigned local_index_2 = p_basal_B_not_A->GetNodeLocalIndex(pNodeB->GetIndex());
                    unsigned next_node_2 = p_basal_B_not_A->GetNodeGlobalIndex(
                        (local_index_2 + 1) % (p_basal_B_not_A->GetNumNodes()));
                    unsigned previous_node_2 = p_basal_B_not_A->GetNodeGlobalIndex(
                        (local_index_2 + p_basal_B_not_A->GetNumNodes() - 1) % (p_basal_B_not_A->GetNumNodes()));

                    if (next_node_1 == previous_node_2 || next_node_2 == previous_node_1)
                    {
                        /*
                         * The node configuration looks like that shown below, and both nodes must be on the boundary.
                         * In this case we remove the void through a call to PerformVoidRemoval().
                         *
                         *    A  C  B                A      B
                         *      /\                 \        /
                         *     /v \                 \  (1) /
                         * (3)o----o (1)  or     (2) o----o (3)    (element number in brackets, v is a void)
                         *   /  (2) \                 \v /
                         *  /        \                 \/
                         *                             C
                         */
                        assert(pNodeA->IsBoundaryNode());
                        assert(pNodeB->IsBoundaryNode());
                        assert(false);

                        // Get the third basal node in the triangular void

                        unsigned nodeC_index;
                        if (next_node_1 == previous_node_2 && next_node_2 != previous_node_1)
                        {
                            nodeC_index = next_node_1;
                        }
                        else if (next_node_2 == previous_node_1 && next_node_1 != previous_node_2)
                        {
                            nodeC_index = next_node_2;
                        }
                        else
                        {
                            assert(next_node_1 == previous_node_2 && next_node_2 == previous_node_1);
                            /**
                             * Here, the triangular element would be along the short edge. Since we
                             * are already checking in CheckForSwapsFromShortEdges() whether the element
                             * is triangular, this exception is redundant for simulations. We leave it in for
                             * clarity.
                             * ///\todo: consider removing the checking for this exception (see #2401)
                             */
                            EXCEPTION("Triangular element next to triangular void, not implemented yet.");
                        }

                        if (p_element_A_not_B->GetNumNodes() == 3u || p_element_B_not_A->GetNumNodes() == 3u)
                        {
                            /**
                             * If this is true then one of the elements adjacent to the triangular void
                             * is triangular. This element will then not share the short edge that is considered
                             * for a swap. Nevertheless, it would loose an edge during the swap. We are currently
                             * not able to deal with this situation.
                             * Related to #2533 and #2401.
                             */
                            EXCEPTION("Triangular element next to triangular void, not implemented yet.");
                        }
                        Node<SPACE_DIM>* pNodeC = this->mNodes[nodeC_index];

                        // Get the two faces that contain node C
                        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_AC = nullptr;
                        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_BC = nullptr;
                        // AC
                        unsigned num_faces_A_not_B = p_element_A_not_B->GetNumFaces();
                        for (unsigned index_face = 0; index_face < num_faces_A_not_B; index_face++)
                        {
                            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_temp_face = p_element_A_not_B->GetFace(index_face);
                            unsigned local_index_A = p_temp_face->GetNodeLocalIndex(pNodeA->GetIndex());
                            unsigned local_index_C = p_temp_face->GetNodeLocalIndex(pNodeC->GetIndex());
                            if (local_index_A != UINT_MAX && local_index_C != UINT_MAX && p_temp_face->GetFaceType() == MonolayerVertexElementType::Lateral)
                            {
                                p_face_AC = p_temp_face;
                                break;
                            }
                        }
                        // BC
                        unsigned num_faces_B_not_A = p_element_B_not_A->GetNumFaces();
                        for (unsigned index_face = 0; index_face < num_faces_B_not_A; index_face++)
                        {
                            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_temp_face = p_element_B_not_A->GetFace(index_face);
                            unsigned local_index_B = p_temp_face->GetNodeLocalIndex(pNodeB->GetIndex());
                            unsigned local_index_C = p_temp_face->GetNodeLocalIndex(pNodeC->GetIndex());
                            if (local_index_B != UINT_MAX && local_index_C != UINT_MAX && p_temp_face->GetFaceType() == MonolayerVertexElementType::Lateral)
                            {
                                p_face_BC = p_temp_face;
                                break;
                            }
                        }
                        if (p_face_AC == nullptr || p_face_BC == nullptr)
                        {
                            EXCEPTION("The nodes that belong to the void do not belong to any faces");
                            // This should not happen
                        }
                        if (activeT1Swap)
                            return;
                        PerformVoidRemoval(pNodeA, pNodeB, pNodeC, pFace, p_face_BC, p_face_AC);
                    }
                    else
                    {
                        /*
                         * The node configuration looks like that below, and both nodes must lie on the boundary.
                         * In this case we perform a T1 swap.
                         *
                         *     A  B                  A  B
                         *   \ empty/              \      /
                         *    \    /                \(1) /
                         * (3) o--o (1)  or      (2) o--o (3)    (element number in brackets)
                         *    / (2)\                /    \
                         *   /      \              /empty \
                         */
                        assert(pNodeA->IsBoundaryNode());
                        assert(pNodeB->IsBoundaryNode());

                        PerformT1Swap(pNodeA, pNodeB, pFace);
                    }
                } // from else if (nodeA_elem_indices.size()==2 && nodeB_elem_indices.size()==2)
                else
                {
                    // In this case, one node must be contained in two elements and the other in three elements.
                    assert((nodeA_elem_indices.size() == 2 && nodeB_elem_indices.size() == 3)
                           || (nodeA_elem_indices.size() == 3 && nodeB_elem_indices.size() == 2));

                    // They can't both be boundary nodes
                    assert(!(pNodeA->IsBoundaryNode() && pNodeB->IsBoundaryNode()));

                    if (pNodeA->IsBoundaryNode() || pNodeB->IsBoundaryNode())
                    {
                        /*
                         * The node configuration looks like that shown below. We perform a T1 swap in this case.
                         *
                         *     A  B                      A  B
                         *   \      /                  \      /
                         *    \ (1)/                    \(1) /
                         * (3) o--o (empty)  or  (empty) o--o (3)    (element number in brackets)
                         *    / (2)\                    /(2) \
                         *   /      \                  /      \
                         */
                        PerformT1Swap(pNodeA, pNodeB, pFace);
                    }
                    else
                    {
                        /*
                         * The node configuration looks like that shown below. We should not be able to reach this case
                         * at present, since we allow only for three-way junctions or boundaries, so we throw an exception.
                         *
                         *     A  B             A  B
                         *   \                       /
                         *    \  (1)           (1)  /
                         * (3) o--o---   or  ---o--o (3)    (element number in brackets)
                         *    /  (2)           (2)  \
                         *   /                       \
                         */
                        EXCEPTION("There are non-boundary nodes contained only in two elements; something has gone wrong.");
                    }
                }
                break;
            }
            case 4:
            {
                /*
                 * The node configuration looks like that shown below. We perform a T1 swap in this case.
                 *
                 *   \(1)/
                 *    \ / Node A
                 * (2) |   (4)      (element number in brackets)
                 *    / \ Node B
                 *   /(3)\
                 */

                /*
                 * This case is handled in a separate method to allow child classes to implement different
                 * functionality for junction remodelling events (see #2664).
                 */
                if (mProtorosetteFormationProbability > RandomNumberGenerator::Instance()->ranf())
                {
                    // First find neighbouring elements of to-be-deleted face
                    std::set<unsigned> intersection_indices, temp_inter_set;
                    std::set_intersection(nodeA_elem_indices.begin(), nodeA_elem_indices.end(),
                                          nodeB_elem_indices.begin(), nodeB_elem_indices.end(),
                                          std::inserter(temp_inter_set, temp_inter_set.begin()));
                    intersection_indices.swap(temp_inter_set); // temp_set will be deleted,

                    // There must be two such elements
                    assert(intersection_indices.size() == 2);

                    std::set<unsigned>::const_iterator temp_iterator = intersection_indices.begin();
                    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_intersection_A = this->mElements[*temp_iterator];
                    temp_iterator++;
                    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_intersection_B = this->mElements[*temp_iterator];

                    // Determine pointers to apical nodes because we need to save these
                    Node<SPACE_DIM>* pApicalNodeA = nullptr;
                    Node<SPACE_DIM>* pApicalNodeB = nullptr;
                    for (unsigned index_node = 0; index_node < pFace->GetNumNodes(); index_node++)
                    {
                        unsigned index_next_node = (index_node + 1) % pFace->GetNumNodes();
                        if (pFace->GetNodeType(index_node) == MonolayerVertexElementType::Apical
                            && pFace->GetNodeType(index_next_node) == MonolayerVertexElementType::Apical)
                        {
                            pApicalNodeA = pFace->GetNode(index_node);
                            pApicalNodeB = pFace->GetNode(index_next_node);
                            break;
                        }
                    }
                    assert(pApicalNodeA != nullptr && pApicalNodeB != nullptr);
                    bool check_intersections = this->mCheckForInternalIntersections && activeT1Swap;

                    bool did_node_merge = this->PerformNodeMerge(pNodeA, pNodeB, pFace, check_intersections);
                    if (!did_node_merge)
                        return;

                    Node<SPACE_DIM>* pProtorosetteNodeBasal = pNodeA->IsDeleted() ? pNodeB : pNodeA;
                    Node<SPACE_DIM>* pProtorosetteNodeApical = pApicalNodeA->IsDeleted() ? pApicalNodeB : pApicalNodeA;
                    assert((!pProtorosetteNodeBasal->IsDeleted()) && (!pProtorosetteNodeApical->IsDeleted()));
                    // Save the protorosette information
                    AddToMapOfProtorosettes(pProtorosetteNodeBasal, pProtorosetteNodeApical, p_element_intersection_A, p_element_intersection_B);
                    // Add the elements that contain the protorosette to the correspinding set
                    for (std::set<unsigned>::const_iterator it = all_indices.begin();
                         it != all_indices.end();
                         ++it)
                    {
                        AddToSetOfElementsWithProtorosette(this->mElements[*it]);
                        mElementsThatUnderwentT1Transitions.push_back(this->mElements[*it]);
                    }

                    if (activeT1Swap)
                        return; // We should not remove faces while in the loop for activeT1Swaps
                    this->RemoveDeletedNodes();
                    // this->RemoveDeletedFaces();
                }
                else
                {
                    this->PerformT1Swap(pNodeA, pNodeB, pFace);
                }
                break;
            }
            default:
                // This can't happen
                NEVER_REACHED;
        }
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool MutableMonolayerVertexMesh<1, 1>::PerformNodeMerge(Node<1>* pBasalNodeA, Node<1>* pBasalNodeB, MonolayerVertexElement<0, 1>* pFace, bool checkIntersection)
{
    EXCEPTION("Node merge only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool MutableMonolayerVertexMesh<1, 2>::PerformNodeMerge(Node<2>* pBasalNodeA, Node<2>* pBasalNodeB, MonolayerVertexElement<0, 2>* pFace, bool checkIntersection)
{
    EXCEPTION("Node merge only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
bool MutableMonolayerVertexMesh<1, 3>::PerformNodeMerge(Node<3>* pBasalNodeA, Node<3>* pBasalNodeB, MonolayerVertexElement<0, 3>* pFace, bool checkIntersection)
{
    EXCEPTION("Node merge only in 3D.");
    return false;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformNodeMerge(Node<SPACE_DIM>* pBasalNodeA,
                                                                          Node<SPACE_DIM>* pBasalNodeB,
                                                                          MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace,
                                                                          bool checkIntersection)
{
    // Get Pointers to the apical nodes
    Node<SPACE_DIM>* pApicalNodeA = nullptr;
    Node<SPACE_DIM>* pApicalNodeB = nullptr;
    unsigned num_nodes_face = pFace->GetNumNodes();
    unsigned number_found_basal_nodes_in_pFace = 0;
    for (unsigned index_node = 0; index_node < num_nodes_face; index_node++)
    {
        if (pFace->GetNodeType(index_node) == MonolayerVertexElementType::Apical)
        {
            if (pApicalNodeA == nullptr)
            {
                pApicalNodeA = pFace->GetNode(index_node);
            }
            else
            {
                pApicalNodeB = pFace->GetNode(index_node);
            }
        }
        else
        {
            // Check if basal node in pFace
            if (pBasalNodeA == pFace->GetNode(index_node))
            {
                number_found_basal_nodes_in_pFace++;
            }
            else if (pBasalNodeB == pFace->GetNode(index_node))
            {
                number_found_basal_nodes_in_pFace++;
            }
        }
    }
    if (pApicalNodeA == nullptr || pApicalNodeB == nullptr)
    {
        EXCEPTION("Tried to merge nodes in a face without apical nodes. This is not implemented.");
    }
    // Check if basal nodes are part of pFace
    if (number_found_basal_nodes_in_pFace != 2)
    {
        EXCEPTION("The nodes given to merge do not belong to pFace");
    }

    // Find the sets of elements containing each of the nodes, sorted by index
    std::set<unsigned> basal_nodeA_elem_indices = pBasalNodeA->rGetContainingElementIndices();
    std::set<unsigned> basal_nodeB_elem_indices = pBasalNodeB->rGetContainingElementIndices();
    std::set<unsigned> apical_nodeA_elem_indices = pApicalNodeA->rGetContainingElementIndices();
    std::set<unsigned> apical_nodeB_elem_indices = pApicalNodeB->rGetContainingElementIndices();

    // Save the old node locations in case we need to reset them.
    c_vector<double, SPACE_DIM> vec_old_basal_A = pBasalNodeA->rGetLocation();
    c_vector<double, SPACE_DIM> vec_old_apical_A = pApicalNodeA->rGetLocation();

    // Move nodes A to the mid-point
    pBasalNodeA->rGetModifiableLocation() += 0.5 * this->GetVectorFromAtoB(pBasalNodeA->rGetLocation(), pBasalNodeB->rGetLocation());
    pApicalNodeA->rGetModifiableLocation() += 0.5 * this->GetVectorFromAtoB(pApicalNodeA->rGetLocation(), pApicalNodeB->rGetLocation());

    // We create data structures that save the changes to faces and elements for possible reversion
    // We save the pointer to the element/face, whether we deleted/replaced the apical or basal node,
    // the index at which we replaced/deleted and whether it was a deletion (an not replacement)
    std::vector<std::tuple<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*, MonolayerVertexElementType, unsigned, bool> > map_faces_deleted_or_replaced_B;
    std::vector<std::tuple<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*, MonolayerVertexElementType, unsigned, bool> > map_elements_deleted_or_replaced_B;

    // Update the elements previously containing basal node B to contain basal node A
    unsigned basal_node_B_index = pBasalNodeB->GetIndex();
    unsigned basal_node_A_index = pBasalNodeA->GetIndex();
    for (std::set<unsigned>::const_iterator it = basal_nodeB_elem_indices.begin(); it != basal_nodeB_elem_indices.end(); ++it)
    {
        // Find the local index of node B in this element
        unsigned node_B_local_index = this->mElements[*it]->GetNodeLocalIndex(basal_node_B_index);
        assert(node_B_local_index < UINT_MAX); // this element contains node B

        /*
         * If this element already contains node A, then just remove node B.
         * Otherwise replace it with node A in the element and remove it from mNodes.
         */
        if (basal_nodeA_elem_indices.count(*it) != 0)
        {
            this->mElements[*it]->DeleteNode(node_B_local_index);
            // save deletion data for possible reversal
            map_elements_deleted_or_replaced_B.push_back(std::make_tuple(this->mElements[*it], MonolayerVertexElementType::Basal, node_B_local_index, true));
        }
        else
        {
            // Replace node B with node A in this element
            this->mElements[*it]->UpdateNode(node_B_local_index, pBasalNodeA, MonolayerVertexElementType::Basal);
            // save deletion data for possible reversal
            map_elements_deleted_or_replaced_B.push_back(std::make_tuple(this->mElements[*it], MonolayerVertexElementType::Basal, node_B_local_index, false));
        }
        // Now we have to do the same procedure for all faces in the element which contain node B
        unsigned num_faces = this->mElements[*it]->GetNumFaces();
        for (unsigned face_in_element_index = 0; face_in_element_index < num_faces; face_in_element_index++)
        {
            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_in_element = this->mElements[*it]->GetFace(face_in_element_index);
            unsigned node_B_local_face_index = p_face_in_element->GetNodeLocalIndex(basal_node_B_index);
            unsigned node_A_local_face_index = p_face_in_element->GetNodeLocalIndex(basal_node_A_index);

            // Only do something for faces which contain node B
            if (node_B_local_face_index == UINT_MAX)
            {
                continue;
            }
            /*
             * If this element already contains node A, then just remove node B.
             * Otherwise replace it with node A in the element and remove it from mNodes.
             */
            if (node_A_local_face_index != UINT_MAX)
            {
                p_face_in_element->DeleteNode(node_B_local_face_index);
                // save deletion data for possible reversal
                map_faces_deleted_or_replaced_B.push_back(std::make_tuple(p_face_in_element, MonolayerVertexElementType::Basal, node_B_local_face_index, true));
            }
            else
            {
                // Replace node B with node A in this face
                p_face_in_element->UpdateNode(node_B_local_face_index, pBasalNodeA, MonolayerVertexElementType::Basal);
                // save deletion data for possible reversal
                map_faces_deleted_or_replaced_B.push_back(std::make_tuple(p_face_in_element, MonolayerVertexElementType::Basal, node_B_local_face_index, false));
            }
        }
    }

    // Update the elements previously containing apical node B to contain apical node A
    unsigned apical_node_B_index = pApicalNodeB->GetIndex();
    unsigned apical_node_A_index = pApicalNodeA->GetIndex();
    for (std::set<unsigned>::const_iterator it = apical_nodeB_elem_indices.begin(); it != apical_nodeB_elem_indices.end(); ++it)
    {
        // Find the local index of node B in this element
        unsigned node_B_local_index = this->mElements[*it]->GetNodeLocalIndex(apical_node_B_index);
        assert(node_B_local_index < UINT_MAX); // this element contains node B

        /*
         * If this element already contains node A, then just remove node B.
         * Otherwise replace it with node A in the element and remove it from mNodes.
         */
        if (apical_nodeA_elem_indices.count(*it) != 0)
        {
            this->mElements[*it]->DeleteNode(node_B_local_index); // save deletion data for possible reversal
            map_elements_deleted_or_replaced_B.push_back(std::make_tuple(this->mElements[*it], MonolayerVertexElementType::Apical, node_B_local_index, true));
        }
        else
        {
            // Replace node B with node A in this element
            this->mElements[*it]->UpdateNode(node_B_local_index, pApicalNodeA, MonolayerVertexElementType::Apical);
            // save deletion data for possible reversal
            map_elements_deleted_or_replaced_B.push_back(std::make_tuple(this->mElements[*it], MonolayerVertexElementType::Apical, node_B_local_index, false));
        }
        // Now we have to do the same procedure for all faces in the element which contain node B
        unsigned num_faces = this->mElements[*it]->GetNumFaces();
        for (unsigned face_in_element_index = 0; face_in_element_index < num_faces; face_in_element_index++)
        {
            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_in_element = this->mElements[*it]->GetFace(face_in_element_index);
            unsigned node_B_local_face_index = p_face_in_element->GetNodeLocalIndex(apical_node_B_index);
            unsigned node_A_local_face_index = p_face_in_element->GetNodeLocalIndex(apical_node_A_index);

            // Only do something for faces which contain node B
            if (node_B_local_face_index == UINT_MAX)
            {
                continue;
            }
            /*
             * If this element already contains node A, then just remove node B.
             * Otherwise replace it with node A in the element and remove it from mNodes.
             */
            if (node_A_local_face_index != UINT_MAX)
            {
                p_face_in_element->DeleteNode(node_B_local_face_index);
                // save deletion data for possible reversal
                map_faces_deleted_or_replaced_B.push_back(std::make_tuple(p_face_in_element, MonolayerVertexElementType::Apical, node_B_local_face_index, true));
            }
            else
            {
                // Replace node B with node A in this face
                p_face_in_element->UpdateNode(node_B_local_face_index, pApicalNodeA, MonolayerVertexElementType::Apical);
                // save deletion data for possible reversal
                map_faces_deleted_or_replaced_B.push_back(std::make_tuple(p_face_in_element, MonolayerVertexElementType::Apical, node_B_local_face_index, false));
            }
        }
    }

    // Now that we have changed the topology, we may check for intersections, if checkIntersection is true
    if (checkIntersection)
    {
        bool found_intersection = false;
        basal_nodeA_elem_indices = pBasalNodeA->rGetContainingElementIndices(); // update the set of elements with A
        for (std::set<unsigned>::const_iterator it = basal_nodeA_elem_indices.begin(); it != basal_nodeA_elem_indices.end(); ++it)
        {
            // Find the local index of node A in this element
            unsigned node_A_local_index = this->mElements[*it]->GetNodeLocalIndex(basal_node_A_index);
            assert(node_A_local_index < UINT_MAX); // this element contains node A

            // now check orientation of apical/basal faces in the element which contain node A
            unsigned num_faces = this->mElements[*it]->GetNumFaces();
            for (unsigned face_in_element_index = 0; face_in_element_index < num_faces; face_in_element_index++)
            {
                MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_in_element = this->mElements[*it]->GetFace(face_in_element_index);
                if (p_face_in_element->GetFaceType() == MonolayerVertexElementType::Basal || p_face_in_element->GetFaceType() == MonolayerVertexElementType::Apical)
                {
                    bool orientation_okay = this->IsFaceConvexAndNonSelfIntersecting(p_face_in_element);
                    if (!orientation_okay)
                    {
                        found_intersection = true;
                        break;
                    }
                }
            }
            if (found_intersection)
                break;
        }
        // If one of the elements has a basal or apical face with intersections or is non-convec,
        // we have to revert to the initial configuration
        if (found_intersection)
        {
            // Reset the positions
            pBasalNodeA->rGetModifiableLocation() = vec_old_basal_A;
            pApicalNodeA->rGetModifiableLocation() = vec_old_apical_A;

            // Re-add the nodes B to the elements:
            // -> deleted or replaced from elements
            // We iterate backwards, as the indices changed during the deletion/replacement
            for (typename std::vector<std::tuple<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*, MonolayerVertexElementType, unsigned, bool> >::reverse_iterator it = map_elements_deleted_or_replaced_B.rbegin();
                 it != map_elements_deleted_or_replaced_B.rend();
                 ++it)
            {
                Node<SPACE_DIM>* p_node_to_add = std::get<1>(*it) == MonolayerVertexElementType::Apical ? pApicalNodeB : pBasalNodeB;
                MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = std::get<0>(*it);
                unsigned old_index = std::get<2>(*it);
                unsigned num_nodes = p_element->GetNumNodes();
                bool was_deleted = std::get<3>(*it);
                if (was_deleted)
                {
                    p_element->AddNode(p_node_to_add, (old_index - 1 + num_nodes) % num_nodes, std::get<1>(*it));
                }
                else
                {
                    p_element->UpdateNode(old_index, p_node_to_add, std::get<1>(*it));
                }
            }
            // -> deleted or replaced from faces
            for (typename std::vector<std::tuple<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*, MonolayerVertexElementType, unsigned, bool> >::reverse_iterator it = map_faces_deleted_or_replaced_B.rbegin();
                 it != map_faces_deleted_or_replaced_B.rend();
                 ++it)
            {
                Node<SPACE_DIM>* p_node_to_add = std::get<1>(*it) == MonolayerVertexElementType::Apical ? pApicalNodeB : pBasalNodeB;
                MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = std::get<0>(*it);
                unsigned old_index = std::get<2>(*it);
                unsigned num_nodes = p_face->GetNumNodes();
                bool was_deleted = std::get<3>(*it);
                if (was_deleted)
                {
                    p_face->AddNode(p_node_to_add, (old_index - 1 + num_nodes) % num_nodes, std::get<1>(*it));
                }
                else
                {
                    p_face->UpdateNode(old_index, p_node_to_add, std::get<1>(*it));
                }
            }
            // We have now reverted to the initial configuration and can return false
            return false;
        }
    }
    // If we do not have intersections/non-convexity or do not check:
    // Now we delete the face from all elements that contain it
    // We assume that we have a single-layered monolayer and thus
    // only need to check the basal network
    // Form the set interesction
    std::set<unsigned> shared_elements;
    std::set_intersection(basal_nodeA_elem_indices.begin(), basal_nodeA_elem_indices.end(),
                          basal_nodeB_elem_indices.begin(), basal_nodeB_elem_indices.end(),
                          std::inserter(shared_elements, shared_elements.begin()));
    // shared_elements now contains all the indices of elements
    // that use this the face

    for (std::set<unsigned>::const_iterator it = shared_elements.begin();
         it != shared_elements.end();
         ++it)
    {
        unsigned num_faces_element = this->mElements[*it]->GetNumFaces();
        for (unsigned local_face_index = 0; local_face_index < num_faces_element; local_face_index++)
        {
            if (this->mElements[*it]->GetFace(local_face_index) == pFace)
            {
                // Delete face from element and update the map from elements to faces
                this->mElements[*it]->DeleteFace(local_face_index);
                this->UpdateElementsFacesMapOfElement(*it);
                break; // The face should be only once in an element
            }
        }
    }

    assert(!(this->mNodes[basal_node_B_index]->IsDeleted()));
    DeleteNodePriorToReMesh(basal_node_B_index);
    // this->mNodes[basal_node_B_index]->MarkAsDeleted();
    assert(!(this->mNodes[apical_node_B_index]->IsDeleted()));
    DeleteNodePriorToReMesh(apical_node_B_index);
    // this->mNodes[apical_node_B_index]->MarkAsDeleted();
    assert(!(pFace->IsDeleted()));
    pFace->MarkAsDeleted();

    mDeletedFaceIndices.push_back(pFace->GetIndex());
    // mDeletedNodeIndices.push_back(basal_node_B_index);
    // mDeletedNodeIndices.push_back(apical_node_B_index);
    //  We have merged the nodes
    return true;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MutableMonolayerVertexMesh<1, 1>::HandleHighOrderJunctions(Node<1>* pBasalNodeA, Node<1>* pBasalNodeB, MonolayerVertexElement<0, 1>* pFace)
{
    EXCEPTION("HandleHighOrderJunction only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MutableMonolayerVertexMesh<1, 2>::HandleHighOrderJunctions(Node<2>* pBasalNodeA, Node<2>* pBasalNodeB, MonolayerVertexElement<0, 2>* pFace)
{
    EXCEPTION("HandleHighOrderJunction only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MutableMonolayerVertexMesh<1, 3>::HandleHighOrderJunctions(Node<3>* pBasalNodeA, Node<3>* pBasalNodeB, MonolayerVertexElement<0, 3>* pFace)
{
    EXCEPTION("HandleHighOrderJunction only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::HandleHighOrderJunctions(Node<SPACE_DIM>* pBasalNodeA,
                                                                                  Node<SPACE_DIM>* pBasalNodeB,
                                                                                  MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace)
{
    // Find the sets of elements containing basal nodes A and B
    std::set<unsigned> nodeA_elem_indices_basal = pBasalNodeA->rGetContainingElementIndices();
    std::set<unsigned> nodeB_elem_indices_basal = pBasalNodeB->rGetContainingElementIndices();

    // We name basal node A as the higher order node
    if (nodeB_elem_indices_basal.size() > nodeA_elem_indices_basal.size())
    {
        Node<SPACE_DIM>* temp = pBasalNodeA;
        pBasalNodeA = pBasalNodeB;
        pBasalNodeB = temp;
    }

    // Determine the apical nodes
    unsigned num_nodes_lateral = pFace->GetNumNodes(); // should be always 4!
    assert(num_nodes_lateral == 4);
    Node<SPACE_DIM>* pApicalNodeA = nullptr;
    Node<SPACE_DIM>* pApicalNodeB = nullptr;
    for (unsigned local_index = 0; local_index < num_nodes_lateral; local_index++)
    {
        Node<SPACE_DIM>* p_current_node = pFace->GetNode(local_index);
        Node<SPACE_DIM>* p_next_node = pFace->GetNode((local_index + 1) % num_nodes_lateral);

        MonolayerVertexElementType current_node_type = pFace->GetNodeType(local_index);
        MonolayerVertexElementType next_node_type = pFace->GetNodeType((local_index + 1) % num_nodes_lateral);

        if (current_node_type == MonolayerVertexElementType::Apical && next_node_type == MonolayerVertexElementType::Apical)
        {
            pApicalNodeA = p_current_node;
            pApicalNodeB = p_next_node;
            break;
        }
    }
    assert(pApicalNodeA != nullptr && pApicalNodeB != nullptr);

    // Find the sets of elements containing apical nodes A and B
    std::set<unsigned> nodeA_elem_indices_apical = pApicalNodeA->rGetContainingElementIndices();
    std::set<unsigned> nodeB_elem_indices_apical = pApicalNodeB->rGetContainingElementIndices();

    // We name apical node A as the higher order node
    if (nodeB_elem_indices_apical.size() > nodeA_elem_indices_apical.size())
    {
        Node<SPACE_DIM>* temp = pApicalNodeA;
        pApicalNodeA = pApicalNodeB;
        pApicalNodeB = temp;
    }

    /*
     * Looks like
     *
     *  \
     *   \ A   B
     * ---o---o---
     *   /
     *  /
     * However if both are protorosettes we determine randomly
     */

    // We now reset the lower order node along the connecting axis such that the distance is above threshold.
    c_vector<double, SPACE_DIM>& basal_nodeA_location = pBasalNodeA->rGetModifiableLocation();
    c_vector<double, SPACE_DIM>& basal_nodeB_location = pBasalNodeB->rGetModifiableLocation();

    c_vector<double, SPACE_DIM> basal_vector_axis = this->GetVectorFromAtoB(basal_nodeA_location, basal_nodeB_location);
    basal_vector_axis *= mCellRearrangementRatio * mCellRearrangementThreshold / norm_2(basal_vector_axis);

    basal_nodeB_location = basal_nodeA_location + basal_vector_axis;

    c_vector<double, SPACE_DIM>& apical_nodeA_location = pApicalNodeA->rGetModifiableLocation();
    c_vector<double, SPACE_DIM>& apical_nodeB_location = pApicalNodeB->rGetModifiableLocation();

    c_vector<double, SPACE_DIM> apical_vector_axis = this->GetVectorFromAtoB(apical_nodeA_location, apical_nodeB_location);
    apical_vector_axis *= mCellRearrangementRatio * mCellRearrangementThreshold / norm_2(apical_vector_axis);

    apical_nodeB_location = apical_nodeA_location + apical_vector_axis;
}

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template <>
void MutableMonolayerVertexMesh<1, 1>::PerformT1Swap(Node<1>* pBasalNodeA, Node<1>* pBasalNodeB, MonolayerVertexElement<0, 1>* pFace)
{
    EXCEPTION("T1 Swaps only implemented in three dimension, for three-dimensional monolayers.");
}
template <>
void MutableMonolayerVertexMesh<1, 2>::PerformT1Swap(Node<2>* pBasalNodeA, Node<2>* pBasalNodeB, MonolayerVertexElement<0, 2>* pFace)
{
    EXCEPTION("T1 Swaps only implemented in three dimension, for three-dimensional monolayers.");
}
template <>
void MutableMonolayerVertexMesh<1, 3>::PerformT1Swap(Node<3>* pBasalNodeA, Node<3>* pBasalNodeB, MonolayerVertexElement<0, 3>* pFace)
{
    EXCEPTION("T1 Swaps only implemented in three dimension, for three-dimensional monolayers.");
}
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformT1Swap(Node<SPACE_DIM>* pBasalNodeA,
                                                                       Node<SPACE_DIM>* pBasalNodeB,
                                                                       MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace)
{
    // We name the node in correct orientation as A->B on both apical and basal sides
    unsigned num_nodes_lateral = pFace->GetNumNodes(); // should be always 4!
    if (num_nodes_lateral != 4)
    {
        EXCEPTION("Found a lateral side which does not have 4 vertices. This should not happen.");
    }

    if ((pFace->GetNodeLocalIndex(pBasalNodeB->GetIndex()) - 1 + num_nodes_lateral) % num_nodes_lateral == pFace->GetNodeLocalIndex(pBasalNodeA->GetIndex()))
    {
        // correct ordering
    }
    else if ((pFace->GetNodeLocalIndex(pBasalNodeB->GetIndex()) + 1) % num_nodes_lateral == pFace->GetNodeLocalIndex(pBasalNodeA->GetIndex()))
    {
        Node<SPACE_DIM>* temp_node = pBasalNodeA;
        pBasalNodeA = pBasalNodeB;
        pBasalNodeB = temp_node;
    }
    else
    {
        EXCEPTION("The basal edge of a lateral face seems to have more than two nodes, has wrong nodes or is not oriented correctly");
    }

    Node<SPACE_DIM>* pApicalNodeA = nullptr;
    Node<SPACE_DIM>* pApicalNodeB = nullptr;
    // unsigned number_found_basal_nodes_in_pFace = 0;
    for (unsigned local_index = 0; local_index < num_nodes_lateral; local_index++)
    {
        Node<SPACE_DIM>* p_current_node = pFace->GetNode(local_index);
        Node<SPACE_DIM>* p_next_node = pFace->GetNode((local_index + 1) % num_nodes_lateral);

        MonolayerVertexElementType current_node_type = pFace->GetNodeType(local_index);
        MonolayerVertexElementType next_node_type = pFace->GetNodeType((local_index + 1) % num_nodes_lateral);

        if (current_node_type == MonolayerVertexElementType::Apical && next_node_type == MonolayerVertexElementType::Apical)
        {
            pApicalNodeA = p_current_node;
            pApicalNodeB = p_next_node;
            break;
        }
    }
    assert(pApicalNodeA != nullptr && pApicalNodeB != nullptr);
    // Find the sets of elements containing each of the nodes, sorted by index
    std::set<unsigned> basal_nodeA_elem_indices = pBasalNodeA->rGetContainingElementIndices();
    std::set<unsigned> basal_nodeB_elem_indices = pBasalNodeB->rGetContainingElementIndices();
    std::set<unsigned> apical_nodeA_elem_indices = pApicalNodeA->rGetContainingElementIndices();
    std::set<unsigned> apical_nodeB_elem_indices = pApicalNodeB->rGetContainingElementIndices();

    // Find the elements which contain the face
    // We assume that we have a single-layered monolayer and thus
    // only need to check the basal network
    // Form the set interesction
    std::set<unsigned> shared_elements;
    std::set_intersection(basal_nodeA_elem_indices.begin(), basal_nodeA_elem_indices.end(),
                          basal_nodeB_elem_indices.begin(), basal_nodeB_elem_indices.end(),
                          std::inserter(shared_elements, shared_elements.begin()));
    // shared_elements now contains all the indices of elements
    // that use this the face

    if (shared_elements.size() > 2)
    {
        EXCEPTION("Found more than two elements which share a basal edge. This is not possible in a monolayer.");
    }

    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_plus = nullptr;
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_minus = nullptr;
    unsigned index_element_plus = 0;
    unsigned index_element_minus = 0;
    for (std::set<unsigned>::const_iterator it = shared_elements.begin();
         it != shared_elements.end();
         ++it)
    {
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_element = this->GetElement(*it);
        unsigned local_face_index = p_temp_element->GetFaceLocalIndex(pFace);
        assert(local_face_index != UINT_MAX);

        // name plus the outward oriented face and minus the other one
        bool false_orientation = p_temp_element->FaceIsOrientatedClockwise(local_face_index);
        if (false_orientation)
        {
            p_element_minus = p_temp_element;
            index_element_minus = *it;
        }
        else
        {
            p_element_plus = p_temp_element;
            index_element_plus = *it;
        }
    }

    /*
     * We are now in the situation, where for the basal network we have:
     *     \   A   /
     *       \   /					On the basal edge we have that the vector A->B
     *         A						points to the left w.r.t. the normal of the face
     *         |						Positively oiented element is denoted as (+) and
     * (-) <---|--- (+)		negatively as (-). Cell A contains basal node A
     *         V						and B basal node B
     *         B
     *       /   \					For the apical side we have that node B is in cell
     *     /   B   \				A and node A in cell B. The  vertical arrow points
     *											in the other direction then
     *
     * We assume that the apical and basal faces are oriented with normals pointing outwards.
     */

    // We want the axis between the moved nodes to be the vector
    // from (-) to (+) projected on the perpedicular space to A->B
    // where we intersect at the middle point with the old axis

    c_vector<double, SPACE_DIM> basal_vector_AB = this->GetVectorFromAtoB(pBasalNodeA->rGetLocation(),
                                                                          pBasalNodeB->rGetLocation());

    c_vector<double, SPACE_DIM> basal_nodeA_location = pBasalNodeA->rGetLocation();
    c_vector<double, SPACE_DIM> basal_nodeB_location = pBasalNodeB->rGetLocation();

    // We want the edge to lie in the surface plane as much as possible and need to take care
    // of the case without (+) or (-) cells.
    c_vector<double, SPACE_DIM> basal_vector_mp = zero_vector<double>(SPACE_DIM);
    if (p_element_plus == nullptr)
    {
        c_vector<double, SPACE_DIM> center = basal_nodeA_location + 0.5 * basal_vector_AB;
        basal_vector_mp = this->GetVectorFromAtoB(this->GetPassiveCenterOfFaceTypeInElement(index_element_minus, MonolayerVertexElementType::Basal),
                                                  center);
        mElementsThatUnderwentT1Transitions.push_back(p_element_minus);
    }
    else if (p_element_minus == nullptr)
    {
        c_vector<double, SPACE_DIM> center = basal_nodeA_location + 0.5 * basal_vector_AB;
        basal_vector_mp = this->GetVectorFromAtoB(center,
                                                  this->GetPassiveCenterOfFaceTypeInElement(index_element_plus, MonolayerVertexElementType::Basal));
        mElementsThatUnderwentT1Transitions.push_back(p_element_plus);
    }
    else
    {
        basal_vector_mp = this->GetVectorFromAtoB(this->GetPassiveCenterOfFaceTypeInElement(index_element_minus, MonolayerVertexElementType::Basal),
                                                  this->GetPassiveCenterOfFaceTypeInElement(index_element_plus, MonolayerVertexElementType::Basal));
        mElementsThatUnderwentT1Transitions.push_back(p_element_plus);
        mElementsThatUnderwentT1Transitions.push_back(p_element_minus);
    }

    // Project for orthogonality and set length
    basal_vector_mp -= inner_prod(basal_vector_AB, basal_vector_mp) / norm_2(basal_vector_AB) * (basal_vector_AB / norm_2(basal_vector_AB));
    basal_vector_mp *= mCellRearrangementRatio * mCellRearrangementThreshold / norm_2(basal_vector_mp);

    basal_nodeA_location = basal_nodeA_location + 0.5 * basal_vector_AB - 0.5 * basal_vector_mp;
    basal_nodeB_location = basal_nodeB_location - 0.5 * basal_vector_AB + 0.5 * basal_vector_mp;

    pBasalNodeA->rGetModifiableLocation() = basal_nodeA_location;
    pBasalNodeB->rGetModifiableLocation() = basal_nodeB_location;

    /*
     * For the apical network we have:
     *    before			#					after:
     *     \   A   /	  	#   \        A     /
     *       \   /			#     \           /   	The apical network is ordered anticlockwise
     *         B    		#       \    |   /		 	seen from above. We must therefore insert
     *         ^			#  (-)    B<---A   (+)	node A after B in face A
     * (-) <---|--- (+)	    #       /  	 |   \
     *         |			#     /		 V    \
     *         A		 	#   /        B     \
     *       /   \			#
     *     /   B   \		#
     *
     */

    // We want the axis between the moved nodes to be the vector
    // from (+) to (-) projected on the perpedicular space to A->B
    // where we intersect at the middle point with the old axis

    c_vector<double, SPACE_DIM> apical_vector_AB = this->GetVectorFromAtoB(pApicalNodeA->rGetLocation(),
                                                                           pApicalNodeB->rGetLocation());

    c_vector<double, SPACE_DIM> apical_nodeA_location = pApicalNodeA->rGetLocation();
    c_vector<double, SPACE_DIM> apical_nodeB_location = pApicalNodeB->rGetLocation();

    // We want the edge to lie in the surface plane as much as possible and need to take care
    // of the case without (+) or (-) cells.
    c_vector<double, SPACE_DIM> apical_vector_mp = zero_vector<double>(SPACE_DIM);
    if (p_element_plus == nullptr)
    {
        c_vector<double, SPACE_DIM> center = apical_nodeA_location + 0.5 * apical_vector_AB;
        apical_vector_mp = this->GetVectorFromAtoB(this->GetPassiveCenterOfFaceTypeInElement(index_element_minus, MonolayerVertexElementType::Apical),
                                                   center);
    }
    else if (p_element_minus == nullptr)
    {
        c_vector<double, SPACE_DIM> center = apical_nodeA_location + 0.5 * apical_vector_AB;
        apical_vector_mp = this->GetVectorFromAtoB(center,
                                                   this->GetPassiveCenterOfFaceTypeInElement(index_element_plus, MonolayerVertexElementType::Apical));
    }
    else
    {
        apical_vector_mp = this->GetVectorFromAtoB(this->GetPassiveCenterOfFaceTypeInElement(index_element_minus, MonolayerVertexElementType::Apical),
                                                   this->GetPassiveCenterOfFaceTypeInElement(index_element_plus, MonolayerVertexElementType::Apical));
    }

    // Project for orthogonality and set length
    apical_vector_mp -= inner_prod(apical_vector_AB, apical_vector_mp) / norm_2(apical_vector_AB) * (apical_vector_AB / norm_2(apical_vector_AB));
    apical_vector_mp *= mCellRearrangementRatio * mCellRearrangementThreshold / norm_2(apical_vector_mp);

    apical_nodeA_location = apical_nodeA_location + 0.5 * apical_vector_AB + 0.5 * apical_vector_mp;
    apical_nodeB_location = apical_nodeB_location - 0.5 * apical_vector_AB - 0.5 * apical_vector_mp;

    pApicalNodeA->rGetModifiableLocation() = apical_nodeA_location;
    pApicalNodeB->rGetModifiableLocation() = apical_nodeB_location;

    // Compute the mean position of pFace around which we do the T1 swap
    c_vector<double, SPACE_DIM> center_of_t1_swap = 0.5 * apical_nodeA_location - 0.25 * apical_vector_mp
        + 0.5 * basal_nodeA_location + 0.25 * basal_vector_mp;
    mLocationsOfT1Swaps.push_back(center_of_t1_swap);

    /*
     * The basal network will look like this:
     *   \       A        /
     *     \            /   	The basal network is ordered clockwise
     *       \   |    /		 	seen from above. We must therefore insert
     *  (-)    A--->B   (+)	node B before A in face A
     *       /   |    \
     *	   /     V      \
     *   /       B        \
     */

    // Add the nodes to basal and apical faces
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_A;
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_B;
    unsigned index_element_A, index_element_B;
    unsigned nodeA_local_index, nodeB_local_index;
    bool element_AorB_found = false;

    // Find the elements which only contain basal node A.
    std::set<unsigned> elements_containing_basal_A;
    std::set_difference(basal_nodeA_elem_indices.begin(), basal_nodeA_elem_indices.end(),
                        basal_nodeB_elem_indices.begin(), basal_nodeB_elem_indices.end(),
                        std::inserter(elements_containing_basal_A, elements_containing_basal_A.begin()));
    // elements_containing_basal_A now contains all the indices
    // of elements that contain A and not B

    for (std::set<unsigned>::const_iterator it = elements_containing_basal_A.begin();
         it != elements_containing_basal_A.end();
         ++it)
    {
        // If this element also contains the apical node B, which it should ...
        if (apical_nodeB_elem_indices.find(*it) != apical_nodeB_elem_indices.end())
        {
            p_element_A = this->GetElement(*it);
            index_element_A = *it;

            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_basal_A = this->GetFaceOfType(index_element_A, MonolayerVertexElementType::Basal);
            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_apical_A = this->GetFaceOfType(index_element_A, MonolayerVertexElementType::Apical);

            // ... add node B to element A just before node A on basal side and ...
            nodeA_local_index = p_face_basal_A->GetNodeLocalIndex(pBasalNodeA->GetIndex());
            assert(nodeA_local_index < UINT_MAX);
            unsigned num_nodes_face_basal_A = p_face_basal_A->GetNumNodes();
            p_face_basal_A->AddNode(pBasalNodeB, (nodeA_local_index - 1 + num_nodes_face_basal_A) % num_nodes_face_basal_A);

            // .. add node A to element A just after node B on apical side.
            nodeB_local_index = p_face_apical_A->GetNodeLocalIndex(pApicalNodeB->GetIndex());
            assert(nodeB_local_index < UINT_MAX);
            p_face_apical_A->AddNode(pApicalNodeA, nodeB_local_index);

            // We also add them to the element itself
            unsigned max_elements = p_element_A->GetNumNodes();
            p_element_A->AddNode(pBasalNodeB, max_elements - 1, MonolayerVertexElementType::Basal);
            p_element_A->AddNode(pApicalNodeA, max_elements, MonolayerVertexElementType::Apical);

            // And we add the face to the element
            p_element_A->AddFace(pFace, MonolayerVertexElementType::Lateral, false); // false because wrong ordering == false

            element_AorB_found = true;
            mElementsThatUnderwentT1Transitions.push_back(p_element_A);
        }
    }

    // Find the elements which only contain basal node B.
    std::set<unsigned> elements_containing_basal_B;
    std::set_difference(basal_nodeB_elem_indices.begin(), basal_nodeB_elem_indices.end(),
                        basal_nodeA_elem_indices.begin(), basal_nodeA_elem_indices.end(),
                        std::inserter(elements_containing_basal_B, elements_containing_basal_B.begin()));
    // elements_containing_basal_A now contains all the indices
    // of elements that contain A and not B

    for (std::set<unsigned>::const_iterator it = elements_containing_basal_B.begin();
         it != elements_containing_basal_B.end();
         ++it)
    {
        // If this element also contains the apical node A, which it should ...
        if (apical_nodeA_elem_indices.find(*it) != apical_nodeA_elem_indices.end())
        {
            p_element_B = this->GetElement(*it);
            index_element_B = *it;

            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_basal_B = this->GetFaceOfType(index_element_B, MonolayerVertexElementType::Basal);
            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_apical_B = this->GetFaceOfType(index_element_B, MonolayerVertexElementType::Apical);

            // Add node A to element B just before node B on basal side
            nodeB_local_index = p_face_basal_B->GetNodeLocalIndex(pBasalNodeB->GetIndex());
            assert(nodeB_local_index < UINT_MAX);
            unsigned num_nodes_face_basal_B = p_face_basal_B->GetNumNodes();
            p_face_basal_B->AddNode(pBasalNodeA, (nodeB_local_index - 1 + num_nodes_face_basal_B) % num_nodes_face_basal_B);

            // Add node B to element B just after node A on apical side
            nodeA_local_index = p_face_apical_B->GetNodeLocalIndex(pApicalNodeA->GetIndex());
            assert(nodeA_local_index < UINT_MAX);
            p_face_apical_B->AddNode(pApicalNodeB, nodeA_local_index);

            // We also add them to the element itself
            unsigned max_elements = p_element_B->GetNumNodes();
            p_element_B->AddNode(pBasalNodeA, max_elements - 1, MonolayerVertexElementType::Basal);
            p_element_B->AddNode(pApicalNodeB, max_elements, MonolayerVertexElementType::Apical);

            // And we add the face to the element
            p_element_B->AddFace(pFace, MonolayerVertexElementType::Lateral, true); // true because wrong ordering != false

            element_AorB_found = true;
            mElementsThatUnderwentT1Transitions.push_back(p_element_B);
        }
    }

    if (p_element_plus != nullptr)
    {
        // Now we remove nodes from the (+) and (-) elements on basal side ...
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_basal_plus = this->GetFaceOfType(index_element_plus, MonolayerVertexElementType::Basal);

        // Remove node A from (+)
        nodeA_local_index = p_face_basal_plus->GetNodeLocalIndex(pBasalNodeA->GetIndex());
        assert(nodeA_local_index < UINT_MAX);
        p_face_basal_plus->DeleteNode(nodeA_local_index);

        // ... and on apical side
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_apical_plus = nullptr;
        p_face_apical_plus = this->GetFaceOfType(index_element_plus, MonolayerVertexElementType::Apical);

        // Remove node B from (+)
        nodeB_local_index = p_face_apical_plus->GetNodeLocalIndex(pApicalNodeB->GetIndex());
        assert(nodeB_local_index < UINT_MAX);
        p_face_apical_plus->DeleteNode(nodeB_local_index);

        // Remove the face from the (+) element
        unsigned local_face_index_in_plus = p_element_plus->GetFaceLocalIndex(pFace);
        assert(local_face_index_in_plus < UINT_MAX);
        p_element_plus->DeleteFace(local_face_index_in_plus);

        // We also remove nodes from the element itself
        unsigned local_node_basalA_index_in_plus = p_element_plus->GetNodeLocalIndex(pBasalNodeA->GetIndex());
        assert(local_node_basalA_index_in_plus < UINT_MAX);
        p_element_plus->DeleteNode(local_node_basalA_index_in_plus);

        unsigned local_node_apicalB_index_in_plus = p_element_plus->GetNodeLocalIndex(pApicalNodeB->GetIndex());
        assert(local_node_apicalB_index_in_plus < UINT_MAX);
        p_element_plus->DeleteNode(local_node_apicalB_index_in_plus);

        // We need to reconnect the lateral sides of (-) and (+)
        unsigned num_faces_plus = p_element_plus->GetNumFaces();
        for (unsigned i_face = 0; i_face < num_faces_plus; i_face++)
        {
            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_in_element = p_element_plus->GetFace(i_face);
            unsigned local_index_apical_B = p_face_in_element->GetNodeLocalIndex(pApicalNodeB->GetIndex());
            if (local_index_apical_B != UINT_MAX)
            {
                // Replace node B with node A on apical side
                p_face_in_element->UpdateNode(local_index_apical_B, pApicalNodeA, MonolayerVertexElementType::Apical);
            }

            unsigned local_index_basal_A = p_face_in_element->GetNodeLocalIndex(pBasalNodeA->GetIndex());
            if (local_index_basal_A != UINT_MAX)
            {
                // Replace node A with node B on basal side
                p_face_in_element->UpdateNode(local_index_basal_A, pBasalNodeB, MonolayerVertexElementType::Basal);
            }
        }
    }
    else // from p_element_plus != nullptr
    {
        // If we do not have the (+) we need to reconnect manually
        for (std::set<unsigned>::const_iterator it = elements_containing_basal_A.begin();
             it != elements_containing_basal_A.end();
             ++it)
        {
            MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(*it);
            unsigned num_faces = p_element->GetNumFaces();
            for (unsigned i_face = 0; i_face < num_faces; i_face++)
            {
                MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_in_element = p_element->GetFace(i_face);
                unsigned local_index_basal_A = p_face_in_element->GetNodeLocalIndex(pBasalNodeA->GetIndex());
                unsigned local_index_apical_B = p_face_in_element->GetNodeLocalIndex(pApicalNodeB->GetIndex());

                unsigned local_face_index_in_minus = p_element_minus->GetFaceLocalIndex(p_face_in_element);

                // If we have a face which contains the basal B node and is not connected to (+)
                // the basal B nodes need to be changed to basal A nodes.
                if (local_index_basal_A != UINT_MAX && local_face_index_in_minus == UINT_MAX
                    && local_index_apical_B != UINT_MAX && p_face_in_element != pFace)
                {
                    // Replace node A with node B on basal side
                    p_face_in_element->UpdateNode(local_index_basal_A, pBasalNodeB, MonolayerVertexElementType::Basal);
                    // Replace node B with node A on apical side
                    p_face_in_element->UpdateNode(local_index_apical_B, pApicalNodeA, MonolayerVertexElementType::Apical);
                }
            }
        }
    }
    if (p_element_minus != nullptr)
    {
        // Now we remove nodes from the (+) and (-) elements on basal side ...
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_basal_minus = this->GetFaceOfType(index_element_minus, MonolayerVertexElementType::Basal);

        // Remove node B from (-)
        nodeB_local_index = p_face_basal_minus->GetNodeLocalIndex(pBasalNodeB->GetIndex());
        assert(nodeB_local_index < UINT_MAX);
        p_face_basal_minus->DeleteNode(nodeB_local_index);

        // ... and on apical side
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_apical_minus = nullptr;
        p_face_apical_minus = this->GetFaceOfType(index_element_minus, MonolayerVertexElementType::Apical);

        // Remove node A from (-)
        nodeA_local_index = p_face_apical_minus->GetNodeLocalIndex(pApicalNodeA->GetIndex());
        assert(nodeA_local_index < UINT_MAX);
        p_face_apical_minus->DeleteNode(nodeA_local_index);

        // Remove the face from the (-) elements
        unsigned local_face_index_in_minus = p_element_minus->GetFaceLocalIndex(pFace);
        assert(local_face_index_in_minus < UINT_MAX);
        p_element_minus->DeleteFace(local_face_index_in_minus);

        // We also remove nodes from the element itself
        unsigned local_node_basalB_index_in_minus = p_element_minus->GetNodeLocalIndex(pBasalNodeB->GetIndex());
        assert(local_node_basalB_index_in_minus < UINT_MAX);
        p_element_minus->DeleteNode(local_node_basalB_index_in_minus);

        unsigned local_node_apicalA_index_in_minus = p_element_minus->GetNodeLocalIndex(pApicalNodeA->GetIndex());
        assert(local_node_apicalA_index_in_minus < UINT_MAX);
        p_element_minus->DeleteNode(local_node_apicalA_index_in_minus);

        // We need to reconnect the lateral sides of (-) and (+)
        unsigned num_faces_minus = p_element_minus->GetNumFaces();
        for (unsigned i_face = 0; i_face < num_faces_minus; i_face++)
        {
            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_in_element = p_element_minus->GetFace(i_face);
            unsigned local_index_apical_A = p_face_in_element->GetNodeLocalIndex(pApicalNodeA->GetIndex());
            if (local_index_apical_A != UINT_MAX)
            {
                // Replace node A with node B on apical side
                p_face_in_element->UpdateNode(local_index_apical_A, pApicalNodeB, MonolayerVertexElementType::Apical);
            }

            unsigned local_index_basal_B = p_face_in_element->GetNodeLocalIndex(pBasalNodeB->GetIndex());
            if (local_index_basal_B != UINT_MAX)
            {
                // Replace node B with node A on basal side
                p_face_in_element->UpdateNode(local_index_basal_B, pBasalNodeA, MonolayerVertexElementType::Basal);
            }
        }
    }
    else // from p_element_minus != nullptr
    {
        // If we do not have the (-) we need to reconnect manually
        for (std::set<unsigned>::const_iterator it = elements_containing_basal_B.begin();
             it != elements_containing_basal_B.end();
             ++it)
        {
            MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(*it);
            unsigned num_faces = p_element->GetNumFaces();
            for (unsigned i_face = 0; i_face < num_faces; i_face++)
            {
                MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face_in_element = p_element->GetFace(i_face);
                unsigned local_index_basal_B = p_face_in_element->GetNodeLocalIndex(pBasalNodeB->GetIndex());
                unsigned local_index_apical_A = p_face_in_element->GetNodeLocalIndex(pApicalNodeA->GetIndex());

                unsigned local_face_index_in_plus = p_element_plus->GetFaceLocalIndex(p_face_in_element);

                // If we have a face which contains the basal B node and is not connected to (+)
                // the basal B nodes need to be changed to basal A nodes.
                if (local_index_basal_B != UINT_MAX && local_face_index_in_plus == UINT_MAX
                    && local_index_apical_A != UINT_MAX && p_face_in_element != pFace)
                {
                    // Replace node B with node A on basal side
                    p_face_in_element->UpdateNode(local_index_basal_B, pBasalNodeA, MonolayerVertexElementType::Basal);
                    // Replace node A with node B on apical side
                    p_face_in_element->UpdateNode(local_index_apical_A, pApicalNodeB, MonolayerVertexElementType::Apical);
                }
            }
        }
    }

    // Sort out boundary nodes
    // If one basal is boundary node, we assume that also the apical is
    if (pBasalNodeA->IsBoundaryNode() || pBasalNodeB->IsBoundaryNode())
    {
        if (pBasalNodeA->GetNumContainingElements() == 3)
        {
            pBasalNodeA->SetAsBoundaryNode(false);
            pApicalNodeB->SetAsBoundaryNode(false);
        }
        else
        {
            pBasalNodeA->SetAsBoundaryNode(true);
            pApicalNodeB->SetAsBoundaryNode(true);
        }
        if (pBasalNodeB->GetNumContainingElements() == 3)
        {
            pBasalNodeB->SetAsBoundaryNode(false);
            pApicalNodeA->SetAsBoundaryNode(false);
        }
        else
        {
            pBasalNodeB->SetAsBoundaryNode(true);
            pApicalNodeA->SetAsBoundaryNode(true);
        }
    }

    // If face is not connected to any element (since A and B were void) delete it
    if (!element_AorB_found)
    {
        pFace->MarkAsDeleted();
        mDeletedFaceIndices.push_back(pFace->GetIndex());
    }
    else // if it exists, check if boundary face
    {
        if (pBasalNodeA->IsBoundaryNode() && pBasalNodeB->IsBoundaryNode() && pApicalNodeA->IsBoundaryNode() && pApicalNodeB->IsBoundaryNode())
        {
            // Update the sets of elements containing each of the nodes, sorted by index
            basal_nodeA_elem_indices = pBasalNodeA->rGetContainingElementIndices();
            basal_nodeB_elem_indices = pBasalNodeB->rGetContainingElementIndices();
            apical_nodeA_elem_indices = pApicalNodeA->rGetContainingElementIndices();
            apical_nodeB_elem_indices = pApicalNodeB->rGetContainingElementIndices();

            // It is a boundary face if additionally only one element contains this face
            // we need this, as single-cell wide situations can have boundary nodes without face
            std::set<unsigned> shared_elements_basal;
            std::set_intersection(basal_nodeA_elem_indices.begin(), basal_nodeA_elem_indices.end(),
                                  basal_nodeB_elem_indices.begin(), basal_nodeB_elem_indices.end(),
                                  std::inserter(shared_elements_basal, shared_elements_basal.begin()));

            std::set<unsigned> shared_elements_apical;
            std::set_intersection(apical_nodeA_elem_indices.begin(), apical_nodeA_elem_indices.end(),
                                  apical_nodeB_elem_indices.begin(), apical_nodeB_elem_indices.end(),
                                  std::inserter(shared_elements_apical, shared_elements_apical.begin()));

            if (shared_elements_basal.size() == 1 && shared_elements_apical.size() == 1)
            {
                pFace->SetAsBoundaryFace(true);
            }
            else
            {
                pFace->SetAsBoundaryFace(false);
            }
        }
        else
        {
            pFace->SetAsBoundaryFace(false);
        }
    }
}

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template <>
void MutableMonolayerVertexMesh<1, 1>::PerformVoidRemoval(Node<1>* pBasalNodeA, Node<1>* pBasalNodeB, Node<1>* pBasalNodeC,
                                                          MonolayerVertexElement<0, 1>* pFaceAB,
                                                          MonolayerVertexElement<0, 1>* pFaceBC,
                                                          MonolayerVertexElement<0, 1>* pFaceAC)
{
    EXCEPTION("Void removals only implemented in three dimension, for three-dimensional monolayers.");
}
template <>
void MutableMonolayerVertexMesh<1, 2>::PerformVoidRemoval(Node<2>* pBasalNodeA, Node<2>* pBasalNodeB, Node<2>* pBasalNodeC,
                                                          MonolayerVertexElement<0, 2>* pFaceAB,
                                                          MonolayerVertexElement<0, 2>* pFaceBC,
                                                          MonolayerVertexElement<0, 2>* pFaceAC)
{
    EXCEPTION("Void removals only implemented in three dimension, for three-dimensional monolayers.");
}
template <>
void MutableMonolayerVertexMesh<1, 3>::PerformVoidRemoval(Node<3>* pBasalNodeA, Node<3>* pBasalNodeB, Node<3>* pBasalNodeC,
                                                          MonolayerVertexElement<0, 3>* pFaceAB,
                                                          MonolayerVertexElement<0, 3>* pFaceBC,
                                                          MonolayerVertexElement<0, 3>* pFaceAC)
{
    EXCEPTION("Void removals only implemented in three dimension, for three-dimensional monolayers.");
}
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformVoidRemoval(Node<SPACE_DIM>* pBasalNodeA, Node<SPACE_DIM>* pBasalNodeB, Node<SPACE_DIM>* pBasalNodeC,
                                                                            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFaceAB,
                                                                            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFaceBC,
                                                                            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFaceAC)
{
    std::cout << "We reach this!";
    // Calculate void centroid
    c_vector<double, SPACE_DIM> basal_nodes_midpoint = pBasalNodeA->rGetLocation()
        + this->GetVectorFromAtoB(pBasalNodeA->rGetLocation(), pBasalNodeB->rGetLocation()) / 3.0
        + this->GetVectorFromAtoB(pBasalNodeA->rGetLocation(), pBasalNodeC->rGetLocation()) / 3.0;

    Node<SPACE_DIM>* pApicalNodeA = nullptr;
    Node<SPACE_DIM>* pApicalNodeB = nullptr;
    Node<SPACE_DIM>* pApicalNodeC = nullptr;
    Node<SPACE_DIM>* p_temp_node = nullptr;
    // Find the corresponding apical nodes, i.e. the lateral edge connectes A with A
    // AB
    unsigned num_nodes_AB = pFaceAB->GetNumNodes();
    bool nodes_AB_ordered_AB = false;
    for (unsigned index_node = 0; index_node < num_nodes_AB; index_node++)
    {
        Node<SPACE_DIM>* p_current_node = pFaceAB->GetNode(index_node);
        Node<SPACE_DIM>* p_next_node = pFaceAB->GetNode((index_node + 1) % num_nodes_AB);
        if (pFaceAB->GetNodeType(index_node) == MonolayerVertexElementType::Basal && pFaceAB->GetNodeType((index_node + 1) % num_nodes_AB) == MonolayerVertexElementType::Basal)
        {
            if (p_current_node == pBasalNodeA)
            {
                nodes_AB_ordered_AB = true;
            }
            else
            {
                nodes_AB_ordered_AB = false;
            }
        }
        else if (pFaceAB->GetNodeType(index_node) == MonolayerVertexElementType::Apical && pFaceAB->GetNodeType((index_node + 1) % num_nodes_AB) == MonolayerVertexElementType::Apical)
        {
            // Opposite order on apical side
            pApicalNodeB = p_current_node;
            pApicalNodeA = p_next_node;
        }
    }
    if (!(nodes_AB_ordered_AB))
    {
        p_temp_node = pApicalNodeA;
        pApicalNodeA = pApicalNodeB;
        pApicalNodeB = p_temp_node;
    }

    // AC
    unsigned num_nodes_AC = pFaceAC->GetNumNodes();
    bool nodes_AC_ordered_AC = false;
    for (unsigned index_node = 0; index_node < num_nodes_AC; index_node++)
    {
        Node<SPACE_DIM>* p_current_node = pFaceAC->GetNode(index_node);
        Node<SPACE_DIM>* p_next_node = pFaceAC->GetNode((index_node + 1) % num_nodes_AC);
        if (pFaceAC->GetNodeType(index_node) == MonolayerVertexElementType::Basal && pFaceAC->GetNodeType((index_node + 1) % num_nodes_AC) == MonolayerVertexElementType::Basal)
        {
            if (p_current_node == pBasalNodeA)
            {
                nodes_AC_ordered_AC = true;
            }
            else
            {
                nodes_AC_ordered_AC = false;
            }
        }
        else if (pFaceAC->GetNodeType(index_node) == MonolayerVertexElementType::Apical && pFaceAC->GetNodeType((index_node + 1) % num_nodes_AC) == MonolayerVertexElementType::Apical)
        {
            // Opposite order on apical side
            pApicalNodeC = p_current_node;
            p_temp_node = p_next_node;
        }
    }
    if (!(nodes_AC_ordered_AC))
    {
        pApicalNodeC = p_temp_node;
    }

    c_vector<double, SPACE_DIM> apical_nodes_midpoint = pApicalNodeA->rGetLocation()
        + this->GetVectorFromAtoB(pApicalNodeA->rGetLocation(), pApicalNodeB->rGetLocation()) / 3.0
        + this->GetVectorFromAtoB(pApicalNodeA->rGetLocation(), pApicalNodeC->rGetLocation()) / 3.0;

    /*
     * In two steps, merge nodes A, B and C into a single node.  This is implemented in such a way that
     * the ordering of their indices does not matter.
     */

    PerformNodeMerge(pBasalNodeA, pBasalNodeB, pFaceAB);

    Node<SPACE_DIM>* p_merged_node_basal = pBasalNodeB;
    Node<SPACE_DIM>* p_merged_node_apical = pApicalNodeB;

    if (pBasalNodeB->IsDeleted())
    {
        p_merged_node_basal = pBasalNodeA;
        p_merged_node_apical = pApicalNodeA;
    }
    // We now have the face pegedNode <-> NodeC double, pFaceAC and pFaceBC
    PerformNodeMerge(pBasalNodeC, p_merged_node_basal, pFaceAC);

    if (p_merged_node_basal->IsDeleted())
    {
        p_merged_node_basal = pBasalNodeC;
        p_merged_node_apical = pApicalNodeC;
    }

    // We need to remove the last face.
    std::set<unsigned> basal_merged_node_elem_indices = p_merged_node_basal->rGetContainingElementIndices();

    for (std::set<unsigned>::const_iterator it = basal_merged_node_elem_indices.begin();
         it != basal_merged_node_elem_indices.end();
         ++it)
    {
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_element = this->GetElement(*it);
        unsigned local_face_index = p_temp_element->GetFaceLocalIndex(pFaceBC);
        if (local_face_index != UINT_MAX)
        {
            p_temp_element->DeleteFace(local_face_index);
            this->UpdateElementsFacesMapOfElement(*it);
            break; // The face should be only once in an element
        }
    }
    assert(!(pFaceBC->IsDeleted()));
    pFaceBC->MarkAsDeleted();
    mDeletedFaceIndices.push_back(pFaceBC->GetIndex());

    p_merged_node_basal->rGetModifiableLocation() = basal_nodes_midpoint;
    p_merged_node_apical->rGetModifiableLocation() = apical_nodes_midpoint;

    // Tag remaining node as non-boundary
    p_merged_node_basal->SetAsBoundaryNode(false);
    p_merged_node_apical->SetAsBoundaryNode(false);

    // Remove the deleted nodes and re-index
    RemoveDeletedNodes();
    // RemoveDeletedFaces();
}

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template <>
void MutableMonolayerVertexMesh<1, 1>::CheckForRosettes()
{
    EXCEPTION("CheckForRosettes only implemented in three dimension, for three-dimensional monolayers.");
}

template <>
void MutableMonolayerVertexMesh<1, 2>::CheckForRosettes()
{
    EXCEPTION("CheckForRosettes only implemented in three dimension, for three-dimensional monolayers.");
}

template <>
void MutableMonolayerVertexMesh<1, 3>::CheckForRosettes()
{
    EXCEPTION("CheckForRosettes only implemented in three dimension, for three-dimensional monolayers.");
}
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::CheckForRosettes()
{
    /**
     * First, we loop over each node and populate vectors of protorosette and rosette nodes which need to undergo
     * resolution for basal nodes
     *
     * We do not perform the resolution events in this initial loop because the resolution events involve changing
     * nodes in the mesh.
     */

    // Sets to store the nodes that need resolution events
    // std::set<Node<SPACE_DIM>* > protorosette_nodes;
    // std::set<Node<SPACE_DIM>* > rosette_nodes;

    // Loop over lateral faces, because then we can access the basal nodes only
    /*
    for(typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator face_iter = this->GetFaceIteratorBegin();
   face_iter != this->GetFaceIteratorEnd();
   ++face_iter)
    {
                    // Only check lateral faces
                    if(face_iter->GetFaceType() != MonolayerVertexElementType::Lateral)
                    {
                            continue;
                    }

                    unsigned num_nodes = face_iter->GetNumNodes();
    assert(num_nodes > 0);

                    Node<SPACE_DIM>* p_basal_node_a = nullptr;
                    Node<SPACE_DIM>* p_basal_node_b = nullptr;

                    // Caveat: We assume that we only have one apical and one basal edge
                    for(unsigned local_index=0; local_index < num_nodes; local_index++)
                    {
                            Node<SPACE_DIM>* p_current_node = face_iter->GetNode(local_index);
                            Node<SPACE_DIM>* p_next_node = face_iter->GetNode((local_index+1)%num_nodes);

                            MonolayerVertexElementType current_node_type = face_iter->GetNodeType(local_index);
                            MonolayerVertexElementType next_node_type = face_iter->GetNodeType((local_index+1)%num_nodes);

                            // Check if neighbouring nodes are either both apical or basal and too short
                            if(current_node_type == MonolayerVertexElementType::Basal && next_node_type == MonolayerVertexElementType::Basal)
                            {
                                    p_basal_node_a = p_current_node;
                                    p_basal_node_b = p_next_node;
                            }
                    }

                    assert(p_basal_node_a!=nullptr && p_basal_node_b!=nullptr);
  unsigned node_rank_a = p_basal_node_a->rGetContainingElementIndices().size();
                    unsigned node_rank_b = p_basal_node_b->rGetContainingElementIndices().size();
                    // A
  if (node_rank_a < 4)
  {
      // Nothing to do if the node is not high-rank
      continue;
  }
  else if (node_rank_a == 4)
  {
      // For protorosette nodes, we check against a random number to decide if resolution is necessary
      if (mProtorosetteResolutionProbabilityPerTimestep >= RandomNumberGenerator::Instance()->ranf())
      {
          protorosette_nodes.insert(p_basal_node_a);
      }
  }
  else // if (node_rank_a > 4)
  {
      // For rosette nodes, we check against a random number to decide if resolution is necessary
      if (mRosetteResolutionProbabilityPerTimestep >= RandomNumberGenerator::Instance()->ranf())
      {
          rosette_nodes.insert(p_basal_node_a);
      }
  }
                    // B
  if (node_rank_b < 4)
  {
      // Nothing to do if the node is not high-rank
      continue;
  }
  else if (node_rank_b == 4)
  {
      // For protorosette nodes, we check against a random number to decide if resolution is necessary
      if (mProtorosetteResolutionProbabilityPerTimestep >= RandomNumberGenerator::Instance()->ranf())
      {
          protorosette_nodes.insert(p_basal_node_b);
      }
  }
  else // if (node_rank_b > 4)
  {
      // For rosette nodes, we check against a random number to decide if resolution is necessary
      if (mRosetteResolutionProbabilityPerTimestep >= RandomNumberGenerator::Instance()->ranf())
      {
          rosette_nodes.insert(p_basal_node_b);
      }
  }
            }
            */

    /**
     * First, we loop over each node and populate vectors of protorosette and rosette nodes which need to undergo
     * resolution for basal nodes
     *
     * We do not perform the resolution events in this initial loop because the resolution events involve changing
     * the MapOfProtorosettes
     */

    // Sets to store the nodes that need resolution events
    std::set<Node<SPACE_DIM>*> protorosette_nodes;

    for (auto it = this->mMapOfProtorosettes.begin();
         it != this->mMapOfProtorosettes.end(); ++it)
    {
        Node<SPACE_DIM>* current_node = it->first;

        // Verify that node has not been marked for deletion, and that it is still contained in four elements
        assert(!(current_node->IsDeleted()));
        assert(current_node->rGetContainingElementIndices().size() == 4);

        // Perform protorosette resolution
        protorosette_nodes.insert(current_node);
    }

    /**
     * Finally, we loop over the contents of each node set and perform the necessary resolution events.
     *
     * Because each resolution event changes nodes, we include several assertions to catch possible unconsidered
     * behaviour.
     */

    for (auto it = protorosette_nodes.begin(); it != protorosette_nodes.end(); ++it)
    {
        Node<SPACE_DIM>* current_node = *it;
        assert(!(current_node->IsDeleted()));
        assert(current_node->rGetContainingElementIndices().size() == 4);
        this->PerformProtorosetteResolution(current_node);
    }

    /*
    // Finally, resolve any rosettes. This should not happen as we do not allow for the formation (yet)
    for (typename std::set<Node<SPACE_DIM>* >::iterator it = rosette_nodes.begin(); it != rosette_nodes.end(); it++)
    {
                                EXCEPTION("We need to resolve rosettes with more than 4 edges. This should not happen!");
        Node<SPACE_DIM>* current_node = *it;

        // Verify that node has not been marked for deletion, and that it is still contained in at least four elements
        assert( !(current_node->IsDeleted()) );
        assert( current_node->rGetContainingElementIndices().size() > 4 );

        // Perform protorosette resolution
        //this->PerformRosetteRankDecrease(current_node);
    }
                */
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::AddToMapOfProtorosettes(Node<SPACE_DIM>* pBasalNode, Node<SPACE_DIM>* pApicalNode,
                                                                                 MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pNeighbourElementA,
                                                                                 MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pNeighbourElementB)
{
    MutableMonolayerVertexMeshProtorosette<ELEMENT_DIM, SPACE_DIM> temp_protorosette{ pApicalNode, pNeighbourElementA, pNeighbourElementB };
    mMapOfProtorosettes[pBasalNode] = temp_protorosette;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MutableMonolayerVertexMeshProtorosette<ELEMENT_DIM, SPACE_DIM> MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetMappedProtorosette(Node<SPACE_DIM>* pBasalNode)
{
    assert(mMapOfProtorosettes.find(pBasalNode) != mMapOfProtorosettes.end());
    return mMapOfProtorosettes[pBasalNode];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::RemoveMappedProtorosette(Node<SPACE_DIM>* pBasalNode)
{
    assert(mMapOfProtorosettes.find(pBasalNode) != mMapOfProtorosettes.end());
    mMapOfProtorosettes.erase(pBasalNode);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::AddToSetOfElementsWithProtorosette(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    mSetOfElementsWithProtorosette.insert(pElement);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::RemoveFromSetOfElementsWithProtorosette(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    assert(mSetOfElementsWithProtorosette.find(pElement) != mSetOfElementsWithProtorosette.end());
    mSetOfElementsWithProtorosette.erase(pElement);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::DoesElementHaveProtorosette(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement)
{
    bool found = false;
    unsigned num_nodes = pElement->GetNumNodes();
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        if (pElement->GetNodeType(node_index) == MonolayerVertexElementType::Basal)
        {
            found = (mMapOfProtorosettes.find(pElement->GetNode(node_index)) != mMapOfProtorosettes.end());
        }
        if (found)
            return found;
    }
    return found;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::ClearSetOFElementsWithProtorosette()
{
    mSetOfElementsWithProtorosette.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::ClearElementsThatUnderwentT1Transitions()
{
    mElementsThatUnderwentT1Transitions.clear();
}

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template <>
void MutableMonolayerVertexMesh<1, 1>::PerformProtorosetteResolution(Node<1>* pProtorosetteNode)
{
    EXCEPTION("PerformProtorosetteResolution only implemented in three dimension, for three-dimensional monolayers.");
}

template <>
void MutableMonolayerVertexMesh<1, 2>::PerformProtorosetteResolution(Node<2>* pProtorosetteNode)
{
    EXCEPTION("PerformProtorosetteResolution only implemented in three dimension, for three-dimensional monolayers.");
}

template <>
void MutableMonolayerVertexMesh<1, 3>::PerformProtorosetteResolution(Node<3>* pProtorosetteNode)
{
    EXCEPTION("PerformProtorosetteResolution only implemented in three dimension, for three-dimensional monolayers.");
}
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::PerformProtorosetteResolution(Node<SPACE_DIM>* pProtorosetteNode)
{
    // Double check we are dealing with a protorosette
    assert(pProtorosetteNode->rGetContainingElementIndices().size() == 4);

    // We assume that we have a mapped protorosette. This might change after/if we implemented rosettes
    MutableMonolayerVertexMeshProtorosette<ELEMENT_DIM, SPACE_DIM> temp_protorosette = this->GetMappedProtorosette(pProtorosetteNode);

    Node<SPACE_DIM>* pProtorosetteNodeApical = temp_protorosette.pApicalNode;
    // Numbered elements which previously shared a face
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElementOne = temp_protorosette.pNeighbourOne;
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElementTwo = temp_protorosette.pNeighbourTwo;

    // Find global indices of elements around the protorosette node
    std::set<unsigned> protorosette_node_containing_elem_indices = pProtorosetteNode->rGetContainingElementIndices();

    // Determine the lettered elements which will share tthe new face
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElementA = nullptr;
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElementB = nullptr;
    for (std::set<unsigned>::const_iterator elem_index_iter = protorosette_node_containing_elem_indices.begin();
         elem_index_iter != protorosette_node_containing_elem_indices.end();
         elem_index_iter++)
    {
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(*elem_index_iter);
        if (p_element != pElementOne && p_element != pElementTwo)
        {
            if (pElementA == nullptr)
            {
                pElementA = p_element;
            }
            else
            {
                pElementB = p_element;
                break;
            }
        }
    }
    assert(pElementA != nullptr && pElementB != nullptr);

    /**
     * Ordering elements as follows, where A<->B are unordered:
     *
     *      \  1  /
     *       \   /
     *        \ /
     *     A   X   B
     *        / \
     *       /   \
     *      /  2  \
     *
     * Element A is the randomly chosen element from {A,B}.
     * These elements will share the new face, while the two elements 1 and
     * 2 which are adjacent to A and B will end up separated by it:
     *    basal network:
     *      \  1  /
     *       \   /
     *        \a/
     *         |
     *    A    |    B
     *         |
     *        /b\
     *       /   \
     *      /  2  \
     *
     * We assume that the apical and basal faces are oriented with normals pointing outwards.
     */

    /**
     * Next, we compute where to place the four nodes which will replace the two protorosette nodes.
     *
     * We place each node along the line joining the protorosette node to the centroid of element which will contain it,
     * and the distance along this line is such that the distance of nodes is larger than the swap distance.
     *
     * To do this, we will move the existing protorosette nodes in to element 1, and create a new node in element 2.  We
     * then need to tidy up the nodes by adding the new node to elements A, B and 2, and removing the protorosette node
     * from element 1.
     *
     */

    double swap_distance_final = (this->mCellRearrangementRatio) * (this->mCellRearrangementThreshold);

    c_vector<double, SPACE_DIM> basal_node_to_elem_1_basal_centre = this->GetPassiveCenterOfFaceTypeInElement(pElementOne->GetIndex(), MonolayerVertexElementType::Basal) - pProtorosetteNode->rGetLocation();
    c_vector<double, SPACE_DIM> basal_node_to_elem_2_basal_centre = this->GetPassiveCenterOfFaceTypeInElement(pElementTwo->GetIndex(), MonolayerVertexElementType::Basal) - pProtorosetteNode->rGetLocation();
    basal_node_to_elem_1_basal_centre /= norm_2(basal_node_to_elem_1_basal_centre);
    basal_node_to_elem_2_basal_centre /= norm_2(basal_node_to_elem_2_basal_centre);

    c_vector<double, SPACE_DIM> apical_node_to_elem_1_apical_centre = this->GetPassiveCenterOfFaceTypeInElement(pElementOne->GetIndex(), MonolayerVertexElementType::Apical) - pProtorosetteNodeApical->rGetLocation();
    c_vector<double, SPACE_DIM> apical_node_to_elem_2_apical_centre = this->GetPassiveCenterOfFaceTypeInElement(pElementTwo->GetIndex(), MonolayerVertexElementType::Apical) - pProtorosetteNodeApical->rGetLocation();
    apical_node_to_elem_1_apical_centre /= norm_2(apical_node_to_elem_1_apical_centre);
    apical_node_to_elem_2_apical_centre /= norm_2(apical_node_to_elem_2_apical_centre);

    double angle_factor_basal = 1.0 / sqrt(0.5 - 0.5 * inner_prod(basal_node_to_elem_1_basal_centre, basal_node_to_elem_2_basal_centre));
    double swap_distance_basal = swap_distance_final * angle_factor_basal;

    double angle_factor_apical = 1.0 / sqrt(0.5 - 0.5 * inner_prod(apical_node_to_elem_1_apical_centre, apical_node_to_elem_2_apical_centre));
    double swap_distance_apical = swap_distance_final * angle_factor_apical;

    // Calculate new node locations
    c_vector<double, SPACE_DIM> new_location_of_protorosette_node_basal = pProtorosetteNode->rGetLocation() + (0.5 * swap_distance_basal) * basal_node_to_elem_1_basal_centre;
    c_vector<double, SPACE_DIM> location_of_new_node_basal = pProtorosetteNode->rGetLocation() + (0.5 * swap_distance_basal) * basal_node_to_elem_2_basal_centre;

    c_vector<double, SPACE_DIM> new_location_of_protorosette_node_apical = pProtorosetteNodeApical->rGetLocation() + (0.5 * swap_distance_apical) * apical_node_to_elem_1_apical_centre;
    c_vector<double, SPACE_DIM> location_of_new_node_apical = pProtorosetteNodeApical->rGetLocation() + (0.5 * swap_distance_apical) * apical_node_to_elem_2_apical_centre;

    // Move protorosette nodes to new location
    pProtorosetteNode->rGetModifiableLocation() = new_location_of_protorosette_node_basal;
    pProtorosetteNodeApical->rGetModifiableLocation() = new_location_of_protorosette_node_apical;

    // Create new nodes in correct location
    unsigned new_node_global_index_basal = this->AddNode(new Node<SPACE_DIM>(this->GetNumNodes(), location_of_new_node_basal, false));
    Node<SPACE_DIM>* p_new_node_basal = this->GetNode(new_node_global_index_basal);

    unsigned new_node_global_index_apical = this->AddNode(new Node<SPACE_DIM>(this->GetNumNodes(), location_of_new_node_apical, false));
    Node<SPACE_DIM>* p_new_node_apical = this->GetNode(new_node_global_index_apical);

    // Create the new face and add it to mesh
    std::vector<Node<SPACE_DIM>*> nodes_new_face = { pProtorosetteNode, pProtorosetteNodeApical, p_new_node_apical, p_new_node_basal };
    std::vector<MonolayerVertexElementType> node_types_new_face = { MonolayerVertexElementType::Basal, MonolayerVertexElementType::Apical,
                                                                    MonolayerVertexElementType::Apical, MonolayerVertexElementType::Basal };
    unsigned new_face_global_index = this->AddFace(new MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>(this->GetNumFaces(), MonolayerVertexElementType::Lateral, nodes_new_face, node_types_new_face));
    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_new_lateral_face = this->GetFace(new_face_global_index);

    /*
     * We now rename the elements A and B such that A is positively oriented and B negatively oriented w.r.t. face
     *    both sides:
     *      \  1  /
     *       \   /
     *        \ /old
     *         |
     *    A    |    B
     *         |
     *        / \new
     *       /   \
     *      /  2  \
     */

    c_vector<double, SPACE_DIM> normal_face;
    this->CalculateUnitNormalToFaceWithArea(p_new_lateral_face, normal_face);
    unsigned global_index_element_A = pElementA->GetIndex();
    c_vector<double, SPACE_DIM> face_to_center_A = this->GetVectorFromAtoB(this->GetPassiveCenterOfFace(p_new_lateral_face),
                                                                           this->GetCentroidOfElement(global_index_element_A));

    if (inner_prod(normal_face, face_to_center_A) > 0.0) // face oriented inwards for A
    {
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element_temp = pElementA;
        pElementA = pElementB;
        pElementB = p_element_temp;
    }

    // Now we add the nodes to the elements
    pElementA->AddNode(p_new_node_basal, pElementA->GetNumNodes() - 1, MonolayerVertexElementType::Basal);
    pElementA->AddNode(p_new_node_apical, pElementA->GetNumNodes() - 1, MonolayerVertexElementType::Apical);

    pElementB->AddNode(p_new_node_basal, pElementB->GetNumNodes() - 1, MonolayerVertexElementType::Basal);
    pElementB->AddNode(p_new_node_apical, pElementB->GetNumNodes() - 1, MonolayerVertexElementType::Apical);

    pElementTwo->AddNode(p_new_node_basal, pElementTwo->GetNumNodes() - 1, MonolayerVertexElementType::Basal);
    pElementTwo->AddNode(p_new_node_apical, pElementTwo->GetNumNodes() - 1, MonolayerVertexElementType::Apical);

    // Add face to elements
    pElementA->AddFace(p_new_lateral_face, MonolayerVertexElementType::Lateral, false);
    // false because clockwise orientation == false
    pElementB->AddFace(p_new_lateral_face, MonolayerVertexElementType::Lateral, true);

    // Remove old nodes from element 2
    pElementTwo->DeleteNode(pElementTwo->GetNodeLocalIndex(pProtorosetteNode->GetIndex()));
    pElementTwo->DeleteNode(pElementTwo->GetNodeLocalIndex(pProtorosetteNodeApical->GetIndex()));

    // Change the nodes in faces for element 2
    // This also corrects for the wrong connections in lateral faces in elements A and B
    // because we loop over ALL faces which belong to element 2 (incl. these ones)
    for (unsigned face_index = 0; face_index < pElementTwo->GetNumFaces(); face_index++)
    {
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = pElementTwo->GetFace(face_index);

        unsigned local_index_basal = p_face->GetNodeLocalIndex(pProtorosetteNode->GetIndex());
        unsigned local_index_apical = p_face->GetNodeLocalIndex(pProtorosetteNodeApical->GetIndex());
        if (local_index_basal != UINT_MAX)
        {
            // Replace old node with new node on basal side
            p_face->UpdateNode(local_index_basal, p_new_node_basal, MonolayerVertexElementType::Basal);
        }
        if (local_index_apical != UINT_MAX)
        {
            // Replace old node with new node on apical side
            p_face->UpdateNode(local_index_apical, p_new_node_apical, MonolayerVertexElementType::Apical);
        }
    }

    // Add the new nodes to the faces of elements A and B
    // First A
    for (unsigned face_index = 0; face_index < pElementA->GetNumFaces(); face_index++)
    {
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = pElementA->GetFace(face_index);
        MonolayerVertexElementType face_type = p_face->GetFaceType();
        if (face_type == MonolayerVertexElementType::Apical)
        {
            unsigned local_index_old_apical = p_face->GetNodeLocalIndex(pProtorosetteNodeApical->GetIndex());
            assert(local_index_old_apical != UINT_MAX);

            // We have to add the new apical element before the old one on the (positively oriented element's) apical face
            unsigned local_index_new_apical = (local_index_old_apical - 1 + p_face->GetNumNodes()) % p_face->GetNumNodes();
            p_face->AddNode(p_new_node_apical, local_index_new_apical, MonolayerVertexElementType::Apical);
        }
        else if (face_type == MonolayerVertexElementType::Basal)
        {
            unsigned local_index_old_basal = p_face->GetNodeLocalIndex(pProtorosetteNode->GetIndex());
            assert(local_index_old_basal != UINT_MAX);

            // We have to add the new basal element after the old one on the (positively oriented element's) basal face
            unsigned local_index_new_basal = local_index_old_basal;
            p_face->AddNode(p_new_node_basal, local_index_new_basal, MonolayerVertexElementType::Apical);
        }
    }
    // Now B
    for (unsigned face_index = 0; face_index < pElementB->GetNumFaces(); face_index++)
    {
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = pElementB->GetFace(face_index);
        MonolayerVertexElementType face_type = p_face->GetFaceType();
        if (face_type == MonolayerVertexElementType::Apical)
        {
            unsigned local_index_old_apical = p_face->GetNodeLocalIndex(pProtorosetteNodeApical->GetIndex());
            assert(local_index_old_apical != UINT_MAX);

            // We have to add the new apical element after the old one on the (negatively oriented element's) apical face
            unsigned local_index_new_apical = local_index_old_apical;
            p_face->AddNode(p_new_node_apical, local_index_new_apical, MonolayerVertexElementType::Apical);
        }
        else if (face_type == MonolayerVertexElementType::Basal)
        {
            unsigned local_index_old_basal = p_face->GetNodeLocalIndex(pProtorosetteNode->GetIndex());
            assert(local_index_old_basal != UINT_MAX);

            // We have to add the new basal element before the old one on the (negatively oriented element's) basal face
            unsigned local_index_new_basal = (local_index_old_basal - 1 + p_face->GetNumNodes()) % p_face->GetNumNodes();
            p_face->AddNode(p_new_node_basal, local_index_new_basal, MonolayerVertexElementType::Apical);
        }
    }

    // For a matter of completeness, remove the entry from the map of protorosettes
    this->RemoveMappedProtorosette(pProtorosetteNode);

    // And delete the elements from the set of elements with protorosettes
    // this->RemoveFromSetOfElementsWithProtorosette(pElementA);
    // this->RemoveFromSetOfElementsWithProtorosette(pElementB);
    // this->RemoveFromSetOfElementsWithProtorosette(pElementOne);
    // this->RemoveFromSetOfElementsWithProtorosette(pElementTwo);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::CalculateActiveT1SwapProbability(double apical_length, double basal_length)
{
    if (mLengthDependentActiveT1BoltzmannParameter > 0.0)
    {
        double mid_plane_length = (apical_length + basal_length) / 2.0;

        return exp(-(mid_plane_length - mCellRearrangementThreshold) / (mLengthDependentActiveT1BoltzmannParameter));
    }

    // if isScalingWithEdges
    double edge_scaling_factor = mIsActiveT1ProbPerEdge ? 1.0 : 1.0 / this->GetNumFacesByType(MonolayerVertexElementType::Lateral);
    double swap_prob = mActiveT1SwapProbability * edge_scaling_factor;

    return swap_prob;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SetActiveT1BoltzmannParameter(double boltzmann_parameter)
{
    mLengthDependentActiveT1BoltzmannParameter = boltzmann_parameter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetAndResetPassiveT1Count()
{
    unsigned count = mPassiveT1TransitionsCounter;

    mPassiveT1TransitionsCounter = 0;

    return count;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MutableMonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetAndResetActiveT1Count()
{
    unsigned count = mActiveT1TransitionsCounter;

    mActiveT1TransitionsCounter = 0;

    return count;
}

// Explicit instantiation
template class MutableMonolayerVertexMesh<1, 1>;
template class MutableMonolayerVertexMesh<1, 2>;
template class MutableMonolayerVertexMesh<1, 3>;
template class MutableMonolayerVertexMesh<2, 2>;
template class MutableMonolayerVertexMesh<2, 3>;
template class MutableMonolayerVertexMesh<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MutableMonolayerVertexMesh)