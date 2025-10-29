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
#include "MonolayerVertexElement.hpp"
#include <algorithm>
#include <cassert>

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElement(unsigned index, MonolayerVertexElementType rFaceType,
                                                                       const std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*>& rFaces,
                                                                       const std::vector<bool>& rOrientations,
                                                                       const std::vector<Node<SPACE_DIM>*>& rNodes,
                                                                       const std::vector<MonolayerVertexElementType>& rNodeTypes)
        : MutableElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes),
          mFaces(rFaces),
          mOrientations(rOrientations),
          mFaceType(rFaceType),
          mNodeTypes(rNodeTypes),
          mBoundaryFace(false)
{
    // This constructor should only be used in 3D
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Each face must have an associated orientation
    assert(mFaces.size() == mOrientations.size());
    assert(this->mNodes.size() == mNodeTypes.size());

    if (SPACE_DIM == ELEMENT_DIM)
    {
        // Register element with nodes
        this->RegisterWithNodes();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElement(unsigned index, MonolayerVertexElementType rFaceType,
                                                                       const std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*>& rFaces,
                                                                       const std::vector<bool>& rOrientations)
        : MutableElement<ELEMENT_DIM, SPACE_DIM>(index),
          mFaces(rFaces),
          mOrientations(rOrientations),
          mFaceType(rFaceType),
          mBoundaryFace(false)
{
    // Each face must have an associated orientation
    assert(mFaces.size() == mOrientations.size());

    // Make a set of nodes and node types with mFaces
    std::set<std::pair<Node<SPACE_DIM>*, MonolayerVertexElementType> > nodes_set;
    for (unsigned face_index = 0; face_index < mFaces.size(); face_index++)
    {
        for (unsigned node_index = 0; node_index < mFaces[face_index]->GetNumNodes(); node_index++)
        {
            nodes_set.insert({ mFaces[face_index]->GetNode(node_index), mFaces[face_index]->GetNodeType(node_index) });
        }
    }

    // Populate mNodes and mNodeTypes
    for (typename std::set<std::pair<Node<SPACE_DIM>*, MonolayerVertexElementType> >::iterator node_iter = nodes_set.begin();
         node_iter != nodes_set.end();
         ++node_iter)
    {
        this->mNodes.push_back(node_iter->first);
        this->mNodeTypes.push_back(node_iter->second);
    }

    // Register element with nodes
    this->RegisterWithNodes();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElement(unsigned index, MonolayerVertexElementType rFaceType,
                                                                       const std::vector<Node<SPACE_DIM>*>& rNodes,
                                                                       const std::vector<MonolayerVertexElementType>& rNodeTypes)
        : MutableElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes),
          mFaceType(rFaceType),
          mNodeTypes(rNodeTypes),
          mBoundaryFace(false)
{
    assert(rNodes.size() == rNodeTypes.size());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElement(unsigned index)
        : MutableElement<ELEMENT_DIM, SPACE_DIM>(index)
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::~MonolayerVertexElement()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNode, const unsigned& rIndex, const MonolayerVertexElementType NodeType)
{
    /**
     * When constructing a VertexMesh as the Voronoi dual to a Delaunay mesh,
     * each MutableElement is initially constructed without nodes. We therefore
     * require the two cases below.
     */
    if (this->mNodes.empty())
    {
        // Populate mNodes with pNode
        this->mNodes.push_back(pNode);
        this->mNodeTypes.push_back(NodeType);

        // Add element to this node
        if (ELEMENT_DIM == SPACE_DIM)
        {
            this->mNodes[0]->AddElement(this->mIndex);
        }
    }
    else
    {
        assert(rIndex < this->mNodes.size());

        // Add pNode to rIndex+1 element of mNodes pushing the others up
        this->mNodes.insert(this->mNodes.begin() + rIndex + 1, pNode);
        this->mNodeTypes.insert(this->mNodeTypes.begin() + rIndex + 1, NodeType);

        // Add element to this node
        if (ELEMENT_DIM == SPACE_DIM)
        {
            this->mNodes[rIndex + 1]->AddElement(this->mIndex);
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNode, const unsigned& rIndex)
{
    this->AddNode(pNode, rIndex, MonolayerVertexElementType::Undetermined);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNode)
{
    this->AddNode(pNode, this->mNodes.size() - 1);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::DeleteNode(const unsigned& rIndex)
{
    assert(rIndex < this->mNodes.size());

    // Remove element from the node at this location
    // if this is a fully-dimensional object
    if (SPACE_DIM == ELEMENT_DIM)
    {
        this->mNodes[rIndex]->RemoveElement(this->mIndex);
    }
    // Remove the node at rIndex (removes node from element)
    this->mNodes.erase(this->mNodes.begin() + rIndex);
    this->mNodeTypes.erase(this->mNodeTypes.begin() + rIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode, MonolayerVertexElementType NodeType)
{
    assert(rIndex < this->mNodes.size());
    assert(rIndex < this->mNodeTypes.size());
    // if this is a fully-dimensional object we can use the standard function
    if (SPACE_DIM == ELEMENT_DIM)
    {
        MutableElement<ELEMENT_DIM, SPACE_DIM>::UpdateNode(rIndex, pNode);
    }
    else // if not we have to do it without registering the element
    {
        // Update the node at this location
        this->mNodes[rIndex] = pNode;
    }

    // Update the node type at this location
    this->mNodeTypes[rIndex] = NodeType;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode, MonolayerVertexElementType NodeType)
{
    assert(pOldNode != pNewNode);
    for (unsigned i = 0; i < this->mNodes.size(); i++)
    {
        if (this->mNodes[i] == pOldNode)
        {
            UpdateNode(i, pNewNode, NodeType);
            return;
        }
    }
    EXCEPTION("You didn't have that node to start with.");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::MarkAsDeleted()
{
    // Mark element as deleted
    this->mIsDeleted = true;

    // Update nodes in the element so they know they are not contained by it
    // if this is a fully-dimensional element and not a face
    if (SPACE_DIM == ELEMENT_DIM)
    {
        for (unsigned i = 0; i < this->GetNumNodes(); i++)
        {
            this->mNodes[i]->RemoveElement(this->mIndex);
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mFaces.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::AddFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, MonolayerVertexElementType FaceType, bool Orientation)
{
    // Add pFace to the end of mFaces
    this->mFaces.push_back(pFace);
    // pFace->mFaceType = FaceType;
    this->mOrientations.push_back(Orientation);

    // Create a set of indices of nodes currently owned by this element
    std::set<unsigned> node_indices;
    for (unsigned local_index = 0; local_index < this->GetNumNodes(); local_index++)
    {
        node_indices.insert(this->GetNodeGlobalIndex(local_index));
    }

    // Loop over nodes owned by pFace
    unsigned end_index = this->GetNumNodes() - 1;
    for (unsigned local_index = 0; local_index < pFace->GetNumNodes(); local_index++)
    {
        // If this node is not already owned by this element...
        unsigned global_index = pFace->GetNodeGlobalIndex(local_index);
        if (node_indices.find(global_index) == node_indices.end())
        {
            // ... then add it to the element (and vice versa)
            this->AddNode(pFace->GetNode(local_index), end_index, pFace->GetNodeType(local_index));
            end_index++;
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::DeleteFace(const unsigned& rIndex)
{
    assert(rIndex < this->mFaces.size());
    assert(rIndex < this->mOrientations.size());

    // Remove the face at rIndex (removes face element from element)
    this->mFaces.erase(this->mFaces.begin() + rIndex);
    this->mOrientations.erase(this->mOrientations.begin() + rIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::GetFace(unsigned index) const
{
    assert(index < mFaces.size());
    return mFaces[index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::GetFaceLocalIndex(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace) const
{
    unsigned local_index = UINT_MAX;
    for (unsigned i = 0; i < this->mFaces.size(); i++)
    {
        if (mFaces[i] == pFace)
        {
            local_index = i;
            break;
        }
    }
    return local_index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElementType MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::GetFaceType() const
{
    return mFaceType;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElementType MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::GetNodeType(unsigned localIndex) const
{
    assert(localIndex < mNodeTypes.size());
    return mNodeTypes[localIndex];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::FaceIsOrientatedClockwise(unsigned index) const
{
    assert(index < mOrientations.size());
    return mOrientations[index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::SwitchOrientation()
{
    std::reverse(this->mNodeTypes.begin(), this->mNodeTypes.end());
    std::reverse(this->mNodes.begin(), this->mNodes.end());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::IsBoundaryFace() const
{
    return mBoundaryFace;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>::SetAsBoundaryFace(bool boundaryFace)
{
    mBoundaryFace = boundaryFace;
}

//////////////////////////////////////////////////////////////////////
//                  Specialization for 1d elements                  //
//                                                                  //
//                 1d elements are just edges (lines)               //
//////////////////////////////////////////////////////////////////////

/**
 * Specialization for 1d elements so we don't get errors from Boost on some
 * compilers.
 */
template <unsigned SPACE_DIM>
MonolayerVertexElement<1, SPACE_DIM>::MonolayerVertexElement(unsigned index, MonolayerVertexElementType rFaceType, const std::vector<Node<SPACE_DIM>*>& rNodes,
                                                             const std::vector<MonolayerVertexElementType>& rNodeTypes)
        : MutableElement<1, SPACE_DIM>(index, rNodes),
          mNodeTypes(rNodeTypes),
          mFaceType(rFaceType)
{
}

template <unsigned SPACE_DIM>
void MonolayerVertexElement<1, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNode, const unsigned& rIndex, const MonolayerVertexElementType NodeType)
{
    /**
     * When constructing a VertexMesh as the Voronoi dual to a Delaunay mesh,
     * each MutableElement is initially constructed without nodes. We therefore
     * require the two cases below.
     */
    if (this->mNodes.empty())
    {
        // Populate mNodes with pNode
        this->mNodes.push_back(pNode);
        this->mNodeTypes.push_back(NodeType);

        // Add element to this node
        this->mNodes[0]->AddElement(this->mIndex);
    }
    else
    {
        assert(rIndex < this->mNodes.size());

        // Add pNode to rIndex+1 element of mNodes pushing the others up
        this->mNodes.insert(this->mNodes.begin() + rIndex + 1, pNode);
        this->mNodeTypes.insert(this->mNodeTypes.begin() + rIndex + 1, NodeType);

        // Add element to this node
        this->mNodes[rIndex + 1]->AddElement(this->mIndex);
    }
}

template <unsigned SPACE_DIM>
void MonolayerVertexElement<1, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNode, const unsigned& rIndex)
{
    this->AddNode(pNode, rIndex, MonolayerVertexElementType::Undetermined);
}

template <unsigned SPACE_DIM>
void MonolayerVertexElement<1, SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNode)
{
    this->AddNode(pNode, this->mNodes.size() - 1);
}

template <unsigned SPACE_DIM>
void MonolayerVertexElement<1, SPACE_DIM>::DeleteNode(const unsigned& rIndex)
{
    assert(rIndex < this->mNodes.size());

    // Remove element from the node at this location
    // if this is a fully-dimensional object
    if (SPACE_DIM == 1)
    {
        this->mNodes[rIndex]->RemoveElement(this->mIndex);
    }

    // Remove the node at rIndex (removes node from element)
    this->mNodes.erase(this->mNodes.begin() + rIndex);
    this->mNodeTypes.erase(this->mNodeTypes.begin() + rIndex);
}

template <unsigned SPACE_DIM>
void MonolayerVertexElement<1, SPACE_DIM>::UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode, MonolayerVertexElementType NodeType)
{
    assert(rIndex < this->mNodes.size());

    MutableElement<1, SPACE_DIM>::UpdateNode(rIndex, pNode);

    // Update the node type at this location
    this->mNodeTypes[rIndex] = NodeType;
}

template <unsigned SPACE_DIM>
void MonolayerVertexElement<1, SPACE_DIM>::MarkAsDeleted()
{
    // Mark element as deleted
    this->mIsDeleted = true;

    // Update nodes in the element so they know they are not contained by it
    // if this is a fully-dimensional element and not a face
    if (SPACE_DIM == 1)
    {
        for (unsigned i = 0; i < this->GetNumNodes(); i++)
        {
            this->mNodes[i]->RemoveElement(this->mIndex);
        }
    }
}

template <unsigned SPACE_DIM>
unsigned MonolayerVertexElement<1, SPACE_DIM>::GetNumFaces() const
{
    return 0;
}

template <unsigned SPACE_DIM>
MonolayerVertexElement<0, SPACE_DIM>* MonolayerVertexElement<1, SPACE_DIM>::GetFace(unsigned index) const
{
    return nullptr;
}

template <unsigned SPACE_DIM>
void MonolayerVertexElement<1, SPACE_DIM>::DeleteFace(const unsigned& rIndex)
{
}

template <unsigned SPACE_DIM>
MonolayerVertexElementType MonolayerVertexElement<1, SPACE_DIM>::GetFaceType() const
{
    return mFaceType;
}

template <unsigned SPACE_DIM>
MonolayerVertexElementType MonolayerVertexElement<1, SPACE_DIM>::GetNodeType(unsigned localIndex) const
{
    assert(localIndex < mNodeTypes.size());
    return mNodeTypes[localIndex];
}

template <unsigned SPACE_DIM>
bool MonolayerVertexElement<1, SPACE_DIM>::IsBoundaryFace() const
{
    return mBoundaryFace;
}

template <unsigned SPACE_DIM>
void MonolayerVertexElement<1, SPACE_DIM>::SetAsBoundaryFace(bool boundaryFace)
{
    mBoundaryFace = boundaryFace;
}

template <unsigned SPACE_DIM>
bool MonolayerVertexElement<1, SPACE_DIM>::FaceIsOrientatedClockwise(unsigned index) const
{
    return false;
}

// Explicit instantiation
template class MonolayerVertexElement<1, 1>;
template class MonolayerVertexElement<1, 2>;
template class MonolayerVertexElement<1, 3>;
template class MonolayerVertexElement<2, 2>;
template class MonolayerVertexElement<2, 3>;
template class MonolayerVertexElement<3, 3>;