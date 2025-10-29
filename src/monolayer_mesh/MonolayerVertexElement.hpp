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
#ifndef MONOLAYERVERTEXELEMENT_HPP_
#define MONOLAYERVERTEXELEMENT_HPP_

#include "MutableElement.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"

/**
 * Type of vertex as datatype, which allows 'Apical', 'Basal, 'Lateral' and 'Undefined'
 */
enum class MonolayerVertexElementType : unsigned
{
    Apical = 1,
    Basal = 2,
    Lateral = 3,
    Undetermined = 0
};

/**
 * An element class for use in the VertexMesh class. The main
 * difference between this and the Element class is that a
 * VertexElement can have a variable number of nodes associated
 * with it.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonolayerVertexElement : public MutableElement<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * Faces of the VertexElement, which should be distinct.
     */
    std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> mFaces;

    /**
     * How each face is oriented. From the perspective of the outside
     * of the element, the vertices of each face should be ordered
     * anti clockwise. If and only if this is false, the order of vertices
     * in the corresponding face should be reversed.
     *
     * N.B. Most faces belong to two elements, but with opposite
     * orientations. This allows us to reuse the face data across the
     * two cells.
     */
    std::vector<bool> mOrientations;

    /**
     * Face types of the MonolayerVertexElement, i.e. apical, basal, lateral or volume
     */
    MonolayerVertexElementType mFaceType;

    /**
     * Node types of the MonolayerVertexElement nodes, i.e. apical, basal or lateral
     */
    std::vector<MonolayerVertexElementType> mNodeTypes;

    /**
     * Whether this is a boundary element (esp. used for boundary faces)
     */
    bool mBoundaryFace;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        // This needs to be first so that MeshBasedCellPopulation::Validate() doesn't go mental.
        archive & mFaces;
        archive & mOrientations;
        archive & mNodeTypes;
        archive & mBoundaryFace;
        archive & mFaceType;
        archive& boost::serialization::base_object<MutableElement<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:
    using MutableElement<ELEMENT_DIM, SPACE_DIM>::UpdateNode;

    /**
     * Constructor.
     *
     * @param index global index of the element
     * @param rFaces vector of faces associated with the element
     * @param rOrientations vector of orientations of the faces associated with the element, where true means counter-clockwise from the outsid
     * @param rFaceTypes vector of face types assiciated with faces
     */
    MonolayerVertexElement(unsigned index, MonolayerVertexElementType FaceType,
                           const std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*>& rFaces,
                           const std::vector<bool>& rOrientations);

    /**
     * Constructor.
     *
     * @param index global index of the element
     * @param rNodes vector of Nodes associated with the element
     * @param rFaceTypes vector of face types assiciated with faces
     */
    MonolayerVertexElement(unsigned index, MonolayerVertexElementType FaceType,
                           const std::vector<Node<SPACE_DIM>*>& rNodes,
                           const std::vector<MonolayerVertexElementType>& rNodeTypes);

    /**
     * Constructor used to specify the element completely. This ensures that
     * the nodes and faces are owned by the element *in a specified order*.
     * See #1076 and #1377 for more details.
     *
     * @param index global index of the element
     * @param rFaces vector of faces associated with the element
     * @param rOrientations vector of orientations of the faces associated with the element, where true means counter-clockwise from the outside
     * @param rNodes vector of Nodes associated with the element
     * @param rFaceTypes vector of face types assiciated with faces
     */
    MonolayerVertexElement(unsigned index, MonolayerVertexElementType FaceType,
                           const std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*>& rFaces,
                           const std::vector<bool>& rOrientations,
                           const std::vector<Node<SPACE_DIM>*>& rNodes,
                           const std::vector<MonolayerVertexElementType>& rNodeTypes);

    /**
     *
     * Alternative constructor.
     *
     * When constructing a VertexMesh as the Voronoi dual to a Delaunay mesh,
     * each VertexElement is initially constructed without nodes.
     *
     * @param index global index of the element
     */
    MonolayerVertexElement(unsigned index);

    /**
     * Destructor.
     */
    ~MonolayerVertexElement();

    /**
     * Add a node to the element between nodes at rIndex and rIndex+1.
     *
     * @param rIndex the local index of the node after which the new node is added
     * @param pNode a pointer to the new node
     * @param NodeType node type of node to be added
     */
    void AddNode(Node<SPACE_DIM>* pNode, const unsigned& rIndex, const MonolayerVertexElementType NodeType);

    /**
     * Add a node to the element with standard class 113
     *
     * @param pNode a pointer to the new node
     */
    void AddNode(Node<SPACE_DIM>* pNode);

    /**
     * Add a node to the element between nodes at rIndex and rIndex+1 of standard type 113
     *
     * @param rIndex the local index of the node after which the new node is added
     * @param pNode a pointer to the new node
     */
    void AddNode(Node<SPACE_DIM>* pNode, const unsigned& rIndex);

    /**
     * Delete a node with given local index.
     *
     * @param rIndex is the local index of the node to remove
     */
    void DeleteNode(const unsigned& rIndex);

    /**
     * Update node at the given index with new node type
     *
     * @param rIndex is an local index to which node to change
     * @param pNode is a pointer to the replacement node
     * @param NodeType
     */
    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode, MonolayerVertexElementType NodeType);

    /**
     * Replace node with new node and corresponding node type
     *
     * @param pOldNode is a pointer to the old node to be repplaced
     * @param pNewNode is a pointer to the replacement node
     * @param NodeType
     */
    void ReplaceNode(Node<SPACE_DIM>* pOldNode, Node<SPACE_DIM>* pNewNode, MonolayerVertexElementType NodeType);

    /**
     * Overridden MarkAsDeleted() method.
     *
     * Mark an element as having been removed from the mesh.
     * Also notify nodes in the element that it has been removed
     * if this node is fully-dimensional
     */
    void MarkAsDeleted();

    /**
     * @return the number of faces owned by this element.
     */
    unsigned GetNumFaces() const;

    /**
     * @return whether this is a boundary face.
     */
    bool IsBoundaryFace() const;

    /**
     * @param boundaryFace whether this is a boundary face.
     */
    void SetAsBoundaryFace(bool boundaryFace = true);

    /**
     * Add a face to the element.
     *
     * @param pFace a pointer to the new face
     * @param FaceType type of the new face
     * @param Orientation orientation of the face w.r.t. the element
     */
    void AddFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, MonolayerVertexElementType FaceType, bool Orientation);

    /**
     * Delete a face with given local index from the element.
     *
     * @param rIndex is the local index of the face to remove
     */
    void DeleteFace(const unsigned& rIndex);

    /**
     * @param index the global index of a specified face
     *
     * @return a pointer to the face
     */
    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* GetFace(unsigned index) const;

    /**
     * Determine the local index of the face in the element. If face not in element return UINT_MAX.
     *
     * @param pointer to the face
     *
     * @return local index of the face
     */
    unsigned GetFaceLocalIndex(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace) const;

    /**
     *
     * @return number of face type of the vertex
     */
    MonolayerVertexElementType GetFaceType() const;

    /**
     * @param index the local node index of the node
     *
     * @return type of node
     */
    MonolayerVertexElementType GetNodeType(unsigned localIndex) const;

    /**
     * @return whether the face with a given index is oriented clockwise.
     *
     * @param index the index of the face
     */
    bool FaceIsOrientatedClockwise(unsigned index) const;

    /**
     * Switch the orientation of the element. This is done via changing the local
     * ordering of the nodes and node types in the internal data structures.
     * This should only be used sparsly as this might be costly.
     *
     * Note: This does NOT change the value for orientation in elements which
     * contain this as a face!
     */
    void SwitchOrientation();
};

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
class MonolayerVertexElement<1, SPACE_DIM> : public MutableElement<1, SPACE_DIM>
{
private:
    /**
     * Node types of the MonolayerVertexElement nodes, i.e. apical, basal or lateral
     */
    std::vector<MonolayerVertexElementType> mNodeTypes;

    /**
     * Face types of the MonolayerVertexElement, i.e. apical, basal, lateral or volume
     */
    MonolayerVertexElementType mFaceType;

    /**
     * Whether this is a boundary element (esp. used for boundary faces)
     */
    bool mBoundaryFace;

public:
    using MutableElement<1, SPACE_DIM>::UpdateNode;

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param FaceType Face type
     * @param index  the index of the element in the mesh
     * @param rNodes the nodes owned by the element
     */
    MonolayerVertexElement(unsigned index, MonolayerVertexElementType FaceType, const std::vector<Node<SPACE_DIM>*>& rNodes, const std::vector<MonolayerVertexElementType>& rNodeTypes);

    /**
     * Add a node to the element between nodes at rIndex and rIndex+1.
     *
     * @param rIndex the local index of the node after which the new node is added
     * @param pNode a pointer to the new node
     * @param NodeType node type of node to be added
     */
    void AddNode(Node<SPACE_DIM>* pNode, const unsigned& rIndex, const MonolayerVertexElementType NodeType);

    /**
     * Add a node to the element of Undetermined type
     *
     * @param pNode a pointer to the new node
     */
    void AddNode(Node<SPACE_DIM>* pNode);

    /**
     * Delete a node with given local index.
     *
     * @param rIndex is the local index of the node to remove
     */
    void DeleteNode(const unsigned& rIndex);

    /**
     * Add a node to the element between nodes at rIndex and rIndex+1 of Undetermined type
     *
     * @param rIndex the local index of the node after which the new node is added
     * @param pNode a pointer to the new node
     */
    void AddNode(Node<SPACE_DIM>* pNode, const unsigned& rIndex);

    /**
     * Update node at the given index with new node type
     *
     * @param rIndex is an local index to which node to change
     * @param pNode is a pointer to the replacement node
     * @param NodeType
     */
    void UpdateNode(const unsigned& rIndex, Node<SPACE_DIM>* pNode, MonolayerVertexElementType NodeType);

    /**
     * Overridden MarkAsDeleted() method.
     *
     * Mark an element as having been removed from the mesh.
     * Also notify nodes in the element that it has been removed
     * if this node is fully-dimensional
     */
    void MarkAsDeleted();

    /**
     * @return the number of faces owned by this element.
     */
    unsigned GetNumFaces() const;

    /**
     * @param index the global index of a specified face
     *
     * @return a pointer to the face
     */
    MonolayerVertexElement<0, SPACE_DIM>* GetFace(unsigned index) const;

    /**
     * Delete a face with given local index from the element.
     *
     * @param rIndex is the local index of the face to remove
     */
    void DeleteFace(const unsigned& rIndex);

    /**
     * @param index the global index of a specified face
     *
     * @return number of face type
     */
    MonolayerVertexElementType GetFaceType() const;

    /**
     * @param index the local node index of the node
     *
     * @return type of node
     */
    MonolayerVertexElementType GetNodeType(unsigned localIndex) const;

    /**
     * @return whether this is a boundary face.
     */
    bool IsBoundaryFace() const;

    /**
     * @param boundaryFace whether this is a boundary face.
     */
    void SetAsBoundaryFace(bool boundaryFace = true);

    /**
     * @return whether the face with a given index is oriented clockwise.
     *
     * @param index the index of the face
     */
    bool FaceIsOrientatedClockwise(unsigned index) const;
};

#endif /*MONOLAYERVERTEXELEMENT_HPP_*/