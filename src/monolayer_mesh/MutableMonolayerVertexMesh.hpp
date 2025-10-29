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

#ifndef MUTABLEMONOLAYERVERTEXMESH_HPP_
#define MUTABLEMONOLAYERVERTEXMESH_HPP_

// Forward declaration prevents circular include chain
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonolayerVertexMeshWriter;

#include <algorithm>
#include <iostream>
#include <map>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"

#include "MonolayerVertexMesh.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * Struct that saves all information on protorosettes, i.e. the apical node and the previously shared faces
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct MutableMonolayerVertexMeshProtorosette
{
    Node<SPACE_DIM>* pApicalNode;
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pNeighbourOne;
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pNeighbourTwo;
};

/**
 * A mutable vertex-based mesh class, which inherits from VertexMesh and allows for local
 * remeshing. This is implemented through simple operations including node merging, neighbour
 * exchange ("T1 swap"), node/edge merging in the case of intersections ("T3 swap") and
 * removal of small triangular elements ("T2 swap").
 *
 * MutableVertexMesh is used as a member of the VertexBasedCellPopulation class to represent
 * the junctional network of cells that forms the basis of simulations of off-lattice
 * vertex-based models.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MutableMonolayerVertexMesh : public MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestMutableMonolayerVertexMesh;
    friend class TestMutableMonolayerVertexMeshReMesh;
    friend class TestMutableMonolayerVertexMeshRosetteMethods;

protected:
    /** The minimum distance apart that two nodes in the mesh can be without causing element rearrangement. */
    double mCellRearrangementThreshold;

    /**
     * The ratio between the minimum distance apart that two nodes in the mesh can be without causing element
     * rearrangement and their separation after remeshing.
     */
    double mCellRearrangementRatio;

    /** The area threshold at which T2 swaps occur in an apoptotic, triangular cell/element. */
    double mT2Threshold;

    /** The probability that a T1swap is done despite the edges beeing too short (active T1) */
    double mActiveT1SwapProbability;

    /** The probability that, instead of a T1 swap, the relevant nodes merge to form a protorosette */
    double mProtorosetteFormationProbability;

    /** The probability that, in a given timestep, a protorosette node resolves into two rank-3 nodes */
    double mProtorosetteResolutionProbabilityPerTimestep;

    /** The probability that, in a given timestep, a rosette node resolves into two lower-rank nodes */
    double mRosetteResolutionProbabilityPerTimestep;

    /** Whether to check for edges intersections (true) or not (false). */
    bool mCheckForInternalIntersections;

    /** Indices of nodes that have been deleted. These indices can be reused when adding new elements/nodes. */
    std::vector<unsigned> mDeletedNodeIndices;

    /** Indices of elements that have been deleted. These indices can be reused when adding new elements. */
    std::vector<unsigned> mDeletedElementIndices;

    /** Indices of faces that have been deleted. These indices can be reused when adding new faces. */
    std::vector<unsigned> mDeletedFaceIndices;

    /**
     * Distance for T3 swap checking. At each time step we check for each boundary node whether
     * it intersects with any boundary elements (cells) whose centroids lie within this distance
     * to the node. Note that T3 swaps may not be resolved correctly if this distance is chosen
     * too small, while large values for mDistanceForT3SwapChecking may slow down the simulation.
     */
    double mDistanceForT3SwapChecking;

    /**
     * Locations of T1 swaps (the mid point of the moving nodes), stored so they can be accessed and output by the cell population.
     * The locations are stored until they are cleared by ClearLocationsOfT1Swaps().
     */
    std::vector<c_vector<double, SPACE_DIM> > mLocationsOfT1Swaps;

    /**
     * The location of the last T2 swap (the centre of the removed triangle), stored so it can be accessed by the T2SwapCellKiller.
     */
    c_vector<double, SPACE_DIM> mLastT2SwapLocation;

    /**
     * Locations of T3 swaps (the location of the intersection with the edge), stored so they can be accessed and output by the cell population.
     * The locations are stored until they are cleared by ClearLocationsOfT3Swaps().
     */
    std::vector<c_vector<double, SPACE_DIM> > mLocationsOfT3Swaps;

    /** Indices of faces that have been swapped. These faces can be ignored when performing active T1 swaps */
    std::set<unsigned> mSwappedFaceIndices;

    /** Map from basal protorosette nodes to saved information on apical node and previously shared faces */
    std::map<Node<SPACE_DIM>*, MutableMonolayerVertexMeshProtorosette<ELEMENT_DIM, SPACE_DIM> > mMapOfProtorosettes;

    /** Set which contains the pointers to elements with protorosette */
    std::set<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> mSetOfElementsWithProtorosette;

    /**
     * Elements that underwent T1 swaps, stored so they can be accessed by the cell population.
     * The locations are stored until they are cleared by ClearElementsThatUnderwentT1Transitions(),
     * which is called at the beginning of each remeshing
     */
    std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> mElementsThatUnderwentT1Transitions;

    /**
     * Active t1 parameter if active t1 swap is length dependent. Then a Boltzmann-probability is used.
     */
    double mLengthDependentActiveT1BoltzmannParameter;

    /**
     * Flag if active T1 rate is per edge (default is false)
     */
    bool mIsActiveT1ProbPerEdge;

    /**
     * Counter for how many passive T1 transitions happen
     */
    unsigned mPassiveT1TransitionsCounter;

    /**
     * Counter for how many active T1 transitions happen
     */
    unsigned mActiveT1TransitionsCounter;

    /**
     * Divide an element at the given nodes, creating a shared lateral face containing the nodes.
     * The nodes in nodeIndices have to be ordered in anti- or clockwise manner. The new element
     * is created in the direction of the normal vector of the new lateral face, if
     * placeOriginalElementOutside is set to false. Otherwise this  is reversed.
     *
     * @param pElement the element to divide
     * @param nodeIndices array of indicies of nodes to divide the element at in an ordered manner
     * @param placeOriginalElementOutside whether the original element is supposed to be outside the new lateral face
     *
     * @return the index of the new element
     */
    unsigned DivideElement(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                           std::array<unsigned, 4> nodeIndices,
                           bool placeOriginalElementOutside = false);

    /**
     * Divide the lateral face such that a new edge between the new nodes exists as a boundary between two
     * lateral faces. We assume that the new nodes both lie inbetween the two nodes (usally along the connection)
     * such that the orientation of the newly created lateral face is identical. We also add the nodes to the
     * apical and basal faces of the elements which share the lateral face.
     *
     * @param pFace lateral face to divide
     * @param pNewApicalNode apical node which we add to the two created lateral faces
     * @param pNewBasalNode basal node which we add to the two created lateral faces
     */
    void DivideLateralFaceWithNodes(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace,
                                    Node<SPACE_DIM>* pNewApicalNode,
                                    Node<SPACE_DIM>* pNewBasalNode);

    /**
     * Helper method for ReMesh().
     *
     * Check if any neighbouring nodes in an element are closer than the mCellRearrangementThreshold
     * and are not contained in any triangular elements. If any such pair of nodes are found, then
     * call IdentifySwapType(), which in turn implements the appropriate local remeshing operation
     * (a T1 swap, void removal, or node merge).
     *
     * @return whether we need to check for, and implement, any further local remeshing operations
     *                   (true if any swaps are performed).
     */
    virtual bool CheckForSwapsFromShortEdges();

    /**
     * Helper method for ReMesh().
     *
     * Go through all the faces, which have not been swapped passively from CheckForSwapsFromShortEdges
     * and perform active T1 swaps with probability mActiveT1SwapProbability for all edges where T1
     * swaps can be performed. The previously swapped faces (marked in mSwappedFaceIndices) are skipped
     *
     */
    virtual void PerformActiveT1Swaps();

    /**
     * Helper method for ReMesh(), called by CheckForSwapsFromShortEdges() when
     * both apical and basal edges of a lateral face have been found to be shorter than the mCellRearrangementThreshold
     * and do not share any triangular elements.
     *
     * Identify the type of local remeshing operation required (T1 swap, void removal, or node merge).
     * If activeT1Swap is true, we perform the operation only if it is a T1 swap.
     *
     * @param pNodeA one of the basal nodes to perform the swap with
     * @param pNodeB the other basal node to perform the swap
     * @param pFace the lateral face which is swapped
     * @param activeT1Swap whether this is a forced T1 swap, in which case we only do it if it is truly a T1
     */
    virtual void IdentifySwapType(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, bool activeT1Swap = false);

    /**
     * Helper method for ReMesh(), called by IdentifySwapType().
     *
     * Merge two given basal nodes in the mesh and their respective apical nodes, which eliminates
     * a lateral face and update node/element ownership, by replacing
     * the nodes contained in the least number of elements with the other node. The merged
     * node is moved to the centre between the two old node positions.
     * If we check for intersections we use IsFaceConvexAndNonSelfIntersecting() on the result and
     * only perform the merge if this is true. This should only be used in active T1s or we get stuck
     * in the ReMesh loop, as we might not be able to solve the situation
     *
     * @param pBasalNodeA one of the basal nodes to perform the merge with
     * @param pBasalNodeB the other basal node to perform the merge with
     * @param pFace the face to be eliminated
     * @param checkIntersections whether to check for intersections
     *
     * @return whether a node merge has happened (only relevant if checkIntersections is true)
     */
    bool PerformNodeMerge(Node<SPACE_DIM>* pBasalNodeA, Node<SPACE_DIM>* pBasalNodeB,
                          MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, bool checkIntersections = false);

    /**
     * Helper method for ReMesh(), called by IdentifySwapType().
     *
     * A high order junction is a short edge between a protorosette node and a lower order node
     * This issue is resolved by moving the node along the axis from the high order to the
     * low order node such a distance, such that they are the minimum distance apart
     *
     * @param pBasalNodeA one of the basal nodes of the junction
     * @param pBasalNodeB the other basal node of the junction
     * @param pFace the face that connects the nodes.
     */
    void HandleHighOrderJunctions(Node<SPACE_DIM>* pBasalNodeA, Node<SPACE_DIM>* pBasalNodeB, MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace);

    /**
     * Helper method for ReMesh(), called by IdentifySwapType().
     *
     * Perform a T1 swap on two given nodes contained in a given set of elements.
     * This involves replacing the two nodes with two new nodes placed on either
     * side of the previous shared edge, such that the edge formed by the two new nodes
     * is the perpendicular bisector of the previous shared edge, and 'just larger' (by a
     * factor mCellRearrangementRatio) than mThresholdDistance.
     *
     * @param pNodeA one of the nodes to perform the swap with
     * @param pNodeB the other node to perform the swap
     * @param pFace pointer to the lateral face that is swapped
     */
    void PerformT1Swap(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace);

    /**
     * Helper method for ReMesh(), called by CheckForRosettes().
     *
     * Split protorosette which was formed as intermediate step in T1 swaps.
     * Create new nodes and redistribute nodes along line joining
     * centres of cells which were previously separated by a face.
     *
     * @param pProtorosetteNode basal node at centre of protorosette
     */
    void PerformProtorosetteResolution(Node<SPACE_DIM>* pProtorosetteNode);

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialize the mesh.
     *
     * Note that if you are calling this method (from subclasses) you should archive your
     * member variables FIRST. So that this method can call a ReMesh
     * (to convert from TrianglesMeshReader input format into your native format).
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        // NOTE - Subclasses must archive their member variables BEFORE calling this method.
        archive & mCellRearrangementThreshold;
        archive & mCellRearrangementRatio;
        archive & mT2Threshold;
        archive & mActiveT1SwapProbability;
        archive & mProtorosetteFormationProbability;
        archive & mProtorosetteResolutionProbabilityPerTimestep;
        archive & mRosetteResolutionProbabilityPerTimestep;
        archive & mCheckForInternalIntersections;
        archive & mDeletedNodeIndices;
        archive & mDeletedElementIndices;
        archive & mDistanceForT3SwapChecking;
        archive & mLengthDependentActiveT1BoltzmannParameter;
        ///\todo: maybe we should archive the mLocationsOfT1Swaps and mDeletedNodeIndices etc. as well?

        archive& boost::serialization::base_object<MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:
    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangement (defaults to 0.01)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     * @param cellRearrangementRatio ratio between the minimum threshold distance for element
     *                                rearrangement node separation after remeshing (defaults to 1.5)
     * @param protorosetteFormationProbability the probability of a protorosette formation event happening instead of
     *                                a T1 swap
     * @param protorosetteResolutionProbabilityPerTimestep the probability that, in a given timestep, a protorosette
     *                                will resolve (similar to the completion of a T1 swap)
     * @param rosetteResolutionProbabilityPerTimestep the probability that, in a given timestep, a rosette will
     *                                resolve (reduce the number of cells sharing a common vertex by 1)
     */
    MutableMonolayerVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                               std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements,
                               double cellRearrangementThreshold = 0.01,
                               double t2Threshold = 0.001,
                               double cellRearrangementRatio = 1.5,
                               double protorosetteFormationProbability = 0.0,
                               double protorosetteResolutionProbabilityPerTimestep = 0.0,
                               double rosetteResolutionProbabilityPerTimestep = 0.0);

    /**
     * Constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param faces vector of pointer to VertexElements
     * @param vertexElements vector of pointers to VertexElement<3,3>
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangement (defaults to 0.01)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     * @param cellRearrangementRatio ratio between the minimum threshold distance for element
     *                                rearrangement node separation after remeshing (defaults to 1.5)
     * @param protorosetteFormationProbability the probability of a protorosette formation event happening instead of
     *                                a T1 swap
     * @param protorosetteResolutionProbabilityPerTimestep the probability that, in a given timestep, a protorosette
     *                                will resolve (similar to the completion of a T1 swap)
     * @param rosetteResolutionProbabilityPerTimestep the probability that, in a given timestep, a rosette will
     *                                resolve (reduce the number of cells sharing a common vertex by 1)
     */
    MutableMonolayerVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                               std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> faces,
                               std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements,
                               double cellRearrangementThreshold = 0.01,
                               double t2Threshold = 0.001,
                               double cellRearrangementRatio = 1.5,
                               double protorosetteFormationProbability = 0.0,
                               double protorosetteResolutionProbabilityPerTimestep = 0.0,
                               double rosetteResolutionProbabilityPerTimestep = 0.0);

    /**
     * Default constructor for use by serializer.
     */
    MutableMonolayerVertexMesh();

    /**
     * Destructor.
     */
    virtual ~MutableMonolayerVertexMesh();

    /**
     * Set method for mCellRearrangementThreshold.
     *
     * @param cellRearrangementThreshold
     */
    void SetCellRearrangementThreshold(double cellRearrangementThreshold);

    /**
     * Set method for mT2Threshold.
     *
     * @param t2Threshold
     */
    void SetT2Threshold(double t2Threshold);

    /**
     * Set method for mCellRearrangementRatio.
     *
     * @param cellRearrangementRatio
     */
    void SetCellRearrangementRatio(double cellRearrangementRatio);

    /**
     * Set method for mProtoRosetteFormationProbability.
     *
     * @param protorosetteFormationProbability the new value of mProtoRosetteFormationProbability
     */
    void SetProtorosetteFormationProbability(double protorosetteFormationProbability);

    /**
     * Set method for mActiveT1SwapProbability.
     *
     * @param mActiveT1SwapProbability the new value of mActiveT1SwapProbability
     */
    void SetActiveT1SwapProbability(double activeT1SwapProbability, bool isActiveT1ProbabilityPerEdge = false);

    /**
     * Set method for mProtoRosetteResolutionProbabilityPerTimestep.
     *
     * @param protorosetteResolutionProbabilityPerTimestep the new value of mProtoRosetteResolutionProbabilityPerTimestep
     */
    void SetProtorosetteResolutionProbabilityPerTimestep(double protorosetteResolutionProbabilityPerTimestep);

    /**
     * Set method for mRosetteResolutionProbabilityPerTimestep.
     *
     * @param rosetteResolutionProbabilityPerTimestep the new value of mRosetteResolutionProbabilityPerTimestep
     */
    void SetRosetteResolutionProbabilityPerTimestep(double rosetteResolutionProbabilityPerTimestep);

    /**
     * Move the node with a particular index to a new point in space.
     *
     * @param nodeIndex the index of the node to be moved
     * @param point the new target location of the node
     */
    virtual void SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point);

    /**
     * Set method for mCheckForInternalIntersections.
     *
     * @param checkForInternalIntersections
     */
    void SetCheckForInternalIntersections(bool checkForInternalIntersections);

    /**
     * @return mCellRearrangementThreshold
     */
    double GetCellRearrangementThreshold() const;

    /**
     * @return mT2Threshold
     */
    double GetT2Threshold() const;

    /**
     * @return mCellRearrangementRatio
     */
    double GetCellRearrangementRatio() const;

    /**
     * Get method for mProtoRosetteFormationProbability.
     *
     * @return mProtoRosetteFormationProbability
     */
    double GetProtorosetteFormationProbability() const;

    /**
     * Get method for mActiveT1SwapProbability.
     *
     * @return mActiveT1SwapProbability
     */
    double GetActiveT1SwapProbability() const;

    /**
     * Get method for mProtoRosetteResolutionProbabilityPerTimestep.
     *
     * @return mProtoRosetteResolutionProbabilityPerTimestep
     */
    double GetProtorosetteResolutionProbabilityPerTimestep() const;

    /**
     * Get method for mRosetteResolutionProbabilityPerTimestep.
     *
     * @return mRosetteResolutionProbabilityPerTimestep
     */
    double GetRosetteResolutionProbabilityPerTimestep() const;

    /**
     * Set distance for T3 swap checking. At each time step we check for each boundary node whether
     * it intersects with any boundary elements (cells) whose centroids lie within this distance
     * to the node. Note that T3 swaps may not be resolved correctly if this distance is chosen
     * too small, while large values for mDistanceForT3SwapChecking may slow down the simulation.
     *
     * @param distanceForT3SwapChecking
     */
    void SetDistanceForT3SwapChecking(double distanceForT3SwapChecking);

    /**
     * Get Distance for T3 swap checking.
     *
     * @return mDistanceForT3SwapChecking
     */
    double GetDistanceForT3SwapChecking() const;

    /**
     * @return the number of Nodes in the mesh.
     */
    unsigned GetNumNodes() const;

    /**
     * @return the number of VertexElements in the mesh.
     */
    unsigned GetNumElements() const;

    /**
     * @return mCheckForInternalIntersections, either to check for edges intersections or not.
     */
    bool GetCheckForInternalIntersections() const;

    /**
     * @return the locations of the T1 swaps
     */
    std::vector<c_vector<double, SPACE_DIM> > GetLocationsOfT1Swaps();

    /**
     * @return the location of the last T2 swap
     */
    c_vector<double, SPACE_DIM> GetLastT2SwapLocation();

    /**
     * @return the locations of the T3 swaps
     */
    std::vector<c_vector<double, SPACE_DIM> > GetLocationsOfT3Swaps();

    /**
     * @return the elements which underwent T1 swaps
     */
    std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> GetElementsThatUnderwentT1Transitions();

    /**
     * Helper method to clear the stored indices of swapped Face indices.
     */
    void ClearSwappedFaceIndices();

    /**
     * Helper method to clear the stored T1 swaps
     */
    void ClearLocationsOfT1Swaps();

    /**
     * Helper method to clear the stored T3 swaps
     */
    void ClearLocationsOfT3Swaps();

    /**
     * Helper method to clear the Protorosette map
     */
    void ClearMapOfProtorosettes();

    /**
     * Helper method to clear the set of elements which contain a protorosette
     */
    void ClearSetOFElementsWithProtorosette();

    /**
     * Helper method to clear the set of elements which underwent a T1 transtion
     */
    void ClearElementsThatUnderwentT1Transitions();

    /**
     * Add a node to the mesh.
     *
     * Note: After calling this one or more times, you must then call ReMesh.
     *
     * @param pNewNode pointer to the new node
     * @return the global index of the new node in the mesh.
     */
    unsigned AddNode(Node<SPACE_DIM>* pNewNode);

    /**
     * Add a face to the mesh.
     *
     * Note: After calling this one or more times, you must then call ReMesh.
     *
     * @param pNewNode pointer to the new face
     * @return the global index of the new face in the mesh.
     */
    unsigned AddFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pNewFace);

    /**
     * Check whether the face is convexish, i.e. the angles with respect to the center are
     * consistent, and non-intersecting.
     * For this we check whether the triangulation has normals pointing in the same direction
     * and whether the scalarproduct of the normals of the triangles spanned by two edges have
     * the same sign
     *
     * @param pFace pointer to the face
     * @return whether the face is convex and non-selfintersecting
     */
    bool IsFaceConvexAndNonSelfIntersecting(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace);

    /**
     * Mark an element as deleted. Note that it DOES NOT deal with the associated
     * nodes and therefore should only be called immediately prior to a ReMesh()
     * being called.
     *
     * @param index  the global index of a specified vertex element
     */
    void DeleteElementPriorToReMesh(unsigned index);

    /**
     * Mark a given node as deleted. Note that this method DOES NOT deal with the
     * associated elements and therefore should only be called immediately prior
     * to a ReMesh() being called.
     *
     * @param index The index of the node to delete
     */
    void DeleteNodePriorToReMesh(unsigned index);

    /**
     * Mark a given face as deleted. Note that this method DOES NOT deal with the
     * associated elements and therefore should only be called immediately prior
     * to a ReMesh() being called.
     *
     * @param index The index of the node to delete
     */
    void DeleteFacePriorToReMesh(unsigned index);

    /**
     * Divide an element along its short axis.
     *
     * @param pElement the element to divide
     * @param placeOriginalElementBelow whether to place the original element below (in the y direction) the new element (defaults to false)
     *
     * @return the index of the new element
     */
    unsigned DivideElementAlongShortAxis(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                         bool placeOriginalElementBelow = false);

    /**
     * Divide an element along a specified axis.
     *
     * If the new nodes (intersections of axis with element) are within
     * mCellRearrangementThreshold of existing nodes then they are
     * moved 2*mCellRearrangementThreshold away.
     *
     * @param pElement the element to divide
     * @param axisOfDivision axis to divide the element by
     * @param placeOriginalElementBelow whether to place the original element below (in the y direction) the new element (defaults to false)
     *
     * @return the index of the new element
     */
    unsigned DivideElementAlongGivenAxis(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement,
                                         c_vector<double, SPACE_DIM> axisOfDivision,
                                         bool placeOriginalElementBelow = false, unsigned iterative_call_counter = 0);

    /**
     * Add an element to the mesh.
     *
     * @param pNewElement the new element
     *
     * @return the index of the new element in the mesh
     */
    unsigned AddElement(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pNewElement);

    /**
     * Helper method for ReMesh().
     *
     * Check for any triangular element whose area is smaller than mT2Threshold
     * and call PerformT2Swap() on any such element.
     *
     * @param rElementMap a VertexElementMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created
     *                   with the correct size, GetNumElements()
     *
     * @return whether we need to check for, and implement, any further local remeshing operations
     */
    bool CheckForT2Swaps(VertexElementMap& rElementMap);

    /**
     * Helper method for ReMesh(), called by IdentifySwapType().
     *
     * Remove a triangular void bounded by three given nodes, in which one of the edges is
     * less than mCellRearrangementThreshold, through calls to PerformNodeMerge().
     *
     * @param pNodeA one of the nodes on the short edge
     * @param pNodeB the other node on the short edge
     * @param pNodeC the other node in the triangular void
     * @param pFaceAB lateral face connected to nodes A and B
     * @param pFaceBC lateral face connected to nodes B and C
     * @param pFaceAC lateral face connected to nodes A and C
     */
    void PerformVoidRemoval(Node<SPACE_DIM>* pBasalNodeA, Node<SPACE_DIM>* pBasalNodeB, Node<SPACE_DIM>* pBasalNodeC,
                            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFaceAB,
                            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFaceBC,
                            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFaceAC);

    /**
     * Helper method for ReMesh().
     *
     * Check whether the mesh contains rosettes or protorosettes, and implement resolution events
     * if necessary.
     */
    void CheckForRosettes();

    /**
     * Add an entry to the map from the basal node to the corresponding apical node
     * and pointers to elements which previously shared an edge
     *
     * @param pBasalNode The basal node of the protorosette
     * @param pApicalNode The basal node of the protorosette
     * @param pNeighbourElementA The first neighbouring element of the deleted face
     * @param pNeighbourElementA The second neighbouring element of the deleted face
     */
    void AddToMapOfProtorosettes(Node<SPACE_DIM>* pBasalNode, Node<SPACE_DIM>* pApicalNode,
                                 MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pNeighbourElementA,
                                 MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pNeighbourElementB);

    /**
     * Add the pointer of the element to the list of elements which have a protorosette
     *
     * @param pElement pointer of element
     */
    void AddToSetOfElementsWithProtorosette(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement);

    /**
     * Retrieve Protorosette object which contains a pointer to the apical node of the protorosette
     * and pointers to the elements which contained the deleted faces when the rosette was formed
     *
     * @param pBasalNode The basal node of the protorosette
     *
     * @return Protorosette object correspinding to basal node
     */
    MutableMonolayerVertexMeshProtorosette<ELEMENT_DIM, SPACE_DIM> GetMappedProtorosette(Node<SPACE_DIM>* pBasalNode);

    /**
     * Check if element does have a node which is part of a Protorosette. For this we check the map
     * of basal nodes to protorosettes, which we maintain for resolution of protorosettes anyways.
     *
     * @param pElement pointer to the element
     *
     * @return boolean whether this element contains a protorosette
     */
    bool DoesElementHaveProtorosette(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement);

    /**
     * Remove map to Protorosette object which contains a pointer to the apical node of the protorosette
     * and pointers to the elements which contained the deleted faces when the rosette was formed
     *
     * @param pBasalNode The basal node of the protorosette to remove from the map
     */
    void RemoveMappedProtorosette(Node<SPACE_DIM>* pBasalNode);

    /**
     * Remove the pointer of the element from the list of elements which have a protorosette
     *
     * @param pElement pointer of element
     */
    void RemoveFromSetOfElementsWithProtorosette(MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement);

    /**
     * Delete mNodes and mElements.
     */
    void Clear();

    /**
     * Add a node on the edge between two nodes.
     *
     * @param pNodeA a pointer to one node
     * @param pNodeB a pointer to the other nodes
     */
    void DivideEdge(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);

    /**
     * Helper method for ReMesh(). Removes the deleted nodes and elements from the mesh and updates the
     * rElementMap accordingly.
     *
     * @param rElementMap a VertexElementMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created
     *                   with the correct size, GetNumElements()
     */
    void RemoveDeletedNodesAndElements(VertexElementMap& rElementMap);

    /**
     * Helper method for ReMesh(). Removes the deleted nodes from the mesh and relabels the node indices.
     */
    void RemoveDeletedNodes();

    /**
     * Helper method for ReMesh(). Removes the deleted faces from the mesh and relabels the face indices.
     */
    void RemoveDeletedFaces();

    /**
     * Update the state of the mesh by implementing any local remeshing operations (node merging,
     * or T1, T2 or T3 swaps) that are required, and store any changes in element indices using
     * the given VertexElementMap.
     *
     * This method calls several other methods, in particular CheckForT2Swaps(), CheckForSwapsFromShortEdges()
     * and CheckForIntersections().
     *
     * @param rElementMap a VertexElementMap which associates the indices of VertexElements in the old mesh
     *                   with indices of VertexElements in the new mesh.  This should be created
     *                   with the correct size, GetNumElements()
     */
    virtual void ReMesh(VertexElementMap& rElementMap);

    /**
     * Alternative version of ReMesh which takes no parameters and does not require a VertexElementMap.
     * Note: inherited classes should overload ReMesh(VertexElementMap&).
     *
     * \todo This method seems to be redundant; remove it? (#2401)
     */
    void ReMesh();

    /**
     * Calculate the probability for a T1 transformation.
     *
     * If the Boltzmann method is used then the averaged mid-plane edge length is used and
     * a Boltzmann factor exp(-(l-l*)/L0)
     * with the  mCellRearrangementThreshold l* and mLengthDependentActiveT1BoltzmannParameter L0
     * Otherwise mActiveT1SwapProbability is used (and normalized per #lateral faces if not mIsActiveT1ProbPerEdge).
     */
    double CalculateActiveT1SwapProbability(double apical_length, double basal_length);

    /**
     * Set if the active T1 swap probability is length dependent.
     */
    void SetActiveT1BoltzmannParameter(double boltzmann_parameter);

    /**
     * Get and reset passive T1 transition count
     */
    unsigned GetAndResetPassiveT1Count();

    /**
     * Get and reset active T1 transition count
     */
    unsigned GetAndResetActiveT1Count();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MutableMonolayerVertexMesh)

#endif /*MUTABLEMONOLAYERVERTEXMESH_HPP_*/