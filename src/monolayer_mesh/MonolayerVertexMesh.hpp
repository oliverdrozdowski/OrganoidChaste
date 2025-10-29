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
#ifndef MONOLAYERVERTEXMESH_HPP_
#define MONOLAYERVERTEXMESH_HPP_

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

#include "AbstractMesh.hpp"
#include "ArchiveLocationInfo.hpp"
#include "MonolayerVertexElement.hpp"
#include "MonolayerVertexMeshReader.hpp"
#include "MonolayerVertexMeshWriter.hpp"
#include "TetrahedralMesh.hpp"
#include "VertexElementMap.hpp"

/**
 * A vertex-based mesh class, in which elements may contain different numbers of nodes.
 * This is facilitated by the VertexElement class.
 *
 * This class has only one application in the vertex model (no voronoi implementation)
 *
 * MonolayerVertexMesh serves as a parent class for MutableMonolayerVertexMesh, which is used as a
 * member of the MonolayerVertexBasedCellPopulation class to represent the junctional network of
 * cells that forms the basis of simulations of finite-thickness off-lattice vertex-based models.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonolayerVertexMesh : public AbstractMesh<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestMonolayerVertexMesh;

protected:
    /** Vector of pointers to VertexElements. */
    std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> mElements;

    /** Vector of pointers to VertexElements. */
    std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> mFaces;

    /** Map of element indices to set of face indices*/
    std::map<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*, std::vector<unsigned> > mElementsFacesMap;

    /**
     * Solve node mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the node
     * @return local index
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Solve element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the element
     * @return local index
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Solve boundary element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the boundary element
     * @return local index
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;

    /**
     * Test whether a given point lies inside a given element.
     *
     * We use a winding number test, which counts the number of times the
     * polygon associated with the element winds around the given point.
     * The point is outside only when this "winding number" vanishes;
     * otherwise, the point is inside.
     *
     * One must decide whether a point on the polygon's boundary is inside
     * or outside: we adopt the standard convention that a point on a left
     * or bottom edge is inside, and a point on a right or top edge is outside.
     * This way, if two distinct polygons share a common boundary segment,
     * then a point on that segment will be in one polygon or the other, but
     * not both at the same time.
     *
     * @param rTestPoint the point to test
     * @param elementIndex global index of the element in the mesh
     *
     * @return if the point is included in the element.
     */
    bool ElementIncludesPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex);

    /**
     * Get the local index of a given element which is the start vertex of the edge
     * of the element that the overlapping point rTestPoint is closest to.
     *
     * @param rTestPoint the point to test
     * @param elementIndex global index of the element in the mesh
     *
     * @return the local index
     */
    unsigned GetLocalIndexForElementEdgeClosestToPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex);

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the VertexMesh and its member variables. Note that this will
     * write out a VertexMeshWriter file to wherever ArchiveLocationInfo has specified.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void save(Archive& archive, const unsigned int version) const
    {
        archive& boost::serialization::base_object<AbstractMesh<ELEMENT_DIM, SPACE_DIM> >(*this);

        // Create a mesh writer pointing to the correct file and directory
        MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
                                                                      ArchiveLocationInfo::GetMeshFilename(),
                                                                      false);
        mesh_writer.WriteFilesUsingMesh(*(const_cast<MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>*>(this)));
    }

    /**
     * Load a mesh by using VertexMeshReader and the location in ArchiveLocationInfo.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void load(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractMesh<ELEMENT_DIM, SPACE_DIM> >(*this);

        MonolayerVertexMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());
        this->ConstructFromMeshReader(mesh_reader);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
    /** Forward declaration of element iterator. */
    class MonolayerVertexElementIterator;

    /** Forward declaration of face iterator. */
    class MonolayerVertexElementFaceIterator;

    /**
     * @return an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline MonolayerVertexElementIterator GetElementIteratorBegin(bool skipDeletedElements = true);

    /**
     * @return an iterator to one past the last element in the mesh.
     */
    inline MonolayerVertexElementIterator GetElementIteratorEnd();

    /**
     * @return an iterator to the first face in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline MonolayerVertexElementFaceIterator GetFaceIteratorBegin(bool skipDeletedElements = true);

    /**
     * @return an iterator to one past the last face in the mesh.
     */
    inline MonolayerVertexElementFaceIterator GetFaceIteratorEnd();

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     */
    MonolayerVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                        std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements);

    /**
     * Constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param faces vector of pointer to VertexElements
     * @param vertexElements vector of pointers to VertexElement<3,3>s
     */
    MonolayerVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                        std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> faces,
                        std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements);

    /**
     * Alternative 3D 'Voronoi' constructor. Creates a Voronoi tessellation of a given tetrahedral mesh,
     * which must be Delaunay (see TetrahedralMesh::CheckIsVoronoi).
     *
     * \todo Merge with 2D Voronoi constructor? (see #1075)
     *
     * @param rMesh a tetrahedral mesh
     */
    // VertexMesh(TetrahedralMesh<3, 3>& rMesh);

    /**
     * Default constructor for use by serializer.
     */
    MonolayerVertexMesh();

    /**
     * Destructor.
     */
    virtual ~MonolayerVertexMesh();

    /**
     * @return the number of Nodes in the mesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * @return the number of VertexElements in the mesh.
     */
    virtual unsigned GetNumElements() const;

    /**
     * @return the number of VertexElements in the mesh, including those marked as deleted.
     */
    unsigned GetNumAllElements() const;

    /**
     * @return the number of Faces in the mesh.
     */
    virtual unsigned GetNumFaces() const;

    /**
     * @return the number of Faces of certain type in the mesh.
     */
    virtual unsigned GetNumFacesByType(MonolayerVertexElementType face_type) const;

    /**
     * @param index  the global index of a specified vertex element.
     *
     * @return a pointer to the vertex element
     */
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /**
     * @param index  the global index of a specified face.
     *
     * @return a pointer to the face
     */
    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* GetFace(unsigned index) const;

    /**
     * Get the first face in the element of a face type
     *
     * @param index  the index of the element for which we want the face.
     * @param faceType the face type of the face we want to find
     *
     * @return a pointer to the face
     */
    MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* GetFaceOfType(unsigned index, MonolayerVertexElementType faceType) const;

    /**
     * Compute the centroid of an element.
     *
     * A formula for the centroid of a plane polygon may be found e.g. in the following reference:
     *
     * Mechanics of Materials
     * James M. Gere (Author), Barry J. Goodno.
     * Cengage Learning; 8th edition (January 1, 2012)
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (centroid_x, centroid_y).
     */
    virtual c_vector<double, SPACE_DIM> GetCentroidOfElement(unsigned index);

    /**
     * Update the map from elements to faces
     *
     */
    void UpdateElementsFacesMap();

    /**
     * Update the map from elements to faces for only one element
     *
     */
    void UpdateElementsFacesMapOfElement(unsigned index);

    /**
     * Get the map from elements to faces
     *
     */
    std::map<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*, std::vector<unsigned> >* GetElementsFacesMap();

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader);

    /**
     * Delete mNodes, mFaces and mElements.
     */
    virtual void Clear();

    /**
     * Get the "rosette rank" of an element.
     *
     * This is defined as the maximum number of elements shared by any node in the specified element.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the rosette rank of the element
     */
    unsigned GetRosetteRankOfElement(unsigned index);

    /**
     * Overridden GetVectorFromAtoB() method. Returns a vector between two points in space.
     *
     * If the mesh is being used to represent a Voronoi tessellation, and mpDelaunayMesh
     * is not NULL, then use that to compute GetVectorFromAtoB.
     *
     * @param rLocationA a c_vector of coordinates
     * @param rLocationB a c_vector of coordinates
     *
     * @return c_vector from location A to location B.
     */
    virtual c_vector<double, SPACE_DIM> GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocationA,
                                                          const c_vector<double, SPACE_DIM>& rLocationB);

    /**
     * Get the volume (or area in 2D, or length in 1D) of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the volume of the element
     */
    virtual double GetVolumeOfElement(unsigned index);

    /**
     * Compute the surface area (or perimeter in 2D) of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the surfacearea of the element
     */
    virtual double GetSurfaceAreaOfElement(unsigned index);

    /**
     * Compute the mid-plane area, which is defined as the area spanned by the mid-points of the
     * apico-basal 2D edges. In 2D this is just the surface area of an element.
     *
     * @param index the global index of a specified vertex element
     *
     * @return the mid plane area of the element
     */
    virtual double GetMidPlaneAreaOfElement(unsigned index);

    /**
     * Compute the unit normal vector to a given face in 3D. This is achieved by calculating scaled normal,
     * which is the effective sum of signed areas of triangle forming the face.
     * Note: this may return the outward or inward normal, depending
     * on the face chirality.
     *
     * @param pFace a face in the mesh
     * @param rNormal vector in which to return the unit normal
     *
     * @return the area
     */
    double CalculateUnitNormalToFaceWithArea(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, c_vector<double, SPACE_DIM>& rNormal);

    /**
     * Get the area of a given face in 3D.  Uses CalculateUnitNormalToFaceWithArea
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param pFace a face in the mesh
     *
     * @return the area
     */
    virtual double CalculateAreaOfFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace);

    /**
     * Get area and normals of triangles of a given face's triangulation with passive center points in 3D.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param pFace a face in the mesh
     * @param pTriangleNormals pointer to a vector in which we save the normals of the triangles
     *
     * @return the areas of the different triangles
     */
    virtual std::vector<double> CalculateAreasAndNormalsOfTriangulationsOfFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, std::vector<c_vector<double, SPACE_DIM> >* pTriangleNormals);

    /**
     * Given a node, find a set containing the indices of its neighbouring nodes.
     *
     * @param nodeIndex global index of the node
     * @return its neighbouring node indices
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);

    /**
     * Given a node and one of its containing elements, find a set containing
     * the indices of those neighbouring node(s) that are NOT also in the element.
     *
     * Note that we allow for more than one such index, since there is no reason
     * a priori to assume that each node is contained by exactly three elements.
     *
     * @param nodeIndex global index of the node
     * @param elemIndex global index of the element
     *
     * @return its neighbouring nodes that are not in the element
     */
    std::set<unsigned> GetNeighbouringNodeNotAlsoInElement(unsigned nodeIndex, unsigned elemIndex);

    /**
     * Given an element, find a set containing the indices of its neighbouring elements.
     *
     * @param elementIndex global index of the element
     * @return its neighbouring element indices
     */
    std::set<unsigned> GetNeighbouringElementIndices(unsigned elementIndex);

    /**
     * Calculate the thickness (height) of an element, as the distance between the average
     * points (passive centers) of apical an basal sides.
     *
     * @param elementIndex global index of the element
     * @return the thickness of the cell
     */
    double GetThicknessOfElement(unsigned elementIndex);

    /**
     * Compute the area gradient of 2D face at one of its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), but cannot be generalized to
     * non-euclidean metrics.
     *
     * @param pElement  pointer to a specified face
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the area of the element, evaluated at this node.
     */
    c_vector<double, SPACE_DIM> GetAreaGradientOfFaceAtNode(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, unsigned localIndex);

    /**
     * Compute the volume gradient of 3D element at one of its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), but cannot be generalized to
     * non-euclidean metrics.
     *
     * @param elementIndex  local index of element
     * @param localNodeIndex  local index of a node in this element
     *
     * @return the gradient of the volume of the element, evaluated at this node.
     */
    c_vector<double, SPACE_DIM> GetVolumeGradientAtNode(unsigned elementIndex, unsigned localNodeIndex);

    /**
     * Return a pointer to the vertex mesh
     *
     * This method may be overridden in daughter classes for non-Euclidean metrics.
     * This can then be used when writing to VTK.
     *
     * @return a pointer to the vertex mesh
     */
    virtual MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>* GetMeshForVtk();

    /**
     * Helper method to determine the centroid of a Monolayer Vertex Element (volume, face or edge)
     * @return vector to position of centroid
     *
     * @param pFace face where to calculate the centroid
     */
    c_vector<double, SPACE_DIM> GetPassiveCenterOfFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace);

    /**
     * Helper method to determine the passive center of the (first) face in an element
     * with the given face type
     * @return vector to position of centroid
     *
     * @param elementIndex index of element
     * @param faceType face type for which we want to determine the first face
     */
    c_vector<double, SPACE_DIM> GetPassiveCenterOfFaceTypeInElement(unsigned elementIndex, MonolayerVertexElementType faceType);

    /**
     * Compute the second moments and product moment of area for a given 2D face
     * projected onto the orthogonal space to the normal about its passive center.
     * The moments are given w.r.t. to the supplied ONB of the orthogonal space
     * These are:
     *
     * I_xx, the second moment of area about an axis through the centroid of the
     * element parallel to the supplied first axis;
     *
     * I_yy, the second moment of area about an axis through the centroid of the
     * element parallel to the supplied first axis;
     *
     * and I_xy, product moment of area through the centroid of the element.
     *
     * Formulae for these quantities may be found e.g. in the following reference:
     *
     * Mechanics of Materials
     * James M. Gere (Author), Barry J. Goodno.
     * Cengage Learning; 8th edition (January 1, 2012)
     *
     * This method is used within GetShortAxisOfElement() to compute the direction
     * of the shortest principal axis passing through the centroid, or 'short axis',
     * of the element.
     *
     * Note that by definition, the second moments of area must be non-negative,
     * while the product moment of area may not be.
     *
     * @param pFace pointer to the face
     * @param normalVector mean normal vector
     * @param otthogonalVector1 x-vector of ONB orthogonal to mean normal
     * @param orthogonalVector2 y-vector of ONB orthogonal to mean normal
     *
     * @return (Ixx,Iyy,Ixy).
     */
    virtual c_vector<double, 3> CalculateMomentsOfFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace,
                                                       c_vector<double, SPACE_DIM> normalVector,
                                                       c_vector<double, SPACE_DIM> orthogonalVector1,
                                                       c_vector<double, SPACE_DIM> orthogonalVector2);

    /**
     * Compute the direction of the mean longest principal axes passing through the passive centers,
     * or 'long axes', of the apical and basal faces of a given element.
     * This is the eigenvector associated with the eigenvalue of smallest magnitude of the inertia
     * matrix for the faces with respect to the basis orthogonal to the mean normal.
     *
     * J = (  I_xx  -I_xy )
     *     ( -I_xy   I_yy )
     *
     * whose entries are computed by calling the method CalculateMomentsOfFace().
     * We consider the moments of the faces projected onto the space orthogonal to the mean normal.
     *
     * Note that if the nodes owned by the element are supplied in clockwise rather than
     * anticlockwise manner, or if this arises when any periodicity is enforced, then the
     * sign of each moment may be incorrect change. This means that we need to consider the eigenvalue
     * of largest magnitude rather than largest value when computing the short axis of the
     * element.
     *
     * If the element has a regular polygon face then the eigenvalues of the inertia tensor are
     * equal: in this case we consider a random unit vector for this face.
     *
     * This method is only implemented in 3D at present.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return a unit vector giving the direction of the long in-plane axis
     */
    c_vector<double, SPACE_DIM> GetMeanApicalBasalLongAxisOfElement(unsigned index);

    /**
     * Compute the direction of the longest principal axes passing through the passive center,
     * or 'long axes', of the basal face of a given element.
     * This is the eigenvector associated with the eigenvalue of smallest magnitude of the inertia
     * matrix for the faces with respect to the basis orthogonal to the mean normal.
     *
     * J = (  I_xx  -I_xy )
     *     ( -I_xy   I_yy )
     *
     * whose entries are computed by calling the method CalculateMomentsOfFace().
     * We consider the moments of the face projected onto the space orthogonal to the mean normal.
     *
     * Note that if the nodes owned by the element are supplied in clockwise rather than
     * anticlockwise manner, or if this arises when any periodicity is enforced, then the
     * sign of each moment may be incorrect change. This means that we need to consider the eigenvalue
     * of largest magnitude rather than largest value when computing the short axis of the
     * element.
     *
     * If the element has a regular polygon face then the eigenvalues of the inertia tensor are
     * equal: in this case we consider a random unit vector for this face.
     *
     * This method is only implemented in 3D at present.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return a unit vector giving the direction of the long in-plane axis
     */
    c_vector<double, SPACE_DIM> GetBasalLongAxisOfElement(unsigned index);

    /**
     * Calculate volume of lumen
     */

    double GetLumenVolume(bool symmetric_content = true);

    /**
     * calculate gradient of lumen volume in respect to one vertex in one face
     */
    c_vector<double, ELEMENT_DIM>
    CalculateLumenVolGradient(unsigned node_index, bool symmetric_content = true);

    /**
     * A smart iterator over the elements in the mesh.
     *
     * \todo This is the same as in AbstractTetrahedralMesh and PottsMesh - merge? (#1379)
     */
    class MonolayerVertexElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         * @return reference
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         * @return true if not equal
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator& rOther);

        /**
         * Prefix increment operator.
         * @return reference to incremented object
         */
        inline MonolayerVertexElementIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * VertexMesh::GetElementIteratorBegin and VertexMesh::GetElementIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        MonolayerVertexElementIterator(MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                       typename std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*>::iterator elementIter,
                                       bool skipDeletedElements = true);

    private:
        /** The mesh we're iterating over. */
        MonolayerVertexMesh& mrMesh;

        /** The actual element iterator. */
        typename std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*>::iterator mElementIter;

        /** Whether to skip deleted elements. */
        bool mSkipDeletedElements;

        /**
         * Helper method to say when we're at the end.
         * @return true if at end
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this element.
         * @return true if allowed
         */
        inline bool IsAllowedElement();
    };

    /**
     * A smart iterator over the faces in the mesh.
     */
    class MonolayerVertexElementFaceIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         * @return reference
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         * @return true if not equal
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator& rOther);

        /**
         * Prefix increment operator.
         * @return reference to incremented object
         */
        inline MonolayerVertexElementFaceIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * VertexMesh::GetElementIteratorBegin and VertexMesh::GetElementIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        MonolayerVertexElementFaceIterator(MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                           typename std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*>::iterator faceIter,
                                           bool skipDeletedElements = true);

    private:
        /** The mesh we're iterating over. */
        MonolayerVertexMesh& mrMesh;

        /** The actual element iterator. */
        typename std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*>::iterator mFaceIter;

        /** Whether to skip deleted elements. */
        bool mSkipDeletedElements;

        /**
         * Helper method to say when we're at the end.
         * @return true if at end
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this element.
         * @return true if allowed
         */
        inline bool IsAllowedElement();
    };
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MonolayerVertexMesh)

//////////////////////////////////////////////////////////////////////////////
// VertexElementIterator class implementation - most methods are inlined    //
//////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorBegin(
    bool skipDeletedElements)
{
    return MonolayerVertexElementIterator(*this, mElements.begin(), skipDeletedElements);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorEnd()
{
    return MonolayerVertexElementIterator(*this, mElements.end());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>& MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator::operator!=(const typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator& MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    } while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator::MonolayerVertexElementIterator(
    MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
    typename std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*>::iterator elementIter,
    bool skipDeletedElements)
        : mrMesh(rMesh),
          mElementIter(elementIter),
          mSkipDeletedElements(skipDeletedElements)
{
    if (mrMesh.mElements.empty())
    {
        // Cope with empty meshes
        mElementIter = mrMesh.mElements.end();
    }
    else
    {
        // Make sure we start at an allowed element
        if (mElementIter == mrMesh.mElements.begin() && !IsAllowedElement())
        {
            ++(*this);
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}

//////////////////////////////////////////////////////////////////////////////////
// VertexElementFaceIterator class implementation - most methods are inlined    //
//////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetFaceIteratorBegin(
    bool skipDeletedElements)
{
    return MonolayerVertexElementFaceIterator(*this, mFaces.begin(), skipDeletedElements);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetFaceIteratorEnd()
{
    return MonolayerVertexElementFaceIterator(*this, mFaces.end());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>& MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator::operator*()
{
    assert(!IsAtEnd());
    return **mFaceIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator::operator->()
{
    assert(!IsAtEnd());
    return *mFaceIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator::operator!=(const typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator& rOther)
{
    return mFaceIter != rOther.mFaceIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator& MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator::operator++()
{
    do
    {
        ++mFaceIter;
    } while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator::MonolayerVertexElementFaceIterator(
    MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
    typename std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*>::iterator faceIter,
    bool skipDeletedElements)
        : mrMesh(rMesh),
          mFaceIter(faceIter),
          mSkipDeletedElements(skipDeletedElements)
{
    if (mrMesh.mFaces.empty())
    {
        // Cope with empty meshes
        mFaceIter = mrMesh.mFaces.end();
    }
    else
    {
        // Make sure we start at an allowed element
        if (mFaceIter == mrMesh.mFaces.begin() && !IsAllowedElement())
        {
            ++(*this);
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator::IsAtEnd()
{
    return mFaceIter == mrMesh.mFaces.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}

#endif /*MONOLAYERVERTEXMESH_HPP_*/