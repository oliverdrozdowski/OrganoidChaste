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

#include "MonolayerVertexMesh.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                                                 std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements)
{

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
        mElements.push_back(p_temp_vertex_element);
    }

    // In 3D, populate mFaces
    if (SPACE_DIM == 3)
    {
        // Use a std::set to keep track of which faces have been added to mFaces
        std::set<unsigned> faces_counted;

        // Loop over mElements
        for (unsigned elem_index = 0; elem_index < mElements.size(); elem_index++)
        {
            // Loop over faces of this element
            for (unsigned face_index = 0; face_index < mElements[elem_index]->GetNumFaces(); face_index++)
            {
                MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = mElements[elem_index]->GetFace(face_index);
                unsigned global_index = p_face->GetIndex();

                // If this face is not already contained in mFaces, add it, and update faces_counted
                if (faces_counted.find(global_index) == faces_counted.end())
                {
                    mFaces.push_back(p_face);
                    faces_counted.insert(global_index);
                }
            }
        }
    }

    // Register elements with nodes
    for (unsigned index = 0; index < mElements.size(); index++)
    {
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = mElements[index];

        unsigned element_index = p_element->GetIndex();
        unsigned num_nodes_in_element = p_element->GetNumNodes();

        for (unsigned node_index = 0; node_index < num_nodes_in_element; node_index++)
        {
            p_element->GetNode(node_index)->AddElement(element_index);
        }
    }

    this->mMeshChangesDuringSimulation = false;
    UpdateElementsFacesMap();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
                                                                 std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> faces,
                                                                 std::vector<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements)
{
    // Reset member variables and clear mNodes, mFaces and mElements
    Clear();

    // Populate mNodes mFaces and mElements
    for (unsigned node_index = 0; node_index < nodes.size(); node_index++)
    {
        Node<SPACE_DIM>* p_temp_node = nodes[node_index];
        this->mNodes.push_back(p_temp_node);
    }

    for (unsigned face_index = 0; face_index < faces.size(); face_index++)
    {
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_temp_face = faces[face_index];
        mFaces.push_back(p_temp_face);
    }

    for (unsigned elem_index = 0; elem_index < vertexElements.size(); elem_index++)
    {
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = vertexElements[elem_index];
        mElements.push_back(p_temp_vertex_element);
    }

    // Register elements with nodes
    for (unsigned index = 0; index < mElements.size(); index++)
    {
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_temp_vertex_element = mElements[index];
        for (unsigned node_index = 0; node_index < p_temp_vertex_element->GetNumNodes(); node_index++)
        {
            Node<SPACE_DIM>* p_temp_node = p_temp_vertex_element->GetNode(node_index);
            p_temp_node->AddElement(p_temp_vertex_element->GetIndex());
        }
    }

    this->mMeshChangesDuringSimulation = false;
    UpdateElementsFacesMap();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexMesh()
{
    this->mMeshChangesDuringSimulation = false;
    Clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::~MonolayerVertexMesh()
{
    Clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveNodeMapping(unsigned index) const
{
    assert(index < this->mNodes.size());
    return index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveElementMapping(unsigned index) const
{
    assert(index < this->mElements.size());
    return index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::SolveBoundaryElementMapping(unsigned index) const
{
    ///\todo sort out boundary elements in a vertex mesh (#1263)
    //    assert(index < this->mBoundaryElements.size() );
    return index;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetRosetteRankOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get pointer to this element
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    // Loop over nodes in the current element and find which is contained in the most elements
    unsigned rosette_rank = 0;
    for (unsigned node_idx = 0; node_idx < p_element->GetNumNodes(); node_idx++)
    {
        unsigned num_elems_this_node = p_element->GetNode(node_idx)->rGetContainingElementIndices().size();

        if (num_elems_this_node > rosette_rank)
        {
            rosette_rank = num_elems_this_node;
        }
    }

    // Return the rosette rank
    return rosette_rank;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::UpdateElementsFacesMap()
{
}

template <>
void MonolayerVertexMesh<3, 3>::UpdateElementsFacesMap()
{
    mElementsFacesMap.clear();
    for (unsigned element_num = 0; element_num < GetNumElements(); element_num++)
    {
        MonolayerVertexElement<3, 3>* pElement = GetElement(element_num);
        std::vector<unsigned> face_indices_in_element;
        for (unsigned face_num = 0; face_num < pElement->GetNumFaces(); face_num++)
        {
            MonolayerVertexElement<2, 3>* pFace = pElement->GetFace(face_num);
            std::vector<MonolayerVertexElement<2, 3>*>::iterator it = std::find(mFaces.begin(), mFaces.end(), pFace);
            face_indices_in_element.push_back(std::distance(mFaces.begin(), it));
        }
        mElementsFacesMap[pElement] = face_indices_in_element;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::UpdateElementsFacesMapOfElement(unsigned index)
{
}

template <>
void MonolayerVertexMesh<3, 3>::UpdateElementsFacesMapOfElement(unsigned index)
{
    MonolayerVertexElement<3, 3>* pElement = GetElement(index);
    std::vector<unsigned> face_indices_in_element;
    for (unsigned face_num = 0; face_num < pElement->GetNumFaces(); face_num++)
    {
        MonolayerVertexElement<2, 3>* pFace = pElement->GetFace(face_num);
        std::vector<MonolayerVertexElement<2, 3>*>::iterator it = std::find(mFaces.begin(), mFaces.end(), pFace);
        face_indices_in_element.push_back(std::distance(mFaces.begin(), it));
    }
    mElementsFacesMap[pElement] = face_indices_in_element;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::map<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*, std::vector<unsigned> >* MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElementsFacesMap()
{
    return &mElementsFacesMap;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::Clear()
{
    // Delete elements
    for (unsigned i = 0; i < mElements.size(); i++)
    {
        delete mElements[i];
    }
    mElements.clear();

    // Delete faces
    for (unsigned i = 0; i < mFaces.size(); i++)
    {
        delete mFaces[i];
    }
    mFaces.clear();

    // Delete nodes
    for (unsigned i = 0; i < this->mNodes.size(); i++)
    {
        delete this->mNodes[i];
    }
    this->mNodes.clear();

    // Clear element faces map
    this->mElementsFacesMap.clear();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumNodes() const
{
    return this->mNodes.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumElements() const
{
    return mElements.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumAllElements() const
{
    return mElements.size();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mFaces.size();
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MonolayerVertexMesh<1, 1>::GetNumFacesByType(MonolayerVertexElementType face_type) const
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,1>::GetNumFacesByType(MonolayerVertexElementType face_type) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MonolayerVertexMesh<1, 2>::GetNumFacesByType(MonolayerVertexElementType face_type) const
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,2>::GetNumFacesByType(MonolayerVertexElementType face_type) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MonolayerVertexMesh<1, 3>::GetNumFacesByType(MonolayerVertexElementType face_type) const
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,3>::GetNumFacesByType(MonolayerVertexElementType face_type) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNumFacesByType(MonolayerVertexElementType face_type) const
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{

    unsigned num_faces = 0;
    for (unsigned i = 0; i < this->GetNumFaces(); i++)
    {
        if (this->GetFace(i)->GetFaceType() == face_type)
        {
            num_faces++;
        }
    }

    return num_faces;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElement(unsigned index) const
{
    assert(index < mElements.size());
    return mElements[index];
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetFace(unsigned index) const
{
    assert(index < mFaces.size());
    return mFaces[index];
}

/**
 * \cond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */
template <>
MonolayerVertexElement<0, 1>* MonolayerVertexMesh<1, 1>::GetFaceOfType(unsigned index, MonolayerVertexElementType faceType) const
{
    return nullptr;
}

template <>
MonolayerVertexElement<0, 2>* MonolayerVertexMesh<1, 2>::GetFaceOfType(unsigned index, MonolayerVertexElementType faceType) const
{
    return nullptr;
}

template <>
MonolayerVertexElement<0, 3>* MonolayerVertexMesh<1, 3>::GetFaceOfType(unsigned index, MonolayerVertexElementType faceType) const
{
    return nullptr;
}
/**
 * \endcond
 * Get Doxygen to ignore, since it's confused by explicit instantiation of templated methods
 */

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetFaceOfType(unsigned index, MonolayerVertexElementType faceType) const
{
    assert(index < mElements.size());
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = mElements[index];
    unsigned num_faces = p_element->GetNumFaces();
    for (unsigned face_index = 0; face_index < num_faces; face_index++)
    {
        if (p_element->GetFace(face_index)->GetFaceType() == faceType)
        {
            return p_element->GetFace(face_index);
        }
    }
    return nullptr;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetCentroidOfElement(unsigned index)
{
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);
    unsigned num_nodes = p_element->GetNumNodes();

    c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);

    switch (SPACE_DIM)
    {
        case 1:
        {
            centroid = 0.5 * (p_element->GetNodeLocation(0) + p_element->GetNodeLocation(1));
        }
        break;
        case 2:
        {
            double centroid_x = 0;
            double centroid_y = 0;

            // Note that we cannot use GetVolumeOfElement() below as it returns the absolute, rather than signed, area
            double element_signed_area = 0.0;

            // Map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
            c_vector<double, SPACE_DIM> first_node_location = p_element->GetNodeLocation(0);
            c_vector<double, SPACE_DIM> pos_1 = zero_vector<double>(SPACE_DIM);

            // Loop over vertices
            for (unsigned local_index = 0; local_index < num_nodes; local_index++)
            {
                c_vector<double, SPACE_DIM> next_node_location = p_element->GetNodeLocation((local_index + 1) % num_nodes);
                c_vector<double, SPACE_DIM> pos_2 = GetVectorFromAtoB(first_node_location, next_node_location);

                double this_x = pos_1[0];
                double this_y = pos_1[1];
                double next_x = pos_2[0];
                double next_y = pos_2[1];

                double signed_area_term = this_x * next_y - this_y * next_x;

                centroid_x += (this_x + next_x) * signed_area_term;
                centroid_y += (this_y + next_y) * signed_area_term;
                element_signed_area += 0.5 * signed_area_term;

                pos_1 = pos_2;
            }

            assert(element_signed_area != 0.0);

            // Finally, map back and employ GetVectorFromAtoB() to allow for periodicity
            centroid = first_node_location;
            centroid(0) += centroid_x / (6.0 * element_signed_area);
            centroid(1) += centroid_y / (6.0 * element_signed_area);
        }
        break;
        case 3:
        {
            ///\todo compute centroid rather than centre of mass (see #1422)
            for (unsigned local_index = 0; local_index < num_nodes; local_index++)
            {
                centroid += p_element->GetNodeLocation(local_index);
            }
            centroid /= ((double)num_nodes);
        }
        break;
        default:
            NEVER_REACHED;
    }
    return centroid;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<1, 1>::GetLumenVolume(bool symmetric_content)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,1>::GetLumenVolume("
              "index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<1, 2>::GetLumenVolume(bool symmetric_content)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,2>::GetLumenVolume("
              "index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<1, 3>::GetLumenVolume(bool symmetric_content)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,3>::GetLumenVolume("
              "index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<2, 2>::GetLumenVolume(bool symmetric_content)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<2,2>::GetLumenVolume("
              "index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<2, 3>::GetLumenVolume(bool symmetric_content)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<2,3>::GetLumenVolume("
              "index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<3, 3>::GetLumenVolume(bool symmetric_content)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    if (symmetric_content)
    {
        double total_lumen_volume{ 0.0 };

        // volume of tetrahedron = 1/6 * abs(r1 * (r2 x r3)); where r1, r2 and r3 are
        // nodes of the triangle at apical cell side and the system origin lies in the
        // lumen
        for (MonolayerVertexElement<2, 3>* p_face : mFaces)
        {
            if (p_face->GetFaceType() == MonolayerVertexElementType::Apical)
            {

                unsigned number_of_nodes = p_face->GetNumNodes();

                // if face has more than three nodes, triangulate face with the center
                // of face every face then consists number_of_nodes + 1 triangles

                c_vector<double, 3> center_of_face = GetPassiveCenterOfFace(p_face);

                for (unsigned local_index = 0; local_index < number_of_nodes; local_index++)
                {

                    c_vector<double, 3> v_cross = VectorProduct(center_of_face, p_face->GetNodeLocation((local_index + 1) % number_of_nodes));

                    total_lumen_volume += inner_prod(p_face->GetNodeLocation(local_index), v_cross);
                }
            }
        }

        return fabs(total_lumen_volume) / 6.0;
    }
    else
    {
        double total_lumen_volume{ 0.0 };

        for (MonolayerVertexElement<2, 3>* p_face : mFaces)
        {
            if (p_face->GetFaceType() == MonolayerVertexElementType::Apical)
            {

                unsigned number_of_nodes = p_face->GetNumNodes();

                c_vector<double, 3> center_of_face = GetPassiveCenterOfFace(p_face);
                c_vector<double, 3> z_normal = zero_vector<double>(3);
                z_normal[2] = 1.0;

                for (unsigned local_index = 0; local_index < number_of_nodes;
                     local_index++)
                {
                    c_vector<double, 3> v_0 = p_face->GetNodeLocation(local_index);
                    c_vector<double, 3> v_1 = p_face->GetNodeLocation((local_index + 1) % number_of_nodes);

                    c_vector<double, 3> s_0 = v_0 - v_1;
                    c_vector<double, 3> s_1 = v_1 - center_of_face;

                    c_vector<double, 3> v_cross = VectorProduct(s_0, s_1);

                    total_lumen_volume += (v_0[2] + v_1[2] + center_of_face[2]) * inner_prod(z_normal, v_cross);
                }
            }
        }

        return fabs(total_lumen_volume) / 6.0;
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
c_vector<double, 1> MonolayerVertexMesh<1, 1>::CalculateLumenVolGradient(
    unsigned node_index, bool symmetric_content)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,1>::CalculateLumenVolGradient("
              "index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
c_vector<double, 1> MonolayerVertexMesh<1, 2>::CalculateLumenVolGradient(
    unsigned node_index, bool symmetric_content)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,2>::CalculateLumenVolGradient("
              "index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
c_vector<double, 1> MonolayerVertexMesh<1, 3>::CalculateLumenVolGradient(
    unsigned node_index, bool symmetric_content)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,3>::CalculateLumenVolGradient("
              "index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
c_vector<double, 2> MonolayerVertexMesh<2, 2>::CalculateLumenVolGradient(
    unsigned node_index, bool symmetric_content)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<2,2>::CalculateLumenVolGradient("
              "index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
c_vector<double, 2> MonolayerVertexMesh<2, 3>::CalculateLumenVolGradient(
    unsigned node_index, bool symmetric_content)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<2,3>::CalculateLumenVolGradient("
              "index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
c_vector<double, 3> MonolayerVertexMesh<3, 3>::CalculateLumenVolGradient(unsigned node_index, bool symmetric_content)
{
    std::set<unsigned> containing_elem_indices = GetNode(node_index)->rGetContainingElementIndices();

    // Initialize total lumen volume gradient for each node
    c_vector<double, 3> v_final_gradient = zero_vector<double>(3);

    // Iterate over these elements
    for (std::set<unsigned>::iterator iter = containing_elem_indices.begin(); iter != containing_elem_indices.end(); ++iter)
    {

        // Get pointer to this element
        MonolayerVertexElement<3, 3>* p_element = GetElement(*iter);
        unsigned local_node_index_in_elem = p_element->GetNodeLocalIndex(node_index);
        /*
         * We assume that the faces are counter-clockwise oriented looking
         *	from the outside for positively oriented faces.
         */
        c_vector<double, 3> this_node_position = p_element->GetNode(local_node_index_in_elem)->rGetLocation();

        unsigned num_faces = p_element->GetNumFaces();

        bool found_apical_face = false;
        int apical_face_index;
        // find apical face
        for (unsigned face_index = 0; face_index < num_faces; face_index++)
        {
            if (p_element->GetFace(face_index)->GetFaceType() == MonolayerVertexElementType::Apical)
            {
                found_apical_face = true;
                apical_face_index = face_index;
                break;
            }
        }

        if (!found_apical_face)
        {
            std::cout << "DIDN'T FIND APICAL FACE IN ELEMENT" << std::endl;
            return v_final_gradient;
        }

        // Is node i in the apical face
        bool is_node_in_apical_face = false;
        for (unsigned face_node_index = 0; face_node_index < p_element->GetFace(apical_face_index)->GetNumNodes(); face_node_index++)
        {
            if (p_element->GetFace(apical_face_index)->GetNode(face_node_index)->GetIndex() == node_index)
            {
                is_node_in_apical_face = true;
                break;
            }
        }

        if (!is_node_in_apical_face)
        {
            return v_final_gradient;
        }

        MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(apical_face_index); // get the apical face
        unsigned num_nodes_in_face = p_face->GetNumNodes();

        c_vector<double, 3> center_of_face = GetPassiveCenterOfFace(p_face);

        // Get local node index in face
        unsigned local_node_index_in_face = p_face->GetNodeLocalIndex(node_index);

        if (symmetric_content)
        {

            // triangle clockwise to v_0
            v_final_gradient += VectorProduct(center_of_face, p_face->GetNodeLocation((local_node_index_in_face + 1) % num_nodes_in_face)) - 1.0 / (double)num_nodes_in_face * VectorProduct(p_face->GetNodeLocation(local_node_index_in_face), p_face->GetNodeLocation((local_node_index_in_face + 1) % num_nodes_in_face));
            // triangle counter clockwise to v_0
            v_final_gradient += VectorProduct(p_face->GetNodeLocation((local_node_index_in_face + num_nodes_in_face - 1) % num_nodes_in_face), center_of_face) + 1.0 / (double)num_nodes_in_face * VectorProduct(p_face->GetNodeLocation(local_node_index_in_face), p_face->GetNodeLocation((local_node_index_in_face + num_nodes_in_face - 1) % num_nodes_in_face));

            for (unsigned node_index = 0; node_index < (num_nodes_in_face - 2); node_index++)
            {
                v_final_gradient -= 1.0 / (double)num_nodes_in_face * VectorProduct(p_face->GetNodeLocation((local_node_index_in_face + node_index + 1) % num_nodes_in_face), p_face->GetNodeLocation((local_node_index_in_face + node_index + 2) % num_nodes_in_face));
            }
        }
        else
        {
            c_vector<double, 3> z_normal = zero_vector<double>(3);
            z_normal[2] = 1.0;
            // triangle clockwise to v_0
            c_vector<double, 3> v_0 = p_face->GetNodeLocation(local_node_index_in_face);
            c_vector<double, 3> v_next = p_face->GetNodeLocation((local_node_index_in_face + 1) % num_nodes_in_face);

            c_vector<double, 3> node_to_center = this->GetVectorFromAtoB(v_0, center_of_face); // s_0
            c_vector<double, 3> center_to_next = this->GetVectorFromAtoB(center_of_face, v_next); // s_1

            c_vector<double, 3> v_cross_1 = VectorProduct(node_to_center, center_to_next); // s_0 x s_1
            c_vector<double, 3> v_cross_2 = VectorProduct(z_normal, center_to_next); // k x s_1
            c_vector<double, 3> v_cross_3 = VectorProduct(node_to_center, z_normal); // s_0 x k

            c_vector<double, 3> gradient_next = (inner_prod(z_normal, v_cross_1) * z_normal * (1.0 + 1.0 / (double)num_nodes_in_face)) + (v_0[2] + v_next[2] + center_of_face[2]) * ((1.0 - 1.0 / (double)num_nodes_in_face) * v_cross_2 + 1.0 / (double)num_nodes_in_face * v_cross_3);

            v_final_gradient += gradient_next;

            // "previous" triangle
            c_vector<double, 3> v_previous = p_face->GetNodeLocation((local_node_index_in_face + num_nodes_in_face - 1) % num_nodes_in_face);

            c_vector<double, 3> node_to_previous = this->GetVectorFromAtoB(v_0, v_previous); // s_0
            c_vector<double, 3> previous_to_center = this->GetVectorFromAtoB(v_previous, center_of_face); // s_1

            c_vector<double, 3> v_cross_4 = VectorProduct(node_to_previous, previous_to_center); // s_0 x s_1
            c_vector<double, 3> v_cross_5 = VectorProduct(z_normal, previous_to_center); // k x s_1
            c_vector<double, 3> v_cross_6 = VectorProduct(z_normal, node_to_previous); // k x s_0

            c_vector<double, 3> gradient_previous = (inner_prod(z_normal, v_cross_4) * z_normal * (1.0 + 1.0 / (double)num_nodes_in_face)) + (v_0[2] + v_previous[2] + center_of_face[2]) * (v_cross_5 + 1.0 / (double)num_nodes_in_face * v_cross_6);

            v_final_gradient += gradient_previous;

            for (unsigned node_index = 0; node_index < (num_nodes_in_face - 2); node_index++)
            {
                c_vector<double, 3> v_1 = p_face->GetNodeLocation((local_node_index_in_face + node_index + 1) % num_nodes_in_face);
                c_vector<double, 3> v_2 = p_face->GetNodeLocation((local_node_index_in_face + node_index + 2) % num_nodes_in_face);

                c_vector<double, 3> one_to_center = this->GetVectorFromAtoB(v_1, center_of_face); // s_0
                c_vector<double, 3> center_to_two = this->GetVectorFromAtoB(center_of_face, v_2); // s_1
                c_vector<double, 3> two_to_one = this->GetVectorFromAtoB(v_2, v_1);

                c_vector<double, 3> v_cross_7 = VectorProduct(one_to_center, center_to_two);
                c_vector<double, 3> v_cross_8 = VectorProduct(z_normal, two_to_one);

                v_final_gradient += 1.0 / (double)num_nodes_in_face * (z_normal * inner_prod(z_normal, v_cross_7) + (v_1[2] + center_of_face[2] + v_2[2]) * v_cross_8);
            }
        }
    }
    v_final_gradient /= 6.0;
    return v_final_gradient;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodeIndices(unsigned nodeIndex)
{
    // Create a set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices;

    // Find the indices of the elements owned by this node
    std::set<unsigned> containing_elem_indices = this->GetNode(nodeIndex)->rGetContainingElementIndices();

    // Iterate over these elements
    for (std::set<unsigned>::iterator elem_iter = containing_elem_indices.begin(); elem_iter != containing_elem_indices.end(); ++elem_iter)
    {
        // Find the local index of this node in this element
        unsigned local_index = GetElement(*elem_iter)->GetNodeLocalIndex(nodeIndex);

        // Find the global indices of the preceding and successive nodes in this element
        unsigned num_nodes = GetElement(*elem_iter)->GetNumNodes();
        unsigned previous_local_index = (local_index + num_nodes - 1) % num_nodes;
        unsigned next_local_index = (local_index + 1) % num_nodes;

        // Add the global indices of these two nodes to the set of neighbouring node indices
        neighbouring_node_indices.insert(GetElement(*elem_iter)->GetNodeGlobalIndex(previous_local_index));
        neighbouring_node_indices.insert(GetElement(*elem_iter)->GetNodeGlobalIndex(next_local_index));
    }

    return neighbouring_node_indices;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringNodeNotAlsoInElement(unsigned nodeIndex, unsigned elemIndex)
{
    // Get a pointer to this element
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(elemIndex);

    // Get the indices of the nodes contained in this element
    std::set<unsigned> node_indices_in_this_element;
    for (unsigned local_index = 0; local_index < p_element->GetNumNodes(); local_index++)
    {
        unsigned global_index = p_element->GetNodeGlobalIndex(local_index);
        node_indices_in_this_element.insert(global_index);
    }

    // Check that the node is in fact contained in the element
    if (node_indices_in_this_element.find(nodeIndex) == node_indices_in_this_element.end())
    {
        EXCEPTION("The given node is not contained in the given element.");
    }

    // Create a set of neighbouring node indices
    std::set<unsigned> neighbouring_node_indices_not_in_this_element;

    // Get the indices of this node's neighbours
    std::set<unsigned> node_neighbours = GetNeighbouringNodeIndices(nodeIndex);

    // Check if each neighbour is also in this element; if not, add it to the set
    for (std::set<unsigned>::iterator iter = node_neighbours.begin();
         iter != node_neighbours.end();
         ++iter)
    {
        if (node_indices_in_this_element.find(*iter) == node_indices_in_this_element.end())
        {
            neighbouring_node_indices_not_in_this_element.insert(*iter);
        }
    }

    return neighbouring_node_indices_not_in_this_element;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetNeighbouringElementIndices(unsigned elementIndex)
{
    // Get a pointer to this element
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = this->GetElement(elementIndex);

    // Create a set of neighbouring element indices
    std::set<unsigned> neighbouring_element_indices;

    // Loop over nodes owned by this element
    for (unsigned local_index = 0; local_index < p_element->GetNumNodes(); local_index++)
    {
        // Get a pointer to this node
        Node<SPACE_DIM>* p_node = p_element->GetNode(local_index);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_node->rGetContainingElementIndices();

        // Form the union of this set with the current set of neighbouring element indices
        std::set<unsigned> all_elements;
        std::set_union(neighbouring_element_indices.begin(), neighbouring_element_indices.end(),
                       containing_elem_indices.begin(), containing_elem_indices.end(),
                       std::inserter(all_elements, all_elements.begin()));

        // Update the set of neighbouring element indices
        neighbouring_element_indices = all_elements;
    }

    // Lastly remove this element's index from the set of neighbouring element indices
    neighbouring_element_indices.erase(elementIndex);

    return neighbouring_element_indices;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>* MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetMeshForVtk()
{
    return this;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MonolayerVertexMesh<1, 1>::ConstructFromMeshReader(AbstractMeshReader<1, 1>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,1>::ConstructFromMeshReader() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MonolayerVertexMesh<1, 2>::ConstructFromMeshReader(AbstractMeshReader<1, 2>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,2>::ConstructFromMeshReader() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MonolayerVertexMesh<1, 3>::ConstructFromMeshReader(AbstractMeshReader<1, 3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,3>::ConstructFromMeshReader() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MonolayerVertexMesh<2, 3>::ConstructFromMeshReader(AbstractMeshReader<2, 3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<2,3>::ConstructFromMeshReader() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MonolayerVertexMesh<2, 2>::ConstructFromMeshReader(AbstractMeshReader<2, 2>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    assert(rMeshReader.HasNodePermutation() == false);
    // Store numbers of nodes and elements
    unsigned num_nodes = rMeshReader.GetNumNodes();
    unsigned num_elements = rMeshReader.GetNumElements();

    // Reserve memory for nodes
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    // Add nodes
    std::vector<double> node_data;
    for (unsigned i = 0; i < num_nodes; i++)
    {
        node_data = rMeshReader.GetNextNode();
        unsigned is_boundary_node = (bool)node_data[2];
        node_data.pop_back();
        this->mNodes.push_back(new Node<2>(i, node_data, is_boundary_node));
    }

    rMeshReader.Reset();

    // Reserve memory for elements
    mElements.reserve(rMeshReader.GetNumElements());

    // Add elements
    for (unsigned elem_index = 0; elem_index < num_elements; elem_index++)
    {
        // Get the data for this element
        ElementData element_data = rMeshReader.GetNextElementData();

        // Get the nodes owned by this element
        std::vector<Node<2>*> nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j = 0; j < num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        // Use nodes and index to construct this element
        // Different node and vertex types not implemented yet
        std::vector<MonolayerVertexElementType> node_types = std::vector<MonolayerVertexElementType>(nodes.size(), MonolayerVertexElementType::Undetermined);
        MonolayerVertexElement<2, 2>* p_element = new MonolayerVertexElement<2, 2>(elem_index, MonolayerVertexElementType::Undetermined, nodes, node_types);
        mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = (unsigned)element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MonolayerVertexMesh<3, 3>::ConstructFromMeshReader(AbstractMeshReader<3, 3>& rMeshReader)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    assert(rMeshReader.HasNodePermutation() == false);

    // Store numbers of nodes and elements
    unsigned num_nodes = rMeshReader.GetNumNodes();
    unsigned num_elements = rMeshReader.GetNumElements();

    // Reserve memory for nodes
    this->mNodes.reserve(num_nodes);

    rMeshReader.Reset();

    // Add nodes
    std::vector<double> node_data;
    for (unsigned i = 0; i < num_nodes; i++)
    {
        node_data = rMeshReader.GetNextNode();
        unsigned is_boundary_node = (bool)node_data[3];
        node_data.pop_back();
        this->mNodes.push_back(new Node<3>(i, node_data, is_boundary_node));
    }

    rMeshReader.Reset();

    // Reserve memory for nodes
    mElements.reserve(rMeshReader.GetNumElements());

    // Use a std::set to keep track of which faces have been added to mFaces
    std::set<unsigned> faces_counted;

    // Add elements
    for (unsigned elem_index = 0; elem_index < num_elements; elem_index++)
    {
        ///\todo Horrible hack! (#1076/#1377)
        typedef MonolayerVertexMeshReader<3, 3> VERTEX_MESH_READER;
        assert(dynamic_cast<VERTEX_MESH_READER*>(&rMeshReader) != nullptr);

        // Get the data for this element
        VertexElementData element_data = static_cast<VERTEX_MESH_READER*>(&rMeshReader)->GetNextElementDataWithFaces();

        // Get the nodes owned by this element
        std::vector<Node<3>*> nodes;
        unsigned num_nodes_in_element = element_data.NodeIndices.size();
        for (unsigned j = 0; j < num_nodes_in_element; j++)
        {
            assert(element_data.NodeIndices[j] < this->mNodes.size());
            nodes.push_back(this->mNodes[element_data.NodeIndices[j]]);
        }

        // Get the faces owned by this element
        std::vector<MonolayerVertexElement<2, 3>*> faces;
        unsigned num_faces_in_element = element_data.Faces.size();
        for (unsigned i = 0; i < num_faces_in_element; i++)
        {
            // Get the data for this face
            ElementData face_data = element_data.Faces[i];

            // Get the face index
            unsigned face_index = (unsigned)face_data.AttributeValue;

            // Get the nodes owned by this face
            std::vector<Node<3>*> nodes_in_face;
            unsigned num_nodes_in_face = face_data.NodeIndices.size();
            for (unsigned j = 0; j < num_nodes_in_face; j++)
            {
                assert(face_data.NodeIndices[j] < this->mNodes.size());
                nodes_in_face.push_back(this->mNodes[face_data.NodeIndices[j]]);
            }

            // If this face index is not already encountered, create a new face and update faces_counted...
            if (faces_counted.find(face_index) == faces_counted.end())
            {
                // Different node and vertex types not implemented yet
                std::vector<MonolayerVertexElementType> face_node_types = std::vector<MonolayerVertexElementType>(nodes_in_face.size(), MonolayerVertexElementType::Undetermined);

                // Use nodes and index to construct this face
                MonolayerVertexElement<2, 3>* p_face = new MonolayerVertexElement<2, 3>(face_index, MonolayerVertexElementType::Undetermined, nodes_in_face, face_node_types);
                mFaces.push_back(p_face);
                faces_counted.insert(face_index);
                faces.push_back(p_face);
            }
            else
            {
                // ... otherwise use the member of mFaces with this index
                bool face_added = false;
                for (unsigned k = 0; k < mFaces.size(); k++)
                {
                    if (mFaces[k]->GetIndex() == face_index)
                    {
                        faces.push_back(mFaces[k]);
                        face_added = true;
                        break;
                    }
                }
                UNUSED_OPT(face_added);
                assert(face_added == true);
            }
        }

        ///\todo Store orientations? (#1076/#1377)
        // Different node and vertex types not implemented yet
        std::vector<bool> orientations = std::vector<bool>(num_faces_in_element, true);
        std::vector<MonolayerVertexElementType> node_types = std::vector<MonolayerVertexElementType>(nodes.size(), MonolayerVertexElementType::Undetermined);

        // Use faces and index to construct this element
        MonolayerVertexElement<3, 3>* p_element = new MonolayerVertexElement<3, 3>(elem_index, MonolayerVertexElementType::Undetermined, faces, orientations, nodes, node_types);
        mElements.push_back(p_element);

        if (rMeshReader.GetNumElementAttributes() > 0)
        {
            assert(rMeshReader.GetNumElementAttributes() == 1);
            unsigned attribute_value = element_data.AttributeValue;
            p_element->SetAttribute(attribute_value);
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(
    const c_vector<double, SPACE_DIM>& rLocationA, const c_vector<double, SPACE_DIM>& rLocationB)
{
    c_vector<double, SPACE_DIM> vector;
    vector = AbstractMesh<ELEMENT_DIM, SPACE_DIM>::GetVectorFromAtoB(rLocationA, rLocationB);

    return vector;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<1, 1>::GetVolumeOfElement(unsigned index)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,1>::GetVolumeOfElement(index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<1, 2>::GetVolumeOfElement(unsigned index)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,2>::GetVolumeOfElement(index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<1, 3>::GetVolumeOfElement(unsigned index)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMesh<1,3>::GetVolumeOfElement(index) is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<2, 2>::GetVolumeOfElement(unsigned index)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    assert(2 == 2 || 2 == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get pointer to this element
    MonolayerVertexElement<2, 2>* p_element = GetElement(index);

    double element_volume = 0.0;

    // We map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
    c_vector<double, 2> first_node_location = p_element->GetNodeLocation(0);
    c_vector<double, 2> pos_1 = zero_vector<double>(2);

    unsigned num_nodes = p_element->GetNumNodes();
    for (unsigned local_index = 0; local_index < num_nodes; local_index++)
    {
        c_vector<double, 2> next_node_location = p_element->GetNodeLocation((local_index + 1) % num_nodes);
        c_vector<double, 2> pos_2 = GetVectorFromAtoB(first_node_location, next_node_location);

        double this_x = pos_1[0];
        double this_y = pos_1[1];
        double next_x = pos_2[0];
        double next_y = pos_2[1];

        element_volume += 0.5 * (this_x * next_y - next_x * this_y);

        pos_1 = pos_2;
    }

    // We take the absolute value just in case the nodes were really oriented clockwise
    return fabs(element_volume);
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<2, 3>::GetVolumeOfElement(unsigned index)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    assert(3 == 2 || 3 == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get pointer to this element
    MonolayerVertexElement<2, 3>* p_element = GetElement(index);

    double element_volume = 0.0;

    assert(2 > 1); // LCOV_EXCL_LINE - code will be removed at compile time
    // Loop over faces and add up pyramid volumes
    // apex is chosen as passive center of first apical face if found, else first face
    // We calculate the volume of the triangulation with passive central nodes as in Krajnc et al. PRE 2018
    unsigned num_faces = p_element->GetNumFaces();
    bool found_apical_face = false;
    c_vector<double, 3> pyramid_apex;
    for (unsigned face_index = 0; face_index < num_faces; face_index++)
    {
        // Chose apex as passive center of face for first apical face
        if (p_element->GetFace(face_index)->GetFaceType() == MonolayerVertexElementType::Apical)
        {
            pyramid_apex = GetPassiveCenterOfFace(p_element->GetFace(face_index));
            found_apical_face = true;
        }
    }
    if (!found_apical_face)
    {
        pyramid_apex = GetPassiveCenterOfFace(p_element->GetFace(0));
    }

    for (unsigned face_index = 0; face_index < num_faces; face_index++)
    {
        // Get pointer to face
        MonolayerVertexElement<1, 3>* p_face = p_element->GetFace(face_index);

        c_vector<double, 3> passive_center_face = GetPassiveCenterOfFace(p_face);

        // Calculate the normals and areas of triangles of face
        std::vector<c_vector<double, 3> > face_triangle_normals;
        std::vector<double> face_triangle_areas = CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &face_triangle_normals);
        assert(face_triangle_areas.size() == face_triangle_normals.size());

        // Distance of passive center (i.e. one triangle corner) to apex
        c_vector<double, 3> passive_center_to_apex = GetVectorFromAtoB(passive_center_face, pyramid_apex);

        for (unsigned triangle_index = 0; triangle_index < face_triangle_normals.size(); triangle_index++)
        {
            // Calculate the perpendicular distance from the triangle to the chosen apex
            double perpendicular_distance = fabs(inner_prod(passive_center_to_apex, face_triangle_normals[triangle_index]));

            // Use triangle area and distance to calculate the volume of the pyramid formed by the triangle
            // and the point pyramid_apex
            element_volume += face_triangle_areas[triangle_index] * perpendicular_distance / 3;
        }
    }

    // We take the absolute value just in case the nodes were really oriented clockwise
    return fabs(element_volume);
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<3, 3>::GetVolumeOfElement(unsigned index)
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    assert(3 == 2 || 3 == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get pointer to this element
    MonolayerVertexElement<3, 3>* p_element = GetElement(index);

    double element_volume = 0.0;

    assert(2 > 1); // LCOV_EXCL_LINE - code will be removed at compile time
    // Loop over faces and add up pyramid volumes
    // apex is chosen as passive center of first apical face if found, else first face
    // We calculate the volume of the triangulation with passive central nodes as in Krajnc et al. PRE 2018
    unsigned num_faces = p_element->GetNumFaces();
    bool found_apical_face = false;
    c_vector<double, 3> pyramid_apex;
    for (unsigned face_index = 0; face_index < num_faces; face_index++)
    {
        // Chose apex as passive center of face for first apical face
        if (p_element->GetFace(face_index)->GetFaceType() == MonolayerVertexElementType::Apical)
        {
            pyramid_apex = GetPassiveCenterOfFace(p_element->GetFace(face_index));
            found_apical_face = true;
        }
    }
    if (!found_apical_face)
    {
        pyramid_apex = GetPassiveCenterOfFace(p_element->GetFace(0));
    }

    for (unsigned face_index = 0; face_index < num_faces; face_index++)
    {
        // Get pointer to face
        MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(face_index);

        c_vector<double, 3> passive_center_face = GetPassiveCenterOfFace(p_face);

        // Calculate the normals and areas of triangles of face
        std::vector<c_vector<double, 3> > face_triangle_normals;
        std::vector<double> face_triangle_areas = CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &face_triangle_normals);
        assert(face_triangle_areas.size() == face_triangle_normals.size());

        // Distance of passive center (i.e. one triangle corner) to apex
        c_vector<double, 3> passive_center_to_apex = GetVectorFromAtoB(passive_center_face, pyramid_apex);

        for (unsigned triangle_index = 0; triangle_index < face_triangle_normals.size(); triangle_index++)
        {
            // Calculate the perpendicular distance from the triangle to the chosen apex
            double perpendicular_distance = fabs(inner_prod(passive_center_to_apex, face_triangle_normals[triangle_index]));

            // Use triangle area and distance to calculate the volume of the pyramid formed by the triangle
            // and the point pyramid_apex
            element_volume += face_triangle_areas[triangle_index] * perpendicular_distance / 3;
        }
    }

    // We take the absolute value just in case the nodes were really oriented clockwise
    return fabs(element_volume);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetVolumeOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get pointer to this element
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    double element_volume = 0.0;
    if (SPACE_DIM == 2)
    {
        // We map the first vertex to the origin and employ GetVectorFromAtoB() to allow for periodicity
        c_vector<double, SPACE_DIM> first_node_location = p_element->GetNodeLocation(0);
        c_vector<double, SPACE_DIM> pos_1 = zero_vector<double>(SPACE_DIM);

        unsigned num_nodes = p_element->GetNumNodes();
        for (unsigned local_index = 0; local_index < num_nodes; local_index++)
        {
            c_vector<double, SPACE_DIM> next_node_location = p_element->GetNodeLocation((local_index + 1) % num_nodes);
            c_vector<double, SPACE_DIM> pos_2 = GetVectorFromAtoB(first_node_location, next_node_location);

            double this_x = pos_1[0];
            double this_y = pos_1[1];
            double next_x = pos_2[0];
            double next_y = pos_2[1];

            element_volume += 0.5 * (this_x * next_y - next_x * this_y);

            pos_1 = pos_2;
        }
    }
    else
    {
        assert(ELEMENT_DIM > 1); // LCOV_EXCL_LINE - code will be removed at compile time
                                 // Loop over faces and add up pyramid volumes
        // apex is chosen as passive center of first apical face if found, else first face
        // We calculate the volume of the triangulation with passive central nodes as in Krajnc et al. PRE 2018
        unsigned num_faces = p_element->GetNumFaces();
        bool found_apical_face = false;
        c_vector<double, SPACE_DIM> pyramid_apex;
        for (unsigned face_index = 0; face_index < num_faces; face_index++)
        {
            // Chose apex as passive center of face for first apical face
            if (p_element->GetFace(face_index)->GetFaceType() == MonolayerVertexElementType::Apical)
            {
                pyramid_apex = GetPassiveCenterOfFace(p_element->GetFace(face_index));
                found_apical_face = true;
                break;
            }
        }
        if (!found_apical_face)
        {
            pyramid_apex = GetPassiveCenterOfFace(p_element->GetFace(0));
        }

        for (unsigned face_index = 0; face_index < num_faces; face_index++)
        {
            // Get pointer to face
            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = p_element->GetFace(face_index);

            c_vector<double, SPACE_DIM> passive_center_face = GetPassiveCenterOfFace(p_face);

            // Calculate the normals and areas of triangles of face
            std::vector<c_vector<double, SPACE_DIM> > face_triangle_normals;
            std::vector<double> face_triangle_areas = CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &face_triangle_normals);
            assert(face_triangle_areas.size() == face_triangle_normals.size());

            // Distance of passive center (i.e. one triangle corner) to apex
            c_vector<double, SPACE_DIM> passive_center_to_apex = GetVectorFromAtoB(passive_center_face, pyramid_apex);

            for (unsigned triangle_index = 0; triangle_index < face_triangle_normals.size(); triangle_index++)
            {
                // Calculate the perpendicular distance from the triangle to the chosen apex
                double perpendicular_distance = fabs(inner_prod(passive_center_to_apex, face_triangle_normals[triangle_index]));

                // Use triangle area and distance to calculate the volume of the pyramid formed by the triangle
                // and the point pyramid_apex
                element_volume += face_triangle_areas[triangle_index] * perpendicular_distance / 3;
            }
        }
    }
    // We take the absolute value just in case the nodes were really oriented clockwise
    return fabs(element_volume);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetSurfaceAreaOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get pointer to this element
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    double surface_area = 0.0;
    if (SPACE_DIM == 2)
    {
        unsigned num_nodes = p_element->GetNumNodes();
        unsigned this_node_index = p_element->GetNodeGlobalIndex(0);
        for (unsigned local_index = 0; local_index < num_nodes; local_index++)
        {
            unsigned next_node_index = p_element->GetNodeGlobalIndex((local_index + 1) % num_nodes);

            surface_area += this->GetDistanceBetweenNodes(this_node_index, next_node_index);
            this_node_index = next_node_index;
        }
    }
    else
    {
        // Loop over faces and add up areas
        for (unsigned face_index = 0; face_index < p_element->GetNumFaces(); face_index++)
        {
            surface_area += CalculateAreaOfFace(p_element->GetFace(face_index));
        }
    }
    return surface_area;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetMidPlaneAreaOfElement(unsigned index)
{
    assert(SPACE_DIM == 2 || SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get pointer to this element
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(index);

    if (SPACE_DIM == 2)
    {
        return GetSurfaceAreaOfElement(index);
    }
    else
    {
        std::vector<std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > > all_mid_points;
        for (unsigned face_index = 0; face_index < p_element->GetNumFaces(); face_index++)
        {
            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = GetFace(face_index);
            if (p_face->GetFaceType() != MonolayerVertexElementType::Lateral)
                continue;

            // For lateral faces we determine the apico-basal mid points

            Node<SPACE_DIM>* p_basal_node_a = nullptr;
            Node<SPACE_DIM>* p_basal_node_b = nullptr;

            Node<SPACE_DIM>* p_apical_node_a = nullptr;
            Node<SPACE_DIM>* p_apical_node_b = nullptr;

            unsigned num_nodes = p_face->GetNumNodes();

            // Caveat: We assume that we only have one apical and one basal edge
            for (unsigned local_index = 0; local_index < num_nodes; local_index++)
            {
                Node<SPACE_DIM>* p_current_node = p_face->GetNode(local_index);
                Node<SPACE_DIM>* p_next_node = p_face->GetNode((local_index + 1) % num_nodes);

                MonolayerVertexElementType current_node_type = p_face->GetNodeType(local_index);
                MonolayerVertexElementType next_node_type = p_face->GetNodeType((local_index + 1) % num_nodes);

                // Get one of the lateral edges
                if (current_node_type == MonolayerVertexElementType::Basal && next_node_type == MonolayerVertexElementType::Apical)
                {
                    p_basal_node_a = p_current_node;
                    p_apical_node_a = p_next_node;
                }

                // Get other lateral edge
                if (current_node_type == MonolayerVertexElementType::Apical && next_node_type == MonolayerVertexElementType::Basal)
                {
                    p_apical_node_b = p_current_node;
                    p_basal_node_b = p_next_node;
                }
            }
            assert(p_basal_node_a != nullptr && p_apical_node_a != nullptr); // LCOV_EXCL_LINE - code will be removed at compile time
            assert(p_basal_node_b != nullptr && p_apical_node_b != nullptr); // LCOV_EXCL_LINE - code will be removed at compile time

            c_vector<double, 3> mid_point_a = (p_basal_node_a->rGetLocation() + p_apical_node_a->rGetLocation()) / 2.0;
            c_vector<double, 3> mid_point_b = (p_basal_node_b->rGetLocation() + p_apical_node_b->rGetLocation()) / 2.0;

            all_mid_points.push_back(std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> >(mid_point_a, mid_point_b));
        }
        // Now determine passive center for mid-plane triangulation
        unsigned num_points = 0.0;
        c_vector<double, SPACE_DIM> passive_center = zero_vector<double>(SPACE_DIM);
        for (auto iter = all_mid_points.begin(); iter != all_mid_points.end(); ++iter)
        {
            num_points += 2;
            passive_center += iter->first;
            passive_center += iter->second;
        }
        passive_center /= num_points;

        // Finally, determine the area by adding triangulation areas
        double mid_plane_area = 0.0;
        for (auto iter = all_mid_points.begin(); iter != all_mid_points.end(); ++iter)
        {
            c_vector<double, SPACE_DIM> ab_1 = GetVectorFromAtoB(iter->first, passive_center);
            c_vector<double, SPACE_DIM> ab_2 = GetVectorFromAtoB(iter->second, passive_center);

            mid_plane_area += norm_2(VectorProduct(ab_1, ab_2)) / 2.0;
        }
        return mid_plane_area;
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMesh<1, 3>::GetMidPlaneAreaOfElement(unsigned index)
{
    EXCEPTION("Mid plane area of elements only in 2D & 3D.");
}

template <>
double MonolayerVertexMesh<1, 2>::GetMidPlaneAreaOfElement(unsigned index)
{
    EXCEPTION("Mid plane area of elements only in 2D & 3D.");
}

template <>
double MonolayerVertexMesh<1, 1>::GetMidPlaneAreaOfElement(unsigned index)
{
    EXCEPTION("Mid plane area of elements only in 2D & 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <>
double MonolayerVertexMesh<1, 1>::GetThicknessOfElement(unsigned elementIndex)
{
    EXCEPTION("GetThicknessOfElement not defined in one dimension");
}

template <>
double MonolayerVertexMesh<1, 2>::GetThicknessOfElement(unsigned elementIndex)
{
    EXCEPTION("GetThicknessOfElement not defined in one dimension");
}

template <>
double MonolayerVertexMesh<1, 3>::GetThicknessOfElement(unsigned elementIndex)
{
    EXCEPTION("GetThicknessOfElement not defined in one dimension");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetThicknessOfElement(unsigned elementIndex)
{
    // Iterate over faces to find apical and basal face
    bool found_apical = false;
    bool found_basal = false;

    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = GetElement(elementIndex);
    unsigned num_faces = p_element->GetNumFaces();

    c_vector<double, SPACE_DIM> apical_position, basal_position;
    for (unsigned index_face = 0; index_face < num_faces; ++index_face)
    {
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = p_element->GetFace(index_face);
        if (p_face->GetFaceType() == MonolayerVertexElementType::Apical)
        {
            apical_position = GetPassiveCenterOfFace(p_face);
            found_apical = true;
        }
        else if (p_face->GetFaceType() == MonolayerVertexElementType::Basal)
        {
            basal_position = GetPassiveCenterOfFace(p_face);
            found_basal = true;
        }
    }
    if (found_apical & found_basal)
    {
        return norm_2(GetVectorFromAtoB(apical_position, basal_position));
    }
    else
    {
        return 0.0;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetAreaGradientOfFaceAtNode(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, unsigned localIndex)
{
    EXCEPTION("GetAreaGradientOfFaceAtNode() only in 3 dimensions.");
}

template <>
c_vector<double, 3> MonolayerVertexMesh<3, 3>::GetAreaGradientOfFaceAtNode(MonolayerVertexElement<2, 3>* pFace, unsigned localIndex)
{
    /** This is a wrong old version
            unsigned num_nodes_in_element = pFace->GetNumNodes();
unsigned next_local_index = (localIndex + 1) % num_nodes_in_element;

// We add an extra num_nodes_in_element in the line below as otherwise this term can be negative, which breaks the % operator
unsigned previous_local_index = (num_nodes_in_element + localIndex - 1) % num_nodes_in_element;

            c_vector<double, 3> this_node_location = pFace->GetNodeLocation(localIndex);
c_vector<double, 3> previous_node_location = pFace->GetNodeLocation(previous_local_index);
c_vector<double, 3> next_node_location = pFace->GetNodeLocation(next_local_index);
            c_vector<double, 3> passive_center_location = GetPassiveCenterOfFace(pFace);

            c_vector<double, 3> vthis_minus_v_pc = this->GetVectorFromAtoB(this_node_location, passive_center_location);

            // triangle with previous
            double a_factor_p = 2.0 * inner_prod(previous_node_location - passive_center_location, previous_node_location - passive_center_location);
            double b_factor_p = 2.0 * ( inner_prod(this_node_location , - previous_node_location + passive_center_location) + inner_prod(passive_center_location , previous_node_location - passive_center_location));
double c_factor_p = 2.0 * ( inner_prod(this_node_location - previous_node_location , previous_node_location) + inner_prod(passive_center_location, previous_node_location - this_node_location));

c_vector<double, 3> vprevious_minus_v_pc = this->GetVectorFromAtoB(previous_node_location, passive_center_location);
            double parallelogram_area_p =  norm_2( VectorProduct(vprevious_minus_v_pc, vthis_minus_v_pc) );

            // triangle with next
            double a_factor_n = 2.0 * inner_prod(next_node_location - passive_center_location, next_node_location - passive_center_location);
            double b_factor_n = 2.0 * ( inner_prod(this_node_location , - next_node_location + passive_center_location) + inner_prod(passive_center_location , next_node_location - passive_center_location));
double c_factor_n = 2.0 * ( inner_prod(this_node_location - next_node_location , next_node_location) + inner_prod(passive_center_location, next_node_location - this_node_location));

c_vector<double, 3> vnext_minus_v_pc = this->GetVectorFromAtoB(next_node_location, passive_center_location);
            double parallelogram_area_n =  norm_2( VectorProduct(vnext_minus_v_pc, vthis_minus_v_pc) );

c_vector<double, 3> area_gradient;

            // grad|*| = grad(|*|^2)/(2 sqrt(|*|^2))
area_gradient += (a_factor_p * this_node_location + b_factor_p * previous_node_location + c_factor_p * passive_center_location)/(2.0 * parallelogram_area_p);
area_gradient += (a_factor_n * this_node_location + b_factor_n * next_node_location + c_factor_n * passive_center_location)/(2.0 * parallelogram_area_n);

            // area of triangle is |*|/2, therefore grad(area)=grad|*|/2
return area_gradient/2.0;
    **/

    // Create empty area gradient for this face
    c_vector<double, 3> area_gradient = zero_vector<double>(3);
    // std::cout << "initial:" << area_gradient[0] << ", " << area_gradient[1] << ", " << area_gradient[2] << ", ";

    // Passive center of face
    unsigned num_nodes_face = pFace->GetNumNodes();
    c_vector<double, 3> this_node_location = pFace->GetNodeLocation(localIndex);
    unsigned this_node_global_index = pFace->GetNode(localIndex)->GetIndex();
    c_vector<double, 3> passive_center_location = GetPassiveCenterOfFace(pFace);

    c_vector<double, 3> v_minus_v_pc = this->GetVectorFromAtoB(passive_center_location, pFace->GetNode(0)->rGetLocation());
    c_vector<double, 3> current_node_location = pFace->GetNodeLocation(0);
    unsigned current_node_global_index = pFace->GetNode(0)->GetIndex();

    // Loop over all triangles and compute the contribution by each triangle
    for (unsigned local_index = 0; local_index < pFace->GetNumNodes(); local_index++)
    {
        c_vector<double, 3> vnext_minus_v_pc = this->GetVectorFromAtoB(passive_center_location, pFace->GetNode((local_index + 1) % pFace->GetNumNodes())->rGetLocation());
        c_vector<double, 3> next_node_location = pFace->GetNodeLocation((local_index + 1) % pFace->GetNumNodes());
        unsigned next_node_global_index = pFace->GetNode((local_index + 1) % pFace->GetNumNodes())->GetIndex();

        // If we now have a triangle with current node being node to be calculated
        if (current_node_global_index == this_node_global_index)
        {
            double a_factor = 2.0 * (inner_prod(next_node_location - passive_center_location, next_node_location - passive_center_location) + 1.0 / double(num_nodes_face) * inner_prod(current_node_location - next_node_location, next_node_location - passive_center_location));

            double b_factor = 2.0 * (inner_prod(current_node_location - passive_center_location, -next_node_location + passive_center_location) + 1.0 / double(num_nodes_face) * inner_prod(passive_center_location - current_node_location, current_node_location - next_node_location));

            double c_factor = 2.0 * (inner_prod(current_node_location - next_node_location, next_node_location - passive_center_location) + 1.0 / double(num_nodes_face) * inner_prod(current_node_location - next_node_location, current_node_location - next_node_location));

            double parallelogram_area = norm_2(VectorProduct(vnext_minus_v_pc, v_minus_v_pc));

            // grad|*| = grad(|*|^2)/(2 sqrt(|*|^2))
            area_gradient[0] += ((a_factor * current_node_location + b_factor * next_node_location + c_factor * passive_center_location) / (2.0 * parallelogram_area))[0];
            area_gradient[1] += ((a_factor * current_node_location + b_factor * next_node_location + c_factor * passive_center_location) / (2.0 * parallelogram_area))[1];
            area_gradient[2] += ((a_factor * current_node_location + b_factor * next_node_location + c_factor * passive_center_location) / (2.0 * parallelogram_area))[2];
        }
        // If we now have a triangle with next node being node to be calculated
        else if (next_node_global_index == this_node_global_index)
        {
            double a_factor = 2.0 * (inner_prod(current_node_location - passive_center_location, current_node_location - passive_center_location) + 1.0 / double(num_nodes_face) * inner_prod(next_node_location - current_node_location, current_node_location - passive_center_location));

            double b_factor = 2.0 * (inner_prod(next_node_location - passive_center_location, -current_node_location + passive_center_location) + 1.0 / double(num_nodes_face) * inner_prod(passive_center_location - next_node_location, next_node_location - current_node_location));

            double c_factor = 2.0 * (inner_prod(next_node_location - current_node_location, current_node_location - passive_center_location) + 1.0 / double(num_nodes_face) * inner_prod(next_node_location - current_node_location, next_node_location - current_node_location));

            double parallelogram_area = norm_2(VectorProduct(vnext_minus_v_pc, v_minus_v_pc));

            // grad|*| = grad(|*|^2)/(2 sqrt(|*|^2))
            area_gradient[0] += ((a_factor * next_node_location + b_factor * current_node_location + c_factor * passive_center_location) / (2.0 * parallelogram_area))[0];
            area_gradient[1] += ((a_factor * next_node_location + b_factor * current_node_location + c_factor * passive_center_location) / (2.0 * parallelogram_area))[1];
            area_gradient[2] += ((a_factor * next_node_location + b_factor * current_node_location + c_factor * passive_center_location) / (2.0 * parallelogram_area))[2];
        }
        // Because the passive center is also moved we have to calculate the energy change from the other triangles from this
        else
        {
            double a_factor = 2.0 / double(num_nodes_face) * inner_prod(current_node_location - next_node_location, next_node_location - passive_center_location);

            double b_factor = 2.0 / double(num_nodes_face) * inner_prod(passive_center_location - current_node_location, current_node_location - next_node_location);

            double c_factor = 2.0 / double(num_nodes_face) * inner_prod(current_node_location - next_node_location, current_node_location - next_node_location);

            double parallelogram_area = norm_2(VectorProduct(vnext_minus_v_pc, v_minus_v_pc));

            // grad|*| = grad(|*|^2)/(2 sqrt(|*|^2))
            area_gradient[0] += ((a_factor * current_node_location + b_factor * next_node_location + c_factor * passive_center_location) / (2.0 * parallelogram_area))[0];
            area_gradient[1] += ((a_factor * current_node_location + b_factor * next_node_location + c_factor * passive_center_location) / (2.0 * parallelogram_area))[1];
            area_gradient[2] += ((a_factor * current_node_location + b_factor * next_node_location + c_factor * passive_center_location) / (2.0 * parallelogram_area))[2];
        }
        v_minus_v_pc = vnext_minus_v_pc;
        current_node_location = next_node_location;
        current_node_global_index = next_node_global_index;
    }
    // area of triangle is |*|/2, therefore grad(area)=grad|*|/\n
    // std::cout << "\nArea Gradient: " << area_gradient[0] << ", " << area_gradient[1] << ", " << area_gradient[2] << ", ";
    area_gradient /= 2.0;
    return area_gradient;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetVolumeGradientAtNode(unsigned elementIndex, unsigned localVertexIndex)
{
    EXCEPTION("GetVolumeGradientAtNode only implemented in 3 dimensions.");
}

template <>
c_vector<double, 3> MonolayerVertexMesh<3, 3>::GetVolumeGradientAtNode(unsigned elementIndex, unsigned localVertexIndex)
{
    // Get pointer to this element
    MonolayerVertexElement<3, 3>* p_element = GetElement(elementIndex);

    c_vector<double, 3> v_final_gradient = zero_vector<double>(3);
    /*
     * We assume that the faces are counter-clockwise oriented looking
     *	from the outside for positively oriented faces.
     */
    c_vector<double, 3> this_node_position = p_element->GetNode(localVertexIndex)->rGetLocation();
    unsigned this_node_global_index = p_element->GetNode(localVertexIndex)->GetIndex();

    // apex is chosen as passive center of first apical face if found, else first face
    // this has to be identical as in the GetVolumeOfElement method!
    unsigned num_faces = p_element->GetNumFaces();
    bool found_apical_face = false;
    c_vector<double, 3> pyramid_apex;
    unsigned apex_face_index;
    for (unsigned face_index = 0; face_index < num_faces; face_index++)
    {
        // Chose apex as passive center of face for first apical face
        if (p_element->GetFace(face_index)->GetFaceType() == MonolayerVertexElementType::Apical)
        {
            pyramid_apex = GetPassiveCenterOfFace(p_element->GetFace(face_index));
            found_apical_face = true;
            apex_face_index = face_index;
            break;
        }
    }
    if (!found_apical_face)
    {
        pyramid_apex = GetPassiveCenterOfFace(p_element->GetFace(0));
        apex_face_index = 0;
    }

    // Is node i in face gamma, which contains pyramid apex
    bool is_node_in_apex_face = false;
    for (unsigned face_node_index = 0; face_node_index < p_element->GetFace(apex_face_index)->GetNumNodes(); face_node_index++)
    {
        if (p_element->GetFace(apex_face_index)->GetNode(face_node_index)->GetIndex() == this_node_global_index)
        {
            is_node_in_apex_face = true;
            break;
        }
    }

    // Loop over faces and triangles
    for (unsigned face_index = 0; face_index < num_faces; face_index++)
    {
        // Get pointer to face
        MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(face_index);
        double orientation = 1.0;
        if (p_element->FaceIsOrientatedClockwise(face_index))
        {
            orientation = -1.0;
        }

        c_vector<double, 3> passive_center_face = GetPassiveCenterOfFace(p_face);

        // Distance of passive center (i.e. one triangle corner) to apex
        c_vector<double, 3> v_passive_center_to_apex = GetVectorFromAtoB(passive_center_face, pyramid_apex);

        // Is node i in face k, which we calculate
        bool is_node_in_face = false;
        for (unsigned face_node_index = 0; face_node_index < p_face->GetNumNodes(); face_node_index++)
        {
            if (p_face->GetNode(face_node_index)->GetIndex() == this_node_global_index)
            {
                is_node_in_face = true;
                break;
            }
        }

        // Loop over all triangles
        unsigned num_nodes_in_face = p_face->GetNumNodes();
        c_vector<double, 3> v_minus_v_pc = this->GetVectorFromAtoB(passive_center_face, p_face->GetNode(0)->rGetLocation());
        unsigned current_node_global_index = p_face->GetNode(0)->GetIndex();

        for (unsigned node_local_index = 0; node_local_index < num_nodes_in_face; node_local_index++)
        {
            c_vector<double, 3> vnext_minus_v_pc = this->GetVectorFromAtoB(passive_center_face, p_face->GetNode((node_local_index + 1) % p_face->GetNumNodes())->rGetLocation());
            unsigned next_node_global_index = p_face->GetNode((node_local_index + 1) % p_face->GetNumNodes())->GetIndex();

            // The three cross products involved in volume calculation
            // (r(l+1)-p) x (r(l)-p)
            // (a-p)			x	(r(l+1)-p)
            // (r(l)-p)		x (a-p)
            c_vector<double, 3> v_next_cross_prev = VectorProduct(vnext_minus_v_pc, v_minus_v_pc);
            c_vector<double, 3> v_apex_cross_next = VectorProduct(v_passive_center_to_apex, vnext_minus_v_pc);
            c_vector<double, 3> v_prev_cross_apex = VectorProduct(v_minus_v_pc, v_passive_center_to_apex);

            // Add the contributions
            if (is_node_in_face && face_index != apex_face_index)
            {
                v_final_gradient -= orientation * 1.0 / double(num_nodes_in_face) * v_next_cross_prev;
                v_final_gradient -= orientation * 1.0 / double(num_nodes_in_face) * v_apex_cross_next;
                v_final_gradient -= orientation * 1.0 / double(num_nodes_in_face) * v_prev_cross_apex;
            }
            // if this_node_global_index is previous/current node index
            if (is_node_in_face && this_node_global_index == current_node_global_index)
            {
                v_final_gradient += orientation * v_apex_cross_next;
            }

            // if this_node_global_index is next node index
            if (is_node_in_face && this_node_global_index == next_node_global_index)
            {
                v_final_gradient += orientation * v_prev_cross_apex;
            }

            // if node is in apex side
            unsigned num_nodes_in_apex_face = p_element->GetFace(apex_face_index)->GetNumNodes();
            if (is_node_in_apex_face && face_index != apex_face_index)
            {
                v_final_gradient += orientation * 1.0 / double(num_nodes_in_apex_face) * v_next_cross_prev;
            }

            v_minus_v_pc = vnext_minus_v_pc;
            current_node_global_index = next_node_global_index;
        }
    }
    // Pyramid volume is area_base*height/3 and triangle size only half the parallelogram
    v_final_gradient /= 6.0;
    return v_final_gradient;
}

//////////////////////////////////////////////////////////////////////
//                        3D-specific methods                       //
//////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetPassiveCenterOfFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace)
{
    unsigned num_nodes = pFace->GetNumNodes();

    c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);

    switch (ELEMENT_DIM - 1)
    {
        case 1:
        {
            centroid = 0.5 * (pFace->GetNodeLocation(0) + pFace->GetNodeLocation(1));
        }
        break;
        case 2:
        {
            // Calculate the mean of the face nodes
            // centroid = pFace->GetNodeLocation(0);

            // Loop over vertices
            for (unsigned local_index = 0; local_index < num_nodes; local_index++)
            {
                centroid += pFace->GetNodeLocation((local_index + 1) % num_nodes);
            }
            centroid /= double(num_nodes);
        }
        break;
        default:
            NEVER_REACHED;
    }
    return centroid;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
c_vector<double, 1> MonolayerVertexMesh<1, 1>::GetPassiveCenterOfFaceTypeInElement(unsigned elementIndex, MonolayerVertexElementType faceType)
{
    EXCEPTION("Passive center of faces only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
c_vector<double, 2> MonolayerVertexMesh<1, 2>::GetPassiveCenterOfFaceTypeInElement(unsigned elementIndex, MonolayerVertexElementType faceType)
{
    EXCEPTION("Passive center of faces only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
c_vector<double, 3> MonolayerVertexMesh<1, 3>::GetPassiveCenterOfFaceTypeInElement(unsigned elementIndex, MonolayerVertexElementType faceType)
{
    EXCEPTION("Passive center of faces only in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/*** Not implemented...
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetVolumeEnclosedByFaces(MonolayerVertexElementType faceType)
{
        // First we determine the apical/basal faces and calculate the centroid, which we will use as apex
        unsigned num_faces = GetNumFaces();
        std::vector<unsigned> faces_of_type;
        unsigned num_faces_of_type = 0;
        c_vector<double, SPACE_DIM> centroid_lumen = zero_vector<double>(SPACE_DIM);
        for(unsigned index = 0; index < num_faces; ++index)
        {
                MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace = GetFace(index);
                if(pFace->GetFaceType() == faceType)
                {
                        faces_of_type.push_back(index);
                        num_faces_of_type++;
                        centroid_lumen += GetPassiveCenterOfFace(pElement);
                }
        }
        centroid_lumen /= num_faces_of_type;

        // We calculate the volume for each face of faceType by calculating the pyramid volume with centroid as apex
        for(auto it = faces_of_type.begin(); it != faces_of_type.end(); ++it)
        {

        }
}
***/

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetPassiveCenterOfFaceTypeInElement(unsigned elementIndex, MonolayerVertexElementType faceType)
{
    c_vector<double, SPACE_DIM> passive_center = zero_vector<double>(SPACE_DIM);
    MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* pElement = GetElement(elementIndex);
    unsigned num_faces = pElement->GetNumFaces();
    for (unsigned index_face = 0; index_face < num_faces; index_face++)
    {
        if (pElement->GetFace(index_face)->GetFaceType() == faceType)
        {
            passive_center = GetPassiveCenterOfFace(pElement->GetFace(index_face));
            break;
        }
    }
    return passive_center;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::CalculateUnitNormalToFaceWithArea(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, c_vector<double, SPACE_DIM>& rNormal)
{
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // As we are in 3D, the face must have at least three vertices
    assert(pFace->GetNumNodes() >= 3u);

    // Reset the answer
    rNormal = zero_vector<double>(SPACE_DIM);
    double surface = 0.0;

    // We calculate the area via the triangulation with passive central nodes as in Krajnc et al. PRE 2018
    c_vector<double, SPACE_DIM> v_passive_center = GetPassiveCenterOfFace(pFace);

    c_vector<double, SPACE_DIM> v_minus_v_pc = this->GetVectorFromAtoB(v_passive_center, pFace->GetNode(0)->rGetLocation());
    for (unsigned local_index = 0; local_index < pFace->GetNumNodes(); local_index++)
    {
        c_vector<double, SPACE_DIM> vnext_minus_v_pc = this->GetVectorFromAtoB(v_passive_center, pFace->GetNode((local_index + 1) % pFace->GetNumNodes())->rGetLocation());

        c_vector<double, SPACE_DIM> localNormal = VectorProduct(v_minus_v_pc, vnext_minus_v_pc);
        surface += norm_2(localNormal) / 2.0;

        // normal direction is implemented as mean normal of triangles (note possible errors from orientation)
        rNormal += localNormal;
        v_minus_v_pc = vnext_minus_v_pc;
    }
    double magnitude = norm_2(rNormal);
    if (magnitude != 0.0)
    {
        // Normalize the normal vector
        rNormal /= magnitude;
        // If all points are co-located, then magnitude==0.0 and there is potential for a floating point exception
        // here if we divide by zero, so we'll move on.
    }
    return surface;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::CalculateAreaOfFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace)
{
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // Get the unit normal to the plane of this face
    c_vector<double, SPACE_DIM> unit_normal;
    return CalculateUnitNormalToFaceWithArea(pFace, unit_normal);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::CalculateAreasAndNormalsOfTriangulationsOfFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, std::vector<c_vector<double, SPACE_DIM> >* pTriangleNormals)
{
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    // As we are in 3D, the face must have at least three vertices
    assert(pFace->GetNumNodes() >= 3u);

    // Create vector to save all triangle areas
    std::vector<double> triangle_surfaces;

    // We calculate the area via the triangulation with passive central nodes as in Krajnc et al. PRE 2018
    c_vector<double, SPACE_DIM> v_passive_center = GetPassiveCenterOfFace(pFace);

    c_vector<double, SPACE_DIM> v_minus_v_pc = this->GetVectorFromAtoB(v_passive_center, pFace->GetNode(0)->rGetLocation());
    for (unsigned local_index = 0; local_index < pFace->GetNumNodes(); local_index++)
    {
        // Calculate normal vector of triangle
        c_vector<double, SPACE_DIM> vnext_minus_v_pc = this->GetVectorFromAtoB(v_passive_center, pFace->GetNode((local_index + 1) % pFace->GetNumNodes())->rGetLocation());
        c_vector<double, SPACE_DIM> localNormal = VectorProduct(v_minus_v_pc, vnext_minus_v_pc);

        // Calculate and save area of triangle
        triangle_surfaces.push_back(norm_2(localNormal) / 2.0);

        v_minus_v_pc = vnext_minus_v_pc;

        // Normalize the local normal vector and push it to the TriangleNormals vector of vectors
        double magnitude = norm_2(localNormal);
        if (magnitude != 0.0)
        {
            // Normalize the normal vector
            localNormal /= magnitude;
            // If all points are co-located, then magnitude==0.0 and there is potential for a floating point exception
            // here if we divide by zero, so we'll move on.
        }
        pTriangleNormals->push_back(localNormal);
    }
    return triangle_surfaces;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, 3> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::CalculateMomentsOfFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace,
                                                                                        c_vector<double, SPACE_DIM> normalVector,
                                                                                        c_vector<double, SPACE_DIM> orthogonalVector1,
                                                                                        c_vector<double, SPACE_DIM> orthogonalVector2)
{
    /* To calculate the Moment we project onto the orthogonal space w.r.t. the mean normal
     * at the passive center and calculate the Moments with respect to a basis of this space,
     * spanned by orthogonalVector1 and orthgonalVector2
     */
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    unsigned num_nodes = pFace->GetNumNodes();
    c_vector<double, 3> passive_center = GetPassiveCenterOfFace(pFace);

    // Project onto orthogonal space at passive center
    std::vector<c_vector<double, 3> > vec_projected_nodes;
    for (unsigned node = 0; node < num_nodes; node++)
    {
        c_vector<double, 3> node_pos = pFace->GetNode(node)->rGetLocation() - passive_center;
        vec_projected_nodes.push_back(node_pos - inner_prod(node_pos, normalVector) * normalVector);
    }
    // Decompose in orthgonal space and calculate Moments in given basis
    c_vector<double, 3> moments = zero_vector<double>(3);

    c_vector<double, 3> curr_proj = vec_projected_nodes[0];
    for (unsigned node = 0; node < num_nodes; node++)
    {
        // Decompose
        c_vector<double, 3> next_proj = vec_projected_nodes[(node + 1) % num_nodes];
        double pos1[2], pos2[2];
        pos1[0] = inner_prod(orthogonalVector1, curr_proj);
        pos1[1] = inner_prod(orthogonalVector2, curr_proj);
        pos2[0] = inner_prod(orthogonalVector1, next_proj);
        pos2[1] = inner_prod(orthogonalVector2, next_proj);

        // Calculate moments
        double signed_area_term = pos1[0] * pos2[1] - pos2[0] * pos1[1];
        // Ixx
        moments(0) += (pos1[1] * pos1[1] + pos1[1] * pos2[1] + pos2[1] * pos2[1]) * signed_area_term;
        // Iyy
        moments(1) += (pos1[0] * pos1[0] + pos1[0] * pos2[0] + pos2[0] * pos2[0]) * signed_area_term;
        // Ixy
        moments(2) += (pos1[0] * pos2[1] + 2.0 * pos1[0] * pos1[1] + 2.0 * pos2[0] * pos2[1] + pos2[0] * pos1[1]) * signed_area_term;

        curr_proj = next_proj;
    }
    moments(0) /= 12.0;
    moments(1) /= 12.0;
    moments(2) /= 24.0;

    /*
     * If the nodes were supplied in a clockwise rather
     * than anticlockwise manner, or if this arose as a result of enforcing
     * periodicity, then our computed quantities will be the wrong sign, so
     * we need to fix this.
     */
    if (moments(0) < 0.0)
    {
        moments(0) = -moments(0);
        moments(1) = -moments(1);
        moments(2) = -moments(2);
    }
    return moments;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetMeanApicalBasalLongAxisOfElement(unsigned index)
{
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    c_vector<double, SPACE_DIM> long_axis = zero_vector<double>(SPACE_DIM);

    // Get Apical and Basal faces
    std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> p_faces = { GetFaceOfType(index, MonolayerVertexElementType::Apical),
                                                                                 GetFaceOfType(index, MonolayerVertexElementType::Basal) };

    // Calculate both faces indepdendently
    std::vector<c_vector<double, SPACE_DIM> > orthogonals;
    std::vector<c_vector<double, 2> > temp_axes;
    bool spherical[2] = { false, false };

    for (unsigned i = 0; i < 2; i++)
    {
        // Calculate normal for the base
        c_vector<double, SPACE_DIM> normal_i;
        CalculateUnitNormalToFaceWithArea(p_faces[i], normal_i);

        // Orthgonals for the base
        c_vector<double, SPACE_DIM> unit = zero_vector<double>(SPACE_DIM);
        unit[0] = 1.0;
        if (inner_prod(normal_i, unit) >= 0.99 || inner_prod(normal_i, unit) <= -0.99)
        {
            // If normal in x-direction, take y-direction
            unit[0] = 0.0;
            unit[1] = 1.0;
        }
        c_vector<double, SPACE_DIM> orthogonal_1 = unit - inner_prod(unit, normal_i) * normal_i;
        orthogonal_1 /= norm_2(orthogonal_1);
        c_vector<double, SPACE_DIM> orthogonal_2 = VectorProduct(normal_i, orthogonal_1);
        orthogonal_2 /= norm_2(orthogonal_2);

        orthogonals.push_back(orthogonal_1);
        orthogonals.push_back(orthogonal_2);

        // Calculate the moments w.r.t. to the basis
        c_vector<double, 3> moments = CalculateMomentsOfFace(p_faces[i], normal_i, orthogonal_1, orthogonal_2);

        // Normalise the moments vector to remove problem of a very small discriminant (see #2874)
        moments /= norm_2(moments);

        // If the principal moments are equal...
        c_vector<double, 2> temp_axis = zero_vector<double>(SPACE_DIM - 1);
        double discriminant = (moments(0) - moments(1)) * (moments(0) - moments(1)) + 4.0 * moments(2) * moments(2);
        if (fabs(discriminant) < DBL_EPSILON)
        {
            // ...then every axis through the centroid is a principal axis, so return a random unit vector
            spherical[i] = true;
            temp_axis(0) = RandomNumberGenerator::Instance()->ranf();
            temp_axis(1) = sqrt(1.0 - temp_axis(0) * temp_axis(0));
        }
        else
        {
            // If the product of inertia is zero, then the coordinate axes are the principal axes
            if (fabs(moments(2)) < DBL_EPSILON)
            {
                if (moments(0) < moments(1))
                {
                    temp_axis(0) = 1.0;
                    temp_axis(1) = 0.0;
                }
                else
                {
                    temp_axis(0) = 0.0;
                    temp_axis(1) = 1.0;
                }
            }
            else
            {
                // Otherwise we find the eigenvector of the inertia matrix corresponding to the smallest eigenvalue
                double lambda = 0.5 * (moments(0) + moments(1) - sqrt(discriminant));

                temp_axis(0) = 1.0;
                temp_axis(1) = (moments(0) - lambda) / moments(2);

                // Normalise the coefficients to get a unit vector
                temp_axis /= norm_2(temp_axis);
            }
        }
        temp_axes.push_back(temp_axis);
    }
    // If the direction is random on both apical an basal side..
    if (spherical[0] && spherical[1])
    {
        // ... choose random direction from apical
        long_axis = temp_axes[0][0] * orthogonals[0] + temp_axes[0][1] * orthogonals[1];
    }
    // if it is random on one side...
    else if (spherical[0] || spherical[1])
    {
        // chose the direction on non-random side
        if (spherical[0])
            long_axis = temp_axes[1][0] * orthogonals[2] + temp_axes[1][1] * orthogonals[3];
        else
            long_axis = temp_axes[0][0] * orthogonals[0] + temp_axes[0][1] * orthogonals[1];
    }
    // if it is not random we take the average
    else
    {
        c_vector<double, SPACE_DIM> long_axis_1 = temp_axes[1][0] * orthogonals[2] + temp_axes[1][1] * orthogonals[3];
        c_vector<double, SPACE_DIM> long_axis_2 = temp_axes[0][0] * orthogonals[0] + temp_axes[0][1] * orthogonals[1];
        // .. but we check that they are not antiparallel.
        if (inner_prod(long_axis_1, long_axis_2) < 0.0)
            long_axis = long_axis_1 - long_axis_2;
        else
            long_axis = long_axis_1 + long_axis_2;

        // std::cout << "1: " << long_axis_1[0] << ", " << long_axis_1[1] << ", " << long_axis_1[2] << "\n";
        // std::cout << "2: " << long_axis_2[0] << ", " << long_axis_2[1] << ", " << long_axis_2[2] << "\n" << std::flush;
    }

    // Normalize to get a unit direction vector
    long_axis /= norm_2(long_axis);
    return long_axis;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::GetBasalLongAxisOfElement(unsigned index)
{
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE - code will be removed at compile time

    c_vector<double, SPACE_DIM> long_axis = zero_vector<double>(SPACE_DIM);

    // Get Basal face
    std::vector<MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> p_faces = { GetFaceOfType(index, MonolayerVertexElementType::Basal) };

    for (unsigned i = 0; i < 1; i++)
    {
        // Calculate normal for the base
        c_vector<double, SPACE_DIM> normal_i;
        CalculateUnitNormalToFaceWithArea(p_faces[i], normal_i);

        // Orthgonals for the base
        c_vector<double, SPACE_DIM> unit = zero_vector<double>(SPACE_DIM);
        unit[0] = 1.0;
        if (inner_prod(normal_i, unit) >= 0.99 || inner_prod(normal_i, unit) <= -0.99)
        {
            // If normal in x-direction, take y-direction
            unit[0] = 0.0;
            unit[1] = 1.0;
        }
        c_vector<double, SPACE_DIM> orthogonal_1 = unit - inner_prod(unit, normal_i) * normal_i;
        orthogonal_1 /= norm_2(orthogonal_1);
        c_vector<double, SPACE_DIM> orthogonal_2 = VectorProduct(normal_i, orthogonal_1);
        orthogonal_2 /= norm_2(orthogonal_2);

        // Calculate the moments w.r.t. to the basis
        c_vector<double, 3> moments = CalculateMomentsOfFace(p_faces[i], normal_i, orthogonal_1, orthogonal_2);

        // Normalise the moments vector to remove problem of a very small discriminant (see #2874)
        moments /= norm_2(moments);

        // If the principal moments are equal...
        c_vector<double, 2> temp_axis = zero_vector<double>(SPACE_DIM - 1);
        double discriminant = (moments(0) - moments(1)) * (moments(0) - moments(1)) + 4.0 * moments(2) * moments(2);
        if (fabs(discriminant) < DBL_EPSILON)
        {
            // ...then every axis through the centroid is a principal axis, so return a random unit vector
            temp_axis(0) = RandomNumberGenerator::Instance()->ranf();
            temp_axis(1) = sqrt(1.0 - temp_axis(0) * temp_axis(0));
        }
        else
        {
            // If the product of inertia is zero, then the coordinate axes are the principal axes
            if (fabs(moments(2)) < DBL_EPSILON)
            {
                if (moments(0) < moments(1))
                {
                    temp_axis(0) = 1.0;
                    temp_axis(1) = 0.0;
                }
                else
                {
                    temp_axis(0) = 0.0;
                    temp_axis(1) = 1.0;
                }
            }
            else
            {
                // Otherwise we find the eigenvector of the inertia matrix corresponding to the smallest eigenvalue
                double lambda = 0.5 * (moments(0) + moments(1) - sqrt(discriminant));

                temp_axis(0) = 1.0;
                temp_axis(1) = (moments(0) - lambda) / moments(2);

                // Normalise the coefficients to get a unit vector
                temp_axis /= norm_2(temp_axis);
            }
        }
        long_axis = temp_axis[0] * orthogonal_1 + temp_axis[1] * orthogonal_2;
    }

    // Normalize to get a unit direction vector
    long_axis /= norm_2(long_axis);
    return long_axis;
}

/// Specialization to avoid compiler error about zero-sized arrays
#if defined(__xlC__)
template <>
double MonolayerVertexMesh<1, 1>::CalculateAreaOfFace(VertexElement<0, 1>* pFace)
{
    NEVER_REACHED;
}
#endif

// Explicit instantiation
template class MonolayerVertexMesh<1, 1>;
template class MonolayerVertexMesh<1, 2>;
template class MonolayerVertexMesh<1, 3>;
template class MonolayerVertexMesh<2, 2>;
template class MonolayerVertexMesh<2, 3>;
template class MonolayerVertexMesh<3, 3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(MonolayerVertexMesh)