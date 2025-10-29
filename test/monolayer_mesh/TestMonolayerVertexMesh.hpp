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

#ifndef TESTVERTEXMESH_HPP_
#define TESTVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "ArchiveOpener.hpp"
#include "MonolayerVertexElement.hpp"
#include "MonolayerVertexMesh.hpp"
#include "MonolayerVertexMeshReader.hpp"
#include "MonolayerVertexMeshWriter.hpp"
#include "MutableMesh.hpp"
// This test is always run sequentially (never in parallel)
#include "UblasCustomFunctions.hpp"
#include "FakePetscSetup.hpp"

class TestMonolayerVertexMesh : public CxxTest::TestSuite
{
private:
    MonolayerVertexMesh<3, 3>* ConstructCubeAndPyramidMesh()
    {
        // Make 8 nodes to assign to a cube and a pyramid element
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(8, false, 0.5, 0.5, 1.5));

        // Set basal/apical nodetypes
        std::vector<MonolayerVertexElementType> node_types;
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);

        // Make six square faces and four triangular faces out of these nodes
        std::vector<std::vector<Node<3>*> > nodes_faces(10);
        std::vector<std::vector<MonolayerVertexElementType> > nodes_faces_types(10);

        unsigned node_indices_face_0[4] = { 0, 2, 4, 1 };
        unsigned node_indices_face_1[4] = { 4, 7, 5, 2 };
        unsigned node_indices_face_2[4] = { 7, 6, 1, 4 };
        unsigned node_indices_face_3[4] = { 0, 3, 5, 2 };
        unsigned node_indices_face_4[4] = { 1, 6, 3, 0 };
        unsigned node_indices_face_5[4] = { 7, 6, 3, 5 };
        unsigned node_indices_face_6[3] = { 6, 7, 8 };
        unsigned node_indices_face_7[3] = { 6, 8, 3 };
        unsigned node_indices_face_8[3] = { 3, 8, 5 };
        unsigned node_indices_face_9[3] = { 5, 8, 7 };

        for (unsigned i = 0; i < 4; i++)
        {
            nodes_faces[0].push_back(nodes[node_indices_face_0[i]]);
            nodes_faces_types[0].push_back(node_types[node_indices_face_0[i]]);

            nodes_faces[1].push_back(nodes[node_indices_face_1[i]]);
            nodes_faces_types[1].push_back(node_types[node_indices_face_1[i]]);

            nodes_faces[2].push_back(nodes[node_indices_face_2[i]]);
            nodes_faces_types[2].push_back(node_types[node_indices_face_2[i]]);

            nodes_faces[3].push_back(nodes[node_indices_face_3[i]]);
            nodes_faces_types[3].push_back(node_types[node_indices_face_3[i]]);

            nodes_faces[4].push_back(nodes[node_indices_face_4[i]]);
            nodes_faces_types[4].push_back(node_types[node_indices_face_4[i]]);

            nodes_faces[5].push_back(nodes[node_indices_face_5[i]]);
            nodes_faces_types[5].push_back(node_types[node_indices_face_5[i]]);

            if (i < 3)
            {
                nodes_faces[6].push_back(nodes[node_indices_face_6[i]]);
                nodes_faces_types[6].push_back(node_types[node_indices_face_6[i]]);

                nodes_faces[7].push_back(nodes[node_indices_face_7[i]]);
                nodes_faces_types[7].push_back(node_types[node_indices_face_7[i]]);

                nodes_faces[8].push_back(nodes[node_indices_face_8[i]]);
                nodes_faces_types[8].push_back(node_types[node_indices_face_8[i]]);

                nodes_faces[9].push_back(nodes[node_indices_face_9[i]]);
                nodes_faces_types[9].push_back(node_types[node_indices_face_9[i]]);
            }
        }

        // Make the faces
        std::vector<MonolayerVertexElement<2, 3>*> faces;

        // Set face types
        std::vector<MonolayerVertexElementType> face_types;
        face_types.push_back(MonolayerVertexElementType::Basal);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Apical);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);

        for (unsigned i = 0; i < 10; i++)
        {
            faces.push_back(new MonolayerVertexElement<2, 3>(i, face_types[i], nodes_faces[i], nodes_faces_types[i]));
        }

        // Make the elements
        std::vector<MonolayerVertexElement<2, 3>*> faces_element_0, faces_element_1;
        std::vector<bool> orientations_0, orientations_1;

        // Cube element
        for (unsigned i = 0; i < 6; i++)
        {
            faces_element_0.push_back(faces[i]);
            orientations_0.push_back(true);
        }

        // Pyramid element
        for (unsigned i = 6; i < 10; i++)
        {
            faces_element_1.push_back(faces[i]);
            orientations_1.push_back(true);
        }
        faces_element_1.push_back(faces[5]);
        orientations_1.push_back(false);

        std::vector<MonolayerVertexElement<3, 3>*> elements;
        elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, faces_element_0, orientations_0));
        elements.push_back(new MonolayerVertexElement<3, 3>(1, MonolayerVertexElementType::Undetermined, faces_element_1, orientations_1));

        return new MonolayerVertexMesh<3, 3>(nodes, faces, elements);
    }

    MonolayerVertexMesh<3, 3>* ConstructPrismMesh()
    {
        // Create nodes
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 1.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 3.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 3.0));

        // Set basal/apical nodetypes
        std::vector<MonolayerVertexElementType> node_types;
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);

        // Make five faces out of these nodes
        std::vector<std::vector<Node<3>*> > nodes_faces(5);
        std::vector<std::vector<MonolayerVertexElementType> > nodes_faces_types(5);

        nodes_faces[0].push_back(nodes[0]);
        nodes_faces[0].push_back(nodes[4]);
        nodes_faces[0].push_back(nodes[5]);
        nodes_faces[0].push_back(nodes[1]);

        nodes_faces_types[0].push_back(node_types[0]);
        nodes_faces_types[0].push_back(node_types[4]);
        nodes_faces_types[0].push_back(node_types[5]);
        nodes_faces_types[0].push_back(node_types[1]);

        nodes_faces[1].push_back(nodes[0]);
        nodes_faces[1].push_back(nodes[3]);
        nodes_faces[1].push_back(nodes[4]);

        nodes_faces_types[1].push_back(node_types[0]);
        nodes_faces_types[1].push_back(node_types[3]);
        nodes_faces_types[1].push_back(node_types[4]);

        nodes_faces[2].push_back(nodes[3]);
        nodes_faces[2].push_back(nodes[2]);
        nodes_faces[2].push_back(nodes[5]);
        nodes_faces[2].push_back(nodes[4]);

        nodes_faces_types[2].push_back(node_types[3]);
        nodes_faces_types[2].push_back(node_types[2]);
        nodes_faces_types[2].push_back(node_types[5]);
        nodes_faces_types[2].push_back(node_types[4]);

        nodes_faces[3].push_back(nodes[1]);
        nodes_faces[3].push_back(nodes[5]);
        nodes_faces[3].push_back(nodes[2]);

        nodes_faces_types[3].push_back(node_types[1]);
        nodes_faces_types[3].push_back(node_types[5]);
        nodes_faces_types[3].push_back(node_types[2]);

        nodes_faces[4].push_back(nodes[3]);
        nodes_faces[4].push_back(nodes[2]);
        nodes_faces[4].push_back(nodes[1]);
        nodes_faces[4].push_back(nodes[0]);

        nodes_faces_types[4].push_back(node_types[3]);
        nodes_faces_types[4].push_back(node_types[2]);
        nodes_faces_types[4].push_back(node_types[1]);
        nodes_faces_types[4].push_back(node_types[0]);

        // Set face types
        std::vector<MonolayerVertexElementType> face_types;
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Basal);

        std::vector<MonolayerVertexElement<2, 3>*> faces(5);
        std::vector<bool> orientations(5);
        for (unsigned i = 0; i < 5; i++)
        {
            faces[i] = new MonolayerVertexElement<2, 3>(i, face_types[i], nodes_faces[i], nodes_faces_types[i]);
            orientations[i] = true;
        }

        // Create cuboidal element
        std::vector<MonolayerVertexElement<3, 3>*> elements;
        elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, faces, orientations));

        // Create mesh
        return new MonolayerVertexMesh<3, 3>(nodes, elements);
    }

    MonolayerVertexMesh<3, 3>* ConstructNonCoplanarCubeMesh()
    {
        // Make 8 nodes to assign to a cube and a pyramid element
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 5.0));

        // Set basal/apical nodetypes
        std::vector<MonolayerVertexElementType> node_types;
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);

        // Make six square faces
        std::vector<std::vector<Node<3>*> > nodes_faces(6);
        std::vector<std::vector<MonolayerVertexElementType> > nodes_faces_types(6);

        unsigned node_indices_face_0[4] = { 0, 2, 4, 1 };
        unsigned node_indices_face_1[4] = { 2, 5, 7, 4 };
        unsigned node_indices_face_2[4] = { 7, 6, 1, 4 };
        unsigned node_indices_face_3[4] = { 0, 3, 5, 2 };
        unsigned node_indices_face_4[4] = { 1, 6, 3, 0 };
        unsigned node_indices_face_5[4] = { 7, 6, 3, 5 };

        // Set face types
        std::vector<MonolayerVertexElementType> face_types;
        face_types.push_back(MonolayerVertexElementType::Basal);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Apical);

        for (unsigned i = 0; i < 4; i++)
        {
            nodes_faces[0].push_back(nodes[node_indices_face_0[i]]);
            nodes_faces_types[0].push_back(node_types[node_indices_face_0[i]]);

            nodes_faces[1].push_back(nodes[node_indices_face_1[i]]);
            nodes_faces_types[1].push_back(node_types[node_indices_face_1[i]]);

            nodes_faces[2].push_back(nodes[node_indices_face_2[i]]);
            nodes_faces_types[2].push_back(node_types[node_indices_face_2[i]]);

            nodes_faces[3].push_back(nodes[node_indices_face_3[i]]);
            nodes_faces_types[3].push_back(node_types[node_indices_face_3[i]]);

            nodes_faces[4].push_back(nodes[node_indices_face_4[i]]);
            nodes_faces_types[4].push_back(node_types[node_indices_face_4[i]]);

            nodes_faces[5].push_back(nodes[node_indices_face_5[i]]);
            nodes_faces_types[5].push_back(node_types[node_indices_face_5[i]]);
        }

        // Make the faces
        std::vector<MonolayerVertexElement<2, 3>*> faces;

        for (unsigned i = 0; i < 6; i++)
        {
            faces.push_back(new MonolayerVertexElement<2, 3>(i, face_types[i], nodes_faces[i], nodes_faces_types[i]));
        }

        // Make the elements
        std::vector<MonolayerVertexElement<2, 3>*> faces_element_0;
        std::vector<bool> orientations_0;

        // Cube element
        for (unsigned i = 0; i < 6; i++)
        {
            faces_element_0.push_back(faces[i]);
            orientations_0.push_back(true);
        }

        std::vector<MonolayerVertexElement<3, 3>*> elements;
        elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, faces_element_0, orientations_0));

        return new MonolayerVertexMesh<3, 3>(nodes, faces, elements);
    }

public:
    void TestBasic1dVertexMesh()
    {
        // Create a 1D mesh comprising four nodes and three elements
        std::vector<Node<1>*> nodes_1d;
        std::vector<MonolayerVertexElementType> node_types_1d;
        for (unsigned i = 0; i < 4; i++)
        {
            nodes_1d.push_back(new Node<1>(i, false, 0.5 * (double)i));
            node_types_1d.push_back(MonolayerVertexElementType::Undetermined);
        }

        std::vector<std::vector<Node<1>*> > nodes_elements_1d(3);
        std::vector<std::vector<MonolayerVertexElementType> > nodes_elements_1d_types(3);
        std::vector<MonolayerVertexElement<1, 1>*> elements_1d;
        for (unsigned i = 0; i < 3; i++)
        {
            nodes_elements_1d[i].push_back(nodes_1d[i]);
            nodes_elements_1d[i].push_back(nodes_1d[i + 1]);
            nodes_elements_1d_types[i].push_back(node_types_1d[i]);
            nodes_elements_1d_types[i].push_back(node_types_1d[i + 1]);

            elements_1d.push_back(new MonolayerVertexElement<1, 1>(i, MonolayerVertexElementType::Undetermined, nodes_elements_1d[i], nodes_elements_1d_types[i]));
        }

        MonolayerVertexMesh<1, 1> mesh_1d(nodes_1d, elements_1d);

        // Test the mesh has the correct number of nodes and elements
        TS_ASSERT_EQUALS(mesh_1d.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(mesh_1d.GetNumElements(), 3u);

        // Test the elements have the correct nodes
        TS_ASSERT_EQUALS(mesh_1d.GetElement(0)->GetNumNodes(), 2u);
        TS_ASSERT_DELTA(mesh_1d.GetElement(0)->GetNodeLocation(0)[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh_1d.GetElement(0)->GetNodeLocation(1)[0], 0.5, 1e-6);

        TS_ASSERT_EQUALS(mesh_1d.GetElement(1)->GetNumNodes(), 2u);
        TS_ASSERT_DELTA(mesh_1d.GetElement(1)->GetNodeLocation(0)[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(mesh_1d.GetElement(1)->GetNodeLocation(1)[0], 1.0, 1e-6);

        TS_ASSERT_EQUALS(mesh_1d.GetElement(2)->GetNumNodes(), 2u);
        TS_ASSERT_DELTA(mesh_1d.GetElement(2)->GetNodeLocation(0)[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(mesh_1d.GetElement(2)->GetNodeLocation(1)[0], 1.5, 1e-6);
    }

    void TestBasic2dVertexMesh()
    {
        // Create a 2D mesh comprising seven nodes and two elements
        std::vector<Node<2>*> nodes_2d;
        nodes_2d.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes_2d.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes_2d.push_back(new Node<2>(2, false, 1.5, 1.0));
        nodes_2d.push_back(new Node<2>(3, false, 1.0, 2.0));
        nodes_2d.push_back(new Node<2>(4, false, 0.0, 1.0));
        nodes_2d.push_back(new Node<2>(5, false, 2.0, 0.0));
        nodes_2d.push_back(new Node<2>(6, false, 2.0, 3.0));

        std::vector<MonolayerVertexElementType> node_types(6, MonolayerVertexElementType::Undetermined);

        std::vector<std::vector<Node<2>*> > nodes_elements_2d(2);
        std::vector<std::vector<MonolayerVertexElementType> > nodes_elements_2d_types(2);
        for (unsigned i = 0; i < 5; i++)
        {
            nodes_elements_2d[0].push_back(nodes_2d[i]);
            nodes_elements_2d_types[0].push_back(node_types[i]);
        }
        nodes_elements_2d[1].push_back(nodes_2d[2]);
        nodes_elements_2d[1].push_back(nodes_2d[5]);
        nodes_elements_2d[1].push_back(nodes_2d[6]);

        nodes_elements_2d_types[1].push_back(node_types[2]);
        nodes_elements_2d_types[1].push_back(node_types[5]);
        nodes_elements_2d_types[1].push_back(node_types[6]);

        std::vector<MonolayerVertexElement<2, 2>*> elements_2d;
        elements_2d.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, nodes_elements_2d[0], nodes_elements_2d_types[0]));
        elements_2d.push_back(new MonolayerVertexElement<2, 2>(1, MonolayerVertexElementType::Undetermined, nodes_elements_2d[1], nodes_elements_2d_types[1]));

        MonolayerVertexMesh<2, 2> mesh_2d(nodes_2d, elements_2d);

        // Test the mesh has the correct numbers of nodes and elements
        TS_ASSERT_EQUALS(mesh_2d.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh_2d.GetNumNodes(), 7u);

        // Further tests of the mesh
        TS_ASSERT_DELTA(mesh_2d.GetNode(2)->rGetLocation()[0], 1.5, 1e-3);
        TS_ASSERT_EQUALS(mesh_2d.GetElement(1)->GetNode(2)->GetIndex(), 6u);

        // Test that the nodes know which elements they are in
        std::set<unsigned> temp_list1;
        temp_list1.insert(0u);

        // Nodes 1 and 4 are only in element 0
        TS_ASSERT_EQUALS(nodes_2d[1]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(nodes_2d[4]->rGetContainingElementIndices(), temp_list1);

        // Node 2 is in elements 0 and 1
        temp_list1.insert(1u);
        TS_ASSERT_EQUALS(nodes_2d[2]->rGetContainingElementIndices(), temp_list1);

        // Node 5 is only in element 1
        std::set<unsigned> temp_list2;
        temp_list2.insert(1u);
        TS_ASSERT_EQUALS(nodes_2d[5]->rGetContainingElementIndices(), temp_list2);

        // Coverage of some methods
        TS_ASSERT_EQUALS(mesh_2d.SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(mesh_2d.SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(mesh_2d.SolveBoundaryElementMapping(0), 0u);
    }

    void TestBasic3dVertexMesh()
    {
        // Create a 3D mesh comprising a cube and pyramid
        MonolayerVertexMesh<3, 3>* p_mesh = ConstructCubeAndPyramidMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 10u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), 2u);

        // Test the location of one of the nodes
        Node<3>* p_node_2 = p_mesh->GetNode(2);
        TS_ASSERT_DELTA(p_node_2->rGetLocation()[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(p_node_2->rGetLocation()[1], 1.0, 1e-3);
        TS_ASSERT_DELTA(p_node_2->rGetLocation()[2], 0.0, 1e-3);

        // Test a couple of the elements
        MonolayerVertexElement<3, 3>* p_element_0 = p_mesh->GetElement(0);
        TS_ASSERT_EQUALS(p_element_0->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(p_element_0->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(p_element_0->GetFaceType(), MonolayerVertexElementType::Undetermined);

        MonolayerVertexElement<3, 3>* p_element_1 = p_mesh->GetElement(1);
        TS_ASSERT_EQUALS(p_element_1->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_element_1->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(p_element_1->GetFaceType(), MonolayerVertexElementType::Undetermined);

        // Check that the nodes know which elements they are in
        std::set<unsigned> temp_list1;
        temp_list1.insert(0);

        // Nodes 0, 1, 2 and 4 are only in element 0
        TS_ASSERT_EQUALS(p_mesh->GetNode(0)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(p_mesh->GetNode(1)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(p_mesh->GetNode(2)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(p_mesh->GetNode(4)->rGetContainingElementIndices(), temp_list1);

        // Node 3, 5, 6 and 7 are in elements 0 and 1
        temp_list1.insert(1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(3)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(p_mesh->GetNode(5)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(p_mesh->GetNode(6)->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(p_mesh->GetNode(7)->rGetContainingElementIndices(), temp_list1);

        // Node 8 is only in element 1
        std::set<unsigned> temp_list2;
        temp_list2.insert(1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(8)->rGetContainingElementIndices(), temp_list2);

        // Coverage of some methods
        TS_ASSERT_EQUALS(p_mesh->SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(p_mesh->SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(p_mesh->SolveBoundaryElementMapping(0), 0u);

        // Avoid memory leaks
        delete p_mesh;
    }

    void TestVertexElementMap()
    {
        VertexElementMap map(10);
        TS_ASSERT_EQUALS(map.Size(), 10u);

        map.ResetToIdentity();
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        map.SetNewIndex(0, 1);
        map.SetNewIndex(1, 0);

        TS_ASSERT_EQUALS(map.GetNewIndex(0), 1u);
        TS_ASSERT_EQUALS(map.GetNewIndex(1), 0u);
        TS_ASSERT_EQUALS(map.GetNewIndex(2), 2u);

        TS_ASSERT_EQUALS(map.IsIdentityMap(), false);

        map.ResetToIdentity();
        map.SetDeleted(4);
        TS_ASSERT_THROWS_THIS(map.GetNewIndex(4), "VertexElement has been deleted");
        TS_ASSERT_EQUALS(map.IsDeleted(4), true);
        TS_ASSERT_EQUALS(map.IsDeleted(5), false);
        TS_ASSERT_EQUALS(map.IsIdentityMap(), false);
    }

    void TestGetCentroidOfElement()
    {
        // Test method with a 1D mesh
        std::vector<Node<1>*> nodes_1d;
        for (unsigned i = 0; i < 4; i++)
        {
            nodes_1d.push_back(new Node<1>(i, false, 0.5 * (double)i));
        }
        std::vector<std::vector<Node<1>*> > nodes_elements_1d(3);
        std::vector<MonolayerVertexElementType> nodes_elements_1d_types(2, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<1, 1>*> elements_1d;
        for (unsigned i = 0; i < 3; i++)
        {
            nodes_elements_1d[i].push_back(nodes_1d[i]);
            nodes_elements_1d[i].push_back(nodes_1d[i + 1]);
            elements_1d.push_back(new MonolayerVertexElement<1, 1>(i, MonolayerVertexElementType::Undetermined, nodes_elements_1d[i], nodes_elements_1d_types));
        }
        MonolayerVertexMesh<1, 1> mesh_1d(nodes_1d, elements_1d);

        TS_ASSERT_DELTA(mesh_1d.GetCentroidOfElement(0)[0], 0.25, 1e-6);
        TS_ASSERT_DELTA(mesh_1d.GetCentroidOfElement(1)[0], 0.75, 1e-6);
        TS_ASSERT_DELTA(mesh_1d.GetCentroidOfElement(2)[0], 1.25, 1e-6);

        // Test method with a single triangular element
        std::vector<Node<2>*> triangle_nodes;
        triangle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        triangle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        std::vector<MonolayerVertexElementType> triangle_node_types(3, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 2>*> triangle_elements;
        triangle_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, triangle_nodes, triangle_node_types));
        MonolayerVertexMesh<2, 2> triangle_mesh(triangle_nodes, triangle_elements);

        c_vector<double, 2> triangle_centroid = triangle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(triangle_centroid[0], 2.0 / 3.0, 1e-4);
        TS_ASSERT_DELTA(triangle_centroid[1], 1.0 / 3.0, 1e-4);

        // Test method with a single square element
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        std::vector<MonolayerVertexElementType> square_node_types(4, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 2>*> square_elements;
        square_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, square_nodes, square_node_types));
        MonolayerVertexMesh<2, 2> square_mesh(square_nodes, square_elements);

        c_vector<double, 2> square_centroid = square_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(square_centroid[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(square_centroid[1], 0.5, 1e-6);

        // Test method with a single rectangular element away from the origin
        std::vector<Node<2>*> far_rectangle_nodes;
        far_rectangle_nodes.push_back(new Node<2>(0, false, 10.0, 10.0));
        far_rectangle_nodes.push_back(new Node<2>(1, false, 11.0, 10.0));
        far_rectangle_nodes.push_back(new Node<2>(2, false, 11.0, 14.0));
        far_rectangle_nodes.push_back(new Node<2>(3, false, 10.0, 14.0));
        std::vector<MonolayerVertexElementType> far_rectangle_node_types(4, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 2>*> far_rectangle_elements;
        far_rectangle_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, far_rectangle_nodes, far_rectangle_node_types));
        MonolayerVertexMesh<2, 2> far_rectangle_mesh(far_rectangle_nodes, far_rectangle_elements);

        c_vector<double, 2> far_rectangle_centroid = far_rectangle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(far_rectangle_centroid[0], 10.5, 1e-4);
        TS_ASSERT_DELTA(far_rectangle_centroid[1], 12.0, 1e-4);

        // Test method with a single rectangular element at a 30 degree angle to the x-axis
        std::vector<Node<2>*> angled_rectangle_nodes;
        angled_rectangle_nodes.push_back(new Node<2>(0, false, 2.0 * 0.5 * sqrt(3.0) - 1.0 * 0.5, 2.0 * 0.5 + 1.0 * 0.5 * sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(1, false, -2.0 * 0.5 * sqrt(3.0) - 1.0 * 0.5, -2.0 * 0.5 + 1.0 * 0.5 * sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(2, false, -2.0 * 0.5 * sqrt(3.0) + 1.0 * 0.5, -2.0 * 0.5 - 1.0 * 0.5 * sqrt(3.0)));
        angled_rectangle_nodes.push_back(new Node<2>(3, false, 2.0 * 0.5 * sqrt(3.0) + 1.0 * 0.5, 2.0 * 0.5 - 1.0 * 0.5 * sqrt(3.0)));
        std::vector<MonolayerVertexElementType> angled_rectangle_node_types(4, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 2>*> angled_rectangle_elements;
        angled_rectangle_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, angled_rectangle_nodes, angled_rectangle_node_types));
        MonolayerVertexMesh<2, 2> angled_rectangle_mesh(angled_rectangle_nodes, angled_rectangle_elements);

        c_vector<double, 2> angled_rectangle_centroid = angled_rectangle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(angled_rectangle_centroid(0), 0.0, 1e-6);
        TS_ASSERT_DELTA(angled_rectangle_centroid(1), 0.0, 1e-6);

        // Test method with a single element that is close to circular
        std::vector<Node<2>*> circle_nodes;
        std::vector<MonolayerVertexElementType> circle_node_types;
        unsigned num_nodes = 1000;
        for (unsigned i = 0; i < num_nodes; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / (double)(num_nodes);
            circle_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
            circle_node_types.push_back(MonolayerVertexElementType::Undetermined);
        }
        std::vector<MonolayerVertexElement<2, 2>*> circle_elements;
        circle_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, circle_nodes, circle_node_types));
        MonolayerVertexMesh<2, 2> circle_mesh(circle_nodes, circle_elements);

        c_vector<double, 2> circle_centroid = circle_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(circle_centroid[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(circle_centroid[1], 0.0, 1e-4);

        // Test method with a single hexagonal element centred at the origin
        std::vector<Node<2>*> hexagon_nodes;
        std::vector<MonolayerVertexElementType> hexagon_node_types;
        for (unsigned i = 0; i < 6; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / 6.0;
            hexagon_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
            hexagon_node_types.push_back(MonolayerVertexElementType::Undetermined);
        }
        std::vector<MonolayerVertexElement<2, 2>*> hexagon_elements;
        hexagon_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, hexagon_nodes, hexagon_node_types));
        MonolayerVertexMesh<2, 2> hexagon_mesh(hexagon_nodes, hexagon_elements);

        c_vector<double, 2> hexagon_centroid = hexagon_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(hexagon_centroid[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(hexagon_centroid[1], 0.0, 1e-6);

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<Node<2>*> irregular_nodes;
        irregular_nodes.push_back(new Node<2>(0, false, 2.7663, 0.1033));
        irregular_nodes.push_back(new Node<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new Node<2>(2, false, 1.9433, 0.5313));
        irregular_nodes.push_back(new Node<2>(3, false, 2.7425, 1.9461));
        irregular_nodes.push_back(new Node<2>(4, false, 3.3158, 1.5588));
        std::vector<MonolayerVertexElementType> irregular_node_types(5, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 2>*> irregular_elements;
        irregular_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, irregular_nodes, irregular_node_types));
        MonolayerVertexMesh<2, 2> irregular_mesh(irregular_nodes, irregular_elements);

        c_vector<double, 2> irregular_centroid = irregular_mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(irregular_centroid[0], 2.6269, 1e-4);
        TS_ASSERT_DELTA(irregular_centroid[1], 0.8930, 1e-4);
        /*
        // Test method with a regular mesh with hexagonal elements of edge length 1/sqrt(3.0)
        MonolayerVertexMeshReader<2, 2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        MonolayerVertexMesh<2, 2> regular_mesh;
        regular_mesh.ConstructFromMeshReader(mesh_reader);

        c_vector<double, 2> regular_centroid_5 = regular_mesh.GetCentroidOfElement(5);
        TS_ASSERT_DELTA(regular_centroid_5[0], 2.0, 1e-4);
        TS_ASSERT_DELTA(regular_centroid_5[1], 2.5 / sqrt(3.0), 1e-4);

        c_vector<double, 2> regular_centroid_7 = regular_mesh.GetCentroidOfElement(7);
        TS_ASSERT_DELTA(regular_centroid_7[0], 4.0, 1e-4);
        TS_ASSERT_DELTA(regular_centroid_7[1], 2.5 / sqrt(3.0), 1e-4);
        */
        // Test method with a 3D mesh
        MonolayerVertexMesh<3, 3>* p_mesh = ConstructPrismMesh();

        // By symmetry, the centroid of the prism should lie in the plane x=0.5
        c_vector<double, 3> centroid = p_mesh->GetCentroidOfElement(0);
        TS_ASSERT_DELTA(centroid(0), 0.5, 1e-5);

        // Avoid memory leaks
        delete p_mesh;
    }

    void TestGetPassiveCenterOfFace()
    {
        // Test method with a single triangular element
        std::vector<Node<3>*> triangle_nodes;
        triangle_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        triangle_nodes.push_back(new Node<3>(1, false, 2.0, 0.0, 0.0));
        triangle_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 3.0));
        std::vector<MonolayerVertexElementType> triangle_node_types(3, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 3>*> triangle_faces;
        triangle_faces.push_back(new MonolayerVertexElement<2, 3>(0, MonolayerVertexElementType::Undetermined, triangle_nodes, triangle_node_types));

        std::vector<bool> triangle_faces_orientation(1, true);
        std::vector<MonolayerVertexElement<3, 3>*> triangle_3d_element;
        triangle_3d_element.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, triangle_faces, triangle_faces_orientation));

        MonolayerVertexMesh<3, 3> triangle_mesh(triangle_nodes, triangle_3d_element);

        MonolayerVertexElement<2, 3>* triangle_face_to_check = triangle_mesh.GetElement(0)->GetFace(0);
        c_vector<double, 3> triangle_passive_center = triangle_mesh.GetPassiveCenterOfFace(triangle_face_to_check);
        TS_ASSERT_DELTA(triangle_passive_center[0], 2.0 / 3.0, 1e-4);
        TS_ASSERT_DELTA(triangle_passive_center[1], 1.0 / 3.0, 1e-4);
        TS_ASSERT_DELTA(triangle_passive_center[2], 1.0, 1e-4);

        // Test method with a single square element
        std::vector<Node<3>*> square_nodes;
        square_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        square_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        square_nodes.push_back(new Node<3>(2, false, 1.0, 1.0, 0.0));
        square_nodes.push_back(new Node<3>(3, false, 0.0, 1.0, 0.0));
        std::vector<MonolayerVertexElementType> square_node_types(4, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 3>*> square_faces;
        square_faces.push_back(new MonolayerVertexElement<2, 3>(0, MonolayerVertexElementType::Undetermined, square_nodes, square_node_types));

        std::vector<bool> square_faces_orientation(1, true);
        std::vector<MonolayerVertexElement<3, 3>*> square_3d_element;
        square_3d_element.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, square_faces, square_faces_orientation));

        MonolayerVertexMesh<3, 3> square_mesh(square_nodes, square_3d_element);

        MonolayerVertexElement<2, 3>* square_face_to_check = square_mesh.GetElement(0)->GetFace(0);
        c_vector<double, 3> square_passive_center = square_mesh.GetPassiveCenterOfFace(square_face_to_check);
        TS_ASSERT_DELTA(square_passive_center[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(square_passive_center[1], 0.5, 1e-6);
        TS_ASSERT_DELTA(square_passive_center[2], 0.0, 1e-6);

        // Test method with a single rectangular element away from the origin
        std::vector<Node<3>*> far_rectangle_nodes;
        far_rectangle_nodes.push_back(new Node<3>(0, false, 10.0, 10.0, 7.0));
        far_rectangle_nodes.push_back(new Node<3>(1, false, 11.0, 10.0, 7.0));
        far_rectangle_nodes.push_back(new Node<3>(2, false, 11.0, 14.0, 7.0));
        far_rectangle_nodes.push_back(new Node<3>(3, false, 10.0, 14.0, 7.0));
        std::vector<MonolayerVertexElementType> far_rectangle_node_types(4, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 3>*> far_rectangle_faces;
        far_rectangle_faces.push_back(new MonolayerVertexElement<2, 3>(0, MonolayerVertexElementType::Undetermined, far_rectangle_nodes, far_rectangle_node_types));
        std::vector<bool> far_rectangle_faces_orientation(1, true);
        std::vector<MonolayerVertexElement<3, 3>*> far_rectangle_3d_element;
        far_rectangle_3d_element.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, far_rectangle_faces, far_rectangle_faces_orientation));

        MonolayerVertexMesh<3, 3> far_rectangle_mesh(far_rectangle_nodes, far_rectangle_3d_element);

        MonolayerVertexElement<2, 3>* far_rectangle_face_to_check = far_rectangle_mesh.GetElement(0)->GetFace(0);
        c_vector<double, 3> far_rectangle_passive_center = far_rectangle_mesh.GetPassiveCenterOfFace(far_rectangle_face_to_check);
        TS_ASSERT_DELTA(far_rectangle_passive_center[0], 10.5, 1e-4);
        TS_ASSERT_DELTA(far_rectangle_passive_center[1], 12.0, 1e-4);
        TS_ASSERT_DELTA(far_rectangle_passive_center[2], 7.0, 1e-4);

        // Test method with a 3D cube and pyramid
        MonolayerVertexMesh<3, 3>* p_mesh = ConstructCubeAndPyramidMesh();

        MonolayerVertexElement<3, 3>* p_cube = p_mesh->GetElement(0);
        MonolayerVertexElement<3, 3>* p_pyramid = p_mesh->GetElement(1);

        // Test cube faces
        MonolayerVertexElement<2, 3>* p_cube_face0 = p_cube->GetFace(0);
        MonolayerVertexElement<2, 3>* p_cube_face2 = p_cube->GetFace(2);

        c_vector<double, 3> cube_face0_passive_center = p_mesh->GetPassiveCenterOfFace(p_cube_face0);
        TS_ASSERT_DELTA(cube_face0_passive_center[0], 0.5, 1e-4);
        TS_ASSERT_DELTA(cube_face0_passive_center[1], 0.5, 1e-4);
        TS_ASSERT_DELTA(cube_face0_passive_center[2], 0.0, 1e-4);

        c_vector<double, 3> cube_face2_passive_center = p_mesh->GetPassiveCenterOfFace(p_cube_face2);
        TS_ASSERT_DELTA(cube_face2_passive_center[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(cube_face2_passive_center[1], 0.5, 1e-4);
        TS_ASSERT_DELTA(cube_face2_passive_center[2], 0.5, 1e-4);

        // Test pyramid faces
        MonolayerVertexElement<2, 3>* p_pyramid_face1 = p_pyramid->GetFace(1);
        MonolayerVertexElement<2, 3>* p_pyramid_face2 = p_pyramid->GetFace(2);

        c_vector<double, 3> pyramid_face1_passive_center = p_mesh->GetPassiveCenterOfFace(p_pyramid_face1);
        TS_ASSERT_DELTA(pyramid_face1_passive_center[0], 1.5 / 3.0, 1e-4);
        TS_ASSERT_DELTA(pyramid_face1_passive_center[1], 0.5 / 3.0, 1e-4);
        TS_ASSERT_DELTA(pyramid_face1_passive_center[2], 3.5 / 3.0, 1e-4);

        c_vector<double, 3> pyramid_face2_passive_center = p_mesh->GetPassiveCenterOfFace(p_pyramid_face2);
        TS_ASSERT_DELTA(pyramid_face2_passive_center[0], 0.5 / 3.0, 1e-4);
        TS_ASSERT_DELTA(pyramid_face2_passive_center[1], 1.5 / 3.0, 1e-4);
        TS_ASSERT_DELTA(pyramid_face2_passive_center[2], 3.5 / 3.0, 1e-4);

        // Test GetPassiveCenterOfFaceTypeInElement
        c_vector<double, 3> cube_face_apical_passive_center = p_mesh->GetPassiveCenterOfFaceTypeInElement(0, MonolayerVertexElementType::Apical);
        c_vector<double, 3> cube_face_basal_passive_center = p_mesh->GetPassiveCenterOfFaceTypeInElement(0, MonolayerVertexElementType::Basal);

        TS_ASSERT_DELTA(cube_face_basal_passive_center[0], 0.5, 1e-4);
        TS_ASSERT_DELTA(cube_face_basal_passive_center[1], 0.5, 1e-4);
        TS_ASSERT_DELTA(cube_face_basal_passive_center[2], 0.0, 1e-4);

        TS_ASSERT_DELTA(cube_face_apical_passive_center[0], 0.5, 1e-4);
        TS_ASSERT_DELTA(cube_face_apical_passive_center[1], 0.5, 1e-4);
        TS_ASSERT_DELTA(cube_face_apical_passive_center[2], 1.0, 1e-4);

        // Test GetThicknessOfElement
        double thickness_cube = p_mesh->GetThicknessOfElement(0);
        TS_ASSERT_DELTA(thickness_cube, 1.0, 1e-4);

        // Avoid memory leaks
        delete p_mesh;
    }

    void TestGetVolumeOfElement()
    {
        // Note that the method GetVolumeOfElement() only works in 2D and 3D

        // Test method with a single triangular element
        std::vector<Node<2>*> triangle_nodes;
        triangle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        triangle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        std::vector<MonolayerVertexElementType> triangle_node_types(3, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 2>*> triangle_elements;
        triangle_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, triangle_nodes, triangle_node_types));
        MonolayerVertexMesh<2, 2> triangle_mesh(triangle_nodes, triangle_elements);

        TS_ASSERT_DELTA(triangle_mesh.GetVolumeOfElement(0), 1.0, 1e-4);

        // Test method with a single square element
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        std::vector<MonolayerVertexElementType> square_node_types(4, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 2>*> square_elements;
        square_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, square_nodes, square_node_types));
        MonolayerVertexMesh<2, 2> square_mesh(square_nodes, square_elements);

        TS_ASSERT_DELTA(square_mesh.GetVolumeOfElement(0), 1.0, 1e-6);

        // Test method with a single element that is close to circular
        std::vector<Node<2>*> circle_nodes;
        std::vector<MonolayerVertexElementType> circle_node_types;
        unsigned num_nodes = 1000;
        for (unsigned i = 0; i < num_nodes; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / (double)(num_nodes);
            circle_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
            circle_node_types.push_back(MonolayerVertexElementType::Undetermined);
        }
        std::vector<MonolayerVertexElement<2, 2>*> circle_elements;
        circle_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, circle_nodes, circle_node_types));
        MonolayerVertexMesh<2, 2> circle_mesh(circle_nodes, circle_elements);

        TS_ASSERT_DELTA(circle_mesh.GetVolumeOfElement(0), M_PI, 1e-4);

        // Test method with a single hexagonal element centred at the origin
        std::vector<Node<2>*> hexagon_nodes;
        std::vector<MonolayerVertexElementType> hexagon_node_types;
        for (unsigned i = 0; i < 6; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / 6.0;
            hexagon_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
            hexagon_node_types.push_back(MonolayerVertexElementType::Undetermined);
        }
        std::vector<MonolayerVertexElement<2, 2>*> hexagon_elements;
        hexagon_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, hexagon_nodes, hexagon_node_types));
        MonolayerVertexMesh<2, 2> hexagon_mesh(hexagon_nodes, hexagon_elements);

        TS_ASSERT_DELTA(hexagon_mesh.GetVolumeOfElement(0), 1.5 * sqrt(3.0), 1e-6);

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<Node<2>*> irregular_nodes;
        irregular_nodes.push_back(new Node<2>(0, false, 2.7663, 0.1033));
        irregular_nodes.push_back(new Node<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new Node<2>(2, false, 1.9433, 0.5313));
        irregular_nodes.push_back(new Node<2>(3, false, 2.7425, 1.9461));
        irregular_nodes.push_back(new Node<2>(4, false, 3.3158, 1.5588));
        std::vector<MonolayerVertexElementType> irregular_node_types(5, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 2>*> irregular_elements;
        irregular_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, irregular_nodes, irregular_node_types));
        MonolayerVertexMesh<2, 2> irregular_mesh(irregular_nodes, irregular_elements);

        TS_ASSERT_DELTA(irregular_mesh.GetVolumeOfElement(0), 1.4684, 1e-3);
        /*
        // Test method with a regular mesh with hexagonal elements of edge length 1/sqrt(3.0)
        MonolayerVertexMeshReader<2, 2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        MonolayerVertexMesh<2, 2> regular_mesh;
        regular_mesh.ConstructFromMeshReader(mesh_reader);

        for (MonolayerVertexMesh<2, 2>::MonolayerVertexElementIterator iter = regular_mesh.GetElementIteratorBegin();
             iter != regular_mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned elem_index = iter->GetIndex();
            TS_ASSERT_DELTA(regular_mesh.GetVolumeOfElement(elem_index), 0.5 * sqrt(3.0), 1e-4);
        }
        */
        // Test method with a 3D mesh
        MonolayerVertexMesh<3, 3>* p_mesh = ConstructPrismMesh();

        // The volume of the prism should be 0.5 * 3 * 2 * 1 = 3
        TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(0), 3.0, 1e-6);

        // Avoid memory leaks
        delete p_mesh;
    }

    void TestGetSurfaceAreaOfElement()
    {
        // Note that the method GetSurfaceAreaOfElement() only works in 2D and 3D

        // Test method with a single triangular element
        std::vector<Node<2>*> triangle_nodes;
        triangle_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        triangle_nodes.push_back(new Node<2>(1, false, 2.0, 0.0));
        triangle_nodes.push_back(new Node<2>(2, false, 0.0, 1.0));
        std::vector<MonolayerVertexElementType> triangle_node_types(3, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 2>*> triangle_elements;
        triangle_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, triangle_nodes, triangle_node_types));
        MonolayerVertexMesh<2, 2> triangle_mesh(triangle_nodes, triangle_elements);

        TS_ASSERT_DELTA(triangle_mesh.GetSurfaceAreaOfElement(0), 3.0 + sqrt(5.0), 1e-4);

        // Test method with a single square element
        std::vector<Node<2>*> square_nodes;
        square_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        square_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        square_nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        square_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        std::vector<MonolayerVertexElementType> square_node_types(4, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 2>*> square_elements;
        square_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, square_nodes, square_node_types));
        MonolayerVertexMesh<2, 2> square_mesh(square_nodes, square_elements);

        TS_ASSERT_DELTA(square_mesh.GetSurfaceAreaOfElement(0), 4.0, 1e-6);

        // Test method with a single element that is close to circular
        std::vector<Node<2>*> circle_nodes;
        std::vector<MonolayerVertexElementType> circle_node_types;
        unsigned num_nodes = 1000;
        for (unsigned i = 0; i < num_nodes; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / (double)(num_nodes);
            circle_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
            circle_node_types.push_back(MonolayerVertexElementType::Undetermined);
        }
        std::vector<MonolayerVertexElement<2, 2>*> circle_elements;
        circle_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, circle_nodes, circle_node_types));
        MonolayerVertexMesh<2, 2> circle_mesh(circle_nodes, circle_elements);

        TS_ASSERT_DELTA(circle_mesh.GetSurfaceAreaOfElement(0), 2.0 * M_PI, 1e-4);

        // Test method with a single hexagonal element centred at the origin
        std::vector<Node<2>*> hexagon_nodes;
        std::vector<MonolayerVertexElementType> hexagon_node_types;
        for (unsigned i = 0; i < 6; i++)
        {
            double theta = 2.0 * M_PI * (double)(i) / 6.0;
            hexagon_nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
            hexagon_node_types.push_back(MonolayerVertexElementType::Undetermined);
        }
        std::vector<MonolayerVertexElement<2, 2>*> hexagon_elements;
        hexagon_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, hexagon_nodes, hexagon_node_types));
        MonolayerVertexMesh<2, 2> hexagon_mesh(hexagon_nodes, hexagon_elements);

        TS_ASSERT_DELTA(hexagon_mesh.GetSurfaceAreaOfElement(0), 6.0, 1e-6);

        // Test method with a single irregular element (compare with pen-and-paper solution)
        std::vector<Node<2>*> irregular_nodes;
        irregular_nodes.push_back(new Node<2>(0, false, 2.7663, 0.1033));
        irregular_nodes.push_back(new Node<2>(1, false, 2.2766, -0.1076));
        irregular_nodes.push_back(new Node<2>(2, false, 1.9433, 0.5313));
        irregular_nodes.push_back(new Node<2>(3, false, 2.7425, 1.9461));
        irregular_nodes.push_back(new Node<2>(4, false, 3.3158, 1.5588));
        std::vector<MonolayerVertexElementType> irregular_node_types(5, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 2>*> irregular_elements;
        irregular_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, irregular_nodes, irregular_node_types));
        MonolayerVertexMesh<2, 2> irregular_mesh(irregular_nodes, irregular_elements);

        TS_ASSERT_DELTA(irregular_mesh.GetSurfaceAreaOfElement(0), 5.1263, 1e-3);
        /*
        // Test method with a regular mesh with hexagonal elements of edge length 1/sqrt(3.0)
        MonolayerVertexMeshReader<2, 2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_4_by_4");
        MonolayerVertexMesh<2, 2> regular_mesh;
        regular_mesh.ConstructFromMeshReader(mesh_reader);

        for (MonolayerVertexMesh<2, 2>::MonolayerVertexElementIterator iter = regular_mesh.GetElementIteratorBegin();
             iter != regular_mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned elem_index = iter->GetIndex();
            TS_ASSERT_DELTA(regular_mesh.GetSurfaceAreaOfElement(elem_index), 2 * sqrt(3.0), 1e-4);
        }
        */
        // Test method with a 3D mesh
        MonolayerVertexMesh<3, 3>* p_mesh = ConstructPrismMesh();

        // The surface area of the prism should be the sum of the face areas
        TS_ASSERT_DELTA(p_mesh->GetSurfaceAreaOfElement(0), 11 + sqrt(13.0), 1e-6);

        // Avoid memory leaks
        delete p_mesh;
    }

    void Test3dMethodsWithPrism()
    {
        // Test method with a 3D mesh
        MonolayerVertexMesh<3, 3>* p_mesh = ConstructPrismMesh();

        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 5u);

        // Face 0 has four vertices, is perpendicular to the y axis, and has area 1*3 = 3
        MonolayerVertexElement<2, 3>* p_face_0 = p_mesh->GetFace(0);
        TS_ASSERT_EQUALS(p_face_0->GetNumNodes(), 4u);
        c_vector<double, 3> unit_normal_0;
        p_mesh->CalculateUnitNormalToFaceWithArea(p_face_0, unit_normal_0);
        TS_ASSERT_DELTA(unit_normal_0[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_0[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_0[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_face_0), 3.0, 1e-6);

        // Face 1 has three vertices, is perpendicular to the x axis, and has area 0.5*2*3 = 3
        MonolayerVertexElement<2, 3>* p_face_1 = p_mesh->GetFace(1);
        TS_ASSERT_EQUALS(p_face_1->GetNumNodes(), 3u);
        c_vector<double, 3> unit_normal_1;
        p_mesh->CalculateUnitNormalToFaceWithArea(p_face_1, unit_normal_1);
        TS_ASSERT_DELTA(unit_normal_1[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_1[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_1[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_face_1), 3.0, 1e-6);

        // Face 2 has four vertices, is at an angle theta to the y axis where tan(theta) = 2/3,
        // and has area 1*sqrt(2^2 + 3^2) = sqrt(13.0)
        MonolayerVertexElement<2, 3>* p_face_2 = p_mesh->GetFace(2);
        TS_ASSERT_EQUALS(p_face_2->GetNumNodes(), 4u);
        c_vector<double, 3> unit_normal_2;
        p_mesh->CalculateUnitNormalToFaceWithArea(p_face_2, unit_normal_2);
        TS_ASSERT_DELTA(unit_normal_2[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_2[1], -sin(atan2(3.0, 2.0)), 1e-6);
        TS_ASSERT_DELTA(unit_normal_2[2], -cos(atan2(3.0, 2.0)), 1e-6);
        TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_face_2), sqrt(13.0), 1e-6);

        // Face 1 has three vertices, is perpendicular to the x axis, and has area 0.5*2*3 = 3
        MonolayerVertexElement<2, 3>* p_face_3 = p_mesh->GetFace(3);
        TS_ASSERT_EQUALS(p_face_3->GetNumNodes(), 3u);
        c_vector<double, 3> unit_normal_3;
        p_mesh->CalculateUnitNormalToFaceWithArea(p_face_3, unit_normal_3);
        TS_ASSERT_DELTA(unit_normal_3[0], -1.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_3[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_3[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_face_3), 3.0, 1e-6);

        // Face 4 has four vertices, is perpendicular to the z axis, and has area 1*2 = 2
        MonolayerVertexElement<2, 3>* p_face_4 = p_mesh->GetFace(4);
        TS_ASSERT_EQUALS(p_face_4->GetNumNodes(), 4u);
        c_vector<double, 3> unit_normal_4;
        p_mesh->CalculateUnitNormalToFaceWithArea(p_face_4, unit_normal_4);
        TS_ASSERT_DELTA(unit_normal_4[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_4[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_4[2], -1.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_face_4), 2.0, 1e-6);

        // Avoid memory leaks
        delete p_mesh;
    }

    void TestNonCoplanarTriangulation()
    {
        MonolayerVertexMesh<3, 3>* p_mesh = ConstructNonCoplanarCubeMesh();

        TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 6u);

        // Face 5 has four vertices, is non-coplanar with passive center at (0.5,0.5,2)
        MonolayerVertexElement<2, 3>* p_face_5 = p_mesh->GetFace(5);
        TS_ASSERT_EQUALS(p_face_5->GetNumNodes(), 4u);

        c_vector<double, 3> passive_center_5;
        passive_center_5 = p_mesh->GetPassiveCenterOfFace(p_face_5);

        TS_ASSERT_DELTA(passive_center_5[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(passive_center_5[1], 0.5, 1e-6);
        TS_ASSERT_DELTA(passive_center_5[2], 2.0, 1e-6);

        // Paper-calculated triangle normals
        c_vector<double, 3> triangle_normal_0;
        c_vector<double, 3> triangle_normal_1;
        c_vector<double, 3> triangle_normal_2;
        c_vector<double, 3> triangle_normal_3;

        triangle_normal_0 <<= 0.0, -2.0 / sqrt(5.0), 1.0 / sqrt(5.0);
        triangle_normal_1 <<= -2.0 / sqrt(5.0), 0.0, 1.0 / sqrt(5.0);
        triangle_normal_2 <<= -4.0 / sqrt(21.0), -2.0 / sqrt(21.0), 1.0 / sqrt(21.0);
        triangle_normal_3 <<= -2.0 / sqrt(21.0), -4.0 / sqrt(21.0), 1.0 / sqrt(21.0);

        std::vector<c_vector<double, 3> > triangle_normals;
        triangle_normals.push_back(triangle_normal_0);
        triangle_normals.push_back(triangle_normal_1);
        triangle_normals.push_back(triangle_normal_2);
        triangle_normals.push_back(triangle_normal_3);

        std::vector<double> triangle_areas;
        triangle_areas.push_back(sqrt(5.0) / 4.0);
        triangle_areas.push_back(sqrt(5.0) / 4.0);
        triangle_areas.push_back(sqrt(21.0) / 4.0);
        triangle_areas.push_back(sqrt(21.0) / 4.0);

        // Check length of paper-calculated triangle normals
        TS_ASSERT_DELTA(norm_2(triangle_normal_0), 1.0, 1e-6);
        TS_ASSERT_DELTA(norm_2(triangle_normal_1), 1.0, 1e-6);
        TS_ASSERT_DELTA(norm_2(triangle_normal_2), 1.0, 1e-6);
        TS_ASSERT_DELTA(norm_2(triangle_normal_3), 1.0, 1e-6);

        // Calculate normal vectors and areas of triangulation triangles
        std::vector<c_vector<double, 3> > calculated_triangle_normals;
        std::vector<double> calculated_triangle_areas = p_mesh->CalculateAreasAndNormalsOfTriangulationsOfFace(p_face_5, &calculated_triangle_normals);

        for (unsigned triangle_index = 0; triangle_index < triangle_normals.size(); triangle_index++)
        {
            bool found_orthogonal_calculated_normal = false;
            for (unsigned calculated_index = 0; calculated_index < calculated_triangle_normals.size(); calculated_index++)
            {
                c_vector<double, 3> cross_product = VectorProduct(triangle_normals[triangle_index], calculated_triangle_normals[calculated_index]);
                double cross_product_norm = norm_2(cross_product);
                if (cross_product_norm <= 1e-6 && cross_product_norm >= -1e-6)
                {
                    double inner_product = 0.0;
                    for (unsigned i = 0; i < 3; i++)
                    {
                        inner_product += triangle_normals[triangle_index][i] * calculated_triangle_normals[calculated_index][i];
                    }
                    // Check if orthogonal normals are projected to 1 on calculated unit normals
                    TS_ASSERT_DELTA(fabs(inner_product), 1.0, 1e-6);
                    found_orthogonal_calculated_normal = true;
                    // Check if calculated area is identical to paper area
                    TS_ASSERT_DELTA(calculated_triangle_areas[calculated_index], triangle_areas[triangle_index], 1e-6)
                    break;
                }
            }
            // Check if all paper-normals exist in calculated normals
            TS_ASSERT(found_orthogonal_calculated_normal);
        }

        // Check if calculated face area is correct
        double correct_face_area_5 = (sqrt(21.0) + sqrt(5.0)) / 2.0;
        TS_ASSERT_DELTA(correct_face_area_5, p_mesh->CalculateAreaOfFace(p_face_5), 1e-6);

        // Check if volume of element is correct
        //  Volume of cube is 1.0, volume of four edges to upper passive center is 1/3=(1/4*1/3)*4
        //  Volume of rest pyramids is (4*1*0.5*0.5/3) * 2pyramids
        TS_ASSERT_DELTA(2.0, p_mesh->GetVolumeOfElement(0), 1e-6)

        // Face 2 has four vertices, is coplanar with passive center at (1.0, 0.5, 1.5)
        MonolayerVertexElement<2, 3>* p_face_2 = p_mesh->GetFace(2);
        TS_ASSERT_EQUALS(p_face_2->GetNumNodes(), 4u);

        c_vector<double, 3> passive_center_2;
        passive_center_2 = p_mesh->GetPassiveCenterOfFace(p_face_2);

        TS_ASSERT_DELTA(passive_center_2[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(passive_center_2[1], 0.5, 1e-6);
        TS_ASSERT_DELTA(passive_center_2[2], 1.5, 1e-6);

        // Calculate normal vectors and areas of triangulation triangles
        std::vector<c_vector<double, 3> > calculated_triangle_normals_2;
        std::vector<double> calculated_triangle_areas_2 = p_mesh->CalculateAreasAndNormalsOfTriangulationsOfFace(p_face_2, &calculated_triangle_normals_2);

        for (unsigned triangle_index = 1; triangle_index < calculated_triangle_normals_2.size(); triangle_index++)
        {
            bool correct_size = false;
            if (calculated_triangle_areas_2[triangle_index] - 0.75 >= -1e-6 && calculated_triangle_areas_2[triangle_index] - 0.75 <= 1e-6)
            {
                correct_size = true;
            }
            else if (calculated_triangle_areas_2[triangle_index] - 0.25 >= -1e-6 && calculated_triangle_areas_2[triangle_index] - 0.25 <= 1e-6)
            {
                correct_size = true;
            }
            else if (calculated_triangle_areas_2[triangle_index] - 1.25 >= -1e-6 && calculated_triangle_areas_2[triangle_index] - 1.25 <= 1e-6)
            {
                correct_size = true;
            }
            // Check if triangle has one of the correct sizes
            TS_ASSERT(correct_size);
            // Check if triangle has correct unit normal
            TS_ASSERT_DELTA(calculated_triangle_normals_2[triangle_index][0], 1.0, 1e-6);
            TS_ASSERT_DELTA(calculated_triangle_normals_2[triangle_index][1], 0.0, 1e-6);
            TS_ASSERT_DELTA(calculated_triangle_normals_2[triangle_index][2], 0.0, 1e-6);
        }
        // Check if total area of face is correct
        TS_ASSERT_DELTA(3.0, p_mesh->CalculateAreaOfFace(p_face_2), 1e-6);

        // Face 1 has four vertices, is coplanar with passive center at (0.5, 1.0, 1.5)
        MonolayerVertexElement<2, 3>* p_face_1 = p_mesh->GetFace(1);
        TS_ASSERT_EQUALS(p_face_1->GetNumNodes(), 4u);

        c_vector<double, 3> passive_center_1;
        passive_center_1 = p_mesh->GetPassiveCenterOfFace(p_face_1);
        TS_ASSERT_DELTA(passive_center_1[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(passive_center_1[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(passive_center_1[2], 1.5, 1e-6);
        TS_ASSERT_DELTA(3.0, p_mesh->CalculateAreaOfFace(p_face_1), 1e-6);
        c_vector<double, 3> unit_normal_1;
        p_mesh->CalculateUnitNormalToFaceWithArea(p_face_1, unit_normal_1);
        TS_ASSERT_DELTA(unit_normal_1[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_1[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_1[2], 0.0, 1e-6);

        // Check Passive Centers, normals and areas of face 0
        MonolayerVertexElement<2, 3>* p_face_0 = p_mesh->GetFace(0);
        c_vector<double, 3> passive_center_0 = p_mesh->GetPassiveCenterOfFace(p_face_0);
        TS_ASSERT_DELTA(passive_center_0[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(passive_center_0[1], 0.5, 1e-6);
        TS_ASSERT_DELTA(passive_center_0[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(1.0, p_mesh->CalculateAreaOfFace(p_face_0), 1e-6);
        std::vector<c_vector<double, 3> > calculated_triangle_normals_0;
        std::vector<double> calculated_triangle_areas_0 = p_mesh->CalculateAreasAndNormalsOfTriangulationsOfFace(p_face_0, &calculated_triangle_normals_0);
        for (unsigned triangle_index = 0; triangle_index < calculated_triangle_normals_0.size(); triangle_index++)
        {
            TS_ASSERT_DELTA(0.25, calculated_triangle_areas_0[triangle_index], 1e-6);
            TS_ASSERT_DELTA(calculated_triangle_normals_0[triangle_index][0], 0.0, 1e-6);
            TS_ASSERT_DELTA(calculated_triangle_normals_0[triangle_index][1], 0.0, 1e-6);
            TS_ASSERT_DELTA(calculated_triangle_normals_0[triangle_index][2], -1.0, 1e-6);
        }

        // Check Passive Centers, normals and areas of face 3
        MonolayerVertexElement<2, 3>* p_face_3 = p_mesh->GetFace(3);
        c_vector<double, 3> passive_center_3 = p_mesh->GetPassiveCenterOfFace(p_face_3);
        TS_ASSERT_DELTA(passive_center_3[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(passive_center_3[1], 0.5, 1e-6);
        TS_ASSERT_DELTA(passive_center_3[2], 0.5, 1e-6);
        TS_ASSERT_DELTA(1.0, p_mesh->CalculateAreaOfFace(p_face_3), 1e-6);
        std::vector<c_vector<double, 3> > calculated_triangle_normals_3;
        std::vector<double> calculated_triangle_areas_3 = p_mesh->CalculateAreasAndNormalsOfTriangulationsOfFace(p_face_3, &calculated_triangle_normals_3);
        for (unsigned triangle_index = 0; triangle_index < calculated_triangle_normals_3.size(); triangle_index++)
        {
            TS_ASSERT_DELTA(0.25, calculated_triangle_areas_3[triangle_index], 1e-6);
            TS_ASSERT_DELTA(calculated_triangle_normals_3[triangle_index][0], -1.0, 1e-6);
            TS_ASSERT_DELTA(calculated_triangle_normals_3[triangle_index][1], 0.0, 1e-6);
            TS_ASSERT_DELTA(calculated_triangle_normals_3[triangle_index][2], 0.0, 1e-6);
        }

        // Check Passive Centers, normals and areas of face 4
        MonolayerVertexElement<2, 3>* p_face_4 = p_mesh->GetFace(4);
        c_vector<double, 3> passive_center_4 = p_mesh->GetPassiveCenterOfFace(p_face_4);
        TS_ASSERT_DELTA(passive_center_4[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(passive_center_4[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(passive_center_4[2], 0.5, 1e-6);
        TS_ASSERT_DELTA(1.0, p_mesh->CalculateAreaOfFace(p_face_4), 1e-6);
        std::vector<c_vector<double, 3> > calculated_triangle_normals_4;
        std::vector<double> calculated_triangle_areas_4 = p_mesh->CalculateAreasAndNormalsOfTriangulationsOfFace(p_face_4, &calculated_triangle_normals_4);
        for (unsigned triangle_index = 0; triangle_index < calculated_triangle_normals_4.size(); triangle_index++)
        {
            TS_ASSERT_DELTA(0.25, calculated_triangle_areas_4[triangle_index], 1e-6);
            TS_ASSERT_DELTA(calculated_triangle_normals_4[triangle_index][0], 0.0, 1e-6);
            TS_ASSERT_DELTA(calculated_triangle_normals_4[triangle_index][1], -1.0, 1e-6);
            TS_ASSERT_DELTA(calculated_triangle_normals_4[triangle_index][2], 0.0, 1e-6);
        }

        // Test GetThicknessOfElement
        double thickness_cube = p_mesh->GetThicknessOfElement(0);
        TS_ASSERT_DELTA(thickness_cube, 2.0, 1e-4);

        // Check if CalculateUnitNormalToFaceWithArea with average unit normal is correct
        c_vector<double, 3> unit_normal_2;
        p_mesh->CalculateUnitNormalToFaceWithArea(p_face_2, unit_normal_2);
        TS_ASSERT_DELTA(unit_normal_2[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_2[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(unit_normal_2[2], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_face_2), 3.0, 1e-6);
    }

    void TestAreaAndVolumeGradients()
    {
        // Create 3D pyramid
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 2.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 2.0, 0.0, 0.0));

        nodes.push_back(new Node<3>(4, false, 0.1, 1.2, 3.0));

        // Construct a cuboid to test non-co-planar volume gradients.

        nodes.push_back(new Node<3>(5, false, 0.0, 0.0, 3.0));
        nodes.push_back(new Node<3>(6, false, 0.0, 1.0, 3.0));
        nodes.push_back(new Node<3>(7, false, 2.0, 1.0, 3.0));
        nodes.push_back(new Node<3>(8, false, 2.0, 0.0, 3.0));

        std::vector<MonolayerVertexElementType> node_types;
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Apical);

        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);

        // Make faces out of these nodes
        std::vector<std::vector<Node<3>*> > nodes_faces(5);
        std::vector<std::vector<MonolayerVertexElementType> > nodes_faces_types(5);

        unsigned node_indices_face_0[4] = { 0, 1, 2, 3 };
        unsigned node_indices_face_1[3] = { 0, 4, 1 };
        unsigned node_indices_face_2[3] = { 1, 4, 2 };
        unsigned node_indices_face_3[3] = { 3, 4, 2 };
        unsigned node_indices_face_4[3] = { 3, 4, 0 };

        for (unsigned i = 0; i < 4; i++)
        {
            nodes_faces[0].push_back(nodes[node_indices_face_0[i]]);
            nodes_faces_types[0].push_back(node_types[node_indices_face_0[i]]);

            if (i < 3)
            {
                nodes_faces[1].push_back(nodes[node_indices_face_1[i]]);
                nodes_faces_types[1].push_back(node_types[node_indices_face_1[i]]);

                nodes_faces[2].push_back(nodes[node_indices_face_2[i]]);
                nodes_faces_types[2].push_back(node_types[node_indices_face_2[i]]);

                nodes_faces[3].push_back(nodes[node_indices_face_3[i]]);
                nodes_faces_types[3].push_back(node_types[node_indices_face_3[i]]);

                nodes_faces[4].push_back(nodes[node_indices_face_4[i]]);
                nodes_faces_types[4].push_back(node_types[node_indices_face_4[i]]);
            }
        }

        // Make the faces
        std::vector<MonolayerVertexElement<2, 3>*> faces;
        std::vector<MonolayerVertexElement<2, 3>*> faces_pyramid;

        // Set face types
        std::vector<MonolayerVertexElementType> face_types;
        face_types.push_back(MonolayerVertexElementType::Basal);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);

        for (unsigned i = 0; i < 5; i++)
        {
            faces.push_back(new MonolayerVertexElement<2, 3>(i, face_types[i], nodes_faces[i], nodes_faces_types[i]));
            faces_pyramid.push_back(faces[i]);
        }

        // Make pyramid elements
        std::vector<bool> orientations(5, false);
        orientations[3] = true;

        std::vector<MonolayerVertexElement<3, 3>*> elements;
        elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, faces_pyramid, orientations, nodes, node_types));

        // Make faces for cuboid
        std::vector<std::vector<Node<3>*> > nodes_faces_cuboid(6);
        std::vector<std::vector<MonolayerVertexElementType> > nodes_faces_types_cuboid(6);

        unsigned node_indices_face_0_cuboid[4] = { 0, 1, 2, 3 };
        unsigned node_indices_face_1_cuboid[4] = { 0, 5, 6, 1 };
        unsigned node_indices_face_2_cuboid[4] = { 1, 6, 7, 2 };
        unsigned node_indices_face_3_cuboid[4] = { 2, 7, 8, 3 };
        unsigned node_indices_face_4_cuboid[4] = { 3, 8, 5, 0 };
        unsigned node_indices_face_5_cuboid[4] = { 5, 8, 7, 6 };

        for (unsigned i = 0; i < 4; i++)
        {
            nodes_faces_cuboid[0].push_back(nodes[node_indices_face_0_cuboid[i]]);
            nodes_faces_types_cuboid[0].push_back(node_types[node_indices_face_0_cuboid[i]]);
            nodes_faces_cuboid[1].push_back(nodes[node_indices_face_1_cuboid[i]]);
            nodes_faces_types_cuboid[1].push_back(node_types[node_indices_face_1_cuboid[i]]);
            nodes_faces_cuboid[2].push_back(nodes[node_indices_face_2_cuboid[i]]);
            nodes_faces_types_cuboid[2].push_back(node_types[node_indices_face_2_cuboid[i]]);
            nodes_faces_cuboid[3].push_back(nodes[node_indices_face_3_cuboid[i]]);
            nodes_faces_types_cuboid[3].push_back(node_types[node_indices_face_3_cuboid[i]]);
            nodes_faces_cuboid[4].push_back(nodes[node_indices_face_4_cuboid[i]]);
            nodes_faces_types_cuboid[4].push_back(node_types[node_indices_face_4_cuboid[i]]);
            nodes_faces_cuboid[5].push_back(nodes[node_indices_face_5_cuboid[i]]);
            nodes_faces_types_cuboid[5].push_back(node_types[node_indices_face_5_cuboid[i]]);
        }

        // Make the faces
        std::vector<MonolayerVertexElement<2, 3>*> faces_cuboid;

        // Set face types
        std::vector<MonolayerVertexElementType> face_types_cuboid;
        face_types_cuboid.push_back(MonolayerVertexElementType::Basal);
        face_types_cuboid.push_back(MonolayerVertexElementType::Lateral);
        face_types_cuboid.push_back(MonolayerVertexElementType::Lateral);
        face_types_cuboid.push_back(MonolayerVertexElementType::Lateral);
        face_types_cuboid.push_back(MonolayerVertexElementType::Lateral);
        face_types_cuboid.push_back(MonolayerVertexElementType::Apical);

        for (unsigned i = 0; i < 6; i++)
        {
            faces.push_back(new MonolayerVertexElement<2, 3>(i + 5, face_types_cuboid[i], nodes_faces_cuboid[i], nodes_faces_types_cuboid[i]));
            faces_cuboid.push_back(faces[i + 5]);
        }

        // Make the elements
        std::vector<bool> orientations_cuboid(6, false);

        elements.push_back(new MonolayerVertexElement<3, 3>(1, MonolayerVertexElementType::Undetermined, faces_cuboid, orientations_cuboid));

        MonolayerVertexMesh<3, 3>* p_mesh = new MonolayerVertexMesh<3, 3>(nodes, faces, elements);

        // Check Area and Volume gradient

        double x_apex = 0.1;
        double y_apex = 1.2;
        double z_apex = 3.0;

        Node<3>* apex_node = p_mesh->GetElement(0)->GetNode(4);

        TS_ASSERT_DELTA(apex_node->rGetLocation()[0], x_apex, 1e-6);
        TS_ASSERT_DELTA(apex_node->rGetLocation()[1], y_apex, 1e-6);
        TS_ASSERT_DELTA(apex_node->rGetLocation()[2], z_apex, 1e-6);

        // Test lateral faces
        MonolayerVertexElement<2, 3>* p_face_0 = p_mesh->GetElement(0)->GetFace(1);
        MonolayerVertexElement<2, 3>* p_face_1 = p_mesh->GetElement(0)->GetFace(2);
        MonolayerVertexElement<2, 3>* p_face_2 = p_mesh->GetElement(0)->GetFace(3);
        MonolayerVertexElement<2, 3>* p_face_3 = p_mesh->GetElement(0)->GetFace(4);

        c_vector<double, 3> grad_0 = p_mesh->GetAreaGradientOfFaceAtNode(p_face_0, 1);
        // std::cout << "\nReturned Gradient: " << grad_0[0] << ", " << grad_0[1] << ", " << grad_0[2] << ", ";

        c_vector<double, 3> grad_1 = p_mesh->GetAreaGradientOfFaceAtNode(p_face_1, 1);
        c_vector<double, 3> grad_2 = p_mesh->GetAreaGradientOfFaceAtNode(p_face_2, 1);
        c_vector<double, 3> grad_3 = p_mesh->GetAreaGradientOfFaceAtNode(p_face_3, 1);

        TS_ASSERT_DELTA(grad_0[0], x_apex / (2.0 * sqrt(x_apex * x_apex + z_apex * z_apex)), 1e-6);
        TS_ASSERT_DELTA(grad_0[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(grad_0[2], z_apex / (2.0 * sqrt(x_apex * x_apex + z_apex * z_apex)), 1e-6);

        TS_ASSERT_DELTA(grad_1[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(grad_1[1], (y_apex - 1.0) / (sqrt((y_apex - 1.0) * (y_apex - 1.0) + z_apex * z_apex)), 1e-6);
        TS_ASSERT_DELTA(grad_1[2], z_apex / (sqrt((y_apex - 1.0) * (y_apex - 1.0) + z_apex * z_apex)), 1e-6);

        TS_ASSERT_DELTA(grad_2[0], (x_apex - 2.0) / (2.0 * sqrt((x_apex - 2.0) * (x_apex - 2.0) + z_apex * z_apex)), 1e-6);
        TS_ASSERT_DELTA(grad_2[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(grad_2[2], z_apex / (2.0 * sqrt((x_apex - 2.0) * (x_apex - 2.0) + z_apex * z_apex)), 1e-6);

        TS_ASSERT_DELTA(grad_3[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(grad_3[1], y_apex / (sqrt(y_apex * y_apex + z_apex * z_apex)), 1e-6);
        TS_ASSERT_DELTA(grad_3[2], z_apex / (sqrt(y_apex * y_apex + z_apex * z_apex)), 1e-6);

        // Test base face

        MonolayerVertexElement<2, 3>* p_face_base = p_mesh->GetElement(0)->GetFace(0);
        c_vector<double, 3> grad_base = p_mesh->GetAreaGradientOfFaceAtNode(p_face_base, 0);

        TS_ASSERT_DELTA(grad_base[0], -0.5, 1e-6);
        TS_ASSERT_DELTA(grad_base[1], -1.0, 1e-6);
        TS_ASSERT_DELTA(grad_base[2], 0.0, 1e-6);

        // Test volume gradient

        c_vector<double, 3> volume_grad_apex = p_mesh->GetVolumeGradientAtNode(0, 4);

        TS_ASSERT_DELTA(volume_grad_apex[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(volume_grad_apex[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(volume_grad_apex[2], 2.0 / 3.0, 1e-6);

        c_vector<double, 3> volume_grad_base = p_mesh->GetVolumeGradientAtNode(0, 0);

        TS_ASSERT_DELTA(volume_grad_base[0], -0.5, 1e-6);
        TS_ASSERT_DELTA(volume_grad_base[1], -1.0, 1e-6);

        /**
         * Check a more complex situation with a non-co-planar base
         **/

        nodes[0]->SetPoint(ChastePoint<3>(0.0, 0.0, -1.0));

        grad_base = p_mesh->GetAreaGradientOfFaceAtNode(p_face_base, 0);
        // Triangle areas
        double A1 = sqrt(1.0 + 1.0 / 16.0) / 2.0;
        double A2 = sqrt(1.0 / 4.0 + 1.0 / 16.0);

        // Coefficients in the cubic area formula
        double zeta_1_3 = (26.0 * 26.0 - 2.0 * 82.0 - 2.0 * 16.0 * 16.0) / (16.0 * 16.0);
        double zeta_2_3 = (26.0 * 14.0 - 228.0) / (32.0);
        double zeta_3_3 = (4.0 * 49.0 - 2.0 * 66.0) / (16.0);
        double A3 = sqrt(zeta_1_3 + zeta_2_3 + zeta_3_3) / 4.0;

        double zeta_1_4 = (26.0 * 26.0 - 2.0 * 82.0 - 2.0 * 16.0 * 16.0) / (16.0 * 16.0);
        double zeta_2_4 = (26.0 * 26.0 - 100.0 - 16.0 * 32.0) / (32.0);
        double zeta_3_4 = (4.0 * 13.0 * 13.0 - 100.0 - 2.0 * 16.0 * 16.0) / (16.0);
        double A4 = sqrt(zeta_1_4 + zeta_2_4 + zeta_3_4) / 4.0;

        double grad_A1 = -1.0 / (64.0 * A1);
        double grad_A2 = -1.0 / (16.0 * A2);
        double grad_A3 = -1.0 / (16.0 * A3) * (2.0 * zeta_1_3 + zeta_2_3);
        double grad_A4 = -1.0 / (16.0 * A4) * (2.0 * zeta_1_4 + zeta_2_4);

        TS_ASSERT_DELTA(grad_base[2], grad_A1 + grad_A2 + grad_A3 + grad_A4, 1e-6);

        // Check volume gradient for cuboid

        c_vector<double, 3> volume_grad_cuboid = p_mesh->GetVolumeGradientAtNode(1, 0);

        TS_ASSERT_DELTA(volume_grad_cuboid[2], 0.5, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(1), 6.5, 1e-6);

        // Avoid memory leaks
        delete p_mesh;
    }

    void TestElementFacesMap()
    {
        // Test method with a single triangular element
        std::vector<Node<3>*> triangle_nodes;
        triangle_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        triangle_nodes.push_back(new Node<3>(1, false, 2.0, 0.0, 0.0));
        triangle_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 3.0));
        std::vector<MonolayerVertexElementType> triangle_node_types(3, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 3>*> triangle_faces;
        triangle_faces.push_back(new MonolayerVertexElement<2, 3>(10, MonolayerVertexElementType::Undetermined, triangle_nodes, triangle_node_types));

        MonolayerVertexElement<2, 3>* p_triangle_face = triangle_faces[0];

        // Test method with a 3D cube and pyramid
        MonolayerVertexMesh<3, 3>* p_mesh = ConstructCubeAndPyramidMesh();

        MonolayerVertexElement<3, 3>* p_cube = p_mesh->GetElement(0);
        MonolayerVertexElement<3, 3>* p_pyramid = p_mesh->GetElement(1);

        // Create the ElementsFacesMap
        p_mesh->UpdateElementsFacesMap();

        std::map<MonolayerVertexElement<3, 3>*, std::vector<unsigned> >* p_map = p_mesh->GetElementsFacesMap();

        TS_ASSERT_EQUALS((*p_map)[p_cube].size(), 6u);
        TS_ASSERT_EQUALS((*p_map)[p_pyramid].size(), 5u);

        std::vector<unsigned> faces_in_cube = (*p_map)[p_cube];
        // Face we will add
        TS_ASSERT(std::find(faces_in_cube.begin(), faces_in_cube.end(), 10) == faces_in_cube.end());
        // Cuboidal faces
        TS_ASSERT(std::find(faces_in_cube.begin(), faces_in_cube.end(), 0) != faces_in_cube.end());
        TS_ASSERT(std::find(faces_in_cube.begin(), faces_in_cube.end(), 1) != faces_in_cube.end());
        // Pyramid faces
        TS_ASSERT(std::find(faces_in_cube.begin(), faces_in_cube.end(), 7) == faces_in_cube.end());
        TS_ASSERT(std::find(faces_in_cube.begin(), faces_in_cube.end(), 8) == faces_in_cube.end());

        // Add a face to the cube
        p_cube->AddFace(p_triangle_face, MonolayerVertexElementType::Undetermined, false);
        p_mesh->UpdateElementsFacesMapOfElement(0);

        TS_ASSERT_EQUALS((*p_map)[p_cube].size(), 7u);
        TS_ASSERT_EQUALS((*p_map)[p_pyramid].size(), 5u);

        faces_in_cube = (*p_map)[p_cube];
        TS_ASSERT(std::find(faces_in_cube.begin(), faces_in_cube.end(), 10) != faces_in_cube.end());

        // Delete a face from pyramid
        unsigned global_index_pyramid_face2 = p_pyramid->GetFace(2)->GetIndex();
        std::vector<unsigned> faces_in_pyramid = (*p_map)[p_pyramid];
        TS_ASSERT(std::find(faces_in_pyramid.begin(), faces_in_pyramid.end(), global_index_pyramid_face2) != faces_in_pyramid.end());

        p_pyramid->DeleteFace(2);
        p_mesh->UpdateElementsFacesMap();

        TS_ASSERT_EQUALS((*p_map)[p_cube].size(), 7u);
        TS_ASSERT_EQUALS((*p_map)[p_pyramid].size(), 4u);

        faces_in_pyramid = (*p_map)[p_pyramid];
        TS_ASSERT(std::find(faces_in_pyramid.begin(), faces_in_pyramid.end(), global_index_pyramid_face2) == faces_in_pyramid.end());

        faces_in_cube = (*p_map)[p_cube];
        TS_ASSERT(std::find(faces_in_cube.begin(), faces_in_cube.end(), 10) != faces_in_cube.end());
    }

    void TestCalculateMomentsOfElement()
    {
        // Useful unit vectors
        c_vector<double, 3> x_unit, y_unit, z_unit;
        x_unit <<= 1.0, 0.0, 0.0;
        y_unit <<= 0.0, 1.0, 0.0;
        z_unit <<= 0.0, 0.0, 1.0;

        // Test method with a single triangular element
        std::vector<Node<3>*> isos_triangle_nodes;
        isos_triangle_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        isos_triangle_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        isos_triangle_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        std::vector<MonolayerVertexElementType> isos_triangle_types(3, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 3>*> isos_triangle_faces;
        std::vector<MonolayerVertexElement<3, 3>*> isos_triangle_elements;
        isos_triangle_faces.push_back(new MonolayerVertexElement<2, 3>(0, MonolayerVertexElementType::Undetermined, isos_triangle_nodes, isos_triangle_types));
        std::vector<bool> isos_triangle_faces_orientations(1, true);
        isos_triangle_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, isos_triangle_faces, isos_triangle_faces_orientations));
        MonolayerVertexMesh<3, 3> isos_triangle_mesh(isos_triangle_nodes, isos_triangle_faces, isos_triangle_elements);

        c_vector<double, 3> isos_triangle_moments = isos_triangle_mesh.CalculateMomentsOfFace(isos_triangle_faces[0], z_unit, x_unit, y_unit);
        TS_ASSERT_DELTA(isos_triangle_moments[0], 0.0277, 1e-4);
        TS_ASSERT_DELTA(isos_triangle_moments[1], 0.0277, 1e-4);
        TS_ASSERT_DELTA(isos_triangle_moments[2], -0.0138, 1e-4);

        // Test method with a single triangular element
        std::vector<Node<3>*> triangle_nodes;
        triangle_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        triangle_nodes.push_back(new Node<3>(1, false, 2.0, 0.0, 0.0));
        triangle_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        std::vector<MonolayerVertexElementType> triangle_types(3, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 3>*> triangle_faces;
        std::vector<MonolayerVertexElement<3, 3>*> triangle_elements;
        triangle_faces.push_back(new MonolayerVertexElement<2, 3>(0, MonolayerVertexElementType::Undetermined, triangle_nodes, triangle_types));
        std::vector<bool> triangle_faces_orientations(1, true);
        triangle_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, triangle_faces, triangle_faces_orientations));
        MonolayerVertexMesh<3, 3> triangle_mesh(triangle_nodes, triangle_faces, triangle_elements);

        c_vector<double, 3> triangle_moments = triangle_mesh.CalculateMomentsOfFace(triangle_faces[0], z_unit, x_unit, y_unit);
        TS_ASSERT_DELTA(triangle_moments(0), 1.0 / 18.0, 1e-6); // Ixx
        TS_ASSERT_DELTA(triangle_moments(1), 2.0 / 9.0, 1e-6); // Iyy
        TS_ASSERT_DELTA(triangle_moments(2), -5.0 / 90.0, 1e-6); // Ixy

        // Test method with a single rectangular element parallel to the x-axis
        std::vector<Node<3>*> horizontal_rectangle_nodes;
        horizontal_rectangle_nodes.push_back(new Node<3>(0, false, 2.0, 1.0, 3.0));
        horizontal_rectangle_nodes.push_back(new Node<3>(1, false, -2.0, 1.0, 3.0));
        horizontal_rectangle_nodes.push_back(new Node<3>(2, false, -2.0, -1.0, 3.0));
        horizontal_rectangle_nodes.push_back(new Node<3>(3, false, 2.0, -1.0, 3.0));
        std::vector<MonolayerVertexElementType> horizontal_rectangle_types(4, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 3>*> horizontal_rectangle_faces;
        std::vector<MonolayerVertexElement<3, 3>*> horizontal_rectangle_elements;
        horizontal_rectangle_faces.push_back(new MonolayerVertexElement<2, 3>(0, MonolayerVertexElementType::Undetermined, horizontal_rectangle_nodes, horizontal_rectangle_types));
        std::vector<bool> horizontal_rectangle_faces_orientations(1, true);
        horizontal_rectangle_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, horizontal_rectangle_faces, horizontal_rectangle_faces_orientations));
        MonolayerVertexMesh<3, 3> horizontal_rectangle_mesh(horizontal_rectangle_nodes, horizontal_rectangle_faces, horizontal_rectangle_elements);

        c_vector<double, 3> horizontal_rectangle_moments = horizontal_rectangle_mesh.CalculateMomentsOfFace(horizontal_rectangle_faces[0], z_unit, x_unit, y_unit);
        TS_ASSERT_DELTA(horizontal_rectangle_moments(0), 8.0 / 3.0, 1e-6); // Ixx
        TS_ASSERT_DELTA(horizontal_rectangle_moments(1), 32.0 / 3.0, 1e-6); // Iyy
        TS_ASSERT_DELTA(horizontal_rectangle_moments(2), 0.0, 1e-6); // Ixy = 0 by symmetry

        // Test method with the same shape, but supply nodes in clockwise manner
        std::vector<Node<3>*> clockwise_rectangle_nodes;
        clockwise_rectangle_nodes.push_back(new Node<3>(0, false, -2.0, -1.0, 5.0));
        clockwise_rectangle_nodes.push_back(new Node<3>(1, false, -2.0, 1.0, 5.0));
        clockwise_rectangle_nodes.push_back(new Node<3>(2, false, 2.0, 1.0, 5.0));
        clockwise_rectangle_nodes.push_back(new Node<3>(3, false, 2.0, -1.0, 5.0));
        std::vector<MonolayerVertexElementType> clockwise_rectangle_types(4, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 3>*> clockwise_rectangle_faces;
        std::vector<MonolayerVertexElement<3, 3>*> clockwise_rectangle_elements;
        clockwise_rectangle_faces.push_back(new MonolayerVertexElement<2, 3>(0, MonolayerVertexElementType::Undetermined, clockwise_rectangle_nodes, clockwise_rectangle_types));
        std::vector<bool> clockwise_rectangle_faces_orientations(1, true);
        clockwise_rectangle_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, clockwise_rectangle_faces, clockwise_rectangle_faces_orientations));
        MonolayerVertexMesh<3, 3> clockwise_rectangle_mesh(clockwise_rectangle_nodes, clockwise_rectangle_faces, clockwise_rectangle_elements);

        c_vector<double, 3> clockwise_rectangle_moments = clockwise_rectangle_mesh.CalculateMomentsOfFace(clockwise_rectangle_faces[0], z_unit, x_unit, y_unit);
        TS_ASSERT_DELTA(clockwise_rectangle_moments(0), 8.0 / 3.0, 1e-6); // Ixx
        TS_ASSERT_DELTA(clockwise_rectangle_moments(1), 32.0 / 3.0, 1e-6); // Iyy
        TS_ASSERT_DELTA(clockwise_rectangle_moments(2), 0.0, 1e-6); // Ixy = 0 by symmetry

        // Test method with a single rectangular element at a 30 degree angle to the x-axis
        std::vector<Node<3>*> angled_rectangle_nodes;
        angled_rectangle_nodes.push_back(new Node<3>(0, false, 2.0 * 0.5 * sqrt(3.0) - 1.0 * 0.5, 2.0 * 0.5 + 1.0 * 0.5 * sqrt(3.0), 2.0));
        angled_rectangle_nodes.push_back(new Node<3>(1, false, -2.0 * 0.5 * sqrt(3.0) - 1.0 * 0.5, -2.0 * 0.5 + 1.0 * 0.5 * sqrt(3.0), 2.0));
        angled_rectangle_nodes.push_back(new Node<3>(2, false, -2.0 * 0.5 * sqrt(3.0) + 1.0 * 0.5, -2.0 * 0.5 - 1.0 * 0.5 * sqrt(3.0), 2.0));
        angled_rectangle_nodes.push_back(new Node<3>(3, false, 2.0 * 0.5 * sqrt(3.0) + 1.0 * 0.5, 2.0 * 0.5 - 1.0 * 0.5 * sqrt(3.0), 2.0));
        std::vector<MonolayerVertexElementType> angled_rectangle_types(4, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 3>*> angled_rectangle_faces;
        std::vector<MonolayerVertexElement<3, 3>*> angled_rectangle_elements;
        angled_rectangle_faces.push_back(new MonolayerVertexElement<2, 3>(0, MonolayerVertexElementType::Undetermined, angled_rectangle_nodes, angled_rectangle_types));
        std::vector<bool> angled_rectangle_faces_orientations(1, true);
        angled_rectangle_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, angled_rectangle_faces, angled_rectangle_faces_orientations));
        MonolayerVertexMesh<3, 3> angled_rectangle_mesh(angled_rectangle_nodes, angled_rectangle_faces, angled_rectangle_elements);

        c_vector<double, 3> angled_rectangle_moments = angled_rectangle_mesh.CalculateMomentsOfFace(angled_rectangle_faces[0], z_unit, x_unit, y_unit);
        TS_ASSERT_DELTA(angled_rectangle_moments[0], 14.0 / 3.0, 1e-4); // Ixx
        TS_ASSERT_DELTA(angled_rectangle_moments[1], 26.0 / 3.0, 1e-4); // Iyy
        TS_ASSERT_DELTA(angled_rectangle_moments[2], 2.0 * sqrt(3.0), 1e-4); // Ixy

        // Test method with a non-coplanar rectangle
        std::vector<Node<3>*> nonplanar_rectangle_nodes;
        nonplanar_rectangle_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 5.0));
        nonplanar_rectangle_nodes.push_back(new Node<3>(1, false, 1.0, -1.0, 5.0));
        nonplanar_rectangle_nodes.push_back(new Node<3>(2, false, 2.0, 0.0, 5.0 + 4.0 / sqrt(2.0)));
        nonplanar_rectangle_nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 5.0));

        std::vector<MonolayerVertexElementType> nonplanar_rectangle_types(4, MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElement<2, 3>*> nonplanar_rectangle_faces;
        std::vector<MonolayerVertexElement<3, 3>*> nonplanar_rectangle_elements;
        nonplanar_rectangle_faces.push_back(new MonolayerVertexElement<2, 3>(0, MonolayerVertexElementType::Undetermined, nonplanar_rectangle_nodes, nonplanar_rectangle_types));
        std::vector<bool> nonplanar_rectangle_faces_orientations(1, true);
        nonplanar_rectangle_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, nonplanar_rectangle_faces, nonplanar_rectangle_faces_orientations));
        MonolayerVertexMesh<3, 3> nonplanar_rectangle_mesh(nonplanar_rectangle_nodes, nonplanar_rectangle_faces, nonplanar_rectangle_elements);

        // (Mean) Normal
        c_vector<double, 3> nrml, orth_1, orth_2;
        nrml <<= -8.0 / sqrt(2.0), 0.0, 4.0;
        nrml /= norm_2(nrml);
        // Orthogonals
        orth_1 <<= 0.0, 1.0, 0.0;
        orth_2 = VectorProduct(nrml, orth_1);
        orth_2 /= norm_2(orth_2);

        // Project the nodes on the orthogonal plane
        c_vector<double, 3> aminp, bminp, cminp, dminp;
        aminp <<= -1.0, 0.0, -1.0 / sqrt(2.0);
        bminp <<= 0.0, -1.0, -1.0 / sqrt(2.0);
        cminp <<= 1.0, 0.0, 3.0 / sqrt(2.0);
        dminp <<= 0.0, 1.0, -1.0 / sqrt(2.0);

        aminp = aminp - inner_prod(aminp, nrml) * nrml;
        bminp = bminp - inner_prod(bminp, nrml) * nrml;
        cminp = cminp - inner_prod(cminp, nrml) * nrml;
        dminp = dminp - inner_prod(dminp, nrml) * nrml;

        // Calculate decomposition into orthogonals
        double acoeff[2], bcoeff[2], ccoeff[2], dcoeff[2];
        acoeff[0] = aminp(1);
        bcoeff[0] = bminp(1);
        ccoeff[0] = cminp(1);
        dcoeff[0] = dminp(1);

        acoeff[1] = inner_prod(aminp, orth_2);
        bcoeff[1] = inner_prod(bminp, orth_2);
        ccoeff[1] = inner_prod(cminp, orth_2);
        dcoeff[1] = inner_prod(dminp, orth_2);

        // Calculate Moments with the formula
        double Iyy_ab = 1.0 / 12.0 * (acoeff[0] * bcoeff[1] - bcoeff[0] * acoeff[1]) * (acoeff[0] * acoeff[0] + acoeff[0] * bcoeff[0] + bcoeff[0] * bcoeff[0]);
        double Ixx_ab = 1.0 / 12.0 * (acoeff[0] * bcoeff[1] - bcoeff[0] * acoeff[1]) * (acoeff[1] * acoeff[1] + acoeff[1] * bcoeff[1] + bcoeff[1] * bcoeff[1]);
        double Ixy_ab = 1.0 / 24.0 * (acoeff[0] * bcoeff[1] - bcoeff[0] * acoeff[1]) * (acoeff[0] * bcoeff[1] + 2.0 * acoeff[0] * acoeff[1] + 2.0 * bcoeff[0] * bcoeff[1] + bcoeff[0] * acoeff[1]);

        double Iyy_bc = 1.0 / 12.0 * (bcoeff[0] * ccoeff[1] - ccoeff[0] * bcoeff[1]) * (bcoeff[0] * bcoeff[0] + bcoeff[0] * ccoeff[0] + ccoeff[0] * ccoeff[0]);
        double Ixx_bc = 1.0 / 12.0 * (bcoeff[0] * ccoeff[1] - ccoeff[0] * bcoeff[1]) * (bcoeff[1] * bcoeff[1] + bcoeff[1] * ccoeff[1] + ccoeff[1] * ccoeff[1]);
        double Ixy_bc = 1.0 / 24.0 * (bcoeff[0] * ccoeff[1] - ccoeff[0] * bcoeff[1]) * (bcoeff[0] * ccoeff[1] + 2.0 * bcoeff[0] * bcoeff[1] + 2.0 * ccoeff[0] * ccoeff[1] + ccoeff[0] * bcoeff[1]);

        double Iyy_cd = 1.0 / 12.0 * (ccoeff[0] * dcoeff[1] - dcoeff[0] * ccoeff[1]) * (ccoeff[0] * ccoeff[0] + ccoeff[0] * dcoeff[0] + dcoeff[0] * dcoeff[0]);
        double Ixx_cd = 1.0 / 12.0 * (ccoeff[0] * dcoeff[1] - dcoeff[0] * ccoeff[1]) * (ccoeff[1] * ccoeff[1] + ccoeff[1] * dcoeff[1] + dcoeff[1] * dcoeff[1]);
        double Ixy_cd = 1.0 / 24.0 * (ccoeff[0] * dcoeff[1] - dcoeff[0] * ccoeff[1]) * (ccoeff[0] * dcoeff[1] + 2.0 * ccoeff[0] * ccoeff[1] + 2.0 * dcoeff[0] * dcoeff[1] + dcoeff[0] * ccoeff[1]);

        double Iyy_da = 1.0 / 12.0 * (dcoeff[0] * acoeff[1] - acoeff[0] * dcoeff[1]) * (dcoeff[0] * dcoeff[0] + dcoeff[0] * acoeff[0] + acoeff[0] * acoeff[0]);
        double Ixx_da = 1.0 / 12.0 * (dcoeff[0] * acoeff[1] - acoeff[0] * dcoeff[1]) * (dcoeff[1] * dcoeff[1] + dcoeff[1] * acoeff[1] + acoeff[1] * acoeff[1]);
        double Ixy_da = 1.0 / 24.0 * (dcoeff[0] * acoeff[1] - acoeff[0] * dcoeff[1]) * (dcoeff[0] * acoeff[1] + 2.0 * dcoeff[0] * dcoeff[1] + 2.0 * acoeff[0] * acoeff[1] + acoeff[0] * dcoeff[1]);

        double Ixx = Ixx_ab + Ixx_bc + Ixx_cd + Ixx_da;
        double Iyy = Iyy_ab + Iyy_bc + Iyy_cd + Iyy_da;
        double Ixy = Ixy_ab + Ixy_bc + Ixy_cd + Ixy_da;

        Ixx = Ixx < 0.0 ? -1.0 * Ixx : Ixx;
        Iyy = Ixx < 0.0 ? -1.0 * Iyy : Iyy;
        Ixy = Ixx < 0.0 ? -1.0 * Ixy : Ixy;

        c_vector<double, 3> nonplanar_rectangle_moments = nonplanar_rectangle_mesh.CalculateMomentsOfFace(nonplanar_rectangle_faces[0], nrml, orth_1, orth_2);
        TS_ASSERT_DELTA(nonplanar_rectangle_moments(0), Ixx, 1e-6); // Ixx
        TS_ASSERT_DELTA(nonplanar_rectangle_moments(1), Iyy, 1e-6); // Iyy
        TS_ASSERT_DELTA(nonplanar_rectangle_moments(2), Ixy, 1e-6); // Ixy
    }

    void TestGetMeanApicalBasalLongAxisOfElement()
    {
        // We test this with an elongated hexagon with long axis parallel to the x-axis
        /*
         *             /o-------o\
         *            / :_ _ _ _: \
         *           o /         \ o
         *           |\           /|
         *            \\o-------o/ /
         *             \|_______| /
         */
        std::vector<Node<3>*> regular_nodes;
        // Basal
        regular_nodes.push_back(new Node<3>(0, false, 3.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(1, false, 2.5, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(2, false, -2.5, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(3, false, -3.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(4, false, -2.5, -1.0, 0.0));
        regular_nodes.push_back(new Node<3>(5, false, 2.5, -1.0, 0.0));
        // Apical
        regular_nodes.push_back(new Node<3>(6, false, 3.0, 0.0, 1.0));
        regular_nodes.push_back(new Node<3>(7, false, 2.5, 1.0, 1.0));
        regular_nodes.push_back(new Node<3>(8, false, -2.5, 1.0, 1.0));
        regular_nodes.push_back(new Node<3>(9, false, -3.0, 0.0, 1.0));
        regular_nodes.push_back(new Node<3>(10, false, -2.5, -1.0, 1.0));
        regular_nodes.push_back(new Node<3>(11, false, 2.5, -1.0, 1.0));

        std::vector<MonolayerVertexElementType> regular_types(6, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> regular_types_2(6, MonolayerVertexElementType::Apical);
        regular_types.insert(regular_types.end(), regular_types_2.begin(), regular_types_2.end());

        std::vector<std::vector<int> > face_indices;
        std::vector<int> face_0 = { 5, 4, 3, 2, 1, 0 };
        std::vector<int> face_1 = { 6, 7, 8, 9, 10, 11 };
        std::vector<int> face_2 = { 0, 5, 11, 6 };
        std::vector<int> face_3 = { 5, 4, 10, 11 };
        std::vector<int> face_4 = { 4, 3, 9, 10 };
        std::vector<int> face_5 = { 3, 2, 8, 9 };
        std::vector<int> face_6 = { 2, 1, 7, 8 };
        std::vector<int> face_7 = { 1, 0, 6, 7 };
        face_indices = { face_0, face_1, face_2, face_3, face_4, face_5, face_6, face_7 };

        std::vector<MonolayerVertexElement<2, 3>*> regular_faces;
        for (unsigned i = 0; i < 8; i++)
        {
            std::vector<Node<3>*> face_nodes;
            std::vector<MonolayerVertexElementType> node_types;
            for (auto it = face_indices[i].begin(); it != face_indices[i].end(); ++it)
            {
                face_nodes.push_back(regular_nodes[*it]);
                node_types.push_back(regular_types[*it]);
            }
            MonolayerVertexElementType facetype;
            if (i == 0)
                facetype = MonolayerVertexElementType::Basal;
            else if (i == 1)
                facetype = MonolayerVertexElementType::Apical;
            else
                facetype = MonolayerVertexElementType::Lateral;

            regular_faces.push_back(new MonolayerVertexElement<2, 3>(i, facetype, face_nodes, node_types));
        }
        std::vector<MonolayerVertexElement<3, 3>*> regular_elements;
        std::vector<bool> regular_faces_orientations(8, false);

        regular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, regular_faces, regular_faces_orientations));
        MonolayerVertexMesh<3, 3> regular_mesh(regular_nodes, regular_faces, regular_elements);

        c_vector<double, 3> long_axis = regular_mesh.GetMeanApicalBasalLongAxisOfElement(0);
        long_axis = long_axis(0) < 0.0 ? -1.0 * long_axis : long_axis;

        TS_ASSERT_DELTA(long_axis(0), 1.0, 1e-6);
        TS_ASSERT_DELTA(long_axis(1), 0.0, 1e-6);
        TS_ASSERT_DELTA(long_axis(2), 0.0, 1e-6);

        // Now check also whether mid-plane area is correct
        double mid_plane_area = regular_mesh.GetMidPlaneAreaOfElement(0);

        TS_ASSERT_DELTA(mid_plane_area, 11.0, 1e-6);

        // Now change one node to get a z-component
        c_vector<double, 3>& node_6 = regular_nodes[6]->rGetModifiableLocation();
        node_6[2] = 2.0;

        c_vector<double, 3> analytical_long_basal;
        analytical_long_basal <<= 1.0, 0.0, 0.0;

        // Check normal calculation
        c_vector<double, 3> nrml;
        regular_mesh.CalculateUnitNormalToFaceWithArea(regular_faces[1], nrml);

        TS_ASSERT_DELTA(nrml(0), -2.0 / sqrt(2.0 * 2.0 + 22.0 * 22.0), 1e-6);
        TS_ASSERT_DELTA(nrml(1), 0.0, 1e-6);
        TS_ASSERT_DELTA(nrml(2), 22.0 / sqrt(2.0 * 2.0 + 22.0 * 22.0), 1e-6);

        c_vector<double, 3> analytical_long_apical;
        // paper calculation
        analytical_long_apical <<= 22.0 / sqrt(2.0 * 2.0 + 22.0 * 22.0), 0, 2.0 / sqrt(2.0 * 2.0 + 22.0 * 22.0);

        c_vector<double, 3> analytical_long = analytical_long_basal + analytical_long_apical;
        analytical_long /= norm_2(analytical_long);

        long_axis = regular_mesh.GetMeanApicalBasalLongAxisOfElement(0);
        long_axis = long_axis(0) < 0.0 ? -1.0 * long_axis : long_axis;

        TS_ASSERT_DELTA(long_axis(0), analytical_long(0), 1e-6);
        TS_ASSERT_DELTA(long_axis(1), analytical_long(1), 1e-6);
        TS_ASSERT_DELTA(long_axis(2), analytical_long(2), 1e-6);

        // We rotate the system with a rotation matrix with two angles phi and psi
        double psi = M_PI / 6.0;
        double phi = M_PI / 3.0;

        // rotate around x
        boost::numeric::ublas::matrix<double> rotation_mat_1(3, 3, 0.0);
        rotation_mat_1(0, 0) = 1.0;
        rotation_mat_1(1, 1) = cos(psi);
        rotation_mat_1(2, 2) = cos(psi);
        rotation_mat_1(1, 2) = -sin(psi);
        rotation_mat_1(2, 1) = sin(psi);

        // rotate around z
        boost::numeric::ublas::matrix<double> rotation_mat_2(3, 3, 0.0);
        rotation_mat_2(2, 2) = 1.0;
        rotation_mat_2(0, 0) = cos(phi);
        rotation_mat_2(1, 1) = cos(phi);
        rotation_mat_2(0, 1) = -sin(phi);
        rotation_mat_2(1, 0) = sin(phi);

        // Rotate nodes
        for (auto it = regular_nodes.begin(); it != regular_nodes.end(); ++it)
        {
            c_vector<double, 3>& rPosition = (*it)->rGetModifiableLocation();
            rPosition = prod(prod(rotation_mat_2, rotation_mat_1), rPosition);
        }

        // Rotate result
        analytical_long = prod(prod(rotation_mat_2, rotation_mat_1), analytical_long);

        // Check if it still works
        long_axis = regular_mesh.GetMeanApicalBasalLongAxisOfElement(0);
        long_axis = long_axis(0) < 0.0 ? -1.0 * long_axis : long_axis;

        TS_ASSERT_DELTA(long_axis(0), analytical_long(0), 1e-6);
        TS_ASSERT_DELTA(long_axis(1), analytical_long(1), 1e-6);
        TS_ASSERT_DELTA(long_axis(2), analytical_long(2), 1e-6);
    }

    void TestGetLumenVolume()
    {
    }
};

#endif /*TESTVERTEXMESH_HPP_*/