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

#ifndef TESTMUTABLEMONOLAYERVERTEXMESH_HPP_
#define TESTMUTABLEMONOLAYERVERTEXMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "ArchiveOpener.hpp"
#include "FiniteThicknessHoneycombVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshReader.hpp"
#include "MonolayerVertexMeshWriter.hpp"
#include "MutableMonolayerVertexMesh.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestMutableMonolayerVertexMesh : public CxxTest::TestSuite
{
private:
    MutableMonolayerVertexMesh<3, 3>* ConstructCubeAndPyramidMesh()
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

        return new MutableMonolayerVertexMesh<3, 3>(nodes, faces, elements);
    }

public:
    void TestMutableMonolayerVertexElementIterator()
    {
        // Create mesh
        MonolayerVertexMeshReader<2, 2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_3_by_3");
        MutableMonolayerVertexMesh<2, 2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 9u);

        unsigned counter = 0;
        for (MutableMonolayerVertexMesh<2, 2>::MonolayerVertexElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter, element_index); // assumes the iterator will give elements 0,1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(mesh.GetNumElements(), counter);

        // For coverage, test with an empty mesh
        MutableMonolayerVertexMesh<2, 2> empty_mesh;

        // Since the mesh is empty, the iterator should be set to mrMesh.mNodes.end() when constructed
        MutableMonolayerVertexMesh<2, 2>::MonolayerVertexElementIterator el_iter = empty_mesh.GetElementIteratorBegin();

        // Check that the iterator is now at the end (we need to check this as a double-negative,
        // as we only have a NOT-equals operator defined on the iterator).
        bool iter_is_not_at_end = (el_iter != empty_mesh.GetElementIteratorEnd());
        TS_ASSERT_EQUALS(iter_is_not_at_end, false);

        // Delete an element from mesh and test the iterator
        mesh.DeleteElementPriorToReMesh(0);

        counter = 0;
        for (MutableMonolayerVertexMesh<2, 2>::MonolayerVertexElementIterator iter = mesh.GetElementIteratorBegin();
             iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            unsigned element_index = iter->GetIndex();
            TS_ASSERT_EQUALS(counter + 1, element_index); // assumes the iterator will give elements 1..,N in that order
            counter++;
        }

        TS_ASSERT_EQUALS(mesh.GetNumElements(), counter);
        TS_ASSERT_EQUALS(mesh.GetNumAllElements(), counter + 1);
        TS_ASSERT_EQUALS(mesh.IsMeshChanging(), true);
    }

    void TestBasic2dMutableMonolayerVertexMesh()
    {
        // Make seven nodes to assign to two elements
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        basic_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        basic_nodes.push_back(new Node<2>(2, false, 1.5, 1.0));
        basic_nodes.push_back(new Node<2>(3, false, 1.0, 2.0));
        basic_nodes.push_back(new Node<2>(4, false, 0.0, 1.0));
        basic_nodes.push_back(new Node<2>(5, false, 2.0, 0.0));
        basic_nodes.push_back(new Node<2>(6, false, 2.0, 3.0));

        // Make two triangular elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[5] = { 0, 1, 2, 3, 4 };
        unsigned node_indices_elem_1[3] = { 2, 5, 6 };
        for (unsigned i = 0; i < 5; i++)
        {
            nodes_elem_0.push_back(basic_nodes[node_indices_elem_0[i]]);
            if (i < 3)
            {
                nodes_elem_1.push_back(basic_nodes[node_indices_elem_1[i]]);
            }
        }

        std::vector<MonolayerVertexElement<2, 2>*> basic_vertex_elements;
        std::vector<MonolayerVertexElementType> node_types_elem_0(nodes_elem_0.size(), MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElementType> node_types_elem_1(nodes_elem_1.size(), MonolayerVertexElementType::Undetermined);
        basic_vertex_elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, nodes_elem_0, node_types_elem_0));
        basic_vertex_elements.push_back(new MonolayerVertexElement<2, 2>(1, MonolayerVertexElementType::Undetermined, nodes_elem_1, node_types_elem_1));

        // Make a vertex mesh
        MutableMonolayerVertexMesh<2, 2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);

        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetNumNodes(), 7u);

        TS_ASSERT_DELTA(basic_vertex_mesh.GetNode(2)->rGetLocation()[0], 1.5, 1e-3);
        TS_ASSERT_EQUALS(basic_vertex_mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);

        // Check that the nodes know which elements they are in
        std::set<unsigned> temp_list1;
        temp_list1.insert(0);

        // Nodes 1 and 4 are only in element 0
        TS_ASSERT_EQUALS(basic_nodes[1]->rGetContainingElementIndices(), temp_list1);
        TS_ASSERT_EQUALS(basic_nodes[4]->rGetContainingElementIndices(), temp_list1);

        // Node 2 is in elements 0 and 1
        temp_list1.insert(1u);
        TS_ASSERT_EQUALS(basic_nodes[2]->rGetContainingElementIndices(), temp_list1);

        // Node 5 is only in element 1
        std::set<unsigned> temp_list2;
        temp_list2.insert(1u);
        TS_ASSERT_EQUALS(basic_nodes[5]->rGetContainingElementIndices(), temp_list2);

        // Test Set and Get methods
        TS_ASSERT_DELTA(basic_vertex_mesh.GetCellRearrangementThreshold(), 0.01, 1e-4); // Default value
        TS_ASSERT_DELTA(basic_vertex_mesh.GetT2Threshold(), 0.001, 1e-4); // Default value

        basic_vertex_mesh.SetCellRearrangementThreshold(0.03);
        basic_vertex_mesh.SetT2Threshold(0.003);

        TS_ASSERT_DELTA(basic_vertex_mesh.GetCellRearrangementThreshold(), 0.03, 1e-4);
        TS_ASSERT_DELTA(basic_vertex_mesh.GetT2Threshold(), 0.003, 1e-4);

        // Coverage
        TS_ASSERT_EQUALS(basic_vertex_mesh.SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(basic_vertex_mesh.SolveBoundaryElementMapping(0), 0u);
    }

    void TestBasic3dMutableMonolayerVertexMesh()
    {
        MutableMonolayerVertexMesh<3, 3>* p_mesh = ConstructCubeAndPyramidMesh();

        // Set/get parameter values
        p_mesh->SetCellRearrangementThreshold(0.1);
        p_mesh->SetT2Threshold(0.01);
        p_mesh->SetCellRearrangementRatio(1.6);

        TS_ASSERT_DELTA(p_mesh->GetCellRearrangementThreshold(), 0.1, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetT2Threshold(), 0.01, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetCellRearrangementRatio(), 1.6, 1e-6);

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

        MonolayerVertexElement<3, 3>* p_element_1 = p_mesh->GetElement(1);
        TS_ASSERT_EQUALS(p_element_1->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(p_element_1->GetNumFaces(), 5u);

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

        // Coverage
        TS_ASSERT_EQUALS(p_mesh->SolveNodeMapping(0), 0u);
        TS_ASSERT_EQUALS(p_mesh->SolveElementMapping(0), 0u);
        TS_ASSERT_EQUALS(p_mesh->SolveBoundaryElementMapping(0), 0u);

        // Tidy up
        delete p_mesh;
    }

    void TestMeshConstructionFromMeshReader()
    {
        // Create mesh
        MonolayerVertexMeshReader<2, 2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");
        MutableMonolayerVertexMesh<2, 2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Test Get methods
        TS_ASSERT_DELTA(mesh.GetCellRearrangementThreshold(), 0.01, 1e-4); // Default value
        TS_ASSERT_DELTA(mesh.GetT2Threshold(), 0.001, 1e-4); // Default value

        // Check we have the right number of nodes and elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        // Check some node co-ordinates
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(0)->GetPoint()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[0], 1.5, 1e-6);
        TS_ASSERT_DELTA(mesh.GetNode(2)->GetPoint()[1], 1.0, 1e-6);

        // Check second element has the right nodes
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNodeGlobalIndex(2), 6u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetNode(1), mesh.GetNode(5));
    }

    void TestSetNode()
    {
        // Create mesh
        MonolayerVertexMeshReader<2, 2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");
        MutableMonolayerVertexMesh<2, 2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        ChastePoint<2> point = mesh.GetNode(3)->GetPoint();
        TS_ASSERT_DELTA(point[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(point[1], 2.0, 1e-6);

        // Nudge node
        point.SetCoordinate(0, 1.1);
        mesh.SetNode(3, point);

        ChastePoint<2> point2 = mesh.GetNode(3)->GetPoint();
        TS_ASSERT_DELTA(point2[0], 1.1, 1e-6);
        TS_ASSERT_DELTA(point2[1], 2.0, 1e-6);

        // Nudge node again
        point.SetCoordinate(1, 1.9);
        mesh.SetNode(3, point);

        ChastePoint<2> point3 = mesh.GetNode(3)->GetPoint();
        TS_ASSERT_DELTA(point3[0], 1.1, 1e-6);
        TS_ASSERT_DELTA(point3[1], 1.9, 1e-6);
    }

    void TestAddNodeAndReMesh()
    {
        // Create mesh
        MutableMonolayerVertexMesh<3, 3>* p_mesh = ConstructCubeAndPyramidMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);

        // Choose a node
        ChastePoint<3> point = p_mesh->GetNode(2)->GetPoint();
        TS_ASSERT_DELTA(point[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(point[1], 1.0, 1e-6);
        TS_ASSERT_DELTA(point[2], 0.0, 1e-6);

        // Create a new node close to this node
        point.SetCoordinate(0, 0.0);
        point.SetCoordinate(1, 1.1);
        point.SetCoordinate(2, 0.1);
        Node<3>* p_node = new Node<3>(p_mesh->GetNumNodes(), point);

        unsigned old_num_nodes = p_mesh->GetNumNodes();

        // Add this new node to the mesh
        unsigned new_index = p_mesh->AddNode(p_node);
        TS_ASSERT_EQUALS(new_index, old_num_nodes);

        // Remesh to update correspondences
        VertexElementMap map(p_mesh->GetNumElements());
        p_mesh->ReMesh(map);

        TS_ASSERT_EQUALS(map.Size(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(map.IsIdentityMap(), true);

        // Check that the mesh is updated
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);

        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[0], 0.0, 1e-7);
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[1], 1.1, 1e-7);
        TS_ASSERT_DELTA(p_mesh->GetNode(new_index)->rGetLocation()[2], 0.1, 1e-7);

        // Now test AddNode() when mDeletedNodeIndices is populated

        // Label node 3 as deleted
        p_mesh->mDeletedNodeIndices.push_back(2);

        // Create a new node close to this node
        ChastePoint<3> point2;
        point2.SetCoordinate(0, 0.0);
        point2.SetCoordinate(1, 0.9);
        point2.SetCoordinate(2, 0.0);
        Node<3>* p_node2 = new Node<3>(p_mesh->GetNumNodes(), point);

        // Add this new node to the mesh
        new_index = p_mesh->AddNode(p_node2);
        TS_ASSERT_EQUALS(new_index, 2u);

        delete p_mesh;
    }

    void TestAddElement()
    {
        // Make four nodes to assign to two elements
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.5, 1.0));
        nodes.push_back(new Node<2>(3, false, 1.0, 2.0));
        nodes.push_back(new Node<2>(4, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(5, false, 2.0, 0.0));
        nodes.push_back(new Node<2>(6, false, 2.0, 3.0));

        // Make two triangular elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[5] = { 0, 1, 2, 3, 4 };
        unsigned node_indices_elem_1[3] = { 2, 5, 6 };
        for (unsigned i = 0; i < 5; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            if (i < 3)
            {
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }
        }

        std::vector<MonolayerVertexElement<2, 2>*> elements;

        std::vector<MonolayerVertexElementType> node_types_elem_0(nodes_elem_0.size(), MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElementType> node_types_elem_1(nodes_elem_1.size(), MonolayerVertexElementType::Undetermined);
        MonolayerVertexElement<2, 2>* p_replaced_vertex_element = new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, nodes_elem_0, node_types_elem_0);
        elements.push_back(p_replaced_vertex_element);
        elements.push_back(new MonolayerVertexElement<2, 2>(1, MonolayerVertexElementType::Undetermined, nodes_elem_1, node_types_elem_1));

        // Make a vertex mesh
        MutableMonolayerVertexMesh<2, 2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        std::vector<Node<2>*> nodes_elem_2, nodes_elem_3;

        // Make two triangular elements out of these nodes
        nodes_elem_2.push_back(nodes[6]);
        nodes_elem_2.push_back(nodes[1]);
        nodes_elem_2.push_back(nodes[2]);

        nodes_elem_3.push_back(nodes[0]);
        nodes_elem_3.push_back(nodes[1]);
        nodes_elem_3.push_back(nodes[2]);
        nodes_elem_3.push_back(nodes[3]);

        // Add a new element to the mesh
        std::vector<MonolayerVertexElementType> node_types_elem_2(nodes_elem_2.size(), MonolayerVertexElementType::Undetermined);
        mesh.AddElement(new MonolayerVertexElement<2, 2>(2, MonolayerVertexElementType::Undetermined, nodes_elem_2, node_types_elem_2));

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);

        // Replace element 0 in the mesh
        std::vector<MonolayerVertexElementType> node_types_elem_3(nodes_elem_3.size(), MonolayerVertexElementType::Undetermined);
        mesh.AddElement(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, nodes_elem_3, node_types_elem_3));

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);

        // Tidy up
        delete p_replaced_vertex_element;
    }

    void TestDeletingNodes()
    {
        // Make a simple vertex mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.5, 1.0));
        nodes.push_back(new Node<2>(3, false, 1.0, 2.0));
        nodes.push_back(new Node<2>(4, false, 0.0, 1.0));
        nodes.push_back(new Node<2>(5, false, 2.0, 0.0));
        nodes.push_back(new Node<2>(6, false, 2.0, 3.0));

        // Make two triangular elements out of these nodes
        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;
        unsigned node_indices_elem_0[5] = { 0, 1, 2, 3, 4 };
        unsigned node_indices_elem_1[3] = { 2, 5, 6 };
        for (unsigned i = 0; i < 5; i++)
        {
            nodes_elem_0.push_back(nodes[node_indices_elem_0[i]]);
            if (i < 3)
            {
                nodes_elem_1.push_back(nodes[node_indices_elem_1[i]]);
            }
        }

        std::vector<MonolayerVertexElement<2, 2>*> elements;
        std::vector<MonolayerVertexElementType> node_types_elem_0(nodes_elem_0.size(), MonolayerVertexElementType::Undetermined);
        std::vector<MonolayerVertexElementType> node_types_elem_1(nodes_elem_1.size(), MonolayerVertexElementType::Undetermined);
        elements.push_back(new MonolayerVertexElement<2, 2>(0, MonolayerVertexElementType::Undetermined, nodes_elem_0, node_types_elem_0));
        elements.push_back(new MonolayerVertexElement<2, 2>(1, MonolayerVertexElementType::Undetermined, nodes_elem_1, node_types_elem_1));

        MutableMonolayerVertexMesh<2, 2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 7u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 2u);

        mesh.DeleteElementPriorToReMesh(0);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 3u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
    }

    void TestArchive2dMutableMonolayerVertexMesh()
    {
        // Set archiving location
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "mutable_vertex_2d.arch";
        ArchiveLocationInfo::SetMeshFilename("mutable_vertex_2d");
        // Create mesh
        MonolayerVertexMeshReader<2, 2> mesh_reader("mesh/test/data/TestVertexMesh/honeycomb_vertex_mesh_5_by_3");
        MutableMonolayerVertexMesh<2, 2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set member variables
        mesh.SetCellRearrangementThreshold(0.54);
        mesh.SetT2Threshold(0.012);
        mesh.SetCellRearrangementRatio(1.6);
        mesh.SetDistanceForT3SwapChecking(7.3);

        AbstractMesh<2, 2>* const p_mesh = &mesh;

        /*
         * You need the const above to stop a BOOST_STATIC_ASSERTION failure.
         * This is because the serialization library only allows you to save tracked
         * objects while the compiler considers them const, to prevent the objects
         * changing during the save, and so object tracking leading to wrong results.
         *
         * E.g. A is saved once via pointer, then changed, then saved again. The second
         * save notes that A was saved before, so doesn't write its data again, and the
         * change is lost.
         */

        // Create an output archive
        {
            TS_ASSERT_EQUALS((static_cast<MutableMonolayerVertexMesh<2, 2>*>(p_mesh))->GetNumNodes(), 46u);
            TS_ASSERT_EQUALS((static_cast<MutableMonolayerVertexMesh<2, 2>*>(p_mesh))->GetNumElements(), 15u);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();
            // We have to serialize via a pointer here, or the derived class information is lost
            (*p_arch) << p_mesh;
        }

        {
            // De-serialize and compare
            AbstractMesh<2, 2>* p_mesh2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_mesh2;

            MutableMonolayerVertexMesh<2, 2>* p_mesh_original = static_cast<MutableMonolayerVertexMesh<2, 2>*>(p_mesh2);
            MutableMonolayerVertexMesh<2, 2>* p_mesh_loaded = static_cast<MutableMonolayerVertexMesh<2, 2>*>(p_mesh);

            // Test member variables were archived correctly
            TS_ASSERT_DELTA(p_mesh_original->GetCellRearrangementThreshold(), 0.54, 1e-6);
            TS_ASSERT_DELTA(p_mesh_loaded->GetCellRearrangementThreshold(), 0.54, 1e-6);
            TS_ASSERT_DELTA(p_mesh_original->GetT2Threshold(), 0.012, 1e-6);
            TS_ASSERT_DELTA(p_mesh_loaded->GetT2Threshold(), 0.012, 1e-6);
            TS_ASSERT_DELTA(p_mesh_loaded->GetCellRearrangementRatio(), 1.6, 1e-6);
            TS_ASSERT_DELTA(p_mesh_original->GetDistanceForT3SwapChecking(), 7.3, 1e-6);
            TS_ASSERT_DELTA(p_mesh_loaded->GetDistanceForT3SwapChecking(), 7.3, 1e-6);

            // Compare the loaded mesh against the original
            TS_ASSERT_EQUALS(p_mesh_original->GetNumNodes(), p_mesh_loaded->GetNumNodes());
            for (unsigned node_index = 0; node_index < p_mesh_original->GetNumNodes(); node_index++)
            {
                Node<2>* p_node = p_mesh_original->GetNode(node_index);
                Node<2>* p_node2 = p_mesh_loaded->GetNode(node_index);

                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());
                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());

                for (unsigned dimension = 0; dimension < 2; dimension++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[dimension], p_node2->rGetLocation()[dimension], 1e-4);
                }
            }
            TS_ASSERT_EQUALS(p_mesh_original->GetNumElements(), p_mesh_loaded->GetNumElements());

            for (unsigned elem_index = 0; elem_index < p_mesh_original->GetNumElements(); elem_index++)
            {
                TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNumNodes(),
                                 p_mesh_loaded->GetElement(elem_index)->GetNumNodes());

                for (unsigned local_index = 0; local_index < p_mesh_original->GetElement(elem_index)->GetNumNodes(); local_index++)
                {
                    TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNodeGlobalIndex(local_index),
                                     p_mesh_loaded->GetElement(elem_index)->GetNodeGlobalIndex(local_index));
                }
            }
            // Tidy up
            delete p_mesh_loaded;
        }
    }

    void TestArchive3dMutableMonolayerVertexMesh()
    {
        // Set archiving location
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "mutable_vertex_3d.arch";
        ArchiveLocationInfo::SetMeshFilename("mutable_vertex_3d");

        // Create mesh
        MutableMonolayerVertexMesh<3, 3>* p_mutable_mesh = ConstructCubeAndPyramidMesh();

        // Set member variables
        p_mutable_mesh->SetCellRearrangementThreshold(0.54);
        p_mutable_mesh->SetT2Threshold(0.012);
        p_mutable_mesh->SetCellRearrangementRatio(1.6);

        AbstractMesh<3, 3>* const p_mesh = p_mutable_mesh;

        /*
         * You need the const above to stop a BOOST_STATIC_ASSERTION failure.
         * This is because the serialization library only allows you to save tracked
         * objects while the compiler considers them const, to prevent the objects
         * changing during the save, and so object tracking leading to wrong results.
         *
         * E.g. A is saved once via pointer, then changed, then saved again. The second
         * save notes that A was saved before, so doesn't write its data again, and the
         * change is lost.
         */

        // Create an output archive
        {
            TS_ASSERT_EQUALS((static_cast<MutableMonolayerVertexMesh<3, 3>*>(p_mesh))->GetNumNodes(), 9u);
            TS_ASSERT_EQUALS((static_cast<MutableMonolayerVertexMesh<3, 3>*>(p_mesh))->GetNumElements(), 2u);
            TS_ASSERT_EQUALS((static_cast<MutableMonolayerVertexMesh<3, 3>*>(p_mesh))->GetNumFaces(), 10u);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // We have to serialize via a pointer here, or the derived class information is lost
            (*p_arch) << p_mesh;
        }

        {
            // De-serialize and compare
            AbstractMesh<3, 3>* p_mesh2;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore from the archive
            (*p_arch) >> p_mesh2;

            MutableMonolayerVertexMesh<3, 3>* p_mesh_original = static_cast<MutableMonolayerVertexMesh<3, 3>*>(p_mesh);
            MutableMonolayerVertexMesh<3, 3>* p_mesh_loaded = static_cast<MutableMonolayerVertexMesh<3, 3>*>(p_mesh2);

            // Compare the loaded mesh against the original
            TS_ASSERT_EQUALS(p_mesh_original->GetNumNodes(), p_mesh_loaded->GetNumNodes());

            for (unsigned node_index = 0; node_index < p_mesh_original->GetNumNodes(); node_index++)
            {
                Node<3>* p_node = p_mesh_original->GetNode(node_index);
                Node<3>* p_node2 = p_mesh_loaded->GetNode(node_index);

                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());

                TS_ASSERT_EQUALS(p_node->IsBoundaryNode(), p_node2->IsBoundaryNode());

                for (unsigned dimension = 0; dimension < 3; dimension++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[dimension], p_node2->rGetLocation()[dimension], 1e-4);
                }
            }

            TS_ASSERT_EQUALS(p_mesh_original->GetNumElements(), p_mesh_loaded->GetNumElements());

            for (unsigned elem_index = 0; elem_index < p_mesh_original->GetNumElements(); elem_index++)
            {
                TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNumNodes(),
                                 p_mesh_loaded->GetElement(elem_index)->GetNumNodes());

                TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNumFaces(),
                                 p_mesh_loaded->GetElement(elem_index)->GetNumFaces());

                for (unsigned local_index = 0; local_index < p_mesh_original->GetElement(elem_index)->GetNumNodes(); local_index++)
                {
                    TS_ASSERT_EQUALS(p_mesh_original->GetElement(elem_index)->GetNodeGlobalIndex(local_index),
                                     p_mesh_loaded->GetElement(elem_index)->GetNodeGlobalIndex(local_index));
                }
            }

            TS_ASSERT_DELTA(p_mesh_loaded->GetCellRearrangementThreshold(), 0.54, 1e-6);
            TS_ASSERT_DELTA(p_mesh_loaded->GetT2Threshold(), 0.012, 1e-6);
            TS_ASSERT_DELTA(p_mesh_loaded->GetCellRearrangementRatio(), 1.6, 1e-6);

            TS_ASSERT_DELTA(p_mesh_original->GetCellRearrangementThreshold(), 0.54, 1e-6);
            TS_ASSERT_DELTA(p_mesh_original->GetT2Threshold(), 0.012, 1e-6);
            TS_ASSERT_DELTA(p_mesh_original->GetCellRearrangementRatio(), 1.6, 1e-6);

            // Tidy up
            delete p_mesh_loaded;
        }

        delete p_mesh;
    }
};

#endif /*TESTMUTABLEMONOLAYERVERTEXMESH_HPP_*/
