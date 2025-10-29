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

#ifndef TESTMUTABLEMONOLAYERVERTEXMESHREMESH_HPP_
#define TESTMUTABLEMONOLAYERVERTEXMESHREMESH_HPP_

#include <cxxtest/TestSuite.h>

#include "FileComparison.hpp"
#include "FiniteThicknessHoneycombVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshWriter.hpp"
#include "MutableMonolayerVertexMesh.hpp"
#include "Warnings.hpp"

#include <cmath> //for M_PI

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestMutableMonolayerVertexMeshReMesh : public CxxTest::TestSuite
{
public:
    void TestPerformNodeMerge() noexcept(false)
    {
        /*
         * Create a mesh comprising a single triangular finite thickness element, as shown below.
         * We will test that the nodes marked with an x are merged correctly.
         *				 /|
         *			  / |
         *      /|  |
         *     / |  |
         *    /  |  |
         *   /...|..|
         *  /: ::| /
         *  --xx- /
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.4, 0.0, 0.0));
        nodes.push_back(new Node<3>(4, true, 0.6, 0.0, 0.0));

        nodes.push_back(new Node<3>(5, true, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(6, true, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, true, 1.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(8, true, 0.4, 0.0, 1.0));
        nodes.push_back(new Node<3>(9, true, 0.6, 0.0, 1.0));

        std::vector<MonolayerVertexElementType> node_types;
        for (unsigned i = 0; i < 10; i++)
        {
            if (i < 5)
            {
                node_types.push_back(MonolayerVertexElementType::Basal);
            }
            else
            {
                node_types.push_back(MonolayerVertexElementType::Apical);
            }
        }

        std::vector<std::vector<unsigned> > node_indices_all_faces;
        std::vector<unsigned> node_indices_face0 = { 0, 2, 1, 4, 3 };
        std::vector<unsigned> node_indices_face1 = { 5, 8, 9, 6, 7 };
        std::vector<unsigned> node_indices_face2 = { 0, 5, 7, 2 };
        std::vector<unsigned> node_indices_face3 = { 2, 7, 6, 1 };
        std::vector<unsigned> node_indices_face4 = { 1, 6, 9, 4 };
        std::vector<unsigned> node_indices_face5 = { 4, 9, 8, 3 };
        std::vector<unsigned> node_indices_face6 = { 3, 8, 5, 0 };

        node_indices_all_faces.push_back(node_indices_face0);
        node_indices_all_faces.push_back(node_indices_face1);
        node_indices_all_faces.push_back(node_indices_face2);
        node_indices_all_faces.push_back(node_indices_face3);
        node_indices_all_faces.push_back(node_indices_face4);
        node_indices_all_faces.push_back(node_indices_face5);
        node_indices_all_faces.push_back(node_indices_face6);

        std::vector<MonolayerVertexElementType> face_types_all;
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);

        std::vector<MonolayerVertexElement<2, 3>*> face_elements;
        std::vector<bool> orientations;
        for (unsigned i = 0; i < 7; i++)
        {
            std::vector<Node<3>*> temp_nodes;
            std::vector<MonolayerVertexElementType> temp_types;
            unsigned max_index = i < 2 ? 5 : 4;
            for (unsigned j = 0; j < max_index; j++)
            {
                temp_nodes.push_back(nodes[node_indices_all_faces[i][j]]);
                temp_types.push_back(node_types[node_indices_all_faces[i][j]]);
            }
            face_elements.push_back(new MonolayerVertexElement<2, 3>(i, face_types_all[i], temp_nodes, temp_types));
            orientations.push_back(true);
            face_elements[i]->SetAsBoundaryFace();
        }

        std::vector<MonolayerVertexElement<3, 3>*> vertex_elements;
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, face_elements, orientations, nodes, node_types));

        // Test the element's area and perimeter are computed correctly
        MutableMonolayerVertexMesh<3, 3> vertex_mesh(nodes, face_elements, vertex_elements);

        // Test if mesh is correctly set up
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 7u);

        MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(0);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(p_element->GetFace(0)), 0.5, 1e-6);

        // Merge nodes 3 and 4
        MonolayerVertexElement<2, 3>* pFace = face_elements[5];
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(3), vertex_mesh.GetNode(4), pFace);

        // Test the mesh is correctly updated
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 6u);

        // Test the correct nodes are boundary nodes
        for (unsigned i = 0; i < 8; i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetFace(4)->GetNodeType(1), MonolayerVertexElementType::Apical);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetFace(4)->GetNodeType(2), MonolayerVertexElementType::Apical);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetFace(0)->GetNodeType(3), MonolayerVertexElementType::Basal);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetFace(1)->GetNodeType(1), MonolayerVertexElementType::Apical);

        // Test the merged nodes are in the correct place
        TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[0], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[1], 0.0, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(3)->rGetLocation()[2], 0.0, 1e-3);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(7)->rGetLocation()[0], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(7)->rGetLocation()[1], 0.0, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(7)->rGetLocation()[2], 1.0, 1e-3);

        // Test the elements own the correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 8u);
        pFace = vertex_mesh.GetElement(0)->GetFace(4);
        TS_ASSERT_EQUALS(pFace->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(pFace->GetNodeGlobalIndex(0), 1u);
        TS_ASSERT_EQUALS(pFace->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(pFace->GetNodeGlobalIndex(2), 7u);
        TS_ASSERT_EQUALS(pFace->GetNodeGlobalIndex(3), 3u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetFace(1)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetFace(1)->GetNodeGlobalIndex(0), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetFace(1)->GetNodeGlobalIndex(1), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetFace(1)->GetNodeGlobalIndex(2), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetFace(1)->GetNodeGlobalIndex(3), 6u);

        // Test the element's area and perimeter are computed correctly
        p_element = vertex_mesh.GetElement(0);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(p_element->GetFace(0)), 0.5, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(p_element->GetFace(1)), 0.5, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(p_element->GetFace(2)), sqrt(2.0), 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(p_element->GetFace(3)), 1.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(p_element->GetFace(4)), 0.5, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(p_element->GetFace(5)), 0.5, 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.5, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 3 + sqrt(2.0), 1e-6);
    }

    void TestPerformNodeMergeWhenLowIndexNodeMustBeAddedToElement()
    {
        /**
         * Create a mesh comprising two square elements, as shown below. We will test that the
         * nodes marked with an x are merged correctly. We will test node merging in the case
         * where, when the elements previously containing the high-index node are updated to
         * contain the low-index node, at least one of these elements did not already contain
         * the low-index node.
         *       basal
         *     6
         *   / .\
         *  |\ . \  5  4
         *  |  \  ' x_x___ 3
         *  |  . \ / /   /|
         *  |  . 'x-x---  |
         *  |  0..|.1...|.|2
         *  | :   |:    |/
         *   -----------
         *    apical
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.00, 0.00, 0.00));
        nodes.push_back(new Node<3>(1, true, 1.00, 0.00, 0.00));
        nodes.push_back(new Node<3>(2, true, 2.00, 0.00, 0.00));
        nodes.push_back(new Node<3>(3, true, 2.00, 1.00, 0.00));
        nodes.push_back(new Node<3>(4, true, 1.01, 1.00, 0.00));
        nodes.push_back(new Node<3>(5, true, 1.00, 1.00, 0.00));
        nodes.push_back(new Node<3>(6, true, 0.00, 2.00, 0.00));

        nodes.push_back(new Node<3>(7, true, 0.00, 0.00, 1.00));
        nodes.push_back(new Node<3>(8, true, 1.00, 0.00, 1.00));
        nodes.push_back(new Node<3>(9, true, 2.00, 0.00, 1.00));
        nodes.push_back(new Node<3>(10, true, 2.00, 1.00, 1.00));
        nodes.push_back(new Node<3>(11, true, 1.01, 1.00, 1.00));
        nodes.push_back(new Node<3>(12, true, 1.00, 1.00, 1.00));
        nodes.push_back(new Node<3>(13, true, 0.00, 2.00, 1.00));

        std::vector<MonolayerVertexElementType> node_types;
        for (unsigned i = 0; i < 14; i++)
        {
            if (i < 7)
            {
                node_types.push_back(MonolayerVertexElementType::Basal);
            }
            else
            {
                node_types.push_back(MonolayerVertexElementType::Apical);
            }
        }

        std::vector<std::vector<unsigned> > node_indices_all_faces;
        // basal
        std::vector<unsigned> node_indices_face0 = { 0, 6, 5, 1 };
        std::vector<unsigned> node_indices_face1 = { 1, 5, 4, 3, 2 };
        // apical
        std::vector<unsigned> node_indices_face2 = { 7, 8, 12, 13 };
        std::vector<unsigned> node_indices_face3 = { 8, 9, 10, 11, 12 };
        // lateral
        std::vector<unsigned> node_indices_face4 = { 0, 7, 13, 6 };
        std::vector<unsigned> node_indices_face5 = { 5, 6, 13, 12 };
        std::vector<unsigned> node_indices_face6 = { 4, 5, 12, 11 };
        std::vector<unsigned> node_indices_face7 = { 3, 4, 11, 10 };
        std::vector<unsigned> node_indices_face8 = { 3, 10, 9, 2 };
        std::vector<unsigned> node_indices_face9 = { 1, 2, 9, 8 };
        std::vector<unsigned> node_indices_face10 = { 0, 1, 8, 7 };
        // shared lateral, correctly oriented to left
        std::vector<unsigned> node_indices_face11 = { 1, 5, 12, 8 };

        node_indices_all_faces.push_back(node_indices_face0);
        node_indices_all_faces.push_back(node_indices_face1);
        node_indices_all_faces.push_back(node_indices_face2);
        node_indices_all_faces.push_back(node_indices_face3);
        node_indices_all_faces.push_back(node_indices_face4);
        node_indices_all_faces.push_back(node_indices_face5);
        node_indices_all_faces.push_back(node_indices_face6);
        node_indices_all_faces.push_back(node_indices_face7);
        node_indices_all_faces.push_back(node_indices_face8);
        node_indices_all_faces.push_back(node_indices_face9);
        node_indices_all_faces.push_back(node_indices_face10);
        node_indices_all_faces.push_back(node_indices_face11);

        std::vector<MonolayerVertexElementType> face_types_all;
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);

        std::vector<MonolayerVertexElement<2, 3>*> face_elements;
        for (unsigned i = 0; i < 12; i++)
        {
            std::vector<Node<3>*> temp_nodes;
            std::vector<MonolayerVertexElementType> temp_types;
            unsigned max_index = (i == 1 || i == 3) ? 5 : 4;
            for (unsigned j = 0; j < max_index; j++)
            {
                temp_nodes.push_back(nodes[node_indices_all_faces[i][j]]);
                temp_types.push_back(node_types[node_indices_all_faces[i][j]]);
            }
            face_elements.push_back(new MonolayerVertexElement<2, 3>(i, face_types_all[i], temp_nodes, temp_types));
            if (i != 11)
            {
                face_elements[i]->SetAsBoundaryFace();
            }
        }

        /*    .   Face numbering
         *   / .\ 0   basal
         *  |\ .  \     1
         *  | \.5  'x_x___
         *  |  .\  /6/ 7 /|
         *  |4 . 'x-x--- 8|
         *  |  ...|11...|.|
         *  | : 10|:  9 |/
         *   -----------
         *    2      3
         *    apical
         */

        std::vector<MonolayerVertexElement<2, 3>*> face_elements_0;
        face_elements_0.push_back(face_elements[0]);
        face_elements_0.push_back(face_elements[2]);
        face_elements_0.push_back(face_elements[4]);
        face_elements_0.push_back(face_elements[5]);
        face_elements_0.push_back(face_elements[11]);
        face_elements_0.push_back(face_elements[10]);
        std::vector<bool> orientations_0(6, false);

        std::vector<MonolayerVertexElement<2, 3>*> face_elements_1;
        face_elements_1.push_back(face_elements[1]);
        face_elements_1.push_back(face_elements[3]);
        face_elements_1.push_back(face_elements[6]);
        face_elements_1.push_back(face_elements[7]);
        face_elements_1.push_back(face_elements[8]);
        face_elements_1.push_back(face_elements[9]);
        face_elements_1.push_back(face_elements[11]);
        std::vector<bool> orientations_1(7, false);
        orientations_1[6] = true;

        std::vector<MonolayerVertexElement<3, 3>*> vertex_elements;
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, face_elements_0, orientations_0));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(1, MonolayerVertexElementType::Undetermined, face_elements_1, orientations_1));

        MutableMonolayerVertexMesh<3, 3> vertex_mesh(nodes, face_elements, vertex_elements);

        // Check if volumes and surface areas correct before merging
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 1.5, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 7.0 + sqrt(2), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 1.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 6.0, 1e-6);

        // Merge nodes 4 and 5
        TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), face_elements[11]), "The nodes given to merge do not belong to pFace");
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), face_elements[6]);

        // Test the mesh is correctly updated
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 12u);

        // Test the correct nodes are boundary nodes
        for (unsigned i = 0; i < vertex_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetNode(i)->IsBoundaryNode(), true);
        }

        // Test that the moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 1.005, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0, 1e-8);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 1.005, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.0, 1e-8);

        // Test the elements own the correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 8u);
        unsigned node_indices_element_0[8] = { 0, 1, 4, 5, 6, 7, 10, 11 };
        unsigned node_indices_element_1[8] = { 1, 2, 3, 4, 7, 8, 9, 10 };
        bool found_a = false;
        bool found_b = false;
        for (unsigned i = 0; i < 8; i++)
        {
            found_a = false;
            found_b = false;
            for (unsigned j = 0; j < 8; j++)
            {
                if (vertex_mesh.GetElement(0)->GetNodeGlobalIndex(i) == node_indices_element_0[j])
                {
                    found_a = true;
                }
                if (vertex_mesh.GetElement(1)->GetNodeGlobalIndex(i) == node_indices_element_1[j])
                {
                    found_b = true;
                }
            }
            if ((!found_a) || (!found_b))
            {
                break;
            }
        }
        TS_ASSERT(found_a);
        TS_ASSERT(found_b);

        // Test the elements own the correct faces
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 6u);
        unsigned node_indices_faces_0[6] = { 0, 2, 4, 5, 11, 10 };
        unsigned node_indices_faces_1[6] = { 1, 3, 7, 8, 9, 11 };
        found_a = false;
        found_b = false;
        for (unsigned i = 0; i < 6; i++)
        {
            found_a = false;
            found_b = false;
            for (unsigned j = 0; j < 6; j++)
            {
                if (vertex_mesh.GetElement(0)->GetFace(i) == face_elements[node_indices_faces_0[j]])
                {
                    found_a = true;
                }
                if (vertex_mesh.GetElement(1)->GetFace(i) == face_elements[node_indices_faces_1[j]])
                {
                    found_b = true;
                }
            }
            if ((!found_a) || (!found_b))
            {
                break;
            }
        }
        TS_ASSERT(found_a);
        TS_ASSERT(found_b);

        // Get index of inner face
        unsigned index_inner_face = face_elements[11]->GetIndex();

        // Check if boundary faces correct
        for (unsigned index_face = 0; index_face < 11; index_face++)
        {
            MonolayerVertexElement<2, 3>* p_face = vertex_mesh.GetFace(index_face);
            if (index_face == index_inner_face)
            {
                TS_ASSERT(!(p_face->IsBoundaryFace()));
            }
            else
            {
                TS_ASSERT(p_face->IsBoundaryFace());
            }
        }

        // Check if volumes and surface areas correct
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 1.0 + 0.0025 + (1.005) / 2.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 5.0 + 0.005 + 1.005 + sqrt(1.0 + (1.005 * 1.005)) + sqrt(1.0 + (0.005 * 0.005)), 1e-6);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.995 + 0.0025, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 2.0 + 0.9975 + 0.9975 + 0.995 + sqrt(1.0 + (0.005 * 0.005)), 1e-6);
    }

    void TestPerformT1SwapAndIdentifySwapType()
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * We will test that that a T1 swap of the two central nodes is correctly implemented.
         *   _____         basal network
         *  /     /|			  3_____2
         * /_____/ |				|\   /|
         * |\   /| |				| \5/ |
         * | \ / | |				|  |  |
         * |  |  | |basal		| /4\ |
         * | / \ | /				|/___\|
         * |/___\|/apical		0			1
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.6, 0.0));

        nodes.push_back(new Node<3>(6, true, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, true, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(8, true, 1.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(9, true, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(10, false, 0.5, 0.4, 1.0));
        nodes.push_back(new Node<3>(11, false, 0.5, 0.6, 1.0));

        std::vector<MonolayerVertexElementType> node_types;
        for (unsigned i = 0; i < 12; i++)
        {
            if (i < 6)
            {
                node_types.push_back(MonolayerVertexElementType::Basal);
            }
            else
            {
                node_types.push_back(MonolayerVertexElementType::Apical);
            }
        }

        std::vector<std::vector<unsigned> > node_indices_all_faces;
        // basal
        std::vector<unsigned> node_indices_face0 = { 0, 4, 1 };
        std::vector<unsigned> node_indices_face1 = { 0, 3, 5, 4 };
        std::vector<unsigned> node_indices_face2 = { 3, 2, 5 };
        std::vector<unsigned> node_indices_face3 = { 5, 2, 1, 4 };
        // apical
        std::vector<unsigned> node_indices_face4 = { 6, 7, 10 };
        std::vector<unsigned> node_indices_face5 = { 6, 10, 11, 9 };
        std::vector<unsigned> node_indices_face6 = { 9, 11, 8 };
        std::vector<unsigned> node_indices_face7 = { 11, 10, 7, 8 };
        // lateral outer
        std::vector<unsigned> node_indices_face8 = { 0, 1, 7, 6 };
        std::vector<unsigned> node_indices_face9 = { 3, 0, 6, 9 };
        std::vector<unsigned> node_indices_face10 = { 2, 3, 9, 8 };
        std::vector<unsigned> node_indices_face11 = { 1, 2, 8, 7 };
        // lateral inner, normals always point to the right
        std::vector<unsigned> node_indices_face12 = { 0, 4, 10, 6 };
        std::vector<unsigned> node_indices_face13 = { 5, 3, 9, 11 };
        std::vector<unsigned> node_indices_face14 = { 5, 2, 8, 11 };
        std::vector<unsigned> node_indices_face15 = { 1, 4, 10, 7 };
        std::vector<unsigned> node_indices_face16 = { 4, 5, 11, 10 };

        node_indices_all_faces.push_back(node_indices_face0);
        node_indices_all_faces.push_back(node_indices_face1);
        node_indices_all_faces.push_back(node_indices_face2);
        node_indices_all_faces.push_back(node_indices_face3);
        node_indices_all_faces.push_back(node_indices_face4);
        node_indices_all_faces.push_back(node_indices_face5);
        node_indices_all_faces.push_back(node_indices_face6);
        node_indices_all_faces.push_back(node_indices_face7);
        node_indices_all_faces.push_back(node_indices_face8);
        node_indices_all_faces.push_back(node_indices_face9);
        node_indices_all_faces.push_back(node_indices_face10);
        node_indices_all_faces.push_back(node_indices_face11);
        node_indices_all_faces.push_back(node_indices_face12);
        node_indices_all_faces.push_back(node_indices_face13);
        node_indices_all_faces.push_back(node_indices_face14);
        node_indices_all_faces.push_back(node_indices_face15);
        node_indices_all_faces.push_back(node_indices_face16);

        /*
         * Face numbering
         *     _____      basal network  apical network
         *    / 10  /|			   _____ 				 _____
         *   /_____/ |				|\ 2 /|				|\ 6 /|
         *   |\   /| |				| \ / |				| \ / |
         *   13\ /14 |11			|1 | 3|				|5 | 7|
         * 9 |  16 | |				| /0\ |				| /4\ |
         *   | / \15 /				|/___\|				|/___\|
         *   |/12_\|/
         *			8
         */

        std::vector<MonolayerVertexElementType> face_types_all;
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);

        std::vector<MonolayerVertexElement<2, 3>*> face_elements;
        for (unsigned i = 0; i < 17; i++)
        {
            std::vector<Node<3>*> temp_nodes;
            std::vector<MonolayerVertexElementType> temp_types;
            unsigned max_index = (i % 2 == 0 && i < 7) ? 3 : 4;
            for (unsigned j = 0; j < max_index; j++)
            {
                temp_nodes.push_back(nodes[node_indices_all_faces[i][j]]);
                temp_types.push_back(node_types[node_indices_all_faces[i][j]]);
            }
            face_elements.push_back(new MonolayerVertexElement<2, 3>(i, face_types_all[i], temp_nodes, temp_types));
            if (i > 7 && i < 12)
            {
                face_elements[i]->SetAsBoundaryFace();
            }
        }

        unsigned face_indices_0[5] = { 0, 8, 12, 15, 4 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_0;
        for (unsigned index = 0; index < 5; index++)
        {
            face_elements_0.push_back(face_elements[face_indices_0[index]]);
        }
        std::vector<bool> orientations_0 = { false, false, true, false, false };

        unsigned face_indices_1[6] = { 1, 12, 13, 16, 9, 5 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_1;
        for (unsigned index = 0; index < 6; index++)
        {
            face_elements_1.push_back(face_elements[face_indices_1[index]]);
        }
        std::vector<bool> orientations_1 = { false, false, false, false, false, false };

        unsigned face_indices_2[5] = { 2, 10, 13, 14, 6 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_2;
        for (unsigned index = 0; index < 5; index++)
        {
            face_elements_2.push_back(face_elements[face_indices_2[index]]);
        }
        std::vector<bool> orientations_2 = { false, false, true, false, false };

        unsigned face_indices_3[6] = { 3, 11, 14, 15, 16, 7 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_3;
        for (unsigned index = 0; index < 6; index++)
        {
            face_elements_3.push_back(face_elements[face_indices_3[index]]);
        }
        std::vector<bool> orientations_3 = { false, false, true, true, true, false };

        unsigned* all_face_indices[4] = { face_indices_0, face_indices_1, face_indices_2, face_indices_3 };

        std::vector<MonolayerVertexElement<3, 3>*> vertex_elements;
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, face_elements_0, orientations_0));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(1, MonolayerVertexElementType::Undetermined, face_elements_1, orientations_1));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(2, MonolayerVertexElementType::Undetermined, face_elements_2, orientations_2));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(3, MonolayerVertexElementType::Undetermined, face_elements_3, orientations_3));

        MutableMonolayerVertexMesh<3, 3> vertex_mesh(nodes, face_elements, vertex_elements);

        // Test that the orientations are correct prior to remeshing
        // Consistency check to see if mesh is correctly defined
        for (unsigned i = 0; i < 4; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1 * 2.0 / 1.5);

        // Perform a T1 swap on nodes 4 and 5
        TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), face_elements[11]), "The basal edge of a lateral face seems to have more than two nodes, has wrong nodes or is not oriented correctly");
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), face_elements[16]);

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0.0, 1e-3);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1.0, 1e-3);

        // Test that each element contains the correct faces following the rearrangement
        bool found_all_faces = true;
        for (unsigned i = 0; i < 4; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            unsigned num_faces = (i == 0 || i == 2) ? 5 : 6;
            for (unsigned j_face = 0; j_face < num_faces; j_face++)
            {
                unsigned global_face_index = *(all_face_indices[i] + j_face);
                MonolayerVertexElement<2, 3>* p_face_expected = face_elements[global_face_index];
                if (global_face_index == 16)
                    continue;
                if (p_element->GetFaceLocalIndex(p_face_expected) == UINT_MAX)
                {
                    found_all_faces = false;
                    break;
                }
            }
            TS_ASSERT(found_all_faces);
        }
        // Check that central face changed elements correctly.
        bool swapped_face_in_correct_elements = vertex_elements[0]->GetFaceLocalIndex(face_elements[16]) != UINT_MAX
            && vertex_elements[1]->GetFaceLocalIndex(face_elements[16]) == UINT_MAX
            && vertex_elements[2]->GetFaceLocalIndex(face_elements[16]) != UINT_MAX
            && vertex_elements[3]->GetFaceLocalIndex(face_elements[16]) == UINT_MAX;

        TS_ASSERT(swapped_face_in_correct_elements);

        // Check that swapped nodes are in correct elements

        bool test = vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Basal)->GetNodeLocalIndex(4) == UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Apical)->GetNodeLocalIndex(10) == UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Apical)->GetNodeLocalIndex(11) != UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Basal)->GetNodeLocalIndex(5) != UINT_MAX;
        TS_ASSERT(test);

        bool swapped_nodes_in_correct_elements_0 = vertex_elements[0]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_1 = vertex_elements[1]->GetNodeLocalIndex(4) == UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(10) == UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_2 = vertex_elements[2]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_3 = vertex_elements[3]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[3]->GetNodeLocalIndex(5) == UINT_MAX
            && vertex_elements[3]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[3]->GetNodeLocalIndex(11) == UINT_MAX;

        TS_ASSERT(swapped_nodes_in_correct_elements_0);
        TS_ASSERT(swapped_nodes_in_correct_elements_1);
        TS_ASSERT(swapped_nodes_in_correct_elements_2);
        TS_ASSERT(swapped_nodes_in_correct_elements_3);

        TS_ASSERT_EQUALS(vertex_elements[0]->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_elements[1]->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[2]->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_elements[3]->GetNumNodes(), 6u);

        TS_ASSERT_EQUALS(vertex_elements[0]->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[1]->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_elements[2]->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[3]->GetNumFaces(), 5u);

        // Test that the orientations are correct
        for (unsigned i = 0; i < 4; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(3), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.2 + 0.2 * sqrt(41.0) + 0.6, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.0 + 0.2 * sqrt(41.0) + 0.4, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.2 + 0.2 * sqrt(41.0) + 0.6, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(3), 1.0 + 0.2 * sqrt(41.0) + 0.4, 1e-6);

        // Test boundary faces
        for (unsigned face_index = 0; face_index < face_elements.size(); face_index++)
        {
            if (face_index > 7 && face_index < 12)
            {
                TS_ASSERT(face_elements[face_index]->IsBoundaryFace());
            }
            else
            {
                TS_ASSERT(!(face_elements[face_index]->IsBoundaryFace()));
            }
        }

        // Test T1 swap location tracking
        std::vector<c_vector<double, 3> > t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 1u);
        TS_ASSERT_DELTA(t1_locations[0][0], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][1], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][2], 0.5, 1e-6);

        // Test T1 swap location clearing
        vertex_mesh.ClearLocationsOfT1Swaps();
        t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 0u);
    }

    void TestPerformT1SwapAndIdentifySwapTypeVaryingThickness()
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * We will test that that a T1 swap of the two central nodes is correctly implemented.
         * Here we assume a thickness profile and check whether the swaps remain in-plane
         *   _____
         *   /     /|        basal network   profile
         *  /     / | 			  	3_____2      --__
         * /_____/  |					 |\   /|      |		--__
         * |\   /|  .						| \5/ |      |       |
         * | \ / | .						|  |  |      |     __|
         * |  |  | : basal			| /4\ |      | __--
         * | / \ | |			  		|/___\|      --
         * |/___\|/apical		  	0			1
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, -1.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, -1.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.4, -0.5));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.6, -0.5));

        nodes.push_back(new Node<3>(6, true, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, true, 1.0, 0.0, 2.0));
        nodes.push_back(new Node<3>(8, true, 1.0, 1.0, 2.0));
        nodes.push_back(new Node<3>(9, true, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(10, false, 0.5, 0.4, 1.5));
        nodes.push_back(new Node<3>(11, false, 0.5, 0.6, 1.5));

        std::vector<MonolayerVertexElementType> node_types;
        for (unsigned i = 0; i < 12; i++)
        {
            if (i < 6)
            {
                node_types.push_back(MonolayerVertexElementType::Basal);
            }
            else
            {
                node_types.push_back(MonolayerVertexElementType::Apical);
            }
        }

        std::vector<std::vector<unsigned> > node_indices_all_faces;
        // basal
        std::vector<unsigned> node_indices_face0 = { 0, 4, 1 };
        std::vector<unsigned> node_indices_face1 = { 0, 3, 5, 4 };
        std::vector<unsigned> node_indices_face2 = { 3, 2, 5 };
        std::vector<unsigned> node_indices_face3 = { 5, 2, 1, 4 };
        // apical
        std::vector<unsigned> node_indices_face4 = { 6, 7, 10 };
        std::vector<unsigned> node_indices_face5 = { 6, 10, 11, 9 };
        std::vector<unsigned> node_indices_face6 = { 9, 11, 8 };
        std::vector<unsigned> node_indices_face7 = { 11, 10, 7, 8 };
        // lateral outer
        std::vector<unsigned> node_indices_face8 = { 0, 1, 7, 6 };
        std::vector<unsigned> node_indices_face9 = { 3, 0, 6, 9 };
        std::vector<unsigned> node_indices_face10 = { 2, 3, 9, 8 };
        std::vector<unsigned> node_indices_face11 = { 1, 2, 8, 7 };
        // lateral inner, normals always point to the right
        std::vector<unsigned> node_indices_face12 = { 0, 4, 10, 6 };
        std::vector<unsigned> node_indices_face13 = { 5, 3, 9, 11 };
        std::vector<unsigned> node_indices_face14 = { 5, 2, 8, 11 };
        std::vector<unsigned> node_indices_face15 = { 1, 4, 10, 7 };
        std::vector<unsigned> node_indices_face16 = { 4, 5, 11, 10 };

        node_indices_all_faces.push_back(node_indices_face0);
        node_indices_all_faces.push_back(node_indices_face1);
        node_indices_all_faces.push_back(node_indices_face2);
        node_indices_all_faces.push_back(node_indices_face3);
        node_indices_all_faces.push_back(node_indices_face4);
        node_indices_all_faces.push_back(node_indices_face5);
        node_indices_all_faces.push_back(node_indices_face6);
        node_indices_all_faces.push_back(node_indices_face7);
        node_indices_all_faces.push_back(node_indices_face8);
        node_indices_all_faces.push_back(node_indices_face9);
        node_indices_all_faces.push_back(node_indices_face10);
        node_indices_all_faces.push_back(node_indices_face11);
        node_indices_all_faces.push_back(node_indices_face12);
        node_indices_all_faces.push_back(node_indices_face13);
        node_indices_all_faces.push_back(node_indices_face14);
        node_indices_all_faces.push_back(node_indices_face15);
        node_indices_all_faces.push_back(node_indices_face16);

        /*
         * Face numbering
         *     _____      basal network  apical network
         *    / 10  /|			   _____ 				 _____
         *   /_____/ |				|\ 2 /|				|\ 6 /|
         *   |\   /| |				| \ / |				| \ / |
         *   13\ /14 |11			|1 | 3|				|5 | 7|
         * 9 |  16 | |				| /0\ |				| /4\ |
         *   | / \15 /				|/___\|				|/___\|
         *   |/12_\|/
         *			8
         */

        std::vector<MonolayerVertexElementType> face_types_all;
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);

        std::vector<MonolayerVertexElement<2, 3>*> face_elements;
        for (unsigned i = 0; i < 17; i++)
        {
            std::vector<Node<3>*> temp_nodes;
            std::vector<MonolayerVertexElementType> temp_types;
            unsigned max_index = (i % 2 == 0 && i < 7) ? 3 : 4;
            for (unsigned j = 0; j < max_index; j++)
            {
                temp_nodes.push_back(nodes[node_indices_all_faces[i][j]]);
                temp_types.push_back(node_types[node_indices_all_faces[i][j]]);
            }
            face_elements.push_back(new MonolayerVertexElement<2, 3>(i, face_types_all[i], temp_nodes, temp_types));
            if (i > 7 && i < 12)
            {
                face_elements[i]->SetAsBoundaryFace();
            }
        }

        unsigned face_indices_0[5] = { 0, 8, 12, 15, 4 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_0;
        for (unsigned index = 0; index < 5; index++)
        {
            face_elements_0.push_back(face_elements[face_indices_0[index]]);
        }
        std::vector<bool> orientations_0 = { false, false, true, false, false };

        unsigned face_indices_1[6] = { 1, 12, 13, 16, 9, 5 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_1;
        for (unsigned index = 0; index < 6; index++)
        {
            face_elements_1.push_back(face_elements[face_indices_1[index]]);
        }
        std::vector<bool> orientations_1 = { false, false, false, false, false, false };

        unsigned face_indices_2[5] = { 2, 10, 13, 14, 6 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_2;
        for (unsigned index = 0; index < 5; index++)
        {
            face_elements_2.push_back(face_elements[face_indices_2[index]]);
        }
        std::vector<bool> orientations_2 = { false, false, true, false, false };

        unsigned face_indices_3[6] = { 3, 11, 14, 15, 16, 7 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_3;
        for (unsigned index = 0; index < 6; index++)
        {
            face_elements_3.push_back(face_elements[face_indices_3[index]]);
        }
        std::vector<bool> orientations_3 = { false, false, true, true, true, false };

        unsigned* all_face_indices[4] = { face_indices_0, face_indices_1, face_indices_2, face_indices_3 };

        std::vector<MonolayerVertexElement<3, 3>*> vertex_elements;
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, face_elements_0, orientations_0));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(1, MonolayerVertexElementType::Undetermined, face_elements_1, orientations_1));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(2, MonolayerVertexElementType::Undetermined, face_elements_2, orientations_2));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(3, MonolayerVertexElementType::Undetermined, face_elements_3, orientations_3));

        MutableMonolayerVertexMesh<3, 3> vertex_mesh(nodes, face_elements, vertex_elements);

        // Test that the orientations are correct prior to remeshing
        // Consistency check to see if mesh is correctly defined
        for (unsigned i = 0; i < 4; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1 * 2.0 / 1.5 * sqrt(2.0));

        // Perform a T1 swap on nodes 4 and 5
        TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), face_elements[11]), "The basal edge of a lateral face seems to have more than two nodes, has wrong nodes or is not oriented correctly");
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), face_elements[16]);

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], -0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], -0.4, 1e-3);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1.4, 1e-3);

        // Test that each element contains the correct faces following the rearrangement
        bool found_all_faces = true;
        for (unsigned i = 0; i < 4; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            unsigned num_faces = (i == 0 || i == 2) ? 5 : 6;
            for (unsigned j_face = 0; j_face < num_faces; j_face++)
            {
                unsigned global_face_index = *(all_face_indices[i] + j_face);
                MonolayerVertexElement<2, 3>* p_face_expected = face_elements[global_face_index];
                if (global_face_index == 16)
                    continue;
                if (p_element->GetFaceLocalIndex(p_face_expected) == UINT_MAX)
                {
                    found_all_faces = false;
                    break;
                }
            }
            TS_ASSERT(found_all_faces);
        }
        // Check that central face changed elements correctly.
        bool swapped_face_in_correct_elements = vertex_elements[0]->GetFaceLocalIndex(face_elements[16]) != UINT_MAX
            && vertex_elements[1]->GetFaceLocalIndex(face_elements[16]) == UINT_MAX
            && vertex_elements[2]->GetFaceLocalIndex(face_elements[16]) != UINT_MAX
            && vertex_elements[3]->GetFaceLocalIndex(face_elements[16]) == UINT_MAX;

        TS_ASSERT(swapped_face_in_correct_elements);

        // Check that swapped nodes are in correct elements

        bool test = vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Basal)->GetNodeLocalIndex(4) == UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Apical)->GetNodeLocalIndex(10) == UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Apical)->GetNodeLocalIndex(11) != UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Basal)->GetNodeLocalIndex(5) != UINT_MAX;
        TS_ASSERT(test);

        bool swapped_nodes_in_correct_elements_0 = vertex_elements[0]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_1 = vertex_elements[1]->GetNodeLocalIndex(4) == UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(10) == UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_2 = vertex_elements[2]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_3 = vertex_elements[3]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[3]->GetNodeLocalIndex(5) == UINT_MAX
            && vertex_elements[3]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[3]->GetNodeLocalIndex(11) == UINT_MAX;

        TS_ASSERT(swapped_nodes_in_correct_elements_0);
        TS_ASSERT(swapped_nodes_in_correct_elements_1);
        TS_ASSERT(swapped_nodes_in_correct_elements_2);
        TS_ASSERT(swapped_nodes_in_correct_elements_3);

        TS_ASSERT_EQUALS(vertex_elements[0]->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_elements[1]->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[2]->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_elements[3]->GetNumNodes(), 6u);

        TS_ASSERT_EQUALS(vertex_elements[0]->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[1]->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_elements[2]->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[3]->GetNumFaces(), 5u);

        // Test that the orientations are correct
        for (unsigned i = 0; i < 4; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Test T1 swap location tracking
        std::vector<c_vector<double, 3> > t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 1u);
        TS_ASSERT_DELTA(t1_locations[0][0], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][1], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][2], 0.5, 1e-6);

        // Test T1 swap location clearing
        vertex_mesh.ClearLocationsOfT1Swaps();
        t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 0u);
    }

    void TestPerformT1SwapOnBoundary()
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * One of the rhomboids is a void. We will test that that a T1 swap of the two central nodes is correctly implemented.
         *   _____         basal network
         *  /     /			  3_____2
         * /_____/ 				|\   /
         * |\   / 				| \5/
         * | \ //\ 				|  |
         * |  |/  \ basal | /4\
         * | / \  /				|/___\
         * |/___\/apical	0			1
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, true, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, true, 0.5, 0.6, 0.0));

        nodes.push_back(new Node<3>(6, true, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, true, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(8, true, 1.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(9, true, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(10, true, 0.5, 0.4, 1.0));
        nodes.push_back(new Node<3>(11, true, 0.5, 0.6, 1.0));

        std::vector<MonolayerVertexElementType> node_types;
        for (unsigned i = 0; i < 12; i++)
        {
            if (i < 6)
            {
                node_types.push_back(MonolayerVertexElementType::Basal);
            }
            else
            {
                node_types.push_back(MonolayerVertexElementType::Apical);
            }
        }

        std::vector<std::vector<unsigned> > node_indices_all_faces;
        // basal
        std::vector<unsigned> node_indices_face0 = { 0, 4, 1 };
        std::vector<unsigned> node_indices_face1 = { 0, 3, 5, 4 };
        std::vector<unsigned> node_indices_face2 = { 3, 2, 5 };
        // apical
        std::vector<unsigned> node_indices_face3 = { 6, 7, 10 };
        std::vector<unsigned> node_indices_face4 = { 6, 10, 11, 9 };
        std::vector<unsigned> node_indices_face5 = { 9, 11, 8 };
        // lateral outer
        std::vector<unsigned> node_indices_face6 = { 0, 1, 7, 6 };
        std::vector<unsigned> node_indices_face7 = { 3, 0, 6, 9 };
        std::vector<unsigned> node_indices_face8 = { 2, 3, 9, 8 };
        // lateral inner, normals always point to the right
        std::vector<unsigned> node_indices_face9 = { 0, 4, 10, 6 };
        std::vector<unsigned> node_indices_face10 = { 5, 3, 9, 11 };
        std::vector<unsigned> node_indices_face11 = { 5, 2, 8, 11 };
        std::vector<unsigned> node_indices_face12 = { 1, 4, 10, 7 };
        std::vector<unsigned> node_indices_face13 = { 4, 5, 11, 10 };

        node_indices_all_faces.push_back(node_indices_face0);
        node_indices_all_faces.push_back(node_indices_face1);
        node_indices_all_faces.push_back(node_indices_face2);
        node_indices_all_faces.push_back(node_indices_face3);
        node_indices_all_faces.push_back(node_indices_face4);
        node_indices_all_faces.push_back(node_indices_face5);
        node_indices_all_faces.push_back(node_indices_face6);
        node_indices_all_faces.push_back(node_indices_face7);
        node_indices_all_faces.push_back(node_indices_face8);
        node_indices_all_faces.push_back(node_indices_face9);
        node_indices_all_faces.push_back(node_indices_face10);
        node_indices_all_faces.push_back(node_indices_face11);
        node_indices_all_faces.push_back(node_indices_face12);
        node_indices_all_faces.push_back(node_indices_face13);

        /*
         * Face numbering
         *     _____      basal network  apical network
         *    / 8   /			   _____ 				 _____
         *   /_____/11 			|\ 2 /				|\ 5 /
         *   |\   /  				| \ / 				| \ /
         *   10\ //\  		  	|1 | 					|4 |
         * 7 |  13  \				| /0\ 				| /3\
         *   | / \12/				|/___\				|/___\
         *   |/9__\/
         *			6
         */

        std::vector<MonolayerVertexElementType> face_types_all;
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);

        std::vector<MonolayerVertexElement<2, 3>*> face_elements;
        for (unsigned i = 0; i < 14; i++)
        {
            std::vector<Node<3>*> temp_nodes;
            std::vector<MonolayerVertexElementType> temp_types;
            unsigned max_index = (i == 0 || i == 2 || i == 3 || i == 5) ? 3 : 4;
            for (unsigned j = 0; j < max_index; j++)
            {
                temp_nodes.push_back(nodes[node_indices_all_faces[i][j]]);
                temp_types.push_back(node_types[node_indices_all_faces[i][j]]);
            }
            face_elements.push_back(new MonolayerVertexElement<2, 3>(i, face_types_all[i], temp_nodes, temp_types));
            if ((i > 5 && i < 8) || (i > 10))
            {
                face_elements[i]->SetAsBoundaryFace();
            }
        }

        unsigned face_indices_0[5] = { 0, 6, 9, 12, 3 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_0;
        for (unsigned index = 0; index < 5; index++)
        {
            face_elements_0.push_back(face_elements[face_indices_0[index]]);
        }
        std::vector<bool> orientations_0 = { false, false, true, false, false };

        unsigned face_indices_1[6] = { 1, 9, 10, 13, 7, 4 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_1;
        for (unsigned index = 0; index < 6; index++)
        {
            face_elements_1.push_back(face_elements[face_indices_1[index]]);
        }
        std::vector<bool> orientations_1 = { false, false, false, false, false, false };

        unsigned face_indices_2[5] = { 2, 8, 10, 11, 5 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_2;
        for (unsigned index = 0; index < 5; index++)
        {
            face_elements_2.push_back(face_elements[face_indices_2[index]]);
        }
        std::vector<bool> orientations_2 = { false, false, true, false, false };

        unsigned* all_face_indices[3] = { face_indices_0, face_indices_1, face_indices_2 };

        std::vector<MonolayerVertexElement<3, 3>*> vertex_elements;
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, face_elements_0, orientations_0));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(1, MonolayerVertexElementType::Undetermined, face_elements_1, orientations_1));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(2, MonolayerVertexElementType::Undetermined, face_elements_2, orientations_2));

        MutableMonolayerVertexMesh<3, 3> vertex_mesh(nodes, face_elements, vertex_elements);

        // Test that the orientations are correct prior to remeshing
        // Consistency check to see if mesh is correctly defined
        for (unsigned i = 0; i < 3; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1 * 2.0 / 1.5);

        // Perform a T1 swap on nodes 4 and 5
        // TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), face_elements[9]), "The basal edge of a lateral face seems to have more than two nodes, has wrong nodes or is not oriented correctly");
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4), face_elements[13]);

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0.0, 1e-3);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1.0, 1e-3);

        // Test that each element contains the correct faces following the rearrangement
        bool found_all_faces = true;
        for (unsigned i = 0; i < 3; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            unsigned num_faces = (i == 0 || i == 2) ? 5 : 6;
            for (unsigned j_face = 0; j_face < num_faces; j_face++)
            {
                unsigned global_face_index = *(all_face_indices[i] + j_face);
                MonolayerVertexElement<2, 3>* p_face_expected = face_elements[global_face_index];
                if (global_face_index == 13)
                    continue;
                if (p_element->GetFaceLocalIndex(p_face_expected) == UINT_MAX)
                {
                    found_all_faces = false;
                    break;
                }
            }
            TS_ASSERT(found_all_faces);
        }
        // Check that central face changed elements correctly.
        bool swapped_face_in_correct_elements = vertex_elements[0]->GetFaceLocalIndex(face_elements[13]) != UINT_MAX
            && vertex_elements[1]->GetFaceLocalIndex(face_elements[13]) == UINT_MAX
            && vertex_elements[2]->GetFaceLocalIndex(face_elements[13]) != UINT_MAX;

        TS_ASSERT(swapped_face_in_correct_elements);

        // Check that swapped nodes are in correct elements

        bool test = vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Basal)->GetNodeLocalIndex(4) == UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Apical)->GetNodeLocalIndex(10) == UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Apical)->GetNodeLocalIndex(11) != UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Basal)->GetNodeLocalIndex(5) != UINT_MAX;
        TS_ASSERT(test);

        bool swapped_nodes_in_correct_elements_0 = vertex_elements[0]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_1 = vertex_elements[1]->GetNodeLocalIndex(4) == UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(10) == UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_2 = vertex_elements[2]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(11) != UINT_MAX;

        TS_ASSERT(swapped_nodes_in_correct_elements_0);
        TS_ASSERT(swapped_nodes_in_correct_elements_1);
        TS_ASSERT(swapped_nodes_in_correct_elements_2);

        TS_ASSERT_EQUALS(vertex_elements[0]->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_elements[1]->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[2]->GetNumNodes(), 8u);

        TS_ASSERT_EQUALS(vertex_elements[0]->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[1]->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_elements[2]->GetNumFaces(), 6u);

        // Test that the orientations are correct
        for (unsigned i = 0; i < 3; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.2 + 0.2 * sqrt(41.0) + 0.6, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.0 + 0.2 * sqrt(41.0) + 0.4, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.2 + 0.2 * sqrt(41.0) + 0.6, 1e-6);

        // Test boundary faces
        for (unsigned face_index = 0; face_index < face_elements.size(); face_index++)
        {
            if ((face_index > 5 && face_index < 8) || (face_index > 10 && face_index < 13))
            {
                TS_ASSERT(face_elements[face_index]->IsBoundaryFace());
            }
            else
            {
                TS_ASSERT(!(face_elements[face_index]->IsBoundaryFace()));
            }
        }

        // Test boundary nodes
        for (unsigned node_index = 0; node_index < nodes.size(); node_index++)
        {
            bool boundary_node = nodes[node_index]->IsBoundaryNode();
            if (node_index == 5 || node_index == 11)
            {
                TS_ASSERT(!boundary_node);
            }
            else
            {
                TS_ASSERT(boundary_node);
            }
        }

        // Test T1 swap location tracking
        std::vector<c_vector<double, 3> > t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 1u);
        TS_ASSERT_DELTA(t1_locations[0][0], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][1], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][2], 0.5, 1e-6);

        // Test T1 swap location clearing
        vertex_mesh.ClearLocationsOfT1Swaps();
        t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 0u);
    }

    void TestPerformT1SwapOnBoundary2()
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * One of the rhomboids is a void. We will test that that a T1 swap of the two central nodes is correctly implemented.
         *   _____         basal network
         *  /     /			  3_____2
         * /_____/ 				|\   /
         * |\   / 				| \5/
         * | \ //\ 				|  |
         * |  |/  \ basal | /4\
         * | / \  /				|/___\
         * |/___\/apical	0			1
         */

        // Same test as before with differently oriented face 13, to check if the mirrore situation works.
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, true, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, true, 0.5, 0.6, 0.0));

        nodes.push_back(new Node<3>(6, true, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, true, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(8, true, 1.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(9, true, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(10, true, 0.5, 0.4, 1.0));
        nodes.push_back(new Node<3>(11, true, 0.5, 0.6, 1.0));

        std::vector<MonolayerVertexElementType> node_types;
        for (unsigned i = 0; i < 12; i++)
        {
            if (i < 6)
            {
                node_types.push_back(MonolayerVertexElementType::Basal);
            }
            else
            {
                node_types.push_back(MonolayerVertexElementType::Apical);
            }
        }

        std::vector<std::vector<unsigned> > node_indices_all_faces;
        // basal
        std::vector<unsigned> node_indices_face0 = { 0, 4, 1 };
        std::vector<unsigned> node_indices_face1 = { 0, 3, 5, 4 };
        std::vector<unsigned> node_indices_face2 = { 3, 2, 5 };
        // apical
        std::vector<unsigned> node_indices_face3 = { 6, 7, 10 };
        std::vector<unsigned> node_indices_face4 = { 6, 10, 11, 9 };
        std::vector<unsigned> node_indices_face5 = { 9, 11, 8 };
        // lateral outer
        std::vector<unsigned> node_indices_face6 = { 0, 1, 7, 6 };
        std::vector<unsigned> node_indices_face7 = { 3, 0, 6, 9 };
        std::vector<unsigned> node_indices_face8 = { 2, 3, 9, 8 };
        // lateral inner, normals always point to the right
        std::vector<unsigned> node_indices_face9 = { 0, 4, 10, 6 };
        std::vector<unsigned> node_indices_face10 = { 5, 3, 9, 11 };
        std::vector<unsigned> node_indices_face11 = { 5, 2, 8, 11 };
        std::vector<unsigned> node_indices_face12 = { 1, 4, 10, 7 };
        std::vector<unsigned> node_indices_face13 = { 4, 10, 11, 5 };

        node_indices_all_faces.push_back(node_indices_face0);
        node_indices_all_faces.push_back(node_indices_face1);
        node_indices_all_faces.push_back(node_indices_face2);
        node_indices_all_faces.push_back(node_indices_face3);
        node_indices_all_faces.push_back(node_indices_face4);
        node_indices_all_faces.push_back(node_indices_face5);
        node_indices_all_faces.push_back(node_indices_face6);
        node_indices_all_faces.push_back(node_indices_face7);
        node_indices_all_faces.push_back(node_indices_face8);
        node_indices_all_faces.push_back(node_indices_face9);
        node_indices_all_faces.push_back(node_indices_face10);
        node_indices_all_faces.push_back(node_indices_face11);
        node_indices_all_faces.push_back(node_indices_face12);
        node_indices_all_faces.push_back(node_indices_face13);

        /*
         * Face numbering
         *     _____      basal network  apical network
         *    / 8   /			   _____ 				 _____
         *   /_____/11 			|\ 2 /				|\ 5 /
         *   |\   /  				| \ / 				| \ /
         *   10\ //\  		  	|1 | 					|4 |
         * 7 |  13  \				| /0\ 				| /3\
         *   | / \12/				|/___\				|/___\
         *   |/9__\/
         *			6
         */

        std::vector<MonolayerVertexElementType> face_types_all;
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);

        std::vector<MonolayerVertexElement<2, 3>*> face_elements;
        for (unsigned i = 0; i < 14; i++)
        {
            std::vector<Node<3>*> temp_nodes;
            std::vector<MonolayerVertexElementType> temp_types;
            unsigned max_index = (i == 0 || i == 2 || i == 3 || i == 5) ? 3 : 4;
            for (unsigned j = 0; j < max_index; j++)
            {
                temp_nodes.push_back(nodes[node_indices_all_faces[i][j]]);
                temp_types.push_back(node_types[node_indices_all_faces[i][j]]);
            }
            face_elements.push_back(new MonolayerVertexElement<2, 3>(i, face_types_all[i], temp_nodes, temp_types));
            if ((i > 5 && i < 8) || (i > 10))
            {
                face_elements[i]->SetAsBoundaryFace();
            }
        }

        unsigned face_indices_0[5] = { 0, 6, 9, 12, 3 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_0;
        for (unsigned index = 0; index < 5; index++)
        {
            face_elements_0.push_back(face_elements[face_indices_0[index]]);
        }
        std::vector<bool> orientations_0 = { false, false, true, false, false };

        unsigned face_indices_1[6] = { 1, 9, 10, 13, 7, 4 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_1;
        for (unsigned index = 0; index < 6; index++)
        {
            face_elements_1.push_back(face_elements[face_indices_1[index]]);
        }
        std::vector<bool> orientations_1 = { false, false, false, true, false, false };

        unsigned face_indices_2[5] = { 2, 8, 10, 11, 5 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_2;
        for (unsigned index = 0; index < 5; index++)
        {
            face_elements_2.push_back(face_elements[face_indices_2[index]]);
        }
        std::vector<bool> orientations_2 = { false, false, true, false, false };

        unsigned* all_face_indices[3] = { face_indices_0, face_indices_1, face_indices_2 };

        std::vector<MonolayerVertexElement<3, 3>*> vertex_elements;
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, face_elements_0, orientations_0));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(1, MonolayerVertexElementType::Undetermined, face_elements_1, orientations_1));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(2, MonolayerVertexElementType::Undetermined, face_elements_2, orientations_2));

        MutableMonolayerVertexMesh<3, 3> vertex_mesh(nodes, face_elements, vertex_elements);

        // Test that the orientations are correct prior to remeshing
        // Consistency check to see if mesh is correctly defined
        for (unsigned i = 0; i < 3; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1 * 2.0 / 1.5);

        // Perform a T1 swap on nodes 4 and 5
        // TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), face_elements[9]), "The basal edge of a lateral face seems to have more than two nodes, has wrong nodes or is not oriented correctly");
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(5), vertex_mesh.GetNode(4), face_elements[13]);

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0.0, 1e-3);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1.0, 1e-3);

        // Test that each element contains the correct faces following the rearrangement
        bool found_all_faces = true;
        for (unsigned i = 0; i < 3; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            unsigned num_faces = (i == 0 || i == 2) ? 5 : 6;
            for (unsigned j_face = 0; j_face < num_faces; j_face++)
            {
                unsigned global_face_index = *(all_face_indices[i] + j_face);
                MonolayerVertexElement<2, 3>* p_face_expected = face_elements[global_face_index];
                if (global_face_index == 13)
                    continue;
                if (p_element->GetFaceLocalIndex(p_face_expected) == UINT_MAX)
                {
                    found_all_faces = false;
                    break;
                }
            }
            TS_ASSERT(found_all_faces);
        }
        // Check that central face changed elements correctly.
        bool swapped_face_in_correct_elements = vertex_elements[0]->GetFaceLocalIndex(face_elements[13]) != UINT_MAX
            && vertex_elements[1]->GetFaceLocalIndex(face_elements[13]) == UINT_MAX
            && vertex_elements[2]->GetFaceLocalIndex(face_elements[13]) != UINT_MAX;

        TS_ASSERT(swapped_face_in_correct_elements);

        // Check that swapped nodes are in correct elements

        bool test = vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Basal)->GetNodeLocalIndex(4) == UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Apical)->GetNodeLocalIndex(10) == UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Apical)->GetNodeLocalIndex(11) != UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Basal)->GetNodeLocalIndex(5) != UINT_MAX;
        TS_ASSERT(test);

        bool swapped_nodes_in_correct_elements_0 = vertex_elements[0]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_1 = vertex_elements[1]->GetNodeLocalIndex(4) == UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(10) == UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_2 = vertex_elements[2]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(11) != UINT_MAX;

        TS_ASSERT(swapped_nodes_in_correct_elements_0);
        TS_ASSERT(swapped_nodes_in_correct_elements_1);
        TS_ASSERT(swapped_nodes_in_correct_elements_2);

        TS_ASSERT_EQUALS(vertex_elements[0]->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_elements[1]->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[2]->GetNumNodes(), 8u);

        TS_ASSERT_EQUALS(vertex_elements[0]->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[1]->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_elements[2]->GetNumFaces(), 6u);

        // Test that the orientations are correct
        for (unsigned i = 0; i < 3; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.2 + 0.2 * sqrt(41.0) + 0.6, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.0 + 0.2 * sqrt(41.0) + 0.4, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.2 + 0.2 * sqrt(41.0) + 0.6, 1e-6);

        // Test boundary faces
        for (unsigned face_index = 0; face_index < face_elements.size(); face_index++)
        {
            if ((face_index > 5 && face_index < 8) || (face_index > 10 && face_index < 13))
            {
                TS_ASSERT(face_elements[face_index]->IsBoundaryFace());
            }
            else
            {
                TS_ASSERT(!(face_elements[face_index]->IsBoundaryFace()));
            }
        }

        // Test boundary nodes
        for (unsigned node_index = 0; node_index < nodes.size(); node_index++)
        {
            bool boundary_node = nodes[node_index]->IsBoundaryNode();
            if (node_index == 5 || node_index == 11)
            {
                TS_ASSERT(!boundary_node);
            }
            else
            {
                TS_ASSERT(boundary_node);
            }
        }

        // Test T1 swap location tracking
        std::vector<c_vector<double, 3> > t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 1u);
        TS_ASSERT_DELTA(t1_locations[0][0], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][1], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][2], 0.5, 1e-6);

        // Test T1 swap location clearing
        vertex_mesh.ClearLocationsOfT1Swaps();
        t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 0u);
    }

    void TestPerformT1SwapOnBoundary3()
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * We will test that that a T1 swap of the two central nodes is correctly implemented.
         *                 basal network
         *  /\    /|			  3     2
         * /  \  / |				|\   /|
         * |\  \/| |				| \5/ |
         * | \ / | |				|  |  |
         * |  |  | |basal		| /4\ |
         * | / \ | /				|/___\|
         * |/___\|/apical		0			1
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, true, 0.5, 0.6, 0.0));

        nodes.push_back(new Node<3>(6, true, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, true, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(8, true, 1.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(9, true, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(10, false, 0.5, 0.4, 1.0));
        nodes.push_back(new Node<3>(11, true, 0.5, 0.6, 1.0));

        std::vector<MonolayerVertexElementType> node_types;
        for (unsigned i = 0; i < 12; i++)
        {
            if (i < 6)
            {
                node_types.push_back(MonolayerVertexElementType::Basal);
            }
            else
            {
                node_types.push_back(MonolayerVertexElementType::Apical);
            }
        }

        std::vector<std::vector<unsigned> > node_indices_all_faces;
        // basal
        std::vector<unsigned> node_indices_face0 = { 0, 4, 1 };
        std::vector<unsigned> node_indices_face1 = { 0, 3, 5, 4 };
        std::vector<unsigned> node_indices_face2 = { 2, 1, 4, 5 };
        // apical
        std::vector<unsigned> node_indices_face3 = { 6, 7, 10 };
        std::vector<unsigned> node_indices_face4 = { 6, 10, 11, 9 };
        std::vector<unsigned> node_indices_face5 = { 11, 10, 7, 8 };
        // lateral outer
        std::vector<unsigned> node_indices_face6 = { 0, 1, 7, 6 };
        std::vector<unsigned> node_indices_face7 = { 3, 0, 6, 9 };
        std::vector<unsigned> node_indices_face8 = { 1, 2, 8, 7 };
        // lateral inner, normals always point to the right
        std::vector<unsigned> node_indices_face9 = { 0, 4, 10, 6 };
        std::vector<unsigned> node_indices_face10 = { 5, 3, 9, 11 };
        std::vector<unsigned> node_indices_face11 = { 5, 2, 8, 11 };
        std::vector<unsigned> node_indices_face12 = { 1, 4, 10, 7 };
        std::vector<unsigned> node_indices_face13 = { 4, 5, 11, 10 };

        node_indices_all_faces.push_back(node_indices_face0);
        node_indices_all_faces.push_back(node_indices_face1);
        node_indices_all_faces.push_back(node_indices_face2);
        node_indices_all_faces.push_back(node_indices_face3);
        node_indices_all_faces.push_back(node_indices_face4);
        node_indices_all_faces.push_back(node_indices_face5);
        node_indices_all_faces.push_back(node_indices_face6);
        node_indices_all_faces.push_back(node_indices_face7);
        node_indices_all_faces.push_back(node_indices_face8);
        node_indices_all_faces.push_back(node_indices_face9);
        node_indices_all_faces.push_back(node_indices_face10);
        node_indices_all_faces.push_back(node_indices_face11);
        node_indices_all_faces.push_back(node_indices_face12);
        node_indices_all_faces.push_back(node_indices_face13);

        /*
         * Face numbering
         *                basal network  apical network
         *    /\    /|
         *   /10\  / |				|\   /|				|\   /|
         *   |\  \/| |				| \ / |				| \ / |
         *   | \ /11 |8			|1 | 2|				|4 | 5|
         * 7 |  13 | |				| /0\ |				| /3\ |
         *   | / \12 /				|/___\|				|/___\|
         *   |/9__\|/
         *			6
         */

        std::vector<MonolayerVertexElementType> face_types_all;
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);

        std::vector<MonolayerVertexElement<2, 3>*> face_elements;
        for (unsigned i = 0; i < 14; i++)
        {
            std::vector<Node<3>*> temp_nodes;
            std::vector<MonolayerVertexElementType> temp_types;
            unsigned max_index = (i == 0 || i == 3) ? 3 : 4;
            for (unsigned j = 0; j < max_index; j++)
            {
                temp_nodes.push_back(nodes[node_indices_all_faces[i][j]]);
                temp_types.push_back(node_types[node_indices_all_faces[i][j]]);
            }
            face_elements.push_back(new MonolayerVertexElement<2, 3>(i, face_types_all[i], temp_nodes, temp_types));
            if (i > 5 && i < 12 && i != 9)
            {
                face_elements[i]->SetAsBoundaryFace();
            }
        }

        unsigned face_indices_0[5] = { 0, 6, 9, 12, 3 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_0;
        for (unsigned index = 0; index < 5; index++)
        {
            face_elements_0.push_back(face_elements[face_indices_0[index]]);
        }
        std::vector<bool> orientations_0 = { false, false, true, false, false };

        unsigned face_indices_1[6] = { 1, 7, 9, 10, 13, 4 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_1;
        for (unsigned index = 0; index < 6; index++)
        {
            face_elements_1.push_back(face_elements[face_indices_1[index]]);
        }
        std::vector<bool> orientations_1 = { false, false, false, false, false, false };

        unsigned face_indices_3[6] = { 2, 8, 11, 12, 13, 5 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_3;
        for (unsigned index = 0; index < 6; index++)
        {
            face_elements_3.push_back(face_elements[face_indices_3[index]]);
        }
        std::vector<bool> orientations_3 = { false, false, true, true, true, false };

        unsigned* all_face_indices[3] = { face_indices_0, face_indices_1, face_indices_3 };

        std::vector<MonolayerVertexElement<3, 3>*> vertex_elements;
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, face_elements_0, orientations_0));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(1, MonolayerVertexElementType::Undetermined, face_elements_1, orientations_1));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(2, MonolayerVertexElementType::Undetermined, face_elements_3, orientations_3));

        MutableMonolayerVertexMesh<3, 3> vertex_mesh(nodes, face_elements, vertex_elements);

        // Test that the orientations are correct prior to remeshing
        // Consistency check to see if mesh is correctly defined
        for (unsigned i = 0; i < 3; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1 * 2.0 / 1.5);

        // Perform a T1 swap on nodes 4 and 5
        TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), face_elements[7]), "The basal edge of a lateral face seems to have more than two nodes, has wrong nodes or is not oriented correctly");
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), face_elements[13]);

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0.0, 1e-3);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.6, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.4, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5, 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1.0, 1e-3);

        // Test that each element contains the correct faces following the rearrangement
        bool found_all_faces = true;
        for (unsigned i = 0; i < 3; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            unsigned num_faces = (i == 0) ? 5 : 6;
            for (unsigned j_face = 0; j_face < num_faces; j_face++)
            {
                unsigned global_face_index = *(all_face_indices[i] + j_face);
                MonolayerVertexElement<2, 3>* p_face_expected = face_elements[global_face_index];
                if (global_face_index == 13)
                    continue;
                if (p_element->GetFaceLocalIndex(p_face_expected) == UINT_MAX)
                {
                    found_all_faces = false;
                    break;
                }
            }
            TS_ASSERT(found_all_faces);
        }
        // Check that central face changed elements correctly.
        bool swapped_face_in_correct_elements = vertex_elements[0]->GetFaceLocalIndex(face_elements[13]) != UINT_MAX
            && vertex_elements[1]->GetFaceLocalIndex(face_elements[13]) == UINT_MAX
            && vertex_elements[2]->GetFaceLocalIndex(face_elements[13]) == UINT_MAX;

        TS_ASSERT(swapped_face_in_correct_elements);

        // Check that swapped nodes are in correct elements

        bool test = vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Basal)->GetNodeLocalIndex(4) == UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Apical)->GetNodeLocalIndex(10) == UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Apical)->GetNodeLocalIndex(11) != UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Basal)->GetNodeLocalIndex(5) != UINT_MAX;
        TS_ASSERT(test);

        bool swapped_nodes_in_correct_elements_0 = vertex_elements[0]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_1 = vertex_elements[1]->GetNodeLocalIndex(4) == UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(10) == UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_3 = vertex_elements[2]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(5) == UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(11) == UINT_MAX;

        TS_ASSERT(swapped_nodes_in_correct_elements_0);
        TS_ASSERT(swapped_nodes_in_correct_elements_1);
        TS_ASSERT(swapped_nodes_in_correct_elements_3);

        TS_ASSERT_EQUALS(vertex_elements[0]->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_elements[1]->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[2]->GetNumNodes(), 6u);

        TS_ASSERT_EQUALS(vertex_elements[0]->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[1]->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_elements[2]->GetNumFaces(), 5u);

        // Test that the orientations are correct
        for (unsigned i = 0; i < 3; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 0.3, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), 0.2, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 1.2 + 0.2 * sqrt(41.0) + 0.6, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 1.0 + 0.2 * sqrt(41.0) + 0.4, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(2), 1.0 + 0.2 * sqrt(41.0) + 0.4, 1e-6);

        // Test boundary faces
        for (unsigned face_index = 0; face_index < face_elements.size(); face_index++)
        {
            if (face_index > 5 && face_index != 9 && face_index != 12)
            {
                TS_ASSERT(face_elements[face_index]->IsBoundaryFace());
            }
            else
            {
                TS_ASSERT(!(face_elements[face_index]->IsBoundaryFace()));
            }
        }

        // Test boundary nodes
        for (unsigned node_index = 0; node_index < nodes.size(); node_index++)
        {
            bool boundary_node = nodes[node_index]->IsBoundaryNode();
            TS_ASSERT(boundary_node);
        }

        // Test T1 swap location tracking
        std::vector<c_vector<double, 3> > t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 1u);
        TS_ASSERT_DELTA(t1_locations[0][0], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][1], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][2], 0.5, 1e-6);

        // Test T1 swap location clearing
        vertex_mesh.ClearLocationsOfT1Swaps();
        t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 0u);
    }

    void TestPerformT1SwapAndIdentifySwapTypeNonOrthogonal()
    {
        /*
         * Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
         * We will test that that a T1 swap of the two central nodes is correctly implemented.
         *   _____         basal network
         *  /     /|			  3_____2
         * /_____/ |				|\   /|
         * |\   /| |				| \5/ |
         * | \ / | |				|  \  |
         * |  \  | |basal		| /4\ |
         * | / \ | /				|/___\|
         * |/___\|/apical		0			1
         */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.6, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.4, 0.6, 0.0));

        nodes.push_back(new Node<3>(6, true, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, true, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(8, true, 1.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(9, true, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(10, false, 0.6, 0.4, 1.0));
        nodes.push_back(new Node<3>(11, false, 0.4, 0.6, 1.0));

        std::vector<MonolayerVertexElementType> node_types;
        for (unsigned i = 0; i < 12; i++)
        {
            if (i < 6)
            {
                node_types.push_back(MonolayerVertexElementType::Basal);
            }
            else
            {
                node_types.push_back(MonolayerVertexElementType::Apical);
            }
        }

        std::vector<std::vector<unsigned> > node_indices_all_faces;
        // basal
        std::vector<unsigned> node_indices_face0 = { 0, 4, 1 };
        std::vector<unsigned> node_indices_face1 = { 0, 3, 5, 4 };
        std::vector<unsigned> node_indices_face2 = { 3, 2, 5 };
        std::vector<unsigned> node_indices_face3 = { 5, 2, 1, 4 };
        // apical
        std::vector<unsigned> node_indices_face4 = { 6, 7, 10 };
        std::vector<unsigned> node_indices_face5 = { 6, 10, 11, 9 };
        std::vector<unsigned> node_indices_face6 = { 9, 11, 8 };
        std::vector<unsigned> node_indices_face7 = { 11, 10, 7, 8 };
        // lateral outer
        std::vector<unsigned> node_indices_face8 = { 0, 1, 7, 6 };
        std::vector<unsigned> node_indices_face9 = { 3, 0, 6, 9 };
        std::vector<unsigned> node_indices_face10 = { 2, 3, 9, 8 };
        std::vector<unsigned> node_indices_face11 = { 1, 2, 8, 7 };
        // lateral inner, normals always point to the right
        std::vector<unsigned> node_indices_face12 = { 0, 4, 10, 6 };
        std::vector<unsigned> node_indices_face13 = { 5, 3, 9, 11 };
        std::vector<unsigned> node_indices_face14 = { 5, 2, 8, 11 };
        std::vector<unsigned> node_indices_face15 = { 1, 4, 10, 7 };
        std::vector<unsigned> node_indices_face16 = { 4, 5, 11, 10 };

        node_indices_all_faces.push_back(node_indices_face0);
        node_indices_all_faces.push_back(node_indices_face1);
        node_indices_all_faces.push_back(node_indices_face2);
        node_indices_all_faces.push_back(node_indices_face3);
        node_indices_all_faces.push_back(node_indices_face4);
        node_indices_all_faces.push_back(node_indices_face5);
        node_indices_all_faces.push_back(node_indices_face6);
        node_indices_all_faces.push_back(node_indices_face7);
        node_indices_all_faces.push_back(node_indices_face8);
        node_indices_all_faces.push_back(node_indices_face9);
        node_indices_all_faces.push_back(node_indices_face10);
        node_indices_all_faces.push_back(node_indices_face11);
        node_indices_all_faces.push_back(node_indices_face12);
        node_indices_all_faces.push_back(node_indices_face13);
        node_indices_all_faces.push_back(node_indices_face14);
        node_indices_all_faces.push_back(node_indices_face15);
        node_indices_all_faces.push_back(node_indices_face16);

        /*
         * Face numbering
         *     _____      basal network  apical network
         *    / 10  /|			   _____ 				 _____
         *   /_____/ |				|\ 2 /|				|\ 6 /|
         *   |\   /| |				| \ / |				| \ / |
         *   13\ /14 |11			|1 | 3|				|5 | 7|
         * 9 |  16 | |				| /0\ |				| /4\ |
         *   | / \15 /				|/___\|				|/___\|
         *   |/12_\|/
         *			8
         */

        std::vector<MonolayerVertexElementType> face_types_all;
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Basal);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Apical);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);
        face_types_all.push_back(MonolayerVertexElementType::Lateral);

        std::vector<MonolayerVertexElement<2, 3>*> face_elements;
        for (unsigned i = 0; i < 17; i++)
        {
            std::vector<Node<3>*> temp_nodes;
            std::vector<MonolayerVertexElementType> temp_types;
            unsigned max_index = (i % 2 == 0 && i < 7) ? 3 : 4;
            for (unsigned j = 0; j < max_index; j++)
            {
                temp_nodes.push_back(nodes[node_indices_all_faces[i][j]]);
                temp_types.push_back(node_types[node_indices_all_faces[i][j]]);
            }
            face_elements.push_back(new MonolayerVertexElement<2, 3>(i, face_types_all[i], temp_nodes, temp_types));
            if (i > 7 && i < 12)
            {
                face_elements[i]->SetAsBoundaryFace();
            }
        }

        unsigned face_indices_0[5] = { 0, 8, 12, 15, 4 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_0;
        for (unsigned index = 0; index < 5; index++)
        {
            face_elements_0.push_back(face_elements[face_indices_0[index]]);
        }
        std::vector<bool> orientations_0 = { false, false, true, false, false };

        unsigned face_indices_1[6] = { 1, 12, 13, 16, 9, 5 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_1;
        for (unsigned index = 0; index < 6; index++)
        {
            face_elements_1.push_back(face_elements[face_indices_1[index]]);
        }
        std::vector<bool> orientations_1 = { false, false, false, false, false, false };

        unsigned face_indices_2[5] = { 2, 10, 13, 14, 6 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_2;
        for (unsigned index = 0; index < 5; index++)
        {
            face_elements_2.push_back(face_elements[face_indices_2[index]]);
        }
        std::vector<bool> orientations_2 = { false, false, true, false, false };

        unsigned face_indices_3[6] = { 3, 11, 14, 15, 16, 7 };
        std::vector<MonolayerVertexElement<2, 3>*> face_elements_3;
        for (unsigned index = 0; index < 6; index++)
        {
            face_elements_3.push_back(face_elements[face_indices_3[index]]);
        }
        std::vector<bool> orientations_3 = { false, false, true, true, true, false };

        unsigned* all_face_indices[4] = { face_indices_0, face_indices_1, face_indices_2, face_indices_3 };

        std::vector<MonolayerVertexElement<3, 3>*> vertex_elements;
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, face_elements_0, orientations_0));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(1, MonolayerVertexElementType::Undetermined, face_elements_1, orientations_1));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(2, MonolayerVertexElementType::Undetermined, face_elements_2, orientations_2));
        vertex_elements.push_back(new MonolayerVertexElement<3, 3>(3, MonolayerVertexElementType::Undetermined, face_elements_3, orientations_3));

        MutableMonolayerVertexMesh<3, 3> vertex_mesh(nodes, face_elements, vertex_elements);

        // Test that the orientations are correct prior to remeshing
        // Consistency check to see if mesh is correctly defined
        for (unsigned i = 0; i < 4; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Set the threshold distance between vertices for a T1 swap as follows, to ease calculations
        vertex_mesh.SetCellRearrangementThreshold(0.1 * 2.0 / 1.5);

        // Perform a T1 swap on nodes 4 and 5
        TS_ASSERT_THROWS_THIS(vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), face_elements[11]), "The basal edge of a lateral face seems to have more than two nodes, has wrong nodes or is not oriented correctly");
        vertex_mesh.IdentifySwapType(vertex_mesh.GetNode(4), vertex_mesh.GetNode(5), face_elements[16]);

        // Test that each moved node has the correct location following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[0], 0.5 + 0.1 / sqrt(2.0), 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[1], 0.5 + 0.1 / sqrt(2.0), 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(4)->rGetLocation()[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[0], 0.5 - 0.1 / sqrt(2.0), 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[1], 0.5 - 0.1 / sqrt(2.0), 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(5)->rGetLocation()[2], 0.0, 1e-3);

        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[0], 0.5 + 0.1 / sqrt(2.0), 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[1], 0.5 + 0.1 / sqrt(2.0), 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(10)->rGetLocation()[2], 1.0, 1e-8);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[0], 0.5 - 0.1 / sqrt(2.0), 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[1], 0.5 - 0.1 / sqrt(2.0), 1e-3);
        TS_ASSERT_DELTA(vertex_mesh.GetNode(11)->rGetLocation()[2], 1.0, 1e-3);

        // Test that each element contains the correct faces following the rearrangement
        bool found_all_faces = true;
        for (unsigned i = 0; i < 4; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            unsigned num_faces = (i == 0 || i == 2) ? 5 : 6;
            for (unsigned j_face = 0; j_face < num_faces; j_face++)
            {
                unsigned global_face_index = *(all_face_indices[i] + j_face);
                MonolayerVertexElement<2, 3>* p_face_expected = face_elements[global_face_index];
                if (global_face_index == 16)
                    continue;
                if (p_element->GetFaceLocalIndex(p_face_expected) == UINT_MAX)
                {
                    found_all_faces = false;
                    break;
                }
            }
            TS_ASSERT(found_all_faces);
        }
        // Check that central face changed elements correctly.
        bool swapped_face_in_correct_elements = vertex_elements[0]->GetFaceLocalIndex(face_elements[16]) != UINT_MAX
            && vertex_elements[1]->GetFaceLocalIndex(face_elements[16]) == UINT_MAX
            && vertex_elements[2]->GetFaceLocalIndex(face_elements[16]) != UINT_MAX
            && vertex_elements[3]->GetFaceLocalIndex(face_elements[16]) == UINT_MAX;

        TS_ASSERT(swapped_face_in_correct_elements);

        // Check that swapped nodes are in correct elements

        bool test = vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Basal)->GetNodeLocalIndex(4) == UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Apical)->GetNodeLocalIndex(10) == UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Apical)->GetNodeLocalIndex(11) != UINT_MAX
            && vertex_mesh.GetFaceOfType(1, MonolayerVertexElementType::Basal)->GetNodeLocalIndex(5) != UINT_MAX;
        TS_ASSERT(test);

        bool swapped_nodes_in_correct_elements_0 = vertex_elements[0]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[0]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_1 = vertex_elements[1]->GetNodeLocalIndex(4) == UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(10) == UINT_MAX
            && vertex_elements[1]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_2 = vertex_elements[2]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(5) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[2]->GetNodeLocalIndex(11) != UINT_MAX;

        bool swapped_nodes_in_correct_elements_3 = vertex_elements[3]->GetNodeLocalIndex(4) != UINT_MAX
            && vertex_elements[3]->GetNodeLocalIndex(5) == UINT_MAX
            && vertex_elements[3]->GetNodeLocalIndex(10) != UINT_MAX
            && vertex_elements[3]->GetNodeLocalIndex(11) == UINT_MAX;

        TS_ASSERT(swapped_nodes_in_correct_elements_0);
        TS_ASSERT(swapped_nodes_in_correct_elements_1);
        TS_ASSERT(swapped_nodes_in_correct_elements_2);
        TS_ASSERT(swapped_nodes_in_correct_elements_3);

        TS_ASSERT_EQUALS(vertex_elements[0]->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_elements[1]->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[2]->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_elements[3]->GetNumNodes(), 6u);

        TS_ASSERT_EQUALS(vertex_elements[0]->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[1]->GetNumFaces(), 5u);
        TS_ASSERT_EQUALS(vertex_elements[2]->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_elements[3]->GetNumFaces(), 5u);

        // Test that the orientations are correct
        for (unsigned i = 0; i < 4; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Test that each element has the correct area and perimeter following the rearrangement
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), (1.0 - 0.5 + 0.1 / sqrt(2.0)) / 2.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), (0.5 - 0.1 / sqrt(2.0)) / 2.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(2), (1.0 - 0.5 + 0.1 / sqrt(2.0)) / 2.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(3), (0.5 - 0.1 / sqrt(2.0)) / 2.0, 1e-6);

        // Test boundary faces
        for (unsigned face_index = 0; face_index < face_elements.size(); face_index++)
        {
            if (face_index > 7 && face_index < 12)
            {
                TS_ASSERT(face_elements[face_index]->IsBoundaryFace());
            }
            else
            {
                TS_ASSERT(!(face_elements[face_index]->IsBoundaryFace()));
            }
        }

        // Test T1 swap location tracking
        std::vector<c_vector<double, 3> > t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 1u);
        TS_ASSERT_DELTA(t1_locations[0][0], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][1], 0.5, 1e-6);
        TS_ASSERT_DELTA(t1_locations[0][2], 0.5, 1e-6);

        // Test T1 swap location clearing
        vertex_mesh.ClearLocationsOfT1Swaps();
        t1_locations = vertex_mesh.GetLocationsOfT1Swaps();
        TS_ASSERT_EQUALS(t1_locations.size(), 0u);
    }

    void TestDivideElement()
    {
        // We test this with an elongated hexagon with long axis parallel to the x-axis
        /*
         *             /o---o---o\         /o---o---o\
         *            / :_ _:_ _: \       / :_ _|_ _: \
         *           o /         \ o =>  o /    |    \ o
         *           |\           /|     |\     |     /|
         *            \\o---o---o/ /      \\o---o---o/ /
         *             \|___|___| /        \|___|___| /
         */
        std::vector<Node<3>*> regular_nodes;
        // Basal
        regular_nodes.push_back(new Node<3>(0, false, 3.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(1, false, 2.5, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(3, false, -2.5, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(4, false, -3.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(5, false, -2.5, -1.0, 0.0));
        regular_nodes.push_back(new Node<3>(6, false, 0.0, -1.0, 0.0));
        regular_nodes.push_back(new Node<3>(7, false, 2.5, -1.0, 0.0));
        // Apical
        regular_nodes.push_back(new Node<3>(8, false, 3.0, 0.0, 1.0));
        regular_nodes.push_back(new Node<3>(9, false, 2.5, 1.0, 1.0));
        regular_nodes.push_back(new Node<3>(10, false, 0.0, 1.0, 1.0));
        regular_nodes.push_back(new Node<3>(11, false, -2.5, 1.0, 1.0));
        regular_nodes.push_back(new Node<3>(12, false, -3.0, 0.0, 1.0));
        regular_nodes.push_back(new Node<3>(13, false, -2.5, -1.0, 1.0));
        regular_nodes.push_back(new Node<3>(14, false, 0.0, -1.0, 1.0));
        regular_nodes.push_back(new Node<3>(15, false, 2.5, -1.0, 1.0));

        std::vector<MonolayerVertexElementType> regular_types(8, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> regular_types_2(8, MonolayerVertexElementType::Apical);
        regular_types.insert(regular_types.end(), regular_types_2.begin(), regular_types_2.end());

        std::vector<std::vector<int> > face_indices;
        std::vector<int> face_0 = { 7, 6, 5, 4, 3, 2, 1, 0 };
        std::vector<int> face_1 = { 8, 9, 10, 11, 12, 13, 14, 15 }; // Begin at inner node
        std::vector<int> face_2 = { 0, 8, 15, 7 };
        std::vector<int> face_3 = { 7, 15, 14, 6 };
        std::vector<int> face_4 = { 6, 14, 13, 5 };
        std::vector<int> face_5 = { 5, 13, 12, 4 };
        std::vector<int> face_6 = { 4, 12, 11, 3 };
        std::vector<int> face_7 = { 3, 11, 10, 2 };
        std::vector<int> face_8 = { 2, 10, 9, 1 };
        std::vector<int> face_9 = { 1, 9, 8, 0 };
        face_indices = { face_0, face_1, face_2, face_3, face_4, face_5, face_6, face_7, face_8, face_9 };

        std::vector<MonolayerVertexElement<2, 3>*> regular_faces;
        for (unsigned i = 0; i < 10; i++)
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
        std::vector<bool> regular_faces_orientations(10, false);

        regular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, regular_faces, regular_faces_orientations, regular_nodes, regular_types));
        MutableMonolayerVertexMesh<3, 3> vertex_mesh(regular_nodes, regular_faces, regular_elements);

        // Test that the orientations are correct prior to dividing
        // Consistency check to see if mesh is correctly defined
        MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(0);
        c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(0);
        for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
        {
            MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
            double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
            c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
            std::vector<c_vector<double, 3> > normal_vectors;
            vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
            passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

            // Check anti-clockwise orientation of nodes
            for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
            {
                TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
            }
            // Check for correct orientation of face in element
            c_vector<double, 3> normal_vector;
            vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
            TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
        }

        // Now we divide the element along the x=0 plane
        std::array<unsigned, 4> indices_division_nodes = { 2, 6, 14, 10 };
        vertex_mesh.DivideElement(regular_elements[0], indices_division_nodes, true); // old element is left

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 13u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 16u)

        // Test that the orientations are correct
        for (unsigned i = 0; i < 2; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);

                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Test if the elements contain the correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 10u);

        std::set<unsigned> old_node_indices = { 2, 3, 4, 5, 6, 10, 11, 12, 13, 14 };
        std::set<unsigned> new_node_indices = { 0, 1, 2, 6, 7, 8, 9, 10, 14, 15 };
        for (unsigned i = 0; i < 2; i++)
        {
            std::set<unsigned> temp_nodes;
            unsigned nodes_in_elt = vertex_mesh.GetElement(i)->GetNumNodes();
            for (unsigned ind_node = 0; ind_node < nodes_in_elt; ind_node++)
            {
                temp_nodes.insert(vertex_mesh.GetElement(i)->GetNode(ind_node)->GetIndex());
            }
            if (i == 0)
            {
                TS_ASSERT(old_node_indices == temp_nodes);
            }
            else
            {
                TS_ASSERT(new_node_indices == temp_nodes);
            }
        }
        // Test if the elmeents contain the correct faces
        std::vector<MonolayerVertexElement<2, 3>*> old_faces = { regular_faces[0], regular_faces[1], regular_faces[4],
                                                                 regular_faces[5], regular_faces[6], regular_faces[7], vertex_mesh.GetFace(10) };

        std::vector<MonolayerVertexElement<2, 3>*> new_faces = { vertex_mesh.GetFace(11), vertex_mesh.GetFace(12), regular_faces[2],
                                                                 regular_faces[3], regular_faces[8], regular_faces[9], vertex_mesh.GetFace(10) };

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 7u);

        for (unsigned i = 0; i < 2; i++)
        {
            unsigned faces_in_elt = vertex_mesh.GetElement(i)->GetNumFaces();
            for (unsigned ind_face = 0; ind_face < faces_in_elt; ind_face++)
            {
                std::vector<MonolayerVertexElement<2, 3>*>& vect = i == 0 ? old_faces : new_faces;
                MonolayerVertexElement<2, 3>* p_face = vertex_mesh.GetElement(i)->GetFace(ind_face);
                auto it = std::find(vect.begin(), vect.end(), p_face);

                TS_ASSERT(it != vect.end());
                vect.erase(it);
            }
        }

        // Test the correct orientation.
        TS_ASSERT_EQUALS(regular_elements[0]->GetIndex(), 0u);
        TS_ASSERT(vertex_mesh.GetFace(10)->GetFaceType() == MonolayerVertexElementType::Lateral);
        unsigned local_ind_face_old = vertex_mesh.GetElement(0)->GetFaceLocalIndex(vertex_mesh.GetFace(10));
        unsigned local_ind_face_new = vertex_mesh.GetElement(1)->GetFaceLocalIndex(vertex_mesh.GetFace(10));
        // old element is outside the face
        TS_ASSERT(vertex_mesh.GetElement(0)->FaceIsOrientatedClockwise(local_ind_face_old) == true);
        TS_ASSERT(vertex_mesh.GetElement(1)->FaceIsOrientatedClockwise(local_ind_face_new) == false);

        /*
         * Now we check if the other orientation works and if we can divide at the zeroth node of the apical face.
         */

        std::vector<Node<3>*> irregular_nodes;
        // Basal
        irregular_nodes.push_back(new Node<3>(0, false, 3.0, 0.0, 0.0));
        irregular_nodes.push_back(new Node<3>(1, false, 2.5, 1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(3, false, -2.5, 1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(4, false, -3.0, 0.0, 0.0));
        irregular_nodes.push_back(new Node<3>(5, false, -2.5, -1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(6, false, 0.0, -1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(7, false, 2.5, -1.0, 0.0));
        // Apical
        irregular_nodes.push_back(new Node<3>(8, false, 3.0, 0.0, 1.0));
        irregular_nodes.push_back(new Node<3>(9, false, 2.5, 1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(10, false, 0.0, 1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(11, false, -2.5, 1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(12, false, -3.0, 0.0, 1.0));
        irregular_nodes.push_back(new Node<3>(13, false, -2.5, -1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(14, false, 0.0, -1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(15, false, 2.5, -1.0, 1.0));

        std::vector<MonolayerVertexElementType> irregular_types(8, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> irregular_types_2(8, MonolayerVertexElementType::Apical);
        irregular_types.insert(irregular_types.end(), irregular_types_2.begin(), irregular_types_2.end());

        std::vector<std::vector<int> > ir_face_indices;
        std::vector<int> ir_face_0 = { 7, 6, 5, 4, 3, 2, 1, 0 };
        std::vector<int> ir_face_1 = { 10, 11, 12, 13, 14, 15, 8, 9 }; // Begin at separation node
        std::vector<int> ir_face_2 = { 0, 8, 15, 7 };
        std::vector<int> ir_face_3 = { 7, 15, 14, 6 };
        std::vector<int> ir_face_4 = { 6, 14, 13, 5 };
        std::vector<int> ir_face_5 = { 5, 13, 12, 4 };
        std::vector<int> ir_face_6 = { 4, 12, 11, 3 };
        std::vector<int> ir_face_7 = { 3, 11, 10, 2 };
        std::vector<int> ir_face_8 = { 2, 10, 9, 1 };
        std::vector<int> ir_face_9 = { 1, 9, 8, 0 };
        ir_face_indices = { ir_face_0, ir_face_1, ir_face_2, ir_face_3, ir_face_4, ir_face_5, ir_face_6, ir_face_7, ir_face_8, ir_face_9 };

        std::vector<MonolayerVertexElement<2, 3>*> irregular_faces;
        for (unsigned i = 0; i < 10; i++)
        {
            std::vector<Node<3>*> face_nodes;
            std::vector<MonolayerVertexElementType> node_types;
            for (auto it = ir_face_indices[i].begin(); it != ir_face_indices[i].end(); ++it)
            {
                face_nodes.push_back(irregular_nodes[*it]);
                node_types.push_back(irregular_types[*it]);
            }
            MonolayerVertexElementType facetype;
            if (i == 0)
                facetype = MonolayerVertexElementType::Basal;
            else if (i == 1)
                facetype = MonolayerVertexElementType::Apical;
            else
                facetype = MonolayerVertexElementType::Lateral;

            irregular_faces.push_back(new MonolayerVertexElement<2, 3>(i, facetype, face_nodes, node_types));
        }
        std::vector<MonolayerVertexElement<3, 3>*> irregular_elements;
        std::vector<bool> irregular_faces_orientations(10, false);

        irregular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, irregular_faces, irregular_faces_orientations, irregular_nodes, irregular_types));
        MutableMonolayerVertexMesh<3, 3> ir_vertex_mesh(irregular_nodes, irregular_faces, irregular_elements);

        // Test that the orientations are correct prior to dividing
        // Consistency check to see if mesh is correctly defined
        MonolayerVertexElement<3, 3>* ir_p_element = ir_vertex_mesh.GetElement(0);
        c_vector<double, 3> ir_centroid = ir_vertex_mesh.GetCentroidOfElement(0);
        for (unsigned j_face = 0; j_face < ir_p_element->GetNumFaces(); j_face++)
        {
            MonolayerVertexElement<2, 3>* p_face = ir_p_element->GetFace(j_face);
            double orientation_factor = ir_p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
            c_vector<double, 3> passive_center = ir_vertex_mesh.GetPassiveCenterOfFace(p_face);
            std::vector<c_vector<double, 3> > normal_vectors;
            ir_vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
            passive_center = ir_vertex_mesh.GetVectorFromAtoB(ir_centroid, passive_center);

            // Check anti-clockwise orientation of nodes
            for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
            {
                TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
            }
            // Check for correct orientation of face in element
            c_vector<double, 3> normal_vector;
            ir_vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
            TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
        }

        // Now we divide the element along the x=0 plane
        std::array<unsigned, 4> ir_indices_division_nodes = { 2, 6, 14, 10 };
        ir_vertex_mesh.DivideElement(irregular_elements[0], ir_indices_division_nodes, false); // new element is left

        // Test that the orientations are correct
        for (unsigned i = 0; i < 2; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = ir_vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = ir_vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);

                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = ir_vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                ir_vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = ir_vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                ir_vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Test if the elements contain the correct nodes
        TS_ASSERT_EQUALS(ir_vertex_mesh.GetElement(0)->GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(ir_vertex_mesh.GetElement(1)->GetNumNodes(), 10u);

        new_node_indices = { 2, 3, 4, 5, 6, 10, 11, 12, 13, 14 };
        old_node_indices = { 0, 1, 2, 6, 7, 8, 9, 10, 14, 15 };
        for (unsigned i = 0; i < 2; i++)
        {
            std::set<unsigned> temp_nodes;
            unsigned nodes_in_elt = ir_vertex_mesh.GetElement(i)->GetNumNodes();
            for (unsigned ind_node = 0; ind_node < nodes_in_elt; ind_node++)
            {
                temp_nodes.insert(ir_vertex_mesh.GetElement(i)->GetNode(ind_node)->GetIndex());
            }
            if (i == 0)
            {
                TS_ASSERT(old_node_indices == temp_nodes);
            }
            else
            {
                TS_ASSERT(new_node_indices == temp_nodes);
            }
        }
        // Test if the elmeents contain the correct faces
        TS_ASSERT_EQUALS(ir_vertex_mesh.GetElement(0)->GetNumFaces(), 7u);
        TS_ASSERT_EQUALS(ir_vertex_mesh.GetElement(1)->GetNumFaces(), 7u);

        old_faces = { irregular_faces[0], irregular_faces[1], irregular_faces[2],
                      irregular_faces[3], irregular_faces[8], irregular_faces[9], ir_vertex_mesh.GetFace(10) };

        new_faces = { ir_vertex_mesh.GetFace(11), ir_vertex_mesh.GetFace(12), irregular_faces[4],
                      irregular_faces[5], irregular_faces[6], irregular_faces[7], ir_vertex_mesh.GetFace(10) };

        for (unsigned i = 0; i < 2; i++)
        {
            unsigned faces_in_elt = ir_vertex_mesh.GetElement(i)->GetNumFaces();
            for (unsigned ind_face = 0; ind_face < faces_in_elt; ind_face++)
            {
                std::vector<MonolayerVertexElement<2, 3>*>& vect = i == 0 ? old_faces : new_faces;
                MonolayerVertexElement<2, 3>* p_face = ir_vertex_mesh.GetElement(i)->GetFace(ind_face);
                auto it = std::find(vect.begin(), vect.end(), p_face);

                TS_ASSERT(it != vect.end());
                vect.erase(it);
            }
        }

        // Test the correct orientation.
        TS_ASSERT_EQUALS(irregular_elements[0]->GetIndex(), 0u);
        TS_ASSERT(ir_vertex_mesh.GetFace(10)->GetFaceType() == MonolayerVertexElementType::Lateral);
        local_ind_face_old = ir_vertex_mesh.GetElement(0)->GetFaceLocalIndex(ir_vertex_mesh.GetFace(10));
        local_ind_face_new = ir_vertex_mesh.GetElement(1)->GetFaceLocalIndex(ir_vertex_mesh.GetFace(10));
        // old element is outside the face
        TS_ASSERT(ir_vertex_mesh.GetElement(0)->FaceIsOrientatedClockwise(local_ind_face_old) == false);
        TS_ASSERT(ir_vertex_mesh.GetElement(1)->FaceIsOrientatedClockwise(local_ind_face_new) == true);
    }

    void TestDivideElementAlternativeApicalSide1()
    {
        // We test this with an elongated hexagon with long axis parallel to the x-axis
        /* However, now we use different starting/ending nodes in apical side
         *             /o---o---o\         /o---o---o\
         *            / :_ _:_ _: \       / :_ _|_ _: \
         *           o /         \ o =>  o /    |    \ o
         *           |\           /|     |\     |     /|
         *            \\o---o---o/ /      \\o---o---o/ /
         *             \|___|___| /        \|___|___| /
         */
        std::vector<Node<3>*> regular_nodes;
        // Basal
        regular_nodes.push_back(new Node<3>(0, false, 3.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(1, false, 2.5, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(3, false, -2.5, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(4, false, -3.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(5, false, -2.5, -1.0, 0.0));
        regular_nodes.push_back(new Node<3>(6, false, 0.0, -1.0, 0.0));
        regular_nodes.push_back(new Node<3>(7, false, 2.5, -1.0, 0.0));
        // Apical
        regular_nodes.push_back(new Node<3>(8, false, 3.0, 0.0, 1.0));
        regular_nodes.push_back(new Node<3>(9, false, 2.5, 1.0, 1.0));
        regular_nodes.push_back(new Node<3>(10, false, 0.0, 1.0, 1.0));
        regular_nodes.push_back(new Node<3>(11, false, -2.5, 1.0, 1.0));
        regular_nodes.push_back(new Node<3>(12, false, -3.0, 0.0, 1.0));
        regular_nodes.push_back(new Node<3>(13, false, -2.5, -1.0, 1.0));
        regular_nodes.push_back(new Node<3>(14, false, 0.0, -1.0, 1.0));
        regular_nodes.push_back(new Node<3>(15, false, 2.5, -1.0, 1.0));

        std::vector<MonolayerVertexElementType> regular_types(8, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> regular_types_2(8, MonolayerVertexElementType::Apical);
        regular_types.insert(regular_types.end(), regular_types_2.begin(), regular_types_2.end());

        std::vector<std::vector<int> > face_indices;
        std::vector<int> face_0 = { 7, 6, 5, 4, 3, 2, 1, 0 };
        std::vector<int> face_1 = { 10, 11, 12, 13, 14, 15, 8, 9 }; // Begin at separtion node
        std::vector<int> face_2 = { 0, 8, 15, 7 };
        std::vector<int> face_3 = { 7, 15, 14, 6 };
        std::vector<int> face_4 = { 6, 14, 13, 5 };
        std::vector<int> face_5 = { 5, 13, 12, 4 };
        std::vector<int> face_6 = { 4, 12, 11, 3 };
        std::vector<int> face_7 = { 3, 11, 10, 2 };
        std::vector<int> face_8 = { 2, 10, 9, 1 };
        std::vector<int> face_9 = { 1, 9, 8, 0 };
        face_indices = { face_0, face_1, face_2, face_3, face_4, face_5, face_6, face_7, face_8, face_9 };

        std::vector<MonolayerVertexElement<2, 3>*> regular_faces;
        for (unsigned i = 0; i < 10; i++)
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
        std::vector<bool> regular_faces_orientations(10, false);

        regular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, regular_faces, regular_faces_orientations, regular_nodes, regular_types));
        MutableMonolayerVertexMesh<3, 3> vertex_mesh(regular_nodes, regular_faces, regular_elements);

        // Test that the orientations are correct prior to dividing
        // Consistency check to see if mesh is correctly defined
        MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(0);
        c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(0);
        for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
        {
            MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
            double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
            c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
            std::vector<c_vector<double, 3> > normal_vectors;
            vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
            passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

            // Check anti-clockwise orientation of nodes
            for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
            {
                TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
            }
            // Check for correct orientation of face in element
            c_vector<double, 3> normal_vector;
            vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
            TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
        }

        // Now we divide the element along the x=0 plane
        std::array<unsigned, 4> indices_division_nodes = { 2, 6, 14, 10 };
        vertex_mesh.DivideElement(regular_elements[0], indices_division_nodes, true); // old element is left

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 13u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 16u)

        // Test that the orientations are correct
        for (unsigned i = 0; i < 2; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);

                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Test if the elements contain the correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 10u);

        std::set<unsigned> old_node_indices = { 2, 3, 4, 5, 6, 10, 11, 12, 13, 14 };
        std::set<unsigned> new_node_indices = { 0, 1, 2, 6, 7, 8, 9, 10, 14, 15 };
        for (unsigned i = 0; i < 2; i++)
        {
            std::set<unsigned> temp_nodes;
            unsigned nodes_in_elt = vertex_mesh.GetElement(i)->GetNumNodes();
            for (unsigned ind_node = 0; ind_node < nodes_in_elt; ind_node++)
            {
                temp_nodes.insert(vertex_mesh.GetElement(i)->GetNode(ind_node)->GetIndex());
            }
            if (i == 0)
            {
                TS_ASSERT(old_node_indices == temp_nodes);
            }
            else
            {
                TS_ASSERT(new_node_indices == temp_nodes);
            }
        }
        // Test if the elmeents contain the correct faces
        std::vector<MonolayerVertexElement<2, 3>*> old_faces = { regular_faces[0], regular_faces[1], regular_faces[4],
                                                                 regular_faces[5], regular_faces[6], regular_faces[7], vertex_mesh.GetFace(10) };

        std::vector<MonolayerVertexElement<2, 3>*> new_faces = { vertex_mesh.GetFace(11), vertex_mesh.GetFace(12), regular_faces[2],
                                                                 regular_faces[3], regular_faces[8], regular_faces[9], vertex_mesh.GetFace(10) };

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 7u);

        for (unsigned i = 0; i < 2; i++)
        {
            unsigned faces_in_elt = vertex_mesh.GetElement(i)->GetNumFaces();
            for (unsigned ind_face = 0; ind_face < faces_in_elt; ind_face++)
            {
                std::vector<MonolayerVertexElement<2, 3>*>& vect = i == 0 ? old_faces : new_faces;
                MonolayerVertexElement<2, 3>* p_face = vertex_mesh.GetElement(i)->GetFace(ind_face);
                auto it = std::find(vect.begin(), vect.end(), p_face);

                TS_ASSERT(it != vect.end());
                vect.erase(it);
            }
        }

        // Test the correct orientation.
        TS_ASSERT_EQUALS(regular_elements[0]->GetIndex(), 0u);
        TS_ASSERT(vertex_mesh.GetFace(10)->GetFaceType() == MonolayerVertexElementType::Lateral);
        unsigned local_ind_face_old = vertex_mesh.GetElement(0)->GetFaceLocalIndex(vertex_mesh.GetFace(10));
        unsigned local_ind_face_new = vertex_mesh.GetElement(1)->GetFaceLocalIndex(vertex_mesh.GetFace(10));
        // old element is outside the face
        TS_ASSERT(vertex_mesh.GetElement(0)->FaceIsOrientatedClockwise(local_ind_face_old) == true);
        TS_ASSERT(vertex_mesh.GetElement(1)->FaceIsOrientatedClockwise(local_ind_face_new) == false);

        /*
         * Now we check if the other orientation works and if we can divide at the zeroth node of the apical face.
         */

        std::vector<Node<3>*> irregular_nodes;
        // Basal
        irregular_nodes.push_back(new Node<3>(0, false, 3.0, 0.0, 0.0));
        irregular_nodes.push_back(new Node<3>(1, false, 2.5, 1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(3, false, -2.5, 1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(4, false, -3.0, 0.0, 0.0));
        irregular_nodes.push_back(new Node<3>(5, false, -2.5, -1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(6, false, 0.0, -1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(7, false, 2.5, -1.0, 0.0));
        // Apical
        irregular_nodes.push_back(new Node<3>(8, false, 3.0, 0.0, 1.0));
        irregular_nodes.push_back(new Node<3>(9, false, 2.5, 1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(10, false, 0.0, 1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(11, false, -2.5, 1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(12, false, -3.0, 0.0, 1.0));
        irregular_nodes.push_back(new Node<3>(13, false, -2.5, -1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(14, false, 0.0, -1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(15, false, 2.5, -1.0, 1.0));

        std::vector<MonolayerVertexElementType> irregular_types(8, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> irregular_types_2(8, MonolayerVertexElementType::Apical);
        irregular_types.insert(irregular_types.end(), irregular_types_2.begin(), irregular_types_2.end());

        std::vector<std::vector<int> > ir_face_indices;
        std::vector<int> ir_face_0 = { 7, 6, 5, 4, 3, 2, 1, 0 };
        std::vector<int> ir_face_1 = { 8, 9, 10, 11, 12, 13, 14, 15 }; // Begin with inner node
        std::vector<int> ir_face_2 = { 0, 8, 15, 7 };
        std::vector<int> ir_face_3 = { 7, 15, 14, 6 };
        std::vector<int> ir_face_4 = { 6, 14, 13, 5 };
        std::vector<int> ir_face_5 = { 5, 13, 12, 4 };
        std::vector<int> ir_face_6 = { 4, 12, 11, 3 };
        std::vector<int> ir_face_7 = { 3, 11, 10, 2 };
        std::vector<int> ir_face_8 = { 2, 10, 9, 1 };
        std::vector<int> ir_face_9 = { 1, 9, 8, 0 };
        ir_face_indices = { ir_face_0, ir_face_1, ir_face_2, ir_face_3, ir_face_4, ir_face_5, ir_face_6, ir_face_7, ir_face_8, ir_face_9 };

        std::vector<MonolayerVertexElement<2, 3>*> irregular_faces;
        for (unsigned i = 0; i < 10; i++)
        {
            std::vector<Node<3>*> face_nodes;
            std::vector<MonolayerVertexElementType> node_types;
            for (auto it = ir_face_indices[i].begin(); it != ir_face_indices[i].end(); ++it)
            {
                face_nodes.push_back(irregular_nodes[*it]);
                node_types.push_back(irregular_types[*it]);
            }
            MonolayerVertexElementType facetype;
            if (i == 0)
                facetype = MonolayerVertexElementType::Basal;
            else if (i == 1)
                facetype = MonolayerVertexElementType::Apical;
            else
                facetype = MonolayerVertexElementType::Lateral;

            irregular_faces.push_back(new MonolayerVertexElement<2, 3>(i, facetype, face_nodes, node_types));
        }
        std::vector<MonolayerVertexElement<3, 3>*> irregular_elements;
        std::vector<bool> irregular_faces_orientations(10, false);

        irregular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, irregular_faces, irregular_faces_orientations, irregular_nodes, irregular_types));
        MutableMonolayerVertexMesh<3, 3> ir_vertex_mesh(irregular_nodes, irregular_faces, irregular_elements);

        // Test that the orientations are correct prior to dividing
        // Consistency check to see if mesh is correctly defined
        MonolayerVertexElement<3, 3>* ir_p_element = ir_vertex_mesh.GetElement(0);
        c_vector<double, 3> ir_centroid = ir_vertex_mesh.GetCentroidOfElement(0);
        for (unsigned j_face = 0; j_face < ir_p_element->GetNumFaces(); j_face++)
        {
            MonolayerVertexElement<2, 3>* p_face = ir_p_element->GetFace(j_face);
            double orientation_factor = ir_p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
            c_vector<double, 3> passive_center = ir_vertex_mesh.GetPassiveCenterOfFace(p_face);
            std::vector<c_vector<double, 3> > normal_vectors;
            ir_vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
            passive_center = ir_vertex_mesh.GetVectorFromAtoB(ir_centroid, passive_center);

            // Check anti-clockwise orientation of nodes
            for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
            {
                TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
            }
            // Check for correct orientation of face in element
            c_vector<double, 3> normal_vector;
            ir_vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
            TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
        }

        // Now we divide the element along the x=0 plane
        std::array<unsigned, 4> ir_indices_division_nodes = { 2, 6, 14, 10 };
        ir_vertex_mesh.DivideElement(irregular_elements[0], ir_indices_division_nodes, false); // new element is left

        // Test that the orientations are correct
        for (unsigned i = 0; i < 2; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = ir_vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = ir_vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);

                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = ir_vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                ir_vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = ir_vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                ir_vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Test if the elements contain the correct nodes
        TS_ASSERT_EQUALS(ir_vertex_mesh.GetElement(0)->GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(ir_vertex_mesh.GetElement(1)->GetNumNodes(), 10u);

        new_node_indices = { 2, 3, 4, 5, 6, 10, 11, 12, 13, 14 };
        old_node_indices = { 0, 1, 2, 6, 7, 8, 9, 10, 14, 15 };
        for (unsigned i = 0; i < 2; i++)
        {
            std::set<unsigned> temp_nodes;
            unsigned nodes_in_elt = ir_vertex_mesh.GetElement(i)->GetNumNodes();
            for (unsigned ind_node = 0; ind_node < nodes_in_elt; ind_node++)
            {
                temp_nodes.insert(ir_vertex_mesh.GetElement(i)->GetNode(ind_node)->GetIndex());
            }
            if (i == 0)
            {
                TS_ASSERT(old_node_indices == temp_nodes);
            }
            else
            {
                TS_ASSERT(new_node_indices == temp_nodes);
            }
        }
        // Test if the elmeents contain the correct faces
        TS_ASSERT_EQUALS(ir_vertex_mesh.GetElement(0)->GetNumFaces(), 7u);
        TS_ASSERT_EQUALS(ir_vertex_mesh.GetElement(1)->GetNumFaces(), 7u);

        old_faces = { irregular_faces[0], irregular_faces[1], irregular_faces[2],
                      irregular_faces[3], irregular_faces[8], irregular_faces[9], ir_vertex_mesh.GetFace(10) };

        new_faces = { ir_vertex_mesh.GetFace(11), ir_vertex_mesh.GetFace(12), irregular_faces[4],
                      irregular_faces[5], irregular_faces[6], irregular_faces[7], ir_vertex_mesh.GetFace(10) };

        for (unsigned i = 0; i < 2; i++)
        {
            unsigned faces_in_elt = ir_vertex_mesh.GetElement(i)->GetNumFaces();
            for (unsigned ind_face = 0; ind_face < faces_in_elt; ind_face++)
            {
                std::vector<MonolayerVertexElement<2, 3>*>& vect = i == 0 ? old_faces : new_faces;
                MonolayerVertexElement<2, 3>* p_face = ir_vertex_mesh.GetElement(i)->GetFace(ind_face);
                auto it = std::find(vect.begin(), vect.end(), p_face);

                TS_ASSERT(it != vect.end());
                vect.erase(it);
            }
        }

        // Test the correct orientation.
        TS_ASSERT_EQUALS(irregular_elements[0]->GetIndex(), 0u);
        TS_ASSERT(ir_vertex_mesh.GetFace(10)->GetFaceType() == MonolayerVertexElementType::Lateral);
        local_ind_face_old = ir_vertex_mesh.GetElement(0)->GetFaceLocalIndex(ir_vertex_mesh.GetFace(10));
        local_ind_face_new = ir_vertex_mesh.GetElement(1)->GetFaceLocalIndex(ir_vertex_mesh.GetFace(10));
        // old element is outside the face
        TS_ASSERT(ir_vertex_mesh.GetElement(0)->FaceIsOrientatedClockwise(local_ind_face_old) == false);
        TS_ASSERT(ir_vertex_mesh.GetElement(1)->FaceIsOrientatedClockwise(local_ind_face_new) == true);
    }

    void TestDivideElementAlternativeApicalSide2()
    {
        // We test this with an elongated hexagon with long axis parallel to the x-axis
        /* However, now we use different starting/ending nodes in apical side
         *             /o---o---o\         /o---o---o\
         *            / :_ _:_ _: \       / :_ _|_ _: \
         *           o /         \ o =>  o /    |    \ o
         *           |\           /|     |\     |     /|
         *            \\o---o---o/ /      \\o---o---o/ /
         *             \|___|___| /        \|___|___| /
         */
        std::vector<Node<3>*> regular_nodes;
        // Basal
        regular_nodes.push_back(new Node<3>(0, false, 3.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(1, false, 2.5, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(3, false, -2.5, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(4, false, -3.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(5, false, -2.5, -1.0, 0.0));
        regular_nodes.push_back(new Node<3>(6, false, 0.0, -1.0, 0.0));
        regular_nodes.push_back(new Node<3>(7, false, 2.5, -1.0, 0.0));
        // Apical
        regular_nodes.push_back(new Node<3>(8, false, 3.0, 0.0, 1.0));
        regular_nodes.push_back(new Node<3>(9, false, 2.5, 1.0, 1.0));
        regular_nodes.push_back(new Node<3>(10, false, 0.0, 1.0, 1.0));
        regular_nodes.push_back(new Node<3>(11, false, -2.5, 1.0, 1.0));
        regular_nodes.push_back(new Node<3>(12, false, -3.0, 0.0, 1.0));
        regular_nodes.push_back(new Node<3>(13, false, -2.5, -1.0, 1.0));
        regular_nodes.push_back(new Node<3>(14, false, 0.0, -1.0, 1.0));
        regular_nodes.push_back(new Node<3>(15, false, 2.5, -1.0, 1.0));

        std::vector<MonolayerVertexElementType> regular_types(8, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> regular_types_2(8, MonolayerVertexElementType::Apical);
        regular_types.insert(regular_types.end(), regular_types_2.begin(), regular_types_2.end());

        std::vector<std::vector<int> > face_indices;
        std::vector<int> face_0 = { 7, 6, 5, 4, 3, 2, 1, 0 };
        std::vector<int> face_1 = { 11, 12, 13, 14, 15, 8, 9, 10 }; // End at separtion node
        std::vector<int> face_2 = { 0, 8, 15, 7 };
        std::vector<int> face_3 = { 7, 15, 14, 6 };
        std::vector<int> face_4 = { 6, 14, 13, 5 };
        std::vector<int> face_5 = { 5, 13, 12, 4 };
        std::vector<int> face_6 = { 4, 12, 11, 3 };
        std::vector<int> face_7 = { 3, 11, 10, 2 };
        std::vector<int> face_8 = { 2, 10, 9, 1 };
        std::vector<int> face_9 = { 1, 9, 8, 0 };
        face_indices = { face_0, face_1, face_2, face_3, face_4, face_5, face_6, face_7, face_8, face_9 };

        std::vector<MonolayerVertexElement<2, 3>*> regular_faces;
        for (unsigned i = 0; i < 10; i++)
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
        std::vector<bool> regular_faces_orientations(10, false);

        regular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, regular_faces, regular_faces_orientations, regular_nodes, regular_types));
        MutableMonolayerVertexMesh<3, 3> vertex_mesh(regular_nodes, regular_faces, regular_elements);

        // Test that the orientations are correct prior to dividing
        // Consistency check to see if mesh is correctly defined
        MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(0);
        c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(0);
        for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
        {
            MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
            double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
            c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
            std::vector<c_vector<double, 3> > normal_vectors;
            vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
            passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

            // Check anti-clockwise orientation of nodes
            for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
            {
                TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
            }
            // Check for correct orientation of face in element
            c_vector<double, 3> normal_vector;
            vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
            TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
        }

        // Now we divide the element along the x=0 plane
        std::array<unsigned, 4> indices_division_nodes = { 2, 6, 14, 10 };
        vertex_mesh.DivideElement(regular_elements[0], indices_division_nodes, true); // old element is left

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 13u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 16u)

        // Test that the orientations are correct
        for (unsigned i = 0; i < 2; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);

                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Test if the elements contain the correct nodes
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 10u);

        std::set<unsigned> old_node_indices = { 2, 3, 4, 5, 6, 10, 11, 12, 13, 14 };
        std::set<unsigned> new_node_indices = { 0, 1, 2, 6, 7, 8, 9, 10, 14, 15 };
        for (unsigned i = 0; i < 2; i++)
        {
            std::set<unsigned> temp_nodes;
            unsigned nodes_in_elt = vertex_mesh.GetElement(i)->GetNumNodes();
            for (unsigned ind_node = 0; ind_node < nodes_in_elt; ind_node++)
            {
                temp_nodes.insert(vertex_mesh.GetElement(i)->GetNode(ind_node)->GetIndex());
            }
            if (i == 0)
            {
                TS_ASSERT(old_node_indices == temp_nodes);
            }
            else
            {
                TS_ASSERT(new_node_indices == temp_nodes);
            }
        }
        // Test if the elmeents contain the correct faces
        std::vector<MonolayerVertexElement<2, 3>*> old_faces = { regular_faces[0], regular_faces[1], regular_faces[4],
                                                                 regular_faces[5], regular_faces[6], regular_faces[7], vertex_mesh.GetFace(10) };

        std::vector<MonolayerVertexElement<2, 3>*> new_faces = { vertex_mesh.GetFace(11), vertex_mesh.GetFace(12), regular_faces[2],
                                                                 regular_faces[3], regular_faces[8], regular_faces[9], vertex_mesh.GetFace(10) };

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumFaces(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 7u);

        for (unsigned i = 0; i < 2; i++)
        {
            unsigned faces_in_elt = vertex_mesh.GetElement(i)->GetNumFaces();
            for (unsigned ind_face = 0; ind_face < faces_in_elt; ind_face++)
            {
                std::vector<MonolayerVertexElement<2, 3>*>& vect = i == 0 ? old_faces : new_faces;
                MonolayerVertexElement<2, 3>* p_face = vertex_mesh.GetElement(i)->GetFace(ind_face);
                auto it = std::find(vect.begin(), vect.end(), p_face);

                TS_ASSERT(it != vect.end());
                vect.erase(it);
            }
        }

        // Test the correct orientation.
        TS_ASSERT_EQUALS(regular_elements[0]->GetIndex(), 0u);
        TS_ASSERT(vertex_mesh.GetFace(10)->GetFaceType() == MonolayerVertexElementType::Lateral);
        unsigned local_ind_face_old = vertex_mesh.GetElement(0)->GetFaceLocalIndex(vertex_mesh.GetFace(10));
        unsigned local_ind_face_new = vertex_mesh.GetElement(1)->GetFaceLocalIndex(vertex_mesh.GetFace(10));
        // old element is outside the face
        TS_ASSERT(vertex_mesh.GetElement(0)->FaceIsOrientatedClockwise(local_ind_face_old) == true);
        TS_ASSERT(vertex_mesh.GetElement(1)->FaceIsOrientatedClockwise(local_ind_face_new) == false);

        /*
         * Now we check if the other orientation works and if we can divide at the zeroth node of the apical face.
         */

        std::vector<Node<3>*> irregular_nodes;
        // Basal
        irregular_nodes.push_back(new Node<3>(0, false, 3.0, 0.0, 0.0));
        irregular_nodes.push_back(new Node<3>(1, false, 2.5, 1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(3, false, -2.5, 1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(4, false, -3.0, 0.0, 0.0));
        irregular_nodes.push_back(new Node<3>(5, false, -2.5, -1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(6, false, 0.0, -1.0, 0.0));
        irregular_nodes.push_back(new Node<3>(7, false, 2.5, -1.0, 0.0));
        // Apical
        irregular_nodes.push_back(new Node<3>(8, false, 3.0, 0.0, 1.0));
        irregular_nodes.push_back(new Node<3>(9, false, 2.5, 1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(10, false, 0.0, 1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(11, false, -2.5, 1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(12, false, -3.0, 0.0, 1.0));
        irregular_nodes.push_back(new Node<3>(13, false, -2.5, -1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(14, false, 0.0, -1.0, 1.0));
        irregular_nodes.push_back(new Node<3>(15, false, 2.5, -1.0, 1.0));

        std::vector<MonolayerVertexElementType> irregular_types(8, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> irregular_types_2(8, MonolayerVertexElementType::Apical);
        irregular_types.insert(irregular_types.end(), irregular_types_2.begin(), irregular_types_2.end());

        std::vector<std::vector<int> > ir_face_indices;
        std::vector<int> ir_face_0 = { 7, 6, 5, 4, 3, 2, 1, 0 };
        std::vector<int> ir_face_1 = { 11, 12, 13, 14, 15, 8, 9, 10 }; // End at separation node
        std::vector<int> ir_face_2 = { 0, 8, 15, 7 };
        std::vector<int> ir_face_3 = { 7, 15, 14, 6 };
        std::vector<int> ir_face_4 = { 6, 14, 13, 5 };
        std::vector<int> ir_face_5 = { 5, 13, 12, 4 };
        std::vector<int> ir_face_6 = { 4, 12, 11, 3 };
        std::vector<int> ir_face_7 = { 3, 11, 10, 2 };
        std::vector<int> ir_face_8 = { 2, 10, 9, 1 };
        std::vector<int> ir_face_9 = { 1, 9, 8, 0 };
        ir_face_indices = { ir_face_0, ir_face_1, ir_face_2, ir_face_3, ir_face_4, ir_face_5, ir_face_6, ir_face_7, ir_face_8, ir_face_9 };

        std::vector<MonolayerVertexElement<2, 3>*> irregular_faces;
        for (unsigned i = 0; i < 10; i++)
        {
            std::vector<Node<3>*> face_nodes;
            std::vector<MonolayerVertexElementType> node_types;
            for (auto it = ir_face_indices[i].begin(); it != ir_face_indices[i].end(); ++it)
            {
                face_nodes.push_back(irregular_nodes[*it]);
                node_types.push_back(irregular_types[*it]);
            }
            MonolayerVertexElementType facetype;
            if (i == 0)
                facetype = MonolayerVertexElementType::Basal;
            else if (i == 1)
                facetype = MonolayerVertexElementType::Apical;
            else
                facetype = MonolayerVertexElementType::Lateral;

            irregular_faces.push_back(new MonolayerVertexElement<2, 3>(i, facetype, face_nodes, node_types));
        }
        std::vector<MonolayerVertexElement<3, 3>*> irregular_elements;
        std::vector<bool> irregular_faces_orientations(10, false);

        irregular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, irregular_faces, irregular_faces_orientations, irregular_nodes, irregular_types));
        MutableMonolayerVertexMesh<3, 3> ir_vertex_mesh(irregular_nodes, irregular_faces, irregular_elements);

        // Test that the orientations are correct prior to dividing
        // Consistency check to see if mesh is correctly defined
        MonolayerVertexElement<3, 3>* ir_p_element = ir_vertex_mesh.GetElement(0);
        c_vector<double, 3> ir_centroid = ir_vertex_mesh.GetCentroidOfElement(0);
        for (unsigned j_face = 0; j_face < ir_p_element->GetNumFaces(); j_face++)
        {
            MonolayerVertexElement<2, 3>* p_face = ir_p_element->GetFace(j_face);
            double orientation_factor = ir_p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
            c_vector<double, 3> passive_center = ir_vertex_mesh.GetPassiveCenterOfFace(p_face);
            std::vector<c_vector<double, 3> > normal_vectors;
            ir_vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
            passive_center = ir_vertex_mesh.GetVectorFromAtoB(ir_centroid, passive_center);

            // Check anti-clockwise orientation of nodes
            for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
            {
                TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
            }
            // Check for correct orientation of face in element
            c_vector<double, 3> normal_vector;
            ir_vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
            TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
        }

        // Now we divide the element along the x=0 plane
        std::array<unsigned, 4> ir_indices_division_nodes = { 2, 6, 14, 10 };
        ir_vertex_mesh.DivideElement(irregular_elements[0], ir_indices_division_nodes, false); // new element is left

        // Test that the orientations are correct
        for (unsigned i = 0; i < 2; i++)
        {
            MonolayerVertexElement<3, 3>* p_element = ir_vertex_mesh.GetElement(i);
            c_vector<double, 3> centroid = ir_vertex_mesh.GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);

                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = ir_vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                ir_vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = ir_vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                ir_vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Test if the elements contain the correct nodes
        TS_ASSERT_EQUALS(ir_vertex_mesh.GetElement(0)->GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(ir_vertex_mesh.GetElement(1)->GetNumNodes(), 10u);

        new_node_indices = { 2, 3, 4, 5, 6, 10, 11, 12, 13, 14 };
        old_node_indices = { 0, 1, 2, 6, 7, 8, 9, 10, 14, 15 };
        for (unsigned i = 0; i < 2; i++)
        {
            std::set<unsigned> temp_nodes;
            unsigned nodes_in_elt = ir_vertex_mesh.GetElement(i)->GetNumNodes();
            for (unsigned ind_node = 0; ind_node < nodes_in_elt; ind_node++)
            {
                temp_nodes.insert(ir_vertex_mesh.GetElement(i)->GetNode(ind_node)->GetIndex());
            }
            if (i == 0)
            {
                TS_ASSERT(old_node_indices == temp_nodes);
            }
            else
            {
                TS_ASSERT(new_node_indices == temp_nodes);
            }
        }
        // Test if the elmeents contain the correct faces
        TS_ASSERT_EQUALS(ir_vertex_mesh.GetElement(0)->GetNumFaces(), 7u);
        TS_ASSERT_EQUALS(ir_vertex_mesh.GetElement(1)->GetNumFaces(), 7u);

        old_faces = { irregular_faces[0], irregular_faces[1], irregular_faces[2],
                      irregular_faces[3], irregular_faces[8], irregular_faces[9], ir_vertex_mesh.GetFace(10) };

        new_faces = { ir_vertex_mesh.GetFace(11), ir_vertex_mesh.GetFace(12), irregular_faces[4],
                      irregular_faces[5], irregular_faces[6], irregular_faces[7], ir_vertex_mesh.GetFace(10) };

        for (unsigned i = 0; i < 2; i++)
        {
            unsigned faces_in_elt = ir_vertex_mesh.GetElement(i)->GetNumFaces();
            for (unsigned ind_face = 0; ind_face < faces_in_elt; ind_face++)
            {
                std::vector<MonolayerVertexElement<2, 3>*>& vect = i == 0 ? old_faces : new_faces;
                MonolayerVertexElement<2, 3>* p_face = ir_vertex_mesh.GetElement(i)->GetFace(ind_face);
                auto it = std::find(vect.begin(), vect.end(), p_face);

                TS_ASSERT(it != vect.end());
                vect.erase(it);
            }
        }

        // Test the correct orientation.
        TS_ASSERT_EQUALS(irregular_elements[0]->GetIndex(), 0u);
        TS_ASSERT(ir_vertex_mesh.GetFace(10)->GetFaceType() == MonolayerVertexElementType::Lateral);
        local_ind_face_old = ir_vertex_mesh.GetElement(0)->GetFaceLocalIndex(ir_vertex_mesh.GetFace(10));
        local_ind_face_new = ir_vertex_mesh.GetElement(1)->GetFaceLocalIndex(ir_vertex_mesh.GetFace(10));
        // old element is outside the face
        TS_ASSERT(ir_vertex_mesh.GetElement(0)->FaceIsOrientatedClockwise(local_ind_face_old) == false);
        TS_ASSERT(ir_vertex_mesh.GetElement(1)->FaceIsOrientatedClockwise(local_ind_face_new) == true);
    }

    void TestDivideLateralFaceWithNodes()
    {
        // We test this with a rectange
        /* The two nodes are not part of the element and then are added as separation nodes
         * Orientation of lateral faces should be conserved
         *                 o---.----o       o---o----o
         *                /| _ . _ /| =>   /|_ _| _ /|
         *               / ,      / /     / ,      / /
         *              o-------o/ /     o-------o/ /
         *              |_______| /      |_______| /
         */
        std::vector<Node<3>*> regular_nodes;
        // Basal
        regular_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(1, false, 2.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(2, false, 2.0, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(3, false, 0.0, 1.0, 0.0));
        // Apical
        regular_nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(5, false, 2.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(6, false, 2.0, 1.0, 3.0));
        regular_nodes.push_back(new Node<3>(7, false, 0.0, 1.0, 3.0));
        // Extra nodes
        Node<3>* p_extra_basal = new Node<3>(0, false, 1.0, 1.0, 0.0);
        Node<3>* p_extra_apical = new Node<3>(0, false, 1.0, 1.0, 3.0);

        std::vector<MonolayerVertexElementType> regular_types(4, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> regular_types_2(4, MonolayerVertexElementType::Apical);
        regular_types.insert(regular_types.end(), regular_types_2.begin(), regular_types_2.end());

        std::vector<std::vector<int> > face_indices;
        std::vector<int> face_0 = { 0, 3, 2, 1 };
        std::vector<int> face_1 = { 4, 5, 6, 7 };
        std::vector<int> face_2 = { 1, 5, 4, 0 };
        std::vector<int> face_3 = { 2, 6, 5, 1 };
        std::vector<int> face_4 = { 3, 7, 6, 2 };
        std::vector<int> face_5 = { 0, 4, 7, 3 };
        face_indices = { face_0, face_1, face_2, face_3, face_4, face_5 };

        std::vector<MonolayerVertexElement<2, 3>*> regular_faces;
        for (unsigned i = 0; i < 6; i++)
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
        std::vector<bool> regular_faces_orientations(6, false);

        regular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, regular_faces, regular_faces_orientations, regular_nodes, regular_types));
        MutableMonolayerVertexMesh<3, 3> vertex_mesh(regular_nodes, regular_faces, regular_elements);

        MonolayerVertexElement<2, 3>* p_face_to_divide = regular_faces[4];

        // Check initial state
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 6u);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 6.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 22.0, 1e-6);

        // divide the cell
        vertex_mesh.AddNode(p_extra_apical);
        vertex_mesh.AddNode(p_extra_basal);
        vertex_mesh.DivideLateralFaceWithNodes(p_face_to_divide, p_extra_apical, p_extra_basal);

        // Check final state
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 7u);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 6.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 22.0, 1e-6);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(0)->GetNumNodes(), 5u);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(vertex_mesh.GetFace(0)), 2.0, 1e-6);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(1)->GetNumNodes(), 5u);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(vertex_mesh.GetFace(1)), 2.0, 1e-6);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(3)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(4)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(5)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(6)->GetNumNodes(), 4u);

        // Test that the orientations are correct
        MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(0);
        c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(0);
        for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
        {
            MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);

            double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
            c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
            std::vector<c_vector<double, 3> > normal_vectors;
            vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
            passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

            // Check anti-clockwise orientation of nodes
            for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
            {
                TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
            }
            // Check for correct orientation of face in element
            c_vector<double, 3> normal_vector;
            vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
            TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
        }
    }

    void TestDivideLateralFaceWithNodesNotPlanar()
    {
        // We test this with a rectangle, but nodes outside the face.
        /* The two nodes are not part of the element and then are added as separation nodes
         * Orientation of lateral faces should be conserved
         *                       .              _ o_
         *                       .            _-  | -_
         *                 o--------o       o    .:. o
         *                /| _ _ _ /| =>   /| .:    /|
         *               / ,      / /     / ,      / /
         *              o-------o/ /     o-------o/ /
         *              |_______| /      |_______| /
         */
        std::vector<Node<3>*> regular_nodes;
        // Basal
        regular_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(1, false, 2.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(2, false, 2.0, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(3, false, 0.0, 1.0, 0.0));
        // Apical
        regular_nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(5, false, 2.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(6, false, 2.0, 1.0, 3.0));
        regular_nodes.push_back(new Node<3>(7, false, 0.0, 1.0, 3.0));
        // Extra nodes
        Node<3>* p_extra_basal = new Node<3>(0, false, 1.0, 2.0, 0.0);
        Node<3>* p_extra_apical = new Node<3>(0, false, 1.0, 2.0, 3.0);

        std::vector<MonolayerVertexElementType> regular_types(4, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> regular_types_2(4, MonolayerVertexElementType::Apical);
        regular_types.insert(regular_types.end(), regular_types_2.begin(), regular_types_2.end());

        std::vector<std::vector<int> > face_indices;
        std::vector<int> face_0 = { 0, 3, 2, 1 };
        std::vector<int> face_1 = { 4, 5, 6, 7 };
        std::vector<int> face_2 = { 1, 5, 4, 0 };
        std::vector<int> face_3 = { 2, 6, 5, 1 };
        std::vector<int> face_4 = { 3, 7, 6, 2 };
        std::vector<int> face_5 = { 0, 4, 7, 3 };
        face_indices = { face_0, face_1, face_2, face_3, face_4, face_5 };

        std::vector<MonolayerVertexElement<2, 3>*> regular_faces;
        for (unsigned i = 0; i < 6; i++)
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
        std::vector<bool> regular_faces_orientations(6, false);

        regular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, regular_faces, regular_faces_orientations, regular_nodes, regular_types));
        MutableMonolayerVertexMesh<3, 3> vertex_mesh(regular_nodes, regular_faces, regular_elements);

        MonolayerVertexElement<2, 3>* p_face_to_divide = regular_faces[4];

        // Check initial state
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 6u);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 6.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 22.0, 1e-6);

        // divide the cell
        vertex_mesh.AddNode(p_extra_apical);
        vertex_mesh.AddNode(p_extra_basal);
        vertex_mesh.DivideLateralFaceWithNodes(p_face_to_divide, p_extra_apical, p_extra_basal);

        // Check final state
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 7u);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 9.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 18.0 + 6.0 * sqrt(2), 1e-6);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(0)->GetNumNodes(), 5u);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(vertex_mesh.GetFace(0)), 3.0, 1e-6);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(1)->GetNumNodes(), 5u);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(vertex_mesh.GetFace(1)), 3.0, 1e-6);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(3)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(4)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(5)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(6)->GetNumNodes(), 4u);

        // Test that the orientations are correct
        MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(0);
        c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(0);
        for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
        {
            MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);

            double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
            c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
            std::vector<c_vector<double, 3> > normal_vectors;
            vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
            passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

            // Check anti-clockwise orientation of nodes
            for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
            {
                TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
            }
            // Check for correct orientation of face in element
            c_vector<double, 3> normal_vector;
            vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
            TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
        }
    }

    void TestDivideLateralFaceWithNodesNotPlanarTwoCells()
    {
        // We test this with a rectangle, but nodes outside the face and with a neighbour
        /* The two nodes are not part of the element and then are added as separation nodes
         * Orientation of lateral faces should be conserved
         *										o--------o        o--------o
         *                   /   .    /|       / _ o_   /|
         *                  /    .   / /      /_-  | -_/ /
         *                 o--------o /      o    .:. o /
         *                /| _ _ _ /|/ =>   /| .:    /|/
         *               / ,      / /      / ,      / /
         *              o-------o/ /      o-------o/ /
         *              |_______| /       |_______| /
         */
        std::vector<Node<3>*> regular_nodes;
        // Basal lower
        regular_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(1, false, 2.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(2, false, 2.0, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(3, false, 0.0, 1.0, 0.0));
        // Apical lower
        regular_nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(5, false, 2.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(6, false, 2.0, 1.0, 3.0));
        regular_nodes.push_back(new Node<3>(7, false, 0.0, 1.0, 3.0));
        // Basal upper
        regular_nodes.push_back(new Node<3>(8, false, 0.0, 4.0, 0.0));
        regular_nodes.push_back(new Node<3>(9, false, 2.0, 4.0, 0.0));
        // Apical upper
        regular_nodes.push_back(new Node<3>(10, false, 0.0, 4.0, 3.0));
        regular_nodes.push_back(new Node<3>(11, false, 2.0, 4.0, 3.0));
        // Extra nodes
        Node<3>* p_extra_basal = new Node<3>(0, false, 1.0, 2.0, 0.0);
        Node<3>* p_extra_apical = new Node<3>(0, false, 1.0, 2.0, 3.0);

        std::vector<MonolayerVertexElementType> regular_types(4, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> regular_types_2(4, MonolayerVertexElementType::Apical);
        regular_types.insert(regular_types.end(), regular_types_2.begin(), regular_types_2.end());
        regular_types.push_back(MonolayerVertexElementType::Basal);
        regular_types.push_back(MonolayerVertexElementType::Basal);
        regular_types.push_back(MonolayerVertexElementType::Apical);
        regular_types.push_back(MonolayerVertexElementType::Apical);

        std::vector<std::vector<int> > face_indices_1;
        std::vector<int> face_0 = { 0, 3, 2, 1 };
        std::vector<int> face_1 = { 4, 5, 6, 7 };
        std::vector<int> face_2 = { 1, 5, 4, 0 };
        std::vector<int> face_3 = { 2, 6, 5, 1 };
        std::vector<int> face_4 = { 3, 7, 6, 2 };
        std::vector<int> face_5 = { 0, 4, 7, 3 };
        face_indices_1 = { face_0, face_1, face_2, face_3, face_4, face_5 };

        std::vector<std::vector<int> > face_indices_2;
        std::vector<int> face_6 = { 3, 8, 9, 2 };
        std::vector<int> face_7 = { 6, 11, 10, 7 };
        std::vector<int> face_8 = { 2, 9, 11, 6 };
        std::vector<int> face_9 = { 10, 11, 9, 8 };
        std::vector<int> face_10 = { 3, 7, 10, 8 };
        face_indices_2 = { face_6, face_7, face_8, face_9, face_10, face_4 };

        // Create lower element
        std::vector<MonolayerVertexElement<2, 3>*> regular_faces_1;
        for (unsigned i = 0; i < 6; i++)
        {
            std::vector<Node<3>*> face_nodes;
            std::vector<MonolayerVertexElementType> node_types;
            for (auto it = face_indices_1[i].begin(); it != face_indices_1[i].end(); ++it)
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

            regular_faces_1.push_back(new MonolayerVertexElement<2, 3>(i, facetype, face_nodes, node_types));
        }
        std::vector<MonolayerVertexElement<3, 3>*> regular_elements;
        std::vector<bool> regular_faces_orientations_1(6, false);

        regular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, regular_faces_1, regular_faces_orientations_1));

        // Create upper element
        std::vector<MonolayerVertexElement<2, 3>*> regular_faces_2;
        for (unsigned i = 0; i < 5; i++)
        {
            std::vector<Node<3>*> face_nodes;
            std::vector<MonolayerVertexElementType> node_types;
            for (auto it = face_indices_2[i].begin(); it != face_indices_2[i].end(); ++it)
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

            regular_faces_2.push_back(new MonolayerVertexElement<2, 3>(i + 6, facetype, face_nodes, node_types));
        }
        regular_faces_2.push_back(regular_faces_1[4]);
        std::vector<bool> regular_faces_orientations_2(6, false);
        regular_faces_orientations_2[5] = true;

        regular_elements.push_back(new MonolayerVertexElement<3, 3>(1, MonolayerVertexElementType::Undetermined, regular_faces_2, regular_faces_orientations_2));

        std::vector<MonolayerVertexElement<2, 3>*> regular_faces = regular_faces_1;
        regular_faces.insert(regular_faces.end(), regular_faces_2.begin(), --regular_faces_2.end());

        MutableMonolayerVertexMesh<3, 3> vertex_mesh(regular_nodes, regular_faces, regular_elements);

        MonolayerVertexElement<2, 3>* p_face_to_divide = regular_faces[4];

        // Check initial state
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 11u);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 6.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 18.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 22.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 42.0, 1e-6);

        // divide the cell
        vertex_mesh.AddNode(p_extra_apical);
        vertex_mesh.AddNode(p_extra_basal);
        vertex_mesh.DivideLateralFaceWithNodes(p_face_to_divide, p_extra_apical, p_extra_basal);

        // Check final state
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 14u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 12u);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 9.0, 1e-6);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(0)->GetNumNodes(), 10u);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 18.0 + 6.0 * sqrt(2), 1e-6);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(0)->GetNumNodes(), 5u);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(vertex_mesh.GetFace(0)), 3.0, 1e-6);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(1)->GetNumNodes(), 5u);
        TS_ASSERT_DELTA(vertex_mesh.CalculateAreaOfFace(vertex_mesh.GetFace(1)), 3.0, 1e-6);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(2)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(3)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(4)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(5)->GetNumNodes(), 4u);

        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 15.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 34.0 + 6.0 * sqrt(2), 1e-6);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1)->GetNumFaces(), 7u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(6)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(7)->GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(8)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(9)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(10)->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(vertex_mesh.GetFace(11)->GetNumNodes(), 4u);

        // Test that the orientations are correct
        for (unsigned index_element = 0; index_element < 2; index_element++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(index_element);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(index_element);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);

                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }
    }

    void TestDivideElementAlongGivenAxis()
    {
        // We test this with a rectangle, where we divide by a vector, which is tilted by 10°
        /*
         *                 o--------o       o--o-----o
         *                /| _ _ _ /| =>   /|_ _\ _ /|
         *               / ,      / /     / ,  :   / /
         *              o-------o/ /     o--o----o/ /
         *              |_______| /      |___\:__| /
         */
        std::vector<Node<3>*> regular_nodes;
        // Basal
        regular_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(1, false, 2.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(2, false, 2.0, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(3, false, 0.0, 1.0, 0.0));
        // Apical
        regular_nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(5, false, 2.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(6, false, 2.0, 1.0, 3.0));
        regular_nodes.push_back(new Node<3>(7, false, 0.0, 1.0, 3.0));

        std::vector<MonolayerVertexElementType> regular_types(4, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> regular_types_2(4, MonolayerVertexElementType::Apical);
        regular_types.insert(regular_types.end(), regular_types_2.begin(), regular_types_2.end());

        std::vector<std::vector<int> > face_indices;
        std::vector<int> face_0 = { 0, 3, 2, 1 };
        std::vector<int> face_1 = { 4, 5, 6, 7 };
        std::vector<int> face_2 = { 1, 5, 4, 0 };
        std::vector<int> face_3 = { 2, 6, 5, 1 };
        std::vector<int> face_4 = { 3, 7, 6, 2 };
        std::vector<int> face_5 = { 0, 4, 7, 3 };
        face_indices = { face_0, face_1, face_2, face_3, face_4, face_5 };

        std::vector<MonolayerVertexElement<2, 3>*> regular_faces;
        for (unsigned i = 0; i < 6; i++)
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
        std::vector<bool> regular_faces_orientations(6, false);

        regular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, regular_faces, regular_faces_orientations, regular_nodes, regular_types));
        MutableMonolayerVertexMesh<3, 3> vertex_mesh(regular_nodes, regular_faces, regular_elements);

        // Check initial state
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 6u);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 6.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 22.0, 1e-6);

        // Create vector for normal of division plane
        c_vector<double, 3> normal_axis = zero_vector<double>(3);
        normal_axis[0] = cos(M_PI / 18.0);
        normal_axis[2] = sin(M_PI / 18.0);
        normal_axis[1] = 0.0;

        vertex_mesh.DivideElementAlongGivenAxis(regular_elements[0], normal_axis);

        // Check final state
        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 11u);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 3.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1), 3.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 11.0 + 3.0 / cos(M_PI / 18.0), 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1), 11.0 + 3.0 / cos(M_PI / 18.0), 1e-6);
        for (unsigned face_index = 0; face_index < 11; face_index++)
        {
            TS_ASSERT_EQUALS(vertex_mesh.GetFace(face_index)->GetNumNodes(), 4u);
        }

        // Test that the orientations are correct
        for (unsigned element_index = 0; element_index < 2; element_index++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(element_index);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(element_index);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);

                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }
    }

    void TestDivideElementAlongGivenAxisDistanceMinimization()
    {
        // We test this with a rectangle, where we divide by a vector, which is tilted by 10°
        /* Here we use an axis, where the basal and apical edges do not belong to the same
         * lateral face, which leads to minimization of the basal distance to the plane such that
         * the node is closest to the plane (with 1% minimal distance) on the same face.
         *
         *                 o----o---o       o--o-o---o
         *                /| _ _|_ /| =>   /|_:_\|_ /|
         *               /,       / /     /, :     / /
         *              o----o--o/ /     o--o-o--o/ /
         *              |____|__| /      |___\|__| /
         */
        std::vector<Node<3>*> regular_nodes;
        // Basal
        regular_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(2, false, 2.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(3, false, 2.0, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(4, false, 1.0, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(5, false, 0.0, 1.0, 0.0));
        // Apical
        regular_nodes.push_back(new Node<3>(6, false, 0.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(7, false, 1.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(8, false, 2.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(9, false, 2.0, 1.0, 3.0));
        regular_nodes.push_back(new Node<3>(10, false, 1.0, 1.0, 3.0));
        regular_nodes.push_back(new Node<3>(11, false, 0.0, 1.0, 3.0));

        std::vector<MonolayerVertexElementType> regular_types(6, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> regular_types_2(6, MonolayerVertexElementType::Apical);
        regular_types.insert(regular_types.end(), regular_types_2.begin(), regular_types_2.end());

        std::vector<std::vector<int> > face_indices;
        std::vector<int> face_0 = { 0, 5, 4, 3, 2, 1 };
        std::vector<int> face_1 = { 6, 7, 8, 9, 10, 11 };
        std::vector<int> face_2 = { 1, 7, 6, 0 };
        std::vector<int> face_3 = { 2, 8, 7, 1 };
        std::vector<int> face_4 = { 3, 9, 8, 2 };
        std::vector<int> face_5 = { 4, 10, 9, 3 };
        std::vector<int> face_6 = { 5, 11, 10, 4 };
        std::vector<int> face_7 = { 0, 6, 11, 5 };
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

        regular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, regular_faces, regular_faces_orientations, regular_nodes, regular_types));
        MutableMonolayerVertexMesh<3, 3> vertex_mesh(regular_nodes, regular_faces, regular_elements);

        // Check initial state
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 8u);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(0), 6.0, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(0), 22.0, 1e-6);

        // Create vector for normal of division plane
        c_vector<double, 3> normal_axis = zero_vector<double>(3);
        normal_axis[0] = cos(M_PI / 18.0);
        normal_axis[2] = sin(M_PI / 18.0);
        normal_axis[1] = 0.0;

        vertex_mesh.DivideElementAlongGivenAxis(regular_elements[0], normal_axis);

        // Check final state
        unsigned left_element_index = vertex_mesh.GetVolumeOfElement(0) < vertex_mesh.GetVolumeOfElement(1) ? 0 : 1;

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 13u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(left_element_index)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(left_element_index)->GetNumNodes(), 8u);
        MonolayerVertexElement<2, 3>* p_face_apical_left = vertex_mesh.GetFaceOfType(left_element_index, MonolayerVertexElementType::Apical);
        MonolayerVertexElement<2, 3>* p_face_basal_left = vertex_mesh.GetFaceOfType(left_element_index, MonolayerVertexElementType::Basal);
        TS_ASSERT_EQUALS(p_face_apical_left->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(p_face_basal_left->GetNumNodes(), 4u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1 - left_element_index)->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1 - left_element_index)->GetNumNodes(), 12u);
        MonolayerVertexElement<2, 3>* p_face_apical_right = vertex_mesh.GetFaceOfType(1 - left_element_index, MonolayerVertexElementType::Apical);
        MonolayerVertexElement<2, 3>* p_face_basal_right = vertex_mesh.GetFaceOfType(1 - left_element_index, MonolayerVertexElementType::Basal);
        TS_ASSERT_EQUALS(p_face_apical_right->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_face_basal_right->GetNumNodes(), 6u);

        double V_left = 3.0 * (1.0 - 1.5 * tan(M_PI / 18.0)) + 1.5 * (1.5 * tan(M_PI / 18.0) - 0.01);
        double A_left = V_left * 2.0 + (1.0 - 1.5 * tan(M_PI / 18.0)) + 0.99 + 3.0 + sqrt((1.5 * tan(M_PI / 18.0) - 0.01) * (1.5 * tan(M_PI / 18.0) - 0.01) + 9.0);
        double A_right = 3.03 * 2.0 + (1.5 * tan(M_PI / 18.0) - 0.01) * 1.5 * 2.0 + 1.01 + (1.0 + 1.5 * tan(M_PI / 18.0))
            + 3.0 + +sqrt((1.5 * tan(M_PI / 18.0) - 0.01) * (1.5 * tan(M_PI / 18.0) - 0.01) + 9.0);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(left_element_index), V_left, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetVolumeOfElement(1 - left_element_index), 6.0 - V_left, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(left_element_index), A_left, 1e-6);
        TS_ASSERT_DELTA(vertex_mesh.GetSurfaceAreaOfElement(1 - left_element_index), A_right, 1e-6);

        // Test that the orientations are correct
        for (unsigned element_index = 0; element_index < 2; element_index++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(element_index);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(element_index);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);

                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }
    }

    void TestDivideElementAlongGivenAxisNonCoPlanarSeparation()
    {
        // We test this with a rectangle, where we divide by a vector, which is tilted by 10°
        /* Here we use an axis, where the basal and apical edges do not belong to the same
         * lateral face, which leads to minimization of the basal distance to the plane such that
         * the node is closest to the plane (with 1% minimal distance) on the same face.
         * Additionally the basal edge is very short such that the average normal of the face does
         * not separate the basal edges anymore. We have to use edge information for that
         *
         *                 o----o---o       o--o-o---o
         *                /| _ _ °./| =>   /|_:_°.°./|
         *               /,       / /     /, :     / /
         *              o----o--o/ /     o--o-o--o/ /
         *              |._°____| /      |.:°____| /
         */
        std::vector<Node<3>*> regular_nodes;
        // Basal
        regular_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(2, false, 2.0, 0.0, 0.0));
        regular_nodes.push_back(new Node<3>(3, false, 2.0, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(4, false, 1.0, 1.0, 0.0));
        regular_nodes.push_back(new Node<3>(5, false, 0.0, 1.0, 0.0));
        // Apical
        regular_nodes.push_back(new Node<3>(6, false, 0.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(7, false, 0.01, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(8, false, 2.0, 0.0, 3.0));
        regular_nodes.push_back(new Node<3>(9, false, 2.0, 1.0, 3.0));
        regular_nodes.push_back(new Node<3>(10, false, 1.99, 1.0, 3.0));
        regular_nodes.push_back(new Node<3>(11, false, 0.0, 1.0, 3.0));

        std::vector<MonolayerVertexElementType> regular_types(6, MonolayerVertexElementType::Basal);
        std::vector<MonolayerVertexElementType> regular_types_2(6, MonolayerVertexElementType::Apical);
        regular_types.insert(regular_types.end(), regular_types_2.begin(), regular_types_2.end());

        std::vector<std::vector<int> > face_indices;
        std::vector<int> face_0 = { 0, 5, 4, 3, 2, 1 };
        std::vector<int> face_1 = { 6, 7, 8, 9, 10, 11 };
        std::vector<int> face_2 = { 1, 7, 6, 0 };
        std::vector<int> face_3 = { 2, 8, 7, 1 };
        std::vector<int> face_4 = { 3, 9, 8, 2 };
        std::vector<int> face_5 = { 4, 10, 9, 3 };
        std::vector<int> face_6 = { 5, 11, 10, 4 };
        std::vector<int> face_7 = { 0, 6, 11, 5 };
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

        regular_elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, regular_faces, regular_faces_orientations, regular_nodes, regular_types));
        MutableMonolayerVertexMesh<3, 3> vertex_mesh(regular_nodes, regular_faces, regular_elements);

        // Check initial state
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 8u);

        // Create vector for normal of division plane
        c_vector<double, 3> normal_axis = zero_vector<double>(3);
        normal_axis[0] = cos(M_PI / 6.0);
        normal_axis[2] = -sin(M_PI / 6.0);
        normal_axis[1] = 0.0;

        vertex_mesh.DivideElementAlongGivenAxis(regular_elements[0], normal_axis);

        // Check final state
        unsigned left_element_index = vertex_mesh.GetVolumeOfElement(0) < vertex_mesh.GetVolumeOfElement(1) ? 0 : 1;

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumFaces(), 13u);

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(left_element_index)->GetNumFaces(), 6u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(left_element_index)->GetNumNodes(), 8u);
        MonolayerVertexElement<2, 3>* p_face_apical_left = vertex_mesh.GetFaceOfType(left_element_index, MonolayerVertexElementType::Apical);
        MonolayerVertexElement<2, 3>* p_face_basal_left = vertex_mesh.GetFaceOfType(left_element_index, MonolayerVertexElementType::Basal);
        TS_ASSERT_EQUALS(p_face_apical_left->GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(p_face_basal_left->GetNumNodes(), 4u);

        for (unsigned face_ind = 0; face_ind < 6; face_ind++)
        {
            if (vertex_mesh.GetElement(left_element_index)->GetFace(face_ind)->GetFaceType() != MonolayerVertexElementType::Lateral)
                continue;
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(left_element_index)->GetFace(face_ind)->GetNumNodes(), 4u);
        }

        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1 - left_element_index)->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(vertex_mesh.GetElement(1 - left_element_index)->GetNumNodes(), 12u);
        MonolayerVertexElement<2, 3>* p_face_apical_right = vertex_mesh.GetFaceOfType(1 - left_element_index, MonolayerVertexElementType::Apical);
        MonolayerVertexElement<2, 3>* p_face_basal_right = vertex_mesh.GetFaceOfType(1 - left_element_index, MonolayerVertexElementType::Basal);
        TS_ASSERT_EQUALS(p_face_apical_right->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_face_basal_right->GetNumNodes(), 6u);

        for (unsigned face_ind = 0; face_ind < 8; face_ind++)
        {
            if (vertex_mesh.GetElement(1 - left_element_index)->GetFace(face_ind)->GetFaceType() != MonolayerVertexElementType::Lateral)
                continue;
            TS_ASSERT_EQUALS(vertex_mesh.GetElement(1 - left_element_index)->GetFace(face_ind)->GetNumNodes(), 4u);
        }

        // Test that the orientations are correct
        for (unsigned element_index = 0; element_index < 2; element_index++)
        {
            MonolayerVertexElement<3, 3>* p_element = vertex_mesh.GetElement(element_index);
            c_vector<double, 3> centroid = vertex_mesh.GetCentroidOfElement(element_index);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);

                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = vertex_mesh.GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                vertex_mesh.CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = vertex_mesh.GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                vertex_mesh.CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }
    }

    void TestDivideElementInHoneycomb()
    {
        FiniteThicknessHoneycombVertexMeshGenerator generator(3, 3, false, 0.1, 0.1, 1.0);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        unsigned num_nodes_old = p_mesh->GetNumNodes();

        MonolayerVertexElement<3, 3>* p_element = p_mesh->GetElement(4);
        c_vector<double, 3> axis;
        axis[0] = cos(M_PI / 180.0) * cos(M_PI * 0.43);
        axis[0] = cos(M_PI / 180.0) * sin(M_PI * 0.43);
        axis[0] = sin(M_PI / 180.0);

        p_mesh->DivideElementAlongGivenAxis(p_element, axis);

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 10u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), num_nodes_old + 4u);
    }
};

#endif /*TESTMUTABLEMONOLAYERVERTEXMESHREMESH_HPP_*/
