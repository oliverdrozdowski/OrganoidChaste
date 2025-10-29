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

#ifndef TESTVERTEXMESH2_HPP_
#define TESTVERTEXMESH2_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "MonolayerVertexElement.hpp"
#include "MonolayerVertexMesh.hpp"
// This test is always run sequentially (never in parallel)
#include "UblasCustomFunctions.hpp"
#include "FakePetscSetup.hpp"

class TestMonolayerVertexMesh2 : public CxxTest::TestSuite
{
private:
    MonolayerVertexMesh<3, 3>* ConstructCubeMesh()
    {
        std::vector<Node<3>*> nodes;

        // construct inner 8 nodes
        nodes.push_back(new Node<3>(0, false, -0.5, -0.5, -0.5));
        nodes.push_back(new Node<3>(1, false, -0.5, 0.5, -0.5));
        nodes.push_back(new Node<3>(2, false, -0.5, 0.5, 0.5));
        nodes.push_back(new Node<3>(3, false, -0.5, -0.5, 0.5));
        nodes.push_back(new Node<3>(4, false, 0.5, -0.5, -0.5));
        nodes.push_back(new Node<3>(5, false, 0.5, -0.5, 0.5));
        nodes.push_back(new Node<3>(6, false, 0.5, 0.5, 0.5));
        nodes.push_back(new Node<3>(7, false, 0.5, 0.5, -0.5));

        // make apical faces
        std::vector<std::vector<Node<3>*> > nodes_apical_faces(6);
        std::vector<std::vector<MonolayerVertexElementType> > nodes_apical_faces_types(6);

        // faces for cube in negative x
        unsigned node_indices_face_0[4] = { 0, 1, 2, 3 }; // apical face -x
        unsigned node_indices_face_1[4] = { 4, 5, 6, 7 }; // apical face +x
        unsigned node_indices_face_2[4] = { 0, 3, 5, 4 }; // apical face -y
        unsigned node_indices_face_3[4] = { 1, 7, 6, 2 }; // apical face +y
        unsigned node_indices_face_4[4] = { 0, 4, 7, 1 }; // apical face -z
        unsigned node_indices_face_5[4] = { 2, 6, 5, 3 }; // apical face +z

        std::vector<MonolayerVertexElementType> node_types;
        for (int i = 0; i < 6; i++)
        {
            node_types.push_back(MonolayerVertexElementType::Apical);
        }

        for (unsigned i = 0; i < 4; i++)
        {
            nodes_apical_faces[0].push_back(nodes[node_indices_face_0[i]]);
            nodes_apical_faces_types[0].push_back(node_types[node_indices_face_0[i]]);

            nodes_apical_faces[1].push_back(nodes[node_indices_face_1[i]]);
            nodes_apical_faces_types[1].push_back(node_types[node_indices_face_1[i]]);

            nodes_apical_faces[2].push_back(nodes[node_indices_face_2[i]]);
            nodes_apical_faces_types[2].push_back(node_types[node_indices_face_2[i]]);

            nodes_apical_faces[3].push_back(nodes[node_indices_face_3[i]]);
            nodes_apical_faces_types[3].push_back(node_types[node_indices_face_3[i]]);

            nodes_apical_faces[4].push_back(nodes[node_indices_face_4[i]]);
            nodes_apical_faces_types[4].push_back(node_types[node_indices_face_4[i]]);

            nodes_apical_faces[5].push_back(nodes[node_indices_face_5[i]]);
            nodes_apical_faces_types[5].push_back(node_types[node_indices_face_5[i]]);
        }

        // Make the faces
        std::vector<MonolayerVertexElement<2, 3>*> faces;

        // Set face types
        std::vector<MonolayerVertexElementType> face_types;
        face_types.push_back(MonolayerVertexElementType::Apical);
        face_types.push_back(MonolayerVertexElementType::Apical);
        face_types.push_back(MonolayerVertexElementType::Apical);
        face_types.push_back(MonolayerVertexElementType::Apical);
        face_types.push_back(MonolayerVertexElementType::Apical);
        face_types.push_back(MonolayerVertexElementType::Apical);

        std::vector<MonolayerVertexElement<3, 3>*> elements;

        for (unsigned i = 0; i < 6; i++)
        {
            faces.push_back(new MonolayerVertexElement<2, 3>(i, face_types[i], nodes_apical_faces[i], nodes_apical_faces_types[i]));

            std::vector<MonolayerVertexElement<2, 3>*> faces_element;
            std::vector<bool> orientations;

            faces_element.push_back(faces[i]);
            orientations.push_back(true);

            elements.push_back(new MonolayerVertexElement<3, 3>(i, MonolayerVertexElementType::Undetermined, faces_element, orientations));
        }

        return new MonolayerVertexMesh<3, 3>(nodes, faces, elements);

        /*

        // cube in negative x direction
        nodes.push_back(new Node<3>(8, true, -1.5, -0.5, -0.5));
        nodes.push_back(new Node<3>(9, true, -1.5, 0.5, -0.5));
        nodes.push_back(new Node<3>(10, true, -1.5, 0.5, 0.5));
        nodes.push_back(new Node<3>(11, true, -1.5, -0.5, 0.5));

        // cube in positive x direction
        nodes.push_back(new Node<3>(12, true, 1.5, -0.5, -0.5));
        nodes.push_back(new Node<3>(13, true, 1.5, 0.5, -0.5));
        nodes.push_back(new Node<3>(14, true, 1.5, 0.5, 0.5));
        nodes.push_back(new Node<3>(15, true, 1.5, -0.5, 0.5));

        // cube in negative y direction
        nodes.push_back(new Node<3>(16, true, -0.5, -1.5, -0.5));
        nodes.push_back(new Node<3>(17, true, 0.5, -1.5, -0.5));
        nodes.push_back(new Node<3>(18, true, 0.5, -1.5, 0.5));
        nodes.push_back(new Node<3>(19, true, -0.5, -1.5, 0.5));

        // cube in positive y direction
        nodes.push_back(new Node<3>(20, true, -0.5, 1.5, -0.5));
        nodes.push_back(new Node<3>(21, true, 0.5, 1.5, -0.5));
        nodes.push_back(new Node<3>(22, true, 0.5, 1.5, 0.5));
        nodes.push_back(new Node<3>(23, true, -0.5, 1.5, 0.5));

        // cube in negative z direction
        nodes.push_back(new Node<3>(24, true, -0.5, -0.5, -1.5));
        nodes.push_back(new Node<3>(25, true, 0.5, -0.5, -1.5));
        nodes.push_back(new Node<3>(26, true, 0.5, 0.5, -1.5));
        nodes.push_back(new Node<3>(27, true, -0.5, 0.5, -1.5));

        // cube in positive z direction
        nodes.push_back(new Node<3>(28, true, -0.5, -0.5, 1.5));
        nodes.push_back(new Node<3>(29, true, 0.5, -0.5, 1.5));
        nodes.push_back(new Node<3>(30, true, 0.5, 0.5, 1.5));
        nodes.push_back(new Node<3>(31, true, -0.5, 0.5, 1.5));

        // Set apical nodetypes
        std::vector<MonolayerVertexElementType> node_types;

        for (int i = 0; i < 8; i++)
        {
            node_types.push_back(MonolayerVertexElementType::Apical);
        }

        for (int i = 0; i < 24; i++)
        {
            node_types.push_back(MonolayerVertexElementType::Basal);
        }

        // make the faces
        std::vector<std::vector<Node<3>*> > nodes_faces(6);
        std::vector<std::vector<MonolayerVertexElementType> > nodes_faces_types(6);

        // faces for cube in negative x
        unsigned node_indices_face_0[4] = { 0, 2, 6, 3 }; // apical face
        unsigned node_indices_face_1[4] = { 8, 11, 10, 9 }; // basal face
        unsigned node_indices_face_2[4] = { 9, 10, 6, 2 }; // lateral face
        unsigned node_indices_face_3[4] = { 0, 3, 11, 8 }; // lateral face
        unsigned node_indices_face_4[4] = { 0, 8, 9, 2 }; // lateral face
        unsigned node_indices_face_5[4] = { 10, 11, 3, 6 }; // lateral face

        // faces for cube in positive x
        unsigned node_indices_face_6[4] = { 1, 4, 7, 5 }; // apical face
        unsigned node_indices_face_7[4] = { 12, 13, 14, 15 }; // basal face
        unsigned node_indices_face_8[4] = { 1, 12, 15, 4 }; // lateral face
        unsigned node_indices_face_9[4] = { 7, 14, 13, 5 }; // lateral face
        unsigned node_indices_face_10[4] = { 0, 8, 9, 2 }; // lateral face
        unsigned node_indices_face_11[4] = { 10, 11, 3, 6 }; // lateral face
        */
    }

public:
    void TestGetLumenVolume()
    {
        MonolayerVertexMesh<3, 3>* p_mesh = ConstructCubeMesh();

        std::cout << "LUMEN VOLUME " << p_mesh->GetLumenVolume() << std::endl;

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetLumenVolume(), 1u);

        TS_ASSERT_EQUALS(p_mesh->GetLumenVolume(true), p_mesh->GetLumenVolume(false));
    }

    void TestGetLumenVolumeGradient()
    {
        MonolayerVertexMesh<3, 3>* p_mesh = ConstructCubeMesh();

        unsigned num_nodes = p_mesh->GetNumNodes();

        c_vector<double, 3> gradient1 = zero_vector<double>(3);
        c_vector<double, 3> gradient2 = zero_vector<double>(3);

        for (unsigned node_index = 0; node_index < num_nodes; node_index++)
        {
            std::cout << "NODE " << node_index << std::endl;
            gradient1 += p_mesh->CalculateLumenVolGradient(node_index, true);
            gradient2 += p_mesh->CalculateLumenVolGradient(node_index, false);
        }

        std::cout << "GRADIENT 1 " << gradient1 << std::endl;
        std::cout << "GRADIENT 2 " << gradient2 << std::endl;

        TS_ASSERT_DELTA(gradient1[0], gradient2[0], 1e-6);
    }
};
#endif /*TESTVERTEXMESH2_HPP_*/