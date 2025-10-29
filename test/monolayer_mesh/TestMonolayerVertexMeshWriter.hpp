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

#ifndef TESTMONOLAYERVERTEXMESHWRITER_HPP_
#define TESTMONOLAYERVERTEXMESHWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <fstream>
#include <string>

#include "MonolayerVertexMeshReader.hpp"
#include "MonolayerVertexMeshWriter.hpp"
#include "MutableMonolayerVertexMesh.hpp"
// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestVertexMeshWriter : public CxxTest::TestSuite
{
public:
    MutableMonolayerVertexMesh<3, 3>* ConstructTwoAdjacentCubes()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(6, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(8, false, 0.0, 2.0, 1.0));
        nodes.push_back(new Node<3>(9, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(10, false, 1.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(11, false, 1.0, 2.0, 1.0));

        // Set basal/apical nodetypes
        std::vector<MonolayerVertexElementType> node_types;
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);

        std::vector<std::vector<Node<3>*> > nodes_faces(11);
        std::vector<std::vector<MonolayerVertexElementType> > nodes_faces_types(11);

        unsigned node_indices_face_0[4] = { 0, 3, 4, 1 };
        unsigned node_indices_face_1[4] = { 0, 6, 9, 3 };
        unsigned node_indices_face_2[4] = { 0, 1, 7, 6 };
        unsigned node_indices_face_3[4] = { 1, 4, 10, 7 };
        unsigned node_indices_face_4[4] = { 4, 3, 9, 10 };
        unsigned node_indices_face_5[4] = { 6, 7, 10, 9 };
        unsigned node_indices_face_6[4] = { 1, 4, 5, 2 };
        unsigned node_indices_face_7[4] = { 1, 2, 8, 7 };
        unsigned node_indices_face_8[4] = { 2, 5, 11, 8 };
        unsigned node_indices_face_9[4] = { 5, 4, 10, 11 };
        unsigned node_indices_face_10[4] = { 7, 8, 11, 10 };

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

            nodes_faces[6].push_back(nodes[node_indices_face_6[i]]);
            nodes_faces_types[6].push_back(node_types[node_indices_face_6[i]]);

            nodes_faces[7].push_back(nodes[node_indices_face_7[i]]);
            nodes_faces_types[7].push_back(node_types[node_indices_face_7[i]]);

            nodes_faces[8].push_back(nodes[node_indices_face_8[i]]);
            nodes_faces_types[8].push_back(node_types[node_indices_face_8[i]]);

            nodes_faces[9].push_back(nodes[node_indices_face_9[i]]);
            nodes_faces_types[9].push_back(node_types[node_indices_face_9[i]]);

            nodes_faces[10].push_back(nodes[node_indices_face_10[i]]);
            nodes_faces_types[10].push_back(node_types[node_indices_face_10[i]]);
        }

        // Make the faces
        std::vector<MonolayerVertexElement<2, 3>*> faces;

        // Set face types
        std::vector<MonolayerVertexElementType> face_types;
        face_types.push_back(MonolayerVertexElementType::Apical);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Basal);
        face_types.push_back(MonolayerVertexElementType::Apical);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Basal);

        for (unsigned i = 0; i < 11; i++)
        {
            faces.push_back(new MonolayerVertexElement<2, 3>(i, face_types[i], nodes_faces[i], nodes_faces_types[i]));
        }

        // Make the elements
        std::vector<MonolayerVertexElement<2, 3>*> faces_element_0, faces_element_1;
        std::vector<bool> orientations_0, orientations_1;

        // 1st Cube element
        for (unsigned i = 0; i < 6; i++)
        {
            faces_element_0.push_back(faces[i]);
            orientations_0.push_back(true);
        }

        // 2nd Cube element
        for (unsigned i = 6; i < 11; i++)
        {
            faces_element_1.push_back(faces[i]);
            orientations_1.push_back(true);
        }
        faces_element_1.push_back(faces[3]);
        orientations_1.push_back(false);

        std::vector<MonolayerVertexElement<3, 3>*> elements;
        elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, faces_element_0, orientations_0));
        elements.push_back(new MonolayerVertexElement<3, 3>(1, MonolayerVertexElementType::Undetermined, faces_element_1, orientations_1));

        return new MutableMonolayerVertexMesh<3, 3>(nodes, faces, elements);
    }

    void TestMonolayerVertexMeshWriterIn3dWithFaces()
    {
        MutableMonolayerVertexMesh<3, 3>* p_mesh = ConstructTwoAdjacentCubes();

        MonolayerVertexMeshWriter<3, 3> writer("TestMonolayerVertexMeshWriter3d", "monolayer_vertex_mesh", false);

        writer.WriteFilesUsingMesh(*p_mesh);
    }
};

#endif /*TESTMONOLAYERVERTEXMESHWRITER_HPP_*/