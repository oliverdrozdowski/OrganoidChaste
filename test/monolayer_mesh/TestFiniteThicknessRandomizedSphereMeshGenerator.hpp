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

#ifndef TESTFINITETHICKNESSRANDOMIZEDSPHEREMESHGENERATOR_HPP_
#define TESTFINITETHICKNESSRANDOMIZEDSPHEREMESHGENERATOR_HPP_

#include <cmath> // for M_PI
#include <cxxtest/TestSuite.h>
#include "FiniteThicknessRandomizedSphereMeshGenerator.hpp"
#include "MonolayerVertexElement.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestFiniteThicknessRandomizedSphereMeshGenerator : public CxxTest::TestSuite
{
public:
    void TestSimpleSphere()
    {
        FiniteThicknessRandomizedSphereMeshGenerator generator(20000, 1.0, 1.0, 1.0, 1.0);

        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 20000u);
        TS_ASSERT_DELTA(p_mesh->GetCellRearrangementThreshold(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetT2Threshold(), 1.0, 1e-12);

        double total_volume = 0.0;
        for (unsigned index = 0; index < 20000; index++)
        {
            TS_ASSERT_DELTA(p_mesh->GetThicknessOfElement(index), 1.0, 1e-3); // Rounding errors from imperfect sphere
            total_volume += p_mesh->GetVolumeOfElement(index);
        }
        TS_ASSERT_DELTA(total_volume, 4.0 / 3.0 * M_PI * 7.0, 1e-1); // Rounding errors from imperfect sphere
    }

    void TestFacesThickSphere()
    {
        FiniteThicknessRandomizedSphereMeshGenerator generator(20000, 1.0, 1.0, 1.0, 0.2);

        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 20000u);
        TS_ASSERT_DELTA(p_mesh->GetCellRearrangementThreshold(), 1.0, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetT2Threshold(), 1.0, 1e-12);

        double total_volume = 0.0;
        for (unsigned index = 0; index < 20000; index++)
        {
            TS_ASSERT_DELTA(p_mesh->GetThicknessOfElement(index), 1.0, 1e-3); // Rounding errors from imperfect sphere
            total_volume += p_mesh->GetVolumeOfElement(index);

            for (unsigned face_index = 0; face_index < p_mesh->GetElement(index)->GetNumFaces(); face_index++)
            {
                MonolayerVertexElement<2, 3>* pFace = p_mesh->GetElement(index)->GetFace(face_index);
                unsigned num_nodes = pFace->GetNumNodes();
                TS_ASSERT(num_nodes >= 3u);
            }
        }
        TS_ASSERT_DELTA(total_volume, 4.0 / 3.0 * M_PI * 1.72, 1e-1); // Rounding errors from imperfect sphere
    }

    void TestBoundaryNodes()
    {
        FiniteThicknessRandomizedSphereMeshGenerator generator(20);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 20u);

        unsigned num_boundary_nodes = 0;
        // Check number of non-boundary nodes in first 4 elements
        for (unsigned element_index = 0; element_index < 20u; element_index++)
        {
            unsigned num_nodes = p_mesh->GetElement(element_index)->GetNumNodes();
            for (unsigned node_index = 0; node_index < num_nodes; node_index++)
            {
                if (p_mesh->GetElement(element_index)->GetNode(node_index)->IsBoundaryNode())
                {
                    num_boundary_nodes++;
                }
            }
        }
        // We should have no boundary nodes, as a sphere is closed
        TS_ASSERT_EQUALS(num_boundary_nodes, 0u);
    }

    void TestFaces()
    {
        FiniteThicknessRandomizedSphereMeshGenerator generator(1000);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 1000u);

        for (unsigned element = 0; element < 1000u; element++)
        {
            unsigned num_faces = p_mesh->GetElement(element)->GetNumFaces();
            unsigned apical_num = 0;
            unsigned basal_num = 0;
            unsigned lateral_num = 0;
            for (unsigned face = 0; face < num_faces; face++)
            {
                switch (p_mesh->GetElement(element)->GetFace(face)->GetFaceType())
                {
                    case MonolayerVertexElementType::Apical:
                        apical_num++;
                        break;
                    case MonolayerVertexElementType::Basal:
                        basal_num++;
                        break;
                    case MonolayerVertexElementType::Lateral:
                        lateral_num++;
                        break;
                    default: // This should NOT happen
                        bool found_wrong_face_type = false;
                        TS_ASSERT(found_wrong_face_type);
                }
            }
            TS_ASSERT_EQUALS(apical_num, 1u);
            TS_ASSERT_EQUALS(basal_num, 1u);
            TS_ASSERT_EQUALS(lateral_num, num_faces - 2);
            unsigned num_nodes = p_mesh->GetElement(element)->GetNumNodes() / 2;
            TS_ASSERT_EQUALS(lateral_num, num_nodes);
        }
    }

    void TestOrientations()
    {
        FiniteThicknessRandomizedSphereMeshGenerator generator(1000);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 1000u);

        // Test that the orientations are correct prior to remeshing
        // Consistency check to see if mesh is correctly defined
        for (unsigned i = 0; i < p_mesh->GetNumElements(); i++)
        {
            MonolayerVertexElement<3, 3>* p_element = p_mesh->GetElement(i);
            c_vector<double, 3> centroid = p_mesh->GetCentroidOfElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                double orientation_factor = p_element->FaceIsOrientatedClockwise(j_face) ? -1.0 : 1.0;
                c_vector<double, 3> passive_center = p_mesh->GetPassiveCenterOfFace(p_face);
                std::vector<c_vector<double, 3> > normal_vectors;
                p_mesh->CalculateAreasAndNormalsOfTriangulationsOfFace(p_face, &normal_vectors);
                passive_center = p_mesh->GetVectorFromAtoB(centroid, passive_center);

                // Check anti-clockwise orientation of nodes
                for (unsigned node_index = 0; node_index < normal_vectors.size() - 1; node_index++)
                {
                    TS_ASSERT(inner_prod(normal_vectors[node_index], normal_vectors[node_index + 1]) > 0.0);
                }
                // Check for correct orientation of face in element
                c_vector<double, 3> normal_vector;
                p_mesh->CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                TS_ASSERT(orientation_factor * inner_prod(normal_vector, passive_center) > 0.0);
            }
        }

        // Apical and basal faces have to be correctly oriented.
        for (unsigned i = 0; i < p_mesh->GetNumElements(); i++)
        {
            MonolayerVertexElement<3, 3>* p_element = p_mesh->GetElement(i);
            for (unsigned j_face = 0; j_face < p_element->GetNumFaces(); j_face++)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(j_face);
                MonolayerVertexElementType face_type = p_face->GetFaceType();
                if (face_type == MonolayerVertexElementType::Apical || face_type == MonolayerVertexElementType::Basal)
                {
                    // Apical and basal faces should be oriented AntiClockwise, viewed from outside
                    TS_ASSERT(!(p_element->FaceIsOrientatedClockwise(j_face)));
                    c_vector<double, 3> normal_vector;
                    p_mesh->CalculateUnitNormalToFaceWithArea(p_face, normal_vector);
                    c_vector<double, 3> passive_center = p_mesh->GetPassiveCenterOfFace(p_face);
                    if (face_type == MonolayerVertexElementType::Apical)
                    {
                        // Apical faces have normals pointing inwards
                        TS_ASSERT(inner_prod(passive_center, normal_vector) < 0.0);
                    }
                    else
                    {
                        // Basal faces have normals pointing outwards
                        TS_ASSERT(inner_prod(passive_center, normal_vector) > 0.0);
                    }
                }
            }
        }
    }

    void TestInternalFunctions()
    {
        FiniteThicknessRandomizedSphereMeshGenerator generator(200, 1.0, 1.0, 1.0, 1.0);

        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        for (unsigned face = 0; face < p_mesh->GetNumFaces(); face++)
        {
            MonolayerVertexElement<2, 3>* p_face = p_mesh->GetFace(face);

            // Check if normal of face matches to result from mesh
            c_vector<double, 3> normal_generator = generator.GetNormalToFace(p_face);
            c_vector<double, 3> normal_mesh;
            p_mesh->CalculateUnitNormalToFaceWithArea(p_face, normal_mesh);

            TS_ASSERT_DELTA(normal_generator[0], normal_mesh[0], 1e-6);
            TS_ASSERT_DELTA(normal_generator[1], normal_mesh[1], 1e-6);
            TS_ASSERT_DELTA(normal_generator[2], normal_mesh[2], 1e-6);

            // Check if passive center of face matches to result from mesh
            c_vector<double, 3> center_generator = generator.GetCenterOfFace(p_face);
            c_vector<double, 3> center_mesh = p_mesh->GetPassiveCenterOfFace(p_face);

            TS_ASSERT_DELTA(center_generator[0], center_mesh[0], 1e-6);
            TS_ASSERT_DELTA(center_generator[1], center_mesh[1], 1e-6);
            TS_ASSERT_DELTA(center_generator[2], center_mesh[2], 1e-6);
        }
    }
};

#endif /*TESTFINITETHICKNESSRANDOMIZEDSPHEREMESHGENERATOR_HPP_*/