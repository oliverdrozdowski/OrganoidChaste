/*
Copyright (c) 2005-2019, University of Oxford.
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

#ifndef TESTFINITETHICKNESSHONEYCOMBVERTEXMESHGENERATOR_HPP_
#define TESTFINITETHICKNESSHONEYCOMBVERTEXMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "FiniteThicknessHoneycombVertexMeshGenerator.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestFiniteThicknessHoneycombVertexMeshGenerator : public CxxTest::TestSuite
{
public:
    void TestSimpleMesh()
    {
        FiniteThicknessHoneycombVertexMeshGenerator generator(2, 2, false, 0.1, 0.1, 1.0);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 32u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 4u);
        TS_ASSERT_DELTA(p_mesh->GetCellRearrangementThreshold(), 0.1, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetT2Threshold(), 0.1, 1e-12);
    }

    void TestBoundaryNodes()
    {
        FiniteThicknessHoneycombVertexMeshGenerator generator(4, 4);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 96u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 16u);

        unsigned num_non_boundary_nodes = 0;
        // Check number of non-boundary nodes in first 4 elements
        for (unsigned element_index = 0; element_index < 4u; element_index++)
        {
            TS_ASSERT_EQUALS(p_mesh->GetElement(element_index)->GetNumNodes(), 12u);
            for (unsigned node_index = 0; node_index < 12u; node_index++)
            {
                if (!p_mesh->GetElement(element_index)->GetNode(node_index)->IsBoundaryNode())
                {
                    num_non_boundary_nodes++;
                }
            }
        }
        // Note that we count interior nodes that belong to two elements twice
        TS_ASSERT_EQUALS(num_non_boundary_nodes, 18u);
    }

    void TestLargeMesh()
    {
        FiniteThicknessHoneycombVertexMeshGenerator generator(100, 100, false, 0.01, 0.001, 1.0, 1.0);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 40800u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 10000u);

        // Check interior element for correct positioning etc.
        unsigned element_number = 322;
        MonolayerVertexElement<3, 3>* p_interior_element = p_mesh->GetElement(element_number);

        // Check positions
        // Apical side
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(0)->rGetLocation()[0], (1.0 + 22.0) * sqrt(2.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(0)->rGetLocation()[1], (9.0) * sqrt(1.0 / 6.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(0)->rGetLocation()[2], 1.0, 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(1)->rGetLocation()[0], (1.5 + 22.0) * sqrt(2.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(1)->rGetLocation()[1], (10.0) * sqrt(1.0 / 6.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(1)->rGetLocation()[2], 1.0, 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(2)->rGetLocation()[0], (1.5 + 22.0) * sqrt(2.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(2)->rGetLocation()[1], (12.0) * sqrt(1.0 / 6.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(2)->rGetLocation()[2], 1.0, 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(3)->rGetLocation()[0], (1.0 + 22.0) * sqrt(2.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(3)->rGetLocation()[1], (13.0) * sqrt(1.0 / 6.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(3)->rGetLocation()[2], 1.0, 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(4)->rGetLocation()[0], (0.5 + 22.0) * sqrt(2.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(4)->rGetLocation()[1], (12.0) * sqrt(1.0 / 6.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(4)->rGetLocation()[2], 1.0, 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(5)->rGetLocation()[0], (0.5 + 22.0) * sqrt(2.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(5)->rGetLocation()[1], (10.0) * sqrt(1.0 / 6.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(0)->GetNode(5)->rGetLocation()[2], 1.0, 1e-10);

        // Basal side
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(0)->rGetLocation()[0], (1.0 + 22.0) * sqrt(2.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(0)->rGetLocation()[1], (9.0) * sqrt(1.0 / 6.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(0)->rGetLocation()[2], 0.0, 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(1)->rGetLocation()[0], (0.5 + 22.0) * sqrt(2.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(1)->rGetLocation()[1], (10.0) * sqrt(1.0 / 6.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(1)->rGetLocation()[2], 0.0, 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(2)->rGetLocation()[0], (0.5 + 22.0) * sqrt(2.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(2)->rGetLocation()[1], (12.0) * sqrt(1.0 / 6.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(2)->rGetLocation()[2], 0.0, 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(3)->rGetLocation()[0], (1.0 + 22.0) * sqrt(2.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(3)->rGetLocation()[1], (13.0) * sqrt(1.0 / 6.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(3)->rGetLocation()[2], 0.0, 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(4)->rGetLocation()[0], (1.5 + 22.0) * sqrt(2.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(4)->rGetLocation()[1], (12.0) * sqrt(1.0 / 6.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(4)->rGetLocation()[2], 0.0, 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(5)->rGetLocation()[0], (1.5 + 22.0) * sqrt(2.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(5)->rGetLocation()[1], (10.0) * sqrt(1.0 / 6.0 / sqrt(3.0)), 1e-10);
        TS_ASSERT_DELTA(p_interior_element->GetFace(1)->GetNode(5)->rGetLocation()[2], 0.0, 1e-10);

        // Check volume
        TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(element_number), 1.0, 1e-10);

        // Check faces
        TS_ASSERT_EQUALS(p_interior_element->GetNumFaces(), 8u);
        TS_ASSERT(p_interior_element->GetFace(0)->GetFaceType() == MonolayerVertexElementType::Apical);
        TS_ASSERT(p_interior_element->GetFace(1)->GetFaceType() == MonolayerVertexElementType::Basal);
        for (unsigned face_index = 2; face_index < 8; face_index++)
        {
            TS_ASSERT(p_interior_element->GetFace(face_index)->GetFaceType() == MonolayerVertexElementType::Lateral);
        }

        // Since mesh is built from bottom up, only the apical, basal and upper right
        // faces are oriented anticlockwise (i.e. normal pointing outside) for interior element
        TS_ASSERT(!p_interior_element->FaceIsOrientatedClockwise(0)); // Apical
        TS_ASSERT(!p_interior_element->FaceIsOrientatedClockwise(1)); // Basal
        TS_ASSERT(p_interior_element->FaceIsOrientatedClockwise(2)); // Lateral 0
        TS_ASSERT(!p_interior_element->FaceIsOrientatedClockwise(3));
        TS_ASSERT(!p_interior_element->FaceIsOrientatedClockwise(4));
        TS_ASSERT(!p_interior_element->FaceIsOrientatedClockwise(5));
        TS_ASSERT(p_interior_element->FaceIsOrientatedClockwise(6));
        TS_ASSERT(p_interior_element->FaceIsOrientatedClockwise(7));

        // Normals point outwards for correctly oriented elements
        c_vector<double, 3> normal_0;
        double area_0 = p_mesh->CalculateUnitNormalToFaceWithArea(p_interior_element->GetFace(0), normal_0);
        TS_ASSERT_DELTA(normal_0[0], 0.0, 1e-10);
        TS_ASSERT_DELTA(normal_0[1], 0.0, 1e-10);
        TS_ASSERT_DELTA(normal_0[2], 1.0, 1e-10);
        TS_ASSERT_DELTA(area_0, 1.0, 1e-10);

        c_vector<double, 3> normal_1;
        double area_1 = p_mesh->CalculateUnitNormalToFaceWithArea(p_interior_element->GetFace(1), normal_1);
        TS_ASSERT_DELTA(normal_1[0], 0.0, 1e-10);
        TS_ASSERT_DELTA(normal_1[1], 0.0, 1e-10);
        TS_ASSERT_DELTA(normal_1[2], -1.0, 1e-10);
        TS_ASSERT_DELTA(area_1, 1.0, 1e-10);

        // lateral face 0 points inwards and is oriented clockwise (from the outside)
        c_vector<double, 3> normal_2;
        double area_2 = p_mesh->CalculateUnitNormalToFaceWithArea(p_interior_element->GetFace(2), normal_2);
        TS_ASSERT_DELTA(normal_2[0], -1.0 / 2.0, 1e-10);
        TS_ASSERT_DELTA(normal_2[1], sqrt(3.0) / 2.0, 1e-10);
        TS_ASSERT_DELTA(normal_2[2], 0.0, 1e-10);
        TS_ASSERT_DELTA(area_2, sqrt(2.0 / 3.0 / sqrt(3.0)), 1e-10);

        // lateral face 0 points outwards and is oriented anticlockwise (from the outside)
        c_vector<double, 3> normal_3;
        double area_3 = p_mesh->CalculateUnitNormalToFaceWithArea(p_interior_element->GetFace(3), normal_3);
        TS_ASSERT_DELTA(normal_3[0], 1.0, 1e-10);
        TS_ASSERT_DELTA(normal_3[1], 0.0, 1e-10);
        TS_ASSERT_DELTA(normal_3[2], 0.0, 1e-10);
        TS_ASSERT_DELTA(area_3, sqrt(2.0 / 3.0 / sqrt(3.0)), 1e-10);

        // Check total surface area
        TS_ASSERT_DELTA(p_mesh->GetSurfaceAreaOfElement(22), 2.0 + 6.0 * sqrt(2.0 / 3.0 / sqrt(3.0)), 1e-10);
    }

    void TestElementArea()
    {
        FiniteThicknessHoneycombVertexMeshGenerator generator(6, 6, false, 0.01, 0.001, 0.88, 2.456);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 192u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 36u);

        for (unsigned elem_index = 0; elem_index < p_mesh->GetNumElements(); elem_index++)
        {
            TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(elem_index), 0.88 * 2.456, 1e-3);
            TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_mesh->GetElement(elem_index)->GetFace(0)), 2.456, 1e-3);
            TS_ASSERT_DELTA(p_mesh->CalculateAreaOfFace(p_mesh->GetElement(elem_index)->GetFace(1)), 2.456, 1e-3);
        }
    }

    void TestOrientations()
    {
        FiniteThicknessHoneycombVertexMeshGenerator generator(100, 100, false, 0.01, 0.001, 1.0, 1.0);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 40800u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 10000u);

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
    }

    void TestNodeAndFaceAndElementNumbers()
    {
        FiniteThicknessHoneycombVertexMeshGenerator generator(7, 7, false, 0.01, 0.001, 0.88, 2.456);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        for (unsigned i_elt = 0; i_elt < 6 * 6; ++i_elt)
        {
            MonolayerVertexElement<3, 3>* p_elt = p_mesh->GetElement(i_elt);

            TS_ASSERT(p_elt->GetNumNodes() == 12);
            TS_ASSERT(p_elt->GetNumFaces() == 8);
            TS_ASSERT(p_elt->GetIndex() == i_elt);
        }
        for (unsigned i_nd = 0; i_nd < p_mesh->GetNumNodes(); ++i_nd)
        {
            TS_ASSERT(p_mesh->GetNode(i_nd)->GetIndex() == i_nd);
        }

        // Check correct node registration
        for (unsigned i_nd = 0; i_nd < p_mesh->GetNumNodes(); ++i_nd)
        {
            auto set_indices = p_mesh->GetNode(i_nd)->rGetContainingElementIndices();
            for (auto it_elt_ind = set_indices.begin(); it_elt_ind != set_indices.end(); ++it_elt_ind)
            {
                MonolayerVertexElement<3, 3>* p_elt = p_mesh->GetElement(*it_elt_ind);
                TS_ASSERT(p_elt->GetNodeLocalIndex(i_nd) != UINT_MAX);
            }
        }
        for (unsigned i_elt = 0; i_elt < p_mesh->GetNumElements(); ++i_elt)
        {
            MonolayerVertexElement<3, 3>* p_elt = p_mesh->GetElement(i_elt);
            for (unsigned i_nd = 0; i_nd < p_elt->GetNumNodes(); ++i_nd)
            {
                auto set_indices = p_elt->GetNode(i_nd)->rGetContainingElementIndices();
                TS_ASSERT(set_indices.find(p_elt->GetIndex()) != set_indices.end());
            }
        }
    }

    void TestNodeAndFaceAndElementNumbersForCylinder()
    {
        FiniteThicknessHoneycombVertexMeshGenerator generator(7, 7, false, 0.01, 0.001, 0.88, 2.456);
        generator.MakeCylindrical();
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        for (unsigned i_elt = 0; i_elt < 6 * 6; ++i_elt)
        {
            MonolayerVertexElement<3, 3>* p_elt = p_mesh->GetElement(i_elt);

            TS_ASSERT(p_elt->GetNumNodes() == 12);
            TS_ASSERT(p_elt->GetNumFaces() == 8);
            TS_ASSERT(p_elt->GetIndex() == i_elt);
        }
        for (unsigned i_nd = 0; i_nd < p_mesh->GetNumNodes(); ++i_nd)
        {
            TS_ASSERT(p_mesh->GetNode(i_nd)->GetIndex() == i_nd);
        }

        // Check correct node registration
        for (unsigned i_nd = 0; i_nd < p_mesh->GetNumNodes(); ++i_nd)
        {
            auto set_indices = p_mesh->GetNode(i_nd)->rGetContainingElementIndices();
            for (auto it_elt_ind = set_indices.begin(); it_elt_ind != set_indices.end(); ++it_elt_ind)
            {
                MonolayerVertexElement<3, 3>* p_elt = p_mesh->GetElement(*it_elt_ind);
                TS_ASSERT(p_elt->GetNodeLocalIndex(i_nd) != UINT_MAX);
            }
        }
        for (unsigned i_elt = 0; i_elt < p_mesh->GetNumElements(); ++i_elt)
        {
            MonolayerVertexElement<3, 3>* p_elt = p_mesh->GetElement(i_elt);
            for (unsigned i_nd = 0; i_nd < p_elt->GetNumNodes(); ++i_nd)
            {
                auto set_indices = p_elt->GetNode(i_nd)->rGetContainingElementIndices();
                TS_ASSERT(set_indices.find(p_elt->GetIndex()) != set_indices.end());
            }
        }

        // Check boundary faces
        for (unsigned i_fc = 0; i_fc < p_mesh->GetNumFaces(); ++i_fc)
        {
            unsigned num_elements = 0;
            for (unsigned i_elt = 0; i_elt < p_mesh->GetNumElements(); ++i_elt)
            {
                unsigned local_index = p_mesh->GetElement(i_elt)->GetFaceLocalIndex(p_mesh->GetFace(i_fc));
                num_elements += local_index != UINT_MAX ? 1 : 0;
            }
            if (p_mesh->GetFace(i_fc)->GetFaceType() != MonolayerVertexElementType::Lateral)
            {
                continue;
            }
            if (num_elements == 1)
            {
                TS_ASSERT(p_mesh->GetFace(i_fc)->IsBoundaryFace());
            }
            if (num_elements != 1)
            {
                TS_ASSERT(!p_mesh->GetFace(i_fc)->IsBoundaryFace());
            }
        }
    }
};

#endif /*TESTFINITETHICKNESSHONEYCOMBVERTEXMESHGENERATOR_HPP_*/