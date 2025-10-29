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

#ifndef TESTVORONOISPHEREGENERATOR_HPP_
#define TESTVORONOISPHEREGENERATOR_HPP_

#include <boost/numeric/ublas/assignment.hpp>
#include <cxxtest/TestSuite.h>
#include <set>
#include "UblasCustomFunctions.hpp"
#include "VoronoiSphereGenerator.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestVoronoiSphereGenerator : public CxxTest::TestSuite
{
public:
    void TestDoublePyramidSphere()
    {
        /* Generate a double pyramid with vertices on the sphere
         *			 North
         *       .o.
         *     .  |  .
         *   .  .:o:.  .
         *   o.   |   .o
         *    . :.o.: .
         *     .  |  .
         *       .o.
         *      South
         */
        std::vector<c_vector<double, 3> > inputPositions;
        inputPositions.push_back(c_vector<double, 3>());
        inputPositions.push_back(c_vector<double, 3>());
        inputPositions.push_back(c_vector<double, 3>());
        inputPositions.push_back(c_vector<double, 3>());
        inputPositions.push_back(c_vector<double, 3>());
        inputPositions.push_back(c_vector<double, 3>());

        inputPositions[0] <<= 0.0, 0.0, 1.0;
        inputPositions[1] <<= -1.0, 0.0, 0.0;
        inputPositions[2] <<= 0.0, -1.0, 0.0;
        inputPositions[3] <<= 1.0, 0.0, 0.0;
        inputPositions[4] <<= 0.0, 1.0, 0.0;
        inputPositions[5] <<= 0.0, 0.0, -1.0;

        VoronoiSphereGenerator spheregenerator(inputPositions);

        c_vector<double, 3> voronoi_0, voronoi_1, voronoi_2, voronoi_3, voronoi_4, voronoi_5, voronoi_6, voronoi_7;
        voronoi_0 <<= -1.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0;
        voronoi_1 <<= 1.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0;
        voronoi_2 <<= 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0;
        voronoi_3 <<= -1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0;
        voronoi_4 <<= -1.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0;
        voronoi_5 <<= 1.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0;
        voronoi_6 <<= 1.0 / 3.0, 1.0 / 3.0, -1.0 / 3.0;
        voronoi_7 <<= -1.0 / 3.0, 1.0 / 3.0, -1.0 / 3.0;

        std::vector<c_vector<double, 3> > voronoi_vertices = { voronoi_0, voronoi_1, voronoi_2, voronoi_3, voronoi_4, voronoi_5,
                                                               voronoi_6, voronoi_7 };

        std::vector<c_vector<double, 3> > sphere_vertices = spheregenerator.GetVoronoiVertices();

        // Check if voronoi vertices match
        for (auto it = sphere_vertices.begin(); it != sphere_vertices.end(); ++it)
        {
            bool found_elt = false;
            for (auto it_inner = voronoi_vertices.begin(); it_inner != voronoi_vertices.end(); ++it_inner)
            {
                double diff_norm = norm_2((*it_inner) - (*it));
                found_elt = diff_norm < 1e-8;
                if (found_elt)
                    break;
            }
            TS_ASSERT(found_elt);
        }
        TS_ASSERT_EQUALS(sphere_vertices.size(), 8u);
        TS_ASSERT_EQUALS(voronoi_vertices.size(), 8u);

        /*****************
         * Check the faces
         */
        // Create set of analytical and computed faces
        std::set<std::set<unsigned> > set_face_indices_sets;
        // side faces
        set_face_indices_sets.insert(std::set<unsigned>{ 0, 1, 4, 5 });
        set_face_indices_sets.insert(std::set<unsigned>{ 1, 2, 5, 6 });
        set_face_indices_sets.insert(std::set<unsigned>{ 2, 3, 6, 7 });
        set_face_indices_sets.insert(std::set<unsigned>{ 0, 3, 4, 7 });
        // top face
        set_face_indices_sets.insert(std::set<unsigned>{ 0, 1, 2, 3 });
        // bottom face
        set_face_indices_sets.insert(std::set<unsigned>{ 4, 5, 6, 7 });

        // Now we remap the indices to the ordering found by the algorithm
        std::map<unsigned, unsigned> map_analytical_to_sphere;
        unsigned index_ana = 0;
        for (auto it_ana = voronoi_vertices.begin(); it_ana != voronoi_vertices.end(); ++it_ana)
        {
            unsigned index_sph = 0;
            for (auto it_sph = sphere_vertices.begin(); it_sph != sphere_vertices.end(); ++it_sph)
            {
                double diff_norm = norm_2((*it_ana) - (*it_sph));
                if (diff_norm < 1e-8)
                {
                    map_analytical_to_sphere[index_ana] = index_sph;
                    break;
                }
                index_sph++;
            }
            index_ana++;
        }
        std::set<std::set<unsigned> > set_face_temp;
        for (auto face = set_face_indices_sets.begin(); face != set_face_indices_sets.end(); ++face)
        {
            std::set<unsigned> face_indices;
            for (auto ind = face->begin(); ind != face->end(); ++ind)
            {
                face_indices.insert(map_analytical_to_sphere[*ind]);
            }
            set_face_temp.insert(face_indices);
        }
        set_face_indices_sets = set_face_temp;

        // Get the computed faces from algorithm
        std::map<unsigned, std::vector<unsigned> > sphere_faces = spheregenerator.GetMapFromFacesToVertices();
        std::set<std::set<unsigned> > set_sphere_face_indices;
        for (std::map<unsigned, std::vector<unsigned> >::iterator it = sphere_faces.begin(); it != sphere_faces.end(); ++it)
        {
            set_sphere_face_indices.insert(std::set<unsigned>(it->second.begin(), it->second.end()));
        }

        // Check if face indices match
        for (auto it = set_sphere_face_indices.begin(); it != set_sphere_face_indices.end(); ++it)
        {
            bool found_face = false;
            for (auto it_inner = set_face_indices_sets.begin(); it_inner != set_face_indices_sets.end(); ++it_inner)
            {
                found_face = ((*it_inner) == (*it));
                if (found_face)
                    break;
            }
            TS_ASSERT(found_face);
        }
        TS_ASSERT_EQUALS(set_face_indices_sets.size(), 6u);
        TS_ASSERT_EQUALS(sphere_faces.size(), 6u);

        /********
         * Check the face ordering
         */
        for (auto face_pair = sphere_faces.begin(); face_pair != sphere_faces.end(); ++face_pair)
        {
            // Calculate passive center
            std::vector<unsigned> face = face_pair->second;
            unsigned num_vertices = face.size();
            c_vector<double, 3> passive_center = zero_vector<double>(3);
            for (auto index = face.begin(); index != face.end(); ++index)
            {
                passive_center += sphere_vertices[*index];
            }
            passive_center /= (double)num_vertices;

            // Calculate normal for each triangle from triangulation
            c_vector<double, 3> to_current = sphere_vertices[face[0]] - passive_center;
            c_vector<double, 3> normal_prev;
            for (unsigned vertex = 0; vertex < num_vertices + 1; vertex++)
            {
                c_vector<double, 3> to_next = sphere_vertices[face[(vertex + 1) % num_vertices]] - passive_center;
                // this normal
                c_vector<double, 3> normal_this = VectorProduct(to_current, to_next);
                if (vertex > 0)
                {
                    // Check if normals point in same direction, which is a good check for ordering
                    TS_ASSERT(inner_prod(normal_prev, normal_this) > 0.0);
                }
                normal_prev = normal_this;
                to_current = to_next;
                // Check triangle area, which should be 4/9*1/4
                TS_ASSERT_DELTA(norm_2(normal_this), 4.0 / 9.0 * 2.0 / 4.0, 1e-6);
            }
        }
    }
};

#endif /*TESTVORONOISPHEREGENERATOR_HPP_*/