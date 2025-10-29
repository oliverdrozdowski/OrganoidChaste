#include "VoronoiSphereGenerator.hpp"
#include <boost/numeric/ublas/assignment.hpp>
#include "UblasCustomFunctions.hpp"

// CGAL library
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/number_utils.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

VoronoiSphereGenerator::VoronoiSphereGenerator(std::vector<c_vector<double, 3> > faceCenters)
        : mInputFaceCenters(faceCenters)
{
    this->ConstructDelaunayTriangulation();
    this->CalculateVoronoiVertices();
    this->OrderVoronoiVertexMap();
}

void VoronoiSphereGenerator::ConstructDelaunayTriangulation()
{
    // We use CGAL to generate the convex hull via the quickhull algorithm
    std::vector<Point_3> points;
    for (typename std::vector<c_vector<double, 3> >::iterator it = this->mInputFaceCenters.begin();
         it != this->mInputFaceCenters.end();
         ++it)
    {
        points.push_back(Point_3((*it)[0], (*it)[1], (*it)[2]));
    }
    // define surface mesh object to hold convex hull
    Surface_mesh poly;
    // compute convex hull of non-collinear points
    CGAL::convex_hull_3(points.begin(), points.end(), poly);

    // Save the convex hull points to a vector of c_vectors
    this->mDelaunayVertices.clear();
    // Save a vector of triangles, which contain the convex hull points
    this->mDelaunayTriangleIndices.clear();

    // First add all the convex hull points
    unsigned index = 0;
    for (Surface_mesh::Vertex_range::iterator vertex = poly.vertices().begin(); vertex != poly.vertices().end(); ++vertex)
    {
        c_vector<double, 3> poly_vertex;
        poly_vertex <<= CGAL::to_double(poly.point(*vertex).x()), CGAL::to_double(poly.point(*vertex).y()), CGAL::to_double(poly.point(*vertex).z());
        this->mDelaunayVertices.push_back(poly_vertex);
        // initialize the map from faces (Delaunay vertices) to voronoi vertices
        this->mMapFaceToVertices[index] = std::vector<unsigned>();
        index++;
    }

    // Now go through all triangles...
    for (Surface_mesh::Face_range::iterator face = poly.faces().begin(); face != poly.faces().end(); ++face)
    {
        Surface_mesh::Halfedge_index halfedge_index = poly.halfedge(*face);
        size_t num_edges = poly.degree(*face);
        std::array<unsigned, 3> face_indices;
        assert(num_edges == 3);
        // ... and save all indices that belong to each one of them
        for (unsigned edge_index = 0; edge_index < num_edges; ++edge_index)
        {
            Surface_mesh::Vertex_index vertex_index = poly.target(halfedge_index);
            face_indices[edge_index] = (unsigned)vertex_index;
            halfedge_index = poly.next(halfedge_index); // go to next edge in face
        }
        this->mDelaunayTriangleIndices.push_back(face_indices);
    }
}

void VoronoiSphereGenerator::CalculateVoronoiVertices()
{
    /**
     * The voronoi vertices are calculated as the circumcenters of
     * the triangles of the Delaunay triangulation
     * The index of the voronoi vertex in mVoronoiVertices and
     * of the delaunay triangle in mDelaunayVertices are identical
     */
    unsigned index = 0; // voronoi vertex index
    for (typename std::vector<std::array<unsigned, 3> >::iterator it = this->mDelaunayTriangleIndices.begin();
         it != this->mDelaunayTriangleIndices.end();
         ++it)
    {
        c_vector<double, 3> vertex_a = this->mDelaunayVertices[(*it)[0]];
        c_vector<double, 3> vertex_b = this->mDelaunayVertices[(*it)[1]];
        c_vector<double, 3> vertex_c = this->mDelaunayVertices[(*it)[2]];

        c_vector<double, 3> line_ab = vertex_b - vertex_a;
        c_vector<double, 3> line_ac = vertex_c - vertex_a;
        c_vector<double, 3> line_bc = vertex_c - vertex_b;

        c_vector<double, 3> ab_cross_ac = VectorProduct(line_ab, line_ac);

        double norm2_ac = inner_prod(line_ac, line_ac);
        double norm2_ab = inner_prod(line_ab, line_ab);
        double norm2_cross = inner_prod(ab_cross_ac, ab_cross_ac);

        c_vector<double, 3> vector_a_circumcenter = (norm2_ac * VectorProduct(ab_cross_ac, line_ab)
                                                     - norm2_ab * VectorProduct(ab_cross_ac, line_ac))
            / (2.0 * norm2_cross);

        mVoronoiVertices.push_back(vertex_a + vector_a_circumcenter);
        mMapFaceToVertices[(*it)[0]].push_back(index);
        mMapFaceToVertices[(*it)[1]].push_back(index);
        mMapFaceToVertices[(*it)[2]].push_back(index);
        index++;
    }
}

void VoronoiSphereGenerator::OrderVoronoiVertexMap()
{
    /**
     * We have a map from face centers (i.e. delaunay vertices) to the voronoi
     * vertex indices. To order them (anti/)clockwise, we start with one random
     * index, look for a random neighbouring Delaunay triangle (by checking shared
     * edges) and then continue this way until we come back.
     */
    // Loop through all (voronoi) faces
    for (typename std::map<unsigned, std::vector<unsigned> >::iterator it = this->mMapFaceToVertices.begin();
         it != this->mMapFaceToVertices.end();
         ++it)
    {
        unsigned face_index = it->first;
        std::vector<unsigned>& voronoi_indices = it->second;
        // We leave the first voronoi index unchanged in position
        std::array<unsigned, 3> delaunay_indices = mDelaunayTriangleIndices[voronoi_indices[0]];
        unsigned current_delaunay_index = (delaunay_indices[0] == face_index) ? delaunay_indices[1] : delaunay_indices[0];
        unsigned num_voronoi_vertices = voronoi_indices.size();
        // Now search for the next voronoi index in this face
        for (unsigned index_voronoi = 1; index_voronoi < num_voronoi_vertices; index_voronoi++)
        {
            // We check all vertices with higher index (not ordered yet)
            for (unsigned larger_index_voronoi = index_voronoi; larger_index_voronoi < num_voronoi_vertices; larger_index_voronoi++)
            {
                delaunay_indices = mDelaunayTriangleIndices[voronoi_indices[larger_index_voronoi]];
                // If we find the current Delaunay index (i.e. the triangle vertex of neighbour) ...
                std::array<unsigned, 3>::iterator found_ind_current = std::find(delaunay_indices.begin(), delaunay_indices.end(), current_delaunay_index);
                if (found_ind_current != delaunay_indices.end())
                {
                    // ... we reassign the current_delaunay_index as the new neighbour
                    std::array<unsigned, 3>::iterator found_ind_face = std::find(delaunay_indices.begin(), delaunay_indices.end(), face_index);

                    unsigned local_face_index = std::distance(delaunay_indices.begin(), found_ind_face);
                    unsigned local_current_index = std::distance(delaunay_indices.begin(), found_ind_current);
                    unsigned local_next_index;
                    if (local_face_index == 0)
                        local_next_index = (local_current_index == 1) ? 2 : 1;
                    else if (local_face_index == 1)
                        local_next_index = (local_current_index == 0) ? 2 : 0;
                    else
                        local_next_index = (local_current_index == 0) ? 1 : 0;

                    current_delaunay_index = delaunay_indices[local_next_index];

                    // ... and we switch the voronoi index to the correct position
                    unsigned temp_index_voronoi = voronoi_indices[index_voronoi];
                    voronoi_indices[index_voronoi] = voronoi_indices[larger_index_voronoi];
                    voronoi_indices[larger_index_voronoi] = temp_index_voronoi;
                    break;
                }
            }
        }
    }
}

std::vector<c_vector<double, 3> > VoronoiSphereGenerator::GetVoronoiVertices()
{
    return mVoronoiVertices;
}

std::map<unsigned, std::vector<unsigned> > VoronoiSphereGenerator::GetMapFromFacesToVertices()
{
    return mMapFaceToVertices;
}

VoronoiSphereGenerator::~VoronoiSphereGenerator()
{
}
