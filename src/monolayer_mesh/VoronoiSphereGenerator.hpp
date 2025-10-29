#ifndef VORONOISPHEREGENERATOR_HPP_
#define VORONOISPHEREGENERATOR_HPP_

#include <cmath>
#include <map>
#include <vector>
#include "UblasVectorInclude.hpp"

/**
 * Generator that takes a vector of positions that will become the centers of
 * a voronoi diagram on a sphere, where the vertices of the voronoi diagram
 * will be saved as a vector of c_vectors.
 *
 * NOTE: the user should delete the object after use to manage memory.
 */

class VoronoiSphereGenerator
{
private:
    unsigned mNumFaces;
    std::vector<c_vector<double, 3> > mInputFaceCenters;
    std::vector<c_vector<double, 3> > mDelaunayVertices;
    std::vector<std::array<unsigned, 3> > mDelaunayTriangleIndices;
    std::vector<c_vector<double, 3> > mVoronoiVertices;

    /**
     * Map from face indices (identical to index in mDelaunayVertices) to
     * voronoi vertices (identical to index od Delaunay triangle in mDelaunayTriangleIndices)
     */
    std::map<unsigned, std::vector<unsigned> > mMapFaceToVertices;

    /**
     * The Delaunay triangulation on a sphere is computed as the convex hull of the points
     *
     */
    void ConstructDelaunayTriangulation();

    /**
     * The Voronoi vertices are the circumcenters of the Delaunay triangles. They are computed
     * and saved to mVoronoiVertices along with the map from the face index to the voronoi
     * vertex index mMapFaceToVertices.
     *
     */
    void CalculateVoronoiVertices();

    /**
     * The Voronoi vertices are unordered after construction. With this function they are
     * ordered such that they go anti- or clockwise, by iterating through neighbouring
     * Delaunay triangles. The direction is random.
     *
     */
    void OrderVoronoiVertexMap();

public:
    /**
     * Constructor. This constructor uses the given face centers to create the voronoi diagram. We assume
     * all vectors to have the same length (i.e. centered around zero with constant radius)
     *
     */
    VoronoiSphereGenerator(std::vector<c_vector<double, 3> > faceCenters);

    /**
     * Null constructor for derived classes to call.
     */
    VoronoiSphereGenerator()
    {
    }

    /**
     * Destructor.
     */
    virtual ~VoronoiSphereGenerator();

    /**
     * @return the vector of vertex positions
     */
    virtual std::vector<c_vector<double, 3> > GetVoronoiVertices();

    /**
     * @return the map from face index to voronoi vertex indices
     */
    virtual std::map<unsigned, std::vector<unsigned> > GetMapFromFacesToVertices();
};

#endif /*VORONOISPHEREGENERATOR_HPP_*/