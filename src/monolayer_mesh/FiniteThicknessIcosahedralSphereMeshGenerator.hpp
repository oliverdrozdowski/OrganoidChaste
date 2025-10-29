#ifndef FINITETHICKNESSICOSAHEDRALSPHEREMESHGENERATOR_HPP_
#define FINITETHICKNESSICOSAHEDRALSPHEREMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "HoneycombVertexMeshGenerator.hpp"
#include "MutableMonolayerVertexMesh.hpp"

/**
 * Honeycomb mesh generator that creates a 3d finite thickness honeycomb mesh (with equal distance
 * between nodes and different thickness) for use in vertex simulations.
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */

class FiniteThicknessIcosahedralSphereMeshGenerator
{
    friend class TestFiniteThicknessIcosahedralSphereMeshGenerator;

private:
    unsigned mNumElements;
    double mInnerRadius;
    double mHeight;

    /**
     * Calculate the normal vector of the face.
     *
     * @param pFace pointer to face
     *
     * @return normal vector
     */
    c_vector<double, 3> GetNormalToFace(MonolayerVertexElement<2, 3>* pFace);

    /**
     * Calculate the center vector of the face.
     *
     * @param pFace pointer to face
     *
     * @return center vector
     */
    c_vector<double, 3> GetCenterOfFace(MonolayerVertexElement<2, 3>* pFace);

    /**
     * Map the nodes, which are assumed to lie on a sphere, onto the icosahedron
     * which circumscribes the sphere. The icosahedron is described by the normal
     * vectors which are normal to each face and point outward.
     *
     * @param vectorOfNodes vector containing pointers to the nodes to remap
     * @param faceNormals vector of normals for icosahedron faces
     */
    void ProjectOntoCircumscribingIcosahedronAndRescale(std::vector<Node<3>*> vectorOfNodes,
                                                        std::vector<c_vector<double, 3> > faceNormals,
                                                        double prescribedOuterRadius);

protected:
    /** A pointer to the mesh this class creates */
    MutableMonolayerVertexMesh<3, 3>* mpMesh;

public:
    /**
     * Constructor.
     *
     * @param casparKlugNumH  H in the Caspar-Klug characterization of hexagonal tilings (H,K)
     * @param casparKlugNumH  K in the Caspar-Klug characterization of hexagonal tilings (H,K) (defaults to 0)
     * @param numElementsUp  The number of rows of elements in the mesh
     * @param isFlatBottom  Whether to enforce a flat bottom to the mesh (defaults to false)
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangement (defaults to 0.01)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     * @param height the element height, which has default value (defaults to 1.0)
     * @param elementArea the element area, which has default value 0.5*sqrt(3.0)
     * @param whether the non-tip-points should be projected onto the sphere, defaults to true
     */
    FiniteThicknessIcosahedralSphereMeshGenerator(unsigned casparKlugNumH,
                                                  unsigned casparKlugNumK = 0,
                                                  double cellRearrangementThreshold = 0.01,
                                                  double t2Threshold = 0.001,
                                                  double height = 1.0,
                                                  double innerRadius = 10.0,
                                                  bool makeSphere = true);

    /**
     * Null constructor for derived classes to call.
     */
    FiniteThicknessIcosahedralSphereMeshGenerator()
            : mNumElements(0),
              mInnerRadius(0.0),
              mHeight(0.0)
    {
    }

    /**
     * Destructor - deletes the mesh object and pointer.
     */
    virtual ~FiniteThicknessIcosahedralSphereMeshGenerator();

    /**
     * @return a 3D spherical mesh
     */
    virtual MutableMonolayerVertexMesh<3, 3>* GetMesh();
};

#endif /*FINITETHICKNESSICOSAHEDRALSPHEREMESHGENERATOR_HPP_*/