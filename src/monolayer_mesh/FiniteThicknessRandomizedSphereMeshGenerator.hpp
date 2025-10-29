#ifndef FINITETHICKNESSRANDOMIZEDSPHEREMESHGENERATOR_HPP_
#define FINITETHICKNESSRANDOMIZEDSPHEREMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "HoneycombVertexMeshGenerator.hpp"
#include "MutableMonolayerVertexMesh.hpp"

/**
 * Mesh generator that creates a spherical 3d finite thickness mesh for use in vertex simulations.
 * Cell centers are chosen arbitrarily, but random sequential adsorption is also supported,
 * which only allows for cell centers with a minimum distance (cf. Roshal et al., Phys. Rev. E 108,
 * 024404 (2023))
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */

class FiniteThicknessRandomizedSphereMeshGenerator
{
    friend class TestFiniteThicknessRandomizedSphereMeshGenerator;

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

protected:
    /** A pointer to the mesh this class creates */
    MutableMonolayerVertexMesh<3, 3>* mpMesh;

public:
    /**
     * Constructor.
     *
     * @param numElementsAcross  The number of columns of elements in the mesh
     * @param numElementsUp  The number of rows of elements in the mesh
     * @param isFlatBottom  Whether to enforce a flat bottom to the mesh (defaults to false)
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangement (defaults to 0.01)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     * @param height the element height, which has default value (defaults to 1.0)
     * @param innerRadius the inner radius of the sphere
     * @param doRandomSequentialAdsorption whether to use RSA to get more regularly spaced cell centers (defaults to False)
     *
     */
    FiniteThicknessRandomizedSphereMeshGenerator(unsigned numElements,
                                                 double cellRearrangementThreshold = 0.01,
                                                 double t2Threshold = 0.001,
                                                 double height = 1.0,
                                                 double innerRadius = 10.0,
                                                 bool doRandomSequentialAdsorption = false);

    /**
     * Null constructor for derived classes to call.
     */
    FiniteThicknessRandomizedSphereMeshGenerator()
            : mNumElements(0),
              mInnerRadius(0.0),
              mHeight(0.0)
    {
    }

    /**
     * Destructor - deletes the mesh object and pointer.
     */
    virtual ~FiniteThicknessRandomizedSphereMeshGenerator();

    /**
     * @return a 3D spherical mesh
     */
    virtual MutableMonolayerVertexMesh<3, 3>* GetMesh();

    /**
     * Check for overlap between new cell center and all old cell centers. If the distance between centers is smaller
     * than the minimum we return true, otherwise false
     *
     * @param new_cell_center
     * @param r_cell_centers reference to old cell centers
     * @param min_distance
     *
     * @return boolean whether new cell center is too close to old centers
     */
    bool DoesOverlap(c_vector<double, 3> new_cell_center, std::vector<c_vector<double, 3> >& r_cell_centers, double min_distance);

    /**
     * Generate a random 3-vector in the sphere.
     *
     * @return random 3-vector sphere
     */
    c_vector<double, 3> GenerateRandomPointOnSphere();
};

#endif /*FINITETHICKNESSRANDOMIZEDSPHEREMESHGENERATOR_HPP_*/