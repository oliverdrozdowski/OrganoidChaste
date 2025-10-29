#ifndef FINITETHICKNESSHONEYCOMBVERTEXMESHGENERATOR_HPP_
#define FINITETHICKNESSHONEYCOMBVERTEXMESHGENERATOR_HPP_

#include <array>
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

class FiniteThicknessHoneycombVertexMeshGenerator
{
private:
    unsigned mNumElementsAcross;
    unsigned mNumElementsUp;
    double mElementArea;
    double mHeight;

    std::vector<Node<3>*> mNodesLeftEdge;
    std::vector<Node<3>*> mNodesRightEdge;
    std::vector<Node<3>*> mNodesTopEdge;
    std::vector<Node<3>*> mNodesBottomEdge;

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
     * @param elementArea the element area, which has default value 0.5*sqrt(3.0)
     */
    FiniteThicknessHoneycombVertexMeshGenerator(unsigned numElementsAcross,
                                                unsigned numElementsUp,
                                                bool isFlatBottom = false,
                                                double cellRearrangementThreshold = 0.01,
                                                double t2Threshold = 0.001,
                                                double height = 1.0,
                                                double elementArea = 0.5 * sqrt(3.0));

    /**
     * Null constructor for derived classes to call.
     */
    FiniteThicknessHoneycombVertexMeshGenerator()
            : mNumElementsAcross(0),
              mNumElementsUp(0),
              mElementArea(0.0)
    {
    }

    /**
     * Function to implement a cylindrically curved initial plane configuration. Note that this does not
     * preserve the cell volumes.
     *
     * @param midPoint mid point of cylinder
     * @param midAxis symmetry axis of cylinder
     * @param radius radius of cylidner (to mid-plane)
     */
    void CurveMonolayerCylinder(c_vector<double, 3> midPoint, c_vector<double, 3> midAxis, double radius);

    /**
     * Function to implement a spherically curved initial plane configuration. Note that this does not
     * preserve the cell volumes and it does not preserve the mid-plane areas.
     *
     * @param midPoint mid point of sphere
     * @param radius radius of sphere (to mid-plane)
     */
    void CurveMonolayerSphere(c_vector<double, 3> midPoint, double radius);

    /**
     * Function to implement a curved initial plane configuration. Note that this does not
     * preserve the cell volumes.
     *
     * @param relativeRadius How much larger the radius is than the monolayer dimensions
     */
    void CurveMonolayer(double relativeRadius);

    /**
     * Make the monolayer cylinderical by glueing the left and right edges together
     * Along the azimuth we have mElementsAcross elements and along the axis mElementsUp.
     * Note that this does not preserve the cell volumes.
     * This function should be run before any remeshing was done, as we assume the initial
     * node, face and element numbering and configuration.
     *
     */
    void MakeCylindrical();

    /**
     * Destructor - deletes the mesh object and pointer.
     */
    virtual ~FiniteThicknessHoneycombVertexMeshGenerator();

    /**
     * @return a 3D honeycomb mesh
     */
    virtual MutableMonolayerVertexMesh<3, 3>* GetMesh();

    /**
     * @return an array containing the left, top, right and bottom edges, given as vectors of edge nodes
     */
    std::array<std::vector<Node<3>*>, 4> GetEdgesWithNodes();
};

#endif /*FINITETHICKNESSHONEYCOMBVERTEXMESHGENERATOR_HPP_*/