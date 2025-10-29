#ifndef RECTANGULAREDGESUBFORCE_HPP_
#define RECTANGULAREDGESUBFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"
#include "Exception.hpp"

#include "AbstractCellBasedSimulation.hpp"
#include "AbstractSubForce.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"

#include <iostream>

/**
 * A force class for use in MonolayerVertex-based simulations. This force applies
 * a constant force onto nodes belonging to a certain edge. It is meant to be used
 * in rectangular grids, where left, top, right and bottom edges can be defined
 * but it is universally usable if constant forces are supposed to be supplied to
 * a list of known nodes.
 */

template <unsigned DIM>
class RectangularEdgeSubForce : public AbstractSubForce<DIM>
{
private:
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractSubForce<DIM> >(*this);
        archive & mLeftEdgeNodes;
        archive & mTopEdgeNodes;
        archive & mRightEdgeNodes;
        archive & mBottomEdgeNodes;
        archive & mLeftEdgeForce;
        archive & mTopEdgeForce;
        archive & mRightEdgeForce;
        archive & mBottomEdgeForce;
        archive & mUseInternalCoordinates;
        archive & mParallelizeOppositeEdgeNormals;
        archive & mpSimulation;
    }

    double CalculateTotalEnergy(AbstractCellPopulation<DIM>& rCellPopulation);

protected:
    /**
     * Vector containing pointers to the nodes belonging to the left edge.
     */
    std::vector<Node<DIM>*> mLeftEdgeNodes;

    /**
     * Vector containing pointers to the nodes belonging to the upper edge.
     */
    std::vector<Node<DIM>*> mTopEdgeNodes;

    /**
     * Vector containing pointers to the nodes belonging to the right edge.
     */
    std::vector<Node<DIM>*> mRightEdgeNodes;

    /**
     * Vector containing pointers to the nodes belonging to the bottom edge.
     */
    std::vector<Node<DIM>*> mBottomEdgeNodes;

    c_vector<double, DIM> mLeftEdgeForce;
    c_vector<double, DIM> mTopEdgeForce;
    c_vector<double, DIM> mRightEdgeForce;
    c_vector<double, DIM> mBottomEdgeForce;

    bool mUseInternalCoordinates = false;

    bool mParallelizeOppositeEdgeNormals = false;

    AbstractCellBasedSimulation<DIM, DIM>* mpSimulation;

public:
    /**
     * Constructor.
     */
    RectangularEdgeSubForce();

    /**
     * Destructor.
     */
    virtual ~RectangularEdgeSubForce();

    /**
     * Overridden GetForceContributions() method.
     *
     * Calculate the force on each node in the vertex-based cell population according to the edges
     *
     *
     * @param rCellPopulation reference to the cell population
     * @param pSimulation pointer to simulation instance
     *
     * @return the map of nodes to forces
     */
    NodeSubForceMap<DIM> GetForceContributions(MonolayerVertexBasedCellPopulation<DIM>& rCellPopulation,
                                               AbstractCellBasedSimulation<DIM, DIM>* pSimulation);

    /**
     * @return mLeftEdgeNodes
     */
    std::vector<Node<DIM>*> GetLeftEdgeNodes();

    /**
     * @return mTopEdgeNodes
     */
    std::vector<Node<DIM>*> GetTopEdgeNodes();

    /**
     * @return mRightEdgeNodes
     */
    std::vector<Node<DIM>*> GetRightEdgeNodes();

    /**
     * @return mBottomEdgeNodes
     */
    std::vector<Node<DIM>*> GetBottomEdgeNodes();

    /**
     * @param leftEdgeNodes
     */
    void SetLeftEdgeNodes(std::vector<Node<DIM>*> leftEdgeNodes);

    /**
     * @param topEdgeNodes
     */
    void SetTopEdgeNodes(std::vector<Node<DIM>*> topEdgeNodes);

    /**
     * @param rightEdgeNodes
     */
    void SetRightEdgeNodes(std::vector<Node<DIM>*> rightEdgeNodes);

    /**
     * @param bottomEdgeNodes
     */
    void SetBottomEdgeNodes(std::vector<Node<DIM>*> bottomEdgeNodes);

    /**
     * @return mLeftEdgeForce
     */
    c_vector<double, DIM> GetLeftEdgeForce();

    /**
     * @return mTopEdgeForce
     */
    c_vector<double, DIM> GetTopEdgeForce();

    /**
     * @return mRightEdgeForce
     */
    c_vector<double, DIM> GetRightEdgeForce();

    /**
     * @return mBottomEdgeForce
     */
    c_vector<double, DIM> GetBottomEdgeForce();

    /**
     * @param leftEdgeForce
     */
    void SetLeftEdgeForce(c_vector<double, DIM> leftEdgeForce);

    /**
     * @param topEdgeForce
     */
    void SetTopEdgeForce(c_vector<double, DIM> topEdgeForce);

    /**
     * @param rightEdgeForce
     */
    void SetRightEdgeForce(c_vector<double, DIM> rightEdgeForce);

    /**
     * @param bottomEdgeForce
     */
    void SetBottomEdgeForce(c_vector<double, DIM> bottomEdgeForce);

    /**
     * Whether to use internal coordinates. Then the first component of the force
     * vectors coorresponds to the edge normal (average vector to population center),
     * the second component to the average normal of the apical faces and the third
     * component to the cross-product, which should be tangential to the edge.
     *
     * @param useInternalCoordinates
     */
    void UseInternalCoordinates(bool useInternalCoordinates = true);

    /**
     * Whether to parallelize opposite edge normals, i.e., use the same coordinate
     * system on opposite edges. This makes sense if the symmetry is broken by hand
     * as the coordinate systems may become incompatible, which will not correct itself.
     * Only makes sense if internal coordinates are used.
     *
     * @param parallelizeOppositeEdgeNormals
     */
    void ParallelizeOppositeEdgeNormals(bool parallelizeOppositeEdgeNormals = true);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RectangularEdgeSubForce)

#endif /*RECTANGULAREDGESUBFORCE_HPP_*/