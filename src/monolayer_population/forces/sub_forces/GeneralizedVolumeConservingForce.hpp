#ifndef GENERALIZEDVOLUMECONSERVINGFORCE_HPP_
#define GENERALIZEDVOLUMECONSERVINGFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"
#include "Exception.hpp"

#include "AbstractCellBasedSimulation.hpp"
#include "AbstractForce.hpp"
#include "AbstractSubForce.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "SmallAreaDampingModifier.hpp"
#include "SurfaceTensionSubForce.hpp"

#include <iostream>

// Forward declaration prevents circular include chain
struct DampingCoefficientPair;

/**
 * A force class for use in MonolayerVertex-based simulations. This force contains different
 * AbstractSubForce objects, which it uses to calculate a total force in each node
 * Then the force is projected onto the manifold of volume conserving forces.
 * Also a correction step is performed if the volume is not the target volume.
 */

template <unsigned DIM>
class GeneralizedVolumeConservingForce : public AbstractForce<DIM>
{
    friend class TestForces;
    friend class TestMonolayerForces;

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
        archive& boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mSubForceCollection;
        archive & mProjectOnVolumeConservingForces;
        archive & mIndividualNodeDampingActivated;
        archive & mMapNodeToDampingThresholds;
        archive & mpSimulation;
    }

protected:
    /**
     * A list of all sub forces which are considered on the nodes and projected
     * onto the volume conserving manifold.
     */
    std::vector<boost::shared_ptr<AbstractSubForce<DIM> > > mSubForceCollection;

    bool mProjectOnVolumeConservingForces = true;

    /**
     * Whether the different nodes have different mobilities. Is set by methods writing the corresponding map
     */
    bool mIndividualNodeDampingActivated;

    /**
     * Map of node pointers to damping thresholds which should be applied to the nodes
     */
    std::map<Node<DIM>*, DampingCoefficientPair> mMapNodeToDampingThresholds;

    /**
     * Pointer to the simulation object in case a subforce needs it.
     */
    AbstractCellBasedSimulation<DIM, DIM>* mpSimulation;

    /**
     * Flag if constant lumen volume is wanted. Default is false
     */
    bool mConstantLumenVolume;

public:
    /**
     * Constructor.
     */
    GeneralizedVolumeConservingForce();

    /**
     * Destructor.
     */
    virtual ~GeneralizedVolumeConservingForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on the provided subforces
     * The resulting total force is then projected onto the manifold of volume conserving forces.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Do we project on volume conserving forces?
     *
     * @return mProjectOnVolumeConservingForces
     */
    bool DoProjectOnVolumeConservingForces();

    /**
     * Set if we project on volume conserving forces.
     *
     * @param doProjection whether to do the projection
     */
    void SetProjectOnVolumeConservingForces(bool doProjection);

    /**
     * Set pointer to Simulation instance
     *
     * @param pSimulation pointer to simulation instance
     */
    void SetSimulationInstance(AbstractCellBasedSimulation<DIM, DIM>* pSimulation);

    /**
     * Get pointer to Simulation instance
     *
     * @return mpSimulation pointer to simulation instance
     */
    AbstractCellBasedSimulation<DIM, DIM>* GetSimulationInstance();

    /**
     * @param MapNodeToDampingThresholds
     */
    void SetMapNodeToDampingThresholds(std::map<Node<DIM>*, DampingCoefficientPair> mapNodeToDampingThresholds);

    /**
     * Add a subforce to consider.
     *
     * @param pForce boost pointer to force
     */
    void AddSubForce(boost::shared_ptr<AbstractSubForce<DIM> > pForce);

    /**
     * Remove a subforce from the consideration.
     *
     * @param pForce boost pointer to force
     */
    void RemoveSubForce(boost::shared_ptr<AbstractSubForce<DIM> > pForce);

    /**
     * Get pointer to SurfaceTensionSubforce if it is in mpSubForceCollection, otherwise
     * return nullptr.
     *
     * @return boost pointer to subforce
     */
    boost::shared_ptr<SurfaceTensionSubForce<DIM> > GetSurfaceTensionSubForce();

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

    /**
     * Set flag for constant lumen volume
     */
    void SetConstantLumenVolumeFlag(bool keepLumenVolumeConstant);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GeneralizedVolumeConservingForce)

#endif /*GENERALIZEDVOLUMECONSERVINGFORCE_HPP_*/