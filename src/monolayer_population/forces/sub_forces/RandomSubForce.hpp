#ifndef RANDOMSUBFORCE_HPP_
#define RANDOMSUBFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"
#include "Exception.hpp"

#include "AbstractCellBasedSimulation.hpp"
#include "AbstractSubForce.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "Temperature.hpp"

#include <iostream>

/**
 * A force class for use in MonolayerVertex-based simulations. This force is based on the
 * Energy function containing cell-specific surface tensions. Each face can have its own tension.
 * Therefore the correpsonding simulation modifier is needed if we have topological changes.
 */

template <unsigned DIM>
class RandomSubForce : public AbstractSubForce<DIM>
{
    friend class TestForces;
    friend class TestMonolayerForces;

private:
    /**
     * Strength of the random force
     */
    double mStrength;

    /**
     * Temperature of the random force
     */
    boost::shared_ptr<Temperature> mTemperature{};

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
        archive & mStrength;
        archive & mTemperature;
    }

public:
    /**
     * Constructor with given Temperature
     */
    RandomSubForce(double strength, boost::shared_ptr<Temperature> temp);

    /**
     * Destructor.
     */
    virtual ~RandomSubForce();

    /**
     * Overridden GetForceContributions() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on given boundary
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
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

    /**
     * Reset the global temperature
     */
    void ResetTemperature();

    /**
     * Set the random strength
     */
    void SetStrength(double strength);

    /**
     * Set the temperature decay time (relative to start time of temperature)
     */
    void SetTemperatureDecay(double decay_time);

    /**
     * Set temperature to constant
     */
    void SetTemperatureToConstant();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomSubForce)

#endif /*RANDOMSUBFORCE_HPP_*/