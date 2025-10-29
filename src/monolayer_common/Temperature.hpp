#ifndef TEMPERATURE_HPP_
#define TEMPERATURE_HPP_

#include <boost/serialization/shared_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include "ChasteSerialization.hpp"
#include "SerializableSingleton.hpp"

/**
 * Simulation temperature object stores the temperature.
 * It uses the singleton pattern to provide a globally consistent temperature.
 *
 * Using a decay time and a start time allows for a simulated annealing with
 * a linearly decaying temperature from 1.0 to 1.
 * Before the start time the temperature is always 1.0
 * If decay time is negative, then the temperature increases linearly from 1 after start time
 * by a factor (t-start_time)/decay_time.
 */
class Temperature
{
public:
    /**
     * Temperature constructor for a constant temperature
     */
    Temperature();

    /**
     * Temperature constructor for a given start time and decay time.
     */
    Temperature(double start_time, double decay_time);

    /*
     * Get the current temperature
     */
    double GetCurrentTemperature();

    /**
     * Reset temperature to start at current timestamp
     */
    void SetStartTimeToNow();

    /**
     * Set start time manually
     */
    void SetStartTime(double start_time);

    /**
     * Set absolute decay time
     */
    void SetAbsoluteDecayTime(double absolute_decay_time);

    /**
     * Set relative decay time
     */
    void SetRelativeDecayTime(double relative_decay_time);

    /**
     * Set temperature to constant
     */
    void SetTemperatureToConstant();

private:
    /**
     * Stores the starting time
     */
    double mStartTime;

    /**
     * Stores the time period over which the temperature decays to zero
     */
    double mDecayTime;

    /**
     * Bool on if temperature is constant
     */
    double mIsConstant;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialization of a Temperature object must be done with care.
     * Do not serialize this singleton directly.  Instead, serialize
     * the object returned by GetSerializationWrapper.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & mDecayTime;
        archive & mStartTime;
        archive & mIsConstant;
    }
};

#endif /*TEMPERATURE_HPP_*/
