#ifndef ABSTRACTMOVEPROFILE_HPP_
#define ABSTRACTMOVEPROFILE_HPP_

#include "AbstractCellPopulation.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A class that implements a temporal move profile where a plate or similar in the
 * AbstractMovingBoundary is situated. It is implemented via a scalar z-position.
 */
template <unsigned DIM>
class AbstractMoveProfile
{
protected:
    double mStartTime;

    double mStopTime;

    double mStartPosition;

    double mFinalPosition;

public:
    /**
     * Constructor for an abstract move profile. Initialize start time and start
     * position directly
     */
    AbstractMoveProfile(double start_time, double stop_time,
                        double start_position, double final_position);

    /**
     * Empty Constructor for derived classes
     */
    AbstractMoveProfile();

    /**
     * Deconstructor
     */
    virtual ~AbstractMoveProfile();

    /**
     * Get the z position at the current simulation time
     */
    virtual double GetZPosition() = 0;
};

#endif // ABSTRACTMOVEPROFILE_HPP_