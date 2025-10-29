#ifndef LINEARMOVEPROFILE_HPP_
#define LINEARMOVEPROFILE_HPP_

#include "AbstractMoveProfile.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

template <unsigned DIM>
class LinearMoveProfile : public AbstractMoveProfile<DIM>
{
private:
    /**
     * Slope of linear movement
     */
    double mSlope;

    unsigned mNumberOfSteps;

    double mTimeStep;

    double mStepRatio;

public:
    /**
     * Constructor for simple plate (can either be on top or bottom of the sphere)
     */

    LinearMoveProfile(double start_time, double stop_time,
                      double start_position, double final_position,
                      unsigned number_of_steps, double step_ratio);

    /**
     * Get the z position at the current simulation time
     */
    virtual double GetZPosition();
};

#endif // LINEARMOVEPROFILE_HPP_