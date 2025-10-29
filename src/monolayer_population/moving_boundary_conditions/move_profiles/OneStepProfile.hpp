#ifndef ONESTEPPROFILE_HPP_
#define ONESTEPPROFILE_HPP_

#include "AbstractMoveProfile.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

template <unsigned DIM>
class OneStepProfile : public AbstractMoveProfile<DIM>
{
private:
    /**
     * Bool on if to keep the move profile constant at the current position
     */
    bool mStatic;

    /**
     * Static position of profile
     */
    double mStaticPosition;

public:
    /**
     * Constructor for simple plate (can either be on top or bottom of the sphere)
     */

    OneStepProfile(double start_time, double stop_time, double start_position, double final_position);

    /**
     * Get the z position at the current simulation time
     */
    virtual double GetZPosition();

    /**
     * Update parameters
     */
    void UpdateParameters(double duration, double distance);

    /**
     * Update parameters if boundary should be stationary
     */
    void UpdateParameters(double duration);

    /**
     * Stop the profile and keep it constant at the current position
     */
    void Stop();
};

#endif // ONESTEPPROFILE_HPP_