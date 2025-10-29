#ifndef STEPFUNCTIONMOVEPROFILE_HPP_
#define STEPFUNCTIONMOVEPROFILE_HPP_

#include "AbstractMoveProfile.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

template <unsigned DIM>
class StepFunctionMoveProfile : public AbstractMoveProfile<DIM>
{
private:
    const unsigned mNumberOfSteps;

    double mTimeStep;

    double mZStep;

public:
    /**
     * Constructor for simple plate (can either be on top or bottom of the sphere)
     */
    StepFunctionMoveProfile(const double start_time, const double stop_time,
                            const double start_position,
                            const double final_position,
                            const unsigned number_of_steps);

    /**
     * Get the z position at the current simulation time
     */
    double GetZPosition();
};

#endif // STEPFUNCTIONMOVEPROFILE_HPP_