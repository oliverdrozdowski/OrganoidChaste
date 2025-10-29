#include "LinearMoveProfile.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
LinearMoveProfile<DIM>::LinearMoveProfile(double start_time, double stop_time,
                                          double start_position,
                                          double final_position,
                                          unsigned number_of_steps,
                                          double step_ratio)
        : AbstractMoveProfile<DIM>(start_time, stop_time, start_position,
                                   final_position),
          mSlope(), mNumberOfSteps(number_of_steps), mTimeStep(),
          mStepRatio(step_ratio)
{

    assert(this->mStartTime < this->mStopTime);

    // calculate slope
    double distance = this->mFinalPosition - this->mStartPosition;
    double duration = (this->mStopTime - this->mStartTime);
    double move_duration = duration / mStepRatio;

    mTimeStep = duration / (mNumberOfSteps * mStepRatio);

    mSlope = distance / move_duration;
}

template <unsigned DIM>
double LinearMoveProfile<DIM>::GetZPosition()
{
    // Get current simulation time
    double current_time = SimulationTime::Instance()->GetTime();

    assert(current_time >= this->mStartTime);
    assert(current_time <= this->mStopTime);

    double relative_time = current_time - this->mStartTime;

    for (unsigned n = 0; n < (mNumberOfSteps * mStepRatio); n++)
    {
        if (relative_time <= (n + 1) * mTimeStep)
        {
            if ((n % static_cast<int>(mStepRatio)) == 0)
            {
                return mSlope * (relative_time - n * mTimeStep) + (mSlope * n / mStepRatio * mTimeStep) + this->mStartPosition;
            }
            else
            {
                return (1 + floor(n / mStepRatio)) * mSlope * mTimeStep + this->mStartPosition;
            }
        }
    }
    return 0.0;
}

// Explicit instantiation
template class LinearMoveProfile<1>;
template class LinearMoveProfile<2>;
template class LinearMoveProfile<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
