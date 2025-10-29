#include "StepFunctionMoveProfile.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
StepFunctionMoveProfile<DIM>::StepFunctionMoveProfile(
    const double start_time, const double stop_time,
    const double start_position, const double final_position,
    const unsigned number_of_steps)
        : AbstractMoveProfile<DIM>(start_time, stop_time, start_position,
                                   final_position),
          mNumberOfSteps(number_of_steps)
{

    assert(this->mStartTime < this->mStopTime);

    double duration = this->mStopTime - this->mStartTime;

    mTimeStep = duration / mNumberOfSteps;

    double distance = this->mFinalPosition - this->mStartPosition;

    mZStep = distance / mNumberOfSteps;
}

template <unsigned DIM>
double StepFunctionMoveProfile<DIM>::GetZPosition()
{

    double current_time = SimulationTime::Instance()->GetTime();
    assert(current_time >= this->mStartTime);
    assert(current_time <= this->mStopTime);

    for (unsigned n = 0; n < mNumberOfSteps; n++)
    {
        if ((((n * mTimeStep) + this->mStartTime) <= current_time) & (current_time < (((n + 1) * mTimeStep) + this->mStartTime)))
        {
            return (n * mZStep) + this->mStartPosition;
        }
    }

    return DOUBLE_UNSET;
}

// Explicit instantiation
template class StepFunctionMoveProfile<1>;
template class StepFunctionMoveProfile<2>;
template class StepFunctionMoveProfile<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
