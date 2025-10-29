#include "OneStepProfile.hpp"
#include "SimulationTime.hpp"

template <unsigned DIM>
OneStepProfile<DIM>::OneStepProfile(double start_time, double stop_time, double start_position, double final_position)
        : AbstractMoveProfile<DIM>(start_time, stop_time, start_position, final_position), mStatic(false), mStaticPosition()
{
    assert(this->mStartTime < this->mStopTime);
}

template <unsigned DIM>
double OneStepProfile<DIM>::GetZPosition()
{
    if (mStatic)
    {
        return mStaticPosition;
    }
    else
    {
        // Get current simulation time
        double current_time = SimulationTime::Instance()->GetTime();

        assert(current_time >= this->mStartTime);
        assert(current_time <= this->mStopTime);

        return (this->mFinalPosition - this->mStartPosition) / (this->mStopTime - this->mStartTime) * (current_time - this->mStartTime) + this->mStartPosition;
    }
}

template <unsigned DIM>
void OneStepProfile<DIM>::UpdateParameters(double duration, double distance)
{

    this->mStartTime = this->mStopTime;
    this->mStopTime = this->mStartTime + duration;
    this->mStartPosition = this->mFinalPosition;
    this->mFinalPosition = this->mStartPosition + distance;
}

template <unsigned DIM>
void OneStepProfile<DIM>::UpdateParameters(double duration)
{

    this->mStartTime = this->mStopTime;
    this->mStopTime = this->mStartTime + duration;
    this->mStartPosition = this->mFinalPosition;
}

template <unsigned DIM>
void OneStepProfile<DIM>::Stop()
{
    mStaticPosition = this->GetZPosition();
    mStatic = true;
}

// Explicit instantiation
template class OneStepProfile<1>;
template class OneStepProfile<2>;
template class OneStepProfile<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
