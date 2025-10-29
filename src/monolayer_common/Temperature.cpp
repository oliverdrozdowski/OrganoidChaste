#include "Temperature.hpp"
#include <iostream>
#include "SimulationTime.hpp"

Temperature::Temperature() : mStartTime(0.0), mDecayTime(0.0), mIsConstant(true)
{
}

Temperature::Temperature(double start_time, double decay_time) : mStartTime(start_time), mDecayTime(decay_time), mIsConstant(false)
{
}

double Temperature::GetCurrentTemperature()
{
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();

    if (mIsConstant || current_time < mStartTime)
    {
        return 1.0;
    }
    else
    {
        if ((current_time - mStartTime) < mDecayTime)
        {
            // Let the temperature decay linearly
            return 1.0 - (current_time - mStartTime) / mDecayTime;
        }
        else
        {
            // Negative DecayTime allows for T-growth
            if (mDecayTime < 0.0)
                return -(current_time - mStartTime) / mDecayTime;
        }
    }

    return 0.0;
}

void Temperature::SetStartTimeToNow()
{
    SimulationTime* p_simulation_time = SimulationTime::Instance();

    mStartTime = p_simulation_time->GetTime();
}

void Temperature::SetStartTime(double start_time)
{
    mStartTime = start_time;
}

void Temperature::SetAbsoluteDecayTime(double absolute_decay_time)
{
    mIsConstant = false;
    mDecayTime = absolute_decay_time - mStartTime;
}

void Temperature::SetRelativeDecayTime(double relative_decay_time)
{
    mIsConstant = false;
    mDecayTime = relative_decay_time;
}

void Temperature::SetTemperatureToConstant()
{
    mIsConstant = true;
}
