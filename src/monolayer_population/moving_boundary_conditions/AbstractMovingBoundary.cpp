#include "AbstractMovingBoundary.hpp"

template <unsigned DIM>
AbstractMovingBoundary<DIM>::AbstractMovingBoundary(
    double force_constant, boost::shared_ptr<AbstractMoveProfile<DIM> > move_profile, bool above)
        : mForceConstant(force_constant), mMoveProfile(move_profile),
          mAbove(above), mForceAppliedOnPopulation{ zero_vector<double>(DIM) } {}

template <unsigned DIM>
AbstractMovingBoundary<DIM>::~AbstractMovingBoundary() {}

template <unsigned DIM>
void AbstractMovingBoundary<DIM>::AddAppliedForceOnPopulation(
    c_vector<double, DIM> force)
{
    mForceAppliedOnPopulation += force;
}

template <unsigned DIM>
void AbstractMovingBoundary<DIM>::ResetAppliedForceOnPopulation()
{
    mForceAppliedOnPopulation = zero_vector<double>(DIM);
}

template <unsigned DIM>
c_vector<double, DIM>
AbstractMovingBoundary<DIM>::GetAppliedForceOnPopulation()
{
    return mForceAppliedOnPopulation;
}

template <unsigned DIM>
double AbstractMovingBoundary<DIM>::GetForceConstant()
{
    return mForceConstant;
}

template <unsigned DIM>
double AbstractMovingBoundary<DIM>::GetZPosition()
{
    return mMoveProfile->GetZPosition();
}

template <unsigned DIM>
void AbstractMovingBoundary<DIM>::SetMoveProfile(
    boost::shared_ptr<AbstractMoveProfile<DIM> > move_profile)
{
    mMoveProfile = move_profile;
}

template <unsigned DIM>
void AbstractMovingBoundary<DIM>::SetForceConstant(double k)
{
    mForceConstant = k;
}

// Explicit instantiation
template class AbstractMovingBoundary<1>;
template class AbstractMovingBoundary<2>;
template class AbstractMovingBoundary<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
