#include "StaticMoveProfile.hpp"

template <unsigned DIM>
StaticMoveProfile<DIM>::StaticMoveProfile(double position)
        : AbstractMoveProfile<DIM>()
{
    this->mStartPosition = position;
}

template <unsigned DIM>
double StaticMoveProfile<DIM>::GetZPosition()
{
    return this->mStartPosition;
}

template <unsigned DIM>
void StaticMoveProfile<DIM>::UpdatePosition(double position)
{
    this->mStartPosition = position;
}

// Explicit instantiation
template class StaticMoveProfile<1>;
template class StaticMoveProfile<2>;
template class StaticMoveProfile<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
