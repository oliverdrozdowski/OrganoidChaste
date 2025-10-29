#include "AbstractMoveProfile.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"

template <unsigned DIM>
AbstractMoveProfile<DIM>::AbstractMoveProfile(double start_time, double stop_time, double start_position, double final_position) : mStartTime(start_time), mStopTime(stop_time), mStartPosition(start_position), mFinalPosition(final_position) {}

template <unsigned DIM>
AbstractMoveProfile<DIM>::AbstractMoveProfile() {}

template <unsigned DIM>
AbstractMoveProfile<DIM>::~AbstractMoveProfile() {}

// Explicit instantiation
template class AbstractMoveProfile<1>;
template class AbstractMoveProfile<2>;
template class AbstractMoveProfile<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
