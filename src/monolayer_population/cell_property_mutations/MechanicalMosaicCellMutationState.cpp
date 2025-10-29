#include "MechanicalMosaicCellMutationState.hpp"

MechanicalMosaicCellMutationState::MechanicalMosaicCellMutationState()
        : AbstractCellMutationState(1)
{
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(MechanicalMosaicCellMutationState)
