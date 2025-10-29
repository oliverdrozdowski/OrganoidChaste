#ifndef MECHANICALMOSAICCELLMUTATIONSTATE_HPP_
#define MECHANICALMOSAICCELLMUTATIONSTATE_HPP_

#include <boost/serialization/base_object.hpp>
#include "AbstractCellMutationState.hpp"
#include "ChasteSerialization.hpp"

/**
 * Subclass of AbstractCellMutationState defining a 'wild type' mutation state.
 */
class MechanicalMosaicCellMutationState : public AbstractCellMutationState
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell mutation state.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellMutationState>(*this);
    }

public:
    /**
     * Constructor.
     */
    MechanicalMosaicCellMutationState();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(MechanicalMosaicCellMutationState)

#endif /* MECHANICALMOSAICCELLMUTATIONSTATE_HPP_ */