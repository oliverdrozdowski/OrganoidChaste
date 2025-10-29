#ifndef POPULATIONSURFACETENSIONENERGYWRITE_HPP_
#define POPULATIONSURFACETENSIONENERGYWRITE_HPP_

#include <boost/serialization/base_object.hpp>
#include "AbstractMonolayerVertexPopulationWriter.hpp"
#include "ChasteSerialization.hpp"
#include "SurfaceTensionSubForce.hpp"

/**
 * A class written using the visitor pattern for writing cell population volume
 * data to file. Used by MeshBasedCellPopulation.
 *
 * The output file is called cellpopulationareas.dat by default.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class PopulationSurfaceTensionEnergyWriter
        : public AbstractMonolayerVertexPopulationWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<
            AbstractMonolayerVertexPopulationWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

    /**
     * shared pointer towards surface tension sub force
     */
    boost::shared_ptr<SurfaceTensionSubForce<SPACE_DIM> > mpForce;

public:
    /**
     * Default constructor.
     */
    PopulationSurfaceTensionEnergyWriter();

    /**
     * Set force pointer.
     *
     * @param pForce
     */
    void SetTensionForcePointer(boost::shared_ptr<SurfaceTensionSubForce<SPACE_DIM> > pForce);

    /**
     *
     * @param pCellPopulation a pointer to the VertexBasedCellPopulation to visit.
     */
    virtual void
    Visit(MonolayerVertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PopulationSurfaceTensionEnergyWriter)

#endif /*POPULATIONSURFACETENSIONENERGYWRITE_HPP_*/