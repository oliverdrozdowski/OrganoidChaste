#ifndef FACETYPEWRITER_HPP_
#define FACETYPEWRITER_HPP_

#include <boost/serialization/base_object.hpp>
#include "AbstractFaceWriter.hpp"
#include "ChasteSerialization.hpp"

/**
 * A class written using the visitor pattern for writing cell volumes (in 3D, or areas in 2D) to file.
 *
 * The output file is called cellareas.dat by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored in
 * the VTK cell data "Cell areas" by default.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class FaceTypeWriter : public AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>
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
        archive& boost::serialization::base_object<AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:
    /**
     * Default constructor.
     */
    FaceTypeWriter();

    /**
     * Overridden GetCellDataForVtkOutput() method.
     *
     * Get a double associated with a face. This method reduces duplication
     * of code between the methods VisitCell() and AddVtkData().
     *
     * @param pFace a face
     * @param pCellPopulation a pointer to the cell population owning the cell
     *
     * @return data associated with the cell
     */
    double GetFaceDataForVtkOutput(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Overridden VisitCell() method.
     *
     * Visit a face and write its face type.
     *
     * Outputs a line of space-separated values of the form:
     * ...[location index] [cell id] [x-pos] [y-pos] [z-pos] [face type] ...
     * with [y-pos] and [z-pos] included for 2 and 3 dimensional simulations, respectively.
     *
     * This is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pFace a face
     * @param pCellPopulation a pointer to the cell population owning the cell
     */
    virtual void VisitFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS2(FaceTypeWriter, 2, 2)
EXPORT_TEMPLATE_CLASS2(FaceTypeWriter, 2, 3)
EXPORT_TEMPLATE_CLASS2(FaceTypeWriter, 3, 3)

#endif /* FACETYPEWRITER_HPP_ */