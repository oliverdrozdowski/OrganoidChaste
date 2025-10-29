#ifndef FACETENSIONWRITER_HPP_
#define FACETENSIONWRITER_HPP_

#include <boost/serialization/base_object.hpp>
#include "AbstractFaceWriter.hpp"
#include "ChasteSerialization.hpp"
#include "SurfaceTensionSubForce.hpp"

/**
 * A class written using the visitor pattern for writing face tensions to file.
 * It must be used with the SurfaceTensionSubForce, which must be given to the
 * writer before using it.
 *
 * If VTK is switched on,
 * then the writer  specifies the VTK output for each cell, which is stored in
 * the VTK cell data "Face tension" by default.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class FaceTensionWriter : public AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>
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
        archive & mpSurfaceTensionForce;
    }

protected:
    /**
     * Pointer for SurfaceTensionSubForce
     */
    SurfaceTensionSubForce<SPACE_DIM>* mpSurfaceTensionForce;

public:
    /**
     * Empty constructor.
     */
    FaceTensionWriter();

    /**
     * Default Constructor.
     */
    FaceTensionWriter(SurfaceTensionSubForce<SPACE_DIM>* pTensionForce);

    /**
     * Sets the pointer to save modifiert
     */
    void SetTensionForce(SurfaceTensionSubForce<SPACE_DIM>* pTensionForce);

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
EXPORT_TEMPLATE_CLASS2(FaceTensionWriter, 2, 2)
EXPORT_TEMPLATE_CLASS2(FaceTensionWriter, 2, 3)
EXPORT_TEMPLATE_CLASS2(FaceTensionWriter, 3, 3)

#endif /* FACETENSIONWRITER_HPP_ */