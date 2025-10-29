#ifndef ABSTRACTFACEWRITER_HPP_
#define ABSTRACTFACEWRITER_HPP_

#include <boost/serialization/base_object.hpp>
#include "AbstractCellBasedWriter.hpp"
#include "Cell.hpp"
#include "ChasteSerialization.hpp"
#include "MonolayerVertexElement.hpp"
#include "UblasVectorInclude.hpp"

// Forward declaration prevents circular include chain
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractCellPopulation;

// Forward declaration prevents circular include chain
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonolayerVertexElement;

/**
 * An abstract class for a writer that visits individual faces of a MonolayerVertexMesh and writes their data.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractFaceWriter : public AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>
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
        archive& boost::serialization::base_object<AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mOutputScalarData;
        archive & mOutputVectorData;
        archive & mVtkFaceDataName;
        archive & mVtkVectorFaceDataName;
    }

protected:
    /** Whether to output scalar data for VTK using GetCellDataForVtkOutput(). Default true. */
    bool mOutputScalarData;

    /** Whether to output scalar data for VTK using GetVectorCellDataForVtkOutput(). Default false. */
    bool mOutputVectorData;

    /** The name of the cell data used in VTK output. */
    std::string mVtkFaceDataName;

    /** The name of the vector cell data used in VTK output. */
    std::string mVtkVectorFaceDataName;

public:
    /**
     * Default constructor.
     * @param rFileName the name of the file to write to.
     */
    AbstractFaceWriter(const std::string& rFileName);

    /**
     * Get a double associated with a face. This method reduces duplication
     * of code between the methods VisitCell() and AddVtkData().
     *
     * By default this method returns a DOUBLE_UNSET, but it may be overridden in subclasses
     *
     * @param pFace a pointer to face
     * @param pCellPopulation a pointer to the cell population owning the cell.
     *
     * @return data associated with the cell
     */
    virtual double GetFaceDataForVtkOutput(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Get a c_vector associated with a cell face. This method reduces duplication
     * of code between the methods VisitCell() and AddVtkData().
     *
     * By default this method returns a c_vector of DOUBLE_UNSET, but it may be overridden in subclasses
     *
     * @param pFace a pointer to face
     * @param pCellPopulation a pointer to the cell population owning the cell.
     *
     * @return data associated with the cell
     */
    virtual c_vector<double, SPACE_DIM> GetVectorFaceDataForVtkOutput(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit a face and write its data.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param pFace a pointer to Fface
     * @param pCellPopulation a pointer to the cell population owning the cell.
     */
    virtual void VisitFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation) = 0;

    /**
     * Get whether to invoke GetCellDataForVtkOutput()
     *
     * @return mOutputScalarData
     */
    bool GetOutputScalarData();

    /**
     * Get whether to invoke GetVectorCellDataForVtkOutput()
     *
     * @return mOutputScalarData
     */
    bool GetOutputVectorData();

    /**
     * Set the name of the scalar cell data used in VTK output.
     * This method allows the user to change mVtkCellDataName from
     * its default value, which is set in each subclass's
     * constructor.
     *
     * @param vtkFaceDataName the name of the VTK field
     */
    void SetVtkFaceDataName(std::string vtkFaceDataName);

    /**
     * Set the name of the vector cell data used in VTK output.
     * This method allows the user to change mVtkCellDataName from
     * its default value, which is set in each subclass's
     * constructor.
     *
     * @param vtkFaceDataName the name of the VTK field
     */
    void SetVtkVectorFaceDataName(std::string vtkFaceDataName);

    /**
     * @return the name of the scalar cell data used in VTK output.
     */
    std::string GetVtkFaceDataName();

    /**
     * @return the name of the vector cell data used in VTK output.
     */
    std::string GetVtkVectorFaceDataName();
};

#endif /*ABSTRACTFACEWRITER_HPP_*/
