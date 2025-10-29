#include "AbstractFaceWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "UblasVectorInclude.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>::AbstractFaceWriter(const std::string& rFileName)
        : AbstractCellBasedWriter<ELEMENT_DIM, SPACE_DIM>(rFileName),
          mOutputScalarData(true),
          mOutputVectorData(false),
          mVtkFaceDataName("DefaultVtkFaceDataName"),
          mVtkVectorFaceDataName("DefaultVtkVectorFaceDataName")
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>::GetOutputScalarData()
{
    return mOutputScalarData;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>::GetOutputVectorData()
{
    return mOutputVectorData;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>::SetVtkFaceDataName(std::string vtkFaceDataName)
{
    mVtkFaceDataName = vtkFaceDataName;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>::SetVtkVectorFaceDataName(std::string vtkFaceDataName)
{
    mVtkVectorFaceDataName = vtkFaceDataName;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>::GetFaceDataForVtkOutput(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace,
                                                                           AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    return DOUBLE_UNSET;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>::GetVectorFaceDataForVtkOutput(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace,
                                                                                                      AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    return scalar_vector<double>(SPACE_DIM, DOUBLE_UNSET);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>::GetVtkFaceDataName()
{
    return mVtkFaceDataName;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>::GetVtkVectorFaceDataName()
{
    return mVtkVectorFaceDataName;
}

// Explicit instantiation
template class AbstractFaceWriter<2, 2>;
template class AbstractFaceWriter<2, 3>;
template class AbstractFaceWriter<3, 3>;
