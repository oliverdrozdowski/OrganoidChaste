#include "FaceIndexWriter.hpp"
#include "MonolayerVertexElement.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FaceIndexWriter<ELEMENT_DIM, SPACE_DIM>::FaceIndexWriter()
        : AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>("faceindices.dat")
{
    this->mVtkFaceDataName = "Face index";
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double FaceIndexWriter<ELEMENT_DIM, SPACE_DIM>::GetFaceDataForVtkOutput(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    return double(pFace->GetIndex());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FaceIndexWriter<ELEMENT_DIM, SPACE_DIM>::VisitFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("VisitFace() is not implemented for FaceIndexWriter");
}

// Explicit instantiation

template class FaceIndexWriter<2, 2>;
template class FaceIndexWriter<2, 3>;
template class FaceIndexWriter<3, 3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS2(FaceIndexWriter, 2, 2)
EXPORT_TEMPLATE_CLASS2(FaceIndexWriter, 2, 3)
EXPORT_TEMPLATE_CLASS2(FaceIndexWriter, 3, 3)