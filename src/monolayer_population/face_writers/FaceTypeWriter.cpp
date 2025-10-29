#include "FaceTypeWriter.hpp"
#include "MonolayerVertexElement.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FaceTypeWriter<ELEMENT_DIM, SPACE_DIM>::FaceTypeWriter()
        : AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>("facetypes.dat")
{
    this->mVtkFaceDataName = "Face types";
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double FaceTypeWriter<ELEMENT_DIM, SPACE_DIM>::GetFaceDataForVtkOutput(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // Map the face type to unsigned ints, using the underlying value
    MonolayerVertexElementType face_type_abstract = pFace->GetFaceType();
    unsigned face_type_unsigned = static_cast<unsigned>(face_type_abstract);

    return double(face_type_unsigned);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FaceTypeWriter<ELEMENT_DIM, SPACE_DIM>::VisitFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("VisitFace() is not implemented for FaceTypeWriter");
}

// Explicit instantiation

template class FaceTypeWriter<2, 2>;
template class FaceTypeWriter<2, 3>;
template class FaceTypeWriter<3, 3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS2(FaceTypeWriter, 2, 2)
EXPORT_TEMPLATE_CLASS2(FaceTypeWriter, 2, 3)
EXPORT_TEMPLATE_CLASS2(FaceTypeWriter, 3, 3)