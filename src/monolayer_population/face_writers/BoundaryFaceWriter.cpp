#include "BoundaryFaceWriter.hpp"
#include "MonolayerVertexElement.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
BoundaryFaceWriter<ELEMENT_DIM, SPACE_DIM>::BoundaryFaceWriter()
        : AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>("facetypes.dat")
{
    this->mVtkFaceDataName = "Boundary Face";
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double BoundaryFaceWriter<ELEMENT_DIM, SPACE_DIM>::GetFaceDataForVtkOutput(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    return pFace->IsBoundaryFace() ? 1.0 : 0.0;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BoundaryFaceWriter<ELEMENT_DIM, SPACE_DIM>::VisitFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("VisitFace() is not implemented for BoundaryFaceWriter");
}

// Explicit instantiation

template class BoundaryFaceWriter<2, 2>;
template class BoundaryFaceWriter<2, 3>;
template class BoundaryFaceWriter<3, 3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS2(BoundaryFaceWriter, 2, 2)
EXPORT_TEMPLATE_CLASS2(BoundaryFaceWriter, 2, 3)
EXPORT_TEMPLATE_CLASS2(BoundaryFaceWriter, 3, 3)