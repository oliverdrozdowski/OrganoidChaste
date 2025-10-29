#include "FaceTensionWriter.hpp"
#include "MonolayerVertexElement.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FaceTensionWriter<ELEMENT_DIM, SPACE_DIM>::FaceTensionWriter()
        : AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>("facetensions.dat")
{
    this->mpSurfaceTensionForce = nullptr;
    this->mVtkFaceDataName = "Face tension";
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FaceTensionWriter<ELEMENT_DIM, SPACE_DIM>::FaceTensionWriter(SurfaceTensionSubForce<SPACE_DIM>* pTensionForce)
        : AbstractFaceWriter<ELEMENT_DIM, SPACE_DIM>("facetensions.dat")
{
    this->mpSurfaceTensionForce = pTensionForce;
    this->mVtkFaceDataName = "Face tension";
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FaceTensionWriter<ELEMENT_DIM, SPACE_DIM>::SetTensionForce(SurfaceTensionSubForce<SPACE_DIM>* pTensionForce)
{
    this->mpSurfaceTensionForce = pTensionForce;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double FaceTensionWriter<2, 2>::GetFaceDataForVtkOutput(MonolayerVertexElement<1, 2>* pFace, AbstractCellPopulation<2, 2>* pCellPopulation)
{
    EXCEPTION("FaceTensionWriter only in 3D with 3D elements.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

/// \cond Get Doxygen to ignore, since it's confused by these templates

template <>
double FaceTensionWriter<2, 3>::GetFaceDataForVtkOutput(MonolayerVertexElement<1, 3>* pFace, AbstractCellPopulation<2, 3>* pCellPopulation)
{
    EXCEPTION("FaceTensionWriter only in 3D with 3D elements.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double FaceTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetFaceDataForVtkOutput(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    if (mpSurfaceTensionForce == nullptr)
    {
        EXCEPTION("SurfaceTensionSubForce must be set initially when FaceTensionWriter is being used!");
    }
    mpSurfaceTensionForce->UpdateSurfaceTensions(pCellPopulation);

    // Get the global index and then look up the tension
    unsigned face_global_index = pFace->GetIndex();
    double tension = mpSurfaceTensionForce->GetSurfaceTensionParameter(face_global_index);
    // double counting as laterals are saved with 1/2 in force
    double count_factor = pFace->GetFaceType() == MonolayerVertexElementType::Lateral ? 2.0 : 1.0;
    return count_factor * tension;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FaceTensionWriter<ELEMENT_DIM, SPACE_DIM>::VisitFace(MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("VisitFace() is not implemented for FaceTensionWriter");
}

// Explicit instantiation

template class FaceTensionWriter<2, 2>;
template class FaceTensionWriter<2, 3>;
template class FaceTensionWriter<3, 3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS2(FaceTensionWriter, 2, 2)
EXPORT_TEMPLATE_CLASS2(FaceTensionWriter, 2, 3)
EXPORT_TEMPLATE_CLASS2(FaceTensionWriter, 3, 3)