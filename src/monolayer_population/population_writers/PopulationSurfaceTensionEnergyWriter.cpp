#include "PopulationSurfaceTensionEnergyWriter.hpp"
#include <iomanip>
#include "MonolayerVertexBasedCellPopulation.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PopulationSurfaceTensionEnergyWriter<ELEMENT_DIM, SPACE_DIM>::PopulationSurfaceTensionEnergyWriter()
        : AbstractMonolayerVertexPopulationWriter<ELEMENT_DIM, SPACE_DIM>("CellPopulationSurfaceTensionEnergy.dat")
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationSurfaceTensionEnergyWriter<ELEMENT_DIM, SPACE_DIM>::SetTensionForcePointer(boost::shared_ptr<SurfaceTensionSubForce<SPACE_DIM> > pForce)
{
    this->mpForce = pForce;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PopulationSurfaceTensionEnergyWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MonolayerVertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE
    if (!mpForce)
    {
        EXCEPTION("To use the PopulationSurfaceTensionEnergyWriter one first needs to set the pointer to the SurfaceTensionSubForce instance!");
    }

    double energy = mpForce->CalculateTotalEnergy(*pCellPopulation);
    *this->mpOutStream << std::setprecision(17) << energy << " ";
}

// Explicit instantiation
template class PopulationSurfaceTensionEnergyWriter<1, 1>;
template class PopulationSurfaceTensionEnergyWriter<1, 2>;
template class PopulationSurfaceTensionEnergyWriter<2, 2>;
template class PopulationSurfaceTensionEnergyWriter<1, 3>;
template class PopulationSurfaceTensionEnergyWriter<2, 3>;
template class PopulationSurfaceTensionEnergyWriter<3, 3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PopulationSurfaceTensionEnergyWriter)