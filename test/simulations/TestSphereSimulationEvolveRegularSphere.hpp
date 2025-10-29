#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "BoundaryHeightWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellThicknessWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FaceTypeWriter.hpp"
#include "FiniteThicknessRandomizedSphereMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"
#include "GeneralizedVolumeConservingForce.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "PopulationLumenVolumeWriter.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SmartPointers.hpp"
#include "SurfaceTensionSubForce.hpp"
#include "FakePetscSetup.hpp"
// Test for memory leaks
#include <cmath> //for M_PI
#include <ctime> // for current timestamp

#include "PetscSetupAndFinalize.hpp"

class TestMonolayerVisualization : public AbstractCellBasedTestSuite
{
private:
    std::string generateOutputPath(int num_elements, double gamma)
    {
        time_t rawtime;
        struct tm* timeinfo;
        char buffer[80];
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(buffer, sizeof(buffer), "%Y%m%d%H%M%S", timeinfo);
        std::string time_string(buffer);
        std::stringstream stream;
        stream << time_string << "_numElements_" << num_elements << "_gamma_"
               << gamma;
        std::string outputPath{ stream.str() };

        return outputPath;
    }

public:
    void Test3DVertexBasedSphereSimulationForVisualizing()
    {

        double gamma = 0.8;
        unsigned numElements = 500;

        // std::string outputPath = generateOutputPath(numElements, gamma);

        double height = cbrt(4.0 * 4.0) * cbrt(3.0 * sqrt(3.0) / 2.0) / cbrt(6.0 * 6.0) * cbrt((gamma + gamma) * (gamma + gamma));
        // double avg_cell_volume = 1.0;
        double inner_radius = sqrt(numElements / 4.0 / M_PI / height) - height / 2.0;

        for (int i = 1; i < 7; i++)
        {
            FiniteThicknessRandomizedSphereMeshGenerator generator(
                numElements, 0.0, 0.001, height, inner_radius, true);

            MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

            unsigned num_elements = p_mesh->GetNumElements();
            double num_hex = 0;

            for (unsigned elem_ind = 0; elem_ind < num_elements; elem_ind++)
            {
                MonolayerVertexElement<3, 3>* elem = p_mesh->GetElement(elem_ind);

                int num_neighbours = (double)elem->GetNumFaces() - 2;
                if (num_neighbours == 6)
                    num_hex++;
            }
            std::cout << "P_6: " << num_hex / (double)num_elements << std::endl;
        }
        /*
                FiniteThicknessRandomizedSphereMeshGenerator generator(
                    numElements, 0.0, 0.001, height, inner_radius, true);

                MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
                // We use protorosettes as intermediate states
                p_mesh->SetProtorosetteFormationProbability(1.0);
                p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);
                // p_mesh -> SetCheckForInternalIntersections(true);
                std::vector<CellPtr> cells;
                CellsGenerator<NoCellCycleModel, 3> cells_generator;
                cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

                MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
                cell_population.SetOutputCellRearrangementLocations(false);
                cell_population.SetRestrictVertexMovementBoolean(false);
                cell_population.SetDoInitialVolumeRelaxation(true);

                cell_population.AddCellWriter<CellVolumesWriter>();
                cell_population.AddCellWriter<CellThicknessWriter>();

                FiniteThicknessSimulation3d simulator(cell_population);
                simulator.SetOutputDirectory(outputPath);
                simulator.SetEndTime(0.1);
                simulator.SetSamplingTimestepMultiple(1);
                simulator.SetDt(0.1);

                // Create Forces
                MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
                MAKE_PTR(SurfaceTensionSubForce<3>, p_tension_subforce);
                p_force->AddSubForce(p_tension_subforce);
                p_force->SetSimulationInstance(&simulator);
                p_tension_subforce->SetSimulationInstance(&simulator);
                p_tension_subforce->CreateSurfaceTensionParametersForCells(gamma, gamma, 1.0, p_mesh); // apical, basal, lateral
                simulator.AddForce(p_force);

                // Add Modifier
                MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_growth_modifier, (&cell_population));
                p_growth_modifier->SetGrowthDuration(0.0);
                p_growth_modifier->SetT1AdaptationDuration(0.100);
                p_growth_modifier->SetReferenceTargetVolume(avg_cell_volume);
                simulator.AddSimulationModifier(p_growth_modifier);

                simulator.Solve();
                */
    }
};