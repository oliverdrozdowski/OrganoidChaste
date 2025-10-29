#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"

#include "CellProliferativeTypesWriter.hpp"
#include "CellThicknessWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FiniteThicknessRandomizedSphereMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

#include "BoundaryHeightWriter.hpp"
#include "FaceTypeWriter.hpp"
#include "PopulationBoundaryForceWriter.hpp"
#include "PopulationLumenVolumeWriter.hpp"

#include "DecayingRandomForce.hpp"
#include "GeneralizedVolumeConservingForce.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SurfaceTensionSubForce.hpp"

#include "SmartPointers.hpp"
#include "SubstrateBoundaryCondition.hpp"

#include "ArchiveOpener.hpp"
#include "MonolayerVertexMeshReader.hpp"
#include "MonolayerVertexMeshWriter.hpp"

#include "FakePetscSetup.hpp"
// Test for memory leaks
#include "PetscSetupAndFinalize.hpp"

#include <cmath> //for M_PI
#include <ctime> // for current timestamp

class TestMonolayerVisualization : public AbstractCellBasedTestSuite
{
public:
    void Test3DVertexBasedSphereSimulationForVisualizing()
    {

        double gamma_a = 0.7;
        double gamma_b = 0.7;
        double height = cbrt(4.0 * 4.0) * cbrt(3.0 * sqrt(3.0) / 2.0) / cbrt(6.0 * 6.0) * cbrt((gamma_a + gamma_b) * (gamma_a + gamma_b));
        unsigned int num_elements{ 500 };
        double inner_radius = (-height / 2.0 + sqrt(num_elements / 4.0 / M_PI / height)) * 1.1;
        double t1_length = 0;

        FiniteThicknessRandomizedSphereMeshGenerator generator(
            true, num_elements, t1_length, 0.001, height, inner_radius);

        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        p_mesh->SetProtorosetteFormationProbability(1.0);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellRearrangementLocations(false);
        cell_population.SetRestrictVertexMovementBoolean(false);
        cell_population.SetDoInitialVolumeRelaxation(true);

        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellThicknessWriter>();
        cell_population.AddFaceWriter<FaceTypeWriter>();
        cell_population.AddMonolayerPopulationWriter<PopulationLumenVolumeWriter>();
        cell_population.AddMonolayerPopulationWriter<BoundaryHeightWriter>();
        cell_population
            .AddMonolayerPopulationWriter<PopulationBoundaryForceWriter>();

        FiniteThicknessSimulation3d simulator(cell_population);
        simulator.SetOutputDirectory("TestWriteReadMesh3");
        simulator.SetEndTime(0.1);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetDt(0.05);

        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_tension_subforce);

        p_force->AddSubForce(p_tension_subforce);
        p_force->SetSimulationInstance(&simulator);
        p_tension_subforce->SetSimulatedAnnealingParameters(0.0, 50.0, 0.0);
        p_tension_subforce->SetSimulationInstance(&simulator);

        p_tension_subforce->CreateSurfaceTensionParametersForCells(
            gamma_a, gamma_b, 1.0, p_mesh); // apical, basal, lateral

        simulator.AddForce(p_force);

        double avg_cell_volume = 1.0;

        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_growth_modifier,
                      (&cell_population));
        p_growth_modifier->SetGrowthDuration(0.0);
        p_growth_modifier->SetT1AdaptationDuration(0.100);
        p_growth_modifier->SetReferenceTargetVolume(avg_cell_volume);
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.Solve();

        MonolayerVertexMeshWriter<3, 3> mesh_writer = MonolayerVertexMeshWriter<3, 3>("write_mesh_test", "random_mesh");

        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }
};