#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"

#include "CellProliferativeTypesWriter.hpp"
#include "CellThicknessWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FiniteThicknessSimulation3d.hpp"
#include "FiniteThicknessVertexMeshFromReaderGenerator.hpp"
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

        // MonolayerVertexMeshReader<3, 3> mesh_reader("/home/chaste/testoutput/write_mesh_test/random_mesh");
        MonolayerVertexMeshReader<3, 3> mesh_reader("/home/chaste/testoutput/final_mesh/evolved_random_sphere");

        FiniteThicknessVertexMeshFromReaderGenerator generator(mesh_reader, 0.0, 0.001);

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
        simulator.SetOutputDirectory("20230818_TestReadMesh");
        simulator.SetEndTime(5.0);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetDt(0.01);

        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_tension_subforce);

        p_force->AddSubForce(p_tension_subforce);
        p_force->SetSimulationInstance(&simulator);
        p_tension_subforce->SetSimulatedAnnealingParameters(0.0, 50.0, 0.0);
        p_tension_subforce->SetSimulationInstance(&simulator);

        p_tension_subforce->CreateSurfaceTensionParametersForCells(
            0.8, 0.8, 1.0, p_mesh); // apical, basal, lateral

        simulator.AddForce(p_force);

        double avg_cell_volume = 1.0;

        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_growth_modifier,
                      (&cell_population));
        p_growth_modifier->SetGrowthDuration(0.0);
        p_growth_modifier->SetT1AdaptationDuration(0.100);
        p_growth_modifier->SetReferenceTargetVolume(avg_cell_volume);
        simulator.AddSimulationModifier(p_growth_modifier);
        std::cout << "Reached point before Solve()" << std::endl;

        std::cout << "Lumen volume before initial relaxation: " << p_mesh->GetLumenVolume() << std::endl;

        double target_vol = 4.0 / 3.0 * M_PI * pow((-1.435 / 2.0 + sqrt(600.0 / 4.0 / M_PI / 1.435)), 3);
        // p_mesh->GetLumenVolume();
        std::cout << "Target Lumen Volume: " << target_vol << std::endl;
        cell_population.SetTargetLumenVolume(target_vol);
        p_force->SetConstantLumenVolumeFlag(true);

        simulator.Solve();
        std::cout << "Reached point after Solve()" << std::endl;
    }
};