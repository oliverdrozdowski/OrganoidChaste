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
#include "LumenPressureSubForce.hpp"
#include "MovingBoundarySubForce.hpp"
#include "RandomSubForce.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SurfaceTensionSubForce.hpp"

#include "SmartPointers.hpp"
#include "SubstrateBoundaryCondition.hpp"

#include "ArchiveOpener.hpp"
#include "MonolayerVertexMeshReader.hpp"
#include "MonolayerVertexMeshWriter.hpp"

#include "OneStepProfile.hpp"
#include "PlateMovingBoundary.hpp"

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
        MonolayerVertexMeshReader<3, 3> mesh_reader("/home/chaste/mesh/N_700_G09/final_mesh/evolved_random_sphere");

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
        cell_population.SetDoInitialVolumeRelaxation(false);

        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellThicknessWriter>();
        cell_population.AddFaceWriter<FaceTypeWriter>();
        cell_population.AddMonolayerPopulationWriter<PopulationLumenVolumeWriter>();
        cell_population.AddMonolayerPopulationWriter<BoundaryHeightWriter>();
        cell_population.AddMonolayerPopulationWriter<PopulationBoundaryForceWriter>();

        FiniteThicknessSimulation3d simulator(cell_population);

        time_t rawtime;
        struct tm* timeinfo;
        char buffer[80];
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(buffer, sizeof(buffer), "%Y%m%d%H%M%S", timeinfo);
        std::string time_string(buffer);
        std::stringstream stream;
        stream << time_string << "_test_pressurized_sphere";
        std::string outputPath{ stream.str() };
        simulator.SetOutputDirectory(outputPath);

        simulator.SetEndTime(5.0);
        simulator.SetSamplingTimestepMultiple(5);
        simulator.SetDt(0.005);

        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_tension_subforce);

        p_force->AddSubForce(p_tension_subforce);
        p_force->SetSimulationInstance(&simulator);
        p_tension_subforce->SetSimulatedAnnealingParameters(0.0, 50.0, 0.0);
        p_tension_subforce->SetSimulationInstance(&simulator);

        p_tension_subforce->CreateSurfaceTensionParametersForCells(
            0.9, 0.9, 1.0, p_mesh); // apical, basal, lateral

        // MAKE_PTR_ARGS(LumenPressureSubForce<3>, p_lumen_pressure_subforce, (2.8));
        // p_force->AddSubForce(p_lumen_pressure_subforce);
        // Temperature temp = Temperature();

        // MAKE_PTR_ARGS(RandomSubForce<3>, p_random_subforce, (1.0));
        //
        // p_force->AddSubForce(p_random_subforce);

        simulator.AddForce(p_force);

        double avg_cell_volume = 1.0;

        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_growth_modifier, (&cell_population));
        p_growth_modifier->SetGrowthDuration(0.0);
        p_growth_modifier->SetT1AdaptationDuration(0.100);
        p_growth_modifier->SetReferenceTargetVolume(avg_cell_volume);
        simulator.AddSimulationModifier(p_growth_modifier);
        std::cout << "Reached point before Solve()" << std::endl;

        std::cout << "Lumen volume before initial relaxation: " << p_mesh->GetLumenVolume() << std::endl;

        double height = cbrt(4.0 * 4.0) * cbrt(3.0 * sqrt(3.0) / 2.0) / cbrt(6.0 * 6.0) * cbrt((0.9 + 0.9) * (0.9 + 0.9));
        int num_cells = 200;
        double inner_radius = sqrt((double)num_cells / 4.0 / M_PI / height) - height / 2.0;
        double target_vol = 4.0 / 3.0 * M_PI * pow(inner_radius, 3);
        std::cout << "Target Lumen Volume: " << target_vol << std::endl;
        cell_population.SetTargetLumenVolume(target_vol);
        p_force->SetConstantLumenVolumeFlag(true);

        // p_mesh->SetLengthDependentActiveT1SwapConstant(0.001);
        // p_mesh->InitializeTemperature();

        simulator.Solve();

        std::cout << "Reached point after Solve()" << std::endl;
    }
};