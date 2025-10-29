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
#include "LumenPressureSubForce.hpp"
#include "MovingBoundarySubForce.hpp"
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
        double gamma_a = 0.8;
        double gamma_b = 0.8;
        int num_cells = 300;

        double height = cbrt(4.0 * 4.0) * cbrt(3.0 * sqrt(3.0) / 2.0) / cbrt(6.0 * 6.0) * cbrt((gamma_a + gamma_b) * (gamma_a + gamma_b));
        double avg_cell_volume = 1.0;

        // for this h and N. make a bit larger because the random distribution of
        // points would fail otherwise
        double true_inner_radius = sqrt(num_cells / 4.0 / M_PI / height); // - height / 2.0;
        double inner_radius = true_inner_radius * 1.2;

        //  length at which T1s become favorable: l^*=0.66/3^(2/3)/gamma^(1/3)
        double t1_length = 0.0;

        // Generate and configure Mesh
        FiniteThicknessRandomizedSphereMeshGenerator generator(
            true, num_cells, t1_length, 0.001, height, inner_radius);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        // We use protorosettes as intermediate states
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

        FiniteThicknessSimulation3d simulator(cell_population);

        time_t rawtime;
        struct tm* timeinfo;
        char buffer[80];
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(buffer, sizeof(buffer), "%Y%m%d%H%M%S", timeinfo);
        std::string time_string(buffer);
        std::stringstream stream;
        stream << time_string << "_test_evolve_sphere";
        std::string outputPath{ stream.str() };
        simulator.SetOutputDirectory(outputPath);

        simulator.SetEndTime(5.0);
        simulator.SetSamplingTimestepMultiple(20);
        simulator.SetDt(0.001);

        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_tension_subforce);

        p_force->AddSubForce(p_tension_subforce);
        p_force->SetSimulationInstance(&simulator);
        p_tension_subforce->SetSimulatedAnnealingParameters(0.0, 50.0, 0.0);
        p_tension_subforce->SetSimulationInstance(&simulator);

        p_tension_subforce->CreateSurfaceTensionParametersForCells(
            gamma_a, gamma_b, 1.0, p_mesh); // apical, basal, lateral

        // MAKE_PTR_ARGS(LumenPressureSubForce<3>, p_lumen_pressure_subforce, (0.9));
        // p_force->AddSubForce(p_lumen_pressure_subforce);

        simulator.AddForce(p_force);

        // Add moving boundaries

        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_growth_modifier,
                      (&cell_population));
        p_growth_modifier->SetGrowthDuration(0.0);
        p_growth_modifier->SetT1AdaptationDuration(0.100);
        p_growth_modifier->SetReferenceTargetVolume(avg_cell_volume);
        simulator.AddSimulationModifier(p_growth_modifier);
        std::cout << "Reached point before Solve()" << std::endl;

        std::cout << "Lumen volume before initial relaxation: " << p_mesh->GetLumenVolume() << std::endl;

        double target_vol = p_mesh->GetLumenVolume();
        std::cout << "Target Lumen Volume: " << target_vol << std::endl;
        cell_population.SetTargetLumenVolume(target_vol);
        p_force->SetConstantLumenVolumeFlag(true);

        simulator.Solve();
    }
};