#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"

#include "CellThicknessWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellsGenerator.hpp"
#include "FaceTypeWriter.hpp"
#include "FiniteThicknessRandomizedSphereMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

#include "CellProliferativeTypesWriter.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"

#include "DecayingRandomForce.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SurfaceTensionForce.hpp"

#include "SmartPointers.hpp"
#include "SubstrateBoundaryCondition.hpp"

#include "FakePetscSetup.hpp"
// Test for memory leaks
#include "PetscSetupAndFinalize.hpp"

#include <cmath> //for M_PI

class TestMonolayerVisualization : public AbstractCellBasedTestSuite
{
public:
    void Test3DVertexBasedSphereSimulationForVisualizing()
    {
        unsigned num_cells = 40;
        double height = 2.0 / 3.0 / sqrt(3.0) * cbrt((9.0 / 2.0) * (9.0 / 2.0)) * cbrt(1.0 / 0.35 * 1.0 / 0.35) * 1.0;
        // double area = 1.0/height;
        //  length at which T1s become favorable: l^*=0.66/3^(2/3)/gamma^(1/3)
        double t1_length = 0.66 / cbrt(3 * 3 * (1.0 + 1.0) / 0.7);
        FiniteThicknessRandomizedSphereMeshGenerator generator(num_cells, t1_length, 0.001, height, 3.0 * height);
        // generator.CurveMonolayer(2.0);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        // We use protorosettes as intermediate states
        p_mesh->SetProtorosetteFormationProbability(1.0);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);
        // p_mesh -> SetCheckForInternalIntersections(true);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 40u);

        // We generate the cell objects
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        for (unsigned cell_index = 0; cell_index < num_cells; cell_index++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel;
            p_model->SetDimension(3);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_model->SetStemCellG1Duration(25.0);
            p_model->SetTransitCellG1Duration(25.0);
            p_model->SetMaxTransitGenerations(1);
            p_model->SetSDuration(25.0);
            p_model->SetG2Duration(0.001);
            p_model->SetMDuration(0.001);
            double birth_time = -RandomNumberGenerator::Instance()->ranf() * 5.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // TS_ASSERT_EQUALS(cells.size(), 400u);

        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellRearrangementLocations(false);
        cell_population.SetRestrictVertexMovementBoolean(false);
        cell_population.SetDoInitialVolumeRelaxation(true);

        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellThicknessWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddFaceWriter<FaceTypeWriter>();

        FiniteThicknessSimulation3d simulator(cell_population);
        simulator.SetOutputDirectory("Test3DSphere_growth_active_t1_1transitgen_t1rate2");
        simulator.SetEndTime(0.300);
        simulator.SetSamplingTimestepMultiple(20);
        simulator.SetDt(0.001);

        MAKE_PTR(SurfaceTensionForce<3>, p_force);
        // MAKE_PTR(DecayingRandomForce<3>, p_random_force);
        p_force->CreateSurfaceTensionParametersForCells(1.0, 1.0, 0.7, p_mesh); // apical, basal, lateral
        p_force->SetSimulatedAnnealingParameters(0.0, 190.0, 10.0);
        p_force->SetSimulationInstance(&simulator);
        // p_force->SetPerformActiveT1Swaps();
        // p_random_force->SetRelativeStrengthParameter(5.0);
        // p_random_force->SetDecayTimeParameter(20.0);
        simulator.AddForce(p_force);
        // simulator.AddForce(p_random_force);

        double in_rad = 3.0 * height;
        double out_rad = 3.0 * height + height;
        double avg_cell_volume = 4.0 / 3.0 * M_PI * (out_rad * out_rad * out_rad - in_rad * in_rad * in_rad) / num_cells;

        // MAKE_PTR(SimpleTargetVolumeModifier<3>, p_growth_modifier);
        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_growth_modifier, (&cell_population));
        // p_growth_modifier->SetGrowthDuration(0.400);
        p_growth_modifier->SetT1AdaptationDuration(0.300);
        p_growth_modifier->SetReferenceTargetVolume(avg_cell_volume);
        // p_growth_modifier->SetT1AdaptationDuration(0.0);
        simulator.AddSimulationModifier(p_growth_modifier);

        // c_vector<double, 3> zero_point = zero_vector<double>(3);
        // c_vector<double, 3> normal_vector = zero_vector<double>(3);
        // normal_vector[2]=1.0;
        // MAKE_PTR_ARGS(SubstrateBoundaryCondition<3>, p_substrate_bc, (&cell_population, zero_point, normal_vector));
        // simulator.AddCellPopulationBoundaryCondition(p_substrate_bc);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 400u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.06, 1e-10);

        simulator.SetEndTime(45.0);
        simulator.SetOutputDirectory("Test3DSphere_growth_active_t1_1transitgen_t1rate2");
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetDt(0.006);

        p_force->SetSimulatedAnnealingParameters(0.001, 50.0, 0.5);
        p_growth_modifier->SetT1AdaptationDuration(0.100);
        p_force->SetSimulationInstance(&simulator);
        p_force->SetPerformActiveT1Swaps(true);

        simulator.Solve();

        for (auto it = cell_population.rGetCells().begin(); it != cell_population.rGetCells().begin(); ++it)
        {
            FixedG1GenerationalCellCycleModel* p_cell_cycle = static_cast<FixedG1GenerationalCellCycleModel*>((*it)->GetCellCycleModel());
            p_cell_cycle->SetStemCellG1Duration(25.0);
            p_cell_cycle->SetTransitCellG1Duration(25.0);
        }

        for (auto it = cell_population.rGetCells().begin(); it != cell_population.rGetCells().begin(); ++it)
        {
            TS_ASSERT_DELTA(static_cast<FixedG1GenerationalCellCycleModel*>((*it)->GetCellCycleModel())->GetStemCellG1Duration(), 0.018, 1e-7);
            TS_ASSERT_DELTA(static_cast<FixedG1GenerationalCellCycleModel*>((*it)->GetCellCycleModel())->GetTransitCellG1Duration(), 0.018, 1e-7);
        }

        simulator.SetEndTime(350.0);
        simulator.SetOutputDirectory("Test3DSphere_growth_active_t1_1transitgen_t1rate2");
        simulator.SetSamplingTimestepMultiple(400);
        simulator.SetDt(0.006);

        p_force->SetSimulatedAnnealingParameters(0.003, 1900000.0, 1.0);
        p_force->SetSimulationInstance(&simulator);
        p_force->SetPerformActiveT1Swaps(true);
        p_force->SetT1TransitionParameters(2.0, false);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 400u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 300.0, 1e-10);
    }
};