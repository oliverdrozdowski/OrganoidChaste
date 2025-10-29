#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"

#include "CellLabelWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellThicknessWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FaceTypeWriter.hpp"
#include "FiniteThicknessRandomizedSphereMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"
#include "MechanicalMosaicCellMutationState.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

#include "DecayingRandomForce.hpp"
#include "GeneralizedVolumeConservingForce.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SurfaceEvolverSaveModifier.hpp"
#include "SurfaceTensionSubForce.hpp"

#include "SmartPointers.hpp"
#include "SubstrateBoundaryCondition.hpp"

#include "FakePetscSetup.hpp"
// Test for memory leaks
#include "PetscSetupAndFinalize.hpp"

#include <cmath> //for M_PI

class TestMonolayerVisualization : public AbstractCellBasedTestSuite
{
public:
    void Test3DVertexBasedSphereSimulationHeterogeneousSurfaceEvolver()
    {
        unsigned num_cells = 400;
        double height = 2.0 / 3.0 / sqrt(3.0) * cbrt((9.0 / 2.0) * (9.0 / 2.0)) * cbrt(1.0 / 0.35 * 1.0 / 0.35) * 1.0;
        // double area = 1.0/height;
        //  length at which T1s become favorable: l^*=0.66/3^(2/3)/gamma^(1/3)
        double t1_length = 0.66 / cbrt(3 * 3 * (1.0 + 1.0) / 0.7);
        FiniteThicknessRandomizedSphereMeshGenerator generator(num_cells, t1_length, 0.001, height, 7.0 * height);
        // generator.CurveMonolayer(2.0);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        // We use protorosettes as intermediate states
        p_mesh->SetProtorosetteFormationProbability(1.0);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);
        // p_mesh -> SetCheckForInternalIntersections(true);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), num_cells);

        // We generate the cell objects
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state_wt);
        MAKE_PTR(MechanicalMosaicCellMutationState, p_state_mosaic);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        for (unsigned cell_index = 0; cell_index < num_cells; cell_index++)
        {
            c_vector<double, 3> centr = p_mesh->GetCentroidOfElement(cell_index);
            NoCellCycleModel* p_model = new NoCellCycleModel;
            p_model->SetDimension(3);

            CellPtr p_cell;
            // Contraction for large z.
            if (centr(2) < 7.5 * 0.8 * height)
                p_cell = CellPtr(new Cell(p_state_wt, p_model));
            else
                p_cell = CellPtr(new Cell(p_state_mosaic, p_model));

            p_cell->SetCellProliferativeType(p_differentiated_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf() * 5.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellRearrangementLocations(false);
        cell_population.SetRestrictVertexMovementBoolean(false);
        cell_population.SetDoInitialVolumeRelaxation(true);

        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellThicknessWriter>();
        cell_population.AddCellWriter<CellLabelWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddFaceWriter<FaceTypeWriter>();

        FiniteThicknessSimulation3d simulator(cell_population);
        simulator.SetOutputDirectory("Test3DSphere_budding_apical_constr_apsftns1.2_surfev");
        simulator.SetEndTime(0.120);
        simulator.SetSamplingTimestepMultiple(5);
        simulator.SetDt(0.003);

        // Make a dictionary for mosaic cell surface tension parameters
        std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3> > dictionary_surface_tensions;
        std::array<double, 3> wild_type_surface_tensions = { 1.0, 1.0, 0.7 }; // apical, basal, lateral
        std::array<double, 3> mosaic_surface_tensions = { 1.2, 1.0, 0.7 };

        dictionary_surface_tensions[p_state_wt] = wild_type_surface_tensions;
        dictionary_surface_tensions[p_state_mosaic] = mosaic_surface_tensions;

        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_subforce);
        p_force->AddSubForce(p_subforce);
        p_force->SetSimulationInstance(&simulator);

        // MAKE_PTR(SurfaceTensionForce<3>, p_force);
        // MAKE_PTR(DecayingRandomForce<3>, p_random_force);
        // p_force->CreateSurfaceTensionParametersForCells(1.0, 1.0, 0.7, p_mesh); // apical, basal, lateral
        p_subforce->SetSurfaceTensionParametersByMutation(dictionary_surface_tensions);
        p_subforce->SetSimulatedAnnealingParameters(0.0, 190.0, 10.0);
        p_subforce->SetSimulationInstance(&simulator);
        // p_force->SetPerformActiveT1Swaps();
        // p_random_force->SetRelativeStrengthParameter(5.0);
        // p_random_force->SetDecayTimeParameter(20.0);
        simulator.AddForce(p_force);
        // simulator.AddForce(p_random_force);

        // Make the map from surface tensions to colors in surface evolver
        std::map<double, int> color_map;
        color_map[1.0] = 1;
        color_map[0.7] = 2;
        color_map[1.2] = 3;

        MAKE_PTR(SurfaceEvolverSaveModifier<3>, p_surface_evolver_modifier);
        p_surface_evolver_modifier->SetSaveAtEnd();
        p_surface_evolver_modifier->SetMapTensionToColor(color_map);
        p_surface_evolver_modifier->SetSurfaceTensionSubForce(p_subforce);
        simulator.AddSimulationModifier(p_surface_evolver_modifier);

        double in_rad = 7.0 * height;
        double out_rad = 7.0 * height + height;
        double avg_cell_volume = 4.0 / 3.0 * M_PI * (out_rad * out_rad * out_rad - in_rad * in_rad * in_rad) / num_cells;

        // MAKE_PTR(SimpleTargetVolumeModifier<3>, p_growth_modifier);
        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_growth_modifier, (&cell_population));
        p_growth_modifier->SetGrowthDuration(0.0);
        p_growth_modifier->SetT1AdaptationDuration(0.100);
        p_growth_modifier->SetReferenceTargetVolume(avg_cell_volume);
        // p_growth_modifier->SetT1AdaptationDuration(0.0);
        simulator.AddSimulationModifier(p_growth_modifier);

        // c_vector<double, 3> zero_point = zero_vector<double>(3);
        // c_vector<double, 3> normal_vector = zero_vector<double>(3);
        // normal_vector[2]=1.0;
        // MAKE_PTR_ARGS(SubstrateBoundaryCondition<3>, p_substrate_bc, (&cell_population, zero_point, normal_vector));
        // simulator.AddCellPopulationBoundaryCondition(p_substrate_bc);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.12, 1e-10);

        simulator.SetEndTime(400.0);
        simulator.SetOutputDirectory("Test3DSphere_budding_apical_constr_apsftns1.2_surfev");
        simulator.SetSamplingTimestepMultiple(400);
        simulator.SetDt(0.02);

        p_subforce->SetSimulatedAnnealingParameters(0.1, 300.0, 5.0);
        p_subforce->SetSimulationInstance(&simulator);
        p_subforce->SetPerformActiveT1Swaps();
        p_subforce->SetT1TransitionParameters(50.0, true);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 1000.0, 1e-10);
    }
};