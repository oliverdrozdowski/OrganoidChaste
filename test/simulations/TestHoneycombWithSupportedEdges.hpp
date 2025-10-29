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
#include "FiniteThicknessHoneycombVertexMeshGenerator.hpp"
#include "FiniteThicknessIcosahedralSphereMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"
#include "MechanicalMosaicCellMutationState.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

#include "DecayingRandomForce.hpp"
#include "GeneralizedVolumeConservingForce.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "RectangularEdgeForce.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SurfaceEvolverSaveModifier.hpp"
#include "SurfaceTensionSubForce.hpp"

#include "SimplySupportedEdgeBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "SubstrateBoundaryCondition.hpp"

#include "FakePetscSetup.hpp"
// Test for memory leaks
#include "PetscSetupAndFinalize.hpp"

#include <cmath> //for M_PI

class TestFlatLatticeSurfaceEvolver : public AbstractCellBasedTestSuite
{
public:
    void Test3DFlatLatticceSupportedEdgeBoundaryCondition()
    {
        // alpha = apical/lateral surface tension
        // beta = basal/lateral surface tension
        double alpha = 1.0;
        double beta = 1.0;

        double height = cbrt(2.0 / sqrt(3.0)) * cbrt((alpha + beta) * (alpha + beta));
        unsigned num_cells_x = 9;
        unsigned num_cells_y = 9;
        unsigned num_cells = num_cells_x * num_cells_y;

        // length at which T1s become favorable: l^*=0.66/3^(2/3)/gamma^(1/3)
        double t1_length = 0.66 / cbrt(3. * 3. * 2. / (alpha + beta));
        t1_length = 0.0;
        FiniteThicknessHoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false,
                                                              t1_length, 0.001, height, 1.0 / height);
        generator.CurveMonolayer(40.0);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        num_cells = p_mesh->GetNumElements();
        // We use protorosettes as intermediate states
        p_mesh->SetProtorosetteFormationProbability(1.0);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);
        // p_mesh -> SetCheckForInternalIntersections(true);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), num_cells);

        // Get the edge nodes
        std::array<std::vector<Node<3>*>, 4> edge_array = generator.GetEdgesWithNodes();
        std::vector<Node<3>*> left_edge_nodes = edge_array[0];
        std::vector<Node<3>*> right_edge_nodes = edge_array[2];
        c_vector<double, 3> left_force = zero_vector<double>(3);
        c_vector<double, 3> right_force = zero_vector<double>(3);

        // We generate the cell objects
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state_wt);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        for (unsigned cell_index = 0; cell_index < num_cells; cell_index++)
        {
            c_vector<double, 3> centr = p_mesh->GetCentroidOfElement(cell_index);
            NoCellCycleModel* p_model = new NoCellCycleModel;
            p_model->SetDimension(3);

            CellPtr p_cell;

            // Randomly assign cell state (i.e. tensions)
            // unsigned rand_int = RandomNumberGenerator::Instance()->randMod(4);
            p_cell = CellPtr(new Cell(p_state_wt, p_model));
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
        simulator.SetOutputDirectory("Test3D_Honeycomb_SupportedEdges_3");
        simulator.SetEndTime(40.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetDt(0.003);

        // Make a dictionary for mosaic cell surface tension parameters
        std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3> > dictionary_surface_tensions;
        std::array<double, 3> wild_type_surface_tensions = { alpha, beta, 1.0 }; // apical, basal, lateral

        dictionary_surface_tensions[p_state_wt] = wild_type_surface_tensions;

        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_subforce);
        p_force->AddSubForce(p_subforce);
        p_force->SetSimulationInstance(&simulator);

        p_subforce->SetSurfaceTensionParametersByMutation(dictionary_surface_tensions);
        p_subforce->SetSimulatedAnnealingParameters(0.1, 50.0, 0.0);
        p_subforce->SetSimulationInstance(&simulator);
        // p_force->SetPerformActiveT1Swaps();
        // p_random_force->SetRelativeStrengthParameter(5.0);
        // p_random_force->SetDecayTimeParameter(20.0);
        simulator.AddForce(p_force);
        // simulator.AddForce(p_random_force);

        // Make the map from surface tensions to colors in surface evolver
        std::map<double, int> color_map;
        color_map[1.0] = 2;
        color_map[0.1] = 3;
        color_map[0.2] = 4;
        color_map[0.3] = 5;
        color_map[0.4] = 6;
        color_map[0.5] = 7;
        color_map[0.6] = 8;
        color_map[0.7] = 9;

        double avg_cell_volume = 1.0;

        // double relative_volume_range = avg_cell_volume * 0.5; // 50% volume variations possible

        MAKE_PTR(SurfaceEvolverSaveModifier<3>, p_surface_evolver_modifier);
        p_surface_evolver_modifier->SetSaveAtEnd();
        p_surface_evolver_modifier->SetMapTensionToColor(color_map);
        p_surface_evolver_modifier->SetSurfaceTensionSubForce(p_subforce);
        // p_surface_evolver_modifier->SetUseRandomizedVolumes(true, relative_volume_range);
        // p_surface_evolver_modifier->SetFixBoundaryNodes(true);
        simulator.AddSimulationModifier(p_surface_evolver_modifier);

        // MAKE_PTR(SimpleTargetVolumeModifier<3>, p_growth_modifier);
        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_growth_modifier, (&cell_population));
        p_growth_modifier->SetGrowthDuration(0.0);
        p_growth_modifier->SetT1AdaptationDuration(0.100);
        p_growth_modifier->SetReferenceTargetVolume(avg_cell_volume);
        // p_growth_modifier->SetT1AdaptationDuration(0.0);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Add simply supported boundary condition at edges
        c_vector<double, 3> mid_zero_point = zero_vector<double>(3);
        c_vector<double, 3> mid_normal_vector = zero_vector<double>(3);
        mid_normal_vector[2] = 1.0;
        mid_zero_point[2] = height / 2.0;
        MAKE_PTR_ARGS(SimplySupportedEdgeBoundaryCondition<3>,
                      p_supported_bc,
                      (&cell_population, mid_zero_point, mid_normal_vector));
        simulator.AddCellPopulationBoundaryCondition(p_supported_bc);

        simulator.Solve();

        simulator.SetEndTime(400.0);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetDt(0.1);

        simulator.Solve();

        simulator.SetEndTime(1000.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetDt(0.2);

        simulator.Solve();

        // Now introduce the force
        left_force[0] = 0.05;
        right_force[0] = 0.05;

        MAKE_PTR(RectangularEdgeForce<3>, p_force_edge);
        p_force_edge->SetLeftEdgeNodes(left_edge_nodes);
        p_force_edge->SetRightEdgeNodes(right_edge_nodes);
        p_force_edge->SetLeftEdgeForce(left_force);
        p_force_edge->SetRightEdgeForce(right_force);
        p_force_edge->UseInternalCoordinates();
        simulator.AddForce(p_force_edge);

        generator.CurveMonolayer(50.0);

        simulator.SetEndTime(4000.0);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetDt(0.1);

        simulator.Solve();

        simulator.SetEndTime(8000.0);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetDt(0.2);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 2000.0, 1e-10);
    }
};