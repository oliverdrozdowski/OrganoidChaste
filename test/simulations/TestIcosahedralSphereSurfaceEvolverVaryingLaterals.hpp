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
#include "FiniteThicknessIcosahedralSphereMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"
#include "MechanicalMosaicCellMutationState.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
#include "NodeTypeWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

#include "DecayingRandomForce.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SurfaceEvolverSaveModifier.hpp"
#include "SurfaceTensionForce.hpp"

#include "SmartPointers.hpp"
#include "SubstrateBoundaryCondition.hpp"

#include "FakePetscSetup.hpp"
// Test for memory leaks
#include "PetscSetupAndFinalize.hpp"

#include <cmath> //for M_PI

class TestIcosahedralSphereSurfaceEvolver : public AbstractCellBasedTestSuite
{
public:
    void TestAndVisualizeSphere()
    {
        // num_cells should be ~ pi*( 2^(1/3) / 3^(1/6) ) * (alpha-beta)^-2
        unsigned num_caspar_klug_h = 5;
        unsigned num_caspar_klug_k = 1;
        unsigned T = (num_caspar_klug_h * num_caspar_klug_h + num_caspar_klug_h * num_caspar_klug_k + num_caspar_klug_k * num_caspar_klug_k);
        unsigned num_cells = 10 * (T - 1) + 12;
        ;
        // alpha = apical/lateral surface tension
        double alpha = 0.75 / 1.0;
        // beta = basal/lateral surface tension
        double beta = 0.65 / 1.0;
        double height = cbrt(2.0 / sqrt(3.0)) * cbrt((alpha + beta) * (alpha + beta));
        double radius = 1.0 / (2.0 * cbrt(alpha + beta) * (alpha - beta)) * 2.0; // make larger to get rid of errors

        std::cout << "Height, radius:" << height << " , " << radius << std::flush;

        // length at which T1s become favorable: l^*=0.66/3^(2/3)/gamma^(1/3)
        double t1_length = 0.66 / cbrt(3 * 3 * (1.0 + 1.0) / 0.7);
        FiniteThicknessIcosahedralSphereMeshGenerator generator(num_caspar_klug_h, num_caspar_klug_k, 0.01, 0.001, height, 10.0);
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
        MAKE_PTR(MechanicalMosaicCellMutationState, p_state_mosaic_1);
        MAKE_PTR(MechanicalMosaicCellMutationState, p_state_mosaic_2);
        MAKE_PTR(MechanicalMosaicCellMutationState, p_state_mosaic_3);
        MAKE_PTR(MechanicalMosaicCellMutationState, p_state_mosaic_4);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        for (unsigned cell_index = 0; cell_index < num_cells; cell_index++)
        {
            NoCellCycleModel* p_model = new NoCellCycleModel;
            p_model->SetDimension(3);

            CellPtr p_cell;
            unsigned rand_int = RandomNumberGenerator::Instance()->randMod(5);

            if (rand_int == 0)
                p_cell = CellPtr(new Cell(p_state_wt, p_model));
            else if (rand_int == 1)
                p_cell = CellPtr(new Cell(p_state_mosaic_1, p_model));
            else if (rand_int == 2)
                p_cell = CellPtr(new Cell(p_state_mosaic_2, p_model));
            else if (rand_int == 3)
                p_cell = CellPtr(new Cell(p_state_mosaic_3, p_model));
            else if (rand_int == 4)
                p_cell = CellPtr(new Cell(p_state_mosaic_4, p_model));
            else // Should not happen
                p_cell = CellPtr(new Cell(p_state_mosaic_4, p_model));

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
        cell_population.AddPopulationWriter<NodeVelocityWriter>();
        cell_population.AddPopulationWriter<NodeTypeWriter>();
        cell_population.AddFaceWriter<FaceTypeWriter>();

        FiniteThicknessSimulation3d simulator(cell_population);
        simulator.SetOutputDirectory("Test3DSphere_icosahedral_" + std::to_string(num_caspar_klug_h) + "_"
                                     + std::to_string(num_caspar_klug_k) + "_a_" + std::to_string(alpha) + "_b_" + std::to_string(beta));
        simulator.SetEndTime(1.20);
        simulator.SetSamplingTimestepMultiple(5);
        simulator.SetDt(0.003);

        // Make a dictionary for mosaic cell surface tension parameters
        std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3> > dictionary_surface_tensions;
        std::array<double, 3> wild_type_surface_tensions = { 0.8, 0.65, 0.8 }; // apical, basal, lateral
        std::array<double, 3> mosaic_surface_tensions_1 = { 0.7, 0.6, 0.9 };
        std::array<double, 3> mosaic_surface_tensions_2 = { 0.8, 0.7, 1.0 };
        std::array<double, 3> mosaic_surface_tensions_3 = { 0.7, 0.6, 1.1 };
        std::array<double, 3> mosaic_surface_tensions_4 = { 0.8, 0.7, 1.2 };

        dictionary_surface_tensions[p_state_wt] = wild_type_surface_tensions;
        dictionary_surface_tensions[p_state_mosaic_1] = mosaic_surface_tensions_1;
        dictionary_surface_tensions[p_state_mosaic_2] = mosaic_surface_tensions_2;
        dictionary_surface_tensions[p_state_mosaic_3] = mosaic_surface_tensions_3;
        dictionary_surface_tensions[p_state_mosaic_4] = mosaic_surface_tensions_4;

        MAKE_PTR(SurfaceTensionForce<3>, p_force);
        // MAKE_PTR(DecayingRandomForce<3>, p_random_force);
        // p_force->CreateSurfaceTensionParametersForCells(1.0, 1.0, 0.7, p_mesh); // apical, basal, lateral
        p_force->SetSurfaceTensionParametersByMutation(dictionary_surface_tensions);
        p_force->SetSimulatedAnnealingParameters(0.0, 190.0, 10.0);
        p_force->SetSimulationInstance(&simulator);
        // p_force->SetPerformActiveT1Swaps();
        // p_random_force->SetRelativeStrengthParameter(5.0);
        // p_random_force->SetDecayTimeParameter(20.0);
        simulator.AddForce(p_force);
        // simulator.AddForce(p_random_force);

        // Make the map from surface tensions to colors in surface evolver
        std::map<double, int> color_map;
        color_map[0.6] = 1;
        color_map[0.65] = 2;
        color_map[0.7] = 3;
        color_map[0.8] = 4;
        color_map[0.85] = 5;
        color_map[0.9] = 6;
        color_map[0.95] = 7;
        color_map[1.0] = 8;
        color_map[1.05] = 9;
        color_map[1.1] = 10;
        color_map[1.15] = 11;
        color_map[1.2] = 12;

        double in_rad = radius - height / 2.0;
        double out_rad = radius + height / 2.0;
        double avg_cell_volume = 4.0 / 3.0 * M_PI * (out_rad * out_rad * out_rad - in_rad * in_rad * in_rad) / num_cells;
        // double avg_cell_volume = 1.0;

        double relative_volume_range = avg_cell_volume * 0.2; // 20% volume variations possible

        MAKE_PTR(SurfaceEvolverSaveModifier<3>, p_surface_evolver_modifier);
        p_surface_evolver_modifier->SetSaveAtEnd();
        p_surface_evolver_modifier->SetMapTensionToColor(color_map);
        p_surface_evolver_modifier->SetSurfaceTensionForce(p_force);
        p_surface_evolver_modifier->SetUseRandomizedVolumes(true, relative_volume_range);
        p_surface_evolver_modifier->SetWriteFaceTypeIntoFile(true);
        // p_surface_evolver_modifier->SetFixBoundaryNodes(true);
        simulator.AddSimulationModifier(p_surface_evolver_modifier);

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
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 1.2, 1e-10);
    }
};