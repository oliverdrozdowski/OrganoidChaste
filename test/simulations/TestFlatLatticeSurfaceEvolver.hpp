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
#include "SimpleTargetVolumeModifier.hpp"
#include "SurfaceEvolverSaveModifier.hpp"
#include "SurfaceTensionSubForce.hpp"

#include "SmartPointers.hpp"
#include "SubstrateBoundaryCondition.hpp"

#include "FakePetscSetup.hpp"
// Test for memory leaks
#include "PetscSetupAndFinalize.hpp"

#include <cmath> //for M_PI

class TestFlatLatticeSurfaceEvolver : public AbstractCellBasedTestSuite
{
public:
    void Test3DFlatLattuceSurfaceEvolver()
    {
        unsigned num_cells = 81;
        double height = 2.0 / 3.0 / sqrt(3.0) * cbrt((9.0 / 2.0) * (9.0 / 2.0)) * cbrt(1.0 / 1.0 * 1.0 / 1.0) * 1.0;
        // double area = 1.0/height;
        //  length at which T1s become favorable: l^*=0.66/3^(2/3)/gamma^(1/3)
        // double t1_length = 0.66/cbrt(3*3 * (1.0+1.0)/0.7);
        FiniteThicknessHoneycombVertexMeshGenerator generator(9, 9, false, 0.01, 0.001, height, 1.0);
        // FiniteThicknessIcosahedralSphereMeshGenerator generator(11, 0, 0.01, 0.001, height, 10.0);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        // num_cells = p_mesh->GetNumElements();
        //  We use protorosettes as intermediate states
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
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        for (unsigned cell_index = 0; cell_index < num_cells; cell_index++)
        {
            c_vector<double, 3> centr = p_mesh->GetCentroidOfElement(cell_index);
            NoCellCycleModel* p_model = new NoCellCycleModel;
            p_model->SetDimension(3);

            CellPtr p_cell;

            // Randomly assign cell state (i.e. tensions)
            unsigned rand_int = RandomNumberGenerator::Instance()->randMod(4);
            if (rand_int == 0)
                p_cell = CellPtr(new Cell(p_state_wt, p_model));
            else if (rand_int == 1)
                p_cell = CellPtr(new Cell(p_state_mosaic_1, p_model));
            else if (rand_int == 2)
                p_cell = CellPtr(new Cell(p_state_mosaic_2, p_model));
            else
                p_cell = CellPtr(new Cell(p_state_mosaic_3, p_model));

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
        simulator.SetOutputDirectory("Test3D_Honeycomb_Surface_Evolver");
        simulator.SetEndTime(0.003);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetDt(0.003);

        // Make a dictionary for mosaic cell surface tension parameters
        std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3> > dictionary_surface_tensions;
        std::array<double, 3> wild_type_surface_tensions = { 1.0, 1.0, 0.7 }; // apical, basal, lateral
        std::array<double, 3> mosaic_surface_tensions_1 = { 1.0, 1.0, 0.5 };
        std::array<double, 3> mosaic_surface_tensions_2 = { 1.0, 1.0, 0.3 };
        std::array<double, 3> mosaic_surface_tensions_3 = { 1.0, 1.0, 0.1 };

        dictionary_surface_tensions[p_state_wt] = wild_type_surface_tensions;
        dictionary_surface_tensions[p_state_mosaic_1] = mosaic_surface_tensions_1;
        dictionary_surface_tensions[p_state_mosaic_2] = mosaic_surface_tensions_2;
        dictionary_surface_tensions[p_state_mosaic_3] = mosaic_surface_tensions_3;

        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_subforce);
        p_force->AddSubForce(p_subforce);
        p_force->SetSimulationInstance(&simulator);
        // MAKE_PTR(DecayingRandomForce<3>, p_random_force);
        // p_force->CreateSurfaceTensionParametersForCells(1.0, 1.0, 0.7, p_mesh); // apical, basal, lateral
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

        double avg_cell_volume = height;

        double relative_volume_range = avg_cell_volume * 0.5; // 50% volume variations possible

        MAKE_PTR(SurfaceEvolverSaveModifier<3>, p_surface_evolver_modifier);
        p_surface_evolver_modifier->SetSaveAtEnd();
        p_surface_evolver_modifier->SetMapTensionToColor(color_map);
        p_surface_evolver_modifier->SetSurfaceTensionSubForce(p_subforce);
        p_surface_evolver_modifier->SetUseRandomizedVolumes(true, relative_volume_range);
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

        // Set the boundary surface tensions to match apical/basal
        p_subforce->UpdateSurfaceTensionsByMutation(&cell_population);
        unsigned num_faces = p_mesh->GetNumFaces();
        for (unsigned index_face = 0; index_face < num_faces; index_face++)
        {
            if (p_mesh->GetFace(index_face)->IsBoundaryFace())
            {
                p_subforce->SetSurfaceTensionParameter(index_face, 0.5); // /2 !
            }
        }

        p_surface_evolver_modifier->UpdateAtEndOfSolve(cell_population);

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.003, 1e-10);
    }
};