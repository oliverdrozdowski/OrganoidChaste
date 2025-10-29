#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"

#include "CellProliferativeTypesWriter.hpp"
#include "CellThicknessWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellsGenerator.hpp"
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
#include "RandomSubForce.hpp"
#include "SimpleTargetVolumeModifier.hpp"
#include "SurfaceTensionSubForce.hpp"

#include "FiniteThicknessHoneycombVertexMeshGenerator.hpp"

#include "SmartPointers.hpp"

#include "FakePetscSetup.hpp"
// Test for memory leaks
#include "PetscSetupAndFinalize.hpp"

#include "UblasCustomFunctions.hpp" // For printing a vector

#include <cmath> //for M_PI
#include <ctime> // for current timestamp

class TestMonolayerVisualization : public AbstractCellBasedTestSuite
{
public:
    void Test3DVertexBasedSphereSimulationForVisualizing()
    {
        unsigned num_elem_up = 10;
        unsigned num_elem_across = 10;

        // double gamma_sum = 2.0;

        /*

        std::cout << "Gamma A   Gamma B   Energy" << std::endl;

        for (int i = 0; i < 10; i++)
        {

            double gamma_a = 0.5 + (double)i / 20.0;
            double gamma_b = gamma_sum - gamma_a;
            double height = cbrt(4.0 * 4.0) * cbrt(3.0 * sqrt(3.0) / 2.0) / cbrt(6.0 * 6.0) * cbrt((gamma_a + gamma_b) * (gamma_a + gamma_b));
            double element_area = 1 / height;

            FiniteThicknessHoneycombVertexMeshGenerator generator(num_elem_across, num_elem_up, false, 0.0, 0.001, height, element_area);
            MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

            std::vector<CellPtr> cells;
            CellsGenerator<NoCellCycleModel, 3> cells_generator;
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

            MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
            cell_population.SetOutputCellRearrangementLocations(false);
            cell_population.SetRestrictVertexMovementBoolean(false);

            cell_population.AddCellWriter<CellVolumesWriter>();
            cell_population.AddCellWriter<CellThicknessWriter>();
            cell_population.AddFaceWriter<FaceTypeWriter>();

            FiniteThicknessSimulation3d simulator(cell_population);

            simulator.SetOutputDirectory("TestWrapHoneycombOverSphere");
            simulator.SetEndTime(0.1 + (double)i / 10.0);
            simulator.SetSamplingTimestepMultiple(1);
            simulator.SetDt(0.1);

            simulator.Solve();

            MAKE_PTR(SurfaceTensionSubForce<3>, p_tension_subforce);

            p_tension_subforce->SetSimulationInstance(&simulator);
            p_tension_subforce->CreateSurfaceTensionParametersForCells(gamma_a, gamma_b, 1.0, p_mesh); // apical, basal, lateral

            std::cout << gamma_a << "   " << gamma_b << "   " << p_tension_subforce->CalculateTotalEnergy(cell_population) << std::endl;
        }
         */
        double gamma_a = 0.8;
        double gamma_b = 1.2;
        double height = cbrt(4.0 * 4.0) * cbrt(3.0 * sqrt(3.0) / 2.0) / cbrt(6.0 * 6.0) * cbrt((gamma_a + gamma_b) * (gamma_a + gamma_b));
        double element_area = 1 / height;

        FiniteThicknessHoneycombVertexMeshGenerator generator(num_elem_across, num_elem_up, false, 0.0, 0.001, height, element_area);

        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellRearrangementLocations(false);
        cell_population.SetRestrictVertexMovementBoolean(false);

        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellThicknessWriter>();
        cell_population.AddFaceWriter<FaceTypeWriter>();

        FiniteThicknessSimulation3d simulator(cell_population);

        // Create Forces
        MAKE_PTR(GeneralizedVolumeConservingForce<3>, p_force);
        MAKE_PTR(SurfaceTensionSubForce<3>, p_tension_subforce);
        p_force->AddSubForce(p_tension_subforce);
        p_force->SetSimulationInstance(&simulator);
        p_tension_subforce->SetSimulationInstance(&simulator);
        p_tension_subforce->CreateSurfaceTensionParametersForCells(gamma_a, gamma_b, 1.0, p_mesh); // apical, basal, lateral
        simulator.AddForce(p_force);

    

        simulator.SetOutputDirectory("Curling_plate");
        simulator.SetEndTime(5);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetDt(0.1);

        simulator.Solve();
    }
};