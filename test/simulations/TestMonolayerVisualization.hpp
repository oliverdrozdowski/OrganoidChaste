#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"

#include "CellProliferativeTypesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FaceTypeWriter.hpp"
#include "FiniteThicknessHoneycombVertexMeshGenerator.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

#include "FakePetscSetup.hpp"

class TestMonolayerVisualization : public AbstractCellBasedTestSuite
{
public:
    void Test2DVertexBasedMonolayerSimulationForVisualizing()
    {
        FiniteThicknessHoneycombVertexMeshGenerator generator(7, 9);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 63u);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        TS_ASSERT_EQUALS(cells.size(), 63u);

        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellRearrangementLocations(false);
        cell_population.SetRestrictVertexMovementBoolean(false);

        CellPtr p_different_cell = cell_population.GetCellUsingLocationIndex(31);
        p_different_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddFaceWriter<FaceTypeWriter>();

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("Test3DMonolayerVisualization");
        simulator.SetEndTime(0.1);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 63u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.1, 1e-10);
    }
};
