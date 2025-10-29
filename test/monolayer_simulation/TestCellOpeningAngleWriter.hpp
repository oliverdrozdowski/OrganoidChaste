#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"

#include "CellOpeningAngleWriter.hpp"
#include "CellsGenerator.hpp"
#include "FiniteThicknessHoneycombVertexMeshGenerator.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MonolayerVertexElement.hpp"
#include "MonolayerVertexMesh.hpp"
#include "MutableMonolayerVertexMesh.hpp"
#include "NoCellCycleModel.hpp"

#include "FakePetscSetup.hpp"
// Test for memory leaks
#include "PetscSetupAndFinalize.hpp"

#include <cmath> //for M_PI

class TestWriteToSurfaceEvolver : public AbstractCellBasedTestSuite
{
private:
    MutableMonolayerVertexMesh<3, 3>* ConstructAngledCubesMesh()
    {
        // Make 8 nodes to assign to a cube with opening upwards
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, -0.5, -0.5, 0.0));
        nodes.push_back(new Node<3>(1, false, 0.5, -0.5, 0.0));
        nodes.push_back(new Node<3>(2, false, -0.5, 0.5, 0.0));
        nodes.push_back(new Node<3>(3, false, -1.5, -1.5, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.0));
        nodes.push_back(new Node<3>(5, false, -1.5, 1.5, 1.0));
        nodes.push_back(new Node<3>(6, false, 1.5, -1.5, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.5, 1.5, 1.0));

        // Set basal/apical nodetypes
        std::vector<MonolayerVertexElementType> node_types;
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Basal);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);
        node_types.push_back(MonolayerVertexElementType::Apical);

        // Make six square faces and four triangular faces out of these nodes
        std::vector<std::vector<Node<3>*> > nodes_faces(10);
        std::vector<std::vector<MonolayerVertexElementType> > nodes_faces_types(10);

        unsigned node_indices_face_0[4] = { 0, 2, 4, 1 };
        unsigned node_indices_face_1[4] = { 4, 7, 5, 2 };
        unsigned node_indices_face_2[4] = { 7, 6, 1, 4 };
        unsigned node_indices_face_3[4] = { 0, 3, 5, 2 };
        unsigned node_indices_face_4[4] = { 1, 6, 3, 0 };
        unsigned node_indices_face_5[4] = { 7, 6, 3, 5 };

        for (unsigned i = 0; i < 4; i++)
        {
            nodes_faces[0].push_back(nodes[node_indices_face_0[i]]);
            nodes_faces_types[0].push_back(node_types[node_indices_face_0[i]]);

            nodes_faces[1].push_back(nodes[node_indices_face_1[i]]);
            nodes_faces_types[1].push_back(node_types[node_indices_face_1[i]]);

            nodes_faces[2].push_back(nodes[node_indices_face_2[i]]);
            nodes_faces_types[2].push_back(node_types[node_indices_face_2[i]]);

            nodes_faces[3].push_back(nodes[node_indices_face_3[i]]);
            nodes_faces_types[3].push_back(node_types[node_indices_face_3[i]]);

            nodes_faces[4].push_back(nodes[node_indices_face_4[i]]);
            nodes_faces_types[4].push_back(node_types[node_indices_face_4[i]]);

            nodes_faces[5].push_back(nodes[node_indices_face_5[i]]);
            nodes_faces_types[5].push_back(node_types[node_indices_face_5[i]]);

            nodes_faces[6].push_back(nodes[node_indices_face_0[i]]);
            nodes_faces_types[6].push_back(node_types[node_indices_face_5[i]]);

            nodes_faces[7].push_back(nodes[node_indices_face_5[i]]);
            nodes_faces_types[7].push_back(node_types[node_indices_face_0[i]]);
        }

        // Make the faces
        std::vector<MonolayerVertexElement<2, 3>*> faces;

        // Set face types
        std::vector<MonolayerVertexElementType> face_types;
        face_types.push_back(MonolayerVertexElementType::Basal);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Lateral);
        face_types.push_back(MonolayerVertexElementType::Apical);

        // Generate opposite apical and basal faces
        face_types.push_back(MonolayerVertexElementType::Apical);
        face_types.push_back(MonolayerVertexElementType::Basal);

        for (unsigned i = 0; i < 8; i++)
        {
            faces.push_back(new MonolayerVertexElement<2, 3>(i, face_types[i], nodes_faces[i], nodes_faces_types[i]));
        }

        // Make the elements
        std::vector<MonolayerVertexElement<2, 3>*> faces_element_0, faces_element_1;
        std::vector<bool> orientations_0, orientations_1;

        // Cube element
        faces_element_1.push_back(faces[6]);
        orientations_1.push_back(true);
        for (unsigned i = 0; i < 6; i++)
        {
            faces_element_0.push_back(faces[i]);
            orientations_0.push_back(true);

            if (i != 0 || i != 5)
            {
                faces_element_1.push_back(faces[i]);
                orientations_1.push_back(true);
            }
        }
        faces_element_1.push_back(faces[7]);
        orientations_1.push_back(true);

        std::vector<MonolayerVertexElement<3, 3>*> elements;
        elements.push_back(new MonolayerVertexElement<3, 3>(0, MonolayerVertexElementType::Undetermined, faces_element_0, orientations_0));
        elements.push_back(new MonolayerVertexElement<3, 3>(1, MonolayerVertexElementType::Undetermined, faces_element_1, orientations_1));

        return new MutableMonolayerVertexMesh<3, 3>(nodes, faces, elements);
    }

public:
    void TestAnglesForPyramidBase()
    {
        MutableMonolayerVertexMesh<3, 3>* p_mesh = ConstructAngledCubesMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        CellOpeningAngleWriter<3, 3> writer;

        double opening_angle = writer.GetCellDataForVtkOutput(*(cell_population.rGetCells().begin()), &cell_population);
        TS_ASSERT_DELTA(opening_angle, M_PI / 2.0, 1e-10);

        // Now switched apical and basal sides should yield the same opening angle (we do not consider direction)
        double opening_angle_2 = writer.GetCellDataForVtkOutput(*(cell_population.rGetCells().begin()++), &cell_population);
        TS_ASSERT_DELTA(opening_angle_2, M_PI / 2.0, 1e-10);
    }

    void TestAngleForHexagon()
    {
        FiniteThicknessHoneycombVertexMeshGenerator generator(1, 1);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        CellOpeningAngleWriter<3, 3> writer;

        double opening_angle = writer.GetCellDataForVtkOutput(*(cell_population.rGetCells().begin()), &cell_population);
        TS_ASSERT_DELTA(opening_angle, 0.0, 1e-10);
    }
};