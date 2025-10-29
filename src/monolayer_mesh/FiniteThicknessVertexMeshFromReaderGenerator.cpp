#include "FiniteThicknessVertexMeshFromReaderGenerator.hpp"

#include <boost/numeric/ublas/assignment.hpp>
#include <cmath> //for M_PI
#include "UblasCustomFunctions.hpp"

#include "MonolayerVertexElement.hpp"

FiniteThicknessVertexMeshFromReaderGenerator::
    FiniteThicknessVertexMeshFromReaderGenerator(MonolayerVertexMeshReader<3, 3>& mesh_reader, double cellRearrangementThreshold,
                                                 double t2Threshold)
{
    assert(cellRearrangementThreshold >= 0.0);
    assert(t2Threshold > 0.0);

    // Store numbers of nodes and elements
    unsigned num_nodes = mesh_reader.GetNumNodes();
    unsigned num_faces = mesh_reader.GetNumFaces();
    unsigned num_elements = mesh_reader.GetNumElements();

    mesh_reader.Reset();

    std::vector<Node<3>*> nodes;

    // Add nodes
    std::vector<double> node_data;
    for (unsigned i = 0; i < num_nodes; i++)
    {
        node_data = mesh_reader.GetNextNode();

        // remove last element of vector (which is isBoundaryNode and we don't use it)
        node_data.pop_back();

        nodes.push_back(new Node<3>(i, false, node_data[0], node_data[1], node_data[2]));
    }

    std::vector<MonolayerVertexElement<2, 3>*> faces;

    for (unsigned face_index = 0; face_index < num_faces; face_index++)
    {
        MonolayerVertexFaceData face_data;
        face_data = mesh_reader.GetNextMonolayerFaceData();

        unsigned num_nodes_in_face = face_data.NodeIndices.size();
        std::vector<Node<3>*> nodes_in_face;
        std::vector<MonolayerVertexElementType> node_types;

        for (unsigned node_index = 0; node_index < num_nodes_in_face; node_index++)
        {
            nodes_in_face.push_back(nodes[face_data.NodeIndices[node_index]]);
            node_types.push_back(static_cast<MonolayerVertexElementType>(face_data.NodeTypes[node_index]));
        }

        MonolayerVertexElement<2, 3>* face = new MonolayerVertexElement<2, 3>(face_index, static_cast<MonolayerVertexElementType>(face_data.FaceType), nodes_in_face, node_types);

        if (face_data.FaceType == 1)
        {
            c_vector<double, 3> center = GetCenterOfFace(face);
            c_vector<double, 3> norml = GetNormalToFace(face);

            if (inner_prod(center, norml) > 0)
                std::cout << "Face " << face_index << " is oriented outward!!" << std::endl;
        }

        faces.push_back(face);
    }

    std::vector<MonolayerVertexElement<3, 3>*> elements;

    // Add elements
    for (unsigned elem_index = 0; elem_index < num_elements; elem_index++)
    {
        MonolayerVertexElementData element_data;
        element_data = mesh_reader.GetNextMonolayerElementData();

        unsigned num_faces_in_element = element_data.FaceIndices.size();
        std::vector<MonolayerVertexElement<2, 3>*> faces_in_element;
        std::vector<bool> orientations;

        for (unsigned face_index = 0; face_index < num_faces_in_element; face_index++)
        {
            faces_in_element.push_back(faces[element_data.FaceIndices[face_index]]);
            orientations.push_back(element_data.FaceOrientations[face_index]);
        }

        MonolayerVertexElement<3, 3>* element = new MonolayerVertexElement<3, 3>(elem_index, MonolayerVertexElementType::Undetermined, faces_in_element, orientations);

        elements.push_back(element);
    }

    mpMesh = new MutableMonolayerVertexMesh<3, 3>(
        nodes, faces, elements, cellRearrangementThreshold, t2Threshold);
}

FiniteThicknessVertexMeshFromReaderGenerator::
    ~FiniteThicknessVertexMeshFromReaderGenerator()
{
    delete mpMesh;
}

MutableMonolayerVertexMesh<3, 3>* FiniteThicknessVertexMeshFromReaderGenerator::GetMesh()
{
    return mpMesh;
}

c_vector<double, 3> FiniteThicknessVertexMeshFromReaderGenerator::GetNormalToFace(MonolayerVertexElement<2, 3>* pFace)
{
    // As we are in 3D, the face must have at least three vertices
    assert(pFace->GetNumNodes() >= 3u);

    // calculate the center
    // unsigned num_nodes = pFace->GetNumNodes();
    c_vector<double, 3> center = GetCenterOfFace(pFace);

    // Calculate the normal as mean value of triangle normals
    c_vector<double, 3> normal_vector = zero_vector<double>(3);

    c_vector<double, 3> v_minus_v_pc = pFace->GetNode(0)->rGetLocation() - center;
    for (unsigned local_index = 0; local_index < pFace->GetNumNodes();
         local_index++)
    {
        c_vector<double, 3> vnext_minus_v_pc = pFace->GetNode((local_index + 1) % pFace->GetNumNodes())
                                                   ->rGetLocation()
            - center;

        c_vector<double, 3> localNormal = VectorProduct(v_minus_v_pc, vnext_minus_v_pc);

        // normal direction is implemented as mean normal of triangles (note
        // possible errors from orientation)
        normal_vector += localNormal;
        v_minus_v_pc = vnext_minus_v_pc;
    }
    double magnitude = norm_2(normal_vector);
    if (magnitude != 0.0)
    {
        // Normalize the normal vector
        normal_vector /= magnitude;
        // If all points are co-located, then magnitude==0.0 and there is potential
        // for a floating point exception here if we divide by zero, so we'll move
        // on.
    }
    return normal_vector;
}

c_vector<double, 3>
FiniteThicknessVertexMeshFromReaderGenerator::GetCenterOfFace(
    MonolayerVertexElement<2, 3>* pFace)
{
    // As we are in 3D, the face must have at least three vertices
    assert(pFace->GetNumNodes() >= 3u);

    // calculate the center
    unsigned num_nodes = pFace->GetNumNodes();
    c_vector<double, 3> center = zero_vector<double>(3);
    for (unsigned ind = 0; ind < num_nodes; ++ind)
    {
        center += pFace->GetNode(ind)->rGetLocation();
    }
    return center / (double)num_nodes;
}