#include "FiniteThicknessHoneycombVertexMeshGenerator.hpp"
#include <cmath>
#include "MonolayerVertexElement.hpp"
#include "UblasCustomFunctions.hpp"

FiniteThicknessHoneycombVertexMeshGenerator::FiniteThicknessHoneycombVertexMeshGenerator(unsigned numElementsAcross,
                                                                                         unsigned numElementsUp,
                                                                                         bool isFlatBottom,
                                                                                         double cellRearrangementThreshold,
                                                                                         double t2Threshold,
                                                                                         double height,
                                                                                         double elementArea)
        : mNumElementsAcross(numElementsAcross),
          mNumElementsUp(numElementsUp),
          mElementArea(elementArea),
          mHeight(height)
{
    assert(numElementsAcross > 0);
    assert(numElementsUp > 0);
    assert(cellRearrangementThreshold >= 0.0);
    assert(t2Threshold > 0.0);
    assert(height > 0.0);
    assert(elementArea > 0.0);

    std::vector<Node<3>*> nodes_apical;
    std::vector<Node<3>*> nodes_basal;
    std::vector<MonolayerVertexElement<2, 3>*> faces_apical;
    std::vector<MonolayerVertexElement<2, 3>*> faces_basal;
    std::vector<MonolayerVertexElement<2, 3>*> faces_lateral;
    std::vector<MonolayerVertexElement<3, 3>*> elements;

    // We create vectors with the edge nodes, if we want to retrieve them later
    std::vector<Node<3>*> nodes_left;
    std::vector<Node<3>*> nodes_top;
    std::vector<Node<3>*> nodes_right;
    std::vector<Node<3>*> nodes_bottom;

    MonolayerVertexElement<2, 3>* p_lateral_faces[6];

    unsigned node_index = 0;
    unsigned node_indices[6];
    unsigned element_index;

    /*
    // Create the apical side first at z=height
    */

    // Create the nodes, row by row, from the bottom up, start with the upper (apical) side

    // On the first row we have numElementsAcross nodes, all of which are boundary nodes
    for (unsigned i = 0; i < numElementsAcross; i++)
    {
        Node<3>* p_node = new Node<3>(node_index, true, i + 0.5, 0, height);
        nodes_apical.push_back(p_node);
        node_index++;
        nodes_bottom.push_back(p_node);
    }

    /*
     * On each interior row we have numElementsAcross+1 nodes. On the second and penultimate
     * row all nodes are boundary nodes. On other rows the first and last nodes only
     * are boundary nodes.
     */
    for (unsigned j = 1; j < 2 * numElementsUp + 1; j++)
    {
        for (unsigned i = 0; i <= numElementsAcross; i++)
        {
            double x_coord = ((j % 4 == 0) || (j % 4 == 3)) ? i + 0.5 : i;
            double y_coord = (1.5 * j - 0.5 * (j % 2)) * 0.5 / sqrt(3.0);
            bool is_boundary_node = (j == 1 || j == 2 * numElementsUp || i == 0 || i == numElementsAcross) ? true : false;

            Node<3>* p_node = new Node<3>(node_index, is_boundary_node, x_coord, y_coord, height);
            nodes_apical.push_back(p_node);
            node_index++;

            if (is_boundary_node)
            {
                if (j == 1)
                {
                    if (i > 0)
                    {
                        nodes_bottom.push_back(p_node);
                    }
                    if (i == 0)
                    {
                        nodes_left.push_back(p_node);
                    }
                    else if (i == numElementsAcross)
                    {
                        nodes_right.push_back(p_node);
                    }
                }
                else if (j == 2 * numElementsUp)
                {
                    if (numElementsUp % 2 == 0)
                    {
                        if (i != numElementsAcross)
                        {
                            nodes_top.push_back(p_node);
                        }
                    }
                    else
                    {
                        if (i != 0)
                        {
                            nodes_top.push_back(p_node);
                        }
                    }
                    if (i == 0)
                    {
                        nodes_left.push_back(p_node);
                    }
                    else if (i == numElementsAcross)
                    {
                        nodes_right.push_back(p_node);
                    }
                }
                else if (i == 0)
                {
                    nodes_left.push_back(p_node);
                }
                else if (i == numElementsAcross)
                {
                    nodes_right.push_back(p_node);
                }
            }
        }
    }

    /*
     * On the last row we have numElementsAcross nodes, all of which are boundary nodes.
     */
    double y_coord = (1.5 * (2 * numElementsUp + 1) - 0.5 * ((2 * numElementsUp + 1) % 2)) * 0.5 / sqrt(3.0);
    if (((2 * numElementsUp + 1) % 4 == 0) || ((2 * numElementsUp + 1) % 4 == 3))
    {
        Node<3>* p_node = new Node<3>(node_index, true, 0.5, y_coord, height);
        nodes_apical.push_back(p_node);
        node_index++;
        nodes_top.push_back(p_node);
    }
    for (unsigned i = 1; i < numElementsAcross; i++)
    {
        double x_coord = (((2 * numElementsUp + 1) % 4 == 0) || ((2 * numElementsUp + 1) % 4 == 3)) ? i + 0.5 : i;

        Node<3>* p_node = new Node<3>(node_index, true, x_coord, y_coord, height);
        nodes_apical.push_back(p_node);
        node_index++;
        nodes_top.push_back(p_node);
    }
    if (((2 * numElementsUp + 1) % 4 == 1) || ((2 * numElementsUp + 1) % 4 == 2))
    {
        Node<3>* p_node = new Node<3>(node_index, true, numElementsAcross, y_coord, height);
        nodes_apical.push_back(p_node);
        node_index++;
        nodes_top.push_back(p_node);
    }

    /*
     * Create the apical elements. The array node_indices contains the
     * global node indices from bottom, going anticlockwise, viewed from above.
     */
    for (unsigned j = 0; j < numElementsUp; j++)
    {
        for (unsigned i = 0; i < numElementsAcross; i++)
        {
            // Create apical faces
            if (j == 0)
            {
                node_indices[0] = i;
            }
            else
            {
                node_indices[0] = 2 * j * (numElementsAcross + 1) - 1 * (j % 2 == 0) + i; // different for even/odd rows
            }
            node_indices[1] = node_indices[0] + numElementsAcross + 1 + 1 * (j % 2 == 0 && j > 0);
            node_indices[2] = node_indices[1] + numElementsAcross + 1;
            node_indices[3] = node_indices[2] + numElementsAcross + 1 * (j % 2 == 1 && j < numElementsUp - 1);
            node_indices[4] = node_indices[2] - 1;
            node_indices[5] = node_indices[1] - 1;

            std::vector<Node<3>*> face_nodes;
            std::vector<MonolayerVertexElementType> face_node_types;
            for (unsigned k = 0; k < 6; k++)
            {
                face_nodes.push_back(nodes_apical[node_indices[k]]);
                face_node_types.push_back(MonolayerVertexElementType::Apical);
            }

            element_index = j * numElementsAcross + i;
            MonolayerVertexElement<2, 3>* p_face_apical = new MonolayerVertexElement<2, 3>(element_index, MonolayerVertexElementType::Apical, face_nodes, face_node_types);
            faces_apical.push_back(p_face_apical);
        }
    }

    /*
    // Create the basal and lateral sides second with basal side at z=0.0
    */

    // Create the nodes, row by row, from the bottom up, now the basal side

    // On the first row we have numElementsAcross nodes, all of which are boundary nodes
    for (unsigned i = 0; i < numElementsAcross; i++)
    {
        Node<3>* p_node = new Node<3>(node_index, true, i + 0.5, 0, 0.0);
        nodes_basal.push_back(p_node);
        node_index++;
        nodes_bottom.push_back(p_node);
    }

    /*
     * On each interior row we have numElementsAcross+1 nodes. On the second and penultimate
     * row all nodes are boundary nodes. On other rows the first and last nodes only
     * are boundary nodes.
     */
    for (unsigned j = 1; j < 2 * numElementsUp + 1; j++)
    {
        for (unsigned i = 0; i <= numElementsAcross; i++)
        {
            double x_coord = ((j % 4 == 0) || (j % 4 == 3)) ? i + 0.5 : i;
            double y_coord = (1.5 * j - 0.5 * (j % 2)) * 0.5 / sqrt(3.0);
            bool is_boundary_node = (j == 1 || j == 2 * numElementsUp || i == 0 || i == numElementsAcross) ? true : false;

            Node<3>* p_node = new Node<3>(node_index, is_boundary_node, x_coord, y_coord, 0.0);
            nodes_basal.push_back(p_node);
            node_index++;

            if (is_boundary_node)
            {
                if (j == 1)
                {
                    if (i > 0)
                    {
                        nodes_bottom.push_back(p_node);
                    }
                    if (i == 0)
                    {
                        nodes_left.push_back(p_node);
                    }
                    else if (i == numElementsAcross)
                    {
                        nodes_right.push_back(p_node);
                    }
                }
                else if (j == 2 * numElementsUp)
                {
                    if (numElementsUp % 2 == 0)
                    {
                        if (i != numElementsAcross)
                        {
                            nodes_top.push_back(p_node);
                        }
                    }
                    else
                    {
                        if (i != 0)
                        {
                            nodes_top.push_back(p_node);
                        }
                    }
                    if (i == 0)
                    {
                        nodes_left.push_back(p_node);
                    }
                    else if (i == numElementsAcross)
                    {
                        nodes_right.push_back(p_node);
                    }
                }
                else if (i == 0)
                {
                    nodes_left.push_back(p_node);
                }
                else
                {
                    nodes_right.push_back(p_node);
                }
            }
        }
    }

    /*
     * On the last row we have numElementsAcross nodes, all of which are boundary nodes.
     */
    y_coord = (1.5 * (2 * numElementsUp + 1) - 0.5 * ((2 * numElementsUp + 1) % 2)) * 0.5 / sqrt(3.0);
    if (((2 * numElementsUp + 1) % 4 == 0) || ((2 * numElementsUp + 1) % 4 == 3))
    {
        Node<3>* p_node = new Node<3>(node_index, true, 0.5, y_coord, 0.0);
        nodes_basal.push_back(p_node);
        node_index++;
        nodes_top.push_back(p_node);
    }
    for (unsigned i = 1; i < numElementsAcross; i++)
    {
        double x_coord = (((2 * numElementsUp + 1) % 4 == 0) || ((2 * numElementsUp + 1) % 4 == 3)) ? i + 0.5 : i;

        Node<3>* p_node = new Node<3>(node_index, true, x_coord, y_coord, 0.0);
        nodes_basal.push_back(p_node);
        node_index++;
        nodes_top.push_back(p_node);
    }
    if (((2 * numElementsUp + 1) % 4 == 1) || ((2 * numElementsUp + 1) % 4 == 2))
    {
        Node<3>* p_node = new Node<3>(node_index, true, numElementsAcross, y_coord, 0.0);
        nodes_basal.push_back(p_node);
        node_index++;
        nodes_top.push_back(p_node);
    }

    /*
     * Create the basal elements. The array node_indices contains the
     * global node indices from bottom, going clockwise, viewed from above.
     */
    for (unsigned j = 0; j < numElementsUp; j++)
    {
        for (unsigned i = 0; i < numElementsAcross; i++)
        {
            // Create basal faces
            if (j == 0)
            {
                node_indices[0] = i;
            }
            else
            {
                node_indices[0] = 2 * j * (numElementsAcross + 1) - 1 * (j % 2 == 0) + i; // different for even/odd rows
            }
            node_indices[5] = node_indices[0] + numElementsAcross + 1 + 1 * (j % 2 == 0 && j > 0);
            node_indices[4] = node_indices[5] + numElementsAcross + 1;
            node_indices[3] = node_indices[4] + numElementsAcross + 1 * (j % 2 == 1 && j < numElementsUp - 1);
            node_indices[2] = node_indices[4] - 1;
            node_indices[1] = node_indices[5] - 1;

            std::vector<Node<3>*> face_nodes;
            std::vector<MonolayerVertexElementType> face_node_types;
            for (unsigned k = 0; k < 6; k++)
            {
                face_nodes.push_back(nodes_basal[node_indices[k]]);
                face_node_types.push_back(MonolayerVertexElementType::Basal);
            }
            // number first apical faces, then basal
            element_index = numElementsAcross * numElementsUp + j * numElementsAcross + i;
            MonolayerVertexElement<2, 3>* p_face_basal = new MonolayerVertexElement<2, 3>(element_index, MonolayerVertexElementType::Basal, face_nodes, face_node_types);
            faces_basal.push_back(p_face_basal);

            // Create lateral faces and element

            /* Face index combinations for lateral faces anti-clockwise oriented, viewed from the outside
            // we start with the bottom-left (basal) node
            // The lateral faces are numbered starting from the bottom going anti-clockwise along the contour (viewed from above)
            //        3   d  2                iv------iii  apical z=height
            //           /\                    | face |
            //        c /  \ e                 |  1   |
            //         |    |                  i------ii  basal z=0
            //     4   |    |   1       face one seen from the right side
            //        b \  / f       (outside) with anti-clockwise numbering
            //           \/
            //       5   a   0
            // node numbering given by a-f (corresponding to 0 to 5) and face numbering by 0-5
            */
            std::vector<std::vector<Node<3>*> > nodes_lateral(6);

            nodes_lateral[0].push_back(nodes_basal[node_indices[0]]);
            nodes_lateral[0].push_back(nodes_basal[node_indices[5]]);
            nodes_lateral[0].push_back(nodes_apical[node_indices[5]]);
            nodes_lateral[0].push_back(nodes_apical[node_indices[0]]);

            nodes_lateral[1].push_back(nodes_basal[node_indices[5]]);
            nodes_lateral[1].push_back(nodes_basal[node_indices[4]]);
            nodes_lateral[1].push_back(nodes_apical[node_indices[4]]);
            nodes_lateral[1].push_back(nodes_apical[node_indices[5]]);

            nodes_lateral[2].push_back(nodes_basal[node_indices[4]]);
            nodes_lateral[2].push_back(nodes_basal[node_indices[3]]);
            nodes_lateral[2].push_back(nodes_apical[node_indices[3]]);
            nodes_lateral[2].push_back(nodes_apical[node_indices[4]]);

            nodes_lateral[3].push_back(nodes_basal[node_indices[3]]);
            nodes_lateral[3].push_back(nodes_basal[node_indices[2]]);
            nodes_lateral[3].push_back(nodes_apical[node_indices[2]]);
            nodes_lateral[3].push_back(nodes_apical[node_indices[3]]);

            nodes_lateral[4].push_back(nodes_basal[node_indices[2]]);
            nodes_lateral[4].push_back(nodes_basal[node_indices[1]]);
            nodes_lateral[4].push_back(nodes_apical[node_indices[1]]);
            nodes_lateral[4].push_back(nodes_apical[node_indices[2]]);

            nodes_lateral[5].push_back(nodes_basal[node_indices[1]]);
            nodes_lateral[5].push_back(nodes_basal[node_indices[0]]);
            nodes_lateral[5].push_back(nodes_apical[node_indices[0]]);
            nodes_lateral[5].push_back(nodes_apical[node_indices[1]]);

            std::vector<MonolayerVertexElementType> lateral_node_types;
            lateral_node_types.push_back(MonolayerVertexElementType::Basal);
            lateral_node_types.push_back(MonolayerVertexElementType::Basal);
            lateral_node_types.push_back(MonolayerVertexElementType::Apical);
            lateral_node_types.push_back(MonolayerVertexElementType::Apical);

            std::vector<bool> face_orientations(6);

            // First row we cannot use faces from previous rows
            if (j == 0)
            {
                unsigned element_index_base = 2 * numElementsAcross * numElementsUp;
                element_index_base += (i == 0 ? 0 : (6 + (i - 1) * 5)); // previous elements in this row

                face_orientations = std::vector<bool>(6, false);

                // For first element create all faces, including left
                if (i == 0)
                {
                    p_lateral_faces[4] = new MonolayerVertexElement<2, 3>(element_index_base + 4, MonolayerVertexElementType::Lateral, nodes_lateral[4], lateral_node_types);
                }
                // Other elements in first row share the left lateral face (#4) with the previous right face (#1)
                else
                {
                    // Copy old pointer
                    p_lateral_faces[4] = p_lateral_faces[1];
                    face_orientations[4] = true;
                }

                p_lateral_faces[0] = new MonolayerVertexElement<2, 3>(element_index_base, MonolayerVertexElementType::Lateral, nodes_lateral[0], lateral_node_types);
                p_lateral_faces[1] = new MonolayerVertexElement<2, 3>(element_index_base + 1, MonolayerVertexElementType::Lateral, nodes_lateral[1], lateral_node_types);
                p_lateral_faces[2] = new MonolayerVertexElement<2, 3>(element_index_base + 2, MonolayerVertexElementType::Lateral, nodes_lateral[2], lateral_node_types);
                p_lateral_faces[3] = new MonolayerVertexElement<2, 3>(element_index_base + 3, MonolayerVertexElementType::Lateral, nodes_lateral[3], lateral_node_types);

                unsigned additionalFactorIfFithFace = (i == 0) ? 1 : 0;
                p_lateral_faces[5] = new MonolayerVertexElement<2, 3>(element_index_base + additionalFactorIfFithFace + 4, MonolayerVertexElementType::Lateral, nodes_lateral[5], lateral_node_types);
            }

            // Elements in uneven rows share lower edges (#0, #5), except for last edge
            //     | X | X | X | uneven
            //    / \ / \ / \ /  <==
            //   | X | X | X |   even
            else if (j % 2 == 1)
            {
                unsigned element_index_base = 2 * numElementsAcross * numElementsUp + 6 + 5 * (numElementsAcross - 1); // until frist lateral row
                element_index_base += (j - 1) * (3 * numElementsAcross + 2); // all previous lateral rows
                element_index_base += (i == 0 ? 0 : (4 + (i - 1) * 3)); // previous elements in this row

                face_orientations = std::vector<bool>(6, false);

                // First element in row does not share left edge, but two bottom edges
                if (i == 0)
                {
                    // Copy old pointer
                    p_lateral_faces[0] = elements[(j - 1) * numElementsAcross + 1]->GetFace(5);
                    p_lateral_faces[5] = elements[(j - 1) * numElementsAcross]->GetFace(4);
                    face_orientations[0] = true;
                    face_orientations[5] = true;
                    // New left edge face
                    p_lateral_faces[4] = new MonolayerVertexElement<2, 3>(element_index_base + 3, MonolayerVertexElementType::Lateral, nodes_lateral[4], lateral_node_types);
                }
                // Last element shares left edge and only one bottom edge
                else if (i == numElementsAcross - 1)
                {
                    // Copy old pointer
                    p_lateral_faces[5] = elements[(j - 1) * numElementsAcross + i]->GetFace(4);
                    p_lateral_faces[4] = p_lateral_faces[1]; // left edge
                    face_orientations[5] = true;
                    face_orientations[4] = true;
                    // New bottom edge face
                    p_lateral_faces[0] = new MonolayerVertexElement<2, 3>(element_index_base, MonolayerVertexElementType::Lateral, nodes_lateral[0], lateral_node_types);
                }
                // Inner elements share left and bottom edges
                else
                {
                    // Copy old pointer
                    p_lateral_faces[0] = elements[(j - 1) * numElementsAcross + i + 1]->GetFace(5);
                    p_lateral_faces[5] = elements[(j - 1) * numElementsAcross + i]->GetFace(4);
                    p_lateral_faces[4] = p_lateral_faces[1]; // left edge
                    face_orientations[0] = true;
                    face_orientations[5] = true;
                    face_orientations[4] = true;
                }
                unsigned additionalFactorIfZerothFace = (i == numElementsAcross - 1) ? 1 : 0;
                // Create new faces
                p_lateral_faces[1] = new MonolayerVertexElement<2, 3>(element_index_base + additionalFactorIfZerothFace, MonolayerVertexElementType::Lateral, nodes_lateral[1], lateral_node_types);
                p_lateral_faces[2] = new MonolayerVertexElement<2, 3>(element_index_base + additionalFactorIfZerothFace + 1, MonolayerVertexElementType::Lateral, nodes_lateral[2], lateral_node_types);
                p_lateral_faces[3] = new MonolayerVertexElement<2, 3>(element_index_base + additionalFactorIfZerothFace + 2, MonolayerVertexElementType::Lateral, nodes_lateral[3], lateral_node_types);
            }

            // Elements in even rows share lower edges (#0, #5), except for first edge
            //   | X | X | X |   even
            //    \ / \ / \ / \  <==
            //     | X | X | X | uneven
            else if (j % 2 == 0)
            {
                unsigned element_index_base = 2 * numElementsAcross * numElementsUp + 6 + 5 * (numElementsAcross - 1); // until frist lateral row
                element_index_base += (j - 1) * (3 * numElementsAcross + 2); // all previous lateral rows
                element_index_base += (i == 0 ? 0 : (5 + (i - 1) * 3)); // previous elements in this row

                face_orientations = std::vector<bool>(6, false);

                // First element in row only shares one edge
                if (i == 0)
                {
                    // Copy old pointer
                    p_lateral_faces[0] = elements[(j - 1) * numElementsAcross]->GetFace(5);
                    face_orientations[0] = true;
                    // New left and bottom edge face
                    p_lateral_faces[4] = new MonolayerVertexElement<2, 3>(element_index_base + 3, MonolayerVertexElementType::Lateral, nodes_lateral[4], lateral_node_types);
                    p_lateral_faces[5] = new MonolayerVertexElement<2, 3>(element_index_base + 4, MonolayerVertexElementType::Lateral, nodes_lateral[5], lateral_node_types);
                }
                // Inner and last elements share left and bottom edges
                else
                {
                    // Copy old pointer
                    p_lateral_faces[0] = elements[(j - 1) * numElementsAcross + i]->GetFace(5);
                    p_lateral_faces[5] = elements[(j - 1) * numElementsAcross + i - 1]->GetFace(4);
                    p_lateral_faces[4] = p_lateral_faces[1]; // left edge
                    face_orientations[0] = true;
                    face_orientations[5] = true;
                    face_orientations[4] = true;
                }
                // Create new faces
                p_lateral_faces[1] = new MonolayerVertexElement<2, 3>(element_index_base, MonolayerVertexElementType::Lateral, nodes_lateral[1], lateral_node_types);
                p_lateral_faces[2] = new MonolayerVertexElement<2, 3>(element_index_base + 1, MonolayerVertexElementType::Lateral, nodes_lateral[2], lateral_node_types);
                p_lateral_faces[3] = new MonolayerVertexElement<2, 3>(element_index_base + 2, MonolayerVertexElementType::Lateral, nodes_lateral[3], lateral_node_types);
            }

            // Set the boundary face property
            // In first row all lower faces (#0,#5) and outer left (#3,#4) and outer right (#1) boundary
            //    | X | X | X | uneven==1
            //  3/ \ / \ / \ /
            // 4| X | X | X |1   even==0
            //   \ / \ / \ /
            //  5  0 5 0 5  0
            if (j == 0)
            {
                if (i == 0)
                {
                    p_lateral_faces[0]->SetAsBoundaryFace();
                    p_lateral_faces[3]->SetAsBoundaryFace();
                    p_lateral_faces[4]->SetAsBoundaryFace();
                    p_lateral_faces[5]->SetAsBoundaryFace();
                }
                if (i == numElementsAcross - 1)
                {
                    p_lateral_faces[0]->SetAsBoundaryFace();
                    p_lateral_faces[1]->SetAsBoundaryFace();
                    p_lateral_faces[5]->SetAsBoundaryFace();
                }
                else
                {
                    p_lateral_faces[0]->SetAsBoundaryFace();
                    p_lateral_faces[5]->SetAsBoundaryFace();
                }
            }
            // Elements in uneven rows have 4 boundary faces
            //   | X | X | X |     even
            //    \ / \ / \ / \2
            //    4| X | X | X |1  uneven
            //    / \ / \ / \ /0
            //   | X | X | X |     even
            else if (j % 2 == 1)
            {
                if (i == 0)
                {
                    p_lateral_faces[4]->SetAsBoundaryFace();
                }
                else if (i == numElementsAcross - 1)
                {
                    p_lateral_faces[0]->SetAsBoundaryFace();
                    p_lateral_faces[1]->SetAsBoundaryFace();
                    p_lateral_faces[2]->SetAsBoundaryFace();
                }
            }
            // Elements in even rows have 4 boundary faces
            //    | X | X | X | uneven
            //  3/ \ / \ / \ /
            // 4| X | X | X |1   even
            //  5\ / \ / \ / \              *
            //    | X | X | X |  uneven
            else if (j % 2 == 0)
            {
                if (i == 0)
                {
                    p_lateral_faces[3]->SetAsBoundaryFace();
                    p_lateral_faces[4]->SetAsBoundaryFace();
                    p_lateral_faces[5]->SetAsBoundaryFace();
                }
                else if (i == numElementsAcross - 1)
                {
                    p_lateral_faces[1]->SetAsBoundaryFace();
                }
            }
            // Elements in last rows have additional upper boundary faces
            // for even we have one left face already
            //  3/ \ / \ / \                *
            // 4| X | X | X |1   even
            //  5\ / \ / \ / \              *
            //    | X | X | X |  uneven
            //
            // for uneven we have one right face already
            //      / \ / \ / \2
            //    4| X | X | X |1  uneven
            //    / \ / \ / \ /0
            //   | X | X | X |     even
            if (j == numElementsUp - 1)
            {
                if (i == 0 and j % 2 == 0)
                {
                    p_lateral_faces[2]->SetAsBoundaryFace();
                    p_lateral_faces[3]->SetAsBoundaryFace();
                }
                else if (i == numElementsAcross - 1 and j % 2 == 1)
                {
                    p_lateral_faces[3]->SetAsBoundaryFace();
                    p_lateral_faces[2]->SetAsBoundaryFace();
                }
                else
                {
                    p_lateral_faces[2]->SetAsBoundaryFace();
                    p_lateral_faces[3]->SetAsBoundaryFace();
                }
            }

            // Now create the element with the correct faces

            std::vector<MonolayerVertexElement<2, 3>*> faces_in_element;
            faces_in_element.push_back(faces_apical[j * numElementsAcross + i]);
            faces_in_element.push_back(faces_basal[j * numElementsAcross + i]);
            faces_in_element.push_back(p_lateral_faces[0]);
            faces_in_element.push_back(p_lateral_faces[1]);
            faces_in_element.push_back(p_lateral_faces[2]);
            faces_in_element.push_back(p_lateral_faces[3]);
            faces_in_element.push_back(p_lateral_faces[4]);
            faces_in_element.push_back(p_lateral_faces[5]);

            std::vector<bool> face_orientations_in_element(8);
            face_orientations_in_element[0] = false;
            face_orientations_in_element[1] = false;
            face_orientations_in_element[2] = face_orientations[0];
            face_orientations_in_element[3] = face_orientations[1];
            face_orientations_in_element[4] = face_orientations[2];
            face_orientations_in_element[5] = face_orientations[3];
            face_orientations_in_element[6] = face_orientations[4];
            face_orientations_in_element[7] = face_orientations[5];

            elements.push_back(new MonolayerVertexElement<3, 3>(j * numElementsAcross + i, MonolayerVertexElementType::Undetermined, faces_in_element, face_orientations_in_element));
        }
    }

    std::vector<Node<3>*> nodes = nodes_apical;
    nodes.insert(nodes.end(), nodes_basal.begin(), nodes_basal.end());

    mpMesh = new MutableMonolayerVertexMesh<3, 3>(nodes, elements, cellRearrangementThreshold, t2Threshold);

    // Save the edges
    mNodesLeftEdge = nodes_left;
    mNodesTopEdge = nodes_top;
    mNodesRightEdge = nodes_right;
    mNodesBottomEdge = nodes_bottom;

    // Scale the mesh so that each element's area takes the value elementArea
    mpMesh->Scale(sqrt(elementArea * 2.0 / sqrt(3.0)), sqrt(elementArea * 2.0 / sqrt(3.0)), 1.0);
}

void FiniteThicknessHoneycombVertexMeshGenerator::CurveMonolayer(double relativeRadius)
{
    double latticeConstant = sqrt(mElementArea * 2.0 / 3.0 / sqrt(3.0));
    double xTotal = (double(mNumElementsAcross) + 0.5) * sqrt(3.0) * latticeConstant;
    double yTotal = (3.0 * double(mNumElementsUp) + 1.0) * latticeConstant / 2.0;

    // Radius is relativeRadius * mean of dimension
    double radius = relativeRadius * (xTotal + yTotal) / 2.0;

    unsigned numNodes = mpMesh->GetNumNodes();

    for (unsigned indexNode = 0; indexNode < numNodes; ++indexNode)
    {
        c_vector<double, 3>& node_position_vector = mpMesh->GetNode(indexNode)->rGetModifiableLocation();
        double x_rel = node_position_vector[0] - xTotal / 2.0;
        double y_rel = node_position_vector[1] - yTotal / 2.0;
        double radius_new = radius + node_position_vector[2];

        double in_plane_radius = sqrt(x_rel * x_rel + y_rel * y_rel);

        node_position_vector[2] -= radius_new * (1.0 - (radius_new / in_plane_radius) / sqrt(1.0 + (radius_new / in_plane_radius) * (radius_new / in_plane_radius)));
    }
}

void FiniteThicknessHoneycombVertexMeshGenerator::CurveMonolayerCylinder(c_vector<double, 3> midPoint, c_vector<double, 3> midAxis, double radius)
{
    unsigned numNodes = mpMesh->GetNumNodes();
    // radius -= mHeight / 2.0;
    radius = radius < 1e-3 ? 1e-2 : radius;

    c_vector<double, 3> z_vector = zero_vector<double>(3);
    z_vector[2] = 1.0;

    midAxis /= norm_2(midAxis);
    c_vector<double, 3> azimuthal_vector = VectorProduct(z_vector, midAxis);
    azimuthal_vector /= norm_2(azimuthal_vector);

    for (unsigned indexNode = 0; indexNode < numNodes; ++indexNode)
    {
        c_vector<double, 3>& node_position_vector = mpMesh->GetNode(indexNode)->rGetModifiableLocation();

        double local_radius = radius + node_position_vector[2] - mHeight / 2.0;
        c_vector<double, 3> in_plane_position = (node_position_vector - midPoint) - inner_prod((node_position_vector - midPoint), z_vector) * z_vector;
        double new_x = inner_prod(in_plane_position, azimuthal_vector);
        double phi = new_x / (2.0 * M_PI * radius);
        double new_y = inner_prod(in_plane_position, midAxis);

        c_vector<double, 3> new_position = zero_vector<double>(3);
        new_position += midPoint + local_radius * cos(phi * 2.0 * M_PI) * z_vector;
        new_position += local_radius * sin(phi * 2.0 * M_PI) * azimuthal_vector;
        new_position += new_y * midAxis;

        node_position_vector = new_position;
    }
}

void FiniteThicknessHoneycombVertexMeshGenerator::CurveMonolayerSphere(c_vector<double, 3> midPoint, double radius)
{
    unsigned numNodes = mpMesh->GetNumNodes();
    radius -= mHeight / 2.0;
    radius = radius < 1e-3 ? 1e-2 : radius;

    c_vector<double, 3> z_vector = zero_vector<double>(3);
    z_vector[2] = 1.0;

    c_vector<double, 3> midAxis = zero_vector<double>(3);
    midAxis[0] = 1.0;
    c_vector<double, 3> azimuthal_vector = VectorProduct(z_vector, midAxis);
    azimuthal_vector /= norm_2(azimuthal_vector);

    for (unsigned indexNode = 0; indexNode < numNodes; ++indexNode)
    {
        c_vector<double, 3>& node_position_vector = mpMesh->GetNode(indexNode)->rGetModifiableLocation();

        double local_radius = radius + node_position_vector[2];

        // project onto upper half-sphere
        // note that this does not conserve the in-plane area,
        // but stereographic projection wouldn't conserve it either.
        c_vector<double, 3> projection_vector = node_position_vector - midPoint;
        projection_vector /= norm_2(projection_vector);

        c_vector<double, 3> new_position = zero_vector<double>(3);
        new_position = midPoint + local_radius * projection_vector;

        node_position_vector = new_position;
    }
}

void FiniteThicknessHoneycombVertexMeshGenerator::MakeCylindrical()
{
    // Fix possible errors in numbering
    mpMesh->RemoveDeletedNodes();
    mpMesh->RemoveDeletedFaces();
    mpMesh->UpdateElementsFacesMap();

    // First we identify the nodes and faces which are glued together
    std::vector<std::pair<Node<3>*, Node<3>*> > old_new_node_pairs;
    std::vector<MonolayerVertexElementType> node_pair_type;
    std::vector<std::pair<MonolayerVertexElement<2, 3>*, MonolayerVertexElement<2, 3>*> > old_new_face_pairs;

    for (unsigned j = 0; j < mNumElementsUp; ++j)
    {
        MonolayerVertexElement<3, 3>* p_right_element = mpMesh->GetElement((j + 1) * mNumElementsAcross - 1);
        MonolayerVertexElement<3, 3>* p_left_element = mpMesh->GetElement((j)*mNumElementsAcross);
        MonolayerVertexElement<3, 3>* p_left_lower_element = nullptr;
        MonolayerVertexElement<3, 3>* p_left_upper_element = nullptr;
        if (j > 0)
        {
            p_left_lower_element = mpMesh->GetElement((j - 1) * mNumElementsAcross);
        }
        if (j < mNumElementsUp - 1)
        {
            p_left_upper_element = mpMesh->GetElement((j + 1) * mNumElementsAcross);
        }

        // First the basal nodes
        MonolayerVertexElement<2, 3>* p_basal_face_right = p_right_element->GetFace(1);
        MonolayerVertexElement<2, 3>* p_basal_face_left = p_left_element->GetFace(1);
        Node<3>* p_old_node_b1 = p_basal_face_right->GetNode(4);
        Node<3>* p_new_node_b1 = p_basal_face_left->GetNode(2);
        Node<3>* p_old_node_b2 = p_basal_face_right->GetNode(5);
        Node<3>* p_new_node_b2 = p_basal_face_left->GetNode(1);

        std::pair<Node<3>*, Node<3>*> pair = std::pair<Node<3>*, Node<3>*>(p_old_node_b1, p_new_node_b1);
        if (std::find(old_new_node_pairs.begin(), old_new_node_pairs.end(), pair) == old_new_node_pairs.end())
        {
            old_new_node_pairs.push_back(pair);
            node_pair_type.push_back(MonolayerVertexElementType::Basal);
        }

        pair = std::pair<Node<3>*, Node<3>*>(p_old_node_b2, p_new_node_b2);
        if (std::find(old_new_node_pairs.begin(), old_new_node_pairs.end(), pair) == old_new_node_pairs.end())
        {
            old_new_node_pairs.push_back(pair);
            node_pair_type.push_back(MonolayerVertexElementType::Basal);
        }

        p_new_node_b1->SetAsBoundaryNode(false);
        p_new_node_b2->SetAsBoundaryNode(false);

        // If the last row is odd we need to treat it differently from other odds
        // We glue together the right and lower right faces
        if (j % 2 != 0)
        {
            MonolayerVertexElement<2, 3>* p_basal_face_left_lower = p_left_lower_element->GetFace(1);
            Node<3>* p_old_node_b3 = p_basal_face_right->GetNode(0);
            Node<3>* p_new_node_b3 = p_basal_face_left_lower->GetNode(2);

            pair = std::pair<Node<3>*, Node<3>*>(p_old_node_b3, p_new_node_b3);
            if (std::find(old_new_node_pairs.begin(), old_new_node_pairs.end(), pair) == old_new_node_pairs.end())
            {
                old_new_node_pairs.push_back(pair);
                node_pair_type.push_back(MonolayerVertexElementType::Basal);
            }

            p_new_node_b3->SetAsBoundaryNode(false);
        }
        // If the row is odd and not last
        // we glue together the right and upper + lower right faces
        if (j % 2 != 0 && j != mNumElementsUp - 1)
        {
            MonolayerVertexElement<2, 3>* p_basal_face_left_upper = p_left_upper_element->GetFace(1);
            Node<3>* p_old_node_b4 = p_basal_face_right->GetNode(3);
            Node<3>* p_new_node_b4 = p_basal_face_left_upper->GetNode(1);

            pair = std::pair<Node<3>*, Node<3>*>(p_old_node_b4, p_new_node_b4);
            if (std::find(old_new_node_pairs.begin(), old_new_node_pairs.end(), pair) == old_new_node_pairs.end())
            {
                old_new_node_pairs.push_back(pair);
                node_pair_type.push_back(MonolayerVertexElementType::Basal);
            }

            p_new_node_b4->SetAsBoundaryNode(false);
        }

        // Now the apical nodes
        MonolayerVertexElement<2, 3>* p_apical_face_right = p_right_element->GetFace(0);
        MonolayerVertexElement<2, 3>* p_apical_face_left = p_left_element->GetFace(0);
        Node<3>* p_old_node_a1 = p_apical_face_right->GetNode(1);
        Node<3>* p_new_node_a1 = p_apical_face_left->GetNode(5);
        Node<3>* p_old_node_a2 = p_apical_face_right->GetNode(2);
        Node<3>* p_new_node_a2 = p_apical_face_left->GetNode(4);

        pair = std::pair<Node<3>*, Node<3>*>(p_old_node_a1, p_new_node_a1);
        if (std::find(old_new_node_pairs.begin(), old_new_node_pairs.end(), pair) == old_new_node_pairs.end())
        {
            old_new_node_pairs.push_back(pair);
            node_pair_type.push_back(MonolayerVertexElementType::Apical);
        }
        pair = std::pair<Node<3>*, Node<3>*>(p_old_node_a2, p_new_node_a2);
        if (std::find(old_new_node_pairs.begin(), old_new_node_pairs.end(), pair) == old_new_node_pairs.end())
        {
            old_new_node_pairs.push_back(pair);
            node_pair_type.push_back(MonolayerVertexElementType::Apical);
        }

        p_new_node_a1->SetAsBoundaryNode(false);
        p_new_node_a2->SetAsBoundaryNode(false);

        if (j % 2 != 0)
        {
            MonolayerVertexElement<2, 3>* p_apical_face_left_lower = p_left_lower_element->GetFace(0);
            Node<3>* p_old_node_a3 = p_apical_face_right->GetNode(0);
            Node<3>* p_new_node_a3 = p_apical_face_left_lower->GetNode(4);

            pair = std::pair<Node<3>*, Node<3>*>(p_old_node_a3, p_new_node_a3);
            if (std::find(old_new_node_pairs.begin(), old_new_node_pairs.end(), pair) == old_new_node_pairs.end())
            {
                old_new_node_pairs.push_back(pair);
                node_pair_type.push_back(MonolayerVertexElementType::Apical);
            }

            p_new_node_a3->SetAsBoundaryNode(false);
        }
        if (j % 2 != 0 && j != mNumElementsUp - 1)
        {
            MonolayerVertexElement<2, 3>* p_apical_face_left_upper = p_left_upper_element->GetFace(0);

            Node<3>* p_old_node_a4 = p_apical_face_right->GetNode(3);
            Node<3>* p_new_node_a4 = p_apical_face_left_upper->GetNode(5);

            pair = std::pair<Node<3>*, Node<3>*>(p_old_node_a4, p_new_node_a4);
            if (std::find(old_new_node_pairs.begin(), old_new_node_pairs.end(), pair) == old_new_node_pairs.end())
            {
                old_new_node_pairs.push_back(pair);
                node_pair_type.push_back(MonolayerVertexElementType::Apical);
            }

            p_new_node_a4->SetAsBoundaryNode(false);
        }
        // Correct the boundary property for top/bottom nodes
        if (j == 0)
        {
            p_new_node_b2->SetAsBoundaryNode(true);
            p_new_node_a1->SetAsBoundaryNode(true);

            auto it = std::find(mNodesBottomEdge.begin(), mNodesBottomEdge.end(), p_old_node_b2);
            if (it != mNodesBottomEdge.end())
            {
                mNodesBottomEdge.erase(it);
            }
            it = std::find(mNodesBottomEdge.begin(), mNodesBottomEdge.end(), p_new_node_b2);
            if (it == mNodesBottomEdge.end())
            {
                mNodesBottomEdge.push_back(p_new_node_b2);
            }

            it = std::find(mNodesBottomEdge.begin(), mNodesBottomEdge.end(), p_old_node_a1);
            if (it != mNodesBottomEdge.end())
            {
                mNodesBottomEdge.erase(it);
            }
            it = std::find(mNodesBottomEdge.begin(), mNodesBottomEdge.end(), p_new_node_a1);
            if (it == mNodesBottomEdge.end())
            {
                mNodesBottomEdge.push_back(p_new_node_a1);
            }
        }
        if (j == mNumElementsUp - 1)
        {
            p_new_node_b1->SetAsBoundaryNode(true);
            p_new_node_a2->SetAsBoundaryNode(true);

            auto it = std::find(mNodesTopEdge.begin(), mNodesTopEdge.end(), p_old_node_b1);
            if (it != mNodesTopEdge.end())
            {
                mNodesTopEdge.erase(it);
            }
            it = std::find(mNodesTopEdge.begin(), mNodesTopEdge.end(), p_new_node_b1);
            if (it == mNodesTopEdge.end())
            {
                mNodesTopEdge.push_back(p_new_node_b1);
            }

            it = std::find(mNodesTopEdge.begin(), mNodesTopEdge.end(), p_old_node_a2);
            if (it != mNodesTopEdge.end())
            {
                mNodesTopEdge.erase(it);
            }
            it = std::find(mNodesTopEdge.begin(), mNodesTopEdge.end(), p_new_node_a2);
            if (it == mNodesTopEdge.end())
            {
                mNodesTopEdge.push_back(p_new_node_a2);
            }
        }

        // Sanity check: No node should be assigned to two new nodes
        std::set<Node<3>*> set_nodes_assigned;
        for (auto it = old_new_node_pairs.begin(); it != old_new_node_pairs.end(); ++it)
        {
            if (set_nodes_assigned.find(it->first) == set_nodes_assigned.end())
            {
                set_nodes_assigned.insert(it->first);
            }
            else
            {
                EXCEPTION("Tried to assign multiple new nodes to one old node when creating cylindrical monolayer.");
            }
        }

        // Now the lateral faces
        MonolayerVertexElement<2, 3>* p_face_to_keep_mid = p_left_element->GetFace(6);
        MonolayerVertexElement<2, 3>* p_face_to_delete_mid = p_right_element->GetFace(3);
        MonolayerVertexElement<2, 3>* p_face_to_keep_lower = nullptr;
        MonolayerVertexElement<2, 3>* p_face_to_delete_lower = nullptr;
        MonolayerVertexElement<2, 3>* p_face_to_keep_upper = nullptr;
        MonolayerVertexElement<2, 3>* p_face_to_delete_upper = nullptr;

        p_face_to_keep_mid->SetAsBoundaryFace(false);

        if (j % 2 != 0)
        {
            p_face_to_keep_lower = p_left_lower_element->GetFace(5);
            p_face_to_delete_lower = p_right_element->GetFace(2);

            // old_new_face_pairs(std::pair<MonolayerVertexElement<2, 3>*, MonolayerVertexElement<2, 3>*>(
            //     p_face_to_delete_lower, p_face_to_keep_lower));

            p_face_to_keep_lower->SetAsBoundaryFace(false);
        }
        if (j % 2 != 0 && j != mNumElementsUp - 1)
        {
            p_face_to_keep_upper = p_left_upper_element->GetFace(7);
            p_face_to_delete_upper = p_right_element->GetFace(4);

            // old_new_face_pairs(std::pair<MonolayerVertexElement<2, 3>*, MonolayerVertexElement<2, 3>*>(
            //     p_face_to_delete_upper, p_face_to_keep_upper));

            p_face_to_keep_upper->SetAsBoundaryFace(false);
        }

        // Now glue together the faces
        p_right_element->DeleteFace(p_right_element->GetFaceLocalIndex(p_face_to_delete_mid));
        p_right_element->AddFace(p_face_to_keep_mid, MonolayerVertexElementType::Lateral, true);
        mpMesh->DeleteFacePriorToReMesh(p_face_to_delete_mid->GetIndex());
        if (p_face_to_delete_lower != nullptr)
        {
            p_right_element->DeleteFace(p_right_element->GetFaceLocalIndex(p_face_to_delete_lower));
            p_right_element->AddFace(p_face_to_keep_lower, MonolayerVertexElementType::Lateral, true);
            mpMesh->DeleteFacePriorToReMesh(p_face_to_delete_lower->GetIndex());
        }
        if (p_face_to_delete_upper != nullptr)
        {
            p_right_element->DeleteFace(p_right_element->GetFaceLocalIndex(p_face_to_delete_upper));
            p_right_element->AddFace(p_face_to_keep_upper, MonolayerVertexElementType::Lateral, true);
            mpMesh->DeleteFacePriorToReMesh(p_face_to_delete_upper->GetIndex());
        }
    }

    // Now we glue together the nodes
    unsigned num_elements = mpMesh->GetNumElements();
    for (unsigned index_element = 0; index_element < num_elements; ++index_element)
    {
        auto it_type = node_pair_type.begin();
        for (auto it_pair = old_new_node_pairs.begin(); it_pair != old_new_node_pairs.end(); ++it_pair)
        {
            MonolayerVertexElement<3, 3>* p_element = mpMesh->GetElement(index_element);
            unsigned local_index = p_element->GetNodeLocalIndex((it_pair)->first->GetIndex());
            if (local_index != UINT_MAX)
            {
                p_element->DeleteNode(local_index);
            }

            unsigned num_faces = p_element->GetNumFaces();
            for (unsigned index_face = 0; index_face < num_faces; ++index_face)
            {
                MonolayerVertexElement<2, 3>* p_face = p_element->GetFace(index_face);
                local_index = p_face->GetNodeLocalIndex((it_pair)->first->GetIndex());
                if (local_index != UINT_MAX)
                {
                    p_face->ReplaceNode((it_pair)->first, (it_pair)->second, (*it_type));
                }
            }
            ++it_type;
        }
    }

    // Now mark the nodes for deletion
    for (auto it_pair = old_new_node_pairs.begin(); it_pair != old_new_node_pairs.end(); ++it_pair)
    {
        unsigned index_delete = (it_pair)->first->GetIndex();
        mpMesh->DeleteNodePriorToReMesh(index_delete);
    }

    // Now remove the unused nodes and faces and update the elements faces map of the mesh
    mpMesh->RemoveDeletedNodes();
    mpMesh->RemoveDeletedFaces();
    mpMesh->UpdateElementsFacesMap();

    // After glueing we now deform into a corresponding cylinder - inside is basal
    double latticeConstant = sqrt(mElementArea * 2.0 / 3.0 / sqrt(3.0));
    double xTotal = (double(mNumElementsAcross)) * sqrt(3.0) * latticeConstant;

    double radius_mid = xTotal / 2.0 / M_PI;
    // If height is too large, we increase the radius to accomodate the height.
    radius_mid = (radius_mid < mHeight / 2.0) ? mHeight / 2.0 + 1e-2 : radius_mid;

    unsigned num_nodes = mpMesh->GetNumNodes();
    for (unsigned index_node = 0; index_node < num_nodes; ++index_node)
    {
        c_vector<double, 3>& node_position_vector = mpMesh->GetNode(index_node)->rGetModifiableLocation();
        double node_radius = radius_mid - mHeight / 2.0 + node_position_vector[2];
        double node_angle = node_position_vector[0] / xTotal * 2.0 * M_PI;

        node_position_vector[0] = node_radius * sin(node_angle);
        node_position_vector[2] = node_radius * cos(node_angle);
    }
    // Now just correct the boundary edge list
    mNodesLeftEdge.clear();
    mNodesRightEdge.clear();
}

std::array<std::vector<Node<3>*>, 4> FiniteThicknessHoneycombVertexMeshGenerator::GetEdgesWithNodes()
{
    std::array<std::vector<Node<3>*>, 4> edge_array;
    edge_array[0] = mNodesLeftEdge;
    edge_array[1] = mNodesTopEdge;
    edge_array[2] = mNodesRightEdge;
    edge_array[3] = mNodesBottomEdge;

    return edge_array;
}

FiniteThicknessHoneycombVertexMeshGenerator::~FiniteThicknessHoneycombVertexMeshGenerator()
{
    delete mpMesh;
}

MutableMonolayerVertexMesh<3, 3>* FiniteThicknessHoneycombVertexMeshGenerator::GetMesh()
{
    return mpMesh;
}