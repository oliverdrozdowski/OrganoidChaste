#include "FiniteThicknessIcosahedralSphereMeshGenerator.hpp"
#include "MonolayerVertexElement.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"
#include "VoronoiSphereGenerator.hpp"

#include <boost/numeric/ublas/assignment.hpp>
#include <cmath> //for M_PI

FiniteThicknessIcosahedralSphereMeshGenerator::FiniteThicknessIcosahedralSphereMeshGenerator(unsigned casparKlugNumH,
                                                                                             unsigned casparKlugNumK,
                                                                                             double cellRearrangementThreshold,
                                                                                             double t2Threshold,
                                                                                             double height,
                                                                                             double innerRadius,
                                                                                             bool makeSphere)
        : mInnerRadius(innerRadius),
          mHeight(height)
{
    unsigned triangulation_number = (casparKlugNumH * casparKlugNumH + casparKlugNumK + casparKlugNumK + casparKlugNumK * casparKlugNumH);

    this->mNumElements = 10 * (triangulation_number - 1);

    // assert(cellRearrangementThreshold >= 0.0);
    assert(t2Threshold > 0.0);
    assert(height > 0.0);
    assert(innerRadius > 0.0);

    // We always consider index H >= K
    if (casparKlugNumH < casparKlugNumK)
    {
        unsigned tmp = casparKlugNumH;
        casparKlugNumH = casparKlugNumK;
        casparKlugNumK = tmp;
    }

    double h = (double)casparKlugNumH;
    double k = (double)casparKlugNumK;

    // First generate pentagonal cell centers on the sphere.
    std::vector<c_vector<double, 3> > cell_centers;
    std::vector<c_vector<double, 3> > icosahedron_face_normals;

    double angles[12][2] = {
        { 0.0, 0.0 },
        { M_PI / 2.0 - 0.4636476, 0.0 * M_PI },
        { M_PI / 2.0 - 0.4636476, 0.4 * M_PI },
        { M_PI / 2.0 - 0.4636476, 0.8 * M_PI },
        { M_PI / 2.0 - 0.4636476, 1.2 * M_PI },
        { M_PI / 2.0 - 0.4636476, 1.6 * M_PI },
        { M_PI / 2.0 + 0.4636476, 0.2 * M_PI },
        { M_PI / 2.0 + 0.4636476, 0.6 * M_PI },
        { M_PI / 2.0 + 0.4636476, 1.0 * M_PI },
        { M_PI / 2.0 + 0.4636476, 1.4 * M_PI },
        { M_PI / 2.0 + 0.4636476, 1.8 * M_PI },
        { M_PI, 0.0 }
    };

    int faces[20][3] = {
        { 0, 1, 2 }, // upper five
        { 0, 2, 3 },
        { 0, 3, 4 },
        { 0, 4, 5 },
        { 0, 5, 1 },
        { 1, 6, 2 }, // upper central five
        { 2, 7, 3 },
        { 3, 8, 4 },
        { 4, 9, 5 },
        { 5, 10, 1 },
        { 6, 7, 2 }, // lower central five
        { 7, 8, 3 },
        { 8, 9, 4 },
        { 9, 10, 5 },
        { 10, 6, 1 },
        { 11, 7, 6 }, // lower five
        { 11, 8, 7 },
        { 11, 9, 8 },
        { 11, 10, 9 },
        { 11, 6, 10 }
    };

    int edges[30][2] = {
        { 0, 1 }, { 0, 2 }, { 0, 3 }, { 0, 4 }, { 0, 5 }, { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 5 }, { 5, 1 }, { 1, 6 }, { 6, 2 }, { 2, 7 }, { 7, 3 }, { 3, 8 }, { 8, 4 }, { 4, 9 }, { 9, 5 }, { 5, 10 }, { 10, 1 }, { 6, 7 }, { 7, 8 }, { 8, 9 }, { 9, 10 }, { 10, 6 }, { 6, 11 }, { 7, 11 }, { 8, 11 }, { 9, 11 }, { 10, 11 }
    };

    std::vector<c_vector<double, 3> > edge_vectors;
    // We calculate the pentagonal centers
    for (unsigned index_center = 0; index_center < 12; ++index_center)
    {
        c_vector<double, 3> point;
        double theta = angles[index_center][0];
        double phi = angles[index_center][1];
        point <<= mInnerRadius * sin(theta) * cos(phi),
            mInnerRadius * sin(theta) * sin(phi),
            mInnerRadius * cos(theta);
        cell_centers.push_back(point);
    }

    // Then we add the edge points
    // for h=k and k=0 we have cells directly on the edges
    if (casparKlugNumK == 0 or casparKlugNumK == casparKlugNumH)
    {
        for (unsigned edge_index = 0; edge_index < 30; ++edge_index)
        {
            if (casparKlugNumH <= 0)
                break;
            unsigned upper_index = casparKlugNumH - 1;
            c_vector<double, 3> pt1 = cell_centers[edges[edge_index][0]];
            c_vector<double, 3> pt2 = cell_centers[edges[edge_index][1]];
            c_vector<double, 3> diff = pt2 - pt1;
            for (unsigned sub_ind = 0; sub_ind < upper_index; ++sub_ind)
            {
                double x_fraction = (double)(sub_ind + 1) / (double)(upper_index + 1);
                c_vector<double, 3> new_pt = pt1 + x_fraction * diff;
                new_pt *= mInnerRadius / norm_2(new_pt);
                cell_centers.push_back(new_pt);
                edge_vectors.push_back(new_pt);
            }
        }
    }

    // First fill the faces inside
    for (unsigned face_index = 0; face_index < 20; ++face_index)
    {
        // Calculate basis vectors of face from pentagon to pentagon
        c_vector<double, 3> pt1 = cell_centers[faces[face_index][0]];
        c_vector<double, 3> pt2 = cell_centers[faces[face_index][1]];
        c_vector<double, 3> pt3 = cell_centers[faces[face_index][2]];
        c_vector<double, 3> A_basis = pt2 - pt1;
        c_vector<double, 3> B_basis = pt3 - pt1;
        c_vector<double, 3> normal_face = VectorProduct(A_basis, B_basis);
        normal_face /= norm_2(normal_face);
        normal_face *= inner_prod((pt1 + pt2 + pt3) / 3.0, normal_face) > 0.0 ? 1.0 : -1.0; // reorient if points inwards
        icosahedron_face_normals.push_back(normal_face);

        // Compute hexagonal basis vectors in pentagonal basis representation,
        // where A = h*a1 + k*a2 and B = h*b + k*a1:
        //
        //        a2     a1
        //        ^     ,
        //     .  |  . /
        //        |   /
        //  .            .
        //            \        :
        //     .     . \       :
        //              \      :
        //                b

        c_vector<double, 3> a1_basis = (k * B_basis + h * A_basis) / (h * k + k * k + h * h);
        c_vector<double, 3> a2_basis = ((h + k) * A_basis - h * B_basis) / (h * k + k * k + h * h);

        for (int i_n = 0; i_n < h + k + 1; ++i_n)
        {
            double n = i_n * 1.0;
            for (int i_m = -h - k - 1; i_m < h + k + 1; ++i_m)
            {
                double m = i_m * 1.0;
                // Check if hexagon is inside the face
                double A_fraction = (n * h + m * (h + k)) / (h * k + k * k + h * h);
                double B_fraction = (n * k - h * m) / (h * k + k * k + h * h);
                bool A_range = A_fraction >= 1.0 or A_fraction <= 0.0;
                bool B_range = B_fraction >= 1.0 or B_fraction <= 0.0;
                bool sum_range = A_fraction + B_fraction >= 1.0;
                if (A_range or B_range or sum_range)
                {
                    continue;
                }
                c_vector<double, 3> new_pt = pt1 + n * a1_basis + m * a2_basis;
                new_pt *= mInnerRadius / norm_2(new_pt);
                cell_centers.push_back(new_pt);
            }
        }
    }

    // Then add the edge cells, which we have not added yet.
    for (unsigned face_index = 0; face_index < 20; ++face_index)
    {
        // Calculate basis vectors of face from pentagon to pentagon
        c_vector<double, 3> pt1 = cell_centers[faces[face_index][0]];
        c_vector<double, 3> pt2 = cell_centers[faces[face_index][1]];
        c_vector<double, 3> pt3 = cell_centers[faces[face_index][2]];
        c_vector<double, 3> A_basis = pt2 - pt1;
        c_vector<double, 3> B_basis = pt3 - pt1;
        c_vector<double, 3> a1_basis = (k * B_basis + h * A_basis) / (h * k + k * k + h * h);
        c_vector<double, 3> a2_basis = ((h + k) * A_basis - h * B_basis) / (h * k + k * k + h * h);
        for (int i_n = 0; i_n < h + k + 1; ++i_n)
        {
            double n = i_n * 1.0;
            for (int i_m = -h - k - 1; i_m < h + k + 1; ++i_m)
            {
                double m = i_m * 1.0;
                double epsilon = norm_2(a1_basis) * 0.8;
                // Check if hexagon is inside the face
                double A_fraction = (n * h + m * (h + k)) / (h * k + k * k + h * h);
                double B_fraction = (n * k - h * m) / (h * k + k * k + h * h);
                bool A_range = A_fraction >= 1.0 or A_fraction <= 0.0;
                bool B_range = B_fraction >= 1.0 or B_fraction <= 0.0;
                bool sum_range = A_fraction + B_fraction >= 1.0;
                if (A_range or B_range or sum_range)
                {
                    // We add the point if it is very close to the edge
                    // but has not been added previously
                    bool B_close = (abs(B_fraction - 1.0) < 1.0 / (h + k) / 20 or abs(B_fraction - 0.0) < 1.0 / (h + k) / 20)
                        and (A_fraction <= 1.0 and A_fraction >= 0.0);
                    bool sum_close = abs(A_fraction + B_fraction - 1.0) < 1.0 / (h + k) / 20
                        and (B_fraction <= 1.0 and B_fraction >= 0.0) and (A_fraction <= 1.0 and A_fraction >= 0.0);
                    if (B_close or sum_close)
                    {
                        c_vector<double, 3> new_pt = pt1 + n * a1_basis + m * a2_basis;
                        new_pt *= mInnerRadius / norm_2(new_pt);

                        const auto findr = std::find_if(
                            cell_centers.begin(),
                            cell_centers.end(),
                            [&new_pt, epsilon](const auto& findr)
                            {
                                return abs(norm_2(findr - new_pt)) < epsilon; // Comparing with the object
                            });
                        if (findr == cell_centers.end())
                            cell_centers.push_back(new_pt);
                    }
                }
            }
        }
    }

    // The VoronoiSphereGenerator gives the voronoi tesselation on the sphere for our vertices.
    // Note that the vertices are not on the sphere anymore and need to be rescaled accordingly!
    VoronoiSphereGenerator voronoi_generator(cell_centers);
    std::vector<c_vector<double, 3> > inner_voronoi_vertices = voronoi_generator.GetVoronoiVertices();

    unsigned node_index = 0;

    // unsigned num_inner_nodes = inner_voronoi_vertices.size();

    // inner is apical
    std::vector<Node<3>*> nodes_apical;
    for (auto vertex = inner_voronoi_vertices.begin(); vertex != inner_voronoi_vertices.end(); ++vertex)
    {
        // No boundary nodes as sphere has no gaps
        double new_radius = mInnerRadius / norm_2(*vertex);
        Node<3>* p_node = new Node<3>(node_index, false, new_radius * (*vertex)[0], new_radius * (*vertex)[1], new_radius * (*vertex)[2]);
        nodes_apical.push_back(p_node);
        node_index++;
    }

    // outer is basal
    std::vector<Node<3>*> nodes_basal;
    for (auto vertex = inner_voronoi_vertices.begin(); vertex != inner_voronoi_vertices.end(); ++vertex)
    {
        // No boundary nodes as sphere has no gaps
        double new_radius = (mHeight + mInnerRadius) / norm_2(*vertex);
        Node<3>* p_node = new Node<3>(node_index, false, new_radius * (*vertex)[0], new_radius * (*vertex)[1], new_radius * (*vertex)[2]);
        nodes_basal.push_back(p_node);
        node_index++;
    }

    std::vector<MonolayerVertexElement<2, 3>*> faces_apical;
    std::vector<bool> orientations_apical;
    std::vector<MonolayerVertexElement<2, 3>*> faces_basal;
    std::vector<bool> orientations_basal;
    std::vector<MonolayerVertexElement<2, 3>*> faces_lateral;
    std::vector<MonolayerVertexElement<3, 3>*> elements;

    unsigned face_index = 0;

    // Create apical faces
    std::map<unsigned, std::vector<unsigned> > map_apical_faces_vertices = voronoi_generator.GetMapFromFacesToVertices();
    unsigned num_apical_faces = map_apical_faces_vertices.size();

    for (auto mapping = map_apical_faces_vertices.begin(); mapping != map_apical_faces_vertices.end(); ++mapping)
    {
        std::vector<unsigned> indices_nodes = mapping->second;
        std::vector<Node<3>*> face_nodes;
        std::vector<MonolayerVertexElementType> face_node_types;
        for (auto node_index = indices_nodes.begin(); node_index != indices_nodes.end(); ++node_index)
        {
            face_nodes.push_back(nodes_apical[*node_index]);
            face_node_types.push_back(MonolayerVertexElementType::Apical);
        }
        MonolayerVertexElement<2, 3>* p_face_apical = new MonolayerVertexElement<2, 3>(face_index, MonolayerVertexElementType::Apical, face_nodes, face_node_types);
        face_index++;
        faces_apical.push_back(p_face_apical);
        // Determine the orientation
        c_vector<double, 3> center = GetCenterOfFace(p_face_apical);
        c_vector<double, 3> norml = GetNormalToFace(p_face_apical);
        bool points_outwards = inner_prod(center, norml) < 0.0 ? false : true;
        // If it points outwards we have to switch the direction of the nodes in the face
        if (points_outwards)
            p_face_apical->SwitchOrientation();
        orientations_apical.push_back(false); // wrong orientation == false
    }

    // Create basal faces
    std::map<unsigned, std::vector<unsigned> > map_basal_faces_vertices = voronoi_generator.GetMapFromFacesToVertices();
    unsigned num_basal_faces = map_basal_faces_vertices.size();
    assert(num_basal_faces == num_apical_faces);

    for (auto mapping = map_basal_faces_vertices.begin(); mapping != map_basal_faces_vertices.end(); ++mapping)
    {
        std::vector<unsigned> indices_nodes = mapping->second;
        std::vector<Node<3>*> face_nodes;
        std::vector<MonolayerVertexElementType> face_node_types;
        for (auto node_index = indices_nodes.begin(); node_index != indices_nodes.end(); ++node_index)
        {
            face_nodes.push_back(nodes_basal[*node_index]);
            face_node_types.push_back(MonolayerVertexElementType::Basal);
        }
        MonolayerVertexElement<2, 3>* p_face_basal = new MonolayerVertexElement<2, 3>(face_index, MonolayerVertexElementType::Basal, face_nodes, face_node_types);
        face_index++;
        faces_basal.push_back(p_face_basal);
        // Determine the orientation
        c_vector<double, 3> center = GetCenterOfFace(p_face_basal);
        c_vector<double, 3> norml = GetNormalToFace(p_face_basal);
        bool points_inwards = inner_prod(center, norml) < 0.0 ? true : false;
        // If it points inwards we have to switch the direction of the nodes in the face
        if (points_inwards)
            p_face_basal->SwitchOrientation();
        orientations_basal.push_back(false); // wrong orientation == false
    }

    // Create lateral faces
    // Map from element index to vector of pairs with lateral index and if orientation wrong
    std::map<unsigned, std::vector<std::pair<unsigned, bool> > > map_lateral_faces;
    // We need to map the indices of apical edges to face pointers to reuse the faces
    // and to check if we already created a lateral face
    std::map<std::set<unsigned>, MonolayerVertexElement<2, 3>*> map_edge_to_face;

    // Loop through all elements
    for (auto mapping = map_apical_faces_vertices.begin(); mapping != map_apical_faces_vertices.end(); ++mapping)
    {
        unsigned element_index = mapping->first;
        std::vector<unsigned> indices_nodes = mapping->second;
        unsigned num_nodes = indices_nodes.size();
        unsigned node_a_index = indices_nodes[0];
        map_lateral_faces[element_index] = std::vector<std::pair<unsigned, bool> >();
        // Loop through edges
        for (unsigned local_ind = 1; local_ind < num_nodes + 1; local_ind++)
        {
            unsigned node_b_index = indices_nodes[(local_ind) % num_nodes];

            // Check if we already created this lateral face
            std::set<unsigned> edge{ node_a_index, node_b_index };
            if (map_edge_to_face.find(edge) == map_edge_to_face.end())
            {
                // If not found create lateral face
                std::vector<Node<3>*> face_nodes;
                std::vector<MonolayerVertexElementType> face_node_types;

                // Apical nodes
                face_nodes.push_back(nodes_apical[node_a_index]);
                face_nodes.push_back(nodes_basal[node_a_index]);
                face_node_types.push_back(MonolayerVertexElementType::Apical);
                face_node_types.push_back(MonolayerVertexElementType::Basal);

                // Basal nodes - reverse order for (anti-)clockwise orientation
                face_nodes.push_back(nodes_basal[node_b_index]);
                face_nodes.push_back(nodes_apical[node_b_index]);
                face_node_types.push_back(MonolayerVertexElementType::Basal);
                face_node_types.push_back(MonolayerVertexElementType::Apical);

                MonolayerVertexElement<2, 3>* p_face_lateral = new MonolayerVertexElement<2, 3>(face_index, MonolayerVertexElementType::Lateral, face_nodes, face_node_types);
                face_index++;
                faces_lateral.push_back(p_face_lateral);

                // For orientation to this element check normal w.r.t. center to face center
                c_vector<double, 3> center_elmt = GetCenterOfFace(faces_apical[element_index]) / 2.0;
                center_elmt += GetCenterOfFace(faces_basal[element_index]) / 2.0;
                c_vector<double, 3> center_lat_face = GetCenterOfFace(p_face_lateral);

                c_vector<double, 3> normal_lat_face = GetNormalToFace(p_face_lateral);
                bool points_inwards = inner_prod((center_lat_face - center_elmt), normal_lat_face) < 0.0 ? true : false;

                map_edge_to_face[edge] = p_face_lateral;
                map_lateral_faces[element_index].push_back(std::pair<unsigned, bool>(faces_lateral.size() - 1, points_inwards));
            }
            else
            {
                // If face already created just add the face to the element with the correct orientation
                MonolayerVertexElement<2, 3>* p_face_lateral = map_edge_to_face[edge];

                // For orientation to this element check normal w.r.t. center to face center
                c_vector<double, 3> center_elmt = GetCenterOfFace(faces_apical[element_index]) / 2.0;
                center_elmt += GetCenterOfFace(faces_basal[element_index]) / 2.0;
                c_vector<double, 3> center_lat_face = GetCenterOfFace(p_face_lateral);

                c_vector<double, 3> normal_lat_face = GetNormalToFace(p_face_lateral);
                bool points_inwards = inner_prod((center_lat_face - center_elmt), normal_lat_face) < 0.0 ? true : false;

                unsigned index_lateral_face = std::distance(faces_lateral.begin(),
                                                            std::find(faces_lateral.begin(), faces_lateral.end(), p_face_lateral));

                map_lateral_faces[element_index].push_back(std::pair<unsigned, bool>(index_lateral_face, points_inwards));
            }
            // Go to next edge
            node_a_index = node_b_index;
        }
    }

    // Now we can create the 3d elements with the faces.

    for (auto mapping = map_lateral_faces.begin(); mapping != map_lateral_faces.end(); ++mapping)
    {
        unsigned element_index = mapping->first;
        std::vector<std::pair<unsigned, bool> > faces_ind_orient = mapping->second;
        std::vector<MonolayerVertexElement<2, 3>*> lateral_faces_in_element;
        std::vector<bool> lateral_orientations_in_element;

        for (auto fc = faces_ind_orient.begin(); fc != faces_ind_orient.end(); ++fc)
        {
            lateral_faces_in_element.push_back(faces_lateral[fc->first]);
            lateral_orientations_in_element.push_back(fc->second);
        }

        // Now create the element with the correct faces

        std::vector<MonolayerVertexElement<2, 3>*> faces_in_element;
        faces_in_element.push_back(faces_apical[element_index]);
        faces_in_element.push_back(faces_basal[element_index]);
        faces_in_element.insert(faces_in_element.end(), lateral_faces_in_element.begin(), lateral_faces_in_element.end());

        std::vector<bool> face_orientations_in_element;
        face_orientations_in_element.push_back(orientations_apical[element_index]);
        face_orientations_in_element.push_back(orientations_basal[element_index]);
        face_orientations_in_element.insert(face_orientations_in_element.end(), lateral_orientations_in_element.begin(), lateral_orientations_in_element.end());

        elements.push_back(new MonolayerVertexElement<3, 3>(element_index, MonolayerVertexElementType::Undetermined, faces_in_element, face_orientations_in_element));
    }

    // Now create the mesh
    std::vector<Node<3>*> nodes = nodes_apical;
    nodes.insert(nodes.end(), nodes_basal.begin(), nodes_basal.end());

    // If we want a real icosahedron, we now map into the circumscribed icosahedron.
    if (!makeSphere)
        ProjectOntoCircumscribingIcosahedronAndRescale(nodes, icosahedron_face_normals, mInnerRadius + mHeight / 2.0);

    mpMesh = new MutableMonolayerVertexMesh<3, 3>(nodes, elements, cellRearrangementThreshold, t2Threshold);
}

c_vector<double, 3> FiniteThicknessIcosahedralSphereMeshGenerator::GetNormalToFace(MonolayerVertexElement<2, 3>* pFace)
{
    // As we are in 3D, the face must have at least three vertices
    assert(pFace->GetNumNodes() >= 3u);

    // calculate the center
    // unsigned num_nodes = pFace->GetNumNodes();
    c_vector<double, 3> center = GetCenterOfFace(pFace);

    // Calculate the normal as mean value of triangle normals
    c_vector<double, 3> normal_vector = zero_vector<double>(3);

    c_vector<double, 3> v_minus_v_pc = pFace->GetNode(0)->rGetLocation() - center;
    for (unsigned local_index = 0; local_index < pFace->GetNumNodes(); local_index++)
    {
        c_vector<double, 3> vnext_minus_v_pc = pFace->GetNode((local_index + 1) % pFace->GetNumNodes())->rGetLocation() - center;

        c_vector<double, 3> localNormal = VectorProduct(v_minus_v_pc, vnext_minus_v_pc);

        // normal direction is implemented as mean normal of triangles (note possible errors from orientation)
        normal_vector += localNormal;
        v_minus_v_pc = vnext_minus_v_pc;
    }
    double magnitude = norm_2(normal_vector);
    if (magnitude != 0.0)
    {
        // Normalize the normal vector
        normal_vector /= magnitude;
        // If all points are co-located, then magnitude==0.0 and there is potential for a floating point exception
        // here if we divide by zero, so we'll move on.
    }
    return normal_vector;
}

void FiniteThicknessIcosahedralSphereMeshGenerator::ProjectOntoCircumscribingIcosahedronAndRescale(
    std::vector<Node<3>*> vectorOfNodes, std::vector<c_vector<double, 3> > faceNormals, double prescribedOuterRadius)
{
    double max_point_distance = 0.0;
    // First we map onto the icosahedron
    for (auto iter_node = vectorOfNodes.begin(); iter_node != vectorOfNodes.end(); ++iter_node)
    {
        /* For each node we check the icosahedral normals to find the normal with largest scalar product to point.
         * Then we map the node to the circumscribing icosahedron by mapping it to the corresponding face.
         */
        c_vector<double, 3>& r_node_vector = (*iter_node)->rGetModifiableLocation();
        double initial_radius = norm_2(r_node_vector);
        double max_product = 0.0;
        for (auto iter_normal = faceNormals.begin(); iter_normal != faceNormals.end(); ++iter_normal)
        {
            double product = inner_prod(r_node_vector, *iter_normal);
            if (product > max_product)
            {
                max_product = product;
            }
        }
        double rescaling_factor = initial_radius / max_product;
        r_node_vector *= rescaling_factor;
        double point_distance = initial_radius * rescaling_factor;
        max_point_distance = point_distance > max_point_distance ? point_distance : max_point_distance;
    }

    // Now we rescale the points as such that the maximum distance is now the prescribedOuterRadius
    for (auto iter_node = vectorOfNodes.begin(); iter_node != vectorOfNodes.end(); ++iter_node)
    {
        c_vector<double, 3>& r_node_vector = (*iter_node)->rGetModifiableLocation();
        r_node_vector *= prescribedOuterRadius / max_point_distance;
    }
}

c_vector<double, 3> FiniteThicknessIcosahedralSphereMeshGenerator::GetCenterOfFace(MonolayerVertexElement<2, 3>* pFace)
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

FiniteThicknessIcosahedralSphereMeshGenerator::~FiniteThicknessIcosahedralSphereMeshGenerator()
{
    delete mpMesh;
}

MutableMonolayerVertexMesh<3, 3>* FiniteThicknessIcosahedralSphereMeshGenerator::GetMesh()
{
    return mpMesh;
}