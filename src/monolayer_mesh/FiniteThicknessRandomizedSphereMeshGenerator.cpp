#include "FiniteThicknessRandomizedSphereMeshGenerator.hpp"

#include <boost/numeric/ublas/assignment.hpp>
#include <cmath> //for M_PI

#include "MonolayerVertexElement.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasCustomFunctions.hpp"
#include "VoronoiSphereGenerator.hpp"

FiniteThicknessRandomizedSphereMeshGenerator::
    FiniteThicknessRandomizedSphereMeshGenerator(
        unsigned numElements, double cellRearrangementThreshold,
        double t2Threshold, double height, double innerRadius, bool doRandomSequentialAdsorption)
        : mNumElements(numElements), mInnerRadius(innerRadius), mHeight(height)
{
    assert(numElements > 0);
    assert(cellRearrangementThreshold >= 0.0);
    assert(t2Threshold > 0.0);
    assert(height > 0.0);
    assert(innerRadius > 0.0);

    // First generate random cell centers on the sphere.
    std::vector<c_vector<double, 3> > random_cell_centers;

    if (doRandomSequentialAdsorption)
    {
        // Magic number from J. Feder, J. Theor. Biol. 87(2), 237 (1980)
        // and Roshal et al., Phys. Rev. E 108, 024404 (2023) which is related to vertex jamming.
        double jamming_number = 0.547;
        double min_distance_between_cell_centers = 1.8 * mInnerRadius / sqrt(mNumElements / (4.0 * jamming_number));
        for (unsigned index_center = 0; index_center < mNumElements;
             ++index_center)
        {

            unsigned max_iterations = 2.0 * (10 + index_center);
            unsigned loop_counter = 0;

            c_vector<double, 3> rand_vector = GenerateRandomPointOnSphere();
            // In RSA we check for overlap, i.e., whether new point is too close to old points.
            bool overlap = DoesOverlap(rand_vector, random_cell_centers, min_distance_between_cell_centers);

            while (overlap && (loop_counter < max_iterations))
            {
                rand_vector = GenerateRandomPointOnSphere();
                overlap = DoesOverlap(rand_vector, random_cell_centers, min_distance_between_cell_centers);
                loop_counter++;
            }

            // If we cannot generate a good random new vector, we just use the old one and print a warning.
            if (loop_counter == max_iterations)
            {
                std::cout << "Warning! Cell center with index " << index_center << " did not find position without overlap." << std::endl
                          << "Random Sequential Adsorption was not fully done and cell centers may be too close." << std::endl;
            }

            random_cell_centers.push_back(rand_vector);
        }
    }
    else
    {
        for (unsigned index_center = 0; index_center < mNumElements;
             ++index_center)
        {
            c_vector<double, 3> rand_vector = GenerateRandomPointOnSphere();
            random_cell_centers.push_back(rand_vector);
        }
    }

    // The VoronoiSphereGenerator gives the voronoi tesselation on the sphere for
    // our vertices. Note that the vertices are not on the sphere anymore and need
    // to be rescaled accordingly!

    VoronoiSphereGenerator voronoi_generator(random_cell_centers);
    std::vector<c_vector<double, 3> > inner_voronoi_vertices = voronoi_generator.GetVoronoiVertices();

    unsigned node_index = 0;

    // unsigned num_inner_nodes = inner_voronoi_vertices.size();

    // inner is apical
    std::vector<Node<3>*> nodes_apical;
    for (auto vertex = inner_voronoi_vertices.begin();
         vertex != inner_voronoi_vertices.end(); ++vertex)
    {
        // No boundary nodes as sphere has no gaps
        double new_radius = mInnerRadius / norm_2(*vertex);
        Node<3>* p_node = new Node<3>(node_index, false, new_radius * (*vertex)[0],
                                      new_radius * (*vertex)[1], new_radius * (*vertex)[2]);
        nodes_apical.push_back(p_node);
        node_index++;
    }

    // outer is basal
    std::vector<Node<3>*> nodes_basal;
    for (auto vertex = inner_voronoi_vertices.begin();
         vertex != inner_voronoi_vertices.end(); ++vertex)
    {
        // No boundary nodes as sphere has no gaps
        double new_radius = (mHeight + mInnerRadius) / norm_2(*vertex);
        Node<3>* p_node = new Node<3>(node_index, false, new_radius * (*vertex)[0],
                                      new_radius * (*vertex)[1], new_radius * (*vertex)[2]);
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

    for (auto mapping = map_apical_faces_vertices.begin();
         mapping != map_apical_faces_vertices.end(); ++mapping)
    {
        std::vector<unsigned> indices_nodes = mapping->second;
        std::vector<Node<3>*> face_nodes;
        std::vector<MonolayerVertexElementType> face_node_types;
        for (auto node_index = indices_nodes.begin();
             node_index != indices_nodes.end(); ++node_index)
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
        // If it points outwards we have to switch the direction of the nodes in the
        // face
        if (points_outwards)
            p_face_apical->SwitchOrientation();
        orientations_apical.push_back(false); // wrong orientation == false
    }

    // Create basal faces
    std::map<unsigned, std::vector<unsigned> > map_basal_faces_vertices = voronoi_generator.GetMapFromFacesToVertices();

    assert(map_basal_faces_vertices.size() == map_apical_faces_vertices.size());

    for (auto mapping = map_basal_faces_vertices.begin();
         mapping != map_basal_faces_vertices.end(); ++mapping)
    {
        std::vector<unsigned> indices_nodes = mapping->second;
        std::vector<Node<3>*> face_nodes;
        std::vector<MonolayerVertexElementType> face_node_types;
        for (auto node_index = indices_nodes.begin();
             node_index != indices_nodes.end(); ++node_index)
        {
            face_nodes.push_back(nodes_basal[*node_index]);
            face_node_types.push_back(MonolayerVertexElementType::Basal);
        }
        MonolayerVertexElement<2, 3>* p_face_basal = new MonolayerVertexElement<2, 3>(face_index,
                                                                                      MonolayerVertexElementType::Basal,
                                                                                      face_nodes, face_node_types);
        face_index++;
        faces_basal.push_back(p_face_basal);
        // Determine the orientation
        c_vector<double, 3> center = GetCenterOfFace(p_face_basal);
        c_vector<double, 3> norml = GetNormalToFace(p_face_basal);
        bool points_inwards = inner_prod(center, norml) < 0.0 ? true : false;
        // If it points inwards we have to switch the direction of the nodes in the
        // face
        if (points_inwards)
            p_face_basal->SwitchOrientation();
        orientations_basal.push_back(false); // wrong orientation == false
    }

    // Create lateral faces
    // Map from element index to vector of pairs with lateral index and if
    // orientation wrong
    std::map<unsigned, std::vector<std::pair<unsigned, bool> > > map_lateral_faces;
    // We need to map the indices of apical edges to face pointers to reuse the
    // faces and to check if we already created a lateral face
    std::map<std::set<unsigned>, MonolayerVertexElement<2, 3>*> map_edge_to_face;

    // Loop through all elements
    for (auto mapping = map_apical_faces_vertices.begin();
         mapping != map_apical_faces_vertices.end(); ++mapping)
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

                MonolayerVertexElement<2, 3>* p_face_lateral = new MonolayerVertexElement<2, 3>(
                    face_index, MonolayerVertexElementType::Lateral, face_nodes,
                    face_node_types);
                face_index++;
                faces_lateral.push_back(p_face_lateral);

                // For orientation to this element check normal w.r.t. center to face
                // center
                c_vector<double, 3> center_elmt = GetCenterOfFace(faces_apical[element_index]) / 2.0;
                center_elmt += GetCenterOfFace(faces_basal[element_index]) / 2.0;
                c_vector<double, 3> center_lat_face = GetCenterOfFace(p_face_lateral);

                c_vector<double, 3> normal_lat_face = GetNormalToFace(p_face_lateral);
                bool points_inwards = inner_prod((center_lat_face - center_elmt), normal_lat_face) < 0.0
                    ? true
                    : false;

                map_edge_to_face[edge] = p_face_lateral;
                map_lateral_faces[element_index].push_back(std::pair<unsigned, bool>(
                    faces_lateral.size() - 1, points_inwards));
            }
            else
            {
                // If face already created just add the face to the element with the
                // correct orientation
                MonolayerVertexElement<2, 3>* p_face_lateral = map_edge_to_face[edge];

                // For orientation to this element check normal w.r.t. center to face
                // center
                c_vector<double, 3> center_elmt = GetCenterOfFace(faces_apical[element_index]) / 2.0;
                center_elmt += GetCenterOfFace(faces_basal[element_index]) / 2.0;
                c_vector<double, 3> center_lat_face = GetCenterOfFace(p_face_lateral);

                c_vector<double, 3> normal_lat_face = GetNormalToFace(p_face_lateral);
                bool points_inwards = inner_prod((center_lat_face - center_elmt), normal_lat_face) < 0.0
                    ? true
                    : false;

                unsigned index_lateral_face = std::distance(faces_lateral.begin(),
                                                            std::find(faces_lateral.begin(), faces_lateral.end(),
                                                                      p_face_lateral));

                map_lateral_faces[element_index].push_back(
                    std::pair<unsigned, bool>(index_lateral_face, points_inwards));
            }
            // Go to next edge
            node_a_index = node_b_index;
        }
    }

    // Now we can create the 3d elements with the faces.

    for (auto mapping = map_lateral_faces.begin();
         mapping != map_lateral_faces.end(); ++mapping)
    {
        unsigned element_index = mapping->first;
        std::vector<std::pair<unsigned, bool> > faces_ind_orient = mapping->second;
        std::vector<MonolayerVertexElement<2, 3>*> lateral_faces_in_element;
        std::vector<bool> lateral_orientations_in_element;

        for (auto fc = faces_ind_orient.begin(); fc != faces_ind_orient.end();
             ++fc)
        {
            lateral_faces_in_element.push_back(faces_lateral[fc->first]);
            lateral_orientations_in_element.push_back(fc->second);
        }

        // Now create the element with the correct faces

        std::vector<MonolayerVertexElement<2, 3>*> faces_in_element;
        faces_in_element.push_back(faces_apical[element_index]);
        faces_in_element.push_back(faces_basal[element_index]);
        faces_in_element.insert(faces_in_element.end(),
                                lateral_faces_in_element.begin(),
                                lateral_faces_in_element.end());

        std::vector<bool> face_orientations_in_element;
        face_orientations_in_element.push_back(orientations_apical[element_index]);
        face_orientations_in_element.push_back(orientations_basal[element_index]);
        face_orientations_in_element.insert(face_orientations_in_element.end(),
                                            lateral_orientations_in_element.begin(),
                                            lateral_orientations_in_element.end());

        elements.push_back(new MonolayerVertexElement<3, 3>(
            element_index, MonolayerVertexElementType::Undetermined,
            faces_in_element, face_orientations_in_element));
    }

    // Now create the mesh
    std::vector<Node<3>*> nodes = nodes_apical;
    nodes.insert(nodes.end(), nodes_basal.begin(), nodes_basal.end());

    mpMesh = new MutableMonolayerVertexMesh<3, 3>(
        nodes, elements, cellRearrangementThreshold, t2Threshold);
}

c_vector<double, 3>
FiniteThicknessRandomizedSphereMeshGenerator::GetNormalToFace(
    MonolayerVertexElement<2, 3>* pFace)
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
FiniteThicknessRandomizedSphereMeshGenerator::GetCenterOfFace(
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

FiniteThicknessRandomizedSphereMeshGenerator::
    ~FiniteThicknessRandomizedSphereMeshGenerator()
{
    delete mpMesh;
}

MutableMonolayerVertexMesh<3, 3>* FiniteThicknessRandomizedSphereMeshGenerator::GetMesh()
{
    return mpMesh;
}

bool FiniteThicknessRandomizedSphereMeshGenerator::DoesOverlap(c_vector<double, 3> new_cell_center, std::vector<c_vector<double, 3> >& cell_centers, double min_distance)
{
    for (c_vector<double, 3> existing_cell_center : cell_centers)
    {
        double distance = norm_2(new_cell_center - existing_cell_center);
        if (distance < min_distance)
        {
            return true;
        }
    }
    return false;
}

c_vector<double, 3> FiniteThicknessRandomizedSphereMeshGenerator::GenerateRandomPointOnSphere()
{
    RandomNumberGenerator* p_random = RandomNumberGenerator::Instance();

    double rand_u1 = p_random->ranf();
    double rand_u2 = p_random->ranf();

    // We get a uniform distribution on the sphere with Inverse Transform
    // Sampling
    double rand_phi = 2.0 * M_PI * rand_u1;
    double rand_theta = acos(1.0 - 2.0 * rand_u2);
    c_vector<double, 3> rand_vector;
    rand_vector <<= mInnerRadius * sin(rand_theta) * cos(rand_phi),
        mInnerRadius * sin(rand_theta) * sin(rand_phi),
        mInnerRadius * cos(rand_theta);

    return rand_vector;
}