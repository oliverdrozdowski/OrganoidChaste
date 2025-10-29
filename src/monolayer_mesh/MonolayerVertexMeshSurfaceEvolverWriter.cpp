#include "MonolayerVertexMeshSurfaceEvolverWriter.hpp"
#include "ChasteSerialization.hpp"
#include "RandomNumberGenerator.hpp"
#include "Version.hpp"

/**
 * Convenience collection of iterators, primarily to get compilation to happen.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct MonolayerMeshSurfaceEvolverWriterIterators
{
    /** Iterator over nodes */
    typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator* pNodeIter;
    /** Iterator over vertex elements */
    typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator* pElemIter;
    /** Iterator over faces */
    typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator* pFaceIter;
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexMeshSurfaceEvolverWriter(const std::string& rDirectory,
                                                                                                         const std::string& rBaseName,
                                                                                                         const bool clearOutputDir)
        : AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName, clearOutputDir),
          mpMesh(nullptr),
          mpIters(new MonolayerMeshSurfaceEvolverWriterIterators<ELEMENT_DIM, SPACE_DIM>),
          mpNodeMap(nullptr),
          mNodeMapCurrentIndex(0),
          mUseRandomizedVolumes(false)
{
    mpIters->pNodeIter = nullptr;
    mpIters->pElemIter = nullptr;
    mpIters->pFaceIter = nullptr;
    mpSurfaceTensionSubForce = nullptr;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::~MonolayerVertexMeshSurfaceEvolverWriter()
{
    if (mpIters->pNodeIter)
    {
        delete mpIters->pNodeIter;
        delete mpIters->pElemIter;
        delete mpIters->pFaceIter;
    }

    delete mpIters;

    if (mpNodeMap)
    {
        delete mpNodeMap;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::SetSurfaceTensionSubForce(boost::shared_ptr<SurfaceTensionSubForce<SPACE_DIM> > p_force)
{
    mpSurfaceTensionSubForce = p_force;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::SetFixBoundaryNodes(bool fixBoundaries)
{
    mFixBoundaryNodes = fixBoundaries;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::SetWriteFaceTypeIntoFile(bool writeFaceType)
{
    mWriteFaceTypeIntoFile = writeFaceType;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::SetConsiderLumenAsCell(bool considerLumenAsCell)
{
    mConsiderLumenAsCell = considerLumenAsCell;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::SetUseRandomizedVolumes(bool useRandomVol, double distributionBoundary)
{
    mUseRandomizedVolumes = useRandomVol;
    mUniformDistributionBoundary = distributionBoundary;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::GetUseRandomizedVolumes()
{
    return mUseRandomizedVolumes;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
double MonolayerVertexMeshSurfaceEvolverWriter<1, 1>::GetSurfaceTensionOfFace(unsigned index)
{
    return 0.0;
}

template <>
double MonolayerVertexMeshSurfaceEvolverWriter<1, 2>::GetSurfaceTensionOfFace(unsigned index)
{
    return 0.0;
}

template <>
double MonolayerVertexMeshSurfaceEvolverWriter<1, 3>::GetSurfaceTensionOfFace(unsigned index)
{
    return 0.0;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::GetSurfaceTensionOfFace(unsigned index)
{
    if (mpSurfaceTensionSubForce)
    {
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = mpMesh->GetFace(index);
        bool boundary_face = p_face->IsBoundaryFace();
        MonolayerVertexElementType type = p_face->GetFaceType();
        double boundary_factor = boundary_face ? 0.5 : 1.0;
        if (type == MonolayerVertexElementType::Lateral)
        {
            // Laterals are saved as half their tensions, except for boundaries, which are saved
            // as full tensions, as they are only looped over once.
            return boundary_factor * mpSurfaceTensionSubForce->GetSurfaceTensionParameter(index) * 2.0;
        }
        else
        {
            return mpSurfaceTensionSubForce->GetSurfaceTensionParameter(index);
        }
    }
    else
    {
        return 1.0;
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
unsigned MonolayerVertexMeshSurfaceEvolverWriter<1, 1>::GetFaceTypeOfFace(unsigned index)
{
    return 0;
}

template <>
unsigned MonolayerVertexMeshSurfaceEvolverWriter<1, 2>::GetFaceTypeOfFace(unsigned index)
{
    return 0;
}

template <>
unsigned MonolayerVertexMeshSurfaceEvolverWriter<1, 3>::GetFaceTypeOfFace(unsigned index)
{
    return 0;
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::GetFaceTypeOfFace(unsigned index)
{
    if (mWriteFaceTypeIntoFile)
    {
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = mpMesh->GetFace(index);
        MonolayerVertexElementType type = p_face->GetFaceType();
        return (unsigned)type;
    }
    else
    {
        return 0;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::SetMapTensionToColor(std::map<double, int> map)
{
    mMapTensionToColor = map;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
int MonolayerVertexMeshSurfaceEvolverWriter<1, 1>::GetColorOfFace(unsigned index)
{
    EXCEPTION("The Surface Evolver Writer only works in 3D.");
}

template <>
int MonolayerVertexMeshSurfaceEvolverWriter<1, 2>::GetColorOfFace(unsigned index)
{
    EXCEPTION("The Surface Evolver Writer only works in 3D.");
}

template <>
int MonolayerVertexMeshSurfaceEvolverWriter<1, 3>::GetColorOfFace(unsigned index)
{
    EXCEPTION("The Surface Evolver Writer only works in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
int MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::GetColorOfFace(unsigned index)
{
    if (mMapTensionToColor.empty())
    {
        return (int)mpMesh->GetFace(index)->GetFaceType() + 1;
    }
    else
    {
        double surface_tension = GetSurfaceTensionOfFace(index);

        const auto findr = std::find_if(mMapTensionToColor.begin(), mMapTensionToColor.end(), [&surface_tension](const auto& findr)
                                        {
                                            return abs(findr.first - surface_tension) < 1e-6; // Comparing with the object
                                        });
        if (findr == mMapTensionToColor.end())
        {
            return 0;
        }
        return findr->second;
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    if (mpMesh)
    {
        // Sanity check
        assert(this->mNumNodes == mpMesh->GetNumNodes());

        std::vector<double> coordinates(SPACE_DIM);

        // Get the node coordinates using the node iterator (thus skipping deleted nodes)
        for (unsigned j = 0; j < SPACE_DIM; j++)
        {
            coordinates[j] = (*(mpIters->pNodeIter))->GetPoint()[j];
        }

        ++(*(mpIters->pNodeIter));

        return coordinates;
    }
    else
    {
        return AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextNode();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::IsNodeOnBoundary()
{
    return (*(mpIters->pNodeIter))->IsBoundaryNode();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::CreateEdgeMapping()
{
    mEdgeIndexPairs.clear();
    mMapEdgeToIndex.clear();

    for (auto face_iterator = mpMesh->GetFaceIteratorBegin(); face_iterator != mpMesh->GetFaceIteratorEnd(); ++face_iterator)
    {
        for (unsigned face_node_index = 0; face_node_index < face_iterator->GetNumNodes(); ++face_node_index)
        {
            unsigned first_index = face_iterator->GetNode(face_node_index)->GetIndex();
            unsigned second_index = face_iterator->GetNode((face_node_index + 1) % face_iterator->GetNumNodes())->GetIndex();
            std::array<unsigned, 2> edge = { first_index, second_index };
            std::array<unsigned, 2> inverse_edge = { second_index, first_index };
            bool edge_known = mMapEdgeToIndex.find(edge) != mMapEdgeToIndex.end();
            bool inverse_edge_known = mMapEdgeToIndex.find(inverse_edge) != mMapEdgeToIndex.end();
            if (!edge_known && !inverse_edge_known)
            {
                mEdgeIndexPairs.push_back(edge);
                mMapEdgeToIndex[edge] = mEdgeIndexPairs.size() - 1;
                mMapFaceIndexToVectorOfEdges[face_iterator->GetIndex()].push_back(mEdgeIndexPairs.size()); // implicit +1!

                // Now determine whether boundary nodes
                bool first_bndry = face_iterator->GetNode(face_node_index)->IsBoundaryNode();
                bool second_bndry = face_iterator->GetNode((face_node_index + 1) % face_iterator->GetNumNodes())->IsBoundaryNode();
                mEdgeOnBoundary.push_back(first_bndry && second_bndry);
            }
            else if (edge_known && !inverse_edge_known)
            {
                unsigned index_edge = mMapEdgeToIndex[edge];
                // We have to add the +1, because otherwise we have problems handling +0 and -0 later on!
                mMapFaceIndexToVectorOfEdges[face_iterator->GetIndex()].push_back((int)index_edge + 1);
            }
            else if (!edge_known && inverse_edge_known)
            {
                unsigned index_edge = mMapEdgeToIndex[inverse_edge];
                mMapFaceIndexToVectorOfEdges[face_iterator->GetIndex()].push_back(-(int)index_edge - 1);
            }
            else
            {
                // This can't happen
                NEVER_REACHED;
            }
        }
    }
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MonolayerVertexMeshSurfaceEvolverWriter<1, 1>::WriteFilesUsingMesh(MonolayerVertexMesh<1, 1>& rMesh)
{
    EXCEPTION("The Surface Evolver Writer only works in 3D.");
}

template <>
void MonolayerVertexMeshSurfaceEvolverWriter<1, 2>::WriteFilesUsingMesh(MonolayerVertexMesh<1, 2>& rMesh)
{
    EXCEPTION("The Surface Evolver Writer only works in 3D.");
}

template <>
void MonolayerVertexMeshSurfaceEvolverWriter<1, 3>::WriteFilesUsingMesh(MonolayerVertexMesh<1, 3>& rMesh)
{
    EXCEPTION("The Surface Evolver Writer only works in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

///\todo Mesh should be const (#1076)
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
    this->mpMeshReader = nullptr;
    mpMesh = &rMesh;

    this->mNumNodes = mpMesh->GetNumNodes();
    this->mNumFaces = mpMesh->GetNumFaces();
    this->mNumElements = mpMesh->GetNumElements();

    typedef typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator NodeIterType;
    mpIters->pNodeIter = new NodeIterType(mpMesh->GetNodeIteratorBegin());

    typedef typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator ElemIterType;
    mpIters->pElemIter = new ElemIterType(mpMesh->GetElementIteratorBegin());

    typedef typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator FaceIterType;
    mpIters->pFaceIter = new FaceIterType(mpMesh->GetFaceIteratorBegin());

    // Set up node map if we might have deleted nodes
    mNodeMapCurrentIndex = 0;
    if (mpMesh->IsMeshChanging())
    {
        mpNodeMap = new NodeMap(mpMesh->GetNumAllNodes());
        for (NodeIterType it = mpMesh->GetNodeIteratorBegin(); it != mpMesh->GetNodeIteratorEnd(); ++it)
        {
            mpNodeMap->SetNewIndex(it->GetIndex(), mNodeMapCurrentIndex++);
        }
    }
    CreateEdgeMapping();
    WriteFiles();
}
/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
void MonolayerVertexMeshSurfaceEvolverWriter<1, 1>::WriteFiles()
{
    EXCEPTION("The Surface Evolver Writer only works in 3D.");
}

template <>
void MonolayerVertexMeshSurfaceEvolverWriter<1, 2>::WriteFiles()
{
    EXCEPTION("The Surface Evolver Writer only works in 3D.");
}

template <>
void MonolayerVertexMeshSurfaceEvolverWriter<1, 3>::WriteFiles()
{
    EXCEPTION("The Surface Evolver Writer only works in 3D.");
}
/// \endcond Get Doxygen to ignore, since it's confused by these templates

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshSurfaceEvolverWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
    std::string comment = "// Surface Evolver starting conditions created with Chaste \n//" + ChasteBuildInfo::GetProvenanceString();

    // Write file
    std::string file_name = this->mBaseName + ".fe";
    out_stream p_file = this->mpOutputFileHandler->OpenOutputFile(file_name);

    // This only works in 3D. Exception throwing would be nice...
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE

    unsigned num_nodes = this->GetNumNodes();

    *p_file << comment << "\n\n";
    *p_file << std::setprecision(6);

    if (mWriteFaceTypeIntoFile)
    {
        *p_file << "define facet attribute facetype integer\n"
                << "define facet attribute faceid integer\n\n";
    }

    *p_file << "vertices\n";
    // Write each node's data
    for (unsigned item_num = 0; item_num < num_nodes; item_num++)
    {
        // First read boundary status, as GetNextNode updates at end of call
        bool boundary_node = IsNodeOnBoundary();
        std::vector<double> current_item = this->GetNextNode();
        *p_file << item_num + 1;
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            *p_file << "\t" << current_item[i];
        }
        if (mFixBoundaryNodes)
        {
            if (boundary_node)
                *p_file << "\t"
                        << "fixed";
        }
        *p_file << "\n";
    }
    *p_file << "\n";

    // Now write each edge's data
    *p_file << "edges\n";
    unsigned edge_num = 1;
    for (auto edge_it = mEdgeIndexPairs.begin(); edge_it != mEdgeIndexPairs.end(); edge_it++)
    {
        *p_file << edge_num << "\t";
        *p_file << (*edge_it)[0] + 1 << "\t";
        *p_file << (*edge_it)[1] + 1;
        if (mFixBoundaryNodes)
        {
            if (mEdgeOnBoundary[edge_num - 1])
                *p_file << "\t"
                        << "fixed";
        }
        *p_file << "\n";
        edge_num++;
    }
    *p_file << "\n";

    // Now write each face's data
    *p_file << "faces\n";
    for (unsigned index_face = 0; index_face < mNumFaces; index_face++)
    {
        *p_file << index_face + 1;
        std::vector<int> edges = mMapFaceIndexToVectorOfEdges[index_face];
        for (auto it = edges.begin(); it != edges.end(); it++)
        {
            *p_file << "\t" << *it; // We already added the +1 in CreateEdgeMapping!
        }
        if (mFixBoundaryNodes)
        {
            if (mpMesh->GetFace(index_face)->IsBoundaryFace())
                *p_file << "\t"
                        << "fixed";
        }
        *p_file << "\t"
                << "tension"
                << "\t";
        *p_file << GetSurfaceTensionOfFace(index_face);
        *p_file << "\t"
                << "color"
                << "\t";
        *p_file << GetColorOfFace(index_face);
        if (mWriteFaceTypeIntoFile)
        {
            *p_file << "\t"
                    << "facetype"
                    << "\t";
            *p_file << GetFaceTypeOfFace(index_face);
        }
        // Write the faceid
        *p_file << "\t"
                << "faceid"
                << "\t";
        *p_file << index_face;
        *p_file << "\n";
    }
    *p_file << "\n";

    // Now write each cell's data
    RandomNumberGenerator* p_random = RandomNumberGenerator::Instance();
    *p_file << "bodies\n";
    for (unsigned index_cell = 0; index_cell < this->mNumElements; index_cell++)
    {
        *p_file << index_cell + 1;
        MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>* p_element = mpMesh->GetElement(index_cell);
        unsigned num_faces_in_cell = p_element->GetNumFaces();
        for (unsigned face_index = 0; face_index < num_faces_in_cell; face_index++)
        {
            // Orientation true means clockwise from outside, i.e. negatively oriented
            if (p_element->FaceIsOrientatedClockwise(face_index))
            {
                *p_file << "\t"
                        << "-" << (p_element->GetFace(face_index)->GetIndex() + 1);
            }
            else
            {
                *p_file << "\t" << (p_element->GetFace(face_index)->GetIndex() + 1);
            }
        }
        *p_file << "\t"
                << "volume";
        if (mUseRandomizedVolumes)
        {
            double rand_uniform = p_random->ranf();
            rand_uniform *= mUniformDistributionBoundary;
            *p_file << "\t" << mpMesh->GetVolumeOfElement(index_cell) + rand_uniform;
        }
        else
        {
            *p_file << "\t" << mpMesh->GetVolumeOfElement(index_cell);
        }
        *p_file << "\n";
    }
    // If we want to consider the lumen as a cell we add it at the end as another cell
    // n.b. We do NOT fix the volume here!
    if (mConsiderLumenAsCell)
    {
        *p_file << this->mNumElements + 1;
        unsigned num_faces = mpMesh->GetNumFaces();
        for (unsigned ind = 0; ind < num_faces; ind++)
        {
            // Add face if apical. We assume apical faces to always be oriented inwards-pointing!
            if (mpMesh->GetFace(ind)->GetFaceType() == MonolayerVertexElementType::Apical)
                *p_file << "\t"
                        << "-" << (ind + 1);
        }
    }
    *p_file << "\n\n";
    p_file->close();
}

///////// Explicit instantiation///////

template class MonolayerVertexMeshSurfaceEvolverWriter<1, 1>;
template class MonolayerVertexMeshSurfaceEvolverWriter<1, 2>;
template class MonolayerVertexMeshSurfaceEvolverWriter<1, 3>;
template class MonolayerVertexMeshSurfaceEvolverWriter<2, 2>;
template class MonolayerVertexMeshSurfaceEvolverWriter<2, 3>;
template class MonolayerVertexMeshSurfaceEvolverWriter<3, 3>;