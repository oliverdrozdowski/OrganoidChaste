/*
Copyright (c) 2005-2019, University of Oxford.
Copyright (c) 2025, Oliver M. Drozdowski and Ulrich S. Schwarz (Heidelberg University)

All rights reserved.
University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.
This file is part of Chaste.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "MonolayerVertexMeshWriter.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "Toroidal2dVertexMesh.hpp"
#include "Version.hpp"

/**
 * Convenience collection of iterators, primarily to get compilation to happen.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct MonolayerMeshWriterIterators
{
    /** Iterator over nodes */
    typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator* pNodeIter;
    /** Iterator over vertex elements */
    typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator* pElemIter;
    /** Iterator over vertex faces */
    typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator* pFaceIter;
};

///////////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexMeshWriter(
    const std::string& rDirectory, const std::string& rBaseName,
    const bool clearOutputDir)
        : AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName,
                                                     clearOutputDir),
          mpMesh(nullptr),
          mpIters(new MonolayerMeshWriterIterators<ELEMENT_DIM, SPACE_DIM>),
          mpNodeMap(nullptr), mNodeMapCurrentIndex(0), mpFaceMap(nullptr), mFaceMapCurrentIndex(0)
{
    mpIters->pNodeIter = nullptr;
    mpIters->pElemIter = nullptr;
    mpIters->pFaceIter = nullptr;

#ifdef CHASTE_VTK
    // Dubious, since we shouldn't yet know what any details of the mesh are.
    mpVtkUnstructedMesh = vtkUnstructuredGrid::New();
    mpVtkUnstructedMeshTriangulation = vtkUnstructuredGrid::New();
#endif // CHASTE_VTK
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::~MonolayerVertexMeshWriter()
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

    if (mpFaceMap)
    {
        delete mpFaceMap;
    }

#ifdef CHASTE_VTK
    // Dubious, since we shouldn't yet know what any details of the mesh are.
    mpVtkUnstructedMesh->Delete(); // Reference counted
    mpVtkUnstructedMeshTriangulation->Delete();
#endif // CHASTE_VTK
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextNode()
{
    if (mpMesh)
    {
        // Sanity check
        assert(this->mNumNodes == mpMesh->GetNumNodes());

        std::vector<double> coordinates(SPACE_DIM + 1);

        // Get the node coordinates using the node iterator (thus skipping deleted nodes)
        for (unsigned j = 0; j < SPACE_DIM; j++)
        {
            coordinates[j] = (*(mpIters->pNodeIter))->GetPoint()[j];
        }
        coordinates[SPACE_DIM] = (*(mpIters->pNodeIter))->IsBoundaryNode();

        ++(*(mpIters->pNodeIter));

        return coordinates;
    }
    else
    {
        return AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextNode();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElementData MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextElementWithFaces()
{
    // This method should only be called in 3D
    assert(SPACE_DIM == 3); // LCOV_EXCL_LINE
    assert(mpMesh);
    assert(this->mNumElements == mpMesh->GetNumElements());

    // Create data structure for this element
    VertexElementData elem_data;

    // Store node indices owned by this element
    elem_data.NodeIndices.resize((*(mpIters->pElemIter))->GetNumNodes());
    for (unsigned j = 0; j < elem_data.NodeIndices.size(); j++)
    {
        unsigned old_index = (*(mpIters->pElemIter))->GetNodeGlobalIndex(j);
        elem_data.NodeIndices[j] = mpMesh->IsMeshChanging()
            ? mpNodeMap->GetNewIndex(old_index)
            : old_index;
    }

    // Store faces owned by this element
    elem_data.Faces.resize((*(mpIters->pElemIter))->GetNumFaces());
    for (unsigned i = 0; i < elem_data.Faces.size(); i++)
    {
        // Get pointer to this face
        MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* p_face = (*(mpIters->pElemIter))->GetFace(i);

        // Create data structure for this face
        ElementData face_data;

        // Store this face's index
        face_data.AttributeValue = p_face->GetIndex();

        // Store node indices owned by this face
        face_data.NodeIndices.resize(p_face->GetNumNodes());
        for (unsigned j = 0; j < face_data.NodeIndices.size(); j++)
        {
            unsigned old_index = p_face->GetNodeGlobalIndex(j);
            face_data.NodeIndices[j] = mpMesh->IsMeshChanging()
                ? mpNodeMap->GetNewIndex(old_index)
                : old_index;
        }

        // Store this face
        elem_data.Faces[i] = face_data;

        ///\todo Store face orientations? (#1076/#1377)
    }

    // Set attribute
    elem_data.AttributeValue = (unsigned)(*(mpIters->pElemIter))->GetAttribute();
    ++(*(mpIters->pElemIter));

    return elem_data;
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
MonolayerVertexFaceData MonolayerVertexMeshWriter<1, 1>::GetNextMonolayerFace()
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMeshWriter<1,1>::GetNextMonolayerFace() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
MonolayerVertexFaceData MonolayerVertexMeshWriter<1, 2>::GetNextMonolayerFace()
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMeshWriter<1,2>::GetNextMonolayerFace() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
MonolayerVertexFaceData MonolayerVertexMeshWriter<1, 3>::GetNextMonolayerFace()
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMeshWriter<1,3>::GetNextMonolayerFace() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
MonolayerVertexFaceData MonolayerVertexMeshWriter<2, 2>::GetNextMonolayerFace()
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMeshWriter<2,2>::GetNextMonolayerFace() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
MonolayerVertexFaceData MonolayerVertexMeshWriter<2, 3>::GetNextMonolayerFace()
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    EXCEPTION("MonolayerVertexMeshWriter<2,3>::GetNextMonolayerFace() is not implemented");
}

/// \cond Get Doxygen to ignore, since it's confused by these templates
template <>
MonolayerVertexFaceData MonolayerVertexMeshWriter<3, 3>::GetNextMonolayerFace()
/// \endcond Get Doxygen to ignore, since it's confused by these templates
{
    assert(mpMesh != nullptr);
    assert(this->mNumFaces == mpMesh->GetNumFaces());

    MonolayerVertexFaceData face_data;
    face_data.NodeIndices.resize((*(mpIters->pFaceIter))->GetNumNodes());
    face_data.NodeTypes.resize((*(mpIters->pFaceIter))->GetNumNodes());

    // Set FaceType
    face_data.FaceType = static_cast<unsigned>((*(mpIters->pFaceIter))->GetFaceType());

    for (unsigned j = 0; j < face_data.NodeIndices.size(); j++)
    {
        unsigned old_index = (*(mpIters->pFaceIter))->GetNodeGlobalIndex(j);
        face_data.NodeIndices[j] = mpMesh->IsMeshChanging()
            ? mpNodeMap->GetNewIndex(old_index)
            : old_index;
        face_data.NodeTypes[j] = (face_data.FaceType < 3) ? face_data.FaceType : static_cast<unsigned>((*(mpIters->pFaceIter))->GetNodeType(j));
    }

    // Set attribute (no use for attributes as of now)
    // face_data.AttributeValue = (*(mpIters->pFaceIter))->GetAttribute();
    ++(*(mpIters->pFaceIter));

    return face_data;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
MonolayerVertexElementData
MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextMonolayerElement()
{
    assert(mpMesh != nullptr);

    assert(this->mNumElements == mpMesh->GetNumElements());

    MonolayerVertexElementData elem_data;
    elem_data.FaceIndices.resize((*(mpIters->pElemIter))->GetNumFaces());
    elem_data.FaceOrientations.resize((*(mpIters->pElemIter))->GetNumFaces());

    for (unsigned j = 0; j < elem_data.FaceIndices.size(); j++)
    {
        // Get face index
        unsigned old_index = (*(mpIters->pElemIter))->GetFace(j)->GetIndex();
        elem_data.FaceIndices[j] = mpMesh->IsMeshChanging()
            ? mpFaceMap->GetNewIndex(old_index)
            : old_index;
        elem_data.FaceOrientations[j] = (*(mpIters->pElemIter))->FaceIsOrientatedClockwise(j);
    }

    // Set attribute
    // elem_data.AttributeValue = (*(mpIters->pElemIter))->GetAttribute();
    ++(*(mpIters->pElemIter));

    return elem_data;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ElementData MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextElement()
{
    ///\todo Assert this method should only be called in 2D? (#1076/#1377)

    if (mpMesh)
    {
        assert(this->mNumElements == mpMesh->GetNumElements());

        ElementData elem_data;
        elem_data.NodeIndices.resize((*(mpIters->pElemIter))->GetNumNodes());
        for (unsigned j = 0; j < elem_data.NodeIndices.size(); j++)
        {
            unsigned old_index = (*(mpIters->pElemIter))->GetNodeGlobalIndex(j);
            elem_data.NodeIndices[j] = mpMesh->IsMeshChanging() ? mpNodeMap->GetNewIndex(old_index) : old_index;
        }

        // Set attribute
        elem_data.AttributeValue = (*(mpIters->pElemIter))->GetAttribute();
        ++(*(mpIters->pElemIter));

        return elem_data;
    }
    else
    {
        return AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>::GetNextElement();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteVtkUsingMesh(MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, std::string stamp)
{
#ifdef CHASTE_VTK
    assert(SPACE_DIM == 3 || SPACE_DIM == 2); // LCOV_EXCL_LINE

    // Create VTK mesh
    // MakeVtkMesh(rMesh);
    MakeVtkMeshTriangulatedVolume(rMesh);
    // Create triangulated mesh
    MakeVtkMeshTriangulated(rMesh);

    // Now write VTK mesh to file
    assert(mpVtkUnstructedMesh->CheckAttributes() == 0);
    assert(mpVtkUnstructedMeshTriangulation->CheckAttributes() == 0);
    vtkXMLUnstructuredGridWriter* p_writer = vtkXMLUnstructuredGridWriter::New();
    vtkXMLUnstructuredGridWriter* p_writer_triangulated = vtkXMLUnstructuredGridWriter::New();
#if VTK_MAJOR_VERSION >= 6
    p_writer->SetInputData(mpVtkUnstructedMesh);
    p_writer_triangulated->SetInputData(mpVtkUnstructedMeshTriangulation);
#else
    p_writer->SetInput(mpVtkUnstructedMesh);
    p_writer_triangulated->SetInput(mpVtkUnstructedMeshTriangulation);
#endif
    // Uninitialised stuff arises (see #1079), but you can remove valgrind
    // problems by removing compression:
    // **** REMOVE WITH CAUTION *****
    // p_writer->SetCompressor(nullptr);
    p_writer_triangulated->SetCompressor(nullptr);
    // **** REMOVE WITH CAUTION *****

    std::string vtk_file_name = this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mBaseName;
    std::string vtk_file_name_triangulated = this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mBaseName + "_tri";
    if (stamp != "")
    {
        vtk_file_name += "_" + stamp;
        vtk_file_name_triangulated += "_" + stamp;
    }
    vtk_file_name += ".vtu";
    vtk_file_name_triangulated += ".vtu";

    p_writer->SetFileName(vtk_file_name.c_str());
    p_writer_triangulated->SetFileName(vtk_file_name_triangulated.c_str());
    // p_writer->PrintSelf(std::cout, vtkIndent());
    p_writer->SetDataModeToAscii();
    p_writer->Update();
    p_writer->Write();
    p_writer->Delete(); // Reference counted

    // p_writer_triangulated->PrintSelf(std::cout, vtkIndent());
    p_writer_triangulated->Write();
    p_writer_triangulated->Delete(); // Reference counted
#endif // CHASTE_VTK
}

/**
 * Write VTK file using a mesh.
 *
 * @param rMesh reference to the vertex-based mesh
 * @param stamp is an optional stamp (like a time-stamp) to put into the name of
 * the file
 */
template <>
void MonolayerVertexMeshWriter<2, 2>::WriteVtkUsingMesh(MonolayerVertexMesh<2, 2>& rMesh, std::string stamp)
{
#ifdef CHASTE_VTK
    // Create VTK mesh
    MonolayerVertexMesh<2, 2>* p_mesh_for_vtk = rMesh.GetMeshForVtk();
    MakeVtkMesh(*p_mesh_for_vtk);

    // Now write VTK mesh to file
    assert(mpVtkUnstructedMesh->CheckAttributes() == 0);
    vtkXMLUnstructuredGridWriter* p_writer = vtkXMLUnstructuredGridWriter::New();
#if VTK_MAJOR_VERSION >= 6
    p_writer->SetInputData(mpVtkUnstructedMesh);
#else
    p_writer->SetInput(mpVtkUnstructedMesh);
#endif
    // Uninitialised stuff arises (see #1079), but you can remove valgrind
    // problems by removing compression:
    // **** REMOVE WITH CAUTION *****
    p_writer->SetCompressor(nullptr);
    // **** REMOVE WITH CAUTION *****

    std::string vtk_file_name = this->mpOutputFileHandler->GetOutputDirectoryFullPath() + this->mBaseName;
    if (stamp != "")
    {
        vtk_file_name += "_" + stamp;
    }
    vtk_file_name += ".vtu";

    p_writer->SetFileName(vtk_file_name.c_str());
    // p_writer->PrintSelf(std::cout, vtkIndent());
    p_writer->Write();
    p_writer->Delete(); // Reference counted
#endif // CHASTE_VTK
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::MakeVtkMesh(MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
#ifdef CHASTE_VTK
    // Make the Vtk mesh
    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");
    for (unsigned node_num = 0; node_num < rMesh.GetNumNodes(); node_num++)
    {
        c_vector<double, SPACE_DIM> position;
        position = rMesh.GetNode(node_num)->rGetLocation();
        if (SPACE_DIM == 2)
        {
            p_pts->InsertPoint(node_num, position[0], position[1], 0.0);
        }
        else
        {
            p_pts->InsertPoint(node_num, position[0], position[1], position[2]);
        }
    }

    mpVtkUnstructedMesh->SetPoints(p_pts);
    p_pts->Delete(); // Reference counted
    for (typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator iter = rMesh.GetElementIteratorBegin();
         iter != rMesh.GetElementIteratorEnd();
         ++iter)
    {
        vtkCell* p_cell;
        if (SPACE_DIM == 2)
        {
            p_cell = vtkPolygon::New();
        }
        else
        {
            p_cell = vtkConvexPointSet::New();
        }
        vtkIdList* p_cell_id_list = p_cell->GetPointIds();
        p_cell_id_list->SetNumberOfIds(iter->GetNumNodes());
        for (unsigned j = 0; j < iter->GetNumNodes(); ++j)
        {
            p_cell_id_list->SetId(j, iter->GetNodeGlobalIndex(j));
        }
        mpVtkUnstructedMesh->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
        p_cell->Delete(); // Reference counted
    }
#endif // CHASTE_VTK
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::MakeVtkMeshTriangulated(MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
#ifdef CHASTE_VTK
    // Make the Vtk mesh
    // Create the points
    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");
    for (unsigned node_num = 0; node_num < rMesh.GetNumNodes() + rMesh.GetNumFaces(); node_num++)
    {
        if (node_num < rMesh.GetNumNodes())
        {
            c_vector<double, SPACE_DIM> position;
            position = rMesh.GetNode(node_num)->rGetLocation();
            if (SPACE_DIM == 2)
            {
                p_pts->InsertPoint(node_num, position[0], position[1], 0.0);
            }
            else
            {
                p_pts->InsertPoint(node_num, position[0], position[1], position[2]);
            }
        }
        else
        {
            c_vector<double, SPACE_DIM> position;
            position = rMesh.GetPassiveCenterOfFace(rMesh.GetFace(node_num - rMesh.GetNumNodes()));
            if (SPACE_DIM == 2)
            {
                p_pts->InsertPoint(node_num, position[0], position[1], 0.0);
            }
            else
            {
                p_pts->InsertPoint(node_num, position[0], position[1], position[2]);
            }
        }
    }

    mpVtkUnstructedMeshTriangulation->SetPoints(p_pts);
    p_pts->Delete(); // Reference counted

    // Create faces in mesh as cells
    unsigned num_face = 0;
    for (typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator iter = rMesh.GetFaceIteratorBegin();
         iter != rMesh.GetFaceIteratorEnd();
         ++iter)
    {
        // Create a 2d cell for every triangle (#triangles = #nodes)
        for (unsigned num_triangle = 0; num_triangle < iter->GetNumNodes(); ++num_triangle)
        {
            vtkCell* p_cell;
            if (SPACE_DIM == 2)
            {
                p_cell = vtkTriangle::New();
            }
            else
            {
                p_cell = vtkTriangle::New();
            }
            vtkIdList* p_cell_id_list = p_cell->GetPointIds();
            p_cell_id_list->SetNumberOfIds(3);

            // passive center
            p_cell_id_list->SetId(0, rMesh.GetNumNodes() + num_face);
            // counter clockwise boundary nodes
            p_cell_id_list->SetId(1, iter->GetNodeGlobalIndex((num_triangle + 1) % iter->GetNumNodes()));
            p_cell_id_list->SetId(2, iter->GetNodeGlobalIndex((num_triangle + 0) % iter->GetNumNodes()));

            mpVtkUnstructedMeshTriangulation->InsertNextCell(p_cell->GetCellType(), p_cell_id_list);
            p_cell->Delete(); // Reference counted
        }
        ++num_face;
    }
#endif // CHASTE_VTK
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::MakeVtkMeshTriangulatedVolume(MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
#ifdef CHASTE_VTK
    // Make the Vtk mesh
    // Create the points

    vtkPoints* p_pts = vtkPoints::New(VTK_DOUBLE);
    p_pts->GetData()->SetName("Vertex positions");
    for (unsigned node_num = 0; node_num < rMesh.GetNumNodes() + rMesh.GetNumFaces(); node_num++)
    {
        if (node_num < rMesh.GetNumNodes())
        {
            c_vector<double, SPACE_DIM> position;
            position = rMesh.GetNode(node_num)->rGetLocation();
            if (SPACE_DIM == 2)
            {
                p_pts->InsertPoint(node_num, position[0], position[1], 0.0);
            }
            else
            {
                p_pts->InsertPoint(node_num, position[0], position[1], position[2]);
            }
        }
        else
        {
            c_vector<double, SPACE_DIM> position;
            position = rMesh.GetPassiveCenterOfFace(rMesh.GetFace(node_num - rMesh.GetNumNodes()));
            if (SPACE_DIM == 2)
            {
                p_pts->InsertPoint(node_num, position[0], position[1], 0.0);
            }
            else
            {
                p_pts->InsertPoint(node_num, position[0], position[1], position[2]);
            }
        }
    }

    mpVtkUnstructedMesh->SetPoints(p_pts);
    p_pts->Delete(); // Reference counted

    // Map to indices of faces when iterating through elements rather than faces
    std::map<MonolayerVertexElement<ELEMENT_DIM, SPACE_DIM>*, std::vector<unsigned> >* pElementsFacesMap;
    pElementsFacesMap = rMesh.GetElementsFacesMap();

    // We iterate through the elements
    for (typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator iter = rMesh.GetElementIteratorBegin();
         iter != rMesh.GetElementIteratorEnd();
         ++iter)
    {
        // Create face stream for polyhedron
        std::vector<unsigned> temporary_face_stream;
        vtkIdType num_faces = 0;
        std::vector<unsigned> face_index_vector = (*pElementsFacesMap)[&*iter];
        if (face_index_vector.size() == 0)
            EXCEPTION("EMPTY");
        // First we create the faces of the element
        for (unsigned num_face = 0; num_face < iter->GetNumFaces(); ++num_face)
        {
            MonolayerVertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace = iter->GetFace(num_face);
            unsigned face_index = face_index_vector[num_face];

            for (unsigned num_triangle = 0; num_triangle < pFace->GetNumNodes(); ++num_triangle)
            {
                temporary_face_stream.push_back(3); // num of points
                // passive center
                temporary_face_stream.push_back(rMesh.GetNumNodes() + face_index);
                // counter clockwise boundary nodes
                temporary_face_stream.push_back(pFace->GetNodeGlobalIndex((num_triangle + 1) % pFace->GetNumNodes()));
                temporary_face_stream.push_back(pFace->GetNodeGlobalIndex((num_triangle + 0) % pFace->GetNumNodes()));
                num_faces++;
            }
        }

        // Now write the face stream, which has to be preallocated
        vtkIdList* face_stream = vtkIdList::New();
        face_stream->SetNumberOfIds(temporary_face_stream.size() + 1);
        face_stream->SetId(0, num_faces);
        for (unsigned face_stream_index = 0; face_stream_index < temporary_face_stream.size(); ++face_stream_index)
        {
            face_stream->SetId(face_stream_index + 1, temporary_face_stream[face_stream_index]);
        }
        // face_stream->PrintSelf(std::cout, vtkIndent(1));
        mpVtkUnstructedMesh->InsertNextCell(VTK_POLYHEDRON, face_stream);
        face_stream->Delete(); // Reference counted
    }
#endif // CHASTE_VTK
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddCellData(std::string dataName, std::vector<double> dataPayload)
{
#ifdef CHASTE_VTK
    vtkDoubleArray* p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i = 0; i < dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkCellData* p_cell_data = mpVtkUnstructedMesh->GetCellData();
    p_cell_data->AddArray(p_scalars);
    p_scalars->Delete(); // Reference counted
#endif // CHASTE_VTK
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddCellVectorData(std::string dataName,
                                                                          std::vector<c_vector<double, SPACE_DIM> > dataPayload)
{
#ifdef CHASTE_VTK
    vtkDoubleArray* p_vectors = vtkDoubleArray::New();
    p_vectors->SetName(dataName.c_str());
    p_vectors->SetNumberOfComponents(SPACE_DIM);
    for (unsigned i = 0; i < dataPayload.size(); i++)
    {
        p_vectors->InsertNextTuple3(dataPayload[i][0], dataPayload[i][1], dataPayload[i][2]);
    }

    vtkCellData* p_cell_data = mpVtkUnstructedMesh->GetCellData();
    p_cell_data->AddArray(p_vectors);
    p_vectors->Delete(); // Reference counted
#endif // CHASTE_VTK
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddFaceData(std::string dataName, std::vector<double> dataPayload)
{
#ifdef CHASTE_VTK
    vtkDoubleArray* p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i = 0; i < dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkCellData* p_face_data = mpVtkUnstructedMeshTriangulation->GetCellData();
    p_face_data->AddArray(p_scalars);
    p_scalars->Delete(); // Reference counted
#endif // CHASTE_VTK
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::AddPointData(std::string dataName, std::vector<double> dataPayload)
{
#ifdef CHASTE_VTK
    vtkDoubleArray* p_scalars = vtkDoubleArray::New();
    p_scalars->SetName(dataName.c_str());
    for (unsigned i = 0; i < dataPayload.size(); i++)
    {
        p_scalars->InsertNextValue(dataPayload[i]);
    }

    vtkPointData* p_point_data = mpVtkUnstructedMesh->GetPointData();
    p_point_data->AddArray(p_scalars);
    p_scalars->Delete(); // Reference counted
#endif // CHASTE_VTK
}

///\todo Mesh should be const (#1076)
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFilesUsingMesh(MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
{
    this->mpMeshReader = nullptr;
    mpMesh = &rMesh;

    this->mNumNodes = mpMesh->GetNumNodes();
    this->mNumElements = mpMesh->GetNumElements();
    this->mNumFaces = mpMesh->GetNumFaces();

    typedef typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator NodeIterType;
    mpIters->pNodeIter = new NodeIterType(mpMesh->GetNodeIteratorBegin());

    typedef typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementFaceIterator FaceIterType;
    mpIters->pFaceIter = new FaceIterType(mpMesh->GetFaceIteratorBegin());

    typedef typename MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>::MonolayerVertexElementIterator ElemIterType;

    mpIters->pElemIter = new ElemIterType(mpMesh->GetElementIteratorBegin());

    // Set up node map if we might have deleted nodes
    mNodeMapCurrentIndex = 0;
    mFaceMapCurrentIndex = 0;
    if (mpMesh->IsMeshChanging())
    {
        mpNodeMap = new NodeMap(mpMesh->GetNumAllNodes());
        for (NodeIterType it = mpMesh->GetNodeIteratorBegin();
             it != mpMesh->GetNodeIteratorEnd(); ++it)
        {
            mpNodeMap->SetNewIndex(it->GetIndex(), mNodeMapCurrentIndex++);
        }

        mpFaceMap = new NodeMap(mpMesh->GetNumFaces());
        for (FaceIterType it = mpMesh->GetFaceIteratorBegin();
             it != mpMesh->GetFaceIteratorEnd(); ++it)
        {
            mpFaceMap->SetNewIndex(it->GetIndex(), mFaceMapCurrentIndex++);
        }
    }
    WriteFiles();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void MonolayerVertexMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
    std::string comment = "# " + ChasteBuildInfo::GetProvenanceString();

    // Write node file
    std::string node_file_name = this->mBaseName + ".node";
    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name);

    // Write the node header
    // unsigned num_attr = 0;
    // unsigned max_bdy_marker = 1; // as we include boundary node information in the node file
    unsigned num_nodes = this->GetNumNodes();

    *p_node_file << num_nodes << "\t";
    *p_node_file << SPACE_DIM << "\n";
    // *p_node_file << num_attr << "\t";
    // *p_node_file << max_bdy_marker << "\n";

    *p_node_file << std::setprecision(6);

    // Write each node's data
    for (unsigned item_num = 0; item_num < num_nodes; item_num++)
    {
        std::vector<double> current_item = this->GetNextNode();
        *p_node_file << item_num;
        for (unsigned i = 0; i < SPACE_DIM + 1; i++)
        {
            *p_node_file << "\t" << current_item[i];
        }
        *p_node_file << "\n";
    }
    *p_node_file << comment << "\n";
    p_node_file->close();

    // Write face file
    std::string face_file_name = this->mBaseName + ".face";
    out_stream p_face_file = this->mpOutputFileHandler->OpenOutputFile(face_file_name);

    unsigned num_faces = mpMesh->GetNumFaces();
    *p_face_file << num_faces << "\n";

    for (unsigned item_num = 0; item_num < num_faces; item_num++)
    {
        MonolayerVertexFaceData face_data = this->GetNextMonolayerFace();

        // Get the node indices owned by this element
        std::vector<unsigned> node_indices = face_data.NodeIndices;
        std::vector<unsigned> node_types = face_data.NodeTypes;

        // Write this face's index and the number of nodes owned by it to file
        *p_face_file << item_num << "\t" << node_indices.size();

        // Write this face's type
        *p_face_file << "\t" << face_data.FaceType;

        // Write the node indices owned by this face to file
        for (unsigned i = 0; i < node_indices.size(); i++)
        {
            *p_face_file << "\t" << node_indices[i];
        }

        // Write the node types owned by this face to file
        for (unsigned i = 0; i < node_indices.size(); i++)
        {
            *p_face_file << "\t" << node_types[i];
        }

        // New line
        *p_face_file << "\n";
    }

    *p_face_file << comment << "\n";
    p_face_file->close();

    // Write element file
    std::string element_file_name = this->mBaseName + ".cell";
    out_stream p_element_file = this->mpOutputFileHandler->OpenOutputFile(element_file_name);

    // Write the element header
    unsigned num_elements = this->GetNumElements();
    *p_element_file << num_elements << "\n";

    // Write each element's data
    for (unsigned item_num = 0; item_num < num_elements; item_num++)
    {
        // Get data for this element
        MonolayerVertexElementData elem_data = this->GetNextMonolayerElement();

        // Get the node indices owned by this element
        std::vector<unsigned> face_indices = elem_data.FaceIndices;

        std::vector<bool> face_orientations = elem_data.FaceOrientations;

        assert(face_indices.size() == face_orientations.size());

        // Write this element's index and the number of nodes owned by it to file
        *p_element_file
            << item_num << "\t" << face_indices.size();

        // Write the face indices owned by this element to file
        for (unsigned i = 0; i < face_indices.size(); i++)
        {
            *p_element_file << "\t" << face_indices[i];
        }
        // Write the face orientations owned by this element to file
        for (unsigned i = 0; i < face_orientations.size(); i++)
        {
            *p_element_file << "\t" << face_orientations[i];
        }

        // New line
        *p_element_file << "\n";
    }

    *p_element_file << comment << "\n";
    p_element_file->close();
}

///////// Explicit instantiation///////

template class MonolayerVertexMeshWriter<1, 1>;
template class MonolayerVertexMeshWriter<1, 2>;
template class MonolayerVertexMeshWriter<1, 3>;
template class MonolayerVertexMeshWriter<2, 2>;
template class MonolayerVertexMeshWriter<2, 3>;
template class MonolayerVertexMeshWriter<3, 3>;