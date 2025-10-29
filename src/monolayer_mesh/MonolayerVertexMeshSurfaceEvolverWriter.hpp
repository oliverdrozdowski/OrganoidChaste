#ifndef MONOLAYERVERTEXMESHSURFACEEVOLVERWRITER_HPP_
#define MONOLAYERVERTEXMESHSURFACEEVOLVERWRITER_HPP_

// Forward declaration prevents circular include chain
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonolayerVertexMesh;

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct MonolayerMeshSurfaceEvolverWriterIterators;

#include "AbstractMeshWriter.hpp"
#include "MonolayerVertexMesh.hpp"
#include "MonolayerVertexMeshWriter.hpp"
#include "NodeMap.hpp"
#include "SurfaceTensionSubForce.hpp"

// Forward declaration prevents circular include chain
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonolayerVertexMesh;

/**
 * A mesh writer class that writes a Monolayer Vertex Mesh into a file which can be used in Surface Evolver
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class MonolayerVertexMeshSurfaceEvolverWriter : public AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * If writing from a mesh object, the mesh to write to disk.
     * Otherwise NULL.
     */
    MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>* mpMesh;

    /** Iterators over the mesh */
    MonolayerMeshSurfaceEvolverWriterIterators<ELEMENT_DIM, SPACE_DIM>* mpIters;

    /** Track deleted nodes so they don't get written */
    NodeMap* mpNodeMap;
    /** What was the last index written to #mpNodeMap ? */
    unsigned mNodeMapCurrentIndex;

    /** Total number of faces in mesh */
    unsigned mNumFaces;

    /** Vector of oriented edges which we save as pairs of indices */
    std::vector<std::array<unsigned, 2> > mEdgeIndexPairs;

    /** Vector containing the information if edge is on boundary (both vertices boundary) */
    std::vector<bool> mEdgeOnBoundary;

    /** Map of oriented edges to index in mEdgeIndexPairs */
    std::map<std::array<unsigned, 2>, unsigned> mMapEdgeToIndex;

    /** Map of face index to vector of (oriented) faces as index in mEdgeIndexPairs */
    std::map<unsigned, std::vector<int> > mMapFaceIndexToVectorOfEdges;

    /** Pointer to a surface tension force, if we want to also include surface tensions */
    boost::shared_ptr<SurfaceTensionSubForce<SPACE_DIM> > mpSurfaceTensionSubForce;

    /** Whether we should add a small random number to the volumes of the cells as target volumes */
    bool mUseRandomizedVolumes;

    /** Upper boudnary for the uniform distribution from which we draw the extra volumes if mUseRandomizedVolumes */
    double mUniformDistributionBoundary = 1.0;

    /** Wheter we write face type (apical/basal/lateral) as extra facet attribute into evolver file */
    bool mWriteFaceTypeIntoFile = false;

    /** Map from surface tension to colors (ints). Note, only 15 colors are available in Surface Evovler! */
    std::map<double, int> mMapTensionToColor;

    /** Whether the boundary nodes of the mesh are set fixed in surface evolver */
    bool mFixBoundaryNodes = false;

    /** Whether the lumen (volume enclosed by all apical faces) is implemented as a cell/body */
    bool mConsiderLumenAsCell = false;

public:
    /**
     * Constructor.
     *
     * @param rDirectory reference to the output directory, relative to where Chaste output is stored
     * @param rBaseName reference to the base name for results files
     * @param clearOutputDir whether to clear the output directory prior to writing files
     */
    MonolayerVertexMeshSurfaceEvolverWriter(const std::string& rDirectory,
                                            const std::string& rBaseName,
                                            const bool clearOutputDir = true);

    /**
     * Destructor.
     */
    ~MonolayerVertexMeshSurfaceEvolverWriter();

    /**
     * Write files using a mesh.
     *
     * @param rMesh reference to the vertex-based mesh
     */
    void WriteFilesUsingMesh(MonolayerVertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

    /**
     * Set the pointer to surface tension force, if we want to print
     * the surface tensions to file
     *
     * @param p_force pointer to the force
     */
    void SetSurfaceTensionSubForce(boost::shared_ptr<SurfaceTensionSubForce<SPACE_DIM> > p_force);

    /**
     * Set the map of surface tensions to colors, which are written into the file
     * for surface evolver. Note that only 15 colors are available in Surface Evolver
     *
     * @param map std::map from tensions (double) to colors (int)
     */
    void SetMapTensionToColor(std::map<double, int> map);

    /**
     * Whether the current node in the iterator is a boundary node
     *
     * @return bool
     */
    bool IsNodeOnBoundary();

    /**
     * Get the color of the face. If mMapTensionToColor is not empty, we use the mapping.
     * If the surface tension is not in the map, we return 0. If the map is empty, we return
     * the color associated to the facetype.
     *
     * @param index index to face
     */
    int GetColorOfFace(unsigned index);

    /**
     * Get the surface tension for face with index
     *
     * @param index index to face
     */
    double GetSurfaceTensionOfFace(unsigned index);

    /**
     * Get the face type of face with index as string
     *
     * @param index index to face
     */
    unsigned GetFaceTypeOfFace(unsigned index);

    /**
     * Set mUseRandomizedVolumes
     *
     * @param use_rand_vol boolean to set
     */
    void SetUseRandomizedVolumes(bool useRandomVol, double distributionBoundary);

    /**
     *
     * @return mUseRandomizedVolumes
     */
    bool GetUseRandomizedVolumes();

    /**
     * Set mFixBoundaryNodes
     */
    void SetFixBoundaryNodes(bool fixBoundaries);

    /**
     * Set mWriteFaceTypeIntoFile
     */
    void SetWriteFaceTypeIntoFile(bool writeFaceType = true);

    /**
     * Set mConsiderLumenAsCell
     */
    void SetConsiderLumenAsCell(bool considerLumenAsCell = true);

    /**
     * @return the coordinates of the next node to be written to file
     */
    std::vector<double> GetNextNode();

    /**
     * Populate mEdgeIndexPairs, mMapEdgeToIndex and mMapFaceIndexToVectorOfEdges.
     * This method creates the mappings from edges to nodes and from faces to edges.
     * Called by WriteFilesUsingMesh.
     */
    void CreateEdgeMapping();

    /**
     * Write mesh data to files.
     * This method must be overridden in concrete classes.
     */
    void WriteFiles();
};

#endif /*MONOLAYERVERTEXMESHSURFACEEVOLVERWRITER_HPP_*/