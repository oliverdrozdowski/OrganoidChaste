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

#ifndef MONOLAYERVERTEXBASEDCELLPOPULATION_HPP_
#define MONOLAYERVERTEXBASEDCELLPOPULATION_HPP_

#include "AbstractCellPopulationWriter.hpp"
#include "AbstractFaceWriter.hpp"
#include "AbstractMonolayerVertexPopulationWriter.hpp"
#include "AbstractMovingBoundary.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "MutableMonolayerVertexMesh.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"

template <unsigned DIM>
class AbstractMonolayerVertexBasedDivisionRule; // Forward declaration to prevent circular include chain

template <unsigned DIM>
class GeometricalTargetVolumeModifier; // Forward declaration to prevent circular include chain

/**
 * A facade class encapsulating a vertex-based cell population.
 *
 * Contains a group of cells and maintains the associations
 * between CellPtrs and elements in the MutableVertexMesh.
 *
 */
template <unsigned DIM>
class MonolayerVertexBasedCellPopulation : public AbstractOffLatticeCellPopulation<DIM>
{
private:
    /**
     * This test uses the private constructor to simplify testing.
     */
    friend class TestVertexBasedDivisionRules;

    /**
     * Whether to delete the mesh when we are destroyed.
     * Needed if this cell population has been de-serialized.
     */
    bool mDeleteMesh;

    /**
     * A static cast of the AbstractMesh from AbstractCellPopulation
     * for use in this class
     */
    MutableMonolayerVertexMesh<DIM, DIM>* mpMutableVertexMesh;

    /** Whether to output the locations of T1 swaps and T3 swaps to files.
     * Defaults to true. */
    bool mOutputCellRearrangementLocations;

    /** A pointer to a division rule that is used to generate the axis when
     * dividing cells. This is a specialisation for Vertex models. */
    boost::shared_ptr<AbstractMonolayerVertexBasedDivisionRule<DIM> >
        mpVertexBasedDivisionRule;

    /**
     * Locations of T2 swaps (the centre of the removed triangle), stored so they
     * can be accessed and output by the cell killer and population writer
     * classes. The locations are stored until they are cleared by
     * ClearLocationsAndCellIdsOfT2Swaps().
     */
    std::vector<c_vector<double, DIM> > mLocationsOfT2Swaps;

    /**
     * The Ids of cells that have undergone T2 swaps, stored so they can be
     * accessed and output by the cell killer and population writer classes. The
     * Ids are stored until they are cleared by
     * ClearLocationsAndCellIdsOfT2Swaps().
     */
    std::vector<unsigned> mCellIdsOfT2Swaps;

    /**
     * Whether to restrict the vertex movement if vertex displacement is larger
     * than the cell rearrangement threshold.
     */
    bool mRestrictVertexMovement;

    /**
     * Whether to throw StepSizeExceptions, which defaults to true but is made
     * false after the first StepSizeException is thrown. In vertex based cell
     * populations a StepSizeException is not considered terminal, so there is no
     * need to throw more than one (as the numerical method uses WARN_ONCE_ONLY).
     */
    bool mThrowStepSizeException = true;

    /*
     * Whether to save the time of the last T1 transition or protorosette formation into the cells' data with
     * the respective volume after the transition. This is only necessary in the combination with
     * GeometricalTargetVolume Modifier, which only slowly adapts the target volume after transitions and
     * it is only set by the modifier's constructor
     */
    bool mSaveT1TimeAndVolumeToCells = false;

    /*
     * Whether the volumes of the cells should be relaxed towards the target
     * volume at the beginning of the simulation. For this the T1 time and initial
     * volume are adapted and the T1 relaxation mechanism is used. Therefore
     * mSaveT1TimeAndVolumeToCells should be set to true. This is especially
     * useful if one starts with randomized cell meshes and wants to get identical
     * volumes.
     */
    bool mDoInitialVolumeRelaxation = false;

    /**
     * Target Lumen Volume
     */
    double mTargetLumenVolume;

    /**
     * Moving boundaries
     */
    std::vector<boost::shared_ptr<AbstractMovingBoundary<DIM> > > mMovingBoundaries;

    /**
     * Overridden WriteVtkResultsToFile() method.
     *
     * @param rDirectory  pathname of the output directory, relative to where
     * Chaste output is stored
     */
    virtual void WriteVtkResultsToFile(const std::string& rDirectory);

    /** A list of face writers. */
    std::vector<boost::shared_ptr<AbstractFaceWriter<DIM, DIM> > > mFaceWriters;

    /** A list of monolayer population writer */
    std::vector<boost::shared_ptr<AbstractMonolayerVertexPopulationWriter<DIM, DIM> > > mMonolayerPopulationWriters;

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by
     * load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractOffLatticeCellPopulation<DIM> >(*this);
        archive & mOutputCellRearrangementLocations;
        archive & mFaceWriters;
        archive & mMonolayerPopulationWriters;
        archive & mpVertexBasedDivisionRule;
        archive & mRestrictVertexMovement;
    }

    /**
     * Check the consistency of internal data structures.
     * Each VertexElement must have a CellPtr associated with it.
     */
    void Validate();

    /**
     * Internal method to save the time of and volume after T1 transitions into
     * the cell data. This is needed by the GeometricalTargetVolumeModifier.
     * The data is only saved if mSaveT1TimeAndVolumeToCells has been set to true.
     */
    void SaveT1TimeAndVolumeToCells();

public:
    /**
     * Create a new cell population facade from a mesh and collection of cells.
     *
     * There must be precisely one CellPtr for each VertexElement in
     * the mesh.
     *
     * @param rMesh reference to a
     * @param rCells reference to a vector of CellPtrs
     * @param deleteMesh set to true if you want the cell population to free the
     * mesh memory on destruction
     * @param validate whether to validate the cell population when it is created
     * (defaults to true)
     * @param locationIndices an optional vector of location indices that
     * correspond to real cells
     */
    MonolayerVertexBasedCellPopulation(MutableMonolayerVertexMesh<DIM, DIM>& rMesh,
                                       std::vector<CellPtr>& rCells,
                                       bool deleteMesh = false,
                                       bool validate = true,
                                       const std::vector<unsigned> locationIndices = std::vector<unsigned>());

    /**
     * Constructor for use by boost serialization ONLY!
     *
     * @param rMesh a vertex mesh.
     */
    MonolayerVertexBasedCellPopulation(
        MutableMonolayerVertexMesh<DIM, DIM>& rMesh);

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~MonolayerVertexBasedCellPopulation();

    /**
     * Overridden GetDampingConstant() method.
     *
     * @param nodeIndex the global index of this node
     * @return the average damping constant of the cells surrounding the node.
     */
    double GetDampingConstant(unsigned nodeIndex);

    /**
     * @return reference to  mrMesh.
     */
    MutableMonolayerVertexMesh<DIM, DIM>& rGetMesh();

    /**
     * @return const reference to mrMesh (used in archiving).
     */
    const MutableMonolayerVertexMesh<DIM, DIM>& rGetMesh() const;

    /**
     * Get a particular VertexElement.
     *
     * @param elementIndex the global index of the VertexElement
     *
     * @return a pointer to the VertexElement.
     */
    MonolayerVertexElement<DIM, DIM>* GetElement(unsigned elementIndex);

    /**
     * @return the number of VertexElements in the cell population.
     */
    unsigned GetNumElements();

    /**
     * Overridden GetNumNodes() method.
     *
     * @return the number of nodes in the cell population.
     */
    unsigned GetNumNodes();

    /**
     * Overridden GetLocationOfCellCentre() method.
     *
     * Find the centre of mass of a given cell (assuming uniform density).
     * Note that, as there is no guarantee of convexity, this may lie
     * outside the VertexElement corresponding to the cell.
     *
     * @param pCell a cell in the population
     *
     * @return the location of the centre of mass of the element corresponding to
     * this cell.
     */
    c_vector<double, DIM> GetLocationOfCellCentre(CellPtr pCell);

    /**
     * Overridden GetNode() method.
     *
     * @param index global index of the specified node
     *
     * @return a pointer to the node.
     */
    Node<DIM>* GetNode(unsigned index);

    /**
     * Overridden GetNeighbouringLocationIndices() method.
     *
     * Given a cell, returns the set of location indices corresponding to
     * neighbouring cells.
     *
     * @param pCell a cell
     * @return the set of neighbouring location indices.
     */
    std::set<unsigned> GetNeighbouringLocationIndices(CellPtr pCell);

    /**
     * Overridden AddNode() method.
     *
     * Add a new node to the cell population.
     *
     * @param pNewNode pointer to the new node
     * @return global index of new node in cell population
     */
    unsigned AddNode(Node<DIM>* pNewNode);

    /**
     * Checks whether a given node displacement violates the movement threshold
     * for this population. If so, a stepSizeException is generated that contains
     * a warning/error message and a suggested smaller dt that should avoid the
     * problem.
     *
     * @param nodeIndex Index of the node in question (allows us to check whether
     * this is a ghost or particle)
     * @param rDisplacement Movement vector of the node at this time step
     * @param dt Current time step size
     */
    virtual void CheckForStepSizeException(unsigned nodeIndex, c_vector<double, DIM>& rDisplacement, double dt);

    /**
     * Overridden SetNode() method.
     *
     * Move the node with a given index to a new point in space.
     *
     * @param index the index of the node to be moved
     * @param rNewLocation the new target location of the node
     */
    void SetNode(unsigned index, ChastePoint<DIM>& rNewLocation);

    /**
     * Get a pointer to the element corresponding to a given CellPtr.
     *
     * @param pCell the cell
     *
     * @return pointer to the element.
     */
    MonolayerVertexElement<DIM, DIM>*
    GetElementCorrespondingToCell(CellPtr pCell);

    /**
     * Overridden AddCell() method.
     *
     * Add a new cell to the cell population.
     *
     * @param pNewCell  the cell to add
     * @param pParentCell pointer to a parent cell (if required)
     * @return address of cell as it appears in the cell list (internal of this
     * method uses a copy constructor along the way)
     */
    CellPtr AddCell(CellPtr pNewCell, CellPtr pParentCell = CellPtr());

    /**
     * Remove all cells labelled as dead.
     *
     * Note that after calling this method the cell population will be in an
     * inconsistent state until the equivalent of a 'remesh' is performed! So
     * don't try iterating over cells or anything like that.
     *
     * @return number of cells removed
     */
    unsigned RemoveDeadCells();

    /**
     * @return whether there is room into which a given cell may divide.
     * We check whether a node has a protorosette and if yes, we return false
     *
     * @param pCell pointer to a cell
     */
    virtual bool IsRoomToDivide(CellPtr pCell);

    /**
     * Overridden IsCellAssociatedWithADeletedLocation() method.
     *
     * @param pCell the cell
     * @return whether a given cell is associated with a deleted element.
     */
    bool IsCellAssociatedWithADeletedLocation(CellPtr pCell);

    /**
     * Remove the VertexElements which have been marked as deleted, perform
     * any cell rearrangements if required, and update the correspondence
     * with CellPtrs.
     *
     * @param hasHadBirthsOrDeaths - a bool saying whether cell population has had
     * Births Or Deaths not needed in this cell population class
     */
    void Update(bool hasHadBirthsOrDeaths = true);

    /**
     * Overridden OpenWritersFiles() method.
     *
     * Open all files in mCellPopulationWriters and mCellWriters for writing (not
     * appending).
     *
     * @param rOutputFileHandler handler for the directory in which to open this file.
     */
    virtual void OpenWritersFiles(OutputFileHandler& rOutputFileHandler);

    /**
     * A virtual method to accept a cell population writer so it can
     * write data from this object to file.
     *
     * @param pPopulationWriter the population writer.
     */
    virtual void AcceptPopulationWriter(
        boost::shared_ptr<AbstractCellPopulationWriter<DIM, DIM> >
            pPopulationWriter);

    /**
     * A virtual method to accept a cell population count writer so it can
     * write data from this object to file.
     *
     * @param pPopulationCountWriter the population count writer.
     */
    virtual void AcceptPopulationCountWriter(
        boost::shared_ptr<AbstractCellPopulationCountWriter<DIM, DIM> >
            pPopulationCountWriter);

    /**
     * A virtual method to accept a cell monolayer population writer so it can
     * write data from this object to file.
     *
     * @param pMonolayerPopulationWriter the monolayer population writer.
     */
    virtual void AcceptMonolayerPopulationWriter(boost::shared_ptr<AbstractMonolayerVertexPopulationWriter<DIM, DIM> > pMonolayerPopulationWriter);

    /**
     * A virtual method to accept a cell writer so it can
     * write data from this object to file.
     *
     * @param pCellWriter the population writer.
     * @param pCell the cell whose data are being written.
     */
    virtual void
    AcceptCellWriter(boost::shared_ptr<AbstractCellWriter<DIM, DIM> > pCellWriter,
                     CellPtr pCell);

    /**
     * Get the "rosette rank" of a cell.
     *
     * This is defined as the maximum number of cells shared by any node in the
     * cell's corresponding element.
     *
     * @param pCell boost shared pointer to a cell
     * @return rosette rank via associated mesh element
     */
    unsigned GetRosetteRankOfCell(CellPtr pCell);

    /**
     * Overridden GetVolumeOfCell() method.
     *
     * @param pCell boost shared pointer to a cell
     * @return volume via associated mesh element
     */
    double GetVolumeOfCell(CellPtr pCell);

    /**
     * Get thickness of cell, i.e. average distance apical and basal sides.
     *
     * @param pCell boost shared pointer to a cell
     * @return thickness via associated mesh element
     */
    double GetThicknessOfCell(CellPtr pCell);

    /**
     * Get force applied on a cell by distributing the force on each vertex onto
     * its elements
     */
    c_vector<double, DIM> GetAppliedForceOnCell(CellPtr pCell);

    /**
     * Add a a moving boundary.
     */
    void AddMovingBoundary(boost::shared_ptr<AbstractMovingBoundary<DIM> > boundary);

    /**
     * Add a a moving boundary.
     */
    void ClearMovingBoundary();

    /**
     * Get moving boundaries
     */

    std::vector<boost::shared_ptr<AbstractMovingBoundary<DIM> > > GetMovingBoundaries();

    /**
     * Return all locations of T2 swaps since the last sampling time step.
     *
     * @return mLocationsOfT2Swaps
     */
    std::vector<c_vector<double, DIM> > GetLocationsOfT2Swaps();

    /**
     * Return all IDs of cells undergoing T2 swaps since the last sampling time step.
     *
     * @return mLocationsOfT2Swaps
     */
    std::vector<unsigned> GetCellIdsOfT2Swaps();

    /**
     * Add location of T2 Swap.
     *
     * @param locationOfT2Swap  Location of the T2 swap
     */
    void AddLocationOfT2Swap(c_vector<double, DIM> locationOfT2Swap);

    /**
     * Add ID of T2 Swap.
     *
     * @param idOfT2Swap ID of the cell undergoing a T2 swap
     */
    void AddCellIdOfT2Swap(unsigned idOfT2Swap);

    /**
     * Clear the locations and IDs of cells undergoing T2 Swaps.
     */
    void ClearLocationsAndCellIdsOfT2Swaps();

    /**
     * @return mSaveT1TimeAndVolumeToCells
     */
    bool GetSaveT1TimeAndVolumeToCells();

    /**
     * @return mDoInitialVolumeRelaxation
     */
    bool GetDoInitialVolumeRelaxation();

    /**
     * @return mOutputCellRearrangementLocations
     */
    bool GetOutputCellRearrangementLocations();

    /**
     * Set mOutputCellRearrangementLocations.
     *
     * @param outputCellRearrangementLocations the new value of
     * mOutputCellRearrangementLocations
     */
    void SetOutputCellRearrangementLocations(bool outputCellRearrangementLocations);

    /**
     * Set mSaveT1TimeAndVolumeToCells.
     *
     * @param mSaveT1TimeAndVolumeToCells the new value of
     * mSaveT1TimeAndVolumeToCells, defaults to true
     */
    void SetSaveT1TimeAndVolumeToCells(bool saveT1TimeAndVolumeToCells = true);

    /**
     * Set mDoInitialVolumeRelaxation.
     *
     * @param mDoInitialVolumeRelaxation the new value of
     * mDoInitialVolumeRelaxation, defaults to true
     */
    void SetDoInitialVolumeRelaxation(bool doInitialVolumeRelaxation = true);

    /**
     * Set target lumen volume
     */
    void SetTargetLumenVolume(double vol);

    /**
     * Get target lumen volume
     */
    double GetTargetLumenVolume();

    /**
     * Initialize the times since last T1 and the volumes of last T1 with the
     * given value. This should be done to guarantee that the values are set,
     * because we assume so when writing the data to file.
     *
     * @param time_last_t1 the time since last t1 to set
     * @param volume_last_t1 the volume since last t1 to set
     */
    void InitializeT1TimeAndVolumeToCells(double time_last_t1, double volume_last_t1);

    /**
     * Overridden OutputCellPopulationParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationParameters(out_stream& rParamsFile);

    /**
     * Overridden GetWidth() method.
     *
     * Calculate the 'width' of any dimension of the cell population by calling
     * GetWidth() on the mesh.
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     */
    double GetWidth(const unsigned& rDimension);

    /**
     * Overridden GetNeighbouringNodeIndices() method.
     *
     * @param index the node index
     * @return the set of neighbouring node indices.
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned index);

    /**
     * Overridden GetTetrahedralMeshForPdeModifier() method.
     *
     * @return a pointer to a tetrahedral mesh using the nodes in the VertexMesh
     * as well as an additional node at the centre of each VertexElement.
     * At present, this method only works in 2D.
     *
     * This method is called by AbstractGrowingDomainPdeModifier.
     */
    virtual TetrahedralMesh<DIM, DIM>* GetTetrahedralMeshForPdeModifier();

    /**
     * Overridden IsPdeNodeAssociatedWithNonApoptoticCell() method.
     *
     * @param pdeNodeIndex inedx of a node in a tetrahedral mesh for use with a
     * PDE modifier
     *
     * @return if a node, specified by its index in a tetrahedral mesh for use
     *         with a PDE modifier, is associated with a non-apoptotic cell.
     * This method can be called by PDE classes.
     */
    virtual bool IsPdeNodeAssociatedWithNonApoptoticCell(unsigned pdeNodeIndex);

    /**
     * Overridden GetCellDataItemAtPdeNode() method.
     *
     * @param pdeNodeIndex index of a node in a tetrahedral mesh for use
     *         with a PDE modifier
     * @param rVariableName the name of the cell data item to get
     * @param dirichletBoundaryConditionApplies where a Dirichlet boundary
     * condition is used (optional; defaults to false)
     * @param dirichletBoundaryValue the value of the Dirichlet boundary
     * condition, if used (optional; defaults to 0.0)
     *
     * @return the value of a CellData item (interpolated if necessary) at a node,
     *         specified by its index in a tetrahedral mesh for use with a PDE
     * modifier. This method can be called by PDE modifier classes.
     */
    virtual double GetCellDataItemAtPdeNode(unsigned pdeNodeIndex,
                                            std::string& rVariableName,
                                            bool dirichletBoundaryConditionApplies = false,
                                            double dirichletBoundaryValue = 0.0);

    /**
     * @return The Vertex division rule that is currently being used.
     */
    boost::shared_ptr<AbstractMonolayerVertexBasedDivisionRule<DIM> >
    GetVertexBasedDivisionRule();

    /**
     * Set the division rule for this population.
     *
     * @param pVertexBasedDivisionRule  pointer to the new division rule
     */
    void SetVertexBasedDivisionRule(
        boost::shared_ptr<AbstractMonolayerVertexBasedDivisionRule<DIM> >
            pVertexBasedDivisionRule);

    /**
     * Overridden GetDefaultTimeStep() method.
     *
     * @return a default value for the time step to use
     * when simulating the cell population.
     *
     * A hard-coded value of 0.002 is returned. However, note that the time
     * step can be reset by calling SetDt() on the simulation object used to
     * simulate the cell population.
     */
    virtual double GetDefaultTimeStep();

    /**
     * Overridden WriteDataToVisualizerSetupFile() method.
     * Write any data necessary to a visualization setup file.
     * Used by AbstractCellBasedSimulation::WriteVisualizerSetupFile().
     *
     * @param pVizSetupFile a visualization setup file
     */
    virtual void WriteDataToVisualizerSetupFile(out_stream& pVizSetupFile);

    /**
     * Overridden SimulationSetupHook() method.
     *
     * Hook method to add a T2SwapCellKiller to a simulation object, which is
     * always required in the case of a VertexBasedCellPopulation. This
     * functionality avoids the need for static or dynamic casts to specific cell
     * population types within simulation methods.
     *
     * Note: In order to inhibit T2 swaps, the user needs to set the threshold for
     * T2 swaps in the MutableVertexMesh object mrMesh to 0, using the
     * SetT2Threshold() method.
     *
     * @param pSimulation pointer to a cell-based simulation object
     */
    virtual void
    SimulationSetupHook(AbstractCellBasedSimulation<DIM, DIM>* pSimulation);

    /**
     * Get the value of the mRestrictVertexMovement boolean.
     *
     * @return True if vertex movement is restricted at each timestep.
     */
    bool GetRestrictVertexMovementBoolean();

    /**
     * Set the value of the mRestrictVertexMovement boolean.
     *
     * @param restrictVertexMovement whether to restrict vertex movement in this
     * simulation.
     */
    void SetRestrictVertexMovementBoolean(bool restrictVertexMovement);

    /**
     * Add a face writer based on its type. Template parameters are inferred from the population.
     * The implementation of this function must be available in the header file.
     *
     * @return This method returns void
     */
    template <template <unsigned, unsigned> class T>
    void AddFaceWriter()
    {
        mFaceWriters.push_back(boost::shared_ptr<T<DIM, DIM> >(new T<DIM, DIM>));
    }

    /**
     * Add a face writer through an input argument.
     * This alternative to the templated AddFaceWriter()
     * method allows the user to, for example, add a writer
     * with a non-default value for its member mFileName.
     *
     * @param pFaceWriter shared pointer to a cell writer
     * @return This method returns void
     */
    void AddFaceWriter(boost::shared_ptr<AbstractFaceWriter<DIM, DIM> > pFaceWriter)
    {
        mFaceWriters.push_back(pFaceWriter);
    }

    /**
     * Get whether the population has a face writer.
     *
     * @return whether the population has this writer
     */
    bool HasFaceWriter() const
    {
        if (!mFaceWriters.empty())
        {
            return true;
        }
        return false;
    }

    /**
     * Close output files associated with any writers in the members
     * mFaceWriters.
     *
     */
    void CloseFaceWritersFiles();

    /**
     * Open all files in mFaceWriters in append mode for writing.
     *
     * @param rOutputFileHandler handler for the directory in which to open this
     * file.
     */
    void OpenRoundRobinFaceWritersFilesForAppend(
        OutputFileHandler& rOutputFileHandler);

    /**
     * Close all files in mFaceWriters.
     *
     */
    void CloseRoundRobinFaceWritersFiles();

    /**
     * Add a monolayer population writer based on its type. Template parameters
     * are inferred from the population. The implementation of this function must
     * be available in the header file.
     */
    template <template <unsigned, unsigned> class T>
    void AddMonolayerPopulationWriter()
    {
        mMonolayerPopulationWriters.push_back(boost::shared_ptr<T<DIM, DIM> >(new T<DIM, DIM>));
    }

    /**
     * Add a MonolayerPopulation writer through an input argument.
     * This alternative to the templated AddMonolayerPopulationWriter()
     * method allows the user to, for example, add a writer
     * with a non-default value for its member mFileName.
     *
     * @param pMonolayerPopulationWriter shared pointer to a cell writer
     */
    void AddMonolayerPopulationWriter(
        boost::shared_ptr<AbstractMonolayerVertexPopulationWriter<DIM, DIM> >
            pMonolayerPopulationWriter)
    {
        mMonolayerPopulationWriters.push_back(pMonolayerPopulationWriter);
    }

    /**
     * Get whether the population has monolayer population writer.
     *
     * @return whether the population has a writer
     */
    bool HasMonolayerPopulationWriter() const
    {
        if (!mMonolayerPopulationWriters.empty())
        {
            return true;
        }
        return false;
    }

    /**
     * write result to file. function copied from source code in order to
     * integrate monolayer population writers.
     */
    void WriteResultsToFiles(const std::string& rDirectory);

    /**
     * Compute the direction of the mean longest principal axes passing through the passive centers,
     * or 'long axes', of the apical and basal faces of a given element.
     * This is the eigenvector associated with the eigenvalue of smallest magnitude of the inertia
     * matrix for the faces with respect to the basis orthogonal to the mean normal.
     *
     * J = (  I_xx  -I_xy )
     *     ( -I_xy   I_yy )
     *
     * whose entries are computed by calling the method CalculateMomentsOfFace().
     * We consider the moments of the faces projected onto the space orthogonal to the mean normal.
     *
     * Note that if the nodes owned by the element are supplied in clockwise rather than
     * anticlockwise manner, or if this arises when any periodicity is enforced, then the
     * sign of each moment may be incorrect change. This means that we need to consider the eigenvalue
     * of largest magnitude rather than largest value when computing the short axis of the
     * element.
     *
     * If the element has a regular polygon face then the eigenvalues of the inertia tensor are
     * equal: in this case we consider a random unit vector for this face.
     *
     * This method is only implemented in 3D at present.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return a unit vector giving the direction of the long in-plane axis
     */
    c_vector<double, DIM> GetMeanApicalBasalLongAxisOfElement(unsigned index);

    /**
     * Get the max z position (of a vertex) of the population
     */
    double GetMaxZPosition();

    /**
     * Get the min z position (of a vertex) of the population
     */
    double GetMinZPosition();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MonolayerVertexBasedCellPopulation)

namespace boost
{
namespace serialization
{
    /**
     * Serialize information required to construct a VertexBasedCellPopulation.
     */
    template <class Archive, unsigned DIM>
    inline void save_construct_data(
        Archive& ar, const MonolayerVertexBasedCellPopulation<DIM>* t, const unsigned int file_version)
    {
        // Save data required to construct instance
        const MutableMonolayerVertexMesh<DIM, DIM>* p_mesh = &(t->rGetMesh());
        ar & p_mesh;
    }

    /**
     * De-serialize constructor parameters and initialise a VertexBasedCellPopulation.
     * Loads the mesh from separate files.
     */
    template <class Archive, unsigned DIM>
    inline void load_construct_data(
        Archive& ar, MonolayerVertexBasedCellPopulation<DIM>* t, const unsigned int file_version)
    {
        // Retrieve data from archive required to construct new instance
        MutableMonolayerVertexMesh<DIM, DIM>* p_mesh;
        ar >> p_mesh;

        // Invoke inplace constructor to initialise instance
        ::new (t) MonolayerVertexBasedCellPopulation<DIM>(*p_mesh);
    }
} // namespace serialization
} // namespace boost

#endif /*MONOLAYERVERTEXBASEDCELLPOPULATION_HPP_*/