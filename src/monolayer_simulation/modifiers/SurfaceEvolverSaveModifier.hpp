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

#ifndef SURFACEEVOLVERSAVEMODIFIER_HPP_
#define SURFACEEVOLVERSAVEMODIFIER_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include "AbstractCellBasedSimulationModifier.hpp"
#include "MonolayerVertexMeshSurfaceEvolverWriter.hpp"
#include "SurfaceTensionSubForce.hpp"

/**
 * A modifier class which at each simulation time step calculates the volume of each cell
 * and stores it in in the CellData property as "volume". To be used in conjunction with
 * contact inhibition cell cycle models.
 */
template <unsigned DIM>
class SurfaceEvolverSaveModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;

    /** Whether we save the surface evolver file at every timestep at which we save data*/
    bool mSaveEveryTime = false;

    /** Whether we save the surface evolver file at the end after the simulation has finished */
    bool mSaveAtEnd = false;

    /** Output directory path */
    std::string mOutputDirectory;

    /** Pointer to surface tension force object */
    boost::shared_ptr<SurfaceTensionSubForce<DIM> > mpSurfaceTensionSubForce = nullptr;

    /** Map from surface tension to colors (ints). Note, only 15 colors are available in Surface Evovler! */
    std::map<double, int> mMapTensionToColor;

    /** Whether we should add a small random number to the volumes of the cells as target volumes */
    bool mUseRandomizedVolumes = false;

    /** Upper boudnary for the uniform distribution from which we draw the extra volumes if mUseRandomizedVolumes */
    double mUniformDistributionBoundary = 1.0;

    /** Whether the boundary nodes of the mesh are set fixed in surface evolver */
    bool mFixBoundaryNodes = false;

    /** Wheter we write face type (apical/basal/lateral) as extra facet attribute into evolver file */
    bool mWriteFaceTypeIntoFile = false;

    /** Whether the lumen (volume enclosed by all apical faces) is implemented as a cell/body */
    bool mConsiderLumenAsCell = false;

    /**
     * Create a MonolayerVertexMeshSurfaceEvolverWriter, change all settings accordingly and return a pointer.
     * This function is called by UpdateAtEndOfOutputTimeStep() and UpdateAtEndOfSolve().
     * Note that the object has to be deleted manually (in these functions)
     *
     * @param filename filename of the file to be saved
     *
     * @return MonolayerVertexMeshSurfaceEvolverWriter<DIM,DIM>* pointer to writer
     */
    MonolayerVertexMeshSurfaceEvolverWriter<DIM, DIM>* CreateSurfaceEvolverWriter(std::string filename);

    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM, DIM> >(*this);
        archive & mSaveEveryTime;
        archive & mSaveAtEnd;
        archive & mOutputDirectory;
        archive & mpSurfaceTensionSubForce;
        archive & mMapTensionToColor;
        archive & mUseRandomizedVolumes;
        archive & mUniformDistributionBoundary;
        archive & mFixBoundaryNodes;
        archive & mWriteFaceTypeIntoFile;
        archive & mConsiderLumenAsCell;
    }

public:
    /**
     * Default constructor.
     */
    SurfaceEvolverSaveModifier();

    /**
     * Destructor.
     */
    virtual ~SurfaceEvolverSaveModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     * Nothing happens.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Overridden UpdateAtEndOfOutputTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step after the output into files.
     * We save the surface evolver file if mSaveEveryTime is true
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     * We save the surface evolver file if mSaveAtEnd is true
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     * The outputDirectory string is saved to mOutputDirectory for subsequent file generation.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Whether we save a file for Surface Evolver at every timestep
     *
     * @param saveEveryTime boolean for mSaveEveryTime
     */
    void SetSaveEveryTime(bool saveEveryTime = true);

    /**
     * Whether we save a file for Surface Evolver after the end of the simulation
     *
     * @param saveAtEnd boolean for mSaveAtEnd
     */
    void SetSaveAtEnd(bool saveAtEnd = true);

    /**
     * Set mUseRandomizedVolumes
     *
     * @param use_rand_vol boolean to set
     */
    void SetUseRandomizedVolumes(bool useRandomVol, double distributionBoundary);

    /**
     * Set mFixBoundaryNodes
     */
    void SetFixBoundaryNodes(bool fixBoundaries = true);

    /**
     * Set mWriteFaceTypeIntoFile
     */
    void SetWriteFaceTypeIntoFile(bool writeFaceType = true);

    /**
     * Set mConsiderLumenAsCell
     */
    void SetConsiderLumenAsCell(bool considerLumenAsCell = true);

    /**
     *
     * @return mUseRandomizedVolumes
     */
    bool GetUseRandomizedVolumes();

    /**
     * Set the internal pointer to the surface tension force object to get the tensions and write them to file
     *
     * @param pForce pointer to force
     */
    void SetSurfaceTensionSubForce(boost::shared_ptr<SurfaceTensionSubForce<DIM> > pForce);

    /**
     * Set the map of surface tensions to colors, which are written into the file
     * for surface evolver. Note that only 15 colors are available in Surface Evolver
     *
     * @param map std::map from tensions (double) to colors (int)
     */
    void SetMapTensionToColor(std::map<double, int> map);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SurfaceEvolverSaveModifier)

#endif /*SURFACEEVOLVERSAVEMODIFIER_HPP_*/
