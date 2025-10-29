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

#ifndef FINITETHICKNESSSIMULATION3D_HPP_
#define FINITETHICKNESSSIMULATION3D_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include "MonolayerVertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"

/**
 * A 3D finite tickness monolayer simulation object.
 * This simulation object hast do be used if one wants to use the GeometricalTargetVolumeModifier
 * because we need to update the cell date after the ReMeshing procedure
 */
class FiniteThicknessSimulation3d : public OffLatticeSimulation<3>
{
    // Allow tests to access private members, in order to test computation of private functions e.g. DoCellBirth()
    friend class TestFiniteThicknessSimulation3d;

protected:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the simulation and member variable.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<OffLatticeSimulation<3> >(*this);
    }

    /**
     * Update the target volume after ReMesh in the population. This is necessary if we use the GeometricalTargetVolumeModifier
     * or the SmallAreaDampingModifier because the target volumes and the node-specific damping coefficients
     * need to be adapted before the forces are applied. This method is called by
     * the overwritten UpdateCellLocationsAndTopolgoy at the beginning and is therefore performed before it.
     */
    virtual void UpdateTargetVolumeAfterReMesh();

    /**
     * Update the simulation modifiers after births and deaths have occured in the UpdateCellPopulation method.
     * This is necessary if we use the GeometricalTargetVolumeModifier because the parameters have to be determined.
     * for newly created faces and cells. This method is called by the overwritten pdateCellPopulation at the end
     * and is therefore performed after it.
     *
     * @param birthsOrDeathsOccured wheter births/deaths occured in this UpdateCellPopulation step
     */
    virtual void UpdateModifiersAfterBirthsAndDeaths(bool birthsOrDeathsOccured);

    /**
     * Overridden SetupSolve() method.
     *
     *  Initialize the time and volume last T1 data in CellData if we use a GeometricalTargetVolumeModifier
     */
    virtual void SetupSolve();

    /**
     * Overridden UpdateCellLocationsAndTopology() method.
     *
     * Perform an UpdateTargetVolumeAfterReMesh at the beginning of this step.
     */
    virtual void UpdateCellLocationsAndTopology();

    /**
     * Overridden UpdateCellPopulation() method.
     *
     * Perform an UpdateModifiersAfterBirthsAndDeaths at the end of this step.
     */
    virtual void UpdateCellPopulation();

public:
    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells whether to initialise cells (defaults to true, set to false when loading from an archive)
     */
    FiniteThicknessSimulation3d(AbstractCellPopulation<3>& rCellPopulation,
                                bool deleteCellPopulationInDestructor = false,
                                bool initialiseCells = true);

    /**
     * Destructor.
     *
     * This frees the CryptSimulationBoundaryCondition.
     */
    virtual ~FiniteThicknessSimulation3d();

    /**
     * Outputs simulation parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(FiniteThicknessSimulation3d)

namespace boost
{
namespace serialization
{
    /**
     * Serialize information required to construct a CryptSimulation2d.
     */
    template <class Archive>
    inline void save_construct_data(
        Archive& ar, const FiniteThicknessSimulation3d* t, const unsigned int file_version)
    {
        // Save data required to construct instance
        const AbstractCellPopulation<3>* p_cell_population = &(t->rGetCellPopulation());
        ar & p_cell_population;
    }

    /**
     * De-serialize constructor parameters and initialise a CryptSimulation2d.
     */
    template <class Archive>
    inline void load_construct_data(
        Archive& ar, FiniteThicknessSimulation3d* t, const unsigned int file_version)
    {
        // Retrieve data from archive required to construct new instance
        AbstractCellPopulation<3>* p_cell_population;
        ar & p_cell_population;

        // Invoke inplace constructor to initialise instance
        ::new (t) FiniteThicknessSimulation3d(*p_cell_population, true, false);
    }
} // namespace serialization
} // namespace boost

#endif /*FINITETHICKNESSSIMULATION3D_HPP_*/