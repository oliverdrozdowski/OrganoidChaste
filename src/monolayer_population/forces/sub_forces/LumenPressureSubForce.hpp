/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef LUMENPRESSURESUBFORCE_HPP_
#define LUMENPRESSURESUBFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"
#include "Exception.hpp"

#include "AbstractCellBasedSimulation.hpp"
#include "AbstractSubForce.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

#include <iostream>

/**
 * A force class for use in MonolayerVertex-based simulations. This force is based on the
 * Energy function containing cell-specific surface tensions. Each face can have its own tension.
 * Therefore the correpsonding simulation modifier is needed if we have topological changes.
 */

template <unsigned DIM>
class LumenPressureSubForce : public AbstractSubForce<DIM>
{
    friend class TestForces;
    friend class TestMonolayerForces;

private:
    /**
     * Lumen pressure
     */
    double mPressureConstant;

    friend class boost::serialization::access;
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
        archive& boost::serialization::base_object<AbstractSubForce<DIM> >(*this);
        archive& mPressureConstant;
    }

public:
    /**
     * Constructor.
     */
    LumenPressureSubForce(double pressure_constant);

    /**
     * Destructor.
     */
    virtual ~LumenPressureSubForce();

    /**
     * Overridden GetForceContributions() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on given boundary
     *
     *
     * @param rCellPopulation reference to the cell population
     * @param pSimulation pointer to simulation instance
     *
     * @return the map of nodes to forces
     */
    NodeSubForceMap<DIM> GetForceContributions(MonolayerVertexBasedCellPopulation<DIM>& rCellPopulation,
                                               AbstractCellBasedSimulation<DIM, DIM>* pSimulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

    /**
     * Set the lumen pressure
     */
    void SetPressureConstant(double pressure_constant);

    /**
     * Get the lumen pressure constant
     */
    double GetPressureConstant();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LumenPressureSubForce)

#endif /*LUMENPRESSURESUBFORCE_HPP_*/