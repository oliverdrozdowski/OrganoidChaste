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

#ifndef ACTIVET1PROBABILITYMODIFIER_HPP_
#define ACTIVET1PROBABILITYMODIFIER_HPP_

#include <boost/serialization/base_object.hpp>
#include "AbstractCellBasedSimulationModifier.hpp"
#include "ChasteSerialization.hpp"
#include "Temperature.hpp"

/**
 * A modifier class in which the T1 rate of of each cell in a MonolayerVertexBasedCellPopulation
 * is updated according to an underlying rule.
 *
 * The possibilities are:
 * (1) active T1 rate with a T1-rate mActiveT1Rate * temperature
 * (2) first-order boltzmann with a Boltzmann parameter in the mesh mBoltzmannParameter * temperature
 *
 */
template <unsigned DIM>
class ActiveT1ProbabilityModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
    /** Needed for serialization. */
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
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM, DIM> >(*this);
        archive & mpTemperature;
        archive & mBoltzmannParameter;
        archive & mActiveT1Rate;
        archive & mIsActiveT1RatePerEdge;
    }

protected:
    /**
     * A pointer to a Temperature object.
     */
    boost::shared_ptr<Temperature> mpTemperature;

    /**
     * The Boltzmann parameter for length dependent active T1 transformations.
     */
    double mBoltzmannParameter;

    /**
     * Constant Active T1 rate.
     */
    double mActiveT1Rate;

    /**
     * Flag if active T1 rate is global or per edge.
     */
    bool mIsActiveT1RatePerEdge;

public:
    /**
     * Default constructor.
     */
    ActiveT1ProbabilityModifier();

    /**
     * Destructor.
     */
    virtual ~ActiveT1ProbabilityModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to update active T1 rate and/or Boltzmann parameter.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateActiveT1Parameters(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Set the pointer to a Temperature object.
     */
    void SetTemperature(boost::shared_ptr<Temperature> temperature);

    /**
     * Set Boltzmann parameter for length dependent active T1 transformations; This also sets a possible constant active T1 rate to zero.
     */
    void SetActiveT1BoltzmannParameter(double boltzmann_parameter);

    /**
     * Set constant active T1 rate. Specify if it is per edge (default false). This also sets a possible length dependent Boltzmann parameter to zero.
     */
    void SetActiveT1Rate(double active_t1_rate, bool is_rate_per_edge = false);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(ActiveT1ProbabilityModifier)

#endif /*ACTIVET1PROBABILITYMODIFIER_HPP_*/