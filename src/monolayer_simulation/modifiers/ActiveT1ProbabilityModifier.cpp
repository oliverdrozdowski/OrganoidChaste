
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

#include "ActiveT1ProbabilityModifier.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MutableMonolayerVertexMesh.hpp"
#include "SimulationTime.hpp"
#include "Temperature.hpp"

template <unsigned DIM>
ActiveT1ProbabilityModifier<DIM>::ActiveT1ProbabilityModifier()
        : AbstractCellBasedSimulationModifier<DIM>(),
          mpTemperature(),
          mBoltzmannParameter(0.0),
          mActiveT1Rate(0.0),
          mIsActiveT1RatePerEdge(false)
{
}

template <unsigned DIM>
ActiveT1ProbabilityModifier<DIM>::~ActiveT1ProbabilityModifier()
{
}

template <unsigned DIM>
void ActiveT1ProbabilityModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    UpdateActiveT1Parameters(rCellPopulation);
}

template <unsigned DIM>
void ActiveT1ProbabilityModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */

    UpdateActiveT1Parameters(rCellPopulation);
}

template <unsigned DIM>
void ActiveT1ProbabilityModifier<DIM>::UpdateActiveT1Parameters(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    // Get the mesh
    MonolayerVertexBasedCellPopulation<DIM>* p_cell_population = static_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    MutableMonolayerVertexMesh<DIM, DIM>& r_mesh = p_cell_population->rGetMesh();

    assert(mBoltzmannParameter * mActiveT1Rate == 0.0);

    // Get current Temperature; Set to 1 if no Temperature object exists.
    double temp = 1.0;
    if (mpTemperature)
    {
        temp = mpTemperature->GetCurrentTemperature();
    }

    if (mBoltzmannParameter > 0.0)
    {
        double boltzmann_temp = mBoltzmannParameter * temp;
        r_mesh.SetActiveT1BoltzmannParameter(boltzmann_temp);
    }
    else if (mActiveT1Rate > 0.0)
    {
        // Rate means that the probability is rate * time_step
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        double time_step = p_simulation_time->GetTimeStep();
        double rate_temp = mActiveT1Rate * temp * time_step;
        r_mesh.SetActiveT1SwapProbability(rate_temp, mIsActiveT1RatePerEdge);
    }
}

template <unsigned DIM>
void ActiveT1ProbabilityModifier<DIM>::SetTemperature(boost::shared_ptr<Temperature> temperature)
{
    mpTemperature = temperature;
}

template <unsigned DIM>
void ActiveT1ProbabilityModifier<DIM>::SetActiveT1BoltzmannParameter(double boltzmann_parameter)
{
    // Check that the new value is greater than 0
    if (boltzmann_parameter < 0.0)
    {
        EXCEPTION("Attempting to assign a negative value.");
    }

    mBoltzmannParameter = boltzmann_parameter;
    mActiveT1Rate = 0.0;
}

template <unsigned DIM>
void ActiveT1ProbabilityModifier<DIM>::SetActiveT1Rate(double active_t1_rate, bool is_rate_per_edge)
{
    // Check that the new value is greater than 0
    if (active_t1_rate < 0.0)
    {
        EXCEPTION("Attempting to assign a negative value.");
    }

    mBoltzmannParameter = 0.0;
    mActiveT1Rate = active_t1_rate;
    mIsActiveT1RatePerEdge = is_rate_per_edge;
}

template <unsigned DIM>
void ActiveT1ProbabilityModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<BoltzmannParameter>" << mBoltzmannParameter << "</BoltzmannParameter>\n";
    *rParamsFile << "\t\t\t<ActiveT1Rate>" << mActiveT1Rate << "</ActiveT1Rate>\n";
    *rParamsFile << "\t\t\t<IsActiveT1RatePerEdge>" << mIsActiveT1RatePerEdge << "</IsActiveT1RatePerEdge>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ActiveT1ProbabilityModifier<1>;
template class ActiveT1ProbabilityModifier<2>;
template class ActiveT1ProbabilityModifier<3>;