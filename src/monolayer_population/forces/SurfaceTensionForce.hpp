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

#ifndef SURFACETENSIONFORCE_HPP_
#define SURFACETENSIONFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"
#include "Exception.hpp"

#include "AbstractCellBasedSimulation.hpp"
#include "AbstractForce.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "SmallAreaDampingModifier.hpp"

#include <iostream>

// Forward declaration prevents circular include chain
struct DampingCoefficientPair;

/**
 * A force class for use in MonolayerVertex-based simulations. This force is based on the
 * Energy function containing cell-specific surface tensions. Each face can have its own tension.
 * Therefore the correpsonding simulation modifier is needed if we have topological changes.
 */

template <unsigned DIM>
class SurfaceTensionForce : public AbstractForce<DIM>
{
    friend class TestForces;
    friend class TestMonolayerForces;

private:
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
        archive& boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mSurfaceTensionParameters;
        archive & mRelativeNoiseStrength;
        archive & mAnnealingDecayTime;
        archive & mInitialTemperature;
        archive & mPerformActiveT1Swaps;
        archive & mpSimulation;
    }

    double CalculateTotalEnergy(AbstractCellPopulation<DIM>& rCellPopulation);

protected:
    /**
     * The strengths of the apical surface tension term in the model. A vector corresponding to faces
     */
    std::vector<double> mSurfaceTensionParameters;
    double mVolumeElasticityParameter;
    bool mProjectOnVolumeConservingForces;
    bool mPerformActiveT1Swaps;

    // Parameters for the simulated annealing
    double mRelativeNoiseStrength;
    double mAnnealingDecayTime;
    double mInitialTemperature;

    // Surface tensions which can be saved into the file, if they are set regularly
    double mApicalStandardTension;
    double mBasalStandardTension;
    double mLateralStandardTension;

    // Parameters for active T1 transitions
    bool mT1AreTemperatureDependent;
    double mT1TransitionRate;

    /*
     * The mode how surface tensions are calculated, is set by functions, which write the tensions
     * 0=not assigned, 1=uniform standard values, 2=matrix, 3=dictionary for mutation states
     */
    unsigned mModeOfSurfaceTensionCalculation;

    /*
     * Whether the different nodes have different mobilities. Is set by methods writing the corresponding map
     */
    bool mIndividualNodeDampingActivated;

    /*
     * Map of node pointers to damping thresholds which should be applied to the nodes
     */
    std::map<Node<DIM>*, DampingCoefficientPair> mMapNodeToDampingThresholds;

    // Dictionary for mutation state-dependent tensions
    std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3> > mMappingMutationTension;

    AbstractCellBasedSimulation<DIM, DIM>* mpSimulation;

public:
    /**
     * Constructor.
     */
    SurfaceTensionForce();

    /**
     * Destructor.
     */
    virtual ~SurfaceTensionForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on the energy function
     *
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Get the surface tension parameter for the face.
     *
     * @param face
     *
     * @return the surface tension parameter for this face.
     */
    virtual double GetSurfaceTensionParameter(unsigned face);

    /**
     * @return mSurfaceTensionParameters
     */
    std::vector<double> GetSurfaceTensionParameters();

    /**
     * @return mApicalStandardTension
     */
    double GetApicalStandardTension();

    /**
     * @return mBasalStandardTension
     */
    double GetBasalStandardTension();

    /**
     * @return mLateralStandardTension
     */
    double GetLateralStandardTension();

    /**
     * Do we project on volume conserving forces?
     *
     * @return mProjectOnVolumeConservingForces
     */
    bool DoProjectOnVolumeConservingForces();

    /**
     * Do we perform active T1 swaps (to be used with simulatead annealing)?
     *
     * @return mPerformActiveT1Swaps
     */
    bool DoPerformActiveT1Swaps();

    /**
     * @return mVolumeElasticityParameters
     */
    double GetVolumeElasticityParameter();

    /**
     * Set surface tension of face to value of surfaceTensionParameter. This also sets the mode how
     * surface tensions are calculated to matrix (=2). Note that this mode cannot be automatically updated!
     *
     * @param face the index of face
     * @param surfaceTensionParameter the new value of mSurfaceTensionParameter at face
     */
    void SetSurfaceTensionParameter(unsigned face, double surfaceTensionParameter);

    /**
     * Set vector mSurfaceTensionParameters. This also sets the mode how surface tensions are calculated
     * to matrix (=2). Note that this mode cannot be automatically updated!
     *
     * @param surfaceTensionParameters the vector of SurfaceTensionParameters
     */
    void SetSurfaceTensionParameters(std::vector<double> surfaceTensionParameters);

    /**
     * Set apical/basal/lateral surface tensions for each cell based on the mutation state
     * For this we supply a map from mutation state to an arrray of surface tensions
     * This also sets the mode how surface tensions are calculated to dictionary (=3)
     * The array is organized as: apical, basal, lateral.
     *
     * @param mappingMutationTension the map from mutation to tension array
     */
    void SetSurfaceTensionParametersByMutation(std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3> > mappingMutationTension);

    /**
     * Update the surface tensions by standard values. This also sets the mode how surface tensions are
     * calculated to uniform standard values (=1)
     *
     * @param pCellPopulation pointer to cell population
     */
    void CreateSurfaceTensionParametersForCells(double apical_surface_tension,
                                                double basal_surface_tension,
                                                double lateral_surface_tension,
                                                MonolayerVertexMesh<DIM, DIM>* p_mesh);

    /**
     * Update the surface tensions by their mutation state. Prior to this SetSurfaceTensionParametersByMutation()
     * must have been called as the last function to set the rule for surface tensions.
     *
     * @param pCellPopulation pointer to cell population
     */
    void UpdateSurfaceTensionsByMutation(AbstractCellPopulation<DIM>* pCellPopulation);

    /**
     * Update the surface tension parameters, based on the rules previously supplied,
     * i.e., the value of mModeOfSurfaceTensionCalculation, prior to calculating the forces
     * This function either calls CreateSurfaceTensionParametersForCells() or
     * SetSurfaceTensionParametersByMutation().	This function does nothing if the mode how
     * surface tensions are calculated is set to matrix (=2) or is not set at all (=0).
     *
     * @param pCellPopulation	pointer to the cell population to update
     */
    void UpdateSurfaceTensions(AbstractCellPopulation<DIM>* pCellPopulation);

    /**
     * Set the parameters for active T1 transitions. The rate is computed as rate*temp*timestep/numfaces
     *
     * @param T1TransitionRate rate for T1 transitions
     * @param T1AreTemperatureDependent whether we multiply the rate with the temperature (defaults to true)
     */
    void SetT1TransitionParameters(double T1TransitionRate, bool T1AreTemperatureDependent = true);

    /**
     * @return mT1TransitionRate
     */
    double GetT1TransitionRate();

    /**
     * @return mT1AreTemperatureDependent
     */
    bool IsT1TransitonRateTemperatureDependent();

    /**
     * Set if we project on volume conserving forces.
     *
     * @param doProjection whether to do the projection
     */
    void SetProjectOnVolumeConservingForces(bool doProjection);

    /**
     * Set if we project on volume conserving forces.
     *
     * @param doProjection whether to do the projection
     */
    void SetPerformActiveT1Swaps(bool doActiveT1Swaps = true);

    /**
     * Set parameters for Metropolis/simulated annealing random forces.
     *
     * @param relativeNoiseStrength the relative strength of the noise
     * @param annealingDecayTime temperature decreases as T=T_0*exp(-t/annealingDecayTime)
     * @param initialTemperature the initial temperature for the Metropolis probability
     */
    void SetSimulatedAnnealingParameters(double relativeNoiseStrength,
                                         double annealingDecayTime,
                                         double initialTemperature = 10.0);

    /**
     * Set pointer to Simulation instance, which is needed for simulated annealing
     *
     * @param pSimulation pointer to simulation instance
     */
    void SetSimulationInstance(AbstractCellBasedSimulation<DIM, DIM>* pSimulation);

    /**
     * @param VolumeElasticityParameter
     */
    void SetVolumeElasticityParameter(double volumeElasticityParameter);

    /**
     * @param MapNodeToDampingThresholds
     */
    void SetMapNodeToDampingThresholds(std::map<Node<DIM>*, DampingCoefficientPair> mapNodeToDampingThresholds);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SurfaceTensionForce)

#endif /*SURFACETENSIONFORCE_HPP_*/