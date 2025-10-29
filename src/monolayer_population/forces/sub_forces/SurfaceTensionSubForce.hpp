#ifndef SURFACETENSIONSUBFORCE_HPP_
#define SURFACETENSIONSUBFORCE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"
#include "Exception.hpp"

#include "AbstractCellBasedSimulation.hpp"
#include "AbstractSubForce.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"

#include <iostream>

// Forward declaration prevents circular include chain
struct DampingCoefficientPair;

/**
 * A force class for use in MonolayerVertex-based simulations. This force is based on the
 * Energy function containing cell-specific surface tensions. Each face can have its own tension.
 * Therefore the correpsonding simulation modifier is needed if we have topological changes.
 */

template <unsigned DIM>
class SurfaceTensionSubForce : public AbstractSubForce<DIM>
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
        archive& boost::serialization::base_object<AbstractSubForce<DIM> >(*this);
        archive & mSurfaceTensionParameters;
        archive & mPerformActiveT1Swaps;
        archive & mAnnealingDecayTime;
        archive & mInitialTemperature;
        archive & mApicalStandardTension;
        archive & mBasalStandardTension;
        archive & mLateralStandardTension;
        archive & mT1AreTemperatureDependent;
        archive & mT1TransitionRate;
        archive & mModeOfSurfaceTensionCalculation;
        archive & mMappingMutationTension;
        archive & mpSimulation;
    }

protected:
    /**
     * The strengths of the apical surface tension term in the model. A vector corresponding to faces
     */
    std::vector<double> mSurfaceTensionParameters;
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

    // Dictionary for mutation state-dependent tensions
    std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3> > mMappingMutationTension;

    AbstractCellBasedSimulation<DIM, DIM>* mpSimulation;

public:
    /**
     * Constructor.
     */
    SurfaceTensionSubForce();

    /**
     * Destructor.
     */
    virtual ~SurfaceTensionSubForce();

    /**
     * Overridden GetForceContributions() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on the energy function
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
     * SetSurfaceTensionParametersByMutation(). This function does nothing if the mode how
     * surface tensions are calculated is set to matrix (=2) or is not set at all (=0).
     *
     * @param pCellPopulation pointer to the cell population to update
     */
    void UpdateSurfaceTensions(AbstractCellPopulation<DIM>* pCellPopulation);

    /**
     * Calculate the total energy based on the surface tensions. Energy of a face is
     * tension times area
     *
     * @return total energy
     **/
    double CalculateTotalEnergy(AbstractCellPopulation<DIM>& rCellPopulation);

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
    void SetPerformActiveT1Swaps(bool doActiveT1Swaps = true);

    /**
     * Set parameters for Metropolis/simulated annealing random forces.
     * Negative annealing decay time will lead to a growth of the temperature prefactor
     * according to the equation T=T_0 - t/annealingDecayTime, i.e. we can fully
     * control the linear factor but use a different equation!
     *
     * @param relativeNoiseStrength the relative strength of the noise
     * @param annealingDecayTime temperature decreases linearly T=T_0*(1-t/annealingDecayTime)
     * @param initialTemperature the initial temperature for the Metropolis probability
     */
    void SetSimulatedAnnealingParameters(double relativeNoiseStrength,
                                         double annealingDecayTime,
                                         double initialTemperature = 10.0);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SurfaceTensionSubForce)

#endif /*SURFACETENSIONSUBFORCE_HPP_*/