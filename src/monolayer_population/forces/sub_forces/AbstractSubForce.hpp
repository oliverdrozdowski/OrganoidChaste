#ifndef ABSTRACTSUBFORCE_HPP_
#define ABSTRACTSUBFORCE_HPP_

#include "AbstractForce.hpp"
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include "MonolayerVertexBasedCellPopulation.hpp"

/**
 * Short name for map from node pointer to force vector
 */
template <unsigned DIM>
using NodeSubForceMap = std::map<Node<DIM>*, c_vector<double, DIM> >;

/**
 * Short name for map from node to map from pointers to vectors
 */
template <unsigned DIM>
using NodeVolumeGradientMap = std::map<Node<DIM>*, std::map<MonolayerVertexElement<DIM, DIM>*, c_vector<double, DIM> > >;

/**
 * An abstract sub force class for use in the Generalized Volume Conserving Forces
 * We only use it for MonolayerVertexBasedCellPopulations.
 * As it is derived from an AbstractForce one may implement an alternative version
 * of AddForceContribution for other populations. For MonolayerVertexBasedCellPopulations
 * the usage as a force (not subforce in Generalized Volume Conserving Forces) leads to the
 * addition of the forces according to the NodeSubForceMap obtained from GetForceContributions.
 *
 * Note that the standard value for the simulation instance is 'nullptr' and the instance has to
 * be set if GetForceContribution needs the real pointer before running AddForceContribution.
 */
template <unsigned DIM>
class AbstractSubForce : public AbstractForce<DIM, DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractForce<DIM> >(*this);
    }

protected:
    /**
     * Pointer to the simulation instance. This can be set manually if we want to use it in implementations of
     * AddForceContribution in derived subclasses.
     */
    AbstractCellBasedSimulation<DIM, DIM>* mpSimulation = nullptr;

public:
    /**
     * Default constructor.
     */
    AbstractSubForce();

    /**
     * Destructor.
     */
    virtual ~AbstractSubForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the cell population based on the provided GetForceContributions
     * method. It is thus only implemented for MonolayerVertexBasedCellPopulations and must be overriden in
     * subclasses if other populations should be handled.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * @param pSimulation pointer to simulation instance
     */
    void SetSimulationInstance(AbstractCellBasedSimulation<DIM, DIM>* pSimulation);

    /**
     * @return mpSimulation pointer to simulation instance
     */
    AbstractCellBasedSimulation<DIM, DIM>* GetSimulationInstance();

    /**
     * Calculates the force on each node and returns a map of nodes to force vectors
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rCellPopulation reference to the cell population
     * @param pSimulation pointer to the simulation instance, if left out is set to nullptr.
     *
     * @return the map from nodes to forces
     */
    virtual NodeSubForceMap<DIM>
    GetForceContributions(MonolayerVertexBasedCellPopulation<DIM>& rCellPopulation,
                          AbstractCellBasedSimulation<DIM, DIM>* pSimulation)
        = 0;

    /**
     * Outputs force parameters to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile) = 0;
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractSubForce)

#endif /*ABSTRACTSUBFORCE_HPP_*/
