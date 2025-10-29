#ifndef STRAIGHTEDGEBOUNDARYCONDITION_HPP_
#define STRAIGHTEDGEBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"

/**
 * A straight edge boundary condition class, which fixes the movement
 * of apico-basal edge centers to variable planes in the domain.
 * The edge centers are constraint to a plane given by a normal vector and
 * the average position along the normal, to allow for global edge movement.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class StraightEdgeBoundaryCondition : public AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * The outward-facing unit normal vector to the boundary plane.
     */
    c_vector<double, SPACE_DIM> mEdgeNormal;

    /**
     * Vector of nodes which are contrained in mid-plane
     */
    std::vector<Node<SPACE_DIM>*> mEdgeNodes;

    /**
     * Whether to jiggle the cells on the basal bottom surface, initialised to false
     * in the constructor.
     */
    bool mUseJiggledNodesOnSubstrate;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM> >(*this);
        // archive & mUseJiggledNodesOnSubstrate;
    }

public:
    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param point a point on the boundary plane
     * @param normal the outward-facing unit normal vector to the boundary plane
     */
    StraightEdgeBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
                                  c_vector<double, SPACE_DIM> edge_normal);

    /**
     * @return #mBasalNormalToSubstrate.
     */
    const c_vector<double, SPACE_DIM>& rGetEdgeNormal() const;

    /**
     * Sets the mEdgeNodes which are constrained in the edge
     * Use nullptr as first element if no nodes are to be constrained.
     * If this function is not called at all (or with empty vector) all nodes in the edges will be constrained.
     *
     * @param #mEdgeNodes
     */
    void SetEdgeNodes(std::vector<Node<SPACE_DIM>*> edgeNodes);

    /**
     * Set method for mUseJiggledNodesOnSubstrate
     *
     * @param useJiggledNodesOnSubstrate whether to jiggle the nodes on the surface of the plane, can help stop overcrowding on plane.
     */
    void SetUseJiggledNodesOnSubstrate(bool useJiggledNodesOnSubstrate);

    /** @return #mUseJiggledNodesOnSubstrate. */
    bool GetUseJiggledNodesOnSubstrate();

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(StraightEdgeBoundaryCondition)

namespace boost
{
namespace serialization
{
    /**
     * Serialize information required to construct a StraightEdgeBoundaryCondition.
     */
    template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    inline void save_construct_data(
        Archive& ar, const StraightEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
    {
        // Save data required to construct instance
        const AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* const p_cell_population = t->GetCellPopulation();
        ar << p_cell_population;

        c_vector<double, SPACE_DIM> edgeNormal = t->rGetEdgeNormal();
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar << edgeNormal[i];
        }
    }

    /**
     * De-serialize constructor parameters and initialize a StraightEdgeBoundaryCondition.
     */
    template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    inline void load_construct_data(
        Archive& ar, StraightEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
    {
        // Retrieve data from archive required to construct new instance
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population;
        ar >> p_cell_population;

        // Archive c_vectors one component at a time
        c_vector<double, SPACE_DIM> edgeNormal;
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar >> edgeNormal[i];
        }

        // Invoke inplace constructor to initialise instance
        ::new (t) StraightEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(p_cell_population, edgeNormal);
    }
} // namespace serialization
} // namespace boost

#endif /*STRAIGHTEDGEBOUNDARYCONDITION_HPP_*/
