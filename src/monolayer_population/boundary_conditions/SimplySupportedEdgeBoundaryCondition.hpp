#ifndef SIMPLYSUPPORTEDEDGEBOUNDARYCONDITION_HPP_
#define SIMPLYSUPPORTEDEDGEBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"

/**
 * A simply supported cell population boundary condition class, which fixes the movement
 * of apical and basal boundary nodes to specified planes in the domain. It thus implements
 * a simply supported boundary condition, where we, however, do not allow for a free height.
 * The height of the lateral face will be fixed by the inclination angle of the cell sheet.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class SimplySupportedEdgeBoundaryCondition : public AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>
{
private:
    /**
     * A point on the basal boundary plane.
     */
    c_vector<double, SPACE_DIM> mMidPlanePoint;

    /**
     * The outward-facing unit normal vector to the basal boundary plane.
     */
    c_vector<double, SPACE_DIM> mMidPlaneNormal;

    /**
     * Vector of nodes which are contrained in mid-plane
     */
    std::vector<Node<SPACE_DIM>*> mMidPlaneNodes;

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
    SimplySupportedEdgeBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation,
                                         c_vector<double, SPACE_DIM> mid_point,
                                         c_vector<double, SPACE_DIM> mid_normal);

    /**
     * @return #mBasalPointOnSubstrate.
     */
    const c_vector<double, SPACE_DIM>& rGetMidPlanePoint() const;

    /**
     * @return #mBasalNormalToSubstrate.
     */
    const c_vector<double, SPACE_DIM>& rGetMidPlaneNormal() const;

    /**
     * Sets the mMidPlaneNodes which are constrained in the mid-plane height
     * Use nullptr as first element if no nodes are to be constrained.
     * If this function is not called at all (or with empty vector) all nodes in the edges will be constrained.
     *
     * @param #mMidPlaneNodes
     */
    void SetMidPlaneNodes(std::vector<Node<SPACE_DIM>*> midPlaneNodes);

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
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SimplySupportedEdgeBoundaryCondition)

namespace boost
{
namespace serialization
{
    /**
     * Serialize information required to construct a SimplySupportedEdgeBoundaryCondition.
     */
    template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    inline void save_construct_data(
        Archive& ar, const SimplySupportedEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
    {
        // Save data required to construct instance
        const AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* const p_cell_population = t->GetCellPopulation();
        ar << p_cell_population;

        // Archive c_vectors one component at a time
        c_vector<double, SPACE_DIM> midPoint = t->rGetMidPlanePoint();
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar << midPoint[i];
        }
        c_vector<double, SPACE_DIM> midNormal = t->rGetMidPlaneNormal();
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar << midNormal[i];
        }
    }

    /**
     * De-serialize constructor parameters and initialize a SimplySupportedEdgeBoundaryCondition.
     */
    template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    inline void load_construct_data(
        Archive& ar, SimplySupportedEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
    {
        // Retrieve data from archive required to construct new instance
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population;
        ar >> p_cell_population;

        // Archive c_vectors one component at a time
        c_vector<double, SPACE_DIM> midPoint;
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar >> midPoint[i];
        }
        c_vector<double, SPACE_DIM> midNormal;
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            ar >> midNormal[i];
        }

        // Invoke inplace constructor to initialise instance
        ::new (t) SimplySupportedEdgeBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(p_cell_population, midPoint, midNormal);
    }
} // namespace serialization
} // namespace boost

#endif /*SIMPLYSUPPORTEDEDGEBOUNDARYCONDITION_HPP_*/
